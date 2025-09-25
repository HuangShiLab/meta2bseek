use crate::cmdline::{ContainArgs, ProfileArgs};
use anyhow::{Result, anyhow, Context};
use std::collections::HashMap;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use rayon::prelude::*;
use std::sync::Mutex;
use std::sync::Arc;
use std::collections::HashSet;
use std::path::PathBuf;
use crate::constants::Hash;
use std::thread;
use std::time::Duration;
use memory_stats::memory_stats;

pub use crate::extract::{SyldbEntry, SylspEntry};


// 内存监控和限制功能 - 采用 sylph 的真正实现
// pub fn check_vram_and_block(max_ram: usize, file: &str) {
//     if let Some(usage) = memory_stats() {
//         let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
//         if (max_ram as f64) < gb_usage_curr {
//             eprintln!("Max memory reached. Blocking processing for {}. Curr memory {:.2} GB, max mem {} GB",
//                      file, gb_usage_curr, max_ram);
//         }
//         while (max_ram as f64) < gb_usage_curr {
//             let five_second = Duration::from_secs(1);
//             thread::sleep(five_second);
//             if let Some(usage) = memory_stats() {
//                 gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
//                 if (max_ram as f64) >= gb_usage_curr {
//                     eprintln!("Processing for {} freed", file);
//                 }
//             } else {
//                 break;
//             }
//         }
//     }
// }

// 定义分类学信息结构体
#[derive(Debug, Clone, Default)]
pub struct TaxonomyInfo {
    pub kingdom: String,    // d__
    pub phylum: String,     // p__
    pub class: String,      // c__
    pub order: String,      // o__
    pub family: String,     // f__
    pub genus: String,      // g__
    pub species: String,    // s__
}

impl TaxonomyInfo {
    pub fn from_gtdb_string(gtdb_str: &str) -> Result<Self> {
        let mut taxonomy = TaxonomyInfo::default();
        
        for part in gtdb_str.split(';') {
            let part = part.trim();
            if part.starts_with("d__") {
                taxonomy.kingdom = part[3..].to_string();
            } else if part.starts_with("p__") {
                taxonomy.phylum = part[3..].to_string();
            } else if part.starts_with("c__") {
                taxonomy.class = part[3..].to_string();
            } else if part.starts_with("o__") {
                taxonomy.order = part[3..].to_string();
            } else if part.starts_with("f__") {
                taxonomy.family = part[3..].to_string();
            } else if part.starts_with("g__") {
                taxonomy.genus = part[3..].to_string();
            } else if part.starts_with("s__") {
                taxonomy.species = part[3..].to_string();
            }
        }
        
        Ok(taxonomy)
    }
    
    pub fn get_species_key(&self) -> String {
        format!("{}|{}|{}|{}|{}|{}|{}", 
                self.kingdom, self.phylum, self.class, 
                self.order, self.family, self.genus, self.species)
    }
}

// 物种级别的丰度结果 - 使用 Arc 共享 TaxonomyInfo
#[derive(Debug, Clone)]
pub struct SpeciesAbundanceResult {
    pub taxonomy: Arc<TaxonomyInfo>,
    pub sample_abundances: FxHashMap<String, f64>,  // sample_id -> abundance
    pub total_tags: usize,
    pub genome_count: usize,  // 该物种包含的genome数量
    pub reads_count: usize,   // 该物种的总 reads 数 (S)
    pub gscore: f64,          // G-score = sqrt(reads_count * total_tags)
}

// 定义比对结果结构
#[derive(Debug, Clone)]
pub struct QueryResult {
    pub sample_file: String,
    pub genome_file: String,
    pub adjusted_ani: f64,
    pub eff_cov: f64,
    pub ani_percentile: (f64, f64),
    pub eff_lambda: f64,
    pub lambda_percentile: (f64, f64),
    pub median_cov: f64,
    pub mean_cov_geq1: f64,
    pub containment_ind: String,
    pub naive_ani: f64,
    pub contig_name: String,
    pub ref_tags: usize,
    pub shared_tags: usize,
    pub query_tags: usize,
    pub taxonomic_abundance: f64,
    pub sequence_abundance: f64,
}

// 新增基因组级别的结果结构体
#[derive(Debug, Clone)]
pub struct GenomeProfileResult {
    pub genome_id: String,
    pub sample_id: String,  // 这个字段将存储实际的样本来源（如 sample1, sample4 等）
    pub file_path: String,  // 新增字段，存储文件路径
    pub adjusted_ani: f64,
    pub taxonomic_abundance: f64,
    pub sequence_abundance: f64,
    pub common_tags: usize,
    pub total_tags: usize,
    pub eff_cov: f64,
}

// 新增：k-mer重新分配相关的结构体和函数

// Winner table条目结构 - 对应sylph中的 (f64, &'a GenomeSketch, bool)
#[derive(Debug, Clone)]
struct WinnerTableEntry {
    pub ani: f64,
    pub genome_id: String,
    pub was_reassigned: bool,
}

// 重新分配统计信息
#[derive(Debug, Clone)]
struct ReassignmentStats {
    pub tags_lost: usize,
    pub total_tags: usize,
    pub reassignment_ratio: f64,
}

// 构建winner table - 借鉴sylph的高效实现
fn build_winner_table<'a>(
    results: &'a [QueryResult], 
    db_entries: &'a [SyldbEntry],
    log_reassign: bool
) -> FxHashMap<Hash, WinnerTableEntry> {
    eprintln!("Building winner table for {} results and {} database entries", results.len(), db_entries.len());
    
    let mut tag_to_genome_map: FxHashMap<Hash, WinnerTableEntry> = FxHashMap::default();
    
    // 关键优化1：构建sequence_id到db_entry的直接映射，避免O(G)查找
    let mut seq_id_to_entry: FxHashMap<String, &SyldbEntry> = FxHashMap::default();
    for entry in db_entries {
        seq_id_to_entry.insert(entry.sequence_id.clone(), entry);
    }
    
    eprintln!("Built seq_id_to_entry mapping with {} entries", seq_id_to_entry.len());
    
    // 关键优化2：借鉴sylph的直接遍历方式，避免复杂的嵌套查找
    for res in results.iter() {
        // O(1)查找，而不是O(G)的find操作
        if let Some(db_entry) = seq_id_to_entry.get(&res.contig_name) {
            // 减少字符串操作：预先计算genome_id
            let genome_id = extract_genome_id_from_path(&db_entry.genome_source);
            
            // 借鉴sylph的简洁遍历方式
            for tag in &db_entry.tags {
                let entry = tag_to_genome_map.entry(*tag).or_insert_with(|| {
                    WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_id: genome_id.to_string(),
                        was_reassigned: false,
                    }
                });
                
                // 关键：选择ANI最高的基因组作为该tag的"赢家"
                if res.adjusted_ani > entry.ani {
                    *entry = WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_id: genome_id.to_string(),
                        was_reassigned: true,
                    };
                }
            }
        } else {
            eprintln!("Warning: No database entry found for contig {}", res.contig_name);
        }
    }
    
    eprintln!("Winner table built with {} unique tags", tag_to_genome_map.len());
    
    // 记录重新分配日志（借鉴sylph的简洁实现）
    if log_reassign {
        eprintln!("------------- Logging tag reassignments -----------------");
        let mut genome_to_index: FxHashMap<String, usize> = FxHashMap::default();
        for (i, res) in results.iter().enumerate() {
            eprintln!("Index\t{}\t{}\t{}", i, res.genome_file, res.contig_name);
            genome_to_index.insert(res.genome_file.clone(), i);
        }
        
        // 关键优化3：借鉴sylph的简洁并行计算，避免复杂的映射查找
        (0..results.len()).into_par_iter().for_each(|i| {
            let res = &results[i];
            let mut reassign_edge_map: FxHashMap<(usize, usize), usize> = FxHashMap::default();
            
            // 使用优化后的seq_id_to_entry映射
            if let Some(db_entry) = seq_id_to_entry.get(&res.contig_name) {
                for tag in &db_entry.tags {
                    if let Some(winner_entry) = tag_to_genome_map.get(tag) {
                        if winner_entry.genome_id != res.genome_file {
                            if let Some(&winner_index) = genome_to_index.get(&winner_entry.genome_id) {
                                let edge_count = reassign_edge_map.entry((winner_index, i)).or_insert(0);
                                *edge_count += 1;
                            }
                        }
                    }
                }
            }
            
            // 直接输出，避免收集到向量中
            for ((from_idx, to_idx), count) in reassign_edge_map {
                if count > 10 {
                    eprintln!("{}->{}\t{}\ttags reassigned", from_idx, to_idx, count);
                }
            }
        });
    }
    
    // 添加调试信息：统计重新分配的情况
    let mut reassigned_tags = 0;
    let mut total_tags = 0;
    for entry in tag_to_genome_map.values() {
        total_tags += 1;
        if entry.was_reassigned {
            reassigned_tags += 1;
        }
    }
    eprintln!("Reassignment statistics: {}/{} tags were reassigned ({:.2}%)", 
              reassigned_tags, total_tags, 
              if total_tags > 0 { reassigned_tags as f64 / total_tags as f64 * 100.0 } else { 0.0 });
    
    // 检查tags重叠情况
    let mut tag_counts: FxHashMap<Hash, usize> = FxHashMap::default();
    for db_entry in db_entries {
        for tag in &db_entry.tags {
            *tag_counts.entry(*tag).or_insert(0) += 1;
        }
    }
    
    let overlapping_tags = tag_counts.values().filter(|&&count| count > 1).count();
    let total_unique_tags = tag_counts.len();
    eprintln!("Tag overlap analysis: {}/{} unique tags are shared between genomes ({:.2}%)", 
              overlapping_tags, total_unique_tags,
              if total_unique_tags > 0 { overlapping_tags as f64 / total_unique_tags as f64 * 100.0 } else { 0.0 });
    
    if overlapping_tags > 0 {
        eprintln!("Shared tags found! This should enable reassignment.");
    } else {
        eprintln!("No shared tags found. Reassignment cannot work without overlapping tags.");
    }
    
    tag_to_genome_map
}

// 使用winner table重新计算统计结果 - 模仿sylph的get_stats函数
fn recalculate_with_winner_table(
    db_entries: &[SyldbEntry],
    sample_entries: &[SylspEntry],
    winner_map: &FxHashMap<Hash, WinnerTableEntry>,
    min_ani: f64,
    log_reassign: bool
) -> Vec<QueryResult> {
    eprintln!("Recalculating with winner table for {} database entries and {} sample entries", 
              db_entries.len(), sample_entries.len());
    
    // 按样本源分组样本条目
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in sample_entries {
        sample_groups.entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }
    
    let mut all_results = Vec::new();
    
    // 为每个样本源分别计算
    for (sample_source, group_entries) in sample_groups {
        let sample_tags: HashSet<Hash> = group_entries.iter()
            .map(|entry| entry.tag.clone())
            .collect();
        
        let total_sample_tags = group_entries.len();
        
        // 为每个数据库条目计算重新分配后的结果
        for db_entry in db_entries {
            // 最小标签数过滤
            if db_entry.tags.len() < MIN_TAGS_FOR_GENOME {
                continue;
            }
            
            let mut owned_tags = 0;
            let mut tags_lost_count = 0;
            
            // 计算属于该基因组的tags（模仿sylph的重新分配逻辑）
            for tag in &db_entry.tags {
                if sample_tags.contains(tag) {
                    if let Some(winner_entry) = winner_map.get(tag) {
                        if winner_entry.genome_id == extract_genome_id_from_path(&db_entry.genome_source) {
                            // 该tag属于当前基因组
                            owned_tags += 1;
                        } else {
                            // 该tag被重新分配给其他基因组
                            tags_lost_count += 1;
                        }
                    } else {
                        // 该tag没有被任何基因组"拥有"（理论上不应该发生）
                        owned_tags += 1;
                    }
                }
            }
            
            let total_ref_tags = db_entry.tags.len();
            
            if log_reassign && tags_lost_count > 0 {
                eprintln!("Genome {} in sample {}: owned_tags={}, lost_tags={}, total_ref_tags={}", 
                         extract_genome_id_from_path(&db_entry.genome_source), 
                         sample_source, owned_tags, tags_lost_count, total_ref_tags);
            }
            
            // 计算统计数据
            let mut result = calculate_statistics(
                owned_tags,
                total_sample_tags,
                total_ref_tags,
            );
            
            // 设置基本信息 - 关键修复：使用正确的样本源
            result.sample_file = sample_source.clone();
            result.genome_file = extract_genome_id_from_path(&db_entry.genome_source).to_string();
            result.contig_name = db_entry.sequence_id.clone();
            result.shared_tags = owned_tags;
            result.query_tags = total_sample_tags;
            result.ref_tags = total_ref_tags;
            
            // 计算平均深度和覆盖度
            if owned_tags > 0 {
                result.mean_cov_geq1 = 1.0;
                result.eff_cov = owned_tags as f64 / total_ref_tags as f64;
                result.median_cov = 1.0;
            }
            
            // 应用profile专用的过滤条件
            if filter_results_for_profile(&result, Some(min_ani)) {
                all_results.push(result);
            }
        }
    }
    
    eprintln!("Recalculation completed: {} results after filtering", all_results.len());
    all_results
}

// 过滤过度重新分配的基因组 - 完全模仿sylph的derep_if_reassign_threshold函数
fn filter_over_reassigned_genomes(
    results_old: &[QueryResult],
    results_new: &[QueryResult],
    ani_thresh: f64,
    k: f64
) -> Vec<QueryResult> {
    eprintln!("Filtering over-reassigned genomes: old_results={}, new_results={}, ani_thresh={}, k={}", 
              results_old.len(), results_new.len(), ani_thresh, k);
    
    let ani_thresh = ani_thresh / 100.0;
    
    // 构建genome_id到旧结果的映射
    let mut genome_to_old_result: FxHashMap<String, &QueryResult> = FxHashMap::default();
    for result in results_old.iter() {
        genome_to_old_result.insert(result.genome_file.clone(), result);
    }
    
    let threshold = f64::powf(ani_thresh, k);
    let mut return_vec = Vec::new();
    let mut filtered_count = 0;
    
    for result in results_new.iter() {
        if let Some(old_res) = genome_to_old_result.get(&result.genome_file) {
            let num_tag_reassign = (old_res.shared_tags - result.shared_tags) as f64;
            let reass_thresh = threshold * result.ref_tags as f64;
            
            if num_tag_reassign < reass_thresh {
                return_vec.push(result.clone());
            } else {
                eprintln!("genome {} had num tags reassigned = {}, threshold was {}, removing.", 
                         result.genome_file, num_tag_reassign, reass_thresh);
                filtered_count += 1;
            }
        } else {
            // 如果没有找到旧结果，保留新结果
            return_vec.push(result.clone());
        }
    }
    
    eprintln!("Filtering completed: {} genomes filtered out, {} genomes retained", 
              filtered_count, return_vec.len());
    
    return_vec
}

// 计算 G-score 的函数
// 输入：某个物种的 reads 数 S（所有 2bRAD markers 的总 reads）和测得的 tag 数目 t
// 输出：G-score = sqrt(S * t)
fn calculate_gscore(reads_count: usize, tag_count: usize) -> f64 {
    ((reads_count as f64) * (tag_count as f64)).sqrt()
}

// 基于 G-score 过滤的函数
// 输入：物种列表（包含 S 和 t 信息）以及一个外部指定的阈值 gscore_threshold
// 输出：过滤后的物种列表，只保留 G-score >= gscore_threshold 的物种
fn filter_species_by_gscore(
    species_results: &mut Vec<SpeciesAbundanceResult>,
    gscore_threshold: f64
) -> Vec<SpeciesAbundanceResult> {
    eprintln!("Filtering species by G-score threshold: {:.2}", gscore_threshold);
    
    let initial_count = species_results.len();
    
    // 首先计算每个物种的 G-score
    for species in species_results.iter_mut() {
        species.gscore = calculate_gscore(species.reads_count, species.total_tags);
    }
    
    // 过滤 G-score >= gscore_threshold 的物种
    let filtered_results: Vec<SpeciesAbundanceResult> = species_results
        .iter()
        .filter(|species| species.gscore >= gscore_threshold)
        .cloned()
        .collect();
    
    let filtered_count = filtered_results.len();
    let removed_count = initial_count - filtered_count;
    
    eprintln!("G-score filtering results: {} species retained, {} species removed (threshold: {:.2})", 
              filtered_count, removed_count, gscore_threshold);
    
    if removed_count > 0 {
        eprintln!("Removed species had G-scores below {:.2}", gscore_threshold);
    }
    
    filtered_results
}

// 重新计算丰度 - 模仿sylph的丰度计算逻辑
fn recalculate_abundances_after_reassignment(
    results: &mut [QueryResult],
    sample_entries: &[SylspEntry]
) {
    eprintln!("Recalculating abundances for {} results", results.len());
    
    // 计算总覆盖度
    let total_cov: f64 = results.iter()
        .map(|r| r.eff_cov)
        .sum();
    
    let total_seq_cov: f64 = results.iter()
        .map(|r| r.eff_cov * r.ref_tags as f64)
        .sum();
    
    eprintln!("Total coverage: {:.6}, Total sequence coverage: {:.6}", total_cov, total_seq_cov);
    
    // 重新计算每个结果的丰度
    for result in results.iter_mut() {
        if total_cov > 0.0 {
            result.taxonomic_abundance = result.eff_cov / total_cov * 100.0;
        } else {
            result.taxonomic_abundance = 0.0;
        }
        
        if total_seq_cov > 0.0 {
            result.sequence_abundance = result.eff_cov * result.ref_tags as f64 / total_seq_cov * 100.0;
        } else {
            result.sequence_abundance = 0.0;
        }
    }
    
    eprintln!("Abundance recalculation completed");
}

// 修改常量
const MIN_COVERAGE: f64 = 0.001;  // 最小覆盖度要求
const MIN_ANI: f64 = 90.0;      // 最小 ANI 要求，调整为与sylph一致
const MIN_SHARED_TAGS: usize = 10;  // 最小共享标签数
const K: f64 = 31.0;  // k-mer 长度
const LAMBDA_THRESHOLD: f64 = 0.05;  // lambda 阈值
const MIN_TAGS_FOR_GENOME: usize = 50;  // 基因组最小标签数要求，调整为与sylph一致
const PROFILE_MIN_ANI: f64 = 95.0;  // profile模式下的最小ANI阈值，调整为与sylph一致
const PROFILE_MIN_COVERAGE: f64 = 0.005;  // profile模式下的最小覆盖率

struct MultiWriter {
    writers: Vec<Box<dyn Write + Send>>,
}

impl MultiWriter {
    fn new() -> Self {
        MultiWriter { writers: Vec::new() }
    }
    fn add_writer(&mut self, writer: Box<dyn Write + Send>) {
        self.writers.push(writer);
    }
}

impl Write for MultiWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        for w in &mut self.writers {
            w.write_all(buf)?;
        }
        Ok(buf.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        for w in &mut self.writers {
            w.flush()?;
        }
        Ok(())
    }
}

pub fn query(args: ContainArgs) -> Result<()> {
    // 首先测试文件格式
    let db_files: Vec<_> = args.files.iter()
        .filter(|f| f.ends_with(".syldb"))
        .collect();
    
    let sample_files: Vec<_> = args.files.iter()
        .filter(|f| f.ends_with(".sylsp"))
        .collect();

    if db_files.is_empty() {
        return Err(anyhow!("No .syldb files found in input files"));
    }

    if sample_files.is_empty() {
        return Err(anyhow!("No .sylsp files found in input files"));
    }

    // 创建输出写入器
    let writer = Arc::new(Mutex::new(create_multi_writer(&args.out_file_name)?));

    // 打印表头（只打印一次）
    print_header(&writer)?;

    // 并行处理所有数据库文件
    for db_path in db_files {
        eprintln!("Processing database file: {}", db_path);
        
        // 读取数据库文件
        let db_file = File::open(db_path)
            .with_context(|| format!("Failed to open database file: {}", db_path))?;
        let db_reader = BufReader::new(db_file);
        let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
            .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;

        eprintln!("Found {} entries in database", db_entries.len());

        // 并行处理所有样本文件
        sample_files.par_iter().try_for_each(|sample_path| -> Result<()> {
            eprintln!("Processing sample file: {}", sample_path);
            
            let sample_file = File::open(sample_path)
                .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
            let sample_reader = BufReader::new(sample_file);
            let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
                .with_context(|| format!("Failed to deserialize sample file: {}", sample_path))?;

            eprintln!("Found {} entries in sample", sample_entries.len());

            // 检查样本数据的有效性
            if sample_entries.is_empty() {
                eprintln!("Warning: Sample {} has no tags", sample_path);
                return Ok(());
            }

                    // 构建样本标签的哈希表
        let sample_tags: HashMap<Hash, usize> = sample_entries.iter()
            .map(|entry| (entry.tag.clone(), 1))
            .collect();

            let total_sample_tags = sample_entries.len();
            eprintln!("Total unique tags in sample: {}", total_sample_tags);

            // 对每个基因组记录进行比对
            for db_entry in &db_entries {
                // 计算共享标签和统计信息
                let mut shared_tags = 0;
                let mut coverages = Vec::new();
                let total_ref_tags = db_entry.tags.len();

                for tag in &db_entry.tags {
                    if sample_tags.contains_key(tag) {
                        shared_tags += 1;
                        coverages.push(1.0); // 简化的覆盖度计算
                    }
                }

                eprintln!("Found {} shared tags between sample and reference {}", 
                         shared_tags, db_entry.sequence_id);

                // 计算统计数据
                let mut result = calculate_statistics(
                    shared_tags,
                    total_sample_tags,
                    total_ref_tags,
                );

                // 设置基本信息
                result.sample_file = sample_path.to_string();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();
                result.shared_tags = shared_tags;
                result.query_tags = total_sample_tags;
                result.ref_tags = total_ref_tags;

                // 计算平均深度和覆盖度
                if shared_tags > 0 {
                    result.mean_cov_geq1 = 1.0; // 简化的深度计算
                    result.eff_cov = shared_tags as f64 / total_ref_tags as f64;
                    
                    // 计算中位数覆盖度
                    if !coverages.is_empty() {
                        coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());
                        result.median_cov = if coverages.len() % 2 == 0 {
                            (coverages[coverages.len()/2 - 1] + coverages[coverages.len()/2]) / 2.0
                        } else {
                            coverages[coverages.len()/2]
                        };
                    }
                }

                // 应用过滤条件
                if filter_results(&result, args.minimum_ani) {
                    eprintln!("Result passed filters: ANI={:.2}, Coverage={:.3}", 
                            result.adjusted_ani, result.eff_cov);
                    // 输出结果
                    print_result(&result, &writer)?;
                } else {
                    eprintln!("Result filtered out: ANI={:.2}, Coverage={:.3}", 
                            result.adjusted_ani, result.eff_cov);
                }
            }
            Ok(())
        })?;
    }

    Ok(())
}

fn create_multi_writer(out_file_name: &Option<String>) -> Result<Box<dyn Write + Send>> {
    let mut mw = MultiWriter::new();
    mw.add_writer(Box::new(BufWriter::new(std::io::stdout())));
    if let Some(path) = out_file_name {
        let file = File::create(path)
            .with_context(|| format!("Failed to create output file: {}", path))?;
        mw.add_writer(Box::new(BufWriter::new(file)));
    }
    Ok(Box::new(mw))
}

fn print_header(writer: &Arc<Mutex<Box<dyn Write + Send>>>) -> Result<()> {
    let mut writer = writer.lock().unwrap();
    writeln!(writer, "{:<20} {:<20} {:<10} {:<10} {:<15} {:<15} {:<10} {:<10} {:<10} {:<15} {:<10} {:<10}",
        "Sample_file", "Genome_file", "ANI(%)", "Eff_cov", "ANI_5-95%", "Eff_lambda", "Lambda_5-95%", "Median_cov", "Mean_cov", "Containment", "Naive_ANI", "Contig_name")?;
    writeln!(writer, "{:-<150}", "")?;
    Ok(())
}

fn print_result(result: &QueryResult, writer: &Arc<Mutex<Box<dyn Write + Send>>>) -> Result<()> {
    let mut writer = writer.lock().unwrap();
    writeln!(writer, "{:<20} {:<20} {:<10.2} {:<10.3} {:<7.2}-{:<7.2} {:<10.3} {:<7.2}-{:<7.2} {:<10.3} {:<10.3} {:<7} {:<10.2} {:<10}",
        result.sample_file,
        result.genome_file,
        result.adjusted_ani,
        result.eff_cov,
        result.ani_percentile.0,
        result.ani_percentile.1,
        result.eff_lambda,
        result.lambda_percentile.0,
        result.lambda_percentile.1,
        result.median_cov,
        result.mean_cov_geq1,
        result.containment_ind,
        result.naive_ani,
        result.contig_name
    )?;
    Ok(())
}

fn calculate_statistics(shared_tags: usize, query_tags: usize, total_ref_tags: usize) -> QueryResult {
    // 避免除零错误
    if query_tags == 0 || total_ref_tags == 0 {
        return QueryResult {
            sample_file: String::new(),
            genome_file: String::new(),
            contig_name: String::new(),
            adjusted_ani: 0.0,
            eff_cov: 0.0,
            ani_percentile: (0.0, 0.0),
            eff_lambda: 0.0,
            lambda_percentile: (0.0, 0.0),
            median_cov: 0.0,
            mean_cov_geq1: 0.0,
            containment_ind: format!("{}/{}", shared_tags, total_ref_tags),
            naive_ani: 0.0,
            ref_tags: total_ref_tags,
            shared_tags: 0,
            query_tags: 0,
            taxonomic_abundance: 0.0,
            sequence_abundance: 0.0,
        };
    }

    // 使用 f64 进行所有计算
    let shared_tags_f64 = shared_tags as f64;
    let total_ref_tags_f64 = total_ref_tags as f64;

    // 计算基础 ANI（使用 sylph 的方法）
    let containment_ratio = shared_tags_f64 / total_ref_tags_f64;
    
    // 只有当共享标签数大于最小要求时才计算 ANI
    let (naive_ani, adjusted_ani) = if shared_tags >= MIN_SHARED_TAGS {
        let naive = f64::powf(containment_ratio, 1.0 / K) * 100.0;
        
        // 计算调整后的 ANI
        let coverage_factor = if containment_ratio < 0.1 {
            1.0 + (0.1 - containment_ratio) * 0.5
        } else {
            1.0
        };
        let adjusted = (naive * coverage_factor).min(100.0);
        (naive, adjusted)
    } else {
        // 当共享标签数太少时，ANI 应该很低但不一定是 0
        let base_ani = (shared_tags_f64 / MIN_SHARED_TAGS as f64) * 80.0; // 使用 80% 作为基准
        (base_ani, base_ani)
    };
    
    // 计算有效覆盖度
    let eff_cov = containment_ratio;
    
    // 计算 Lambda 值
    let eff_lambda = if eff_cov < LAMBDA_THRESHOLD {
        eff_cov * 1.2
    } else {
        eff_cov
    };

    // 计算置信区间
    let base_uncertainty = 1.0;
    let coverage_uncertainty = (1.0 - eff_cov) * 1.5;
    let total_uncertainty = base_uncertainty + coverage_uncertainty;
    
    let ani_low = (adjusted_ani - total_uncertainty).max(0.0);
    let ani_high = (adjusted_ani + total_uncertainty).min(100.0);
    
    // Lambda 置信区间
    let lambda_uncertainty = 0.02 + (1.0 - eff_lambda) * 0.04;
    let lambda_low = (eff_lambda - lambda_uncertainty).max(0.0);
    let lambda_high = (eff_lambda + lambda_uncertainty).min(1.0);

    QueryResult {
        sample_file: String::new(),
        genome_file: String::new(),
        contig_name: String::new(),
        adjusted_ani,
        eff_cov,
        ani_percentile: (ani_low, ani_high),
        eff_lambda,
        lambda_percentile: (lambda_low, lambda_high),
        median_cov: 1.0,
        mean_cov_geq1: 1.0,
        containment_ind: format!("{}/{}", shared_tags, total_ref_tags),
        naive_ani,
        ref_tags: total_ref_tags,
        shared_tags,
        query_tags,
        taxonomic_abundance: 0.0,
        sequence_abundance: 0.0,
    }
}

fn filter_results(result: &QueryResult, min_ani: Option<f64>) -> bool {
    // 只有在有共享标签时才进行过滤
    if result.shared_tags == 0 {
        return false;
    }

    // 当计算丰度时，只要有共享标签就包含在结果中
    if result.shared_tags > 0 {
        return true;
    }

    // 基本过滤条件
    if result.eff_cov < MIN_COVERAGE {
        return false;
    }

    // ANI 过滤
    if let Some(min_ani) = min_ani {
        if result.adjusted_ani < min_ani {
            return false;
        }
    } else if result.adjusted_ani < MIN_ANI {
        return false;
    }

    true
}

// 新增profile专用的过滤函数
fn filter_results_for_profile(result: &QueryResult, min_ani: Option<f64>) -> bool {
    // 只有在有共享标签时才进行过滤
    if result.shared_tags == 0 {
        return false;
    }

    // profile模式下的最小共享标签数过滤（更严格）
    if result.shared_tags < MIN_SHARED_TAGS {
        return false;
    }

    // profile模式下的最小覆盖率过滤（更严格）
    if result.eff_cov < PROFILE_MIN_COVERAGE {
        return false;
    }

    // profile模式下的ANI过滤（更严格）
    let effective_min_ani = min_ani.unwrap_or(PROFILE_MIN_ANI);
    if result.adjusted_ani < effective_min_ani {
        return false;
    }

    // 最小标签数过滤（确保genome有足够的标签）
    if result.ref_tags < MIN_TAGS_FOR_GENOME {
        return false;
    }

    true
}

// 内部函数：使用缓存的数据库数据进行查询 - 优化大文件读取
fn query_single_file_with_cached_db(
    sample_path: &str, 
    db_path: &str, 
    cached_db_entries: &[SyldbEntry], 
    cached_sample_entries: &FxHashMap<String, Vec<SylspEntry>>,
    min_ani: f64
) -> Result<Vec<QueryResult>> {
    eprintln!("Processing sample file with cached database: {}", sample_path);
    
    // 从缓存中获取样本数据
    let sample_entries = cached_sample_entries.get(sample_path)
        .ok_or_else(|| anyhow!("Sample file not found in cache: {}", sample_path))?;

    eprintln!("Found {} entries in sample", sample_entries.len());

    // 检查样本数据的有效性
    if sample_entries.is_empty() {
        eprintln!("Warning: Sample {} has no tags", sample_path);
        return Ok(Vec::new());
    }

    // 按样本源分组 - 这是关键：处理合并文件中的多个样本
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in sample_entries {
        sample_groups.entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }

    eprintln!("Found {} different sample sources in file: {:?}", 
              sample_groups.len(), 
              sample_groups.keys().collect::<Vec<_>>());

    // 并行处理每个样本组，然后合并结果
    let mut all_results: Vec<QueryResult> = sample_groups.par_iter()
        .flat_map(|(sample_source, entries)| {
            eprintln!("Processing sample source: {} with {} entries", sample_source, entries.len());
            
            // 构建样本标签的哈希表 - 使用更高效的HashSet
            let sample_tags: HashSet<Hash> = entries.iter()
                .map(|entry| entry.tag.clone())
                .collect();

            let total_sample_tags = entries.len();

            // 并行处理每个基因组记录进行比对
            cached_db_entries.par_iter().filter_map(|db_entry| {
                // 最小标签数过滤（参考sylph的min_number_kmers）
                if db_entry.tags.len() < MIN_TAGS_FOR_GENOME {
                    return None;
                }

                // 计算共享标签和统计信息 - 优化计算方式
                let shared_tags = db_entry.tags.iter()
                    .filter(|tag| sample_tags.contains(tag))
                    .count();

                let total_ref_tags = db_entry.tags.len();

                // 计算统计数据
                let mut result = calculate_statistics(
                    shared_tags,
                    total_sample_tags,
                    total_ref_tags,
                );

                // 设置基本信息 - 关键：使用实际的样本源ID
                result.sample_file = sample_source.clone();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();
                result.shared_tags = shared_tags;
                result.query_tags = total_sample_tags;
                result.ref_tags = total_ref_tags;

                // 计算平均深度和覆盖度
                if shared_tags > 0 {
                    result.mean_cov_geq1 = 1.0;
                    result.eff_cov = shared_tags as f64 / total_ref_tags as f64;
                    result.median_cov = 1.0;
                }

                // 应用profile专用的过滤条件
                if filter_results_for_profile(&result, Some(min_ani)) {
                    Some(result)
                } else {
                    None
                }
            }).collect::<Vec<QueryResult>>()
        })
        .collect();

    // 按ANI排序（参考sylph的排序机制）
    all_results.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());

    eprintln!("Generated {} results for file {}", all_results.len(), sample_path);
    Ok(all_results)
}

// 添加新的公共函数用于单个文件的查询（保持原有接口不变）
pub fn query_single_file(sample_path: &str, db_path: &str, min_ani: f64) -> Result<Vec<QueryResult>> {
    eprintln!("Processing database file: {}", db_path);
    
    // 读取数据库文件
    let db_file = File::open(db_path)
        .with_context(|| format!("Failed to open database file: {}", db_path))?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;

    eprintln!("Found {} entries in database", db_entries.len());

    // 读取样本文件 - 优化大文件读取
    let sample_file = File::open(sample_path)
        .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
    let sample_reader = BufReader::with_capacity(100_000_000, sample_file); // 100MB 缓冲区
    let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
        .with_context(|| format!("Failed to deserialize sample file: {}", sample_path))?;

    eprintln!("Found {} entries in sample", sample_entries.len());

    // 检查样本数据的有效性
    if sample_entries.is_empty() {
        eprintln!("Warning: Sample {} has no tags", sample_path);
        return Ok(Vec::new());
    }

    // 按样本源分组
    let mut sample_groups: FxHashMap<String, Vec<&SylspEntry>> = FxHashMap::default();
    for entry in &sample_entries {
        sample_groups.entry(entry.sample_source.clone())
            .or_default()
            .push(entry);
    }

    // 并行处理每个样本组，然后合并结果
    let all_results: Vec<QueryResult> = sample_groups.par_iter()
        .flat_map(|(sample_source, entries)| {
            // 构建样本标签的哈希表 - 使用更高效的HashSet
            let sample_tags: HashSet<Hash> = entries.iter()
                .map(|entry| entry.tag.clone())
                .collect();

            let total_sample_tags = entries.len();
            eprintln!("Total unique tags in sample {}: {}", sample_source, total_sample_tags);

            // 并行处理每个基因组记录进行比对
            db_entries.par_iter().filter_map(|db_entry| {
                // 计算共享标签和统计信息
                let shared_tags = db_entry.tags.par_iter()
                    .map(|tag| if sample_tags.contains(tag) { 1 } else { 0 })
                    .sum::<usize>();

                let total_ref_tags = db_entry.tags.len();

                eprintln!("Found {} shared tags between sample {} and reference {}", 
                         shared_tags, sample_source, db_entry.sequence_id);

                // 计算统计数据
                let mut result = calculate_statistics(
                    shared_tags,
                    total_sample_tags,
                    total_ref_tags,
                );

                // 设置基本信息
                result.sample_file = sample_source.clone();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();
                result.shared_tags = shared_tags;
                result.query_tags = total_sample_tags;
                result.ref_tags = total_ref_tags;

                // 计算平均深度和覆盖度
                if shared_tags > 0 {
                    result.mean_cov_geq1 = 1.0;
                    result.eff_cov = shared_tags as f64 / total_ref_tags as f64;
                    
                    // 简化的中位数覆盖度计算（避免生成大量向量）
                    result.median_cov = 1.0;
                }

                // 应用过滤条件
                if filter_results(&result, Some(min_ani)) {
                    eprintln!("Result passed filters: ANI={:.2}, Coverage={:.3}", 
                            result.adjusted_ani, result.eff_cov);
                    Some(result)
                } else {
                    eprintln!("Result filtered out: ANI={:.2}, Coverage={:.3}", 
                            result.adjusted_ani, result.eff_cov);
                    None
                }
            }).collect::<Vec<QueryResult>>()
        })
        .collect();

    Ok(all_results)
}

// 从文件路径中提取基因组ID的函数
fn extract_genome_id(file_path: &str) -> String {
    let path = std::path::Path::new(file_path);
    
    // 获取文件名（不含路径）
    let file_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    // 移除.fasta.gz或.fasta扩展名
    let name_without_ext = file_name
        .strip_suffix(".fasta.gz")
        .or_else(|| file_name.strip_suffix(".fasta"))
        .unwrap_or(file_name);
    
    // 从test_files/目录中提取
    let clean_name = if let Some(stripped) = name_without_ext.strip_prefix("test_files/") {
        stripped
    } else {
        name_without_ext
    };
    
    clean_name.to_string()
}

// 读取taxonomy文件并建立genome到分类信息的映射
fn read_taxonomy_file(taxonomy_file: &str) -> Result<FxHashMap<String, Arc<TaxonomyInfo>>> {
    use std::io::BufRead;
    
    let file = File::open(taxonomy_file)
        .with_context(|| format!("Failed to open taxonomy file: {}", taxonomy_file))?;
    let reader = BufReader::new(file);
    
    // 预分配 HashMap 容量以提高性能
    let mut taxonomy_map = FxHashMap::default(); // 使用默认容量
    let mut line_count = 0;
    
    for line in reader.lines() {
        line_count += 1;
        let line = line.context("Failed to read line from taxonomy file")?;
        let line = line.trim();
        
        // 跳过标题行和空行
        if line_count == 1 || line.is_empty() || line.starts_with('#') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            eprintln!("Warning: Invalid line format at line {}: {}", line_count, line);
            continue;
        }
        
        let accession = parts[0].trim();
        let gtdb_taxonomy = parts[1].trim();
        
        // 处理genome ID：去掉前缀 (RS_, GB_, 等)
        let genome_id = if accession.starts_with("RS_") {
            &accession[3..]
        } else if accession.starts_with("GB_") {
            &accession[3..]
        } else if accession.starts_with("GS_") {
            &accession[3..]
        } else {
            accession
        };
        
        // 解析GTDB分类信息并包装在 Arc 中
        match TaxonomyInfo::from_gtdb_string(gtdb_taxonomy) {
            Ok(taxonomy) => {
                let taxonomy_arc = Arc::new(taxonomy);
                // 添加基本ID（例如：GCF_000006685.1）
                taxonomy_map.insert(genome_id.to_string(), Arc::clone(&taxonomy_arc));
                
                // 添加带_genomic后缀的ID（例如：GCF_000006685.1_genomic）
                let genomic_id = format!("{}_genomic", genome_id);
                taxonomy_map.insert(genomic_id, taxonomy_arc);
            }
            Err(e) => {
                eprintln!("Warning: Failed to parse taxonomy for {}: {}", genome_id, e);
            }
        }
    }
    
    eprintln!("Loaded taxonomy information for {} genome variants", taxonomy_map.len());
    Ok(taxonomy_map)
}

// 从基因组级别聚合到物种级别
fn aggregate_to_species_level(
    sample_groups: &HashMap<String, Vec<GenomeProfileResult>>,
    taxonomy_map: &FxHashMap<String, Arc<TaxonomyInfo>>,
    effective_min_ani: f64,
) -> Result<Vec<SpeciesAbundanceResult>> {
    use std::sync::Mutex;
    
    let species_map = Arc::new(Mutex::new(FxHashMap::<String, SpeciesAbundanceResult>::default()));
    
    // 采用 sylph 的分层并行策略进行物种聚合
    let sample_groups_arc = Arc::new(sample_groups);
    
    // 外层并行：处理样本组
    sample_groups_arc.par_iter()
        .for_each(|(sample_id, genome_results)| {
            // 为每个样本组创建局部聚合结果 - 预分配容量
            let mut local_species_map: FxHashMap<String, SpeciesAbundanceResult> = FxHashMap::default();
            
            for genome_result in genome_results {
                // 额外的过滤条件：确保只有高质量的genome参与物种聚合
                if genome_result.adjusted_ani < effective_min_ani || 
                   genome_result.eff_cov < PROFILE_MIN_COVERAGE ||
                   genome_result.common_tags < MIN_SHARED_TAGS {
                    continue;
                }
                
                // 从genome_id中提取标准化的genome标识符
                let genome_id = extract_genome_id_from_path(&genome_result.genome_id);
                
                // 查找对应的分类信息 - 使用字符串切片
                if let Some(taxonomy_arc) = taxonomy_map.get(genome_id) {
                    let species_key = taxonomy_arc.get_species_key();
                    
                    // 获取或创建物种条目 - 使用 Arc 共享而不是克隆
                    let species_result = local_species_map.entry(species_key).or_insert_with(|| {
                        SpeciesAbundanceResult {
                            taxonomy: Arc::clone(taxonomy_arc),
                            sample_abundances: FxHashMap::default(),
                            total_tags: 0,
                            genome_count: 0,
                            reads_count: 0,
                            gscore: 0.0,
                        }
                    });
                    
                    // 累加样本丰度
                    *species_result.sample_abundances.entry(sample_id.clone()).or_insert(0.0) += 
                        genome_result.taxonomic_abundance;
                    
                    // 累加标签数、基因组计数和 reads 数
                    species_result.total_tags += genome_result.total_tags;
                    species_result.genome_count += 1;
                    // 使用 common_tags 作为该基因组在该样本中的 reads 数代理
                    species_result.reads_count += genome_result.common_tags;
                } else {
                    eprintln!("Warning: No taxonomy information found for genome: {}", genome_id);
                }
            }
            
            // 将局部结果合并到全局结果中
            let mut global_map = species_map.lock().unwrap();
            for (species_key, local_result) in local_species_map {
                let global_result = global_map.entry(species_key).or_insert_with(|| {
                    SpeciesAbundanceResult {
                        taxonomy: Arc::clone(&local_result.taxonomy),
                        sample_abundances: FxHashMap::default(),
                        total_tags: 0,
                        genome_count: 0,
                        reads_count: 0,
                        gscore: 0.0,
                    }
                });
                
                // 合并样本丰度
                for (sample_id, abundance) in local_result.sample_abundances {
                    *global_result.sample_abundances.entry(sample_id).or_insert(0.0) += abundance;
                }
                
                // 合并标签数、基因组计数和 reads 数
                global_result.total_tags += local_result.total_tags;
                global_result.genome_count += local_result.genome_count;
                global_result.reads_count += local_result.reads_count;
            }
        });
    
    let species_map = Arc::try_unwrap(species_map).unwrap().into_inner().unwrap();
    let mut results: Vec<SpeciesAbundanceResult> = species_map.into_values().collect();
    
    // 按分类学层级排序
    results.sort_by(|a, b| {
        a.taxonomy.kingdom.cmp(&b.taxonomy.kingdom)
            .then_with(|| a.taxonomy.phylum.cmp(&b.taxonomy.phylum))
            .then_with(|| a.taxonomy.class.cmp(&b.taxonomy.class))
            .then_with(|| a.taxonomy.order.cmp(&b.taxonomy.order))
            .then_with(|| a.taxonomy.family.cmp(&b.taxonomy.family))
            .then_with(|| a.taxonomy.genus.cmp(&b.taxonomy.genus))
            .then_with(|| a.taxonomy.species.cmp(&b.taxonomy.species))
    });
    
    Ok(results)
}

// 从文件路径或genome_id中提取标准化的genome标识符
fn extract_genome_id_from_path(input: &str) -> &str {
    // 如果输入包含路径分隔符，提取文件名
    let file_name = if input.contains('/') {
        std::path::Path::new(input)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(input)
    } else {
        input
    };
    
    // 移除可能的文件扩展名 - 返回字符串切片而不是 String
    file_name
        .strip_suffix(".fasta.gz").or_else(|| file_name.strip_suffix(".fasta"))
        .or_else(|| file_name.strip_suffix(".fa.gz")).or_else(|| file_name.strip_suffix(".fa"))
        .or_else(|| file_name.strip_suffix(".fna.gz")).or_else(|| file_name.strip_suffix(".fna"))
        .unwrap_or(file_name)
}

// 从syldb文件中读取基因组映射关系
fn read_genome_mapping(db_path: &str) -> Result<FxHashMap<String, (String, String)>> {
    let db_file = File::open(db_path)?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;
    
    // 并行处理数据库条目生成映射
    let genome_map: FxHashMap<String, (String, String)> = db_entries.par_iter()
        .map(|entry| {
            // 获取原始基因组文件路径
            let genome_source = entry.genome_source.clone();
            let genome_id = if let Some(file_name) = std::path::Path::new(&genome_source)
                .file_name()
                .and_then(|s| s.to_str()) 
            {
                // 移除.fasta.gz或.fasta扩展名
                file_name.strip_suffix(".fasta.gz")
                    .or_else(|| file_name.strip_suffix(".fasta"))
                    .unwrap_or(file_name)
                    .to_string()
            } else {
                genome_source.clone()
            };
            
            // 返回(序列ID, (基因组ID, 基因组源文件))
            (entry.sequence_id.clone(), (genome_id, genome_source))
        })
        .collect();
    
    Ok(genome_map)
}

// 读取样本文件列表
fn read_sample_list(list_file: &str) -> Result<Vec<String>> {
    let content = std::fs::read_to_string(list_file)?;
    // 预分配 Vec 容量并减少字符串克隆
    let lines: Vec<String> = content.lines()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect();
    Ok(lines)
}

// 生成TSV格式的丰度矩阵
fn write_abundance_matrix(
    sample_groups: &HashMap<String, Vec<GenomeProfileResult>>,
    all_genomes: &HashSet<String>,
    log_path: Option<String>,
    tsv_name: &str,
    writer: &mut Box<dyn Write + Send>,
) -> Result<()> {
    // 如果指定了log_path，使用它，否则使用当前目录
    let output_dir = if let Some(path) = log_path {
        PathBuf::from(path)
    } else {
        std::env::current_dir()?
    };

    // 确保目录存在
    std::fs::create_dir_all(&output_dir)?;

    // 构建TSV文件路径
    let tsv_path = output_dir.join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(tsv_path)?);

    // 获取所有样本ID并排序
    let mut sample_ids: Vec<_> = sample_groups.keys().collect();
    sample_ids.sort();

    // 写入表头
    write!(tsv_writer, "Genome")?;
    write!(writer, "\nAbundance Matrix:\n")?;
    write!(writer, "Genome")?;
    for sample_id in &sample_ids {
        write!(tsv_writer, "\t{}", sample_id)?;
        write!(writer, "\t{}", sample_id)?;
    }
    writeln!(tsv_writer)?;
    writeln!(writer)?;

    // 采用 sylph 的高效并行数据收集策略
    let genome_data: Vec<(String, Vec<f64>)> = all_genomes.par_iter()
        .map(|genome_id| {
            let abundances: Vec<f64> = sample_ids.iter()
                .map(|sample_id| {
                    sample_groups.get(sample_id.as_str())
                        .and_then(|results| results.iter()
                            .find(|r| r.genome_id == *genome_id))
                        .map(|r| r.taxonomic_abundance)
                        .unwrap_or(0.0)
                })
                .collect();
            (genome_id.clone(), abundances)
        })
        .collect();

    // 写入每个基因组的丰度数据
    for (genome_id, abundances) in genome_data {
        write!(tsv_writer, "{}", genome_id)?;
        write!(writer, "{}", genome_id)?;
        for abundance in abundances {
            write!(tsv_writer, "\t{:.4}", abundance)?;
            write!(writer, "\t{:.4}", abundance)?;
        }
        writeln!(tsv_writer)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;

    Ok(())
}

// 生成物种级别的TSV格式丰度矩阵
fn write_species_abundance_matrix(
    species_results: &[SpeciesAbundanceResult],
    all_samples: &HashSet<String>,
    log_path: Option<String>,
    tsv_name: &str,
    writer: &mut Box<dyn Write + Send>,
) -> Result<()> {
    // 如果指定了log_path，使用它，否则使用当前目录
    let output_dir = if let Some(path) = log_path {
        PathBuf::from(path)
    } else {
        std::env::current_dir()?
    };

    // 确保目录存在
    std::fs::create_dir_all(&output_dir)?;

    // 构建TSV文件路径
    let tsv_path = output_dir.join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(tsv_path)?);

    // 获取所有样本ID并排序
    let mut sample_ids: Vec<_> = all_samples.iter().collect();
    sample_ids.sort();

    // 写入表头 (参考Abundance_Stat.all.xls格式)
    write!(tsv_writer, "#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")?;
    write!(writer, "\nSpecies-level Abundance Matrix:\n")?;
    write!(writer, "#Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies")?;
    for sample_id in &sample_ids {
        write!(tsv_writer, "\t{}", sample_id)?;
        write!(writer, "\t{}", sample_id)?;
    }
    writeln!(tsv_writer)?;
    writeln!(writer)?;

    // 采用 sylph 的高效并行数据收集策略
    let species_data: Vec<(Arc<TaxonomyInfo>, Vec<f64>)> = species_results.par_iter()
        .map(|species_result| {
            let abundances: Vec<f64> = sample_ids.iter()
                .map(|sample_id| {
                    species_result.sample_abundances
                        .get(sample_id.as_str())
                        .copied()
                        .unwrap_or(0.0)
                })
                .collect();
            (Arc::clone(&species_result.taxonomy), abundances)
        })
        .collect();

    // 写入每个物种的丰度数据
    for (taxonomy_arc, abundances) in species_data {
        // 写入分类学信息（7列）
        write!(tsv_writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
               taxonomy_arc.kingdom, taxonomy_arc.phylum, taxonomy_arc.class,
               taxonomy_arc.order, taxonomy_arc.family, taxonomy_arc.genus, taxonomy_arc.species)?;
        write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
               taxonomy_arc.kingdom, taxonomy_arc.phylum, taxonomy_arc.class,
               taxonomy_arc.order, taxonomy_arc.family, taxonomy_arc.genus, taxonomy_arc.species)?;
        
        // 写入各个样本的丰度值
        for abundance in abundances {
            write!(tsv_writer, "\t{:.6}", abundance)?;
            write!(writer, "\t{:.6}", abundance)?;
        }
        writeln!(tsv_writer)?;
        writeln!(writer)?;
    }
    writeln!(writer)?;

    Ok(())
}

// 从缓存的数据库条目中构建基因组映射关系
fn build_genome_mapping_from_cache(cached_db_entries: &[SyldbEntry]) -> FxHashMap<String, (String, String)> {
    // 预分配 HashMap 容量以提高性能
    let mut genome_map = FxHashMap::default();
    
    // 并行处理数据库条目生成映射
    for entry in cached_db_entries {
        // 获取原始基因组文件路径
        let genome_source = entry.genome_source.clone();
        let genome_id = if let Some(file_name) = std::path::Path::new(&genome_source)
            .file_name()
            .and_then(|s| s.to_str()) 
        {
            // 移除.fasta.gz或.fasta扩展名
            file_name.strip_suffix(".fasta.gz")
                .or_else(|| file_name.strip_suffix(".fasta"))
                .unwrap_or(file_name)
                .to_string()
        } else {
            genome_source.clone()
        };
        
        // 返回(序列ID, (基因组ID, 基因组源文件))
        genome_map.insert(entry.sequence_id.clone(), (genome_id, genome_source));
    }
    
    genome_map
}

// 更新profile函数
pub fn profile(args: ProfileArgs) -> Result<()> {
    // 处理minimum_ani参数：如果没有传入参数，使用默认值
    let effective_min_ani = args.minimum_ani.unwrap_or(PROFILE_MIN_ANI);
    eprintln!("Using minimum ANI threshold: {:.1}%", effective_min_ani);
    
    // 优化线程池配置 - 采用 sylph 的策略
    let _max_ram = args.threads * 2; // 简单的内存限制，每线程2GB
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Failed to initialize thread pool")?;

    // 一次性读取并缓存数据库文件 - 优化大文件读取
    eprintln!("Loading database file: {}", args.db_file);
    
    let db_file = File::open(&args.db_file)
        .with_context(|| format!("Failed to open database file: {}", args.db_file))?;
    let db_reader = BufReader::with_capacity(100_000_000, db_file); // 100MB 缓冲区
    let cached_db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", args.db_file))?;
    
    eprintln!("Cached {} entries from database", cached_db_entries.len());

    // 一次性读取并缓存所有样本文件 - 优化大文件读取
    eprintln!("Loading sample files: {}", args.sample_file);
    let sample_files: Vec<String> = if args.sample_file.ends_with(".txt") {
        read_sample_list(&args.sample_file)?
    } else {
        vec![args.sample_file.clone()]
    };

    let mut cached_sample_entries: FxHashMap<String, Vec<SylspEntry>> = FxHashMap::default();
    for sample_path in &sample_files {
        let sample_file = File::open(sample_path)
            .with_context(|| format!("Failed to open sample file: {}", sample_path))?;
        let sample_reader = BufReader::with_capacity(100_000_000, sample_file); // 100MB 缓冲区
        let sample_entries: Vec<SylspEntry> = bincode::deserialize_from(sample_reader)
            .with_context(|| format!("Failed to deserialize sample file: {}", sample_path))?;
        cached_sample_entries.insert(sample_path.clone(), sample_entries);
    }
    eprintln!("Cached {} sample files", cached_sample_entries.len());

    // 从缓存的数据库构建基因组映射关系
    let genome_mapping = build_genome_mapping_from_cache(&cached_db_entries);
    
    // 创建输出写入器
    let mut writer = create_multi_writer(&args.out_file_name)?;



    // 存储所有样本的结果 - 预分配容量，使用 Mutex 保护
    let all_results = Arc::new(Mutex::new(FxHashMap::<(String, String), GenomeProfileResult>::default()));

        // 采用 sylph 的简化并行处理策略
    let step = usize::max(args.threads/3 + 1, usize::min(sample_files.len(), args.threads));
    let chunks: Vec<Vec<String>> = sample_files.chunks(step).map(|chunk| chunk.to_vec()).collect();
    
    // 使用 sylph 风格的分块处理，集成k-mer重新分配机制
    chunks.into_iter().for_each(|chunk| {
        chunk.into_par_iter().for_each(|sample_file| {
            // 第一阶段：计算初步结果（不使用重新分配）
            if let Ok(initial_results) = query_single_file_with_cached_db(&sample_file, &args.db_file, &cached_db_entries, &cached_sample_entries, effective_min_ani) {
                // 按ANI排序
                let mut initial_results = initial_results;
                initial_results.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());
                
                // 第二阶段：构建winner table并重新分配（模仿sylph的两阶段处理）
                eprintln!("{} taxonomic profiling; reassigning tags for {} genomes...", &sample_file, initial_results.len());
                
                // 构建winner table
                let winner_map = build_winner_table(&initial_results, &cached_db_entries, true); // 启用日志
                
                // 使用winner table重新计算结果
                if let Some(sample_entries) = cached_sample_entries.get(&sample_file) {
                    let mut reassigned_results = recalculate_with_winner_table(
                        &cached_db_entries,
                        sample_entries,
                        &winner_map,
                        effective_min_ani,
                        false
                    );
                    
                    // 第三阶段：过滤过度重新分配的基因组
                    reassigned_results = filter_over_reassigned_genomes(
                        &initial_results,
                        &reassigned_results,
                        effective_min_ani,
                        K
                    );
                    
                    // 第四阶段：重新计算丰度
                    recalculate_abundances_after_reassignment(&mut reassigned_results, sample_entries);
                    
                    eprintln!("{} has {} genomes passing profiling threshold after reassignment.", &sample_file, reassigned_results.len());
                    
                    // 按基因组ID分组结果 - 修复：确保每个样本源都被正确处理
                    for result in reassigned_results {
                        if let Some((genome_id, _)) = genome_mapping.get(&result.contig_name) {
                            // 关键修复：使用实际的样本源ID作为key的一部分
                            let key = (genome_id.clone(), result.sample_file.clone());
                            let mut all_results = all_results.lock().unwrap();
                            let entry = all_results.entry(key)
                                .or_insert_with(|| {
                                    GenomeProfileResult {
                                        genome_id: genome_id.clone(),
                                        sample_id: result.sample_file.clone(), // 这里保存的是实际的样本源ID
                                        file_path: sample_file.clone(),
                                        adjusted_ani: 0.0,
                                        taxonomic_abundance: 0.0,
                                        sequence_abundance: 0.0,
                                        common_tags: 0,
                                        total_tags: 0,
                                        eff_cov: 0.0,
                                    }
                                });
                            
                            // 累加标签数
                            entry.common_tags += result.shared_tags;
                            entry.total_tags += result.ref_tags;
                            entry.eff_cov += result.eff_cov;
                            
                            // 使用共享标签数作为权重计算加权平均ANI
                            if entry.common_tags > 0 {
                                entry.adjusted_ani = (entry.adjusted_ani * (entry.common_tags - result.shared_tags) as f64 
                                    + result.adjusted_ani * result.shared_tags as f64) / entry.common_tags as f64;
                            }
                        }
                    }
                }
            }
        });
    });
    
    // 收集所有基因组ID
    let mut all_genomes: HashSet<String> = HashSet::new();
    for entry in genome_mapping.values() {
        all_genomes.insert(entry.0.clone());
    }

    // 转换为向量以便排序和分组
    let results: Vec<_> = all_results.lock().unwrap().values().cloned().collect();
    
    // 按样本分组计算丰度
    let mut sample_groups: HashMap<String, Vec<GenomeProfileResult>> = HashMap::new();
    for result in results {
        // 使用实际的样本来源而不是文件名
        sample_groups.entry(result.sample_id.clone())
            .or_default()
            .push(result);
    }
    
    // 采用 sylph 的简单策略 - 顺序计算丰度，避免复杂的并行迭代器组合
    for (_sample_id, group) in sample_groups.iter_mut() {
        // 按ANI排序（参考sylph的排序机制）
        group.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());
        
        // 过滤掉不符合profile要求的genome
        group.retain(|r| {
            r.common_tags >= MIN_SHARED_TAGS && 
            r.eff_cov >= PROFILE_MIN_COVERAGE && 
            r.adjusted_ani >= effective_min_ani &&
            r.total_tags >= MIN_TAGS_FOR_GENOME
        });
        
        // 计算总覆盖度，包括所有检测到的标签
        let total_genome_cov: f64 = group.iter()
            .map(|r| if r.common_tags > 0 { r.eff_cov } else { 0.0 })
            .sum();
        
        let total_seq_cov: f64 = group.iter()
            .map(|r| if r.common_tags > 0 { 
                r.eff_cov * r.total_tags as f64 
            } else { 
                0.0 
            })
            .sum();
        
        // 计算每个结果的丰度 - 采用 sylph 的顺序处理方式
        for result in group.iter_mut() {
            // 只要有共享标签就计算丰度
            if result.common_tags > 0 {
                result.taxonomic_abundance = if total_genome_cov > 0.0 {
                    result.eff_cov / total_genome_cov * 100.0
                } else {
                    0.0
                };
                
                result.sequence_abundance = if total_seq_cov > 0.0 {
                    result.eff_cov * result.total_tags as f64 / total_seq_cov * 100.0
                } else {
                    0.0
                };
            } else {
                result.taxonomic_abundance = 0.0;
                result.sequence_abundance = 0.0;
            }
        }
    }

    // 检查是否提供了taxonomy文件以进行物种级别聚合
    if let Some(taxonomy_file) = &args.taxonomy_file {
        eprintln!("Loading taxonomy information from: {}", taxonomy_file);
        
        // 读取分类学信息
        let taxonomy_map = read_taxonomy_file(taxonomy_file)?;
        
        // 聚合到物种级别
        let mut species_results = aggregate_to_species_level(&sample_groups, &taxonomy_map, effective_min_ani)?;
        
        // 获取所有样本ID
        let all_samples: HashSet<String> = sample_groups.keys().cloned().collect();
        
        // 生成过滤前的物种级别TSV格式丰度矩阵
        let pre_filter_tsv_name = format!("pre_gscore_filter_{}", args.tsv_name);
        eprintln!("Writing pre-filter species abundance matrix: {}", pre_filter_tsv_name);
        write_species_abundance_matrix(&species_results, &all_samples, args.log_path.clone(), &pre_filter_tsv_name, &mut writer)?;
        
        // 应用 G-score 过滤
        eprintln!("Applying G-score filtering with threshold: {:.2}", args.gscore_threshold);
        species_results = filter_species_by_gscore(&mut species_results, args.gscore_threshold);
        
        // 生成过滤后的物种级别TSV格式丰度矩阵
        eprintln!("Writing post-filter species abundance matrix: {}", args.tsv_name);
        write_species_abundance_matrix(&species_results, &all_samples, args.log_path.clone(), &args.tsv_name, &mut writer)?;
        
        // 输出物种级别的统计信息
        writeln!(writer, "Species-level Profile Results:")?;
        writeln!(writer, "------------------------------")?;
        writeln!(writer, "Sample files: {} files processed", sample_files.len())?;
        writeln!(writer, "Database file: {}", args.db_file)?;
        writeln!(writer, "Taxonomy file: {}", taxonomy_file)?;
        writeln!(writer, "Total species detected: {}", species_results.len())?;
        writeln!(writer, "\nSpecies composition summary:")?;
        writeln!(writer, "{:<50} {:<15} {:<15} {:<15} {:<10} {:<10}", 
            "Species", "Genomes", "Total_Tags", "Reads_Count", "G-score", "Avg_Abundance")?;
        writeln!(writer, "{:-<120}", "")?;
        
        for species_result in &species_results {
            let avg_abundance: f64 = species_result.sample_abundances.values().sum::<f64>() 
                / species_result.sample_abundances.len() as f64;
            let species_name = if species_result.taxonomy.species.is_empty() {
                format!("{}_sp", species_result.taxonomy.genus)
            } else {
                species_result.taxonomy.species.clone()
            };
            
            writeln!(writer, "{:<50} {:<15} {:<15} {:<15} {:<10.2} {:<10.4}", 
                species_name,
                species_result.genome_count,
                species_result.total_tags,
                species_result.reads_count,
                species_result.gscore,
                avg_abundance)?;
        }
        
    } else {
        // 原始的基因组级别输出
        // 生成TSV格式的丰度矩阵
        write_abundance_matrix(&sample_groups, &all_genomes, args.log_path.clone(), &args.tsv_name, &mut writer)?;

        // 将所有结果收集到一个新的向量中
        let mut final_results: Vec<GenomeProfileResult> = sample_groups.into_values().flatten().collect();
        
        // 按基因组ID和ANI值排序
        final_results.sort_by(|a, b| {
            a.genome_id.cmp(&b.genome_id)
                .then_with(|| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap())
        });
        
        // 输出结果
        writeln!(writer, "Genome-level Profile Results:")?;
        writeln!(writer, "-----------------------------")?;
        writeln!(writer, "Sample files: {} files processed", sample_files.len())?;
        writeln!(writer, "Database file: {}", args.db_file)?;
        writeln!(writer, "\nGenome composition:")?;
        writeln!(writer, "{:<30} {:<20} {:<10} {:<12} {:<12} {:<12} {:<12} {:<10}", 
            "Genome_ID", "Sample_ID", "ANI(%)", "Tax_Abund(%)", "Seq_Abund(%)", "Common_Tags", "Total_Tags", "Eff_cov")?;
        writeln!(writer, "{:-<110}", "")?;
        
        let mut current_genome = String::new();
        for result in final_results {
            if current_genome != result.genome_id {
                if !current_genome.is_empty() {
                    writeln!(writer)?;
                }
                current_genome = result.genome_id.clone();
            }
            
            writeln!(writer, "{:<30} {:<20} {:<10.2} {:<12.2} {:<12.2} {:<12} {:<12} {:<10.3}", 
                result.genome_id,
                result.sample_id,  // 使用实际的样本来源
                result.adjusted_ani,
                result.taxonomic_abundance,
                result.sequence_abundance,
                result.common_tags,
                result.total_tags,
                result.eff_cov)?;
        }
    }
    
    Ok(())
}