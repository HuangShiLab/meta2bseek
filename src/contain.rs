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
use crate::constants::{Hash, hash_bytes};
use std::thread;
use std::time::Duration;
use memory_stats::memory_stats;

pub use crate::extract::{SyldbEntry, SylspEntry, TagHash, one_mismatch_canonical_hashes, ani_from_containment_one_mismatch};


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
// genome 用 u32 索引（通过 GenomeInterner 驻留），而不是每个 tag 都存一份 String，
// winner map 可能有上千万条目，String→u32 大幅降低内存并把比较变成整数比较。
#[derive(Debug, Clone)]
struct WinnerTableEntry {
    pub ani: f64,
    pub genome_idx: u32,
    pub was_reassigned: bool,
}

/// 把 genome_id 字符串驻留成 u32 索引，供 winner table 复用。
struct GenomeInterner {
    ids: Vec<String>,
    idx: FxHashMap<String, u32>,
}

impl GenomeInterner {
    fn new() -> Self {
        Self { ids: Vec::new(), idx: FxHashMap::default() }
    }
    fn intern(&mut self, id: &str) -> u32 {
        if let Some(&i) = self.idx.get(id) {
            return i;
        }
        let i = self.ids.len() as u32;
        self.ids.push(id.to_string());
        self.idx.insert(id.to_string(), i);
        i
    }
    fn get(&self, id: &str) -> Option<u32> {
        self.idx.get(id).copied()
    }
    fn resolve(&self, i: u32) -> &str {
        &self.ids[i as usize]
    }
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
) -> (FxHashMap<Hash, WinnerTableEntry>, GenomeInterner) {
    eprintln!("Building winner table for {} results and {} database entries", results.len(), db_entries.len());

    let mut tag_to_genome_map: FxHashMap<Hash, WinnerTableEntry> = FxHashMap::default();
    let mut interner = GenomeInterner::new();

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
            // 每个基因组只驻留一次，tag 循环里只用 u32，避免逐 tag 的 String 分配
            let genome_idx = interner.intern(extract_genome_id_from_path(&db_entry.genome_source));

            // 借鉴sylph的简洁遍历方式
            for tag in &db_entry.tags {
                let entry = tag_to_genome_map.entry(*tag).or_insert_with(|| {
                    WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_idx,
                        was_reassigned: false,
                    }
                });

                // 关键：选择ANI最高的基因组作为该tag的"赢家"
                if res.adjusted_ani > entry.ani {
                    *entry = WinnerTableEntry {
                        ani: res.adjusted_ani,
                        genome_idx,
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
                        let winner_genome = interner.resolve(winner_entry.genome_idx);
                        if winner_genome != res.genome_file {
                            if let Some(&winner_index) = genome_to_index.get(winner_genome) {
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

    // tag 重叠分析纯属诊断信息，却要对所有基因组的所有 tag 再建一次 FxHashMap，
    // 代价高。仅在 debug 日志开启时执行。
    if log::log_enabled!(log::Level::Debug) {
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
    }

    (tag_to_genome_map, interner)
}

// 使用winner table重新计算统计结果 - 模仿sylph的get_stats函数
fn recalculate_with_winner_table(
    db_entries: &[SyldbEntry],
    sample_entries: &[SylspEntry],
    winner_map: &FxHashMap<Hash, WinnerTableEntry>,
    interner: &GenomeInterner,
    min_ani: f64,
    log_reassign: bool,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
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
        let sample_counts = sample_tag_counts_ref(&group_entries);
        let total_sample_tags = group_entries.len();

        // 为每个基因组（已聚合）计算重新分配后的结果
        for db_entry in db_entries {
            // 最小标签数过滤
            if db_entry.tags.len() < min_tags_genome {
                continue;
            }

            let genome_id = extract_genome_id_from_path(&db_entry.genome_source);
            // 当前基因组在 winner table 中的索引（若从未赢得任何 tag 则为 None）
            let my_idx = interner.get(genome_id);
            let mut covs: Vec<u32> = Vec::new();
            let mut tags_lost_count = 0;

            // 只统计被 winner table 判给当前基因组的 tag 的覆盖度
            for (i, tag) in db_entry.tags.iter().enumerate() {
                let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                let c = match lookup_tag_coverage(*tag, seq, &sample_counts, mismatch) {
                    Some(c) if c > 0 => c,
                    _ => continue,
                };
                // winner table 目前基于 exact tag hash 构建；mismatch 命中时若 exact hash 不在 winner_map 中，
                // 则保守地视为归属于当前基因组（不扣除），避免过度惩罚真实但带 1 mismatch 的 tag。
                match winner_map.get(tag) {
                    Some(winner_entry) if Some(winner_entry.genome_idx) != my_idx => {
                        tags_lost_count += 1; // 该 tag 被分配给了其他基因组
                    }
                    _ => {
                        covs.push(c); // 归当前基因组所有
                    }
                }
            }

            let total_ref_tags = db_entry.tags.len();

            if log_reassign && tags_lost_count > 0 {
                eprintln!("Genome {} in sample {}: owned_tags={}, lost_tags={}, total_ref_tags={}",
                         genome_id, sample_source, covs.len(), tags_lost_count, total_ref_tags);
            }

            // 覆盖度校正后的统计
            let mut result = calculate_statistics(covs, &db_entry.tag_lengths, total_sample_tags, total_ref_tags, mismatch);
            result.sample_file = sample_source.clone();
            result.genome_file = genome_id.to_string();
            result.contig_name = db_entry.sequence_id.clone();

            // 应用profile专用的过滤条件
            if filter_results_for_profile(&result, Some(min_ani), min_shared_tags, min_tags_genome) {
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

// 阈值（对齐 sylph 默认值）
const MIN_ANI: f64 = 90.0;                // sylph query 默认 (was 95)
const PROFILE_MIN_ANI: f64 = 95.0;        // sylph profile 默认 (was 97)
pub const MIN_TAGS_FOR_GENOME: usize = 50;    // 基因组最小标签数 (sylph min_number_kmers)
pub const MIN_SHARED_TAGS_FOR_PROFILE: usize = 100; // 覆盖度校正后报告一个基因组所需的最小共享 tag 数；防止极低 tag 数下的虚假优势调用
// 2bRAD tag 长度，用作 containment->ANI 的指数。BcgI=32；理想情况下应随酶/数据库存储。
const K: f64 = 32.0;

// sylph 覆盖度校正参数（来自 sylph constants.rs / inference.rs）
const SAMPLE_SIZE_CUTOFF: usize = 25;     // 估计 lambda 所需的最小非零覆盖点数
const MEDIAN_ANI_THRESHOLD: f64 = 2.0;    // 中位覆盖 > 此值时朴素 ANI 已足够准确
const MIN_COUNT_CORRECT: f64 = 3.0;       // sylph ratio_lambda 的 min_count_correct 默认值

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
        .filter(|f| f.ends_with(".db"))
        .collect();
    
    let sample_files: Vec<_> = args.files.iter()
        .filter(|f| f.ends_with(".sp"))
        .collect();

    if db_files.is_empty() {
        return Err(anyhow!("No .db files found in input files"));
    }

    if sample_files.is_empty() {
        return Err(anyhow!("No .sp files found in input files"));
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
        // 关键：按基因组聚合所有 contig（去重 tag），避免逐 contig 碎片化
        let db_entries = aggregate_db_by_genome(db_entries);

        eprintln!("Found {} genomes in database (after aggregating contigs)", db_entries.len());
        if let Some(enzyme) = detect_db_enzyme(&db_entries) {
            eprintln!("Database enzyme: {}", enzyme);
        }
        ensure_tag_seqs_for_mismatch(&db_entries, args.mismatch)?;

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

            // 统计样本中每个 tag 的覆盖度（重数），用于覆盖度校正
            let sample_counts = sample_tag_counts(&sample_entries);
            let total_sample_tags = sample_entries.len();
            eprintln!("Total tags in sample: {}", total_sample_tags);

            // 对每个基因组进行比对
            for db_entry in &db_entries {
                // 命中 tag 的覆盖度向量（支持 ≤1 mismatch）
                let covs: Vec<u32> = db_entry.tags.iter().enumerate()
                    .filter_map(|(i, t)| {
                        let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                        lookup_tag_coverage(*t, seq, &sample_counts, args.mismatch)
                    })
                    .collect();

                // 覆盖度校正后的 ANI 统计
                let mut result = calculate_statistics(covs, &db_entry.tag_lengths, total_sample_tags, db_entry.tags.len(), args.mismatch);
                result.sample_file = sample_path.to_string();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();

                // 应用过滤条件
                if filter_results(&result, args.minimum_ani, args.min_tags_genome) {
                    print_result(&result, &writer)?;
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

fn mean_u32(v: &[u32]) -> f64 {
    if v.is_empty() { 0.0 } else { v.iter().map(|&x| x as f64).sum::<f64>() / v.len() as f64 }
}

/// sylph 的 ratio 估计器：从覆盖度直方图估计泊松均值 lambda。
fn ratio_lambda(full_covs: &[u32], min_count_correct: f64) -> Option<f64> {
    let mut num_zero = 0usize;
    let mut count_map: FxHashMap<usize, usize> = FxHashMap::default();
    for &x in full_covs {
        if x == 0 { num_zero += 1; }
        else { *count_map.entry(x as usize).or_insert(0) += 1; }
    }
    if count_map.len() <= 1 { return None; }
    if full_covs.len() - num_zero < SAMPLE_SIZE_CUTOFF { return None; }
    let mut sort_vec: Vec<(usize, usize)> = count_map.iter().map(|(&k, &v)| (v, k)).collect();
    sort_vec.sort_by(|a, b| b.0.cmp(&a.0));
    let most_ind = sort_vec[0].1;
    let count_p1 = *count_map.get(&(most_ind + 1))? as f64;
    let count = count_map[&most_ind] as f64;
    if count_p1 < min_count_correct || count < min_count_correct { return None; }
    Some(count_p1 / count * (most_ind + 1) as f64)
}

/// sylph 的覆盖度校正 ANI：用零截断泊松的检出概率 (1 - e^-lambda) 还原因低覆盖
/// 而漏掉的 tag，再映射回 ANI。这是 meta2bseek 之前缺失、导致大量假阴性的核心。
fn ani_from_lambda(lambda: f64, k: f64, full_covs: &[u32]) -> Option<f64> {
    if full_covs.is_empty() { return None; }
    let contain_count = full_covs.iter().filter(|&&x| x != 0).count() as f64;
    let denom = 1.0 - (-lambda).exp();
    if denom <= 0.0 { return None; }
    let adj_index = contain_count / denom / full_covs.len() as f64;
    let ani = adj_index.powf(1.0 / k);
    if ani.is_nan() || ani < 0.0 { None } else { Some(ani) }
}

fn empty_result(shared: usize, total_ref_tags: usize, query_tags: usize) -> QueryResult {
    QueryResult {
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
        containment_ind: format!("{}/{}", shared, total_ref_tags),
        naive_ani: 0.0,
        ref_tags: total_ref_tags,
        shared_tags: shared,
        query_tags,
        taxonomic_abundance: 0.0,
        sequence_abundance: 0.0,
    }
}

/// 多酶联合 ANI：按 tag 长度分区，信息加权最小二乘估计。
/// `classes` 映射 length -> (N_l, D_l)；`lambda` 为覆盖度校正参数，None 表示不做校正。
fn joint_ani_from_classes(classes: &FxHashMap<u8, (usize, usize)>, lambda: Option<f64>, mismatch: usize) -> Option<f64> {
    if mismatch == 0 {
        // exact matching：c_adj = a^ℓ，对 log 空间做加权最小二乘。
        let mut num = 0.0;
        let mut den = 0.0;
        for (&len, &(n_l, d_l)) in classes {
            if n_l == 0 || d_l == 0 {
                continue;
            }
            let c = d_l as f64 / n_l as f64;
            let c_adj = if let Some(lam) = lambda {
                let p = 1.0 - (-lam).exp();
                if p <= 0.0 {
                    continue;
                }
                (c / p).min(1.0)
            } else {
                c
            };
            let l = len as f64;
            num += n_l as f64 * l * c_adj.ln();
            den += n_l as f64 * l * l;
        }
        if den == 0.0 {
            return None;
        }
        Some((num / den).exp())
    } else {
        // ≤1 mismatch：先除以覆盖度校正因子得到 P_detect(a,ℓ)，再逐长度类反解 a，最后加权平均。
        let mut weighted_a = 0.0;
        let mut total_weight = 0.0;
        for (&len, &(n_l, d_l)) in classes {
            if n_l == 0 || d_l == 0 {
                continue;
            }
            let c = d_l as f64 / n_l as f64;
            let c_adj = if let Some(lam) = lambda {
                let p = 1.0 - (-lam).exp();
                if p <= 0.0 {
                    continue;
                }
                (c / p).min(1.0)
            } else {
                c
            };
            let a_est = ani_from_containment_one_mismatch(c_adj, len as usize);
            let weight = n_l as f64 * len as f64;
            weighted_a += weight * a_est;
            total_weight += weight;
        }
        if total_weight == 0.0 {
            return None;
        }
        Some((weighted_a / total_weight).min(1.0))
    }
}

/// 计算单个基因组的朴素 + 覆盖度校正 ANI（sylph get_stats 风格）。
/// `covs` = 该基因组在样本中命中的 tag 的覆盖度（重数），长度 = 共享 tag 数；
/// `tag_lengths` = 该基因组每个参考 tag 的长度（与完整 tag 集一一对应）；
/// `query_tags` = 样本中总 tag 数；
/// `total_ref_tags` = 该基因组去重后的 tag 总数。
fn calculate_statistics(
    mut covs: Vec<u32>,
    tag_lengths: &[u8],
    query_tags: usize,
    total_ref_tags: usize,
    mismatch: usize,
) -> QueryResult {
    let shared = covs.len();
    if shared == 0 || total_ref_tags == 0 || tag_lengths.len() != total_ref_tags {
        return empty_result(shared, total_ref_tags, query_tags);
    }

    let containment = shared as f64 / total_ref_tags as f64;

    covs.sort_unstable();
    let median_cov = covs[covs.len() / 2] as f64;
    let mean_geq1 = mean_u32(&covs);

    // 完整覆盖向量：命中 tag 的覆盖度 + 漏掉的 tag 记为 0
    let mut full_covs = covs;
    full_covs.extend(std::iter::repeat(0u32).take(total_ref_tags - shared));

    // 按长度统计 (N_l, D_l)
    let mut length_classes: FxHashMap<u8, (usize, usize)> = FxHashMap::default();
    for (i, &len) in tag_lengths.iter().enumerate() {
        let entry = length_classes.entry(len).or_insert((0, 0));
        entry.0 += 1;
        if full_covs[i] > 0 {
            entry.1 += 1;
        }
    }

    // 朴素 ANI：不做覆盖度校正
    let naive_ani = joint_ani_from_classes(&length_classes, None, mismatch)
        .unwrap_or_else(|| {
            if mismatch == 1 {
                // 单长度类时退到主导长度；无长度信息时用 K=32 近似。
                let dominant_len = length_classes.iter()
                    .max_by_key(|(_, (n, _))| *n)
                    .map(|(&l, _)| l as usize)
                    .unwrap_or(K as usize);
                ani_from_containment_one_mismatch(containment, dominant_len)
            } else {
                containment.powf(1.0 / K)
            }
        });

    // 是否需要覆盖度校正
    let (adjusted_ani, eff_cov) = if median_cov > MEDIAN_ANI_THRESHOLD {
        // 覆盖足够高，朴素 containment 已能反映一致性
        (naive_ani, median_cov)
    } else {
        match ratio_lambda(&full_covs, MIN_COUNT_CORRECT) {
            Some(lam) => (
                joint_ani_from_classes(&length_classes, Some(lam), mismatch).unwrap_or(naive_ani),
                lam,
            ),
            None => (naive_ani, mean_geq1), // 点数不足，无法估计 lambda：回退到朴素 ANI
        }
    };
    let adjusted_ani = adjusted_ani.min(1.0); // 截断到 100%

    let mut result = empty_result(shared, total_ref_tags, query_tags);
    result.naive_ani = naive_ani * 100.0;
    result.adjusted_ani = adjusted_ani * 100.0;
    result.eff_cov = eff_cov;
    result.eff_lambda = eff_cov;
    result.median_cov = median_cov;
    result.mean_cov_geq1 = mean_geq1;

    // 简单置信带（仅用于显示）
    let unc = 1.0 + (1.0 - containment) * 1.5;
    result.ani_percentile = ((result.adjusted_ani - unc).max(0.0), (result.adjusted_ani + unc).min(100.0));
    result.lambda_percentile = ((eff_cov - 0.1).max(0.0), eff_cov + 0.1);

    result
}

/// 把每条 contig 的 SyldbEntry 按基因组（genome_source 推导出的 id）合并成一个
/// SyldbEntry，tag 去重。这样 query/profile 都按"整基因组"而非"逐 contig"统计，
/// 这是降低假阴性最关键的一步。
/// 从 DB 条目中推断构建时使用的酶（取非空众数）。
fn detect_db_enzyme(entries: &[SyldbEntry]) -> Option<String> {
    let mut counts: FxHashMap<String, usize> = FxHashMap::default();
    for e in entries {
        if !e.enzyme.is_empty() {
            *counts.entry(e.enzyme.clone()).or_insert(0) += 1;
        }
    }
    counts.into_iter().max_by_key(|(_, c)| *c).map(|(e, _)| e)
}

/// mismatch 模式需要数据库中存储的 tag 序列来生成 1-mismatch 邻居；
/// 用 --no-tag-seqs 构建的（或更旧版本的）数据库不含序列，应明确报错而非静默退化为 exact match。
fn ensure_tag_seqs_for_mismatch(entries: &[SyldbEntry], mismatch: usize) -> Result<()> {
    if mismatch > 0 && entries.iter().any(|e| e.tag_seqs.is_none()) {
        return Err(anyhow!(
            "mismatch mode requires a database built with tag sequences; this database was built with --no-tag-seqs (or predates tag sequence storage). Rebuild the database without --no-tag-seqs or use --mismatch 0"
        ));
    }
    Ok(())
}

fn aggregate_db_by_genome(entries: Vec<SyldbEntry>) -> Vec<SyldbEntry> {
    let mut map: FxHashMap<String, SyldbEntry> = FxHashMap::default();
    for e in entries {
        let gid = extract_genome_id_from_path(&e.genome_source).to_string();
        let agg = map.entry(gid.clone()).or_insert_with(|| SyldbEntry {
            sequence_id: gid.clone(),
            tags: Vec::new(),
            tag_lengths: Vec::new(),
            genome_source: e.genome_source.clone(),
            tag_uniqueness: None,
            tag_seqs: None,
            enzyme: e.enzyme.clone(),
        });
        agg.tags.extend(e.tags);
        agg.tag_lengths.extend(e.tag_lengths);
        if let Some(e_seqs) = e.tag_seqs.as_ref() {
            agg.tag_seqs.get_or_insert_with(Vec::new).extend(e_seqs.iter().cloned());
        }
    }
    let mut out: Vec<SyldbEntry> = map.into_values().collect();
    for e in &mut out {
        // 按 (hash, length, seq) 去重，保证多酶模式下同一 hash 的不同长度 tag 不被错误合并
        let n = e.tags.len();
        let seqs = e.tag_seqs.as_ref().map(|s| s.as_slice());
        let mut triples: Vec<(Hash, u8, TagHash)> = Vec::with_capacity(n);
        for i in 0..n {
            let seq = seqs.map(|s| s.get(i).cloned()).flatten().unwrap_or_default();
            triples.push((e.tags[i], e.tag_lengths[i], seq));
        }
        triples.sort_unstable();
        triples.dedup();
        e.tags = triples.iter().map(|(h, _, _)| *h).collect();
        e.tag_lengths = triples.iter().map(|(_, l, _)| *l).collect();
        if e.tag_seqs.is_some() {
            e.tag_seqs = Some(triples.iter().map(|(_, _, s)| s.clone()).collect());
        }
    }
    out
}

/// 统计样本中每个 tag 的出现次数（覆盖度/重数）。之前的实现把重数压成 1，
/// 导致无法估计覆盖度 lambda。
fn sample_tag_counts(sample_entries: &[SylspEntry]) -> FxHashMap<Hash, u32> {
    let mut counts: FxHashMap<Hash, u32> = FxHashMap::default();
    for e in sample_entries {
        *counts.entry(e.tag).or_insert(0) += 1;
    }
    counts
}

fn sample_tag_counts_ref(sample_entries: &[&SylspEntry]) -> FxHashMap<Hash, u32> {
    let mut counts: FxHashMap<Hash, u32> = FxHashMap::default();
    for e in sample_entries {
        *counts.entry(e.tag).or_insert(0) += 1;
    }
    counts
}

/// 查询一个参考 tag 在样本中的覆盖度。`mismatch=0` 时精确匹配；`mismatch=1` 时额外检查
/// 该参考 tag 的所有 canonical 1-mismatch 邻居。优先返回 exact match 的计数，若不存在则
/// 返回第一个命中的 neighbor 计数（避免同一参考 tag 被多个样本 tag 重复计数）。
fn lookup_tag_coverage(
    tag_hash: Hash,
    tag_seq: Option<&TagHash>,
    sample_counts: &FxHashMap<Hash, u32>,
    mismatch: usize,
) -> Option<u32> {
    if let Some(&c) = sample_counts.get(&tag_hash) {
        if c > 0 {
            return Some(c);
        }
    }
    if mismatch == 0 {
        return None;
    }
    let seq = tag_seq?;
    for neighbor_hash in one_mismatch_canonical_hashes(seq) {
        if let Some(&c) = sample_counts.get(&neighbor_hash) {
            if c > 0 {
                return Some(c);
            }
        }
    }
    None
}

fn filter_results(result: &QueryResult, min_ani: Option<f64>, min_tags_genome: usize) -> bool {
    // 只按 sylph 的方式过滤：有命中、基因组足够大、且（覆盖度校正后的）ANI 达标。
    if result.shared_tags == 0 {
        return false;
    }
    if result.ref_tags < min_tags_genome {
        return false;
    }
    let effective_min_ani = min_ani.unwrap_or(MIN_ANI);
    result.adjusted_ani >= effective_min_ani
}

// profile 专用过滤（阈值更严，但同样依赖覆盖度校正后的 ANI）
fn filter_results_for_profile(result: &QueryResult, min_ani: Option<f64>, min_shared_tags: usize, min_tags_genome: usize) -> bool {
    if result.shared_tags == 0 {
        return false;
    }
    if result.ref_tags < min_tags_genome {
        return false;
    }
    if result.shared_tags < min_shared_tags {
        return false;
    }
    let effective_min_ani = min_ani.unwrap_or(PROFILE_MIN_ANI);
    result.adjusted_ani >= effective_min_ani
}

// 内部函数：使用缓存的数据库数据进行查询 - 优化大文件读取
fn query_single_file_with_cached_db(
    sample_path: &str,
    db_path: &str,
    cached_db_entries: &[SyldbEntry],
    cached_sample_entries: &FxHashMap<String, Vec<SylspEntry>>,
    min_ani: f64,
    mismatch: usize,
    min_shared_tags: usize,
    min_tags_genome: usize,
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
            
            // 统计样本中每个 tag 的覆盖度（重数）
            let sample_counts = sample_tag_counts_ref(entries);
            let total_sample_tags = entries.len();

            // 并行处理每个基因组（已按基因组聚合）进行比对
            cached_db_entries.par_iter().filter_map(|db_entry| {
                // 基因组太小/太碎则跳过（sylph min_number_kmers）
                if db_entry.tags.len() < min_tags_genome {
                    return None;
                }

                let covs: Vec<u32> = db_entry.tags.iter().enumerate()
                    .filter_map(|(i, t)| {
                        let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                        lookup_tag_coverage(*t, seq, &sample_counts, mismatch)
                    })
                    .collect();

                let mut result = calculate_statistics(covs, &db_entry.tag_lengths, total_sample_tags, db_entry.tags.len(), mismatch);
                result.sample_file = sample_source.clone();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();

                // 应用profile专用的过滤条件
                if filter_results_for_profile(&result, Some(min_ani), min_shared_tags, min_tags_genome) {
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
pub fn query_single_file(sample_path: &str, db_path: &str, min_ani: f64, mismatch: usize, min_tags_genome: usize) -> Result<Vec<QueryResult>> {
    eprintln!("Processing database file: {}", db_path);
    
    // 读取数据库文件
    let db_file = File::open(db_path)
        .with_context(|| format!("Failed to open database file: {}", db_path))?;
    let db_reader = BufReader::new(db_file);
    let db_entries: Vec<SyldbEntry> = bincode::deserialize_from(db_reader)
        .with_context(|| format!("Failed to deserialize database file: {}", db_path))?;

    eprintln!("Found {} entries in database", db_entries.len());
    if let Some(enzyme) = detect_db_enzyme(&db_entries) {
        eprintln!("Database enzyme: {}", enzyme);
    }

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
            // 统计样本中每个 tag 的覆盖度（重数）
            let sample_counts = sample_tag_counts_ref(entries);
            let total_sample_tags = entries.len();

            // 并行处理每个基因组进行比对
            db_entries.par_iter().filter_map(|db_entry| {
                let covs: Vec<u32> = db_entry.tags.iter().enumerate()
                    .filter_map(|(i, t)| {
                        let seq = db_entry.tag_seqs.as_ref().and_then(|v| v.get(i));
                        lookup_tag_coverage(*t, seq, &sample_counts, mismatch)
                    })
                    .collect();

                let mut result = calculate_statistics(covs, &db_entry.tag_lengths, total_sample_tags, db_entry.tags.len(), mismatch);
                result.sample_file = sample_source.clone();
                result.genome_file = db_path.to_string();
                result.contig_name = db_entry.sequence_id.clone();

                if filter_results(&result, Some(min_ani), min_tags_genome) {
                    Some(result)
                } else {
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
                   genome_result.common_tags == 0 {
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

// 从db文件中读取基因组映射关系
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
    let min_eff_coverage = args.min_eff_coverage;
    let min_shared_tags = args.min_shared_tags;
    let min_tags_genome = args.min_tags_genome;
    eprintln!("Using minimum ANI threshold: {:.1}%", effective_min_ani);
    eprintln!("Using minimum shared tags threshold: {}", min_shared_tags);
    eprintln!("Using minimum tags per genome threshold: {}", min_tags_genome);
    if min_eff_coverage > 0.0 {
        eprintln!("Using minimum effective coverage threshold: {:.2}", min_eff_coverage);
    }
    
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
    // 关键：按基因组聚合所有 contig（去重 tag），整个 profile 流程按"整基因组"统计
    let cached_db_entries = aggregate_db_by_genome(cached_db_entries);

    eprintln!("Cached {} genomes from database (after aggregating contigs)", cached_db_entries.len());
    if let Some(enzyme) = detect_db_enzyme(&cached_db_entries) {
        eprintln!("Database enzyme: {}", enzyme);
    }
    ensure_tag_seqs_for_mismatch(&cached_db_entries, args.mismatch)?;

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
            if let Ok(initial_results) = query_single_file_with_cached_db(&sample_file, &args.db_file, &cached_db_entries, &cached_sample_entries, effective_min_ani, args.mismatch, min_shared_tags, min_tags_genome) {
                // 按ANI排序
                let mut initial_results = initial_results;
                initial_results.sort_by(|a, b| b.adjusted_ani.partial_cmp(&a.adjusted_ani).unwrap());
                
                // 第二阶段：构建winner table并重新分配（模仿sylph的两阶段处理）
                eprintln!("{} taxonomic profiling; reassigning tags for {} genomes...", &sample_file, initial_results.len());
                
                // 构建winner table
                let (winner_map, interner) = build_winner_table(&initial_results, &cached_db_entries, true); // 启用日志

                // 使用winner table重新计算结果
                if let Some(sample_entries) = cached_sample_entries.get(&sample_file) {
                    let mut reassigned_results = recalculate_with_winner_table(
                        &cached_db_entries,
                        sample_entries,
                        &winner_map,
                        &interner,
                        effective_min_ani,
                        false,
                        args.mismatch,
                        min_shared_tags,
                        min_tags_genome
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
        
        // 过滤掉不符合profile要求的genome（依赖覆盖度校正后的 ANI + 基因组大小 + 最小共享 tag 数 + effective coverage）
        group.retain(|r| {
            r.common_tags > 0 &&
            r.common_tags >= min_shared_tags &&
            r.adjusted_ani >= effective_min_ani &&
            r.total_tags >= min_tags_genome &&
            r.eff_cov >= min_eff_coverage
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

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_result(shared: usize, ref_tags: usize, ani: f64) -> QueryResult {
        QueryResult {
            sample_file: String::new(),
            genome_file: String::new(),
            adjusted_ani: ani,
            eff_cov: 0.0,
            ani_percentile: (0.0, 0.0),
            eff_lambda: 0.0,
            lambda_percentile: (0.0, 0.0),
            median_cov: 0.0,
            mean_cov_geq1: 0.0,
            containment_ind: String::new(),
            naive_ani: 0.0,
            contig_name: String::new(),
            ref_tags,
            shared_tags: shared,
            query_tags: 0,
            taxonomic_abundance: 0.0,
            sequence_abundance: 0.0,
        }
    }

    #[test]
    fn test_min_shared_tag_guard() {
        // A result with enough shared tags and high enough ANI should pass.
        let ok_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(filter_results_for_profile(&ok_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));

        // A result with too few shared tags should be filtered out, even if ANI is high.
        let low_tag_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE - 1, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&low_tag_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));

        // Lowering the threshold should admit the same low-tag result.
        assert!(filter_results_for_profile(&low_tag_result, None, 1, MIN_TAGS_FOR_GENOME));

        // A small genome (few reference tags) is filtered at the default threshold but
        // admitted when min_tags_genome is lowered.
        let small_genome_result = dummy_result(MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME - 1, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&small_genome_result, None, MIN_SHARED_TAGS_FOR_PROFILE, MIN_TAGS_FOR_GENOME));
        assert!(filter_results_for_profile(&small_genome_result, None, MIN_SHARED_TAGS_FOR_PROFILE, 20));

        // A result with zero shared tags should always be filtered out.
        let zero_result = dummy_result(0, MIN_TAGS_FOR_GENOME, PROFILE_MIN_ANI + 1.0);
        assert!(!filter_results_for_profile(&zero_result, None, 1, MIN_TAGS_FOR_GENOME));
    }

    #[test]
    fn test_lookup_tag_coverage_exact_and_mismatch() {
        let mut sample_counts: FxHashMap<Hash, u32> = FxHashMap::default();
        let ref_seq = b"AAACCC";
        let exact_hash = hash_bytes(ref_seq);
        let neighbor_seq = b"AAACCG"; // 1 mismatch
        let neighbor_hash = hash_bytes(neighbor_seq);

        sample_counts.insert(exact_hash, 5);
        sample_counts.insert(neighbor_hash, 3);

        // exact mode：只命中 exact
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 0), Some(5));
        // mismatch mode：优先 exact（5 > 3）
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 1), Some(5));

        // 把 exact 计数设为 0，mismatch mode 应命中 neighbor
        sample_counts.insert(exact_hash, 0);
        assert_eq!(lookup_tag_coverage(exact_hash, Some(&ref_seq.to_vec()), &sample_counts, 1), Some(3));

        // 无 tag_seqs 时 mismatch mode 无法生成邻居，只能 exact match
        assert_eq!(lookup_tag_coverage(exact_hash, None, &sample_counts, 1), None);
    }

    #[test]
    fn test_calculate_statistics_mismatch_increases_shared() {
        // 构造 100 个长度 32 的参考 tag，其中 70 个在样本中被观察到。
        // 在 mismatch 模型下，同样的 observed fraction 对应更低的真实 ANI。
        let total_ref = 100usize;
        let covs = vec![1u32; 70];
        let tag_lengths = vec![32u8; total_ref];

        let exact_result = calculate_statistics(covs.clone(), &tag_lengths, 1000, total_ref, 0);
        let mm_result = calculate_statistics(covs, &tag_lengths, 1000, total_ref, 1);

        assert_eq!(exact_result.shared_tags, mm_result.shared_tags);
        assert!(exact_result.adjusted_ani > mm_result.adjusted_ani,
            "exact ANI {} should be > mismatch ANI {}", exact_result.adjusted_ani, mm_result.adjusted_ani);
    }
}
