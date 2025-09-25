use crate::cmdline::SketchArgs;
use crate::extract::{
    GenomeSketch, get_memory_usage,
};
use anyhow::{Result, Context, anyhow};
use fxhash::{FxHashMap, FxHashSet, FxHasher};
use log::{info, warn, debug};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, BufReader, BufRead};
use std::path::Path;
use std::sync::Mutex;
use std::thread;
use std::time::Duration;
use needletail::parse_fastx_file;
use scalable_cuckoo_filter::ScalableCuckooFilterBuilder;
use bincode;

pub type Hash = u64;
pub type Kmer = u64;
type Marker = u32;

// 从sylph复制的常量和函数
const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

const SAMPLE_FILE_SUFFIX: &str = ".sylsp";
const QUERY_FILE_SUFFIX: &str = ".syldb";
const MAX_DEDUP_COUNT: u32 = 10000;

// 文件格式检查函数
pub fn is_fastq(file: &str) -> bool {
    file.ends_with(".fq")
        || file.ends_with(".fnq")
        || file.ends_with(".fastq")
        || file.ends_with(".fq.gz")
        || file.ends_with(".fnq.gz")
        || file.ends_with(".fastq.gz")
}

pub fn is_fasta(file: &str) -> bool {
    file.ends_with(".fa")
        || file.ends_with(".fna")
        || file.ends_with(".fasta")
        || file.ends_with(".fa.gz")
        || file.ends_with(".fna.gz")
        || file.ends_with(".fasta.gz")
}

// 内存检查和阻塞函数
pub fn check_vram_and_block(max_ram: usize, file: &str) {
    if let Some(usage) = get_memory_usage() {
        let mut gb_usage_curr = usage;
        if (max_ram as f64) < gb_usage_curr {
            debug!(
                "Max memory reached. Blocking sketch for {}. Curr memory {}, max mem {}",
                file, gb_usage_curr, max_ram
            );
        }
        while (max_ram as f64) < gb_usage_curr {
            let one_second = Duration::from_secs(1);
            thread::sleep(one_second);
            if let Some(usage) = get_memory_usage() {
                gb_usage_curr = usage;
                if (max_ram as f64) >= gb_usage_curr {
                    debug!("Sketching for {} freed", file);
                }
            } else {
                break;
            }
        }
    }
}

// Hash函数 (从sylph seeding.rs复制)
#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    let mut key = kmer;
    key = !key.wrapping_add(key << 21);
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8);
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4);
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

// K-mer提取函数
pub fn extract_kmers(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    if string.len() < k {
        return;
    }

    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    let reverse_shift_dist = 2 * (k - 1);
    let mask = u64::MAX >> (std::mem::size_of::<u64>() * 8 - 2 * k);
    let rev_mask = !(3 << (2 * k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);

    // 初始化前 k-1 个核苷酸
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        if nuc_f >= 4 {
            return; // 跳过包含无效核苷酸的序列
        }
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
    }

    // 滑动窗口提取k-mers
    for i in k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        if nuc_f >= 4 {
            continue; // 跳过无效核苷酸
        }
        let nuc_r = 3 - nuc_f;
        
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= mask;
        
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        // 选择canonical k-mer
        let canonical_kmer = if rolling_kmer_f < rolling_kmer_r {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };
        
        let hash_value = mm_hash64(canonical_kmer);
        
        // 基于hash值进行采样
        if hash_value < threshold {
            kmer_vec.push(hash_value);
        }
    }
}

// 提取k-mers及其位置信息
pub fn extract_kmers_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    if string.len() < k {
        return;
    }

    let mut rolling_kmer_f: u64 = 0;
    let mut rolling_kmer_r: u64 = 0;

    let reverse_shift_dist = 2 * (k - 1);
    let mask = u64::MAX >> (std::mem::size_of::<u64>() * 8 - 2 * k);
    let rev_mask = !(3 << (2 * k - 2));
    let len = string.len();
    let threshold = u64::MAX / (c as u64);

    // 初始化前 k-1 个核苷酸
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        if nuc_f >= 4 {
            return; // 跳过包含无效核苷酸的序列
        }
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
    }

    // 滑动窗口提取k-mers
    for i in k - 1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        if nuc_f >= 4 {
            continue; // 跳过无效核苷酸
        }
        let nuc_r = 3 - nuc_f;
        
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= mask;
        
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;

        // 选择canonical k-mer
        let canonical_kmer = if rolling_kmer_f < rolling_kmer_r {
            rolling_kmer_f
        } else {
            rolling_kmer_r
        };
        
        let hash_value = mm_hash64(canonical_kmer);
        
        // 基于hash值进行采样
        if hash_value < threshold {
            kmer_vec.push((contig_number, i, hash_value));
        }
    }
}

// Sketch数据结构
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SequencesSketch {
    pub kmer_counts: FxHashMap<Kmer, u32>,
    pub file_name: String,
    pub c: usize,
    pub k: usize,
    pub paired: bool,
    pub sample_name: Option<String>,
    pub mean_read_length: f64,
}

impl SequencesSketch {
    pub fn new(
        file_name: String,
        c: usize,
        k: usize,
        paired: bool,
        sample_name: Option<String>,
        mean_read_length: f64,
    ) -> Self {
        Self {
            kmer_counts: FxHashMap::default(),
            file_name,
            c,
            k,
            paired,
            sample_name,
            mean_read_length,
        }
    }
}

// 重复去除函数
fn dup_removal_exact(
    kmer_counts: &mut FxHashMap<Kmer, u32>,
    kmer_to_pair_set: &mut FxHashSet<(u64, [Marker; 2])>,
    km: &u64,
    kmer_pair: Option<([Marker; 2], [Marker; 2])>,
    num_dup_removed: &mut usize,
    no_dedup: bool,
    threshold: Option<u32>,
) {
    let c = kmer_counts.entry(*km).or_insert(0);
    let mut c_threshold = u32::MAX;
    if let Some(t) = threshold {
        c_threshold = t;
    }
    
    if !no_dedup && *c < c_threshold {
        if let Some(doublepairs) = kmer_pair {
            let mut ret = false;
            if kmer_to_pair_set.contains(&(*km, doublepairs.0)) {
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert((*km, doublepairs.0));
            }
            if kmer_to_pair_set.contains(&(*km, doublepairs.1)) {
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert((*km, doublepairs.1));
            }
            if ret {
                *num_dup_removed += 1;
                return;
            }
        }
    }
    *c += 1;
}

// pair kmer 函数（用于配对reads去重）
#[inline]
fn pair_kmer_single(s1: &[u8]) -> Option<([Marker; 2], [Marker; 2])> {
    let k = std::mem::size_of::<Marker>() * 4;
    if s1.len() < 4 * k + 2 {
        return None;
    } else {
        let mut kmer_f = 0;
        let mut kmer_g = 0;
        let mut kmer_r = 0;
        let mut kmer_t = 0;
        let halfway = s1.len() / 2;
        
        for i in 0..k {
            let nuc_1 = BYTE_TO_SEQ[s1[2 * i] as usize] as Marker;
            let nuc_2 = BYTE_TO_SEQ[s1[2 * i + halfway] as usize] as Marker;
            let nuc_3 = BYTE_TO_SEQ[s1[1 + 2 * i] as usize] as Marker;
            let nuc_4 = BYTE_TO_SEQ[s1[1 + 2 * i + halfway] as usize] as Marker;

            kmer_f <<= 2;
            kmer_f |= nuc_1;

            kmer_r <<= 2;
            kmer_r |= nuc_2;

            kmer_g <<= 2;
            kmer_g |= nuc_3;

            kmer_t <<= 2;
            kmer_t |= nuc_4;
        }
        return Some(([kmer_f, kmer_r], [kmer_g, kmer_t]));
    }
}

// pair kmer 配对reads
#[inline]
fn pair_kmer(s1: &[u8], s2: &[u8]) -> Option<([Marker; 2], [Marker; 2])> {
    let k = std::mem::size_of::<Marker>() * 4;
    if s1.len() < 2 * k + 1 || s2.len() < 2 * k + 1 {
        return None;
    } else {
        let mut kmer_f = 0;
        let mut kmer_g = 0;
        let mut kmer_r = 0;
        let mut kmer_t = 0;
        
        for i in 0..k {
            let nuc_1 = BYTE_TO_SEQ[s1[2 * i] as usize] as Marker;
            let nuc_2 = BYTE_TO_SEQ[s2[2 * i] as usize] as Marker;
            let nuc_3 = BYTE_TO_SEQ[s1[1 + 2 * i] as usize] as Marker;
            let nuc_4 = BYTE_TO_SEQ[s2[1 + 2 * i] as usize] as Marker;

            kmer_f <<= 2;
            kmer_f |= nuc_1;

            kmer_r <<= 2;
            kmer_r |= nuc_2;

            kmer_g <<= 2;
            kmer_g |= nuc_3;

            kmer_t <<= 2;
            kmer_t |= nuc_4;
        }
        return Some(([kmer_f, kmer_r], [kmer_g, kmer_t]));
    }
}

// 检查参数有效性
fn check_args_valid(args: &SketchArgs) -> Result<()> {
    let level = if args.trace {
        log::LevelFilter::Trace
    } else if args.debug {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .with_context(|| "Failed to initialize thread pool")?;

    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .with_context(|| "Failed to initialize logger")?;

    // 检查是否有输入文件
    if args.genomes.is_none()
        && args.genome_list.is_none()
        && args.reads.is_none()
        && args.sample_list.is_none()
        && args.first_pair.is_empty()
        && args.second_pair.is_empty()
        && args.first_pair_list.is_none()
        && args.second_pair_list.is_none()
        && args.list_sequence.is_none()
        && args.files.is_empty()
    {
        return Err(anyhow!("No input sequences found; see meta2bseek sketch -h for help"));
    }

    if args.fpr < 0. || args.fpr >= 1. {
        return Err(anyhow!("Invalid FPR for sketching. Must be in [0,1)"));
    }

    Ok(())
}

// 解析文件列表
fn parse_line_file(file_name: &str) -> Result<Vec<String>> {
    let file = File::open(file_name)
        .with_context(|| format!("Failed to open file list: {}", file_name))?;
    let reader = BufReader::new(file);
    let mut result = Vec::new();
    
    for line in reader.lines() {
        let line = line.with_context(|| "Failed to read line from file list")?;
        if !line.trim().is_empty() {
            result.push(line.trim().to_string());
        }
    }
    
    Ok(result)
}

// 解析模糊输入文件
fn parse_ambiguous_files(
    args: &SketchArgs,
    read_inputs: &mut Vec<String>,
    genome_inputs: &mut Vec<String>,
) -> Result<()> {
    let mut all_files = vec![];
    
    if let Some(list_file) = &args.list_sequence {
        let files = parse_line_file(list_file)?;
        all_files.extend(files);
    }

    all_files.extend(args.files.clone());

    for file in all_files {
        if is_fastq(&file) {
            read_inputs.push(file);
        } else if is_fasta(&file) {
            genome_inputs.push(file);
        } else {
            warn!(
                "{} does not have a fasta/fastq/gzip type extension; skipping",
                file
            );
        }
    }
    
    Ok(())
}

// 解析reads和genomes
fn parse_reads_and_genomes(
    args: &SketchArgs,
    read_inputs: &mut Vec<String>,
    genome_inputs: &mut Vec<String>,
) -> Result<()> {
    if let Some(genomes) = &args.genomes {
        genome_inputs.extend(genomes.clone());
    }
    
    if let Some(reads) = &args.reads {
        read_inputs.extend(reads.clone());
    }

    if let Some(sample_list) = &args.sample_list {
        let files = parse_line_file(sample_list)?;
        read_inputs.extend(files);
    }

    if let Some(genome_list) = &args.genome_list {
        let files = parse_line_file(genome_list)?;
        genome_inputs.extend(files);
    }
    
    Ok(())
}

// 解析配对reads
fn parse_paired_end_reads(
    args: &SketchArgs,
    first_pairs: &mut Vec<String>,
    second_pairs: &mut Vec<String>,
) -> Result<()> {
    if args.first_pair.len() != args.second_pair.len() {
        return Err(anyhow!("Different number of paired sequences"));
    }

    first_pairs.extend(args.first_pair.clone());
    second_pairs.extend(args.second_pair.clone());

    if let Some(first_list) = &args.first_pair_list {
        let files = parse_line_file(first_list)?;
        first_pairs.extend(files);
    }

    if let Some(second_list) = &args.second_pair_list {
        let files = parse_line_file(second_list)?;
        second_pairs.extend(files);
    }

    if first_pairs.len() != second_pairs.len() {
        return Err(anyhow!("Different number of paired sequences after parsing lists"));
    }
    
    Ok(())
}

// 解析样本名称
fn parse_sample_names(args: &SketchArgs) -> Result<Option<Vec<String>>> {
    if args.list_sample_names.is_none() && args.sample_names.is_none() {
        return Ok(None);
    }
    
    let mut sample_names = Vec::new();
    
    if let Some(file) = &args.list_sample_names {
        let names = parse_line_file(file)?;
        sample_names.extend(names);
    }
    
    if let Some(names) = &args.sample_names {
        sample_names.extend(names.clone());
    }
    
    Ok(Some(sample_names))
}

// sketch单个序列文件
pub fn sketch_sequences_needle(
    read_file: &str,
    c: usize,
    k: usize,
    sample_name: Option<String>,
    no_dedup: bool,
) -> Result<SequencesSketch> {
    let mut kmer_map = HashMap::default();
    let reader = parse_fastx_file(read_file)
        .with_context(|| format!("Failed to parse fastx file: {}", read_file))?;
    
    let mut mean_read_length = 0.;
    let mut counter = 0.;
    let mut kmer_to_pair_table = FxHashSet::default();
    let mut num_dup_removed = 0;

    let mut reader = reader;
    while let Some(record) = reader.next() {
        let record = record.with_context(|| "Failed to read record")?;
        let mut vec = vec![];
        let seq = record.seq();
        
        let kmer_pair = if seq.len() > 400 {
            None
        } else {
            pair_kmer_single(&seq)
        };
        
        extract_kmers(&seq, &mut vec, c, k);
        
        for km in vec {
            dup_removal_exact(
                &mut kmer_map,
                &mut kmer_to_pair_table,
                &km,
                kmer_pair,
                &mut num_dup_removed,
                no_dedup,
                Some(MAX_DEDUP_COUNT),
            );
        }
        
        // 计算移动平均
        counter += 1.;
        mean_read_length = mean_read_length + ((seq.len() as f64) - mean_read_length) / counter;
    }

    let percent = (num_dup_removed as f64) / 
        ((kmer_map.values().sum::<u32>() as f64) + num_dup_removed as f64) * 100.;
    debug!(
        "Number of sketched k-mers removed due to read duplication for {}: {}. Percentage: {:.2}%",
        read_file, num_dup_removed, percent
    );

    Ok(SequencesSketch {
        kmer_counts: kmer_map,
        file_name: read_file.to_string(),
        c,
        k,
        paired: false,
        sample_name,
        mean_read_length,
    })
}

// sketch配对reads
pub fn sketch_pair_sequences(
    read_file1: &str,
    read_file2: &str,
    c: usize,
    k: usize,
    sample_name: Option<String>,
    no_dedup: bool,
    dedup_fpr: f64,
) -> Result<SequencesSketch> {
    let r1o = parse_fastx_file(read_file1)
        .with_context(|| format!("Failed to parse first pair file: {}", read_file1))?;
    let r2o = parse_fastx_file(read_file2)
        .with_context(|| format!("Failed to parse second pair file: {}", read_file2))?;
    
    let mut read_sketch = SequencesSketch::new(
        read_file1.to_string(),
        c,
        k,
        true,
        sample_name,
        0.0
    );

    let mut num_dup_removed = 0;
    let mut reader1 = r1o;
    let mut reader2 = r2o;

    let mut kmer_pair_set = FxHashSet::default();
    let mut fpr = 0.001;
    if dedup_fpr != 0. {
        fpr = dedup_fpr;
    }
    
    let mut _kmer_pair_set_approx = ScalableCuckooFilterBuilder::new()
        .initial_capacity(10_000_000)
        .false_positive_probability(fpr)
        .hasher(FxHasher::default())
        .finish::<(u64, [Marker; 2])>();

    let mut mean_read_length: f64 = 0.;
    let mut counter: f64 = 0.;

    loop {
        let n1 = reader1.next();
        let n2 = reader2.next();
        
        if let (Some(rec1_o), Some(rec2_o)) = (n1, n2) {
            let rec1 = rec1_o.with_context(|| "Failed to read first pair record")?;
            let rec2 = rec2_o.with_context(|| "Failed to read second pair record")?;
            
            let mut temp_vec1 = vec![];
            let mut temp_vec2 = vec![];

            extract_kmers(&rec1.seq(), &mut temp_vec1, c, k);
            extract_kmers(&rec2.seq(), &mut temp_vec2, c, k);
            let kmer_pair = pair_kmer(&rec1.seq(), &rec2.seq());

            // 计算移动平均
            counter += 1.;
            mean_read_length = mean_read_length + 
                ((rec1.seq().len() as f64) - mean_read_length) / counter;

            for km in temp_vec1.iter() {
                dup_removal_exact(
                    &mut read_sketch.kmer_counts,
                    &mut kmer_pair_set,
                    km,
                    kmer_pair,
                    &mut num_dup_removed,
                    no_dedup,
                    None,
                );
            }
            
            for km in temp_vec2.iter() {
                if temp_vec1.contains(km) {
                    continue;
                }
                dup_removal_exact(
                    &mut read_sketch.kmer_counts,
                    &mut kmer_pair_set,
                    km,
                    kmer_pair,
                    &mut num_dup_removed,
                    no_dedup,
                    None,
                );
            }
        } else {
            break;
        }
    }
    
    let percent = (num_dup_removed as f64) / 
        ((read_sketch.kmer_counts.values().sum::<u32>() as f64) + num_dup_removed as f64) * 100.;
    debug!(
        "Number of sketched k-mers removed due to read duplication for {}: {}. Percentage: {:.2}%",
        read_sketch.file_name, num_dup_removed, percent
    );
    
    read_sketch.mean_read_length = mean_read_length;
    Ok(read_sketch)
}

// sketch基因组文件
pub fn sketch_genome(
    c: usize,
    k: usize,
    ref_file: &str,
    min_spacing: usize,
    pseudotax: bool,
) -> Result<GenomeSketch> {
    let reader = parse_fastx_file(ref_file)
        .with_context(|| format!("Failed to parse genome file: {}", ref_file))?;
    
    let mut vec = vec![];
    let mut pseudotax_track_kmers = vec![];
    let mut reader = reader;
    let mut first = true;
    let mut return_genome_sketch = GenomeSketch {
        file_name: ref_file.to_string(),
        first_contig_name: String::new(),
        gn_size: 0,
        c,
        k,
        min_spacing,
        genome_kmers: Vec::new(),
    };
    
    let mut contig_number = 0;
    
    while let Some(record) = reader.next() {
        let record = record.with_context(|| "Failed to read genome record")?;
        
        if first {
            let contig_name = String::from_utf8_lossy(record.id()).to_string();
            return_genome_sketch.first_contig_name = contig_name;
            first = false;
        }
        
        let seq = record.seq();
        return_genome_sketch.gn_size += seq.len();
        extract_kmers_positions(&seq, &mut vec, c, k, contig_number);
        contig_number += 1;
    }
    
    let mut kmer_set = FxHashSet::default();
    let mut duplicate_set = FxHashSet::default();
    let mut new_vec = Vec::with_capacity(vec.len());
    vec.sort();
    
    for (_, _, km) in vec.iter() {
        if !kmer_set.contains(km) {
            kmer_set.insert(*km);
        } else {
            duplicate_set.insert(*km);
        }
    }

    let mut last_pos = 0;
    let mut last_contig = 0;
    for (contig, pos, km) in vec.iter() {
        if !duplicate_set.contains(km) {
            if last_pos == 0 || last_contig != *contig || pos - last_pos > min_spacing {
                new_vec.push(*km);
                last_contig = *contig;
                last_pos = *pos;
            } else if pseudotax {
                pseudotax_track_kmers.push(*km);
            }
        }
    }
    
    return_genome_sketch.genome_kmers = new_vec;
    Ok(return_genome_sketch)
}

// sketch基因组每个contig单独处理
pub fn sketch_genome_individual(
    c: usize,
    k: usize,
    ref_file: &str,
    min_spacing: usize,
    pseudotax: bool,
) -> Result<Vec<GenomeSketch>> {
    let reader = parse_fastx_file(ref_file)
        .with_context(|| format!("Failed to parse genome file: {}", ref_file))?;
    
    let mut reader = reader;
    let mut return_vec = vec![];
    
    while let Some(record) = reader.next() {
        let record = record.with_context(|| "Failed to read genome record")?;
        let mut return_genome_sketch = GenomeSketch {
            file_name: ref_file.to_string(),
            first_contig_name: String::from_utf8_lossy(record.id()).to_string(),
            gn_size: record.seq().len(),
            c,
            k,
            min_spacing,
            genome_kmers: Vec::new(),
        };
        
        let mut pseudotax_track_kmers = vec![];
        let mut kmer_vec = vec![];
        let seq = record.seq();

        extract_kmers_positions(&seq, &mut kmer_vec, c, k, 0);

        let mut kmer_set = FxHashSet::default();
        let mut duplicate_set = FxHashSet::default();
        let mut new_vec = Vec::with_capacity(kmer_vec.len());
        kmer_vec.sort();
        
        for (_, _pos, km) in kmer_vec.iter() {
            if !kmer_set.contains(km) {
                kmer_set.insert(*km);
            } else {
                duplicate_set.insert(*km);
            }
        }
        
        let mut last_pos = 0;
        for (_, pos, km) in kmer_vec.iter() {
            if !duplicate_set.contains(km) {
                if last_pos == 0 || pos - last_pos > min_spacing {
                    new_vec.push(*km);
                    last_pos = *pos;
                } else if pseudotax {
                    pseudotax_track_kmers.push(*km);
                }
            }
        }

        return_genome_sketch.genome_kmers = new_vec;
        return_vec.push(return_genome_sketch);
    }
    
    Ok(return_vec)
}

// 生成合并的样本文件
fn generate_merged_sample_file(
    args: &SketchArgs,
    read_inputs: &[String],
    first_pairs: &[String],
    sample_names: &Option<Vec<String>>,
) -> Result<()> {
    let mut all_sketches = Vec::new();
    
    // 读取所有单端reads的sketch文件
    for (i, read_file) in read_inputs.iter().enumerate() {
        let mut sample_name = None;
        if let Some(names) = sample_names {
            sample_name = Some(names[i + first_pairs.len()].clone());
        }
        
        let sketch_name = if sample_name.is_some() {
            sample_name.as_ref().unwrap()
        } else {
            read_file
        };
        
        let read_file_path = Path::new(sketch_name).file_name().unwrap();
        let file_path = Path::new(&args.sample_output_dir).join(read_file_path);
        let file_path_str = format!("{}{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);
        
        if Path::new(&file_path_str).exists() {
            let file = File::open(&file_path_str)
                .with_context(|| format!("Failed to open sketch file: {}", file_path_str))?;
            let reader = BufReader::new(file);
            let sketch: SequencesSketch = bincode::deserialize_from(reader)
                .with_context(|| format!("Failed to deserialize sketch from: {}", file_path_str))?;
            all_sketches.push(sketch);
        }
    }
    
    // 读取所有配对reads的sketch文件
    for (i, read_file1) in first_pairs.iter().enumerate() {
        let mut sample_name = None;
        if let Some(names) = sample_names {
            sample_name = Some(names[i].clone());
        }
        
        let sketch_name = if sample_name.is_some() {
            sample_name.as_ref().unwrap()
        } else {
            read_file1
        };
        
        let read_file_path = Path::new(sketch_name).file_name().unwrap();
        let file_path = Path::new(&args.sample_output_dir).join(read_file_path);
        let file_path_str = format!("{}.paired{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);
        
        if Path::new(&file_path_str).exists() {
            let file = File::open(&file_path_str)
                .with_context(|| format!("Failed to open paired sketch file: {}", file_path_str))?;
            let reader = BufReader::new(file);
            let sketch: SequencesSketch = bincode::deserialize_from(reader)
                .with_context(|| format!("Failed to deserialize paired sketch from: {}", file_path_str))?;
            all_sketches.push(sketch);
        }
    }
    
    if !all_sketches.is_empty() {
        // 创建合并的sketch文件
        let merged_name = args.out_name.as_deref().unwrap_or("merged_samples");
        let merged_file_path = Path::new(&args.sample_output_dir)
            .join(format!("{}{}", merged_name, SAMPLE_FILE_SUFFIX));
        
        let merged_file = File::create(&merged_file_path)
            .with_context(|| format!("Failed to create merged sample file: {}", merged_file_path.display()))?;
        let mut writer = BufWriter::new(merged_file);
        
        bincode::serialize_into(&mut writer, &all_sketches)
            .with_context(|| "Failed to serialize merged sample sketches")?;
        
        info!("Merged sample file created: {}", merged_file_path.display());
    }
    
    Ok(())
}

// 生成合并的基因组数据库文件
fn generate_merged_genome_file(
    args: &SketchArgs,
    genome_inputs: &[String],
) -> Result<()> {
    let mut all_sketches = Vec::new();
    
    // 读取所有基因组的sketch文件
    for genome_file in genome_inputs.iter() {
        let genome_path = Path::new(genome_file);
        let file_stem = genome_path.file_stem().unwrap().to_str().unwrap();
        
        if args.individual {
            // 对于individual模式，可能有多个文件
            let mut file_index = 0;
            loop {
                let individual_path = Path::new(&args.output_dir)
                    .join(format!("{}_{}{}", file_stem, file_index, QUERY_FILE_SUFFIX));
                
                if individual_path.exists() {
                    let file = File::open(&individual_path)
                        .with_context(|| format!("Failed to open individual genome file: {}", individual_path.display()))?;
                    let reader = BufReader::new(file);
                    let sketches: Vec<GenomeSketch> = bincode::deserialize_from(reader)
                        .with_context(|| format!("Failed to deserialize individual genome sketch from: {}", individual_path.display()))?;
                    all_sketches.extend(sketches);
                    file_index += 1;
                } else {
                    break;
                }
            }
        } else {
            let individual_path = Path::new(&args.output_dir)
                .join(format!("{}{}", file_stem, QUERY_FILE_SUFFIX));
            
            if individual_path.exists() {
                let file = File::open(&individual_path)
                    .with_context(|| format!("Failed to open genome file: {}", individual_path.display()))?;
                let reader = BufReader::new(file);
                let sketches: Vec<GenomeSketch> = bincode::deserialize_from(reader)
                    .with_context(|| format!("Failed to deserialize genome sketch from: {}", individual_path.display()))?;
                all_sketches.extend(sketches);
            }
        }
    }
    
    if !all_sketches.is_empty() {
        // 创建合并的数据库文件
        let merged_name = args.out_name.as_deref().unwrap_or("merged_database");
        let merged_file_path = Path::new(&args.output_dir)
            .join(format!("{}{}", merged_name, QUERY_FILE_SUFFIX));
        
        let merged_file = File::create(&merged_file_path)
            .with_context(|| format!("Failed to create merged genome database file: {}", merged_file_path.display()))?;
        let mut writer = BufWriter::new(merged_file);
        
        bincode::serialize_into(&mut writer, &all_sketches)
            .with_context(|| "Failed to serialize merged genome sketches")?;
        
        info!("Merged genome database file created: {}", merged_file_path.display());
    }
    
    Ok(())
}

// 主sketch函数
pub fn sketch(args: SketchArgs) -> Result<()> {
    let mut read_inputs = vec![];
    let mut genome_inputs = vec![];
    let mut first_pairs = vec![];
    let mut second_pairs = vec![];

    check_args_valid(&args)?;
    parse_ambiguous_files(&args, &mut read_inputs, &mut genome_inputs)?;
    parse_reads_and_genomes(&args, &mut read_inputs, &mut genome_inputs)?;
    parse_paired_end_reads(&args, &mut first_pairs, &mut second_pairs)?;

    let sample_names = parse_sample_names(&args)?;
    if let Some(names) = &sample_names {
        if names.len() != first_pairs.len() + read_inputs.len() {
            return Err(anyhow!(
                "Sample name length is not equal to the number of reads"
            ));
        }
    }

    let mut max_ram = usize::MAX;
    if let Some(ram) = args.max_ram {
        max_ram = ram;
        if max_ram < 7 {
            return Err(anyhow!("Max ram must be >= 7"));
        }
    }

    if genome_inputs.is_empty() && args.db_out_name != "database" {
        warn!("-o is set but no genomes are present. -o only applies to genomes; see -d for reads");
    }

    // 处理配对reads - 只生成单个子文件，不合并
    if !first_pairs.is_empty() && !second_pairs.is_empty() {
        info!("Sketching paired sequences...");
        let iter_vec: Vec<usize> = (0..first_pairs.len()).collect();
        
        iter_vec.into_par_iter().try_for_each(|i| -> Result<()> {
            let read_file1 = &first_pairs[i];
            let read_file2 = &second_pairs[i];

            let mut sample_name = None;
            if let Some(names) = &sample_names {
                sample_name = Some(names[i].clone());
            }
            
            let read_sketch = sketch_pair_sequences(
                read_file1,
                read_file2,
                args.c,
                args.k,
                sample_name.clone(),
                args.no_dedup,
                args.fpr,
            )?;

            // 创建输出目录
            fs::create_dir_all(&args.sample_output_dir)
                .with_context(|| format!("Could not create directory at {}", args.sample_output_dir))?;
            
            let pref = Path::new(&args.sample_output_dir);
            let sketch_name = if sample_name.is_some() {
                read_sketch.sample_name.as_ref().unwrap()
            } else {
                &read_sketch.file_name
            };

            // 生成单个配对文件的子文件
            let read_file_path = Path::new(sketch_name).file_name().unwrap();
            let file_path = pref.join(read_file_path);
            let file_path_str = format!("{}.paired{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);

            let mut read_sk_file = BufWriter::new(
                File::create(&file_path_str)
                    .with_context(|| format!("Failed to create file: {}", file_path_str))?
            );

            bincode::serialize_into(&mut read_sk_file, &read_sketch)
                .with_context(|| "Failed to serialize paired read sketch")?;
            info!("Individual sketching {} complete.", file_path_str);
            
            Ok(())
        })?;
    }

    // 处理单端reads - 只生成单个子文件，不合并
    if !read_inputs.is_empty() {
        info!("Sketching non-paired sequences...");
        let iter_vec: Vec<usize> = (0..read_inputs.len()).collect();
        
        iter_vec.into_par_iter().try_for_each(|i| -> Result<()> {
            let pref = Path::new(&args.sample_output_dir);
            fs::create_dir_all(pref)
                .with_context(|| "Could not create directory for output sample files (-d)")?;

            let read_file = &read_inputs[i];
            check_vram_and_block(max_ram, read_file);
            
            let mut sample_name = None;
            if let Some(names) = &sample_names {
                sample_name = Some(names[i + first_pairs.len()].clone());
            }

            let read_sketch = sketch_sequences_needle(
                read_file,
                args.c,
                args.k,
                sample_name.clone(),
                args.no_dedup,
            )?;

            let sketch_name = if sample_name.is_some() {
                read_sketch.sample_name.as_ref().unwrap()
            } else {
                &read_sketch.file_name
            };
            
            // 生成单个文件的子文件
            let read_file_path = Path::new(sketch_name).file_name().unwrap();
            let file_path = pref.join(read_file_path);
            let file_path_str = format!("{}{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);

            let mut read_sk_file = BufWriter::new(
                File::create(&file_path_str)
                    .with_context(|| format!("Failed to create file: {}", file_path_str))?
            );

            bincode::serialize_into(&mut read_sk_file, &read_sketch)
                .with_context(|| "Failed to serialize read sketch")?;
            info!("Individual sketching {} complete.", file_path_str);
            
            Ok(())
        })?;
    }

    // 处理基因组文件 - 只生成单个子文件，不合并
    if !genome_inputs.is_empty() {
        info!("Sketching genomes...");
        let iter_vec: Vec<usize> = (0..genome_inputs.len()).collect();
        let counter: Mutex<usize> = Mutex::new(0);

        // 创建输出目录
        let output_dir = Path::new(&args.output_dir);
        fs::create_dir_all(output_dir)
            .with_context(|| "Could not create directory for output database files (-o)")?;

        iter_vec.into_par_iter().try_for_each(|i| -> Result<()> {
            let genome_file = &genome_inputs[i];
            
            if args.individual {
                let indiv_gn_sketches = sketch_genome_individual(
                    args.c,
                    args.k,
                    genome_file,
                    args.min_spacing_kmer,
                    !args.no_pseudotax,
                )?;
                
                // 生成单个基因组文件的子文件
                for (j, sketch) in indiv_gn_sketches.iter().enumerate() {
                    let genome_path = Path::new(genome_file);
                    let file_stem = genome_path.file_stem().unwrap().to_str().unwrap();
                    let individual_path = output_dir.join(format!("{}_{}{}", file_stem, j, QUERY_FILE_SUFFIX));
                    
                    let mut individual_file = BufWriter::new(
                        File::create(&individual_path)
                            .with_context(|| format!("Failed to create individual genome file: {}", individual_path.display()))?
                    );
                    bincode::serialize_into(&mut individual_file, &vec![sketch.clone()])
                        .with_context(|| "Failed to serialize individual genome sketch")?;
                    info!("Individual genome sketch {} complete.", individual_path.display());
                }
            } else {
                let genome_sketch = sketch_genome(
                    args.c,
                    args.k,
                    genome_file,
                    args.min_spacing_kmer,
                    !args.no_pseudotax,
                )?;
                
                // 生成单个基因组文件的子文件
                let genome_path = Path::new(genome_file);
                let file_stem = genome_path.file_stem().unwrap().to_str().unwrap();
                let individual_path = output_dir.join(format!("{}{}", file_stem, QUERY_FILE_SUFFIX));
                
                let mut individual_file = BufWriter::new(
                    File::create(&individual_path)
                        .with_context(|| format!("Failed to create individual genome file: {}", individual_path.display()))?
                );
                bincode::serialize_into(&mut individual_file, &vec![genome_sketch.clone()])
                    .with_context(|| "Failed to serialize individual genome sketch")?;
                info!("Individual genome sketch {} complete.", individual_path.display());
            }
            
            let mut c = counter.lock().unwrap();
            *c += 1;
            if *c % 100 == 0 && *c != 0 {
                info!("{} genomes processed.", *c);
            }
            
            Ok(())
        })?;
    }

    // 生成合并的样本文件
    if !read_inputs.is_empty() || !first_pairs.is_empty() {
        info!("Generating merged sample file...");
        generate_merged_sample_file(&args, &read_inputs, &first_pairs, &sample_names)?;
    }

    // 生成合并的基因组数据库文件  
    if !genome_inputs.is_empty() {
        info!("Generating merged genome database file...");
        generate_merged_genome_file(&args, &genome_inputs)?;
    }

    info!("Finished.");
    Ok(())
}
