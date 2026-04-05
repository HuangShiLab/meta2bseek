use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use needletail::parse_fastx_file;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use regex::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    thread,
    time::Duration,
};

use crate::cmdline::ExtractArgs;
use serde::{Serialize, Deserialize};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use crate::constants::{Hash, hash_bytes};
// 添加fxhash导入
use fxhash::{FxHashMap, FxHashSet};

// 添加内存统计导入
use memory_stats::memory_stats;
use log::*;

// 类型别名，与sylph保持一致
pub type TagHash = Vec<u8>;
pub type TagCount = u32;
pub type SampleId = String;
pub type TagFrequencyMap = FxHashMap<TagHash, TagCount>;
pub type SampleStatsMap = FxHashMap<SampleId, ExtractionStats>;

// AVX2相关导入
#[cfg(any(target_arch = "x86_64"))]
use std::arch::x86_64::*;

// AVX2优化的DNA序列匹配函数
#[cfg(any(target_arch = "x86_64"))]
unsafe fn extract_tags_avx2(seq: &[u8], enzyme: &EnzymeSpec) -> Result<Vec<TagHash>> {
    if !is_x86_feature_detected!("avx2") {
        return extract_and_validate_tags(seq, enzyme);
    }

    let seq_str = String::from_utf8_lossy(seq);
    let mut tags = Vec::with_capacity(64);
    let mut seen_tags = FxHashSet::default();

    // 获取酶的标签长度
    let tag_length = ENZYME_TAG_LENGTHS
        .iter()
        .find(|(name, _)| *name == enzyme.name)
        .map(|(_, len)| *len)
        .ok_or_else(|| anyhow::anyhow!("Unknown enzyme: {}", enzyme.name))?;

    // 使用AVX2优化的模式匹配
    for pattern in &enzyme.patterns {
        for m in pattern.find_iter(&seq_str) {
            let matched = m.as_str().as_bytes();
            
            // 使用AVX2优化的序列处理
            let tag = if matched.len() > tag_length {
                let start = (matched.len() - tag_length) / 2;
                let tag_slice = &matched[start..start + tag_length];
                
                // AVX2优化的序列验证
                if is_valid_dna_avx2(tag_slice) {
                    tag_slice.to_vec()
                } else {
                    continue;
                }
            } else {
                // AVX2优化的序列验证
                if is_valid_dna_avx2(matched) {
                    matched.to_vec()
                } else {
                    continue;
                }
            };
            
            // 获取 canonical 版本的 tag
            let canonical_tag = get_canonical_sequence(&tag);
            
            // 使用FxHashSet进行去重
            if seen_tags.insert(canonical_tag.clone()) {
                tags.push(canonical_tag);
            }
        }
    }

    Ok(tags)
}

// AVX2优化的DNA序列验证函数
#[cfg(any(target_arch = "x86_64"))]
unsafe fn is_valid_dna_avx2(seq: &[u8]) -> bool {
    if seq.len() < 32 {
        // 对于短序列，使用标准方法
        return seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'));
    }

    // 使用AVX2进行向量化验证
    let mut i = 0;
    let len = seq.len();
    
    // 处理32字节对齐的部分
    while i + 32 <= len {
        let chunk = _mm256_loadu_si256(seq.as_ptr().add(i) as *const __m256i);
        
        // 检查是否为ACGT
        let a_mask = _mm256_cmpeq_epi8(chunk, _mm256_set1_epi8(b'A' as i8));
        let c_mask = _mm256_cmpeq_epi8(chunk, _mm256_set1_epi8(b'C' as i8));
        let g_mask = _mm256_cmpeq_epi8(chunk, _mm256_set1_epi8(b'G' as i8));
        let t_mask = _mm256_cmpeq_epi8(chunk, _mm256_set1_epi8(b'T' as i8));
        
        let valid_mask = _mm256_or_si256(
            _mm256_or_si256(a_mask, c_mask),
            _mm256_or_si256(g_mask, t_mask)
        );
        
        let result = _mm256_movemask_epi8(valid_mask);
        // 修复：使用正确的常量，0xFFFFFFFFu32 as i32 等于 -1
        if result != -1i32 {
            return false;
        }
        
        i += 32;
    }
    
    // 处理剩余部分
    for &b in &seq[i..] {
        if !matches!(b, b'A' | b'C' | b'G' | b'T') {
            return false;
        }
    }
    
    true
}



// 优化的压缩设置
fn get_optimal_compression() -> Compression {
    // 根据系统性能调整压缩级别
    if std::env::var("FAST_COMPRESSION").is_ok() {
        Compression::fast()
    } else if std::env::var("BEST_COMPRESSION").is_ok() {
        Compression::best()
    } else {
        Compression::default()
    }
}

// 优化的文件大小检测
fn get_file_size_optimized(path: &Path) -> Result<u64> {
    let metadata = std::fs::metadata(path)?;
    Ok(metadata.len())
}

// 优化的缓冲区大小计算
fn calculate_optimal_buffer_size(file_size: u64, is_compressed: bool) -> usize {
    let base_size = if file_size > 1024 * 1024 * 1024 {
        // 大文件：256KB
        256 * 1024
    } else if file_size > 100 * 1024 * 1024 {
        // 中等文件：128KB
        128 * 1024
    } else {
        // 小文件：64KB
        64 * 1024
    };
    
    if is_compressed {
        base_size * 2
    } else {
        base_size
    }
}



// 内存监控函数，参考sketch中的check_vram_and_block
pub fn check_vram_and_block(max_ram: usize, file: &str) {
    if let Some(usage) = memory_stats() {
        let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
        if (max_ram as f64) < gb_usage_curr {
            log::debug!(
                "Max memory reached. Blocking extract for {}. Curr memory {}, max mem {}",
                file,
                gb_usage_curr,
                max_ram
            );
        }
        while (max_ram as f64) < gb_usage_curr {
            let five_second = Duration::from_secs(1);
            thread::sleep(five_second);
            if let Some(usage) = memory_stats() {
                gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
                if (max_ram as f64) >= gb_usage_curr {
                    log::debug!("Extract for {} freed", file);
                }
            } else {
                break;
            }
        }
    }
}

// 动态内存管理函数
pub fn get_memory_usage() -> Option<f64> {
    memory_stats().map(|usage| usage.virtual_mem as f64 / 1_000_000_000 as f64)
}







// 内存安全的处理函数
pub fn safe_process_with_memory_check<F, T>(
    max_ram: usize,
    file: &str,
    process_fn: F,
) -> Result<T>
where
    F: FnOnce() -> Result<T>,
{
    // 检查当前内存使用
    if let Some(current_memory) = get_memory_usage() {
        if current_memory > max_ram as f64 {
            check_vram_and_block(max_ram, file);
        }
    }
    
    // 执行处理函数
    let result = process_fn()?;
    
    // 处理完成后再次检查内存
    if let Some(current_memory) = get_memory_usage() {
        if current_memory > max_ram as f64 * 0.8 {
            // 如果内存使用超过80%，等待一下让系统回收内存
            thread::sleep(Duration::from_millis(100));
        }
    }
    
    Ok(result)
}

pub const ENZYME_DEFINITIONS: &[(&str, &[&str])] = &[
    ("CspCI", &[
        r"[ACGT]{11}CAA[ACGT]{5}GTGG[ACGT]{10}",
        r"[ACGT]{10}CCAC[ACGT]{5}TTG[ACGT]{11}",
    ]),
    ("AloI", &[
        r"[ACGT]{7}GAAC[ACGT]{6}TCC[ACGT]{7}",
        r"[ACGT]{7}GGA[ACGT]{6}GTTC[ACGT]{7}",
    ]),
    ("BsaXI", &[
        r"[ACGT]{9}AC[ACGT]{5}CTCC[ACGT]{7}",
        r"[ACGT]{7}GGAG[ACGT]{5}GT[ACGT]{9}",
    ]),
    ("BaeI", &[
        r"[ACGT]{10}AC[ACGT]{4}GTA[CT]C[ACGT]{7}",
        r"[ACGT]{7}G[AG]TAC[ACGT]{4}GT[ACGT]{10}",
    ]),
    ("BcgI", &[
        r"[ACGT]{10}CGA[ACGT]{6}TGC[ACGT]{10}",
        r"[ACGT]{10}GCA[ACGT]{6}TCG[ACGT]{10}",
    ]),
    ("CjeI", &[
        r"[ACGT]{8}CCA[ACGT]{6}GT[ACGT]{9}",
        r"[ACGT]{9}AC[ACGT]{6}TGG[ACGT]{8}",
    ]),
    ("PpiI", &[
        r"[ACGT]{7}GAAC[ACGT]{5}CTC[ACGT]{8}",
        r"[ACGT]{8}GAG[ACGT]{5}GTTC[ACGT]{7}",
    ]),
    ("PsrI", &[
        r"[ACGT]{7}GAAC[ACGT]{6}TAC[ACGT]{7}",
        r"[ACGT]{7}GTA[ACGT]{6}GTTC[ACGT]{7}",
    ]),
    ("BplI", &[
        r"[ACGT]{8}GAG[ACGT]{5}CTC[ACGT]{8}",
    ]),
    ("FalI", &[
        r"[ACGT]{8}AAG[ACGT]{5}CTT[ACGT]{8}",
    ]),
    ("Bsp24I", &[
        r"[ACGT]{8}GAC[ACGT]{6}TGG[ACGT]{7}",
        r"[ACGT]{7}CCA[ACGT]{6}GTC[ACGT]{8}",
    ]),
    ("HaeIV", &[
        r"[ACGT]{7}GA[CT][ACGT]{5}[AG]TC[ACGT]{9}",
        r"[ACGT]{9}GA[CT][ACGT]{5}[AG]TC[ACGT]{7}",
    ]),
    ("CjePI", &[
        r"[ACGT]{7}CCA[ACGT]{7}TC[ACGT]{8}",
        r"[ACGT]{8}GA[ACGT]{7}TGG[ACGT]{7}",
    ]),
    ("Hin4I", &[
        r"[ACGT]{8}GA[CT][ACGT]{5}[GAC]TC[ACGT]{8}",
        r"[ACGT]{8}GA[CTG][ACGT]{5}[AG]TC[ACGT]{8}",
    ]),
    ("AlfI", &[
        r"[ACGT]{10}GCA[ACGT]{6}TGC[ACGT]{10}",
    ]),
    ("BslFI", &[
        r"[ACGT]{6}GGGAC[ACGT]{14}",
        r"[ACGT]{14}GTCCC[ACGT]{6}",
    ]),
];

// 定义每个内切酶的标签长度（固定匹配碱基数 + 自由匹配碱基数）
pub const ENZYME_TAG_LENGTHS: &[(&str, usize)] = &[
    ("CspCI", 33),  // 11 + 3 + 5 + 4 + 10 = 33
    ("AloI", 20),   // 7 + 4 + 6 + 3 = 20
    ("BsaXI", 27),  // 9 + 2 + 5 + 4 + 7 = 27
    ("BaeI", 27),   // 10 + 2 + 4 + 4 + 7 = 27
    ("BcgI", 32),   // 10 + 3 + 6 + 3 + 10 = 32
    ("CjeI", 28),   // 8 + 3 + 6 + 2 + 9 = 28
    ("PpiI", 27),   // 7 + 4 + 5 + 3 + 8 = 27
    ("PsrI", 27),   // 7 + 4 + 6 + 3 + 7 = 27
    ("BplI", 27),   // 8 + 3 + 5 + 3 + 8 = 27
    ("FalI", 27),   // 8 + 3 + 5 + 3 + 8 = 27
    ("Bsp24I", 27), // 8 + 3 + 6 + 3 + 7 = 27
    ("HaeIV", 25),  // 7 + 2 + 5 + 2 + 9 = 25
    ("CjePI", 27),  // 7 + 3 + 7 + 2 + 8 = 27
    ("Hin4I", 25),  // 8 + 2 + 5 + 2 + 8 = 25
    ("AlfI", 32),   // 10 + 3 + 6 + 3 + 10 = 32
    ("BslFI", 25),  // 6 + 5 + 14 = 25
];

#[derive(Debug)]
pub struct EnzymeSpec {
    pub name: String,
    pub patterns: Vec<Regex>,
}

impl EnzymeSpec {
    pub fn new(name: &str) -> Result<Self> {
        let def = ENZYME_DEFINITIONS
            .iter()
            .find(|(e, _)| *e == name)
            .ok_or_else(|| anyhow::anyhow!("Unsupported enzyme: {}", name))?;

        let patterns = def.1
            .iter()
            .map(|p| Regex::new(p).context(format!("Invalid regex pattern: {}", p)))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            name: def.0.to_string(),
            patterns,
        })
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SyldbEntry {
    pub sequence_id: String,
    pub tags: Vec<Hash>,
    pub positions: Vec<usize>,
    pub genome_source: String,
    // 新增字段：标记每个tag是否为unique（taxa-specific）
    pub tag_uniqueness: Option<Vec<bool>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SylspEntry {
    pub sequence_id: String,
    pub tag: Hash,
    pub quality: Option<String>,
    pub sample_source: String,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketch {
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub c: usize,
    pub k: usize,
    pub min_spacing: usize,
    pub genome_kmers: Vec<Hash>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketchInspect {
    pub file_name: String,
    pub genome_kmers_num: usize,
    pub first_contig_name: String,
    pub genome_size: usize,
}

impl From<GenomeSketch> for GenomeSketchInspect {
    fn from(sk: GenomeSketch) -> Self {
        GenomeSketchInspect {
            genome_kmers_num: sk.genome_kmers.len(),
            file_name: sk.file_name,
            first_contig_name: sk.first_contig_name,
            genome_size: sk.gn_size,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct DatabaseSketch {
    pub database_file: String,
    pub c: usize,
    pub k: usize,
    pub min_spacing_parameter: usize,
    pub genome_files: Vec<GenomeSketchInspect>,
}

#[allow(dead_code)]
pub fn process_input(
    input_files: Vec<PathBuf>,
    sample_output_dir: &Path,
    enzyme_name: &str,
    _threads: usize,
    format: &str,
) -> Result<()> {
    let enzyme = EnzymeSpec::new(enzyme_name)
        .context(format!("Unsupported enzyme: {}", enzyme_name))?;

    for input_path in &input_files {
        // 确定输入文件类型
        let is_fasta = is_fasta_file(input_path)
            .context("Failed to determine if file is FASTA")?;
        let is_fastq = is_fastq_file(input_path)
            .context("Failed to determine if file is FASTQ")?;

        if !is_fasta && !is_fastq {
            return Err(anyhow::anyhow!("Unsupported file format: {}", input_path.display()));
        }

        let file_stem = input_path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("output");
            
        let output_name = if format == "fq" {
            format!("{}.fq", file_stem)
        } else {
            format!("{}.fa", file_stem)
        };
        
        let mut output_path = PathBuf::from(sample_output_dir);
        output_path.push(output_name);

        // 根据文件类型处理
        if is_fasta {
            process_fasta(input_path, &output_path, &enzyme, format, input_path.to_string_lossy().ends_with(".gz"))?;
        } else {
            process_fastq(input_path, &output_path, &enzyme, format, input_path.to_string_lossy().ends_with(".gz"))?;
        }
    }

    Ok(())
}



#[allow(dead_code)]
fn is_fasta_file(path: &Path) -> Result<bool> {
    // 检查文件扩展名
    let ext = path.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase();

    // 如果是压缩文件，获取原始扩展名
    let base_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    let is_fasta_ext = if base_name.ends_with(".gz") {
        // 移除 .gz 后缀并检查原始扩展名
        let without_gz = base_name.trim_end_matches(".gz");
        without_gz.ends_with(".fa") || 
        without_gz.ends_with(".fasta") || 
        without_gz.ends_with(".fna") || 
        without_gz.ends_with(".ffn") || 
        without_gz.ends_with(".faa") || 
        without_gz.ends_with(".frn")
    } else {
        matches!(ext.as_str(), 
            "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn"
        )
    };

    // 如果扩展名不明确，检查文件内容
    if !is_fasta_ext {
        let mut reader = create_reader(path)?;
        let mut first_char = [0u8; 1];
        if reader.read_exact(&mut first_char).is_ok() {
            return Ok(first_char[0] == b'>');
        }
    }

    Ok(is_fasta_ext)
}

#[allow(dead_code)]
fn is_fastq_file(path: &Path) -> Result<bool> {
    // 检查文件扩展名
    let base_name = path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    
    let is_fastq_ext = if base_name.ends_with(".gz") {
        // 移除 .gz 后缀并检查原始扩展名
        let without_gz = base_name.trim_end_matches(".gz");
        without_gz.ends_with(".fq") || 
        without_gz.ends_with(".fastq")
    } else {
        base_name.ends_with(".fq") || 
        base_name.ends_with(".fastq")
    };

    // 如果扩展名不明确，检查文件内容
    if !is_fastq_ext {
        let mut reader = create_reader(path)?;
        let mut first_char = [0u8; 1];
        if reader.read_exact(&mut first_char).is_ok() {
            return Ok(first_char[0] == b'@');
        }
    }

    Ok(is_fastq_ext)
}

// 完全按照sylph方式处理FASTA文件
fn process_fasta_sylph_style(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    
    // 完全按照sylph的模式
    let reader = parse_fastx_file(input);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", input.display());
        return Ok(());
    }
    
    let mut reader = reader.unwrap();
    let mut kmer_to_tag_table = FxHashSet::default();
    
    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());
            
            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();
            
            // 使用sylph风格的标签提取（现在包含canonical处理）
            let tags = extract_and_validate_tags(&seq, enzyme)
                .context(format!("Failed to process sequence: {}", seq_id))?;
            
            // 按照sylph的去重模式（现在使用canonical tags）
            for tag in tags {
                if kmer_to_tag_table.insert(tag.clone()) {
                    stats.total_tags += 1;
                    write_tags(&mut *writer, &seq_id, &[tag], format)?;
                }
            }
        } else {
            warn!("Invalid record in file {}", input.display());
        }
    }
    
    log_stats(stats, enzyme);
    Ok(())
}



fn process_fasta(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    // 直接使用sylph风格的处理
    process_fasta_sylph_style(input, output, enzyme, format, compress)
}

// 按照sylph风格处理FASTQ文件
fn process_fastq_sylph_style(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    
    // 完全按照sylph的模式
    let reader = parse_fastx_file(input);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", input.display());
        return Ok(());
    }
    
    let mut reader = reader.unwrap();
    let mut kmer_to_tag_table = FxHashSet::default();
    
    while let Some(record) = reader.next() {
        if record.is_ok() {
            let record = record.expect(&format!("Invalid record for file {} ", input.display()));
            let seq = record.seq();
            let seq_id = String::from_utf8_lossy(record.id());
            
            stats.total_sequences += 1;
            stats.total_sequence_length += seq.len();
            
            // 使用sylph风格的标签提取（现在包含canonical处理）
            let tags = extract_and_validate_tags(&seq, enzyme)
                .context(format!("Failed to process sequence: {}", seq_id))?;
            
            // 按照sylph的去重模式（现在使用canonical tags）
            for tag in tags {
                if kmer_to_tag_table.insert(tag.clone()) {
                    stats.total_tags += 1;
                    write_tags(&mut *writer, &seq_id, &[tag], format)?;
                }
            }
        } else {
            warn!("Invalid record in file {}", input.display());
        }
    }
    
    log_stats(stats, enzyme);
    Ok(())
}

fn process_fastq(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    // 直接使用sylph风格的处理
    process_fastq_sylph_style(input, output, enzyme, format, compress)
}

fn extract_and_validate_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Result<Vec<TagHash>> {
    #[cfg(any(target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe {
                return extract_tags_avx2(seq, enzyme);
            }
        }
    }
    
    // 标准实现（非AVX2或非x86_64架构）
    let seq_str = String::from_utf8_lossy(seq);
    // 预估每个序列可能产生的标签数量，减少重新分配
    let mut tags = Vec::with_capacity(64);
    // 使用FxHashSet进行去重，提高性能
    let mut seen_tags = FxHashSet::default();

    // 获取酶的标签长度
    let tag_length = ENZYME_TAG_LENGTHS
        .iter()
        .find(|(name, _)| *name == enzyme.name)
        .map(|(_, len)| *len)
        .ok_or_else(|| anyhow::anyhow!("Unknown enzyme: {}", enzyme.name))?;

    for pattern in &enzyme.patterns {
        for m in pattern.find_iter(&seq_str) {
            let matched = m.as_str().as_bytes();
            // 只保留酶切位点之间的序列
            let tag = if matched.len() > tag_length {
                let start = (matched.len() - tag_length) / 2;
                matched[start..start + tag_length].to_vec()
            } else {
                matched.to_vec()
            };
            
            // 获取 canonical 版本的 tag
            let canonical_tag = get_canonical_sequence(&tag);
            
            // 使用FxHashSet进行去重
            if seen_tags.insert(canonical_tag.clone()) {
                tags.push(canonical_tag);
            }
        }
    }

    Ok(tags)
}

fn write_tags(
    writer: &mut dyn Write,
    seq_id: &str,
    tags: &[TagHash],
    format: &str,
) -> Result<()> {
    for (i, tag) in tags.iter().enumerate() {
        let header = format!("{}_tag{}", seq_id, i + 1);
        match format {
            "fq" => writeln!(writer, "@{}\n{}\n+\n{}", 
                           header, 
                           String::from_utf8_lossy(tag),
                           "~".repeat(tag.len()))
                .context("Failed to write FASTQ record")?,
            _ => writeln!(writer, ">{}\n{}", 
                        header, 
                        String::from_utf8_lossy(tag))
                .context("Failed to write FASTA record")?,
        }
    }
    Ok(())
}

fn create_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .context(format!("Failed to open input file: {}", path.display()))?;

    // 使用优化的文件大小检测
    let file_size = get_file_size_optimized(path)?;
    let is_compressed = path.to_string_lossy().ends_with(".gz");
    
    // 使用优化的缓冲区大小计算
    let buffer_size = calculate_optimal_buffer_size(file_size, is_compressed);

    Ok(if is_compressed {
        Box::new(BufReader::with_capacity(buffer_size, GzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(buffer_size, file))
    })
}

fn create_writer(path: &Path, compress: bool) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .context(format!("Failed to create output file: {}", path.display()))?;

    // 使用优化的缓冲区大小和压缩设置
    let buffer_size = if compress {
        256 * 1024
    } else {
        128 * 1024
    };

    Ok(if compress {
        let compression = get_optimal_compression();
        Box::new(BufWriter::with_capacity(buffer_size, GzEncoder::new(file, compression)))
    } else {
        Box::new(BufWriter::with_capacity(buffer_size, file))
    })
}

#[derive(Debug, Clone)]
pub struct ExtractionStats {
    total_sequences: usize,
    total_tags: usize,
    total_sequence_length: usize,
}



impl ExtractionStats {
    fn new() -> Self {
        Self {
            total_sequences: 0,
            total_tags: 0,
            total_sequence_length: 0,
        }
    }
}

fn log_stats(stats: ExtractionStats, enzyme: &EnzymeSpec) {
    let k = enzyme.patterns[0].as_str().len();
    let total_kmers = if stats.total_sequence_length >= (k - 1) * stats.total_sequences {
        stats.total_sequence_length - (k - 1) * stats.total_sequences
    } else {
        0
    };
    let percentage = calculate_tag_percentage(stats.total_tags, total_kmers);
    
    // 获取酶的标签长度
    let tag_length = ENZYME_TAG_LENGTHS
        .iter()
        .find(|(name, _)| *name == enzyme.name)
        .map(|(_, len)| *len)
        .unwrap_or(k); // 如果找不到对应的长度，使用模式长度作为后备
    
    let tag_bases_percentage = (stats.total_tags * tag_length) as f64 / stats.total_sequence_length as f64 * 100.0;
    
    println!(
        "\nProcessing complete for {}:\n\
        =============================\n\
        - Total sequences processed: {}\n\
        - Total sequence length: {}\n\
        - Average sequence length: {:.2}\n\
        - Total tags extracted: {}\n\
        - Average tags per sequence: {:.2}\n\
        - Extractable k-mers: {}\n\
        - 2bRAD tag percentage: {:.4}%\n\
        - 2bRAD tag bases percentage: {:.4}%\n\
        - Recognition patterns used: {}",
        enzyme.name,
        stats.total_sequences,
        stats.total_sequence_length,
        stats.total_sequence_length as f32 / stats.total_sequences.max(1) as f32,
        stats.total_tags,
        stats.total_tags as f32 / stats.total_sequences.max(1) as f32,
        total_kmers,
        percentage,
        tag_bases_percentage,
        enzyme.patterns
            .iter()
            .map(|r| r.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    );
}

// 新增函数：处理单对双端测序文件
fn process_paired_fastq_files(
    first_file: &str,
    second_file: &str,
    enzyme: &EnzymeSpec,
    _sample_output_dir: &Path,
    _out_name: Option<&str>,
) -> Result<()> {
    // 从文件名中提取样本名
    let file_stem = Path::new(first_file)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .split('.')
        .next()
        .unwrap_or("unknown")
        .to_string();

    // 处理一对文件
    let fa_entries = process_paired_fastq_to_sylsp(
        first_file,
        second_file,
        enzyme,
        &file_stem,
    )?;

    // 注释掉生成单个文件的代码 - 只保留合并后的文件
    // let output_base = Path::new(sample_output_dir).join(&file_stem);
    // let fa_path = output_base.with_extension("fa");
    // let mut fa_writer = create_writer(&fa_path, false)?;

    let mut sylsp_entries = Vec::new();
    for (id, tag, sample_source) in &fa_entries {
        // 注释掉单个FASTA文件写入
        // writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(tag))
        //     .context("Failed to write FASTA record")?;

        let entry = SylspEntry {
            sequence_id: id.clone(),
            tag: hash_bytes(tag),
            quality: None,
            sample_source: sample_source.clone(),
        };
        sylsp_entries.push(entry.clone());
    }

    // 注释掉生成单个sylsp文件的代码
    // let sylsp_path = if let Some(name) = out_name {
    //     Path::new(sample_output_dir).join(format!("{}.sylsp", name))
    // } else {
    //     output_base.with_extension("sylsp")
    // };
    // 
    // let sylsp_file = File::create(&sylsp_path)
    //     .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
    // let sylsp_writer = BufWriter::new(sylsp_file);
    // bincode::serialize_into(sylsp_writer, &sylsp_entries)
    //     .context("Failed to serialize sylsp data")?;

    Ok(())
}

pub fn extract(args: ExtractArgs) -> Result<()> {
    // 初始化线程池
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // 创建输出目录
    std::fs::create_dir_all(&args.sample_output_dir)
        .context("Failed to create output directory")?;

    // 设置内存限制，如果没有指定则使用默认值
    let max_ram = args.max_ram.unwrap_or(16); // 默认16GB内存限制
    if max_ram < 7 {
        return Err(anyhow::anyhow!("Max ram must be >= 7. Exiting."));
    }

    // 处理单对双端测序文件（-1 和 -2 参数）
    if !args.first_pair.is_empty() && !args.second_pair.is_empty() {
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        for (first_file, second_file) in args.first_pair.iter().zip(args.second_pair.iter()) {
            safe_process_with_memory_check(max_ram, first_file, || {
                process_paired_fastq_files(
                    first_file,
                    second_file,
                    &enzyme,
                    Path::new(&args.sample_output_dir),
                    args.out_name.as_deref(),
                )
            })?;
        }
    }

    // 处理批处理双端测序文件（--l1 和 --l2 参数）
    if let (Some(first_pair_list), Some(second_pair_list)) = (&args.first_pair_list, &args.second_pair_list) {
        // 读取文件列表
        let first_pairs = read_file_list(first_pair_list)
            .context("Failed to read first pair list")?;
        let second_pairs = read_file_list(second_pair_list)
            .context("Failed to read second pair list")?;

        if first_pairs.len() != second_pairs.len() {
            return Err(anyhow::anyhow!("Number of files in first pair list and second pair list do not match"));
        }

        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_sylsp_entries = Vec::new();

        // 并行处理所有配对文件，添加内存监控
        let results: Vec<Result<(String, Vec<SylspEntry>)>> = first_pairs.par_iter()
            .zip(second_pairs.par_iter())
            .map(|(first_file, second_file)| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, first_file);
                    }
                }
                
                let input_path1 = PathBuf::from(first_file);
                let file_stem = input_path1.file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .split('.')
                    .next()
                    .unwrap_or("unknown")
                    .to_string();

                // 处理一对文件
                let fa_entries = process_paired_fastq_to_sylsp(
                    first_file,
                    second_file,
                    &enzyme,
                    &file_stem,
                )?;

                // 注释掉生成单个文件的代码 - 只保留合并后的文件
                // let output_base = Path::new(&args.sample_output_dir).join(&file_stem);
                // let fa_path = output_base.with_extension("fa");
                // let mut fa_writer = create_writer(&fa_path, false)?;

                let mut sylsp_entries = Vec::new();
                for (id, tag, sample_source) in &fa_entries {
                    // 注释掉单个FASTA文件写入
                    // writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(tag))
                    //     .context("Failed to write FASTA record")?;

                    let entry = SylspEntry {
                        sequence_id: id.clone(),
                        tag: hash_bytes(tag),
                        quality: None,
                        sample_source: sample_source.clone(),
                    };
                    sylsp_entries.push(entry.clone());
                }

                // 注释掉生成单个sylsp文件的代码
                // let sylsp_path = output_base.with_extension("sylsp");
                // let sylsp_file = File::create(&sylsp_path)
                //     .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
                // let sylsp_writer = BufWriter::new(sylsp_file);
                // bincode::serialize_into(sylsp_writer, &sylsp_entries)
                //     .context("Failed to serialize sylsp data")?;

                Ok((file_stem, sylsp_entries))
            })
            .collect();

        // 处理结果并收集所有 sylsp 条目
        for result in results {
            match result {
                Ok((_, entries)) => {
                    all_sylsp_entries.extend(entries);
                },
                Err(e) => eprintln!("Error processing paired files: {}", e),
            }
        }

        // 生成合并的 sylsp 文件
        if !all_sylsp_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
            let combined_sylsp_file = File::create(&combined_sylsp_path)
                .context(format!("Failed to create combined sylsp file: {}", combined_sylsp_path.display()))?;
            let combined_sylsp_writer = BufWriter::new(combined_sylsp_file);
            
            bincode::serialize_into(combined_sylsp_writer, &all_sylsp_entries)
                .context("Failed to serialize combined sylsp data")?;
        }
    }

    // 处理单端测序文件
    if let Some(read_files) = args.reads {
        // 存储所有 FASTQ 文件的 sylsp 条目
        let mut all_sylsp_entries = Vec::new();
        let mut all_fa_entries = Vec::new();
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        
        for file in read_files {
            // 检查内存使用
            if let Some(current_memory) = get_memory_usage() {
                if current_memory > max_ram as f64 {
                    check_vram_and_block(max_ram, &file);
                }
            }
            
            let input_path = PathBuf::from(&file);
            let file_stem = input_path.file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .split('.')
                .next()
                .unwrap_or("unknown")
                .to_string();
                
            let reader = fastq::Reader::new(create_reader(&input_path)?);
            let mut stats = ExtractionStats::new();

            for result in reader.records() {
                let record = result.context("Failed to read FASTQ record")?;
                stats.total_sequences += 1;
                stats.total_sequence_length += record.seq().len();
                
                let tags = extract_and_validate_tags(record.seq(), &enzyme)
                    .context(format!("Failed to process read: {}", record.id()))?;
                    
                for (i, tag) in tags.iter().enumerate() {
                    let id = format!("{}_tag{}", record.id(), i + 1);
                    all_fa_entries.push((id.clone(), tag.clone()));
                    
                    let entry = SylspEntry {
                        sequence_id: id,
                        tag: hash_bytes(tag),
                        quality: Some(String::from_utf8_lossy(record.qual()).to_string()),
                        sample_source: file_stem.clone(),
                    };
                    all_sylsp_entries.push(entry);
                }
                
                stats.total_tags += tags.len();
            }
            
            log_stats(stats, &enzyme);
        }
        
        // 生成合并的输出文件
        let output_name = args.out_name.as_ref().map_or_else(|| "reads".to_string(), |s| s.clone());
        
        // 生成 FASTA 文件
        let fa_path = Path::new(&args.sample_output_dir).join(format!("{}.fasta", output_name));
        let mut fa_writer = create_writer(&fa_path, false)?;
        
        for (id, tag) in all_fa_entries {
            writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(&tag))
                .context("Failed to write FASTA record")?;
        }

        // 生成 .sylsp 文件
        let sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
        let sylsp_file = File::create(&sylsp_path)
            .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
        let sylsp_writer = BufWriter::new(sylsp_file);
        
        bincode::serialize_into(sylsp_writer, &all_sylsp_entries)
            .context("Failed to serialize sylsp data")?;
    }

    // 处理基因组列表文件
    if let Some(genome_list) = &args.genome_list {
        let file = File::open(genome_list)
            .context(format!("Failed to open genome list file: {}", genome_list))?;
        let reader = BufReader::new(file);
        let genome_files: Vec<String> = reader.lines()
            .filter_map(|line| line.ok())
            .collect();

        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_syldb_entries = Vec::new();
        
        // 并行处理所有 FASTA 文件，添加内存监控
        let results: Vec<Result<Vec<SyldbEntry>>> = genome_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = Path::new(file);
                let output_base = Path::new(&args.sample_output_dir).join(input_path.file_stem().unwrap_or_default());
                process_fasta_to_syldb(
                    input_path,
                    &output_base,
                    &enzyme,
                    &args.format,
                    file.ends_with(".gz"),
                )
            })
            .collect();

        // 收集所有结果
        for (file, result) in genome_files.iter().zip(results) {
            match result {
                Ok(mut entries) => {
                    // 为每个条目添加基因组来源信息
                    for entry in &mut entries {
                        entry.genome_source = file.clone();
                    }
                    all_syldb_entries.extend(entries);
                },
                Err(e) => {
                    eprintln!("Error processing FASTA file: {}", e);
                }
            }
        }

        // 生成合并的 .syldb 文件
        if !all_syldb_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_syldb_path = Path::new(&args.sample_output_dir).join(format!("{}.syldb", output_name));
            let combined_syldb_file = File::create(&combined_syldb_path)
                .context(format!("Failed to create combined syldb file: {}", combined_syldb_path.display()))?;
            let combined_syldb_writer = BufWriter::new(combined_syldb_file);
            
            bincode::serialize_into(combined_syldb_writer, &all_syldb_entries)
                .context("Failed to serialize combined syldb data")?;
        }
    }

    // 处理基因组文件
    if let Some(genome_files) = &args.genomes {
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        let mut all_syldb_entries = Vec::new();
        
        // 并行处理所有 FASTA 文件，添加内存监控
        let results: Vec<Result<Vec<SyldbEntry>>> = genome_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = Path::new(file);
                let output_base = Path::new(&args.sample_output_dir).join(input_path.file_stem().unwrap_or_default());
                process_fasta_to_syldb(
                    input_path,
                    &output_base,
                    &enzyme,
                    &args.format,
                    file.ends_with(".gz"),
                )
            })
            .collect();

        // 收集所有结果
        for (file, result) in genome_files.iter().zip(results) {
            match result {
                Ok(mut entries) => {
                    // 为每个条目添加基因组来源信息
                    for entry in &mut entries {
                        entry.genome_source = file.clone();
                    }
                    all_syldb_entries.extend(entries);
                },
                Err(e) => {
                    eprintln!("Error processing FASTA file: {}", e);
                }
            }
        }

        // 生成合并的 .syldb 文件
        if !all_syldb_entries.is_empty() {
            let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
            let combined_syldb_path = Path::new(&args.sample_output_dir).join(format!("{}.syldb", output_name));
            let combined_syldb_file = File::create(&combined_syldb_path)
                .context(format!("Failed to create combined syldb file: {}", combined_syldb_path.display()))?;
            let combined_syldb_writer = BufWriter::new(combined_syldb_file);
            
            bincode::serialize_into(combined_syldb_writer, &all_syldb_entries)
                .context("Failed to serialize combined syldb data")?;
        }
    }

    // 处理样本列表文件
    if let Some(sample_list) = &args.sample_list {
        let mut all_sylsp_entries = Vec::new();
        let enzyme = EnzymeSpec::new(&args.enzyme)?;
        
        // 读取样本列表文件
        let file = File::open(sample_list)
            .context(format!("Failed to open sample list file: {}", sample_list))?;
        let reader = BufReader::new(file);
        
        // 并行处理所有样本文件
        let sample_files: Vec<String> = reader.lines()
            .filter_map(|line| line.ok())
            .collect();
            
        // 使用FxHashMap优化样本处理
        let sample_stats = Arc::new(Mutex::new(SampleStatsMap::default()));
        
        let results: Vec<Result<(String, Vec<(String, TagHash)>, Vec<SylspEntry>)>> = sample_files.par_iter()
            .map(|file| {
                // 检查内存使用
                if let Some(current_memory) = get_memory_usage() {
                    if current_memory > max_ram as f64 {
                        check_vram_and_block(max_ram, file);
                    }
                }
                
                let input_path = PathBuf::from(file);
                // 修正样本名提取逻辑
                let file_name = input_path.file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                let file_stem = file_name.split('.').next().unwrap_or("unknown").to_string();
                
                let reader = fastq::Reader::new(create_reader(&input_path)?);
                let mut fa_entries = Vec::new();
                let mut sylsp_entries = Vec::new();
                let mut stats = ExtractionStats::new();
                // 使用FxHashMap优化标签统计
                let mut tag_frequency = TagFrequencyMap::default();

                for result in reader.records() {
                    let record = result.context("Failed to read FASTQ record")?;
                    stats.total_sequences += 1;
                    stats.total_sequence_length += record.seq().len();
                    
                    let tags = extract_and_validate_tags(record.seq(), &enzyme)
                        .context(format!("Failed to process read: {}", record.id()))?;
                        
                    for (i, tag) in tags.iter().enumerate() {
                        let id = format!("{}_tag{}", record.id(), i + 1);
                        fa_entries.push((id.clone(), tag.clone()));
                        
                        // 统计标签频率
                        *tag_frequency.entry(tag.clone()).or_insert(0) += 1;
                        
                        let entry = SylspEntry {
                            sequence_id: id,
                            tag: hash_bytes(tag),
                            quality: Some(String::from_utf8_lossy(record.qual()).to_string()),
                            sample_source: file_stem.clone(), // 用文件名去除扩展名作为样本名
                        };
                        sylsp_entries.push(entry);
                    }
                    
                    stats.total_tags += tags.len();
                }
                
                // 更新全局统计
                let mut global_stats = sample_stats.lock().unwrap();
                global_stats.insert(file_stem.clone(), stats.clone());
                
                log_stats(stats, &enzyme);
                Ok((file_stem, fa_entries, sylsp_entries))
            })
            .collect();
            
        // 处理每个样本的结果
        for result in results {
            match result {
                Ok((_file_stem, _fa_entries, sylsp_entries)) => {
                    // 注释掉为每个样本生成独立文件的代码
                    // let fa_path = Path::new(&args.sample_output_dir)
                    //     .join(format!("{}.fasta", file_stem));
                    // let mut fa_writer = create_writer(&fa_path, false)?;
                    // 
                    // for (id, tag) in fa_entries {
                    //     writeln!(fa_writer, ">{}\n{}", id, String::from_utf8_lossy(&tag))
                    //         .context("Failed to write FASTA record")?;
                    // }
                    // 
                    // // 为每个样本生成独立的 sylsp 文件
                    // let sample_sylsp_path = Path::new(&args.sample_output_dir)
                    //     .join(format!("{}.sylsp", file_stem));
                    // let sample_sylsp_file = File::create(&sample_sylsp_path)
                    //     .context(format!("Failed to create sylsp file: {}", sample_sylsp_path.display()))?;
                    // let sample_sylsp_writer = BufWriter::new(sample_sylsp_file);
                    // 
                    // bincode::serialize_into(sample_sylsp_writer, &sylsp_entries)
                    //     .context(format!("Failed to serialize sylsp data for sample: {}", file_stem))?;
                    
                    // 收集所有 sylsp 条目用于合并
                    all_sylsp_entries.extend(sylsp_entries);
                },
                Err(e) => eprintln!("Error processing sample file: {}", e),
            }
        }
        
        // 生成合并的 .sylsp 文件
        let output_name = args.out_name.as_ref().map_or_else(|| "combined".to_string(), |s| s.clone());
        let sylsp_path = Path::new(&args.sample_output_dir).join(format!("{}.sylsp", output_name));
        let sylsp_file = File::create(&sylsp_path)
            .context(format!("Failed to create combined sylsp file: {}", sylsp_path.display()))?;
        let sylsp_writer = BufWriter::new(sylsp_file);
        
        bincode::serialize_into(sylsp_writer, &all_sylsp_entries)
            .context("Failed to serialize combined sylsp data")?;
    }

    Ok(())
}

fn process_fasta_to_syldb(
    input: &Path,
    _output_base: &Path,
    enzyme: &EnzymeSpec,
    _format: &str,
    _compress: bool,
) -> Result<Vec<SyldbEntry>> {
    // 注释掉生成单个.fa文件的代码
    // let fa_path = output_base.with_extension("fa");
    // let mut fa_writer = BufWriter::with_capacity(64 * 1024, File::create(&fa_path)?);
    
    let mut stats = ExtractionStats::new();
    // 预分配容量 - 估计每个序列平均产生50个标签
    let mut syldb_entries = Vec::with_capacity(100);
    // 使用FxHashMap优化标签去重和统计
    let mut tag_frequency = TagFrequencyMap::default();

    // 读取和处理 FASTA 记录
    let reader = create_reader(input)?;
    for record in fasta::Reader::new(reader).records() {
        let record = record.context("Failed to read FASTA record")?;
        let seq_len = record.seq().len();
        stats.total_sequences += 1;
        stats.total_sequence_length += seq_len;
        
        // 使用包含canonical处理的标签提取
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process sequence: {}", record.id()))?;
            
        // 预分配位置向量容量
        let mut positions = Vec::with_capacity(tags.len());
        for (i, tag) in tags.iter().enumerate() {
            // 注释掉单个FASTA文件写入
            // writeln!(fa_writer, ">{}_{}\n{}", 
            //     record.id(), 
            //     i + 1,
            //     String::from_utf8_lossy(tag))
            //     .context("Failed to write FASTA record")?;
            
            positions.push(i);
            
            // 统计标签频率（现在使用canonical tags）
            *tag_frequency.entry(tag.clone()).or_insert(0) += 1;
        }

        // 创建 syldb 条目 - 直接使用hash_bytes（现在处理canonical tags）
        let entry = SyldbEntry {
            sequence_id: record.id().to_string(),
            tags: tags.iter().map(|t| hash_bytes(t)).collect(),
            positions,
            genome_source: input.to_string_lossy().to_string(),
            tag_uniqueness: None, // 初始时未标记，将由mark命令处理
        };
        syldb_entries.push(entry);
            
        stats.total_tags += tags.len();
    }

    // 注释掉生成单个.syldb文件的代码
    // let syldb_path = output_base.with_extension("syldb");
    // let syldb_file = File::create(&syldb_path)
    //     .context(format!("Failed to create syldb file: {}", syldb_path.display()))?;
    // let syldb_writer = BufWriter::with_capacity(64 * 1024, syldb_file);
    // 
    // // 使用标准序列化 API
    // bincode::serialize_into(syldb_writer, &syldb_entries)
    //     .context("Failed to serialize syldb data")?;

    log_stats(stats, enzyme);
    Ok(syldb_entries)
}



fn process_paired_fastq_to_sylsp(
    input1: &str,
    input2: &str,
    enzyme: &EnzymeSpec,
    sample_source: &str,
) -> Result<Vec<(String, TagHash, String)>> {
    let reader1 = fastq::Reader::new(create_reader(Path::new(input1))?);
    let reader2 = fastq::Reader::new(create_reader(Path::new(input2))?);
    let mut stats = ExtractionStats::new();
    let mut fa_entries = Vec::new();
    // 使用FxHashSet优化双端测序的去重
    let mut seen_pairs = FxHashSet::default();

    let mut iter1 = reader1.records();
    let mut iter2 = reader2.records();

    loop {
        let record1 = match iter1.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading first pair: {}", e)),
            None => break,
        };

        let record2 = match iter2.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading second pair: {}", e)),
            None => break,
        };

        let seq_len1 = record1.seq().len();
        let seq_len2 = record2.seq().len();
        stats.total_sequences += 1;
        stats.total_sequence_length += seq_len1 + seq_len2;
        
        // 处理第一条序列（使用canonical处理）
        let tags1 = extract_and_validate_tags(record1.seq(), enzyme)
            .context(format!("Failed to process read: {}", record1.id()))?;
            
        // 处理第二条序列（使用canonical处理）
        let tags2 = extract_and_validate_tags(record2.seq(), enzyme)
            .context(format!("Failed to process read: {}", record2.id()))?;
            
        stats.total_tags += tags1.len() + tags2.len();
            
        // 收集 FASTA 条目，使用FxHashSet去重（现在使用canonical tags）
        for (i, tag) in tags1.iter().enumerate() {
            let entry_key = (record1.id().to_string(), i, tag.clone());
            if seen_pairs.insert(entry_key.clone()) {
                fa_entries.push((format!("{}_{}", record1.id(), i + 1), tag.clone(), sample_source.to_string()));
            }
        }

        for (i, tag) in tags2.iter().enumerate() {
            let entry_key = (record2.id().to_string(), i, tag.clone());
            if seen_pairs.insert(entry_key.clone()) {
                fa_entries.push((format!("{}_{}", record2.id(), i + 1), tag.clone(), sample_source.to_string()));
            }
        }
    }

    log_stats(stats, enzyme);
    Ok(fa_entries)
}

fn calculate_tag_percentage(tag_count: usize, total_kmers: usize) -> f64 {
    if total_kmers == 0 {
        0.0
    } else {
        (tag_count as f64 / total_kmers as f64) * 100.0
    }
}



fn read_file_list(path: &str) -> Result<Vec<String>> {
    let file = File::open(path)
        .context(format!("Failed to open file list: {}", path))?;
    let reader = BufReader::new(file);
    let files: Vec<String> = reader.lines()
        .filter_map(|line| line.ok())
        .collect();
    Ok(files)
}

// 优化的标签去重函数，使用FxHashSet提高性能





// 添加反向互补序列计算函数
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        let complement = match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b, // 保持其他字符不变
        };
        rc.push(complement);
    }
    rc
}

// 获取 canonical 版本的序列（字典序较小的）
fn get_canonical_sequence(seq: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(seq);
    
    // 比较正向和反向互补序列的字典序
    if seq <= rc.as_slice() {
        seq.to_vec()
    } else {
        rc
    }
}
