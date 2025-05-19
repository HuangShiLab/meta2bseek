use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
//use rayon::prelude::*;
use regex::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};
use glob::glob;
use crate::cmdline::ExtractArgs;
use serde::{Serialize, Deserialize};
// use std::collections::HashSet;

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

#[derive(Serialize, Deserialize)]
pub struct SyldbEntry {
    pub sequence_id: String,
    pub tags: Vec<String>,
    pub positions: Vec<usize>,
}

#[derive(Serialize, Deserialize)]
pub struct SylspEntry {
    pub sequence_id: String,
    pub tag: String,
    pub quality: Option<String>,
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

// Helper function to expand glob patterns
pub fn expand_input_files(patterns: &Vec<String>) -> Result<Vec<PathBuf>> {
    let mut all_paths = Vec::new();
    for pattern in patterns {
        let paths: Result<Vec<_>, _> = glob(pattern)?.collect();
        let mut pattern_paths = paths.context("Failed to read glob pattern")?;
        all_paths.append(&mut pattern_paths);
    }
    Ok(all_paths)
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

#[allow(dead_code)]
fn process_fasta(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let reader = create_reader(input)?;
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    let mut syldb_entries = Vec::new();

    for record in fasta::Reader::new(reader).records() {
        let record = record.context("Failed to read FASTA record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process sequence: {}", record.id()))?;
            
        // 记录 tags 和它们的位置信息
        let mut positions = Vec::new();
        for (i, tag) in tags.iter().enumerate() {
            write_tags(&mut writer, record.id(), &[tag.clone()], format)
                .context("Failed to write tags")?;
            positions.push(i);
        }

        // 创建 syldb 条目
        let entry = SyldbEntry {
            sequence_id: record.id().to_string(),
            tags: tags.iter().map(|t| String::from_utf8_lossy(t).to_string()).collect(),
            positions,
        };
        syldb_entries.push(entry);
            
        stats.total_tags += tags.len();
    }

    // 生成 .syldb 文件
    let syldb_path = output.with_extension("syldb");
    let syldb_file = File::create(&syldb_path)
        .context(format!("Failed to create syldb file: {}", syldb_path.display()))?;
    let syldb_writer = BufWriter::new(syldb_file);
    
    bincode::serialize_into(syldb_writer, &syldb_entries)
        .context("Failed to serialize syldb data")?;

    log_stats(stats, enzyme);
    Ok(())
}

#[allow(dead_code)]
fn process_fastq(
    input: &Path,
    output: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    let reader = fastq::Reader::new(create_reader(input)?);
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();
    let mut syldb_entries = Vec::new();

    for result in reader.records() {
        let record = result.context("Failed to read FASTQ record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process read: {}", record.id()))?;
            
        // 记录 tags 和它们的位置信息
        let mut positions = Vec::new();
        for (i, tag) in tags.iter().enumerate() {
            write_tags(&mut writer, record.id(), &[tag.clone()], format)
                .context("Failed to write tags")?;
            positions.push(i);
        }

        // 创建 syldb 条目
        let entry = SyldbEntry {
            sequence_id: record.id().to_string(),
            tags: tags.iter().map(|t| String::from_utf8_lossy(t).to_string()).collect(),
            positions,
        };
        syldb_entries.push(entry);
            
        stats.total_tags += tags.len();
    }

    // 生成 .syldb 文件
    let syldb_path = output.with_extension("syldb");
    let syldb_file = File::create(&syldb_path)
        .context(format!("Failed to create syldb file: {}", syldb_path.display()))?;
    let syldb_writer = BufWriter::new(syldb_file);
    
    bincode::serialize_into(syldb_writer, &syldb_entries)
        .context("Failed to serialize syldb data")?;

    log_stats(stats, enzyme);
    Ok(())
}

/*
fn extract_and_validate_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Result<Vec<Vec<u8>>> {
    let seq_str = String::from_utf8_lossy(seq);
    let mut tags = Vec::new();

    for pattern in &enzyme.patterns {
        for captures in pattern.captures_iter(&seq_str) {
            if let Some(m) = captures.get(0) {
                // Directly use the matched pattern as the tag
                tags.push(m.as_str().as_bytes().to_vec());
            }
        }
    }

    Ok(tags)
}
*/

fn extract_and_validate_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Result<Vec<Vec<u8>>> {
    let seq_str = String::from_utf8_lossy(seq);
    let mut tags = Vec::new();

    for pattern in &enzyme.patterns {
        for m in pattern.find_iter(&seq_str) {
            let matched = m.as_str().as_bytes().to_vec();
            tags.push(matched);
        }
    }

    Ok(tags)
}


fn write_tags(
    writer: &mut dyn Write,
    seq_id: &str,
    tags: &[Vec<u8>],
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

    Ok(if path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

fn create_writer(path: &Path, compress: bool) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .context(format!("Failed to create output file: {}", path.display()))?;

    Ok(if compress {
        Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
    } else {
        Box::new(BufWriter::new(file))
    })
}

struct ExtractionStats {
    total_sequences: usize,
    total_tags: usize,
}

impl ExtractionStats {
    fn new() -> Self {
        Self {
            total_sequences: 0,
            total_tags: 0,
        }
    }
}

fn log_stats(stats: ExtractionStats, enzyme: &EnzymeSpec) {
    println!(
        "\nProcessing complete for {}:\n\
        =============================\n\
        - Total sequences processed: {}\n\
        - Total tags extracted: {}\n\
        - Average tags per sequence: {:.2}\n\
        - Recognition patterns used: {}",
        enzyme.name,
        stats.total_sequences,
        stats.total_tags,
        stats.total_tags as f32 / stats.total_sequences.max(1) as f32,
        enzyme.patterns
            .iter()
            .map(|r| r.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    );
}

pub fn extract(args: ExtractArgs) -> Result<()> {
    // 创建输出目录
    std::fs::create_dir_all(&args.sample_output_dir)
        .context("Failed to create output directory")?;

    // 处理基因组文件（使用 glob 模式）
    if let Some(genome_pattern) = &args.genomes {
        let genome_files = expand_input_files(genome_pattern)
            .context("Failed to expand genome file pattern")?;
        
        for file in genome_files {
            let output_base = Path::new(&args.sample_output_dir).join(file.file_stem().unwrap_or_default());
            
            // 处理 FASTA 文件，生成 .fa 和 .syldb
            process_fasta_to_syldb(
                &file,
                &output_base,
                &EnzymeSpec::new(&args.enzyme)?,
                &args.format,
                file.to_string_lossy().ends_with(".gz"),
            )?;
        }
    }

    // 处理读取文件，生成 .sylsp
    if let Some(read_files) = args.reads {
        for file in read_files {
            let input_path = PathBuf::from(&file);
            let output_base = Path::new(&args.sample_output_dir).join(input_path.file_stem().unwrap_or_default());
            
            // 处理 FASTQ 文件，只生成 .sylsp
            process_fastq_to_sylsp(
                &input_path,
                &output_base,
                &EnzymeSpec::new(&args.enzyme)?,
            )?;
        }
    }

    Ok(())
}

fn process_fasta_to_syldb(
    input: &Path,
    output_base: &Path,
    enzyme: &EnzymeSpec,
    format: &str,
    compress: bool,
) -> Result<()> {
    // 创建 .fa 文件的写入器
    let fa_path = output_base.with_extension(format);
    let mut fa_writer = create_writer(&fa_path, compress)?;
    
    let mut stats = ExtractionStats::new();
    let mut syldb_entries = Vec::new();

    // 读取和处理 FASTA 记录
    let reader = create_reader(input)?;
    for record in fasta::Reader::new(reader).records() {
        let record = record.context("Failed to read FASTA record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process sequence: {}", record.id()))?;
            
        // 记录 tags 和它们的位置信息
        let mut positions = Vec::new();
        for (i, tag) in tags.iter().enumerate() {
            // 写入到 .fa 文件
            write_tags(&mut fa_writer, record.id(), &[tag.clone()], format)
                .context("Failed to write tags")?;
            positions.push(i);
        }

        // 创建 syldb 条目
        let entry = SyldbEntry {
            sequence_id: record.id().to_string(),
            tags: tags.iter().map(|t| String::from_utf8_lossy(t).to_string()).collect(),
            positions,
        };
        syldb_entries.push(entry);
            
        stats.total_tags += tags.len();
    }

    // 生成 .syldb 文件
    let syldb_path = output_base.with_extension("syldb");
    let syldb_file = File::create(&syldb_path)
        .context(format!("Failed to create syldb file: {}", syldb_path.display()))?;
    let syldb_writer = BufWriter::new(syldb_file);
    
    bincode::serialize_into(syldb_writer, &syldb_entries)
        .context("Failed to serialize syldb data")?;

    log_stats(stats, enzyme);
    Ok(())
}

fn process_fastq_to_sylsp(
    input: &Path,
    output_base: &Path,
    enzyme: &EnzymeSpec,
) -> Result<()> {
    let reader = fastq::Reader::new(create_reader(input)?);
    let mut stats = ExtractionStats::new();
    let mut sylsp_entries = Vec::new();

    for result in reader.records() {
        let record = result.context("Failed to read FASTQ record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process read: {}", record.id()))?;
            
        stats.total_tags += tags.len();
            
        // 为每个标签创建一个 sylsp 条目
        for tag in &tags {
            let entry = SylspEntry {
                sequence_id: record.id().to_string(),
                tag: String::from_utf8_lossy(tag).to_string(),
                quality: Some(String::from_utf8_lossy(record.qual()).to_string()),
            };
            sylsp_entries.push(entry);
        }
    }

    // 生成 .sylsp 文件
    let sylsp_path = output_base.with_extension("sylsp");
    let sylsp_file = File::create(&sylsp_path)
        .context(format!("Failed to create sylsp file: {}", sylsp_path.display()))?;
    let sylsp_writer = BufWriter::new(sylsp_file);
    
    bincode::serialize_into(sylsp_writer, &sylsp_entries)
        .context("Failed to serialize sylsp data")?;

    log_stats(stats, enzyme);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bcgi_pattern_comparison() {
        let enzyme = EnzymeSpec::new("BcgI").unwrap();
        let test_seq = b"CCCGGATGGTAGCGTCAACGTACGTCGCTTCGGCGGCATGAAAATCGAGCGCACCTGGTTCGCCGCCGATAAGACCGGCT";
        
        // Test captures_iter version
        let mut tags1 = Vec::new();
        let seq_str = String::from_utf8_lossy(test_seq);
        for pattern in &enzyme.patterns {
            for captures in pattern.captures_iter(&seq_str) {
                if let Some(m) = captures.get(0) {
                    tags1.push(m.as_str().as_bytes().to_vec());
                }
            }
        }
        
        // Test find_iter version
        let mut tags2 = Vec::new();
        for pattern in &enzyme.patterns {
            for m in pattern.find_iter(&seq_str) {
                tags2.push(m.as_str().as_bytes().to_vec());
            }
        }
        
        println!("Number of tags found by captures_iter: {}", tags1.len());
        println!("Number of tags found by find_iter: {}", tags2.len());
        
        // Print the actual tags found
        println!("\nTags found by captures_iter:");
        for (i, tag) in tags1.iter().enumerate() {
            println!("Tag {}: {}", i + 1, String::from_utf8_lossy(tag));
        }
        
        println!("\nTags found by find_iter:");
        for (i, tag) in tags2.iter().enumerate() {
            println!("Tag {}: {}", i + 1, String::from_utf8_lossy(tag));
        }
        
        // Assert that both methods find the same number of tags
        assert_eq!(tags1.len(), tags2.len());
    }
}

