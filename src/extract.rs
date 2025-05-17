// extract.rs
use anyhow::{Context, Result};
use bincode;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;
use regex::Regex;
use bio::io::{fasta, fastq};
use std::{
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};
use crate::cmdline::ExtractArgs;


// 酶切位点定义
const ENZYME_DEFINITIONS: &[(&str, &[&str])] = &[
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
struct EnzymeSpec {
    name: String,
    patterns: Vec<Regex>,
    tag_lengths: HashSet<usize>,
}

impl EnzymeSpec {
    fn new(name: &str) -> Result<Self> {
        let def = ENZYME_DEFINITIONS
            .iter()
            .find(|(e, _)| *e == name)
            .ok_or_else(|| anyhow::anyhow!("Unsupported enzyme: {}", name))?;

        let patterns = def.1
            .iter()
            .map(|p| Regex::new(p))
            .collect::<Result<Vec<_>, _>>()
            .context("Failed to compile regex patterns")?;

        // 计算预期标签长度
        let mut tag_lengths = HashSet::new();
        for p in def.1 {
            if let Some(caps) = Regex::new(r"\^([ACGT]+)\$")?.captures(p) {
                if let Some(m) = caps.get(1) {
                    tag_lengths.insert(m.as_str().len());
                }
            }
        }

        Ok(Self {
            name: def.0.to_string(),
            patterns,
            tag_lengths,
        })
    }

    fn validate_tag(&self, tag: &[u8]) -> bool {
        self.tag_lengths.contains(&tag.len())
    }
}

// 基因组数据库结构
#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct GenomeDatabase {
    enzyme: String,
    records: Vec<GenomeRecord>,
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct GenomeRecord {
    file_name: String,
    contigs: Vec<ContigTags>,
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]  
struct ContigTags {
    name: String,
    tags: Vec<Vec<u8>>,
    positions: Vec<usize>,
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct SampleProfile {
    enzyme: String,
    file_name: String,
    total_tags: usize,
    unique_tags: HashMap<Vec<u8>, u32>,
    mean_read_length: f64,
}

pub fn extract(args: ExtractArgs) -> Result<()> {
    let enzyme = EnzymeSpec::new(&args.enzyme)?;

    // 创建输出目录
    fs::create_dir_all(&args.sample_output_dir)?;

    // 并行处理所有输入
    rayon::scope(|s| {
        if let Some(genomes) = &args.genomes {
            s.spawn(|_| {
                let _ = process_genomes(genomes, &args.db_out_name, &enzyme);
            });
        }

        if let Some(reads) = &args.reads {
            reads.par_iter().for_each(|file| {
                let _ = process_single_sample(file, &args.sample_output_dir, &enzyme);
            });
        }

        if !args.first_pair.is_empty() {
            args.first_pair.par_iter()
                .zip(&args.second_pair)
                .for_each(|(f1, f2)| {
                    let _ = process_paired_sample(f1, f2, &args.sample_output_dir, &enzyme);
                });
        }
    });

    Ok(())
}

fn process_genomes(files: &[String], db_name: &str, enzyme: &EnzymeSpec) -> Result<()> {
    let records: Vec<GenomeRecord> = files
        .par_iter()
        .map(|file| {
            let reader = create_reader(Path::new(file))?;
            let fasta_reader = fasta::Reader::new(reader);
            let mut contigs = Vec::new();

            // 一次性获取所有记录
            for record in fasta_reader.records() {
                let record = record?;
                let (tags, positions) = extract_and_validate_tags(record.seq(), enzyme)?;
                contigs.push(ContigTags {
                    name: String::from_utf8_lossy(record.id().as_bytes()).to_string(),
                    tags,
                    positions,
                });
            }

            Ok(GenomeRecord {
                file_name: file.clone(),
                contigs,
            })
        })
        .collect::<Result<_>>()?;

    let db = GenomeDatabase {
        enzyme: enzyme.name.clone(),
        records,
    };

    let db_path = PathBuf::from(db_name).with_extension("syldb");
    let file = BufWriter::new(File::create(&db_path)?);
    bincode::serialize_into(file, &db)?;

    Ok(()) // 明确返回Result
}

fn process_single_sample(file: &str, output_dir: &str, enzyme: &EnzymeSpec) -> Result<()> {
    // 创建文件读取器
    let reader = create_reader(Path::new(file))?;
    let fastq_reader = fastq::Reader::new(reader);
    
    // 初始化统计信息
    let mut profile = SampleProfile {
        enzyme: enzyme.name.clone(),
        file_name: file.to_string(),
        total_tags: 0,
        unique_tags: HashMap::new(),
        mean_read_length: 0.0,
    };

    // 初始化统计计数器
    let mut total_length = 0u64;
    let mut read_count = 0u64;

    // 遍历所有fastq记录
    for result in fastq_reader.records() {
        let record = result?;
        let seq = record.seq();
        
        // 统计读长
        total_length += seq.len() as u64;
        read_count += 1;

        // 提取验证标签
        let (tags, _positions) = extract_and_validate_tags(seq, enzyme)?;
        
        // 更新标签统计
        profile.total_tags += tags.len();
        for tag in tags {
            *profile.unique_tags.entry(tag).or_insert(0) += 1;
        }
    }

    // 计算平均读长（至少处理1条read时）
    if read_count > 0 {
        profile.mean_read_length = total_length as f64 / read_count as f64;
    }

    // 写入结果文件
    write_profile(output_dir, profile)
}

fn process_paired_sample(f1: &str, f2: &str, output_dir: &str, enzyme: &EnzymeSpec) -> Result<()> {
    // 创建双端文件读取器
    let r1_reader = fastq::Reader::new(create_reader(Path::new(f1))?);
    let r2_reader = fastq::Reader::new(create_reader(Path::new(f2))?);

    // 初始化统计信息
    let mut profile = SampleProfile {
        enzyme: enzyme.name.clone(),
        file_name: format!("{}_&_{}", f1, f2),
        total_tags: 0,
        unique_tags: HashMap::new(),
        mean_read_length: 0.0,
    };

    // 初始化统计计数器
    let mut total_length = 0u64;
    let mut read_count = 0u64;

    // 并行迭代双端reads
    for (r1_result, r2_result) in r1_reader.records().zip(r2_reader.records()) {
        // 处理可能的读取错误
        let r1 = r1_result?;
        let r2 = r2_result?;

        // 处理R1
        let seq1 = r1.seq();
        total_length += seq1.len() as u64;
        
        // 处理R2
        let seq2 = r2.seq();
        total_length += seq2.len() as u64;
        
        read_count += 2; // 每个循环处理一对reads

        // 提取双端tags
        let (tags1, _) = extract_and_validate_tags(seq1, enzyme)?;
        let (tags2, _) = extract_and_validate_tags(seq2, enzyme)?;

        // 合并统计tags
        profile.total_tags += tags1.len() + tags2.len();
        for tag in tags1.into_iter().chain(tags2) {
            *profile.unique_tags.entry(tag).or_insert(0) += 1;
        }
    }

    // 计算平均读长（至少处理1对reads时）
    if read_count > 0 {
        profile.mean_read_length = total_length as f64 / read_count as f64;
    }

    // 写入结果文件
    write_profile(output_dir, profile)
}

fn extract_and_validate_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Result<(Vec<Vec<u8>>, Vec<usize>)> {
    let seq_str = String::from_utf8_lossy(seq);
    let mut tags = Vec::new();
    let mut positions = Vec::new();

    for pattern in &enzyme.patterns {
        for m in pattern.find_iter(&seq_str) {
            let start = m.start();
            let end = m.end();
            let tag = &seq[start..end];

            if enzyme.validate_tag(tag) {
                tags.push(tag.to_vec());
                positions.push(start);
            }
        }
    }

    Ok((tags, positions))
}

fn create_reader(path: &Path) -> Result<Box<dyn Read + Send>, std::io::Error> {
    let file = File::open(path)?;
    if path.to_string_lossy().ends_with(".gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn create_writer(path: &Path, compress: bool) -> Result<Box<dyn Write>, std::io::Error> {
    let file = File::create(path)?;
    if compress {
        Ok(Box::new(GzEncoder::new(file, Compression::default())))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_profile(output_dir: &str, profile: SampleProfile) -> Result<()> {
    let output_path = PathBuf::from(output_dir)
        .join(Path::new(&profile.file_name).file_name().unwrap())
        .with_extension("sylsp");

    let file = BufWriter::new(File::create(&output_path)?);
    bincode::serialize_into(file, &profile)?;

    Ok(())
}

// 判断是否为fasta文件
fn is_fasta_file(path: &Path) -> Result<bool> {
    let ext = path.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase();

    if matches!(ext.as_str(), "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn") {
        return Ok(true);
    }

    // 内容验证
    let mut file = File::open(path)?;
    let mut buffer = [0; 1];
    file.read_exact(&mut buffer)?;
    Ok(buffer[0] == b'>')
}

fn process_read(file_path: &str, args: &ExtractArgs) -> Result<()> {
    let enzyme = EnzymeSpec::new(&args.enzyme)
        .context(format!("Failed to initialize enzyme {}", args.enzyme))?;

    // 安全处理文件名
    let input_path = Path::new(file_path);
    let file_name = input_path.file_name()
        .ok_or_else(|| anyhow!("Invalid file path: {}", file_path))?;
    
    // 构建输出路径
    let output_path = Path::new(&args.sample_output_dir)
        .join(file_name)
        .with_extension("sylsp");

    // 创建带缓冲的写入器
    let mut writer = create_writer(&output_path, true)
        .context(format!("Failed to create output file: {}", output_path.display()))?;

    // 创建带错误处理的读取器
    let reader = create_reader(input_path)
        .context(format!("Failed to open input file: {}", file_path))?;
    let fastq_reader = fastq::Reader::new(reader);

    // 初始化统计计数器
    let mut total_reads = 0u64;
    let mut total_tags = 0u64;

    // 处理每条记录
    for result in fastq_reader.records() {
        let record = result
            .context(format!("Failed to read record #{}", total_reads + 1))?;
        
        total_reads += 1;

        // 清理可能包含非法字符的ID
        let clean_id = record.id().replace(|c: char| c.is_whitespace(), "_");
        
        // 提取并验证标签
        let (tags, _) = extract_and_validate_tags(record.seq(), &enzyme)
            .context(format!("Failed to process read {}: {}", clean_id, record.seq().escape_ascii()))?;

        total_tags += tags.len() as u64;

        // 写入所有标签
        for (i, tag) in tags.into_iter().enumerate() {
            write_fastq_tag(&mut writer, &clean_id, i, &tag)
                .context(format!("Failed to write tag for read {}", clean_id))?;
        }
    }

    // 记录处理统计
    log::info!(
        "Processed {} reads from {}, extracted {} tags",
        total_reads,
        file_path,
        total_tags
    );

    Ok(())
}

// 写入FASTA格式的tag
fn write_fasta_tag(writer: &mut dyn Write, seq_id: &str, tag: &[u8]) -> std::io::Result<()> {
    writeln!(writer, ">{}_2brad_tag\n{}", seq_id, String::from_utf8_lossy(tag))?;
    Ok(())
}

// 写入FASTQ格式的tag
fn write_fastq_tag(writer: &mut dyn Write, seq_id: &str, index: usize, tag: &[u8]) -> std::io::Result<()> {
    let qual = "~".repeat(tag.len());
    writeln!(writer, "@{}_tag{}\n{}\n+\n{}", seq_id, index+1, 
        String::from_utf8_lossy(tag), qual)?;
    Ok(())
}

// 日志统计
fn log_processing_stats(stats: &ProcessingStats, enzyme: &EnzymeSpec) {
    println!(
        "\nProcessed {} sequences with {}:\n\
        - Total tags extracted: {}\n\
        - Average tags per sequence: {:.2}",
        stats.total_sequences,
        enzyme.name,
        stats.total_tags,
        stats.total_tags as f32 / stats.total_sequences.max(1) as f32
    );
}

// 处理列表输入
fn process_list_inputs(args: &ExtractArgs) -> Result<()> {
    if let Some(list_path) = &args.list_sequence {
        let file = File::open(list_path)?;
        let reader = BufReader::new(file);
        
        for line in reader.lines() {
            let path = line?;
            if is_fasta_file(Path::new(&path))? {
                process_genome(&path, args)?;
            } else {
                process_read(&path, args)?;
            }
        }
    }
    Ok(())
}

// 处理列表中的双端文件
fn process_paired_list(list1: &str, list2: &str, args: &ExtractArgs) -> Result<()> {
    let file1 = File::open(list1)?;
    let file2 = File::open(list2)?;
    let list1 = BufReader::new(file1);
    let list2 = BufReader::new(file2);

    for (line1, line2) in list1.lines().zip(list2.lines()) {
        let f1 = line1?;
        let f2 = line2?;
        process_paired_end(&f1, &f2, args)?;
    }

    Ok(())
}

// 处理单个双端样本
// Add this at the top of extract.rs
use anyhow::anyhow;

// Fix the extract_and_validate_tags usage in process_paired_end
fn process_paired_end(f1: &str, f2: &str, args: &ExtractArgs) -> Result<()> {
    let enzyme = EnzymeSpec::new(&args.enzyme)?;
    let output_dir = Path::new(&args.sample_output_dir);
    let base_name = format!("{}_paired", Path::new(f1).file_stem().unwrap().to_str().unwrap());
    let output_path = output_dir.join(format!("{}.sylsp", base_name));

    let mut writer: Box<dyn Write + 'static> = create_writer(&output_path, true)?;
    let r1_reader = fastq::Reader::new(create_reader(Path::new(f1))?);
    let r2_reader = fastq::Reader::new(create_reader(Path::new(f2))?);

    for (rec1, rec2) in r1_reader.records().zip(r2_reader.records()) {
        let rec1 = rec1?;
        let rec2 = rec2?;

        let (tags1, _) = extract_and_validate_tags(rec1.seq(), &enzyme)?;
        let (tags2, _) = extract_and_validate_tags(rec2.seq(), &enzyme)?;

        for (i, tag) in tags1.into_iter().chain(tags2).enumerate() {
            write_fastq_tag(&mut writer, rec1.id(), i, &tag)?;
        }
    }

    Ok(())
}

// Fix the fasta reader loop in process_genome
fn process_genome(path: &str, args: &ExtractArgs) -> Result<()> {
    let enzyme = EnzymeSpec::new(&args.enzyme)?;
    let reader = bio::io::fasta::Reader::new(create_reader(Path::new(path))?);
    
    let mut tags_collection = Vec::new();
    
    // Collect records first to avoid move issues
    let records: Vec<_> = reader.records().collect::<Result<_, _>>()?;
    
    for record in records {
        let (tags, _) = extract_and_validate_tags(record.seq(), &enzyme)?;
        tags_collection.extend(tags);
    }
    
    let output_path = Path::new(&args.db_out_name).with_extension("syldb");
    let writer = create_writer(&output_path, false)?;
    bincode::serialize_into(writer, &tags_collection)?;
    
    Ok(())
}

// 统计数据结构
struct ProcessingStats {
    total_sequences: usize,
    total_tags: usize,
}

/* 
use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
//use rayon::prelude::*;
use regex::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
};

const ENZYME_DEFINITIONS: &[(&str, &[&str])] = &[
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
struct EnzymeSpec {
    name: String,
    patterns: Vec<Regex>,
}

impl EnzymeSpec {
    fn new(name: &str) -> Result<Self> {
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

pub fn process_input(
    input_path: &Path,
    output_path: &Path,
    enzyme_name: &str,
    threads: usize,
    format: &str,
    compress: bool,
) -> Result<()> {
    let enzyme = EnzymeSpec::new(enzyme_name)
        .context(format!("Unsupported enzyme: {}", enzyme_name))?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Failed to create thread pool")?;

    let is_fasta = is_fasta_file(input_path)
        .context("Failed to determine input file type")?;

    pool.install(|| {
        if is_fasta {
            process_fasta(input_path, output_path, &enzyme, format, compress)
        } else {
            process_fastq(input_path, output_path, &enzyme, format, compress)
        }
    })
}

fn is_fasta_file(path: &Path) -> Result<bool> {
    // Check common FASTA extensions
    let ext = path.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase();

    let is_fasta_ext = matches!(ext.as_str(), 
        "fa" | "fasta" | "fna" | "ffn" | "faa" | "frn"
    );

    // Content-based validation for ambiguous cases
    if !is_fasta_ext {
        let mut file = File::open(path)
            .context(format!("Failed to open file for content check: {}", path.display()))?;
        let mut buffer = [0; 1];
        file.read_exact(&mut buffer)
            .context("Failed to read first byte for content check")?;
        return Ok(buffer[0] == b'>');
    }

    Ok(true)
}

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

    for record in fasta::Reader::new(reader).records() {
        let record = record.context("Failed to read FASTA record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process sequence: {}", record.id()))?;
            
        write_tags(&mut writer, record.id(), &tags, format)
            .context("Failed to write tags")?;
            
        stats.total_tags += tags.len();
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
    let reader = fastq::Reader::new(create_reader(input)?);
    let mut writer = create_writer(output, compress)?;
    let mut stats = ExtractionStats::new();

    for result in reader.records() {
        let record = result.context("Failed to read FASTQ record")?;
        stats.total_sequences += 1;
        
        let tags = extract_and_validate_tags(record.seq(), enzyme)
            .context(format!("Failed to process read: {}", record.id()))?;
            
        write_tags(&mut writer, record.id(), &tags, format)
            .context("Failed to write tags")?;
            
        stats.total_tags += tags.len();
    }

    log_stats(stats, enzyme);
    Ok(())
}

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
*/