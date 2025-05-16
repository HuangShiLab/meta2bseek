// src/extract.rs
use crate::cmdline::ExtractArgs;
use anyhow::Result;

pub fn extract(args: ExtractArgs) -> Result<()> {
    // 打印基础参数
    println!("Extract command parameters:");
    println!("- Input files: {:?}", args.files);
    println!("- DB output name: {}", args.db_out_name);
    println!("- Sample output dir: {}", args.sample_output_dir);
    println!("- Threads: {}", args.threads);
    println!("- Restriction enzyme: {}", args.enzyme); 
    println!("- Subsampling rate (c): {}", args.c);

    // 可选参数打印
    if let Some(reads) = &args.reads {
        println!("- Reads files: {:?}", reads);
    }
    if let Some(genomes) = &args.genomes {
        println!("- Genome files: {:?}", genomes);
    }
    if args.individual {
        println!("- Individual records mode");
    }
    if args.no_dedup {
        println!("- Disabled deduplication");
    }

    Ok(())
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
        
        let tags = extract_tags(record.seq(), enzyme)
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
        
        let tags = extract_tags(record.seq(), enzyme)
            .context(format!("Failed to process read: {}", record.id()))?;
            
        write_tags(&mut writer, record.id(), &tags, format)
            .context("Failed to write tags")?;
            
        stats.total_tags += tags.len();
    }

    log_stats(stats, enzyme);
    Ok(())
}

fn extract_tags(seq: &[u8], enzyme: &EnzymeSpec) -> Result<Vec<Vec<u8>>> {
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