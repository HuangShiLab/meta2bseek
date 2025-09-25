// src/view.rs

use crate::cmdline::ViewArgs;
use crate::sketch::SequencesSketch;
use crate::extract::GenomeSketch;
use anyhow::{Context, Result};
use bincode;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::collections::HashMap;
use crate::constants::Hash;

#[derive(Serialize, Deserialize, Debug)]
struct ViewResult {
    file_type: String,
    file_name: String,
    c: usize,
    k: usize,
    num_records: usize,
    total_kmers: usize,
    unique_kmers: usize,
    kmer_frequency_stats: Vec<(Hash, u32)>,
    mean_read_length: Option<f64>,
    first_contig_name: Option<String>,
    genome_sources: Option<Vec<String>>,
    sample_sources: Option<Vec<SampleStats>>,
    // 每个样本对应的 kmer -> 计数
    per_sample_kmer_counts: Option<std::collections::HashMap<String, std::collections::HashMap<Hash, u32>>>,
    kmer_lengths: Vec<usize>,
    kmer_length_distribution: Vec<(usize, usize, f64)>,
    min_spacing: Option<usize>,
    genome_stats: Option<Vec<GenomeStats>>,
}

#[derive(Serialize, Deserialize, Debug)]
struct GenomeStats {
    source: String,
    num_records: usize,
    total_kmers: usize,
    kmer_length_distribution: Vec<(usize, usize, f64)>,
}

#[derive(Serialize, Deserialize, Debug)]
struct SampleStats {
    source: String,
    num_records: usize,
    total_kmers: usize,
    kmer_length_distribution: Vec<(usize, usize, f64)>,
}

#[derive(Debug)]
struct KmerMatrix {
    samples: Vec<String>,
    kmers: Vec<Hash>,
    matrix: HashMap<(String, Hash), u32>,
}

pub fn view(args: ViewArgs) -> Result<()> {
    let mut writer = match args.out_file_name {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)) as Box<dyn Write>,
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write>,
    };

    // 收集所有样本的k-mer信息用于生成矩阵
    let mut kmer_matrix = KmerMatrix {
        samples: Vec::new(),
        kmers: Vec::new(),
        matrix: HashMap::new(),
    };

    for file in &args.files {
        match view_file(file) {
            Ok(result) => {
                // 输出文件信息
                writeln!(writer, "File Information:")?;
                writeln!(writer, "----------------")?;
                writeln!(writer, "File: {}", result.file_name)?;
                writeln!(writer, "Type: {}", result.file_type)?;
                writeln!(writer, "C-value: {}", result.c)?;
                writeln!(writer, "K-size: {}", result.k)?;
                
                // 样本有无与数量（若存在多个sample则按sample显示数量）
                match &result.sample_sources {
                    Some(samples) if !samples.is_empty() => {
                        writeln!(writer, "Samples: yes ({} samples)", samples.len())?;
                        writeln!(writer, "Per-sample total k-mers:")?;
                        for s in samples {
                            writeln!(writer, "  {}: {}", s.source, s.total_kmers)?;
                        }
                    }
                    _ => {
                        writeln!(writer, "Samples: no")?;
                    }
                }
                writeln!(writer, "Total k-mers: {}", result.total_kmers)?;
                writeln!(writer, "Unique k-mers: {}", result.unique_kmers)?;
                
                if let Some(mean_length) = result.mean_read_length {
                    writeln!(writer, "Mean read length: {:.2}", mean_length)?;
                }
                
                writeln!(writer, "\nK-mer Frequency Statistics:")?;
                writeln!(writer, "---------------------------")?;
                writeln!(writer, "{:<20} {:<10}", "K-mer Hash", "Count")?;
                writeln!(writer, "{:-<30}", "")?;
                
                // 显示前20个最常见的k-mer
                let display_count = std::cmp::min(20, result.kmer_frequency_stats.len());
                for (kmer_hash, count) in result.kmer_frequency_stats.iter().take(display_count) {
                    writeln!(writer, "{:<20} {:<10}", format!("{:016x}", kmer_hash), count)?;
                }
                if result.kmer_frequency_stats.len() > display_count {
                    writeln!(writer, "... and {} more unique k-mers", result.kmer_frequency_stats.len() - display_count)?;
                }
                
                writeln!(writer, "\nK-mer Length Distribution:")?;
                writeln!(writer, "-------------------------")?;
                writeln!(writer, "{:<10} {:<10}", "Length", "Count")?;
                writeln!(writer, "{:-<20}", "")?;
                
                for (length, count, _) in &result.kmer_length_distribution {
                    writeln!(writer, "{:<10} {:<10}", length, count)?;
                }

                if let Some(sample_stats) = &result.sample_sources {
                    writeln!(writer, "\nSample-specific statistics:")?;
                    writeln!(writer, "------------------------")?;
                    for sample in sample_stats {
                        writeln!(writer, "\nSample: {}", sample.source)?;
                        writeln!(writer, "  Records: {}", sample.num_records)?;
                        writeln!(writer, "  Total k-mers: {}", sample.total_kmers)?;
                        writeln!(writer, "  K-mer length distribution:")?;
                        for (length, count, percentage) in &sample.kmer_length_distribution {
                            writeln!(writer, "    Length {}: {} k-mers ({:.2}%)", length, count, percentage)?;
                        }
                    }
                }

                if let Some(stats) = &result.genome_stats {
                    writeln!(writer, "\nGenome-specific statistics:")?;
                    writeln!(writer, "------------------------")?;
                    for genome in stats {
                        writeln!(writer, "\nGenome: {}", genome.source)?;
                        writeln!(writer, "  Records: {}", genome.num_records)?;
                        writeln!(writer, "  Total k-mers: {}", genome.total_kmers)?;
                        writeln!(writer, "  K-mer length distribution:")?;
                        for (length, count, _) in &genome.kmer_length_distribution {
                            writeln!(writer, "    Length {}: {} k-mers", length, count)?;
                        }
                    }
                }
                writeln!(writer, "\n")?;

                // 收集k-mer矩阵信息
                collect_kmer_matrix_data(&result, &mut kmer_matrix);
            }
            Err(e) => {
                eprintln!("Failed to view {}: {}", file, e);
            }
        }
    }

    // 如果指定了输出路径，生成TSV矩阵
    if let Some(log_path) = &args.log_path {
        generate_tsv_matrix(&kmer_matrix, log_path, &args.tsv_name)?;
    }

    Ok(())
}

fn view_file(file_path: &str) -> Result<ViewResult> {
    let path = Path::new(file_path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match path.extension().and_then(|s| s.to_str()) {
        Some("syldb") => view_syldb(reader, file_path),
        Some("sylsp") => view_sylsp(reader, file_path),
        _ => Err(anyhow::anyhow!("Unknown file extension, expected .syldb or .sylsp")),
    }
}

fn view_syldb(reader: BufReader<File>, file_path: &str) -> Result<ViewResult> {
    println!("Attempting to deserialize {} as genome sketches...", file_path);
    
    let entries: Vec<GenomeSketch> = bincode::deserialize_from(reader)
        .with_context(|| format!("Failed to deserialize .syldb file: {}", file_path))?;

    if entries.is_empty() {
        return Err(anyhow::anyhow!("Empty .syldb file"));
    }

    println!("Successfully deserialized {} genome sketches from {}", entries.len(), file_path);

    let mut kmer_lengths = Vec::new();
    let mut genome_sources = std::collections::HashSet::new();
    let mut genome_stats = std::collections::HashMap::new();
    let mut kmer_frequency = std::collections::HashMap::new();

    // 获取第一个条目的参数作为参考
    let first_entry = &entries[0];
    let c = first_entry.c;
    let k = first_entry.k;
    let min_spacing = first_entry.min_spacing;

    for (i, entry) in entries.iter().enumerate() {
        if i % 100 == 0 {
            println!("Processing genome sketch {}/{}", i + 1, entries.len());
        }
        
        let source = &entry.file_name;
        genome_sources.insert(source.clone());
        
        for kmer in &entry.genome_kmers {
            kmer_lengths.push(k); // K-mer size is always k
            *kmer_frequency.entry(*kmer).or_insert(0) += 1;
        }

        let stats = genome_stats.entry(source.clone()).or_insert(GenomeStats {
            source: source.clone(),
            num_records: 1, // Each entry represents one genome
            total_kmers: 0,
            kmer_length_distribution: Vec::new(),
        });
        
        stats.total_kmers += entry.genome_kmers.len();
    }

    for stats in genome_stats.values_mut() {
        let mut lengths = Vec::new();
        for entry in &entries {
            if entry.file_name == stats.source {
                for _kmer in &entry.genome_kmers {
                    lengths.push(k); // K-mer size is always k
                }
            }
        }
        stats.kmer_length_distribution = calculate_kmer_distribution(&lengths);
    }

    let distribution = calculate_kmer_distribution(&kmer_lengths);

    // 计算k-mer统计信息
    let unique_kmers = kmer_frequency.len();
    let mut kmer_frequency_stats: Vec<(Hash, u32)> = kmer_frequency.into_iter().collect();
    kmer_frequency_stats.sort_by(|a, b| b.1.cmp(&a.1)); // 按频率降序排序

    Ok(ViewResult {
        file_type: "GenomeSketch".to_string(),
        file_name: file_path.to_string(),
        c,
        k,
        num_records: entries.len(),
        total_kmers: kmer_lengths.len(),
        unique_kmers,
        kmer_frequency_stats,
        mean_read_length: None,
        first_contig_name: entries.first().map(|e| e.first_contig_name.clone()),
        genome_sources: if genome_sources.is_empty() {
            None
        } else {
            Some(genome_sources.into_iter().collect())
        },
        sample_sources: None,
        per_sample_kmer_counts: None,
        kmer_lengths,
        kmer_length_distribution: distribution,
        min_spacing: Some(min_spacing),
        genome_stats: Some(genome_stats.into_values().collect()),
    })
}

fn view_sylsp(reader: BufReader<File>, file_path: &str) -> Result<ViewResult> {
    // 尝试反序列化为单个SequencesSketch
    let single_sketch: Result<SequencesSketch, _> = bincode::deserialize_from(reader);
    
    if let Ok(sketch) = single_sketch {
        // 单个样本文件
        return view_single_sylsp(sketch, file_path);
    }
    
    // 如果单个反序列化失败，尝试作为多个sketch的列表
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    
    // 添加错误处理和日志
    println!("Attempting to deserialize {} as multiple sketches...", file_path);
    
    let sketches: Result<Vec<SequencesSketch>, _> = bincode::deserialize_from(reader);
    
    if let Ok(sketches) = sketches {
        if sketches.is_empty() {
            return Err(anyhow::anyhow!("Empty .sylsp file"));
        }

        println!("Successfully deserialized {} sketches from {}", sketches.len(), file_path);

        let mut kmer_lengths = Vec::new();
        let mut sample_stats = std::collections::HashMap::new();
        let mut kmer_frequency = std::collections::HashMap::new();
        let mut per_sample_kmer_counts: std::collections::HashMap<String, std::collections::HashMap<Hash, u32>> = std::collections::HashMap::new();

        // 获取第一个条目的参数作为参考
        let first_sketch = &sketches[0];
        let c = first_sketch.c;
        let k = first_sketch.k;

        for (i, sketch) in sketches.iter().enumerate() {
            if i % 100 == 0 {
                println!("Processing sketch {}/{}", i + 1, sketches.len());
            }
            
            let sample_name = sketch.sample_name.as_ref().unwrap_or(&sketch.file_name);
            
            for (kmer, count) in &sketch.kmer_counts {
                for _ in 0..*count {
                    kmer_lengths.push(k); // K-mer size is always k
                }
                *kmer_frequency.entry(*kmer).or_insert(0) += count;

                // 累积每个样本的 k-mer 计数
                let sample_entry = per_sample_kmer_counts
                    .entry(sample_name.clone())
                    .or_insert_with(std::collections::HashMap::new);
                *sample_entry.entry(*kmer).or_insert(0) += count;
            }

            let stats = sample_stats.entry(sample_name.clone()).or_insert(SampleStats {
                source: sample_name.clone(),
                num_records: 1, // Each sketch represents one sample
                total_kmers: 0,
                kmer_length_distribution: Vec::new(),
            });
            stats.total_kmers += sketch.kmer_counts.values().sum::<u32>() as usize;
        }

        for stats in sample_stats.values_mut() {
            let mut sample_lengths = Vec::new();
            for sketch in &sketches {
                let sample_name = sketch.sample_name.as_ref().unwrap_or(&sketch.file_name);
                if sample_name == &stats.source {
                    for (_, count) in &sketch.kmer_counts {
                        for _ in 0..*count {
                            sample_lengths.push(k); // K-mer size is always k
                        }
                    }
                }
            }
            stats.kmer_length_distribution = calculate_kmer_distribution(&sample_lengths);
        }

        let distribution = calculate_kmer_distribution(&kmer_lengths);

        // 计算k-mer统计信息
        let unique_kmers = kmer_frequency.len();
        let mut kmer_frequency_stats: Vec<(Hash, u32)> = kmer_frequency.into_iter().collect();
        kmer_frequency_stats.sort_by(|a, b| b.1.cmp(&a.1)); // 按频率降序排序

        // 计算平均read长度
        let total_mean_length: f64 = sketches.iter()
            .map(|s| s.mean_read_length)
            .sum::<f64>() / sketches.len() as f64;

        return Ok(ViewResult {
            file_type: "SampleSketch".to_string(),
            file_name: file_path.to_string(),
            c,
            k,
            num_records: sketches.len(),
            total_kmers: kmer_lengths.len(),
            unique_kmers,
            kmer_frequency_stats,
            mean_read_length: Some(total_mean_length),
            first_contig_name: None,
            genome_sources: None,
            sample_sources: Some(sample_stats.into_values().collect()),
            per_sample_kmer_counts: Some(per_sample_kmer_counts),
            kmer_lengths,
            kmer_length_distribution: distribution,
            min_spacing: None,
            genome_stats: None,
        });
    }
    
    // 如果Meta2bseek格式失败，尝试sylph格式
    println!("Meta2bseek format failed, attempting sylph format...");
    return Err(anyhow::anyhow!("File format not recognized. This file may be in sylph format or corrupted."));
}

fn view_single_sylsp(sketch: SequencesSketch, file_path: &str) -> Result<ViewResult> {
    let mut kmer_lengths = Vec::new();
    let mut kmer_frequency = std::collections::HashMap::new();
    let mut per_sample_kmer_counts: std::collections::HashMap<String, std::collections::HashMap<Hash, u32>> = std::collections::HashMap::new();

    let sample_name = sketch.sample_name.as_ref().unwrap_or(&sketch.file_name);
    
    for (kmer, count) in &sketch.kmer_counts {
        for _ in 0..*count {
            kmer_lengths.push(sketch.k); // K-mer size is always k
        }
        *kmer_frequency.entry(*kmer).or_insert(0) += count;

        // 累积样本的 k-mer 计数
        let sample_entry = per_sample_kmer_counts
            .entry(sample_name.clone())
            .or_insert_with(std::collections::HashMap::new);
        *sample_entry.entry(*kmer).or_insert(0) += count;
    }

    let distribution = calculate_kmer_distribution(&kmer_lengths);

    // 计算k-mer统计信息
    let unique_kmers = kmer_frequency.len();
    let mut kmer_frequency_stats: Vec<(Hash, u32)> = kmer_frequency.into_iter().collect();
    kmer_frequency_stats.sort_by(|a, b| b.1.cmp(&a.1)); // 按频率降序排序

    let sample_stats = vec![SampleStats {
        source: sample_name.clone(),
        num_records: 1,
        total_kmers: kmer_lengths.len(),
        kmer_length_distribution: distribution.clone(),
    }];

    Ok(ViewResult {
        file_type: "SampleSketch".to_string(),
        file_name: file_path.to_string(),
        c: sketch.c,
        k: sketch.k,
        num_records: 1,
        total_kmers: kmer_lengths.len(),
        unique_kmers,
        kmer_frequency_stats,
        mean_read_length: Some(sketch.mean_read_length),
        first_contig_name: None,
        genome_sources: None,
        sample_sources: Some(sample_stats),
        per_sample_kmer_counts: Some(per_sample_kmer_counts),
        kmer_lengths,
        kmer_length_distribution: distribution,
        min_spacing: None,
        genome_stats: None,
    })
}

fn calculate_kmer_distribution(kmer_lengths: &[usize]) -> Vec<(usize, usize, f64)> {
    use std::collections::HashMap;
    
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for &length in kmer_lengths {
        *length_counts.entry(length).or_insert(0) += 1;
    }

    let total = kmer_lengths.len() as f64;
    let mut distribution: Vec<_> = length_counts
        .into_iter()
        .map(|(length, count)| {
            let percentage = (count as f64 / total) * 100.0;
            (length, count, percentage)
        })
        .collect();

    distribution.sort_by_key(|&(length, _, _)| length);
    distribution
}

fn collect_kmer_matrix_data(result: &ViewResult, kmer_matrix: &mut KmerMatrix) {
    match result.file_type.as_str() {
        "SampleSketch" => {
            if let Some(per_sample) = &result.per_sample_kmer_counts {
                for (sample_name, kmer_counts) in per_sample {
                    if !kmer_matrix.samples.contains(sample_name) {
                        kmer_matrix.samples.push(sample_name.clone());
                    }
                    for (kmer_hash, count) in kmer_counts {
                        if !kmer_matrix.kmers.contains(kmer_hash) {
                            kmer_matrix.kmers.push(*kmer_hash);
                        }
                        let key = (sample_name.clone(), *kmer_hash);
                        let entry = kmer_matrix.matrix.entry(key).or_insert(0);
                        *entry += *count;
                    }
                }
            }
        }
        "GenomeSketch" => {
            // 对于基因组sketch文件，使用文件名作为样本名
            let sample_name = Path::new(&result.file_name)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            
            if !kmer_matrix.samples.contains(&sample_name) {
                kmer_matrix.samples.push(sample_name.clone());
            }
            
            // 收集该数据库中的所有k-mer
            for (kmer_hash, count) in &result.kmer_frequency_stats {
                if !kmer_matrix.kmers.contains(kmer_hash) {
                    kmer_matrix.kmers.push(*kmer_hash);
                }
                let key = (sample_name.clone(), *kmer_hash);
                let entry = kmer_matrix.matrix.entry(key).or_insert(0);
                *entry += *count as u32;
            }
        }
        _ => {}
    }
}

fn generate_tsv_matrix(kmer_matrix: &KmerMatrix, log_path: &str, tsv_name: &str) -> Result<()> {
    // 确保输出目录存在
    std::fs::create_dir_all(log_path)?;
    
    // 构建完整的文件路径
    let tsv_path = Path::new(log_path).join(tsv_name);
    let mut tsv_writer = BufWriter::new(File::create(&tsv_path)?);
    
    // 排序样本和k-mer以确保输出的一致性
    let mut sorted_samples = kmer_matrix.samples.clone();
    sorted_samples.sort();
    let mut sorted_kmers = kmer_matrix.kmers.clone();
    sorted_kmers.sort();
    
    // 写入表头
    write!(tsv_writer, "K-mer")?; // 目前保存的是哈希，无法还原原始序列，这里用16进制hash展示
    for sample in &sorted_samples {
        write!(tsv_writer, "\t{}", sample)?;
    }
    writeln!(tsv_writer)?;
    
    // 写入数据行
    for kmer in &sorted_kmers {
        write!(tsv_writer, "{:016x}", kmer)?;
        for sample in &sorted_samples {
            let count = kmer_matrix.matrix.get(&(sample.clone(), *kmer)).unwrap_or(&0);
            write!(tsv_writer, "\t{}", count)?;
        }
        writeln!(tsv_writer)?;
    }
    
    println!("K-mer count matrix saved to: {}", tsv_path.display());
    Ok(())
}
