// src/inspect.rs

use crate::cmdline::InspectArgs;
use anyhow::{Context, Result};
use bincode;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use regex::Regex;

#[derive(Serialize, Deserialize, Debug)]
struct InspectResult {
    file_type: String,
    file_name: String,
    total_tags: usize,
    tag_lengths: Vec<usize>,
    tag_length_distribution: Vec<(usize, usize, f64)>, // (length, count, percentage)
    enzyme: String,
    patterns: Vec<String>,
}

pub fn inspect(args: InspectArgs) -> Result<()> {
    let mut writer = match args.out_file_name {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)) as Box<dyn Write>,
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write>,
    };

    for file in &args.files {
        match inspect_file(file) {
            Ok(result) => {
                writeln!(writer, "File Information:")?;
                writeln!(writer, "----------------")?;
                writeln!(writer, "File: {}", result.file_name)?;
                writeln!(writer, "Type: {}", result.file_type)?;
                writeln!(writer, "Total tags: {}", result.total_tags)?;
                
                writeln!(writer, "\nDetected Enzyme Information:")?;
                writeln!(writer, "-------------------------")?;
                writeln!(writer, "Enzyme: {}", result.enzyme)?;
                writeln!(writer, "\nRecognition patterns:")?;
                for (i, pattern) in result.patterns.iter().enumerate() {
                    writeln!(writer, "  {}. {}", i + 1, pattern)?;
                }
                
                writeln!(writer, "\nTag Length Distribution:")?;
                writeln!(writer, "----------------------")?;
                writeln!(writer, "{:<10} {:<10}", "Length", "Count")?;
                writeln!(writer, "{:-<20}", "")?;
                
                for (length, count, _) in result.tag_length_distribution {
                    writeln!(writer, "{:<10} {:<10}", length, count)?;
                }
                writeln!(writer, "\n")?;
            }
            Err(e) => {
                eprintln!("Failed to inspect {}: {}", file, e);
            }
        }
    }

    Ok(())
}

fn inspect_file(file_path: &str) -> Result<InspectResult> {
    let path = Path::new(file_path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    match path.extension().and_then(|s| s.to_str()) {
        Some("syldb") => inspect_syldb(reader, file_path),
        Some("sylsp") => inspect_sylsp(reader, file_path),
        _ => Err(anyhow::anyhow!("Unknown file extension, expected .syldb or .sylsp")),
    }
}

fn inspect_syldb(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let entries: Vec<crate::extract::SyldbEntry> = bincode::deserialize_from(reader)
        .context("Failed to deserialize .syldb file")?;

    let mut tag_lengths = Vec::new();
    let mut tags = Vec::new();
    for entry in &entries {
        for tag in &entry.tags {
            tag_lengths.push(tag.len());
            tags.push(tag.clone());
        }
    }

    let distribution = calculate_tag_distribution(&tag_lengths);

    // 更智能的酶识别
    let (enzyme, patterns, matched_count, matched_ratio) = determine_enzyme_by_regex(&tags);

    println!("Total tags: {}", tags.len());
    println!("Best matched enzyme: {}", enzyme);
    println!("Recognition patterns:");
    for (i, pattern) in patterns.iter().enumerate() {
        println!("  {}. {}", i + 1, pattern);
    }
    println!("Matched tags: {} ({:.2}%)", matched_count, matched_ratio * 100.0);

    Ok(InspectResult {
        file_type: "GenomeDatabase".to_string(),
        file_name: file_path.to_string(),
        total_tags: tag_lengths.len(),
        tag_lengths,
        tag_length_distribution: distribution,
        enzyme,
        patterns,
    })
}

fn inspect_sylsp(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let entries: Vec<crate::extract::SylspEntry> = bincode::deserialize_from(reader)
        .context("Failed to deserialize .sylsp file")?;

    let mut tag_lengths = Vec::new();
    let mut tags = Vec::new();
    for entry in &entries {
        tag_lengths.push(entry.tag.len());
        tags.push(entry.tag.clone());
    }

    let distribution = calculate_tag_distribution(&tag_lengths);

    // 更智能的酶识别
    let (enzyme, patterns, matched_count, matched_ratio) = determine_enzyme_by_regex(&tags);

    println!("Total tags: {}", tags.len());
    println!("Best matched enzyme: {}", enzyme);
    println!("Recognition patterns:");
    for (i, pattern) in patterns.iter().enumerate() {
        println!("  {}. {}", i + 1, pattern);
    }
    println!("Matched tags: {} ({:.2}%)", matched_count, matched_ratio * 100.0);

    Ok(InspectResult {
        file_type: "SampleProfile".to_string(),
        file_name: file_path.to_string(),
        total_tags: tag_lengths.len(),
        tag_lengths,
        tag_length_distribution: distribution,
        enzyme,
        patterns,
    })
}

fn determine_enzyme_by_regex(tags: &[String]) -> (String, Vec<String>, usize, f64) {
    if tags.is_empty() {
        return ("unknown".to_string(), Vec::new(), 0, 0.0);
    }
    let mut best_enzyme = "unknown".to_string();
    let mut best_patterns = Vec::new();
    let mut best_count = 0;
    let mut best_ratio = 0.0;
    for (name, patterns) in crate::extract::ENZYME_DEFINITIONS {
        let regexes: Vec<Regex> = patterns.iter().filter_map(|p| Regex::new(p).ok()).collect();
        let mut count = 0;
        for tag in tags {
            if regexes.iter().any(|re| re.is_match(tag)) {
                count += 1;
            }
        }
        let ratio = count as f64 / tags.len() as f64;
        if count > best_count || (count == best_count && ratio > best_ratio) {
            best_enzyme = name.to_string();
            best_patterns = patterns.iter().map(|p| p.to_string()).collect();
            best_count = count;
            best_ratio = ratio;
        }
    }
    (best_enzyme, best_patterns, best_count, best_ratio)
}

fn calculate_tag_distribution(tag_lengths: &[usize]) -> Vec<(usize, usize, f64)> {
    use std::collections::HashMap;
    
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for &length in tag_lengths {
        *length_counts.entry(length).or_insert(0) += 1;
    }

    let total = tag_lengths.len() as f64;
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

/*
// src/inspect.rs
use crate::cmdline::InspectArgs;
use anyhow::{Context, Result};
use bincode;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Serialize, Deserialize, Debug)]
struct InspectResult {
    file_type: String,
    file_name: String,
    enzyme: String,
    num_records: usize,
    total_tags: usize,
    mean_read_length: Option<f64>,
    first_contig_name: Option<String>,
}

pub fn inspect(args: InspectArgs) -> Result<()> {
    let mut writer = match args.out_file_name {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)) as Box<dyn Write>,
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write>,
    };

    for file in &args.files {
        match inspect_file(file) {
            Ok(result) => {
                let yaml = serde_yaml::to_string(&result)?;
                writeln!(writer, "{}", yaml)?;
            }
            Err(e) => {
                warn!("Failed to inspect {}: {}", file, e);
            }
        }
    }

    Ok(())
}

fn inspect_file(file_path: &str) -> Result<InspectResult> {
    let path = Path::new(file_path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    if path.extension().and_then(|s| s.to_str()) == Some("syldb") {
        inspect_syldb(reader, file_path)
    } else if path.extension().and_then(|s| s.to_str()) == Some("sylsp") {
        inspect_sylsp(reader, file_path)
    } else {
        Err(anyhow::anyhow!("Unknown file extension, expected .syldb or .sylsp"))
    }
}

fn inspect_syldb(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let db: crate::extract::GenomeDatabase = bincode::deserialize_from(reader)
        .context("Failed to deserialize .syldb file")?;

    let first_contig_name = db.records.first()
        .and_then(|r| r.contigs.first())
        .map(|c| c.name.clone());

    let total_tags = db.records.iter()
        .flat_map(|r| r.contigs.iter())
        .map(|c| c.tags.len())
        .sum();

    Ok(InspectResult {
        file_type: "GenomeDatabase".to_string(),
        file_name: file_path.to_string(),
        enzyme: db.enzyme,
        num_records: db.records.len(),
        total_tags,
        mean_read_length: None,
        first_contig_name,
    })
}

fn inspect_sylsp(reader: BufReader<File>, file_path: &str) -> Result<InspectResult> {
    let profile: crate::extract::SampleProfile = bincode::deserialize_from(reader)
        .context("Failed to deserialize .sylsp file")?;

    Ok(InspectResult {
        file_type: "SampleProfile".to_string(),
        file_name: profile.file_name.clone(),
        enzyme: profile.enzyme,
        num_records: 1,
        total_tags: profile.total_tags,
        mean_read_length: Some(profile.mean_read_length),
        first_contig_name: None,
    })
}
*/