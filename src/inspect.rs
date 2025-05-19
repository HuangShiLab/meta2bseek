// src/inspect.rs

use crate::cmdline::InspectArgs;
use anyhow::Result;

pub fn inspect(args: InspectArgs) -> Result<()> {
    println!("Inspect command parameters:");
    
    // 文件列表
    println!("- Files to inspect ({}):", args.files.len());
    for (i, file) in args.files.iter().enumerate() {
        println!("  {}. {}", i + 1, file);
    }
    
    // 输出路径
    match &args.out_file_name {
        Some(path) => println!("- Output file: {}", path),
        None => println!("- Output to stdout"),
    }

    Ok(())
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