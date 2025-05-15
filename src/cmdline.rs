use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "Meta2bSeek")]
#[command(version, about, long_about = r#"
Species-resolved microbiome analysis using 2bRAD technology

EXAMPLES:
  Extract 2bRAD tags from FASTA:
  Meta2bSeek extract -i genome.fasta -o tags.fa -e BcgI
  
  Process FASTQ data:
  Meta2bSeek extract -i reads.fastq.gz -o processed_tags.fq -e CjeI -t 8
  
  Show full help:
  Meta2bSeek extract --help
"#)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Extract 2bRAD tags from input sequences
    #[command(arg_required_else_help = true)]
    Extract {
        /// Input file path (FASTA/FASTQ)
        #[arg(short, long)]
        input: PathBuf,
        
        /// Output file path
        #[arg(short, long)]
        output: PathBuf,
        
        /// Restriction enzyme to use
        #[arg(short, long, default_value = "BcgI")]
        enzyme: String,
        
        /// Number of processing threads
        #[arg(short, long, default_value_t = 4)]
        threads: usize,
        
        /// Output format (fa/fq)
        #[arg(long, default_value = "fa")]
        format: String,
        
        /// Enable gzip compression
        #[arg(long, default_value_t = false)]
        compress: bool,
    },
    /// Estimate the containment average nucleotide identity (ANI) of a reference genome to the genomes in your sample.
    #[command(arg_required_else_help = true)]
    Query {
        #[clap(short, long)]
        sample: PathBuf,
        #[clap(short, long)]
        database: PathBuf,
        #[clap(short, long, default_value = "4")]
        threads: usize,
        #[clap(short, long)]
        output: PathBuf,
    },
    /// Profile metagenome sample
    Profile {
        #[clap(short, long)]
        sample: PathBuf,
        #[clap(short, long)]
        database: PathBuf,
        #[clap(short, long, default_value = "4")]
        threads: usize,
        #[clap(short, long)]
        output: PathBuf,
    },
    ///  Inspect sketched .syldb and .sylsp files
    Inspect {
        #[clap(short, long)]
        sample: PathBuf,
        #[clap(short, long)]
        database: PathBuf,
        #[clap(short, long, default_value = "4")]
        threads: usize,
        #[clap(short, long)]
        output: PathBuf,
    }
}