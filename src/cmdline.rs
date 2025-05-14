use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "speci2bseek")]
#[command(version, about, long_about = r#"
Species-resolved microbiome analysis using 2bRAD technology

EXAMPLES:
  Extract 2bRAD tags from FASTA:
  speci2bseek extract -i genome.fasta -o tags.fa -e BcgI
  
  Process FASTQ data:
  speci2bseek extract -i reads.fastq.gz -o processed_tags.fq -e CjeI -t 8
  
  Show full help:
  speci2bseek extract --help
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
}