use clap::{Args, Parser, Subcommand};
pub(crate) use crate::constants::*;

#[derive(Parser)]
#[clap(author, version, about = "Ultrafast genome ANI queries and taxonomic profiling for metagenomic shotgun samples.\n\n--- Preparing inputs by extracting (indexing) 2bRAD tags\n## fastq (reads) and fasta (genomes all at once)\n## *.sylsp found in -d; *.syldb given by -o\nmeta2bseek extract -t 5 sample1.fq sample2.fq genome1.fa genome2.fa -o genome1+genome2 -d sample_dir\n\n## paired-end reads\nmeta2bseek extract -1 a_1.fq b_1.fq -2 b_2.fq b_2.fq -d paired_extracts\n\n--- Nearest neighbour containment ANI\nmeta2bseek query *.syldb *.sylsp > all-to-all-query.tsv\n\n--- Taxonomic profiling with relative abundances and ANI\nmeta2bseek profile *.syldb *.sylsp > all-to-all-profile.tsv", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// extract sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb. 
    #[clap(display_order = 1)]
    Extract(ExtractArgs),
    /// Coverage-adjusted ANI querying between databases and samples.
    #[clap(display_order = 3)]
    Query(ContainArgs),
    ///Species-level taxonomic profiling with abundances and ANIs. 
    #[clap(display_order = 2)]
    Profile(ContainArgs),
    ///Inspect extracted .syldb and .sylsp files.
    #[clap(arg_required_else_help = true, display_order = 4)]
    Inspect(InspectArgs),
}


#[derive(Args, Default)]
pub struct ExtractArgs {
    #[clap(num_args=1.., help_heading = "INPUT", help = "fasta/fastq files; gzip optional. Default: fastq file produces a sample extract (*.sylsp) while fasta files are combined into a database (*.syldb).")]
    pub files: Vec<String>,
    #[clap(short='o',long="out-name-db", default_value = "database", help_heading = "OUTPUT", help = "Output name for database extract (with .syldb appended)")]
    pub db_out_name: String,
    #[clap(short='d',long="sample-output-directory", default_value = "./", help_heading = "OUTPUT", help = "Output directory for sample extractes")]
    pub sample_output_dir: String,
    #[clap(short,long="individual-records", help_heading = "GENOME INPUT", help = "Use individual records (contigs) for database construction")]
    pub individual: bool,
    #[clap(num_args=1..,short,long="reads", help_heading = "SINGLE-END INPUT", help = "Single-end fasta/fastq reads")]
    pub reads: Option<Vec<String>>,
    #[clap(num_args=1..,short='g', long="genomes", help_heading = "GENOME INPUT", help = "Genomes in fasta format")]
    pub genomes: Option<Vec<String>>,
    #[clap(short,long="list", help_heading = "INPUT", help = "Newline delimited file with inputs; fastas -> database, fastq -> sample")]
    pub list_sequence: Option<String>,
    #[clap(long="rl", hide=true, help_heading = "SINGLE-END INPUT", help = "Newline delimited file; inputs assumed reads")]
    pub list_reads: Option<String>,
    #[clap(long="gl", help_heading = "GENOME INPUT", help = "Newline delimited file; inputs assumed genomes")]
    pub list_genomes: Option<String>,
    #[clap(long="l1", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are first pair of PE reads")]
    pub list_first_pair: Option<String>,
    #[clap(long="l2", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are second pair of PE reads")]
    pub list_second_pair: Option<String>,
    #[clap(long="lS", help_heading = "INPUT", help = "Newline delimited file; read extractes are renamed to given sample names")]
    pub list_sample_names: Option<String>,
    #[clap(num_args=1.., short='S', long="sample-names", help_heading = "INPUT", help = "Read extractes are renamed to given sample names")]
    pub sample_names: Option<Vec<String>>,

    #[clap(short, long, default_value = "BcgI", help_heading = "ALGORITHM", help = "Restriction enzyme to use. Available options: CspCI, AloI, BsaXI, BaeI, BcgI, CjeI, PpiI, PsrI, BplI, FalI, Bsp24I, HaeIV, CjePI, Hin4I, AlfI, BslFI")]
    pub enzyme: String,
    #[clap(short, long, help_heading = "ALGORITHM", help = "The format of the output file")]
    pub format: String,
    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="ram-barrier", help = "Stop multi-threaded read extracting when (virtual) RAM is past this value (in GB). Does NOT guarantee max RAM limit", hide=true)]
    pub max_ram: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,


    #[clap(long="no-dedup", help_heading = "ALGORITHM", help = "Disable read deduplication procedure. Reduces memory; not recommended for illumina data")]
    pub no_dedup: bool,
    #[clap(long="disable-profiling", help_heading = "ALGORITHM", help = "Disable meta2bseek profile usage for databases; may decrease size and make meta2bseek query slightly faster", hide=true)]
    pub no_pseudotax: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "ALGORITHM", help = "Minimum spacing between selected 2bRAD tags on the genomes")]
    pub min_spacing_kmer: usize,
    #[clap(long="fpr", default_value_t = DEFAULT_FPR, help_heading = "ALGORITHM", help = "False positive rate for read deduplicate hashing; valid values in [0,1).")]
    pub fpr: f64,
    #[clap(short='1',long="first-pairs", num_args=1.., help_heading = "PAIRED-END INPUT", help = "First pairs for paired end reads")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pairs", num_args=1.., help_heading = "PAIRED-END INPUT", help = "Second pairs for paired end reads")]
    pub second_pair: Vec<String>,
}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(num_args=1.., help = "Pre-extracted *.syldb/*.sylsp files. Raw single-end fastq/fasta are allowed and will be automatically extracted to .sylsp/.syldb")]
    pub files: Vec<String>,

    #[clap(short='l',long="list", help = "Newline delimited file of file inputs",help_heading = "INPUT/OUTPUT")]
    pub file_list: Option<String>,

    #[clap(long,default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum 2bRAD tag multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(short='M',long,default_value_t = 50., help_heading = "ALGORITHM", help = "Exclude genomes with less than this number of extracted 2bRAD tags")]
    pub min_number_kmers: f64,
    #[clap(short, long="minimum-ani", help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 90 for query and 95 for profile. Smaller than 95 for profile will give inaccurate results." )]
    pub minimum_ani: Option<f64>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='s', long="sample-threads", help = "Number of samples to be processed concurrently. Default: (# of total threads / 3) + 1 for profile, 1 for query")]
    pub sample_threads: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(short='u', long="estimate-unknown", help_heading = "ALGORITHM", help = "Estimate true coverage and scale sequence abundance in `profile` by estimated unknown sequence percentage" )]
    pub estimate_unknown: bool,

    
    #[clap(short='I',long="read-seq-id", help_heading = "ALGORITHM", help = "Sequence identity (%) of reads. Only used in -u option and overrides automatic detection. ")]
    pub seq_id: Option<f64>,

    //#[clap(short='l', long="read-length", help_heading = "ALGORITHM", help = "Read length (single-end length for pairs). Only necessary for short-read coverages when using --estimate-unknown. Not needed for long-reads" )]
    //pub read_length: Option<usize>,

    #[clap(short='R', long="redundancy-threshold", help_heading = "ALGORITHM", help = "Removes redundant genomes up to a rough ANI percentile when profiling", default_value_t = 99.0, hide=true)]
    pub redundant_ani: f64,

    #[clap(short='r',long="reads", num_args=1.., help = "Single-end raw reads (fastx/gzip)", display_order = 1, help_heading = "extracting")]
    pub reads: Vec<String>,

    #[clap(short='1', long="first-pairs", num_args=1.., help = "First pairs for raw paired-end reads (fastx/gzip)", help_heading = "extracting")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second-pairs", num_args=1.., help = "Second pairs for raw paired-end reads (fastx/gzip)", help_heading = "extracting")]
    pub second_pair: Vec<String>,

    #[clap(short, default_value_t = 200, help_heading = "extracting", help = "Subsampling rate. Does nothing for pre-extracted files")]
    pub c: usize,
    #[clap(short,long="individual-records", help_heading = "extracting", help = "Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-extracted files")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "extracting", help = "Minimum spacing between selected 2bRAD tags on the database genomes. Does nothing for pre-extracted files")]
    pub min_spacing_kmer: usize,

    #[clap(short='o',long="output-file", help = "Output to this file (TSV format). [default: stdout]", help_heading="INPUT/OUTPUT")]
    pub out_file_name: Option<String>,
    #[clap(long="log-reassignments", help = "Output information for how 2bRAD tags for genomes are reassigned during `profile`. Caution: can be verbose and slows down computation.")]
    pub log_reassignments: bool,


    //Hidden options that are embedded in the args but no longer used... 
    #[clap(short, hide=true, long="pseudotax", help_heading = "ALGORITHM", help = "Pseudo taxonomic classification mode. This removes shared 2bRAD tags between species by assigning 2bRAD tags to the highest ANI species. Requires extractes with --enable-pseudotax option" )]
    pub pseudotax: bool,
    #[clap(long="ratio", hide=true)]
    pub ratio: bool,
    #[clap(long="mme", hide=true)]
    pub mme: bool,
    #[clap(long="mle", hide=true)]
    pub mle: bool,
    #[clap(long="nb", hide=true)]
    pub nb: bool,
    #[clap(long="no-ci", help = "Do not output confidence intervals", hide=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hide=true)]
    pub no_adj: bool,
    #[clap(long="mean-coverage", help_heading = "ALGORITHM", help = "Use the robust mean coverage estimator instead of median estimator", hide=true )]
    pub mean_coverage: bool,

}

#[derive(Args)]
pub struct InspectArgs {
    #[clap(num_args=1.., help = "Pre-extracted *.syldb/*.sylsp files.")]
    pub files: Vec<String>,
    #[clap(short='o',long="output-file", help = "Output to this file (YAML format). [default: stdout]")]
    pub out_file_name: Option<String>,

}
    