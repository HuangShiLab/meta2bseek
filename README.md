# meta2bseek

Meta2bSeek is composed of tools in three categories: 
* taxonomic profiling tool; 
* strain and functional profiling tool; 
* statistical analysis tools. 

**Taxonomic profiling tool**: We engineered a Rust-optimized pipeline for ultrafast species identification using k-mers selected by type 2b sites, achieving algorithmic acceleration through parallel computing and memory-efficient indexing. A comprehensive database spanning 394921 bacterial, 7777 archaeal, and 11184 fungal genomes was curated. Cross-kingdom phylogenetic relationships were reconstructed using 2b-site mash distance-based hierarchical clustering. 

**Strain and functional profiling tool**: A novel k-mer hashing strategy was implemented to identify both known and novel strains at a high ANI resolution, enabling to map uncharacterized strains to their closest phylogenetic neighbors in our expanded reference database.

## Installation

Meta2bSeek is a Rust-based application that requires compilation from source.

### Supported Platforms

- **macOS** (10.15+ recommended)
- **Linux Distributions**:
  - Ubuntu (16.04+, including WSL)
  - CentOS (7+)
  - Other Linux distributions with glibc 2.17+
- **Windows Subsystem for Linux (WSL)**

### Prerequisites

#### 1. Install Rust Toolchain

##### MacOS
```bash
# Install Xcode command line tools (if not already installed)
xcode-select --install

# Install Rust using rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

##### Ubuntu/Debian (including WSL)
```bash
# Update package list
sudo apt update

# Install build dependencies
sudo apt install build-essential curl

# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

##### CentOS/RHEL
```bash
# Install development tools
sudo yum groupinstall "Development Tools"
sudo yum install curl

# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

#### 2. Verify Rust Installation
```bash
rustc --version
cargo --version
```
Note: Ensure you have Rust 1.70.0 or higher

### Installation Methods

#### Method 1: Install from Source (Recommended)

Clone the repository:

```bash
git clone https://github.com/HuangShiLab/meta2bseek.git
cd meta2bseek
```

Build in release mode:

```bash
cargo build --release
```

After compilation, the executable binary file meta2bseek will be generated in the "target/release" directory.


## How to run?

### Overview of Meta2bSeek:

```
meta2bseek -h

Ultrafast genome ANI queries and taxonomic profiling for metagenomic shotgun samples.

--- Preparing inputs by extracting (indexing) 2bRAD tags
## fastq (reads) and fasta (genomes all at once)
## *.sylsp found in -d; *.syldb given by -o
meta2bseek extract -t 5 sample1.fq sample2.fq genome1.fa genome2.fa -o genome1+genome2 -d sample_dir

## paired-end reads
meta2bseek extract -1 a_1.fq b_1.fq -2 b_2.fq b_2.fq -d paired_extracts

## batch process reads from a list file
meta2bseek extract -s reads_list.txt -d batch_output --out-name batch_process

--- Nearest neighbour containment ANI
meta2bseek query *.syldb *.sylsp > all-to-all-query.tsv

--- Taxonomic profiling with relative abundances and ANI
meta2bseek profile *.syldb *.sylsp > all-to-all-profile.tsv

Usage: meta2bseek <COMMAND>

Commands:
  extract  extract sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb
  sketch   sketch sequences using k-mer sampling (similar to sylph). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb
  profile  Species-level taxonomic profiling with abundances and ANIs
  query    Coverage-adjusted ANI querying between databases and samples
  inspect  Inspect extracted .syldb and .sylsp files
  view     View sketched .syldb and .sylsp files (Meta2bseek sketch format)
  mark     Mark unique (taxa-specific) tags in .syldb files

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### `extract`: Extract 2bRAD tags from genome file(s) or sample file(s)

**For Reference Genomes**

•	**Input**: A list of paths to reference genome FASTA files (compressed or uncompressed).

•	**Example input file** (genome_list.txt):

```
/lustre1/g/aos_shihuang/databases/GTDBtk_r220/gtdb_genomes_reps_r220/database/GCA/902/513/955/GCA_902513955.1_genomic.fna.gz
/lustre1/g/aos_shihuang/databases/GTDBtk_r220/gtdb_genomes_reps_r220/database/GCA/902/513/415/GCA_902513415.1_genomic.fna.gz
```

•	**Example Command:**

```
meta2bseek extract -t 8 -k genome_list.txt --sample-output-dir /path/to/output --out-name reference_db
# Note: 8 is the number of threads you want to used
```

•	**Output:** A .syldb file (e.g., reference_db.syldb) in the specified output directory.

**For Sample Sequences**

- For single-end sequencing samples:

•	**Input:** A list of paths to sample FASTQ files (compressed or uncompressed).

•	**Example input file** (sample_list.txt):

```
/lustre1/g/aos_shihuang/data/BCB_12hps_monomicrobial_20241001/Enterobacteriaceae_EB/24M9200249_B_12hps.fastq.gz
/lustre1/g/aos_shihuang/data/BCB_12hps_monomicrobial_20241001/Enterobacteriaceae_EB/24M9203369_A_12hps.fastq.gz
```

•	**Example Command:**

```
meta2bseek extract -t 8 -s sample_list.txt --sample-output-dir /path/to/output --out-name samples
```

•	**Output:** A .sylsp file (e.g., samples.sylsp) in the specified output directory.

- For paired-end sequencing samples:

•	**Input:** Two separate lists of paths to sample FASTQ files for the left and right reads (compressed or uncompressed).

•	**Example input file**:

(sample_left_list.txt)

```
/lustre1/g/aos_shihuang/Strain2b/data/MSA1002_1.clean.fq.gz
/lustre1/g/aos_shihuang/Strain2b/data/MSA1003_1.clean.fq.gz
```

(sample_right_list.txt):

```
/lustre1/g/aos_shihuang/Strain2b/data/MSA1002_2.clean.fq.gz
/lustre1/g/aos_shihuang/Strain2b/data/MSA1003_2.clean.fq.gz
```

•	**Example Command:**

```
meta2bseek extract -t 20 --l1 sample_left_list.txt  --l2 sample_right_list.txt --sample-output-dir /path/to/output --out-name samples)
```

•	**Output:** A .sylsp file (e.g., samples.sylsp) in the specified output directory.

**Usages:**
```
meta2bseek extract -h

extract sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb

Usage: meta2bseek extract [OPTIONS] --sample-output-dir <SAMPLE_OUTPUT_DIR>

Options:
  -t, --threads <THREADS>  Number of threads [default: 3]
  -f, --format <FORMAT>    Output format (fa or fq) [default: fa]
      --debug              Debug output
  -h, --help               Print help

GENOME INPUT:
  -g, --genomes <GENOMES>...       One or more genome files in fasta format
  -k, --genome-list <GENOME_LIST>  Text file containing paths to genome files (one per line)

READ INPUT:
  -r, --reads <READS>...           One or more fastq files for reads
  -s, --sample-list <SAMPLE_LIST>  Text file containing paths to fastq files (one per line)

OUTPUT:
  -o, --output <OUTPUT_DIR>                    Output directory for extracted tags [default: .]
  -d, --sample-output-dir <SAMPLE_OUTPUT_DIR>  Output directory for sample files
  -n, --out-name <OUT_NAME>                    Output name for generated files

ALGORITHM:
  -e, --enzyme <ENZYME>  Restriction enzyme to use [default: BcgI]

PAIRED READ INPUT:
  -1, --first-pair <FIRST_PAIR>    First pair of paired-end reads
  -2, --second-pair <SECOND_PAIR>  Second pair of paired-end reads

BATCH PAIRED READ INPUT:
      --l1 <FIRST_PAIR_LIST>   Text file containing paths to first pair of paired-end reads (one per line)
      --l2 <SECOND_PAIR_LIST>  Text file containing paths to second pair of paired-end reads (one per line)

INPUT:
  -l, --list-sequence <LIST_SEQUENCE>  File containing list of input sequences

MEMORY:
      --max-ram <MAX_RAM>  Maximum RAM usage in GB (default: 16)
```

### `inspect`: Inspect extracted .syldb and .sylsp files

**Usage:**
```
meta2bseek inspect -h

Inspect extracted .syldb and .sylsp files

Usage: meta2bseek inspect [OPTIONS] [FILES]...

Arguments:
  [FILES]...  Pre-extracted *.syldb/*.sylsp files.

Options:
  -o, --output-file <OUT_FILE_NAME>  Output to this file (YAML format). [default: stdout]
      --log-path <LOG_PATH>          Path to store TSV output files
      --tsv-name <TSV_NAME>          Name of the TSV file for tag count matrix [default: tag_matrix.tsv]
  -h, --help                         Print help
```

### `query`: Coverage-adjusted ANI querying between databases and samples
```
meta2bseek query -h

Coverage-adjusted ANI querying between databases and samples

Usage: meta2bseek query [OPTIONS] [FILES]...

Arguments:
  [FILES]...  Pre-extracted *.syldb/*.sylsp files. Raw single-end fastq/fasta are allowed and will be automatically extracted to .sylsp/.syldb

Options:
  -t <THREADS>
          Number of threads [default: 3]
  -s, --sample-threads <SAMPLE_THREADS>
          Number of samples to be processed concurrently. Default: (# of total threads / 3) + 1 for profile, 1 for query
      --trace
          Trace output (caution: very verbose)
      --debug
          Debug output
      --log-reassignments
          Output information for how 2bRAD tags for genomes are reassigned during `profile`. Caution: can be verbose and slows down computation.
  -h, --help
          Print help

INPUT/OUTPUT:
  -l, --list <FILE_LIST>             Newline delimited file of file inputs
  -o, --output-file <OUT_FILE_NAME>  Output to this file (TSV format). [default: stdout]

ALGORITHM:
      --min-count-correct <MIN_COUNT_CORRECT>
          Minimum 2bRAD tag multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity [default: 3]
  -M, --min-number-kmers <MIN_NUMBER_KMERS>
          Exclude genomes with less than this number of extracted 2bRAD tags [default: 50]
  -m, --minimum-ani <MINIMUM_ANI>
          Minimum adjusted ANI to consider (0-100). Default is 90 for query and 95 for profile. Smaller than 95 for profile will give inaccurate results.
  -u, --estimate-unknown
          Estimate true coverage and scale sequence abundance in `profile` by estimated unknown sequence percentage
  -I, --read-seq-id <SEQ_ID>
          Sequence identity (%) of reads. Only used in -u option and overrides automatic detection.

extracting:
  -r, --reads <READS>...                Single-end raw reads (fastx/gzip)
  -1, --first-pairs <FIRST_PAIR>...     First pairs for raw paired-end reads (fastx/gzip)
  -2, --second-pairs <SECOND_PAIR>...   Second pairs for raw paired-end reads (fastx/gzip)
  -c <C>                                Subsampling rate. Does nothing for pre-extracted files [default: 200]
  -i, --individual-records              Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-extracted files
      --min-spacing <MIN_SPACING_KMER>  Minimum spacing between selected 2bRAD tags on the database genomes. Does nothing for pre-extracted files [default: 30]
```

### `profile`: Species-level taxonomic profiling with abundances and ANIs

**Required Inputs**  
•	`.sylsp` file (from sample extraction).  
•	`.syldb` file (from reference genome extraction).  

**Optional Input**  

•	**Taxonomy File:** A TSV file with genome accessions and their corresponding taxonomy information. If not provided, results will display genome file names instead of species names.    

•	**Example Taxonomy File Format** (taxonomy.tsv):  

```
accession	taxonomy
RS_GCF_000657795.2	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Bordetella;s__Bordetella pseudohinzii
RS_GCF_001072555.1	d__Bacteria;p__Bacillota;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus epidermidis
```

**Example Command:**

```
meta2bseek profile \  
  --sample-file samples.sylsp \  
  --db-file reference_db.syldb \  
  --log-path /path/to/output \  
  --tsv-name profiling_results \  
  --threads 8 \  
  --taxonomy-file taxonomy.tsv \
  --minimum-ani 0.95 \
  --gscore-threshold 10
```

**Output:** Two TSV files will be generated in the output directory:

1.	`profiling_results.tsv`: Filtered results with gscore >= 10.
   
2.	`pre_gscore_filter_profiling_results.tsv`: Unfiltered results before applying the gscore filter.

**Usage:**

```
meta2bseek profile -h
Species-level taxonomic profiling with abundances and ANIs

Usage: meta2bseek profile [OPTIONS] --sample-file <SAMPLE_FILE> --db-file <DB_FILE>

Options:
      --sample-file <SAMPLE_FILE>
      --db-file <DB_FILE>
      --threads <THREADS>              [default: 1]
      --out-file-name <OUT_FILE_NAME>
      --log-path <LOG_PATH>
      --tsv-name <TSV_NAME>            [default: abundance_matrix.tsv]
      --taxonomy-file <TAXONOMY_FILE>  Taxonomy annotation file (e.g., taxonomy.txt) for species-level aggregation
  -h, --help                           Print help

ALGORITHM:
      --minimum-ani <MINIMUM_ANI>
          Minimum adjusted ANI to consider (0-100). Default is 95 for profile. Smaller than 95 for profile will give inaccurate results.
      --gscore-threshold <GSCORE_THRESHOLD>
          Minimum G-score threshold for species filtering. G-score = sqrt(reads_count * tag_count). Default is 10.0 [default: 10]
```

## How to interpret the results?
