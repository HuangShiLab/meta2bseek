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

## How to run?

## How to interpret the results?
