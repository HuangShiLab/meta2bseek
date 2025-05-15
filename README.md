# meta2bseek

Meta2bSeek is composed of tools in three categories: 
* taxonomic profiling tool (2bRAD-M2); 
* strain and functional profiling tool (Strain2bFunc); 
* statistical analysis tools. 

2bRAD-M2: We engineered a Rust-optimized pipeline for ultrafast species identification using k-mers selected by type 2b sites, achieving algorithmic acceleration through parallel computing and memory-efficient indexing. A comprehensive database spanning 394921 bacterial, 7777 archaeal, and 11184 fungal genomes was curated. Cross-kingdom phylogenetic relationships were reconstructed using 2b-site mash distance-based hierarchical clustering. 

Strain2bFunc: A novel k-mer hashing strategy was implemented to identify both known and novel strains at a high ANI resolution, enabling to map uncharacterized strains to their closest phylogenetic neighbors in our expanded reference database.

