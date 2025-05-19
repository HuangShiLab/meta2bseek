use crate::cmdline::ContainArgs;
use anyhow::Result;

// 定义 ANI 计算结果的结构体
#[allow(dead_code)]
#[derive(Debug)]
pub struct AniResult {
    pub adjusted_ani: f64,        // 调整后的 ANI 值
    pub effective_coverage: f64,  // 有效覆盖度
    pub ani_confidence: (f64, f64), // ANI 的置信区间
    pub tag_containment: f64,     // 2bRAD 标签包含度
    pub naive_ani: f64,          // 原始 ANI 值
    pub g_score: f64,            // G 分数
    pub theoretical_tags: usize,  // 理论标签数
    pub sequenced_tags: usize,    // 测序标签数
}

impl AniResult {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self {
            adjusted_ani: 0.0,
            effective_coverage: 0.0,
            ani_confidence: (0.0, 0.0),
            tag_containment: 0.0,
            naive_ani: 0.0,
            g_score: 0.0,
            theoretical_tags: 0,
            sequenced_tags: 0,
        }
    }
}

// 计算 ANI 的主要函数
#[allow(dead_code)]
pub fn calculate_ani(
    shared_tags: usize,     // 共享的 2bRAD 标签数量
    query_tags: usize,      // 查询序列的标签总数
    _ref_tags: usize,        // 参考序列的标签总数
    coverage: f64,          // 覆盖度
    theoretical_tags: usize, // 理论标签数
) -> AniResult {
    let mut result = AniResult::new();
    
    // 计算 2bRAD 标签包含度
    result.tag_containment = shared_tags as f64 / query_tags as f64;
    
    // 计算原始 ANI（基于标签包含度）
    result.naive_ani = result.tag_containment * 100.0;
    
    // 计算有效覆盖度
    result.effective_coverage = coverage;
    
    // 计算 G 分数（几何平均数）
    result.g_score = (shared_tags as f64 * query_tags as f64).sqrt();
    
    // 记录标签统计
    result.theoretical_tags = theoretical_tags;
    result.sequenced_tags = shared_tags;
    
    // 计算调整后的 ANI
    // 使用统计方法调整 ANI，考虑覆盖度和 G 分数的影响
    result.adjusted_ani = adjust_ani_for_coverage_and_gscore(
        result.naive_ani,
        result.effective_coverage,
        result.g_score
    );
    
    // 计算 ANI 的置信区间
    result.ani_confidence = calculate_ani_confidence(
        result.adjusted_ani,
        result.effective_coverage,
        result.g_score
    );
    
    result
}

// 根据覆盖度和 G 分数调整 ANI 值
#[allow(dead_code)]
fn adjust_ani_for_coverage_and_gscore(naive_ani: f64, coverage: f64, g_score: f64) -> f64 {
    // 当覆盖度较低时，ANI 值会被低估
    // 使用经验公式进行调整，同时考虑 G 分数的影响
    let coverage_factor = 1.0 + (1.0 - coverage).powf(2.0) * 0.1;
    let g_score_factor = (g_score / 1000.0).min(1.0); // 归一化 G 分数的影响
    
    naive_ani * coverage_factor * (1.0 + g_score_factor * 0.05)
}

// 计算 ANI 的置信区间
#[allow(dead_code)]
fn calculate_ani_confidence(ani: f64, coverage: f64, g_score: f64) -> (f64, f64) {
    // 基于覆盖度和 G 分数计算置信区间
    // 覆盖度越低，G 分数越低，置信区间越宽
    let coverage_margin = (1.0 - coverage) * 0.5;
    let g_score_margin = (1.0 - (g_score / 1000.0).min(1.0)) * 0.3;
    let total_margin = coverage_margin + g_score_margin;
    
    (ani - total_margin, ani + total_margin)
}

pub fn contain(args: ContainArgs, is_profile: bool) -> Result<()> {
    // 打印基础参数
    println!("{} command parameters:", if is_profile { "Profile" } else { "Query" });
    println!("- Input files: {:?}", args.files);
    println!("- Threads: {}", args.threads);

    // 处理可选参数
    if let Some(file_list) = &args.file_list {
        println!("- File list: {}", file_list);
    }
    if let Some(min_ani) = args.minimum_ani {
        println!("- Minimum ANI: {}%", min_ani);
    }
    if args.estimate_unknown {
        println!("- Estimating unknown sequences");
    }

    // 打印调试选项
    if args.trace {
        println!("- Trace logging enabled");
    }
    if args.debug {
        println!("- Debug logging enabled");
    }

    Ok(())
}

/*
use crate::cmdline::*;
use std::path::Path;
use std::io::prelude::*;
use std::io;
use std::io::BufWriter;
use fxhash::FxHashMap;
use crate::constants::*;
use crate::inference::*;
use crate::extract::*;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult, pseudotax: bool, writer: &mut Box<dyn Write + Send>) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High {
        lambda_print = format!("HIGH");
    } else {
        lambda_print = format!("LOW");
    }
    let low_ani = ani_result.ani_ci.0;
    let high_ani = ani_result.ani_ci.1;
    let low_lambda = ani_result.lambda_ci.0;
    let high_lambda = ani_result.lambda_ci.1;

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none() {
        ci_ani = "NA-NA".to_string();
    } else {
        ci_ani = format!(
            "{:.2}-{:.2}",
            low_ani.unwrap() * 100.,
            high_ani.unwrap() * 100.
        );
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none() {
        ci_lambda = "NA-NA".to_string();
    } else {
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }


    //"Sample_file\tQuery_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",

    if !pseudotax{
        writeln!(writer, 
            "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.contig_name,
        ).expect("Error writing to file");
    }
    else{
        writeln!(writer,
            "{}\t{}\t{:.4}\t{:.4}\t{}\t{:.3}\t{}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{:.2}\t{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            ani_result.rel_abund.unwrap(),
            ani_result.seq_abund.unwrap(),
            print_final_ani,
            ani_result.final_est_cov,
            ci_ani,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.naive_ani * 100.,
            ani_result.kmers_lost.unwrap(),
            ani_result.contig_name,
        ).expect("Error writing to file");

    }
}

fn get_chunks(indices: &Vec<usize>, steps: usize) -> Vec<Vec<usize>>{
    let mut start = 0;
    let mut end = steps;
    let len = indices.len();
    let mut return_chunks = vec![];

    while start < len {
        if end > len {
            end = len;
        }

        let chunk: Vec<usize> = (start..end).collect();
        start = end;
        end += steps;
        return_chunks.push(chunk);
    }
    return_chunks
}

pub fn contain(mut args: ContainArgs, pseudotax_in: bool) {

    if pseudotax_in{
        args.pseudotax = true;
    }

    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    }
    else{
        level = log::LevelFilter::Info;
    }
    
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let out_writer = match args.out_file_name {
        Some(ref x) => {
            let path = Path::new(&x);
            Box::new(BufWriter::new(File::create(&path).unwrap())) as Box<dyn Write + Send>
        }
        None => Box::new(BufWriter::new(io::stdout())) as Box<dyn Write + Send>,
    };

    log::info!("Obtaining sketches...");
    let mut genome_sketch_files = vec![];
    let mut genome_files = vec![];
    let mut read_sketch_files = vec![];
    let mut read_files = vec![];

    let mut all_files = args.files.clone();

    if let Some(ref newline_file) = args.file_list{
        let file = File::open(newline_file).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            all_files.push(line.unwrap());
        }

    }


    for file in all_files.iter(){

        let mut genome_sketch_good_suffix = false;
        for suff in QUERY_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                genome_sketch_good_suffix = true;
                break
            }
        }

        let mut sample_sketch_good_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                sample_sketch_good_suffix = true;
                break
            }
        }

        if genome_sketch_good_suffix{
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix{
            read_sketch_files.push(file);
        } else if is_fasta(&file) {
            genome_files.push(file);
        } else if is_fastq(&file) {
            read_files.push(vec![file]);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    if args.first_pair.len() != args.second_pair.len() {
        error!("Different number of paired sequences (-1, -2) for sketching. Exiting.");
        std::process::exit(1);
    }

    // zip together the first and second pair files, push them to read_files
    for (first, second) in args.first_pair.iter().zip(args.second_pair.iter()) {
        read_files.push(vec![first,second]);
    }

    for read in args.reads.iter() {
        read_files.push(vec![read]);
    }

    if genome_sketch_files.is_empty() && genome_files.is_empty(){
        log::error!("No genome files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    if read_sketch_files.is_empty() && read_files.is_empty(){
        log::error!("No read files found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    let genome_sketches = get_genome_sketches(&args, &genome_sketch_files, &genome_files);
    let genome_index_vec = (0..genome_sketches.len()).collect::<Vec<usize>>();
    log::info!("Finished obtaining genome sketches.");

    if genome_sketches.is_empty() {
        log::error!("No genome sketches found; see sylph query/profile -h for help. Exiting");
        std::process::exit(1);
    }

    if genome_sketches.first().unwrap().pseudotax_tracked_nonused_kmers.is_none() && args.pseudotax{
        log::error!("Attempting profiling, but *.syldb was sketched with the --disable-profiling option. Exiting");
        std::process::exit(1);
    }

    let num_raw_read_files = read_files.len();
    let step;
    if let Some(sample_threads) = args.sample_threads{
        if sample_threads > 0{
            step = sample_threads;
        }
        else{
            step = 1;
        }
    }
    else{
        if args.pseudotax{
            step = usize::max(args.threads/3 + 1, usize::min(num_raw_read_files, args.threads))
        }
        else{
            step = usize::max(1, usize::min(num_raw_read_files, args.threads))
        }
    }

    let read_sketch_files_as_vec = read_sketch_files.clone().into_iter().map(|x| vec![x]).collect::<Vec<Vec<&String>>>();
    read_files.extend(read_sketch_files_as_vec);
    let sequence_index_vec = (0..read_files.len()).collect::<Vec<usize>>();
    let out_writer:Mutex<Box<dyn Write + Send>> = Mutex::new(out_writer);

    let chunks = get_chunks(&sequence_index_vec, step);

    print_header(args.pseudotax,&mut *out_writer.lock().unwrap(), args.estimate_unknown);
    chunks.into_iter().for_each(|chunk| {
        chunk.into_par_iter().for_each(|j|{
            let is_sketch = j >= read_files.len() - read_sketch_files.len();
            let sequence_sketch = get_seq_sketch(&args, &read_files[j], is_sketch, genome_sketches[0].c, genome_sketches[0].k);
            if sequence_sketch.is_some(){
                let first_read_file = read_files[j][0];
                let sequence_sketch = sequence_sketch.unwrap();
                
                let kmer_id_opt;
                if args.seq_id.is_some(){
                    kmer_id_opt = Some((args.seq_id.unwrap()/100.).powf(sequence_sketch.k as f64));
                }
                else{
                    kmer_id_opt = get_kmer_identity(&sequence_sketch, args.estimate_unknown);
                    log::debug!("{} has estimated identity {:.3}.", &first_read_file, kmer_id_opt.unwrap().powf(1./sequence_sketch.k as f64) * 100.);
                }
                
                let stats_vec_seq: Mutex<Vec<AniResult>> = Mutex::new(vec![]);
                genome_index_vec.par_iter().for_each(|i| {
                    let genome_sketch = &genome_sketches[*i];
                    let res = get_stats(&args, &genome_sketch, &sequence_sketch, None, args.log_reassignments);
                    if res.is_some() {
                        //res.as_mut().unwrap().genome_sketch_index = *i;
                        stats_vec_seq.lock().unwrap().push(res.unwrap());
                        
                    }
                });

                let mut stats_vec_seq = stats_vec_seq.into_inner().unwrap();
                estimate_true_cov(&mut stats_vec_seq, kmer_id_opt, args.estimate_unknown, sequence_sketch.mean_read_length, sequence_sketch.k);

                if args.pseudotax{
                    log::info!("{} taxonomic profiling; reassigning k-mers for {} genomes...", &first_read_file, stats_vec_seq.len());
                    let winner_map = winner_table(&stats_vec_seq, args.log_reassignments);
                    let remaining_genomes = stats_vec_seq.iter().map(|x| x.genome_sketch).collect::<Vec<&GenomeSketch>>();
                    let stats_vec_seq_2 = Mutex::new(vec![]);
                    remaining_genomes.into_par_iter().for_each(|genome_sketch|{
                        let res = get_stats(&args, &genome_sketch, &sequence_sketch, Some(&winner_map), args.log_reassignments);
                        if res.is_some() {
                            stats_vec_seq_2.lock().unwrap().push(res.unwrap());
                        }
                    });
                    stats_vec_seq = derep_if_reassign_threshold(&stats_vec_seq, stats_vec_seq_2.into_inner().unwrap(), args.redundant_ani, sequence_sketch.k);
                    //stats_vec_seq = stats_vec_seq_2.into_inner().unwrap();
                    estimate_true_cov(&mut stats_vec_seq, kmer_id_opt, args.estimate_unknown, sequence_sketch.mean_read_length, sequence_sketch.k);
                    log::info!("{} has {} genomes passing profiling threshold. ", &first_read_file, stats_vec_seq.len());

                    let mut bases_explained = 1.;
                    if args.estimate_unknown{
                        bases_explained = estimate_covered_bases(&stats_vec_seq, &sequence_sketch, sequence_sketch.mean_read_length, sequence_sketch.k);
                        log::info!("{} has {:.2}% of reads detected in database by profile", &first_read_file, bases_explained * 100.);
                    }

                    let total_cov = stats_vec_seq.iter().map(|x| x.final_est_cov).sum::<f64>();
                    let total_seq_cov = stats_vec_seq.iter().map(|x| x.final_est_cov * x.genome_sketch.gn_size as f64).sum::<f64>();
                    for thing in stats_vec_seq.iter_mut(){
                        thing.rel_abund = Some(thing.final_est_cov/total_cov * 100.);
                    }
                    for thing in stats_vec_seq.iter_mut(){
                        thing.seq_abund = Some(thing.final_est_cov * thing.genome_sketch.gn_size as f64 / total_seq_cov * 100. * bases_explained);
                    }
                }

                if args.pseudotax{
                    stats_vec_seq.sort_by(|x,y| y.rel_abund.unwrap().partial_cmp(&x.rel_abund.unwrap()).unwrap());
                }
                else{
                    stats_vec_seq.sort_by(|x,y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
                }
                
                let mut out_writer = out_writer.lock().unwrap();
                for res in stats_vec_seq{
                    print_ani_result(&res, args.pseudotax, &mut *out_writer);
                }
            }
            if read_files[j].len() > 1{
                log::info!("Finished paired sample {}.", &read_files[j][0]);
            }
            else{
                log::info!("Finished sample {}.", &read_files[j][0]);
            }
        });
    });

    log::info!("sylph finished.");
}

fn derep_if_reassign_threshold<'a>(results_old: &Vec<AniResult>, results_new: Vec<AniResult<'a>>, ani_thresh: f64, k: usize) -> Vec<AniResult<'a>>{
    let ani_thresh = ani_thresh/100.;

    let mut gn_sketch_to_contain = FxHashMap::default();
    for result in results_old.iter(){
        gn_sketch_to_contain.insert(result.genome_sketch, result);
    }

    let threshold = f64::powf(ani_thresh, k as f64);
    let mut return_vec = vec![];
    for result in results_new.into_iter(){
        let old_res = &gn_sketch_to_contain[result.genome_sketch];
        let num_kmer_reassign = (old_res.containment_index.0 - result.containment_index.0) as f64;
        let reass_thresh = threshold * result.containment_index.1 as f64;
        if num_kmer_reassign < reass_thresh{
            return_vec.push(result);
        }
        else{
            log::debug!("genome {} had num k-mers reassigned = {}, threshold was {}, removing.", result.gn_name, num_kmer_reassign, reass_thresh);
        }
    }
    return return_vec;
}

fn estimate_true_cov(results: &mut Vec<AniResult>, kmer_id_opt: Option<f64>, 
                     estimate_unknown: bool, read_length: f64, k: usize){
    let mut multiplier = 1.;
    if estimate_unknown{
        multiplier = read_length / (read_length - k as f64 + 1.);
    }
    if estimate_unknown && kmer_id_opt.is_some(){
        let id = kmer_id_opt.unwrap();
        for res in results.iter_mut(){
            res.final_est_cov = res.final_est_cov / id * multiplier ;
        }
    }
}

fn estimate_covered_bases(results: &Vec<AniResult>, sequence_sketch: &SequencesSketch, read_length: f64, k: usize) -> f64{
    let multiplier = read_length / (read_length - (k as f64) + 1.);

    let mut num_covered_bases = 0.;
    for res in results.iter(){
        num_covered_bases += (res.genome_sketch.gn_size as f64) * res.final_est_cov
    }
    let mut num_total_counts = 0;
    for count in sequence_sketch.kmer_counts.values(){
        num_total_counts += *count as usize;
    }
    let num_tentative_bases = sequence_sketch.c * num_total_counts;
    let num_tentative_bases = num_tentative_bases as f64 * multiplier;
    if num_tentative_bases == 0.{
        return 0.;
    }
    return f64::min(num_covered_bases as f64 / num_tentative_bases, 1.);
}

fn winner_table<'a>(results : &'a Vec<AniResult>, log_reassign: bool) -> FxHashMap<Kmer, (f64,&'a GenomeSketch, bool)> {
    let mut kmer_to_genome_map : FxHashMap<_,_> = FxHashMap::default();
    for res in results.iter(){
        //let gn_sketch = &genome_sketches[res.genome_sketch_index];
        let gn_sketch = res.genome_sketch;
        for kmer in gn_sketch.genome_kmers.iter(){
            let v = kmer_to_genome_map.entry(*kmer).or_insert((res.final_est_ani, res.genome_sketch, false));
            if res.final_est_ani > v.0{
                *v = (res.final_est_ani, gn_sketch, true);
            }
        }
        
        if gn_sketch.pseudotax_tracked_nonused_kmers.is_some(){
            for kmer in gn_sketch.pseudotax_tracked_nonused_kmers.as_ref().unwrap().iter(){
                let v = kmer_to_genome_map.entry(*kmer).or_insert((res.final_est_ani, res.genome_sketch, false));
                if res.final_est_ani > v.0{
                    *v = (res.final_est_ani, gn_sketch, true);
                }
            }
        }
    }

    //log reassigned kmers
    if log_reassign{
        log::info!("------------- Logging k-mer reassignments -----------------");
        let mut sketch_to_index = FxHashMap::default();
        for (i,res) in results.iter().enumerate(){
            log::info!("Index\t{}\t{}\t{}", i, res.genome_sketch.file_name, res.genome_sketch.first_contig_name);
            sketch_to_index.insert(res.genome_sketch, i);
        }
        (0..results.len()).into_par_iter().for_each(|i|{
            let res = &results[i];
            let mut reassign_edge_map = FxHashMap::default();
            for kmer in res.genome_sketch.genome_kmers.iter(){
                let value = kmer_to_genome_map[kmer].1;
                if value != res.genome_sketch{
                    let edge_count = reassign_edge_map.entry((sketch_to_index[value],i)).or_insert(0);
                    *edge_count += 1;
                }
            }
            for (key,val) in reassign_edge_map{
                if val > 10{
                    log::info!("{}->{}\t{}\tkmers reassigned", key.0, key.1, val);
                }
            }
        });
    }

    return kmer_to_genome_map;
}

fn print_header(pseudotax: bool, writer: &mut Box<dyn Write + Send>, estimate_unknown: bool) {
    if !pseudotax{
        writeln!(writer,
            //"Sample_file\tQuery_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            "Sample_file\tGenome_file\tAdjusted_ANI\tEff_cov\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tContig_name",
            ).expect("Error writing to file.");
    }
    else{
        let cov_head;
        if estimate_unknown{
            cov_head = "True_cov";
        }
        else{
            cov_head = "Eff_cov";
        }
        writeln!(writer,
            "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\t{}\tANI_5-95_percentile\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tNaive_ANI\tkmers_reassigned\tContig_name", cov_head
            ).expect("Error writing to file.");
    }
}

fn get_genome_sketches(
    args: &ContainArgs,
    genome_sketch_files: &Vec<&String>,
    genome_files: &Vec<&String>,
) -> Vec<GenomeSketch> {
    let mut lowest_genome_c = None;
    let mut current_k = None;

    let genome_sketches = Mutex::new(vec![]);

    for genome_sketch_file in genome_sketch_files {
        let file = File::open(genome_sketch_file).expect(&format!("The sketch `{}` could not be opened. Exiting", genome_sketch_file));
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
                &genome_sketch_file
            ));
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    return genome_sketches.into_inner().unwrap();
}

fn get_seq_sketch(
    args: &ContainArgs,
    read_file: &Vec<&String>,
    is_sketch_file: bool,
    genome_c: usize,
    genome_k: usize,
) -> Option<SequencesSketch> {
    if is_sketch_file {
        let read_file = read_file[0];
        let read_sketch_file = read_file;
        let file = File::open(read_sketch_file).expect(&format!(
            "The sketch `{}` could not be opened",
            &read_sketch_file
        ));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch: SequencesSketch = bincode::deserialize_from(read_reader).expect(
            &format!("The sketch `{}` is not a valid sketch. Perhaps it is an older incompatible version ", read_sketch_file),
        );
        if read_sketch.c > genome_c {
            error!("{} value of -c is {}; this is greater than the smallest value of -c = {} for a genome sketch. Exiting.", read_file, read_sketch.c, genome_c);
            return None;
        }
        else if read_sketch.c < genome_c{
            info!("{} value of -c for reads is {}; this is smaller than the -c for a genome sketch. Using the larger -c {} instead.", read_file, read_sketch.c,  genome_c);
        }

        return Some(read_sketch);
    } else {
        if args.c > genome_c{
            info!("{} value of -c for reads is {}; this is smaller than the -c for a genome sketch. Using the larger -c {} instead.", read_file[0], args.c,  genome_c);
        }
        if genome_c < args.c {
            error!("{} error: value of -c for contain = {} -- greater than the smallest value of -c for a genome sketch = {}. Continuing without sketching.", read_file[0], args.c, genome_c);
            return None;
        } 
        else if genome_k != args.k {
            error!(
                "{} -k {} is not equal to -k {} found in sketches. Continuing without sketching.",
                read_file[0], args.k, genome_k
            );
            return None;
        } else {
            if read_file.len() == 1{
                let read_sketch_opt = sketch_sequences_needle(&read_file[0], args.c, args.k, None, false);
                return read_sketch_opt;
            }
            else if read_file.len() == 2{
                let read_sketch_opt = sketch_pair_sequences(&read_file[0], &read_file[1], args.c, args.k, None, false, DEFAULT_FPR);
                return read_sketch_opt;
            }
            else{
                panic!("Internal Error: read_file has length {}. Something went wrong...", read_file.len());
            }
        }
    }
}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &SequencesSketch,
    winner_map: Option<&FxHashMap<Kmer, (f64,& GenomeSketch, bool)>>,
    log_reassign: bool
) -> Option<AniResult<'a>> {
    if genome_sketch.k != sequence_sketch.k {
        log::error!(
            "k parameter for reads {} != k parameter for genome {}",
            sequence_sketch.k,
            genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c < sequence_sketch.c {
        log::error!(
            "c parameter for reads {} > c parameter for genome {}",
            sequence_sketch.c,
            genome_sketch.c
        );
        std::process::exit(1);
    }
    let mut contain_count = 0;
    let mut covs = vec![];
    let gn_kmers = &genome_sketch.genome_kmers;
    if (gn_kmers.len() as f64) < args.min_number_kmers{
        return None
    }

    let mut kmers_lost_count = 0;
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            if sequence_sketch.kmer_counts[kmer] == 0{
                continue
            }
            if winner_map.is_some(){
                let map = &winner_map.unwrap();
                if map[kmer].1 != genome_sketch{
                    kmers_lost_count += 1;
                    continue
                }
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);

            }
            else{
                contain_count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            }
        }
    }

    if covs.is_empty() {
        return None;
    }
    let naive_ani = f64::powf(
        contain_count as f64 / gn_kmers.len() as f64,
        1. / genome_sketch.k as f64,
    );
    covs.sort();
    //let covs = &covs[0..covs.len() * 99 / 100];
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov).unwrap();
    let mut max_cov = f64::MAX;
    if median_cov < 30.{
        for i in covs.len() / 2..covs.len(){
            let cov = covs[i];
            if pois.cdf(cov.into()) < CUTOFF_PVALUE {
                max_cov = cov as f64;
            } else {
                break;
            }
        }
    }
    log::trace!("COV VECTOR for {}/{}: {:?}, MAX_COV_THRESHOLD: {}", sequence_sketch.file_name, genome_sketch.file_name ,covs, max_cov);


    let mut full_covs = vec![0; gn_kmers.len() - contain_count];
    for cov in covs.iter() {
        if (*cov as f64) <= max_cov {
            full_covs.push(*cov);
        }
    }
    let var = var(&full_covs);
    if var.is_some(){
        log::trace!("VAR {} {}", var.unwrap(), genome_sketch.file_name);
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    let geq1_mean_cov = full_covs.iter().sum::<u32>() as f64 / covs.len() as f64;

    let use_lambda;
    if median_cov > MEDIAN_ANI_THRESHOLD {
        use_lambda = AdjustStatus::High
    } else {
        let test_lambda;
        if args.ratio {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        } else if args.mme {
            test_lambda = mme_lambda(&full_covs)
        } else if args.nb {
            test_lambda = binary_search_lambda(&full_covs)
        } else if args.mle {
            test_lambda = mle_zip(&full_covs, sequence_sketch.k as f64)
        } else {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        };
        if test_lambda.is_none() {
            use_lambda = AdjustStatus::Low
        } else {
            use_lambda = AdjustStatus::Lambda(test_lambda.unwrap());
        }
    }

    let final_est_cov;

    if let AdjustStatus::Lambda(lam) = use_lambda {
        final_est_cov = lam
    } else if median_cov < MAX_MEDIAN_FOR_MEAN_FINAL_EST{
        final_est_cov = geq1_mean_cov;
    } else{
        if args.mean_coverage{
            final_est_cov = geq1_mean_cov;
        }
        else{
            final_est_cov = median_cov;
        }
    }

    let opt_lambda;
    if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {
        opt_lambda = None
    } else {
        opt_lambda = Some(final_est_cov)
    };

    let opt_est_ani = ani_from_lambda(opt_lambda, mean_cov, sequence_sketch.k as f64, &full_covs);
    
    let final_est_ani;
    if opt_lambda.is_none() || opt_est_ani.is_none() || args.no_adj {
        final_est_ani = naive_ani;
    } else {
        final_est_ani = opt_est_ani.unwrap();
    }

    let min_ani = if args.minimum_ani.is_some() {args.minimum_ani.unwrap()/100. }
        else if args.pseudotax { MIN_ANI_P_DEF } 
        else { MIN_ANI_DEF };
    if final_est_ani < min_ani {
        if winner_map.is_some(){
            //Used to be > min ani, now it is not after reassignment
            if log_reassign{
                log::info!("Genome/contig {}/{} has ANI = {} < {} after reassigning {} k-mers ({} contained k-mers after reassign)", 
                    genome_sketch.file_name,
                    genome_sketch.first_contig_name,
                    final_est_ani * 100.,
                    min_ani * 100.,
                    kmers_lost_count,
                    contain_count)
            }

        }
        return None;
    }

    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if !args.no_ci && opt_lambda.is_some() {
        let bootstrap = bootstrap_interval(&full_covs, sequence_sketch.k as f64, &args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    }

    
    let seq_name;
    if let Some(sample) = &sequence_sketch.sample_name{
        seq_name = sample.clone();
    }
    else{
        seq_name = sequence_sketch.file_name.clone();
    }

    let kmers_lost;
    if winner_map.is_some(){
        kmers_lost = Some(kmers_lost_count)
    }
    else{
        kmers_lost = None;
    }

    let ani_result = AniResult {
        naive_ani,
        final_est_ani,
        final_est_cov,
        seq_name: seq_name,
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov: geq1_mean_cov,
        median_cov,
        containment_index: (contain_count, gn_kmers.len()),
        lambda: use_lambda,
        ani_ci: (low_ani, high_ani),
        lambda_ci: (low_lambda, high_lambda),
        genome_sketch: genome_sketch,
        rel_abund: None,
        seq_abund: None,
        kmers_lost: kmers_lost,

    };
    //log::trace!("Other time {:?}", Instant::now() - start_t_initial);

    return Some(ani_result);
}


fn ani_from_lambda(lambda: Option<f64>, _mean: f64, k: f64, full_cov: &[u32]) -> Option<f64> {
    if lambda.is_none() {
        return None;
    }
    let mut contain_count = 0;
    let mut _zero_count = 0;
    for x in full_cov {
        if *x != 0 {
            contain_count += 1;
        } else {
            _zero_count += 1;
        }
    }

    let lambda = lambda.unwrap();
    let adj_index =
        contain_count as f64 / (1. - f64::exp(-lambda)) / full_cov.len() as f64;
    let ret_ani;
    //let ani = f64::powf(1. - pi, 1./k);
    let ani = f64::powf(adj_index, 1. / k);
    if ani < 0. || ani.is_nan() {
        ret_ani = None;
    } else {
        if ani > 1. {
            ret_ani = Some(ani)
        } else {
            ret_ani = Some(ani);
        }
    }
    return ret_ani;
}

fn bootstrap_interval(
    covs_full: &Vec<u32>,
    k: f64,
    args: &ContainArgs,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    fastrand::seed(7);
    let num_samp = covs_full.len();
    let iters = 100;
    let mut res_ani = vec![];
    let mut res_lambda = vec![];

    for _ in 0..iters {
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp {
            rand_vec.push(covs_full[fastrand::usize(..covs_full.len())]);
        }
        let lambda;
        if args.ratio {
            lambda = ratio_lambda(&rand_vec, args.min_count_correct);
        } else if args.mme {
            lambda = mme_lambda(&rand_vec);
        } else if args.nb {
            lambda = binary_search_lambda(&rand_vec);
        } else if args.mle {
            lambda = mle_zip(&rand_vec, k);
        } else {
            lambda = ratio_lambda(&rand_vec,args.min_count_correct);
        }
        let ani = ani_from_lambda(lambda, mean(&rand_vec).unwrap().into(), k, &rand_vec);
        if ani.is_some() && lambda.is_some() {
            if !ani.unwrap().is_nan() && !lambda.unwrap().is_nan() {
                res_ani.push(ani);
                res_lambda.push(lambda);
            }
        }
    }
    res_ani.sort_by(|x, y| x.partial_cmp(y).unwrap());
    res_lambda.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if res_ani.len() < 50 {
        return (None, None, None, None);
    }
    let suc = res_ani.len();
    let low_ani = res_ani[suc * 5 / 100 - 1];
    let high_ani = res_ani[suc * 95 / 100 - 1];
    let low_lambda = res_lambda[suc * 5 / 100 - 1];
    let high_lambda = res_lambda[suc * 95 / 100 - 1];

    return (low_ani, high_ani, low_lambda, high_lambda);
}


fn get_kmer_identity(seq_sketch: &SequencesSketch, estimate_unknown: bool) -> Option<f64>{

    if !estimate_unknown{
        return None
    }

    let mut median = 0;
    let mut mov_avg_median = 0.;
    let mut n = 1.;
    for count in seq_sketch.kmer_counts.values(){
        if *count > 1{
            if *count > median{
                median += 1;
            }
            else{
                median -= 1;
            }
            mov_avg_median += median as f64;
            n += 1.;
        }
    }

    mov_avg_median /= n;
    log::debug!("Estimated continuous median k-mer count for {} is {:.3}", &seq_sketch.file_name, mov_avg_median);
    
    let mut num_1s = 0;
    let mut num_not1s = 0;
    for count in seq_sketch.kmer_counts.values(){
        if *count == 1{
            num_1s += 1;
        }
        else{
            num_not1s += *count;
        }
    }
    //0.1 so no div by 0 error
    let eps = num_not1s as f64 / (num_not1s as f64 + num_1s as f64 + 0.1);
    //dbg!("Automatic id est, 1-to-2 ratio, 2-to-3", eps.powf(1./31.), num_1s as f64 / num_2s as f64, two_to_three);

    if mov_avg_median < MED_KMER_FOR_ID_EST && seq_sketch.mean_read_length < 400.{
        log::info!("{} short-read sample has high diversity compared to sequencing depth (approx. avg depth < 3). Using 99.5% as read accuracy estimate instead of automatic detection for --estimate-unknown.", &seq_sketch.file_name);
        return Some(0.995f64.powf(seq_sketch.k as f64));
    }

    if eps < 1.{
        return Some(eps)
    }
    else{
        return Some(1.)
    }
}
*/