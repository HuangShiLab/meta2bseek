use crate::cmdline::ContainArgs;
use anyhow::Result;

// 定义 ANI 计算结果的结构体
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
pub fn calculate_ani(
    shared_tags: usize,     // 共享的 2bRAD 标签数量
    query_tags: usize,      // 查询序列的标签总数
    ref_tags: usize,        // 参考序列的标签总数
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
fn adjust_ani_for_coverage_and_gscore(naive_ani: f64, coverage: f64, g_score: f64) -> f64 {
    // 当覆盖度较低时，ANI 值会被低估
    // 使用经验公式进行调整，同时考虑 G 分数的影响
    let coverage_factor = 1.0 + (1.0 - coverage).powf(2.0) * 0.1;
    let g_score_factor = (g_score / 1000.0).min(1.0); // 归一化 G 分数的影响
    
    naive_ani * coverage_factor * (1.0 + g_score_factor * 0.05)
}

// 计算 ANI 的置信区间
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