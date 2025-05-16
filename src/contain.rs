use crate::cmdline::ContainArgs;
use anyhow::Result;

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