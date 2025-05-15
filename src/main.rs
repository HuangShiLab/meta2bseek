use anyhow::Result;
use clap::Parser;

// 仅在musl静态编译时启用jemalloc
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod cmdline;
mod extract_2brad;
mod query_module;

fn main() -> Result<()> {
    // 初始化日志系统
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    
    let cli = cmdline::Cli::parse();

    match cli.command {
        cmdline::Commands::Extract { 
            input,
            output,
            enzyme,
            threads,
            format,
            compress,
        } => extract_2brad::process_input(
            &input,
            &output,
            &enzyme,
            threads,
            &format,
            compress
        )?,

        cmdline::Commands::Query { 
        sample,
        database,
        threads,
        output
        } => {
            query_module::process_query(
                &sample,
                &database,
                threads,
                &output
            )?
        }
    }

    Ok(())
}