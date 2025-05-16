// src/inspect.rs

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. #[cfg(target_env = "musl")]

use clap::Parser;

mod cmdline;
mod constants;
mod extract;
mod contain;
mod inspect;

use crate::cmdline::{Cli, Mode};
use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc; //use std::panic::set_hook;

fn main() -> anyhow::Result<()> {

    let cli = Cli::parse();

    match cli.mode {
        Mode::Extract(extract_args) => extract::extract(extract_args),
        Mode::Query(contain_args) => contain::contain(contain_args, false),
        Mode::Profile(contain_args) => contain::contain(contain_args, true),
        Mode::Inspect(inspect_args) => inspect::inspect(inspect_args),
    }?;

    Ok(())
}