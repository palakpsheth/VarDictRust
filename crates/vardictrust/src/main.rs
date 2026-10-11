mod cli;

use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    env_logger::init();
    let raw_args: Vec<String> = std::env::args().collect();
    let args = cli::normalize_args(&raw_args);
    cli::check_unsupported(&args)?;
    let args = cli::Cli::parse_from(args);
    let cfg = args.to_config()?;
    vardictrust_core::run(&cfg)?;
    Ok(())
}
