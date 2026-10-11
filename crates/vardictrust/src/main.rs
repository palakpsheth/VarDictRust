mod cli;

use anyhow::Result;
use clap::Parser;

fn is_broken_pipe(err: &vardictrust_core::types::VardictError) -> bool {
    match err {
        vardictrust_core::types::VardictError::Io(e) => {
            e.kind() == std::io::ErrorKind::BrokenPipe || e.raw_os_error() == Some(32)
        }
        _ => false,
    }
}

fn main() -> Result<()> {
    env_logger::init();
    let raw_args: Vec<String> = std::env::args().collect();
    let args = cli::normalize_args(&raw_args);
    cli::check_unsupported(&args)?;
    let args = cli::Cli::parse_from(args);
    let cfg = args.to_config()?;
    if let Err(err) = vardictrust_core::run(&cfg) {
        if is_broken_pipe(&err) {
            return Ok(());
        }
        return Err(err.into());
    }
    Ok(())
}
