mod cli;

use anyhow::Result;
use clap::Parser;

// Task 4 allocator experiment: jemalloc global allocator. glibc ptmalloc
// retains freed per-read churn spans (measured: 267 KB live vs 46 MB HWM on
// the chrM slice, 173x); jemalloc's background purging + decay returns them.
// Baked-in config: background purging on, fast dirty decay, no muzzy
// retention, arenas = host CPUs (box-tuned: revisit before merge).
// Toggle off with --no-default-features. No application code changes;
// allocator swaps cannot change output (verified by gates).
#[cfg(feature = "jemalloc")]
#[global_allocator]
static ALLOC: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(feature = "jemalloc")]
#[allow(non_upper_case_globals)]
#[export_name = "_rjem_malloc_conf"]
pub static malloc_conf: &[u8] =
    b"background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0,narenas:20\0";

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
