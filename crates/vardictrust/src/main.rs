mod cli;
mod profile;

use anyhow::Result;
use clap::Parser;

// Priority-based allocator selection: dhat wins when profile-dhat is enabled;
// otherwise jemalloc (if jemalloc feature enabled); otherwise system allocator.
#[cfg(feature = "profile-dhat")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;
#[cfg(all(feature = "jemalloc", not(feature = "profile-dhat")))]
#[global_allocator]
static ALLOC: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

// Thin-pointer _rjem_malloc_conf:
// &[u8; N] coerces to a thin pointer (8 bytes on 64-bit), matching jemalloc's
// expected C `const char *` ABI, unlike a fat slice `&[u8]` (16 bytes).
#[cfg(all(feature = "jemalloc", not(feature = "profile-heap")))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 60] =
    b"background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0";
#[cfg(feature = "profile-heap")]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 105] = b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0";

#[cfg(feature = "jemalloc")]
#[allow(clippy::size_of_ref)]
const _: () = assert!(core::mem::size_of_val(&MALLOC_CONF) == core::mem::size_of::<*const u8>());

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
    let _session = profile::Session::start()?;
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
