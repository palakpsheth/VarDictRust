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
// M-DOCUMENTED-MAGIC: jemalloc's `background_thread:true` relies on Linux pthread
// background decay threads; on non-Linux targets (e.g. macOS Darwin aarch64/x86_64)
// jemalloc rejects `background_thread:true`, prints `<jemalloc>: option
// background_thread currently supports pthread only` to stderr on every invocation,
// and aborts parsing the remainder of `_rjem_malloc_conf`. Target-splitting to
// `background_thread:false` on non-Linux preserves `dirty_decay_ms`, `muzzy_decay_ms`,
// `narenas:8`, and heap-profiling options without stderr noise or fallback to defaults.
const MALLOC_CONF_LINUX_DEFAULT: &[u8; 69] =
    b"background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0";
const MALLOC_CONF_NON_LINUX_DEFAULT: &[u8; 70] =
    b"background_thread:false,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0";
const MALLOC_CONF_LINUX_PROF: &[u8; 105] = b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0";
const MALLOC_CONF_NON_LINUX_PROF: &[u8; 106] = b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:false,dirty_decay_ms:5000,muzzy_decay_ms:0\0";

#[cfg(all(
    feature = "jemalloc",
    not(feature = "profile-heap"),
    target_os = "linux"
))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 69] = MALLOC_CONF_LINUX_DEFAULT;

#[cfg(all(
    feature = "jemalloc",
    not(feature = "profile-heap"),
    not(target_os = "linux")
))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 70] = MALLOC_CONF_NON_LINUX_DEFAULT;

#[cfg(all(feature = "profile-heap", target_os = "linux"))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 105] = MALLOC_CONF_LINUX_PROF;

#[cfg(all(feature = "profile-heap", not(target_os = "linux")))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 106] = MALLOC_CONF_NON_LINUX_PROF;

const fn malloc_conf_bytes_for_target(is_linux: bool, profile_heap: bool) -> &'static [u8] {
    match (is_linux, profile_heap) {
        (true, false) => MALLOC_CONF_LINUX_DEFAULT,
        (false, false) => MALLOC_CONF_NON_LINUX_DEFAULT,
        (true, true) => MALLOC_CONF_LINUX_PROF,
        (false, true) => MALLOC_CONF_NON_LINUX_PROF,
    }
}

const _: () = {
    assert!(malloc_conf_bytes_for_target(true, false).len() == 69);
    assert!(malloc_conf_bytes_for_target(false, false).len() == 70);
    assert!(malloc_conf_bytes_for_target(true, true).len() == 105);
    assert!(malloc_conf_bytes_for_target(false, true).len() == 106);
};

#[cfg(feature = "jemalloc")]
// Lint false positive: size_of_val(&MALLOC_CONF) intentionally asserts the size of the static reference itself (thin 8-byte pointer vs 16-byte slice).
#[allow(clippy::size_of_ref)]
const _: () = assert!(core::mem::size_of_val(&MALLOC_CONF) == core::mem::size_of::<*const u8>());

fn should_ignore_broken_pipe(
    cfg: &vardictrust_core::config::Config,
    err: &vardictrust_core::types::VardictError,
) -> bool {
    cfg.output_path.is_none() && err.is_broken_pipe()
}

fn collect_utf8_args<I: IntoIterator<Item = std::ffi::OsString>>(
    iter: I,
) -> std::result::Result<Vec<String>, vardictrust_core::types::VardictError> {
    iter.into_iter()
        .enumerate()
        .map(|(idx, os)| {
            os.into_string().map_err(|bad| {
                vardictrust_core::types::VardictError::Parse(format!(
                    "CLI argument #{idx} ({}) is not valid UTF-8",
                    bad.to_string_lossy()
                ))
            })
        })
        .collect()
}

fn main() -> Result<()> {
    env_logger::init();
    let _session = profile::Session::start()?;
    let raw_args = collect_utf8_args(std::env::args_os())?;
    let args = cli::normalize_args(&raw_args);
    cli::check_unsupported(&args)?;
    let args = cli::Cli::parse_from(args);
    let cfg = args.to_config()?;
    if let Err(err) = vardictrust_core::run(&cfg) {
        if should_ignore_broken_pipe(&cfg, &err) {
            return Ok(());
        }
        return Err(err.into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malloc_conf_target_split_and_nul_termination() {
        let cases: &[(bool, bool, usize, &[u8])] = &[
            (
                true,
                false,
                69,
                b"background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0",
            ),
            (
                false,
                false,
                70,
                b"background_thread:false,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0",
            ),
            (
                true,
                true,
                105,
                b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0",
            ),
            (
                false,
                true,
                106,
                b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:false,dirty_decay_ms:5000,muzzy_decay_ms:0\0",
            ),
        ];

        for &(is_linux, profile_heap, expected_len, expected) in cases {
            let actual = malloc_conf_bytes_for_target(is_linux, profile_heap);
            assert_eq!(
                actual, expected,
                "is_linux={is_linux}, profile_heap={profile_heap}"
            );
            assert_eq!(actual.len(), expected_len);
            assert_eq!(actual.last(), Some(&0));
            assert!(!actual[..actual.len() - 1].contains(&0));
            let s = core::str::from_utf8(&actual[..actual.len() - 1]).expect("valid ASCII");
            if is_linux {
                assert!(s.contains("background_thread:true"));
                assert!(!s.contains("background_thread:false"));
            } else {
                assert!(s.contains("background_thread:false"));
                assert!(!s.contains("background_thread:true"));
            }
        }

        #[cfg(feature = "jemalloc")]
        {
            #[allow(clippy::size_of_ref)]
            let ptr_size = core::mem::size_of_val(&MALLOC_CONF);
            assert_eq!(ptr_size, core::mem::size_of::<*const u8>());
            assert_eq!(
                &MALLOC_CONF[..],
                malloc_conf_bytes_for_target(
                    cfg!(target_os = "linux"),
                    cfg!(feature = "profile-heap")
                )
            );
            assert_eq!(MALLOC_CONF.last(), Some(&0));
            assert!(!MALLOC_CONF[..MALLOC_CONF.len() - 1].contains(&0));
        }

        #[cfg(unix)]
        {
            use std::ffi::OsString;
            use std::os::unix::ffi::OsStringExt as _;
            let bad_os = vec![
                OsString::from("vardictrust"),
                OsString::from_vec(vec![0xff, b'x']),
            ];
            let err = collect_utf8_args(bad_os).expect_err("non-UTF-8 argv must fail cleanly");
            assert!(err.to_string().contains("not valid UTF-8"));
        }

        let bp = vardictrust_core::types::VardictError::Io(std::io::Error::from(
            std::io::ErrorKind::BrokenPipe,
        ));
        let mut cfg = vardictrust_core::config::Config::default();
        assert!(should_ignore_broken_pipe(&cfg, &bp));
        cfg.output_path = Some("out.fifo".to_string());
        assert!(!should_ignore_broken_pipe(&cfg, &bp));
    }
}
