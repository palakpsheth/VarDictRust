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
#[expect(
    clippy::size_of_ref,
    reason = "intentionally asserts the size of the static reference itself (thin 8-byte pointer vs 16-byte slice)"
)]
const _: () = assert!(core::mem::size_of_val(&MALLOC_CONF) == core::mem::size_of::<*const u8>());

fn should_ignore_broken_pipe(
    cfg: &vardictrust_core::config::Config,
    err: &vardictrust_core::types::VardictError,
) -> bool {
    matches!(cfg.output_path.as_deref(), None | Some("-")) && err.is_broken_pipe()
}

fn finish_stdout_utility(res: std::io::Result<()>) -> Result<()> {
    match res {
        Ok(()) => Ok(()),
        Err(err) if err.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(err) => Err(err.into()),
    }
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

fn sanitize_for_terminal(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        let code = c as u32;
        if c != '\n' && c != '\t' && (code < 0x20 || c == '\x7f' || (0x80..=0x9f).contains(&code)) {
            for esc in c.escape_default() {
                out.push(esc);
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn write_error_ignoring_broken_pipe<W: std::io::Write>(
    err: &anyhow::Error,
    mut out: W,
) -> std::io::Result<()> {
    let safe_msg = sanitize_for_terminal(&format!("{err:?}"));
    let line = format!("Error: {safe_msg}\n");
    match out.write_all(line.as_bytes()).and_then(|()| out.flush()) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(e),
    }
}

fn try_main() -> Result<()> {
    env_logger::init();
    let _session = profile::Session::start()?;
    let raw_args = collect_utf8_args(std::env::args_os())?;
    let args = cli::normalize_args(&raw_args);
    cli::check_unsupported(&args)?;
    let args = cli::Cli::parse_from(args);
    if let Some(shell) = args.completions.first().copied() {
        let mut stdout = std::io::BufWriter::new(std::io::stdout().lock());
        let res = cli::completions_man::write_completions(shell, &mut stdout)
            .and_then(|()| std::io::Write::flush(&mut stdout));
        return finish_stdout_utility(res);
    }
    if args.manpage {
        let mut stdout = std::io::BufWriter::new(std::io::stdout().lock());
        let res = cli::completions_man::write_manpage(&mut stdout)
            .and_then(|()| std::io::Write::flush(&mut stdout));
        return finish_stdout_utility(res);
    }
    let cfg = match std::env::var("RAYON_NUM_THREADS").ok() {
        Some(rayon_env) => args.to_config_with_rayon_env(Some(rayon_env.as_str()))?,
        None => args.to_config()?,
    };
    if let Err(err) = vardictrust_core::run(&cfg) {
        if should_ignore_broken_pipe(&cfg, &err) {
            return Ok(());
        }
        return Err(err.into());
    }
    Ok(())
}

fn main() -> std::process::ExitCode {
    match try_main() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            let _ = write_error_ignoring_broken_pipe(&err, std::io::stderr().lock());
            std::process::ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch12_sanitize_for_terminal_escapes_control_sequences() {
        let hostile =
            "line1\tcol2\n\x1b]52;c;SGVsbG8=\x07\u{9b}2J\roverwrite\0nul\x7fdel λ café 🧬";
        let got = sanitize_for_terminal(hostile);
        assert_eq!(
            got,
            "line1\tcol2\n\\u{1b}]52;c;SGVsbG8=\\u{7}\\u{9b}2J\\roverwrite\\u{0}nul\\u{7f}del λ café 🧬"
        );
        assert!(got.contains('\n') && got.contains('\t'));
        assert!(!got.chars().any(|c| {
            c != '\n'
                && c != '\t'
                && ((c as u32) < 0x20 || c == '\x7f' || (0x80..=0x9f).contains(&(c as u32)))
        }));

        let err = anyhow::anyhow!("inner \x1b[31mred\x07\r\0").context("outer\tctx");
        let mut sink = Vec::new();
        write_error_ignoring_broken_pipe(&err, &mut sink).unwrap();
        let rendered = String::from_utf8(sink).unwrap();
        assert!(rendered.starts_with("Error: outer\tctx\n"));
        assert!(rendered.ends_with('\n'));
        assert!(rendered.contains("inner \\u{1b}[31mred\\u{7}\\r\\u{0}"));
        assert!(!rendered.contains('\x1b'));

        struct BrokenPipeWriter;
        impl std::io::Write for BrokenPipeWriter {
            fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::from(std::io::ErrorKind::BrokenPipe))
            }
        }
        assert!(write_error_ignoring_broken_pipe(&err, BrokenPipeWriter).is_ok());
    }

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
        cfg.output_path = Some("-".to_string());
        assert!(should_ignore_broken_pipe(&cfg, &bp));
        cfg.output_path = Some("out.fifo".to_string());
        assert!(!should_ignore_broken_pipe(&cfg, &bp));
    }
}
