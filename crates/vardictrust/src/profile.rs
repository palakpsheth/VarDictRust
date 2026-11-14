#![deny(unsafe_op_in_unsafe_fn)]

use anyhow::Result;
use std::path::PathBuf;

pub struct Session {
    #[cfg(feature = "profile-cpu")]
    cpu: Option<CpuSession>,
    #[cfg(feature = "profile-dhat")]
    _dhat: Option<dhat::Profiler>,
    #[cfg(feature = "profile-dhat")]
    dhat_paths: Option<(PathBuf, PathBuf)>,
    #[cfg(feature = "profile-heap")]
    heap_dir: Option<PathBuf>,
}

#[cfg(feature = "profile-cpu")]
struct CpuSession {
    guard: pprof::ProfilerGuard<'static>,
    dir: PathBuf,
    old_act: Option<libc::sigaction>,
}

#[cfg(feature = "profile-cpu")]
static TICKS_TOTAL: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(feature = "profile-cpu")]
static OLD_ACT: std::sync::atomic::AtomicPtr<libc::sigaction> =
    std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());

#[cfg(feature = "profile-cpu")]
#[inline]
fn is_chainable_sigprof_handler(addr: usize) -> bool {
    addr != 0
        && addr != libc::SIG_IGN
        && addr != libc::SIG_DFL
        && addr != sigprof_handler as *const () as usize
}

// `sigprof_handler` is an `extern "C"` signal handler invoked by the OS kernel on `SIGPROF`.
// `TICKS_TOTAL` and `OLD_ACT` use lock-free atomics (async-signal-safe), and `OLD_ACT` points
// exclusively to permanently leaked `Box<libc::sigaction>` allocations (`'static` lifetime).
#[cfg(feature = "profile-cpu")]
unsafe extern "C" fn sigprof_handler(
    sig: libc::c_int,
    info: *mut libc::siginfo_t,
    ucontext: *mut libc::c_void,
) {
    TICKS_TOTAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let old_ptr = OLD_ACT.load(std::sync::atomic::Ordering::Acquire);
    if !old_ptr.is_null() {
        // SAFETY: `old_ptr` was published via `Ordering::Release` from `Box::leak(Box::new(prev))`,
        // so it is non-null, properly aligned, initialized for `'static`, and never freed or mutated.
        let old = unsafe { &*old_ptr };
        let addr = old.sa_sigaction;
        if addr != 0
            && addr != libc::SIG_IGN
            && addr != libc::SIG_DFL
            && addr != sigprof_handler as *const () as usize
        {
            if old.sa_flags & libc::SA_SIGINFO != 0 {
                // SAFETY: `addr` came from `libc::sigaction` and passed the non-null, non-`SIG_IGN`,
                // non-`SIG_DFL`, and non-self check above; with `SA_SIGINFO` set, POSIX guarantees
                // `sa_sigaction` holds a valid 3-argument `extern "C"` signal handler function pointer.
                let act = unsafe {
                    std::mem::transmute::<
                        usize,
                        extern "C" fn(libc::c_int, *mut libc::siginfo_t, *mut libc::c_void),
                    >(addr)
                };
                act(sig, info, ucontext);
            } else {
                // SAFETY: `addr` came from `libc::sigaction` and passed the non-null, non-`SIG_IGN`,
                // non-`SIG_DFL`, and non-self check above; without `SA_SIGINFO`, POSIX guarantees
                // `sa_sigaction` (`sa_handler`) holds a valid 1-argument `extern "C"` handler pointer.
                let act = unsafe { std::mem::transmute::<usize, extern "C" fn(libc::c_int)>(addr) };
                act(sig);
            }
        }
    }
}

#[cfg(any(
    test,
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]
fn normalize_dir_for_lstat(dir: &std::path::Path) -> std::borrow::Cow<'_, std::path::Path> {
    #[cfg(unix)]
    {
        use std::ffi::OsStr;
        use std::os::unix::ffi::OsStrExt as _;

        let orig = dir.as_os_str().as_bytes();
        let mut bytes = orig;
        while bytes.len() > 1 {
            if bytes.ends_with(b"/") {
                bytes = &bytes[..bytes.len() - 1];
            } else if bytes.ends_with(b"/.") {
                if bytes == b"/." {
                    bytes = b"/";
                } else {
                    bytes = &bytes[..bytes.len() - 2];
                }
            } else {
                break;
            }
        }
        if bytes.len() == orig.len() {
            std::borrow::Cow::Borrowed(dir)
        } else {
            std::borrow::Cow::Owned(PathBuf::from(OsStr::from_bytes(bytes)))
        }
    }
    #[cfg(not(unix))]
    {
        let orig = dir.to_string_lossy();
        let trimmed = orig.trim_end_matches(['/', '\\']);
        if trimmed.is_empty() {
            std::borrow::Cow::Owned(PathBuf::from("/"))
        } else if trimmed.len() == orig.len() {
            std::borrow::Cow::Borrowed(dir)
        } else {
            std::borrow::Cow::Owned(PathBuf::from(trimmed))
        }
    }
}

#[cfg(any(
    test,
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]
fn ensure_non_symlink_dir(dir: &std::path::Path) -> std::io::Result<()> {
    let normalized = normalize_dir_for_lstat(dir);
    let meta = std::fs::symlink_metadata(normalized.as_ref())?;
    if meta.file_type().is_symlink() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "profile directory '{}' must not be a symlink",
                dir.display()
            ),
        ));
    }
    if !meta.is_dir() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("profile path '{}' is not a directory", dir.display()),
        ));
    }
    Ok(())
}

#[cfg(any(
    test,
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]
fn create_unique_temp_file(
    dir: &std::path::Path,
    prefix: &str,
) -> std::io::Result<(std::fs::File, PathBuf)> {
    use std::collections::hash_map::RandomState;
    use std::hash::{BuildHasher, Hasher};
    use std::sync::atomic::{AtomicU64, Ordering};

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    ensure_non_symlink_dir(dir)?;

    let pid = std::process::id();
    for attempt in 0..64u32 {
        let seq = COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut hasher = RandomState::new().build_hasher();
        hasher.write_u32(pid);
        hasher.write_u64(seq);
        hasher.write_u32(attempt);
        let nonce = hasher.finish();
        let candidate = dir.join(format!("{prefix}.{pid}.{nonce:016x}.{seq}"));

        let mut opts = std::fs::OpenOptions::new();
        opts.write(true).create_new(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt as _;
            opts.mode(0o600);
        }
        match opts.open(&candidate) {
            Ok(file) => return Ok((file, candidate)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }

    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!(
            "failed to create unique temporary file in '{}' after 64 attempts",
            dir.display()
        ),
    ))
}

#[cfg(any(test, feature = "profile-cpu", feature = "profile-heap"))]
fn write_profile_file_atomic(
    dir: &std::path::Path,
    filename: &str,
    bytes: &[u8],
) -> std::io::Result<()> {
    use std::io::Write as _;

    let prefix = format!(".{filename}.tmp");
    let (mut file, tmp_path) = create_unique_temp_file(dir, &prefix)?;
    if let Err(e) = file.write_all(bytes).and_then(|()| file.sync_all()) {
        drop(file);
        let _ = std::fs::remove_file(&tmp_path);
        return Err(e);
    }
    drop(file);

    if let Err(e) =
        ensure_non_symlink_dir(dir).and_then(|()| std::fs::rename(&tmp_path, dir.join(filename)))
    {
        let _ = std::fs::remove_file(&tmp_path);
        return Err(e);
    }
    Ok(())
}

#[cfg(any(test, feature = "profile-dhat"))]
fn finalize_dhat_temp_file(
    dir: &std::path::Path,
    dhat_tmp_path: &std::path::Path,
) -> std::io::Result<()> {
    if let Err(e) = ensure_non_symlink_dir(dir)
        .and_then(|()| std::fs::rename(dhat_tmp_path, dir.join("dhat-heap.json")))
    {
        let _ = std::fs::remove_file(dhat_tmp_path);
        return Err(e);
    }
    Ok(())
}

#[cfg(any(
    test,
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]
fn probe_dir(dir: &std::path::Path, var_name: &str) -> Result<()> {
    let normalized = normalize_dir_for_lstat(dir);
    if let Ok(meta) = std::fs::symlink_metadata(normalized.as_ref()) {
        if meta.file_type().is_symlink() {
            anyhow::bail!(
                "profile directory '{}' for {var_name} must not be a symlink",
                dir.display()
            );
        }
    }
    std::fs::create_dir_all(normalized.as_ref()).map_err(|e| {
        anyhow::anyhow!(
            "cannot create profile directory '{}' for {var_name}: {e}",
            dir.display()
        )
    })?;
    let (probe_file, probe_path) =
        create_unique_temp_file(normalized.as_ref(), ".vardict_profile_probe").map_err(|e| {
            anyhow::anyhow!(
                "profile directory '{}' for {var_name} is not writable: {e}",
                dir.display()
            )
        })?;
    drop(probe_file);
    let _ = std::fs::remove_file(&probe_path);
    Ok(())
}

impl Session {
    pub fn start() -> Result<Self> {
        let cpu_env = std::env::var_os("VARDICT_PROFILE_CPU");
        let heap_env = std::env::var_os("VARDICT_PROFILE_HEAP");
        let dhat_env = std::env::var_os("VARDICT_PROFILE_DHAT");

        let validate_dir =
            |var_name: &str, val: Option<std::ffi::OsString>| -> Result<Option<PathBuf>> {
                match val {
                    Some(v) => {
                        if v.is_empty() {
                            anyhow::bail!("{var_name} must be a non-empty directory path");
                        }
                        Ok(Some(PathBuf::from(v)))
                    }
                    None => Ok(None),
                }
            };

        let cpu_dir = validate_dir("VARDICT_PROFILE_CPU", cpu_env)?;
        #[allow(unused_mut)]
        let mut heap_dir = validate_dir("VARDICT_PROFILE_HEAP", heap_env)?;
        let dhat_dir = validate_dir("VARDICT_PROFILE_DHAT", dhat_env)?;

        #[cfg(not(feature = "profile-cpu"))]
        if cpu_dir.is_some() {
            vardictrust_core::verbose_eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_CPU is set, but the feature is not compiled in"
            );
        }

        #[cfg(not(feature = "profile-heap"))]
        if heap_dir.is_some() {
            vardictrust_core::verbose_eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_HEAP is set, but the feature is not compiled in"
            );
        }

        #[cfg(not(feature = "profile-dhat"))]
        if dhat_dir.is_some() {
            vardictrust_core::verbose_eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_DHAT is set, but the feature is not compiled in"
            );
        }

        #[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
        if heap_dir.is_some() {
            vardictrust_core::verbose_eprintln!("vardictrust: warning: VARDICT_PROFILE_HEAP requested, but global allocator does not match; disabling heap profiling");
            heap_dir = None;
        }

        // Minor 3: Probe ALL requested active directories before starting any profiler or signal handler.
        #[cfg(feature = "profile-cpu")]
        if let Some(dir) = &cpu_dir {
            probe_dir(dir, "VARDICT_PROFILE_CPU")?;
        }
        #[cfg(feature = "profile-dhat")]
        if let Some(dir) = &dhat_dir {
            probe_dir(dir, "VARDICT_PROFILE_DHAT")?;
        }
        #[cfg(feature = "profile-heap")]
        if let Some(dir) = &heap_dir {
            probe_dir(dir, "VARDICT_PROFILE_HEAP")?;
        }

        #[cfg(feature = "profile-cpu")]
        let cpu = if let Some(dir) = cpu_dir {
            let guard = pprof::ProfilerGuardBuilder::default()
                .frequency(499)
                .blocklist(&["libc", "libgcc", "pthread", "vdso", "ld-linux"])
                .build()?;

            // SAFETY: `libc::sigaction` is a C POD struct of integer/pointer fields and `sigset_t`;
            // an all-zero value is valid as an out-parameter for `libc::sigaction`.
            let mut prev: libc::sigaction = unsafe { std::mem::zeroed() };
            // SAFETY: Passing NULL `act` queries the current `SIGPROF` handler into `prev` without
            // altering signal disposition.
            let queried = unsafe { libc::sigaction(libc::SIGPROF, std::ptr::null(), &mut prev) };
            let old_act_opt = if queried == 0 {
                TICKS_TOTAL.store(0, std::sync::atomic::Ordering::Relaxed);
                let self_addr = sigprof_handler as *const () as usize;
                if is_chainable_sigprof_handler(prev.sa_sigaction) {
                    // Minor 1: Box::leak ensures the memory backing OLD_ACT is never freed during
                    // process life, eliminating use-after-free if SIGPROF fires during teardown.
                    let old_box = Box::leak(Box::new(prev));
                    OLD_ACT.store(old_box, std::sync::atomic::Ordering::Release);
                } else if prev.sa_sigaction != self_addr {
                    OLD_ACT.store(std::ptr::null_mut(), std::sync::atomic::Ordering::Release);
                }

                if prev.sa_sigaction != self_addr {
                    let mut new_act: libc::sigaction = prev;
                    new_act.sa_sigaction = self_addr;
                    new_act.sa_flags |= libc::SA_SIGINFO;
                    // SAFETY: Installing chained `SIGPROF` handler `sigprof_handler` with `SA_SIGINFO`
                    // matching its 3-argument `extern "C"` signature. If installation fails, `pprof`'s
                    // existing handler continues running.
                    let _ =
                        unsafe { libc::sigaction(libc::SIGPROF, &new_act, std::ptr::null_mut()) };
                    Some(prev)
                } else {
                    None
                }
            } else {
                None
            };

            Some(CpuSession {
                guard,
                dir,
                old_act: old_act_opt,
            })
        } else {
            None
        };

        #[cfg(feature = "profile-dhat")]
        let (dhat, dhat_paths) = if let Some(dir) = dhat_dir {
            let (tmp_file, dhat_tmp_path) = create_unique_temp_file(&dir, ".dhat-heap.tmp")
                .map_err(|e| {
                    anyhow::anyhow!(
                        "cannot create temporary DHAT profile file in '{}': {e}",
                        dir.display()
                    )
                })?;
            // Close the descriptor while keeping the owned `0600` regular file reserved on disk
            // so `dhat` writes into this exact inode on drop without symlink-follow risk (`AR-18`).
            drop(tmp_file);
            let profiler = dhat::Profiler::builder()
                .file_name(&dhat_tmp_path)
                .trim_backtraces(Some(32))
                .build();
            (Some(profiler), Some((dir, dhat_tmp_path)))
        } else {
            (None, None)
        };

        #[cfg(feature = "profile-heap")]
        let heap_dir_active = if let Some(dir) = heap_dir {
            if let Some(ctl) = jemalloc_pprof::PROF_CTL.as_ref() {
                let mut ctl = ctl.blocking_lock();
                if !ctl.activated() {
                    if let Err(e) = ctl.activate() {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile heap activate failed: {e}"
                        );
                    }
                }
            } else {
                vardictrust_core::verbose_eprintln!("vardictrust: warning: profile heap controller uninitialized; jemalloc profiling unavailable");
            }
            Some(dir)
        } else {
            None
        };

        Ok(Self {
            #[cfg(feature = "profile-cpu")]
            cpu,
            #[cfg(feature = "profile-dhat")]
            _dhat: dhat,
            #[cfg(feature = "profile-dhat")]
            dhat_paths,
            #[cfg(feature = "profile-heap")]
            heap_dir: heap_dir_active,
        })
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        #[cfg(feature = "profile-cpu")]
        if let Some(cpu) = self.cpu.take() {
            if let Some(old_act) = cpu.old_act {
                // SAFETY: `old_act` is a valid `libc::sigaction` previously returned by `libc::sigaction`;
                // restoring it upon profiler session end resets `SIGPROF` before `cpu.guard` drops.
                unsafe {
                    libc::sigaction(libc::SIGPROF, &old_act, std::ptr::null_mut());
                }
                // Note: OLD_ACT static pointer is intentionally not freed (leaked permanently),
                // ensuring concurrent or lingering SIGPROF signals never experience use-after-free.
            }

            match cpu.guard.report().build() {
                Ok(report) => {
                    let samples_kept: u64 =
                        report.data.values().copied().map(|c| c.max(0) as u64).sum();
                    // Minor 4: ticks_total is incremented on every SIGPROF signal delivery before passing to pprof.
                    // It therefore accounts for both samples kept in pprof's hash map and ticks dropped by pprof
                    // (including address blocklist matches and any internal try_write contention drops).
                    // ticks_total.max(samples_kept) guards against any potential signal-receipt skew.
                    let ticks_total = TICKS_TOTAL
                        .load(std::sync::atomic::Ordering::Relaxed)
                        .max(samples_kept);

                    let meta = format!(
                        "{{\"ticks_total\":{},\"samples_kept\":{}}}\n",
                        ticks_total, samples_kept
                    );
                    if let Err(e) =
                        write_profile_file_atomic(&cpu.dir, "cpu-meta.json", meta.as_bytes())
                    {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile write cpu-meta.json failed: {e}"
                        );
                    }

                    match report.pprof() {
                        Ok(profile) => {
                            use pprof::protos::Message;
                            let mut buf = Vec::new();
                            match profile.encode(&mut buf) {
                                Ok(()) => {
                                    if let Err(e) =
                                        write_profile_file_atomic(&cpu.dir, "cpu-profile.pb", &buf)
                                    {
                                        vardictrust_core::verbose_eprintln!(
                                            "vardictrust: warning: profile write cpu-profile.pb failed: {e}"
                                        );
                                    }
                                }
                                Err(e) => {
                                    vardictrust_core::verbose_eprintln!(
                                        "vardictrust: warning: profile pprof encode failed: {e}"
                                    );
                                }
                            }
                        }
                        Err(e) => {
                            vardictrust_core::verbose_eprintln!(
                                "vardictrust: warning: profile pprof proto generation failed: {e}"
                            );
                        }
                    }

                    use std::fmt::Write;
                    let mut lines: Vec<String> = report
                        .data
                        .iter()
                        .map(|(key, value)| {
                            let mut line = key.thread_name_or_id();
                            line.push(';');

                            for frame in key.frames.iter().rev() {
                                for symbol in frame.iter().rev() {
                                    write!(&mut line, "{};", symbol.name()).ok();
                                }
                            }

                            if line.ends_with(';') {
                                line.pop();
                            }
                            write!(&mut line, " {}", value).ok();

                            line
                        })
                        .collect();
                    lines.sort();

                    let collapsed_content = if lines.is_empty() {
                        String::new()
                    } else {
                        let mut s = lines.join("\n");
                        s.push('\n');
                        s
                    };
                    if let Err(e) = write_profile_file_atomic(
                        &cpu.dir,
                        "cpu-collapsed.txt",
                        collapsed_content.as_bytes(),
                    ) {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile write cpu-collapsed.txt failed: {e}"
                        );
                    }

                    let mut svg_bytes = Vec::new();
                    if !lines.is_empty() {
                        let mut opts = inferno::flamegraph::Options::default();
                        if let Err(e) = inferno::flamegraph::from_lines(
                            &mut opts,
                            lines.iter().map(|s| s.as_str()),
                            &mut svg_bytes,
                        ) {
                            vardictrust_core::verbose_eprintln!(
                                "vardictrust: warning: profile render flamegraph SVG failed: {e}"
                            );
                        }
                    }
                    if let Err(e) =
                        write_profile_file_atomic(&cpu.dir, "cpu-flamegraph.svg", &svg_bytes)
                    {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile write cpu-flamegraph.svg failed: {e}"
                        );
                    }
                }
                Err(e) => {
                    vardictrust_core::verbose_eprintln!(
                        "vardictrust: warning: profile CPU report failed: {e}"
                    );
                }
            }
        }

        #[cfg(feature = "profile-heap")]
        if let Some(dir) = self.heap_dir.take() {
            if let Some(ctl) = jemalloc_pprof::PROF_CTL.as_ref() {
                let mut ctl = ctl.blocking_lock();
                match ctl.dump_pprof() {
                    Ok(pprof_data) => {
                        let filename = if pprof_data.starts_with(&[0x1f, 0x8b]) {
                            "heap.pb.gz"
                        } else {
                            "heap.pb"
                        };
                        if let Err(e) = write_profile_file_atomic(&dir, filename, &pprof_data) {
                            vardictrust_core::verbose_eprintln!(
                                "vardictrust: warning: profile write heap pprof failed: {e}"
                            );
                        }
                    }
                    Err(e) => {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile dump heap pprof failed: {e}"
                        );
                    }
                }
                match ctl.dump_flamegraph() {
                    Ok(flamegraph_svg) => {
                        if let Err(e) =
                            write_profile_file_atomic(&dir, "heap-flamegraph.svg", &flamegraph_svg)
                        {
                            vardictrust_core::verbose_eprintln!(
                                "vardictrust: warning: profile write heap flamegraph failed: {e}"
                            );
                        }
                    }
                    Err(e) => {
                        vardictrust_core::verbose_eprintln!(
                            "vardictrust: warning: profile dump heap flamegraph failed: {e}"
                        );
                    }
                }
            } else {
                vardictrust_core::verbose_eprintln!(
                    "vardictrust: warning: profile heap controller uninitialized"
                );
            }
        }

        #[cfg(feature = "profile-dhat")]
        {
            drop(self._dhat.take());
            if let Some((dir, dhat_tmp_path)) = self.dhat_paths.take() {
                if let Err(e) = finalize_dhat_temp_file(&dir, &dhat_tmp_path) {
                    vardictrust_core::verbose_eprintln!(
                        "vardictrust: warning: profile write dhat-heap.json failed: {e}"
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch12_profile_temp_and_atomic_write_refuse_symlinks_and_set_0600() {
        let base = tempfile::tempdir().unwrap();
        let real_dir = base.path().join("real_prof");
        probe_dir(&real_dir, "VARDICT_PROFILE_CPU").unwrap();
        assert_eq!(std::fs::read_dir(&real_dir).unwrap().count(), 0);

        let (f1, p1) = create_unique_temp_file(&real_dir, ".dhat-heap.tmp").unwrap();
        let (f2, p2) = create_unique_temp_file(&real_dir, ".dhat-heap.tmp").unwrap();
        assert_ne!(p1, p2);
        drop(f1);
        drop(f2);

        #[cfg(unix)]
        {
            use std::os::unix::fs::{symlink, PermissionsExt as _};

            let mode1 = std::fs::symlink_metadata(&p1).unwrap().permissions().mode() & 0o777;
            assert_eq!(mode1, 0o600);

            // Symlinked profile directory (both existing target and dangling target) must be refused,
            // including trailing slashes (`path_resolution(7)` `lstat` bypass), while a real
            // subdirectory inside a symlinked parent mount (`symlink_dir/real_sub`) must succeed
            // so HPC/cloud/container symlinked `/scratch` or `$TMPDIR` mounts work cleanly.
            assert_eq!(
                normalize_dir_for_lstat(std::path::Path::new("/")).as_ref(),
                std::path::Path::new("/")
            );
            assert_eq!(
                normalize_dir_for_lstat(std::path::Path::new("///")).as_ref(),
                std::path::Path::new("/")
            );
            assert_eq!(
                normalize_dir_for_lstat(std::path::Path::new("/.")).as_ref(),
                std::path::Path::new("/")
            );
            assert_eq!(
                normalize_dir_for_lstat(std::path::Path::new(".")).as_ref(),
                std::path::Path::new(".")
            );
            assert_eq!(
                normalize_dir_for_lstat(std::path::Path::new("./")).as_ref(),
                std::path::Path::new(".")
            );
            let real_trailing_one = PathBuf::from(format!("{}/", real_dir.display()));
            let real_trailing_three = PathBuf::from(format!("{}///", real_dir.display()));
            probe_dir(&real_trailing_one, "VARDICT_PROFILE_CPU").unwrap();
            ensure_non_symlink_dir(&real_trailing_three).unwrap();

            let symlink_dir = base.path().join("symlink_prof");
            symlink(&real_dir, &symlink_dir).unwrap();
            let symlink_slash_1 = PathBuf::from(format!("{}/", symlink_dir.display()));
            let symlink_slash_3 = PathBuf::from(format!("{}///", symlink_dir.display()));
            let symlink_dot = PathBuf::from(format!("{}/.", symlink_dir.display()));

            for bad_sym_dir in [
                &symlink_dir,
                &symlink_slash_1,
                &symlink_slash_3,
                &symlink_dot,
            ] {
                assert!(
                    probe_dir(bad_sym_dir, "VARDICT_PROFILE_CPU").is_err(),
                    "probe_dir must refuse symlink/trailing-slash path '{}'",
                    bad_sym_dir.display()
                );
                assert!(
                    ensure_non_symlink_dir(bad_sym_dir).is_err(),
                    "ensure_non_symlink_dir must refuse '{}'",
                    bad_sym_dir.display()
                );
                assert!(
                    create_unique_temp_file(bad_sym_dir, ".probe").is_err(),
                    "create_unique_temp_file must refuse '{}'",
                    bad_sym_dir.display()
                );
                assert!(
                    write_profile_file_atomic(bad_sym_dir, "cpu-meta.json", b"{}\n").is_err(),
                    "write_profile_file_atomic must refuse '{}'",
                    bad_sym_dir.display()
                );
            }

            // A real subdirectory inside a symlinked parent mount (`symlink_dir/real_sub`) must succeed.
            let real_sub_under_sym = symlink_dir.join("real_sub");
            probe_dir(&real_sub_under_sym, "VARDICT_PROFILE_CPU").unwrap();
            ensure_non_symlink_dir(&real_sub_under_sym).unwrap();
            write_profile_file_atomic(&real_sub_under_sym, "cpu-meta.json", b"{}\n").unwrap();
            assert_eq!(
                std::fs::read(real_dir.join("real_sub/cpu-meta.json")).unwrap(),
                b"{}\n"
            );

            let nonexistent_target = base.path().join("never_create_me");
            let dangling_dir = base.path().join("dangling_prof");
            symlink(&nonexistent_target, &dangling_dir).unwrap();
            let dangling_slash = PathBuf::from(format!("{}/", dangling_dir.display()));
            let dangling_child = dangling_dir.join("child");
            for bad_dangling in [&dangling_dir, &dangling_slash, &dangling_child] {
                assert!(
                    probe_dir(bad_dangling, "VARDICT_PROFILE_CPU").is_err(),
                    "probe_dir must refuse dangling symlink path '{}'",
                    bad_dangling.display()
                );
            }
            assert!(
                !nonexistent_target.exists(),
                "probe_dir must not create target of dangling symlink"
            );

            // Pre-existing symlink at target output filename must be atomically replaced, never written through
            let victim = base.path().join("victim.txt");
            std::fs::write(&victim, b"precious").unwrap();
            let target_meta = real_dir.join("cpu-meta.json");
            symlink(&victim, &target_meta).unwrap();

            write_profile_file_atomic(&real_dir, "cpu-meta.json", b"{\"ticks_total\":1}\n")
                .unwrap();
            assert_eq!(std::fs::read(&victim).unwrap(), b"precious");
            let out_meta = std::fs::symlink_metadata(&target_meta).unwrap();
            assert!(!out_meta.file_type().is_symlink() && out_meta.is_file());
            assert_eq!(out_meta.permissions().mode() & 0o777, 0o600);
            assert_eq!(
                std::fs::read(&target_meta).unwrap(),
                b"{\"ticks_total\":1}\n"
            );

            // DHAT temp finalization replaces pre-existing symlink at dhat-heap.json and unlinks temp on error
            let target_dhat = real_dir.join("dhat-heap.json");
            symlink(&victim, &target_dhat).unwrap();
            std::fs::write(&p1, b"{\"dhat\":true}\n").unwrap();
            finalize_dhat_temp_file(&real_dir, &p1).unwrap();
            assert!(!p1.exists());
            assert_eq!(std::fs::read(&victim).unwrap(), b"precious");
            let dhat_meta = std::fs::symlink_metadata(&target_dhat).unwrap();
            assert!(!dhat_meta.file_type().is_symlink() && dhat_meta.is_file());
            assert_eq!(dhat_meta.permissions().mode() & 0o777, 0o600);
            assert_eq!(std::fs::read(&target_dhat).unwrap(), b"{\"dhat\":true}\n");

            // If dir is a symlink (including trailing-slash symlink) at DHAT finalization time,
            // refuse and unlink the staged temp file.
            assert!(finalize_dhat_temp_file(&symlink_dir, &p2).is_err());
            assert!(
                !p2.exists(),
                "finalize_dhat_temp_file must unlink staged temp file on error"
            );
            let (f3, p3) = create_unique_temp_file(&real_dir, ".dhat-heap.tmp").unwrap();
            drop(f3);
            assert!(finalize_dhat_temp_file(&symlink_slash_1, &p3).is_err());
            assert!(
                !p3.exists(),
                "finalize_dhat_temp_file must unlink staged temp file on trailing-slash symlink error"
            );
        }
    }

    #[cfg(feature = "profile-cpu")]
    #[test]
    fn batch12_sigprof_handler_guards_against_self_recursion_and_sentinels() {
        static DUMMY_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        extern "C" fn dummy_siginfo(
            _sig: libc::c_int,
            _info: *mut libc::siginfo_t,
            _ucontext: *mut libc::c_void,
        ) {
            DUMMY_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
        extern "C" fn dummy_simple(_sig: libc::c_int) {
            DUMMY_HITS.fetch_add(10, std::sync::atomic::Ordering::Relaxed);
        }

        let prev_ptr = OLD_ACT.load(std::sync::atomic::Ordering::Acquire);
        for bad_addr in [
            0usize,
            libc::SIG_DFL,
            libc::SIG_IGN,
            sigprof_handler as *const () as usize,
        ] {
            assert!(!is_chainable_sigprof_handler(bad_addr));
            for flags in [0, libc::SA_SIGINFO] {
                // SAFETY: `libc::sigaction` is a C POD struct valid when zero-initialized.
                let mut fake: libc::sigaction = unsafe { std::mem::zeroed() };
                fake.sa_sigaction = bad_addr;
                fake.sa_flags = flags;
                let leaked = Box::leak(Box::new(fake));
                OLD_ACT.store(leaked, std::sync::atomic::Ordering::Release);
                let before = TICKS_TOTAL.load(std::sync::atomic::Ordering::Relaxed);
                // SAFETY: Invoking `sigprof_handler` directly with null info/ucontext while
                // `OLD_ACT` holds a non-chainable sentinel/self address must not dereference
                // info/ucontext or recurse.
                unsafe {
                    sigprof_handler(libc::SIGPROF, std::ptr::null_mut(), std::ptr::null_mut());
                }
                let after = TICKS_TOTAL.load(std::sync::atomic::Ordering::Relaxed);
                assert_eq!(after, before + 1);
            }
        }

        assert!(is_chainable_sigprof_handler(
            dummy_siginfo as *const () as usize
        ));
        assert!(is_chainable_sigprof_handler(
            dummy_simple as *const () as usize
        ));
        DUMMY_HITS.store(0, std::sync::atomic::Ordering::Relaxed);
        for (addr, flags) in [
            (dummy_siginfo as *const () as usize, libc::SA_SIGINFO),
            (dummy_simple as *const () as usize, 0),
        ] {
            // SAFETY: `libc::sigaction` is a C POD struct valid when zero-initialized.
            let mut valid: libc::sigaction = unsafe { std::mem::zeroed() };
            valid.sa_sigaction = addr;
            valid.sa_flags = flags;
            let leaked = Box::leak(Box::new(valid));
            OLD_ACT.store(leaked, std::sync::atomic::Ordering::Release);
            // SAFETY: `dummy_siginfo` and `dummy_simple` match `flags` and ignore null pointers.
            unsafe {
                sigprof_handler(libc::SIGPROF, std::ptr::null_mut(), std::ptr::null_mut());
            }
        }
        assert_eq!(DUMMY_HITS.load(std::sync::atomic::Ordering::Relaxed), 11);

        OLD_ACT.store(prev_ptr, std::sync::atomic::Ordering::Release);
    }
}
