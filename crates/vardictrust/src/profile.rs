use anyhow::Result;
use std::path::PathBuf;

pub struct Session {
    #[cfg(feature = "profile-cpu")]
    cpu: Option<CpuSession>,
    #[cfg(feature = "profile-dhat")]
    _dhat: Option<dhat::Profiler>,
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

// SAFETY: sigprof_handler is an extern "C" signal handler called by the OS kernel on SIGPROF.
// TICKS_TOTAL uses atomic fetch_add (lock addq on x86_64), which is async-signal-safe.
// OLD_ACT points to a permanently-leaked Box<libc::sigaction> that is never deallocated,
// preventing use-after-free even if signals arrive concurrently or during process teardown.
#[cfg(feature = "profile-cpu")]
unsafe extern "C" fn sigprof_handler(
    sig: libc::c_int,
    info: *mut libc::siginfo_t,
    ucontext: *mut libc::c_void,
) {
    TICKS_TOTAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let old_ptr = OLD_ACT.load(std::sync::atomic::Ordering::Relaxed);
    if !old_ptr.is_null() {
        // SAFETY: old_ptr is valid for 'static and never deallocated.
        let old = &*old_ptr;
        if old.sa_flags & libc::SA_SIGINFO != 0 {
            // SAFETY: With SA_SIGINFO set, sa_sigaction holds a valid function pointer with the siginfo signature.
            let act: extern "C" fn(libc::c_int, *mut libc::siginfo_t, *mut libc::c_void) =
                std::mem::transmute(old.sa_sigaction);
            act(sig, info, ucontext);
        } else if old.sa_sigaction != 0
            && old.sa_sigaction != libc::SIG_IGN
            && old.sa_sigaction != libc::SIG_DFL
        {
            // SAFETY: Without SA_SIGINFO, sa_sigaction/sa_handler holds a valid 1-argument function pointer.
            let act: extern "C" fn(libc::c_int) = std::mem::transmute(old.sa_sigaction);
            act(sig);
        }
    }
}

#[cfg(any(
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]
fn probe_dir(dir: &std::path::Path, var_name: &str) -> Result<()> {
    std::fs::create_dir_all(dir).map_err(|e| {
        anyhow::anyhow!(
            "cannot create profile directory '{}' for {var_name}: {e}",
            dir.display()
        )
    })?;
    let probe_file = dir.join(".vardict_profile_write_probe");
    std::fs::write(&probe_file, b"").map_err(|e| {
        anyhow::anyhow!(
            "profile directory '{}' for {var_name} is not writable: {e}",
            dir.display()
        )
    })?;
    let _ = std::fs::remove_file(&probe_file);
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
            eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_CPU is set, but the feature is not compiled in"
            );
        }

        #[cfg(not(feature = "profile-heap"))]
        if heap_dir.is_some() {
            eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_HEAP is set, but the feature is not compiled in"
            );
        }

        #[cfg(not(feature = "profile-dhat"))]
        if dhat_dir.is_some() {
            eprintln!(
                "vardictrust: warning: VARDICT_PROFILE_DHAT is set, but the feature is not compiled in"
            );
        }

        #[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
        if heap_dir.is_some() {
            eprintln!("vardictrust: warning: VARDICT_PROFILE_HEAP requested, but global allocator does not match; disabling heap profiling");
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
                .blocklist(&["libc", "libgcc", "pthread", "vdso"])
                .build()?;

            // SAFETY: Initializing zeroed libc::sigaction structure for querying current handler.
            let mut old_act: libc::sigaction = unsafe { std::mem::zeroed() };
            // SAFETY: Passing NULL new_act queries the current SIGPROF handler without modifying it.
            let queried = unsafe { libc::sigaction(libc::SIGPROF, std::ptr::null(), &mut old_act) };
            let old_act_opt = if queried == 0 {
                // Minor 1: Box::leak ensures the memory backing OLD_ACT is never freed during process life,
                // completely eliminating use-after-free if SIGPROF fires concurrently during teardown.
                let old_box = Box::leak(Box::new(old_act));
                OLD_ACT.store(old_box, std::sync::atomic::Ordering::Release);
                TICKS_TOTAL.store(0, std::sync::atomic::Ordering::Relaxed);

                let mut new_act: libc::sigaction = old_act;
                new_act.sa_sigaction = sigprof_handler as *const () as usize;
                new_act.sa_flags |= libc::SA_SIGINFO;
                // SAFETY: Installing chained SIGPROF handler with SA_SIGINFO so sigprof_handler can delegate.
                // Minor 4: If sigaction chaining installation fails (e.g. unsupported platform), it is safely
                // ignored and pprof's default handler continues running.
                let _ = unsafe { libc::sigaction(libc::SIGPROF, &new_act, std::ptr::null_mut()) };
                Some(old_act)
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
        let dhat = dhat_dir.map(|dir| {
            dhat::Profiler::builder()
                .file_name(dir.join("dhat-heap.json"))
                .trim_backtraces(Some(32))
                .build()
        });

        #[cfg(feature = "profile-heap")]
        let heap_dir_active = if let Some(dir) = heap_dir {
            if let Some(ctl) = jemalloc_pprof::PROF_CTL.as_ref() {
                let mut ctl = ctl.blocking_lock();
                if !ctl.activated() {
                    if let Err(e) = ctl.activate() {
                        eprintln!("vardictrust: warning: profile heap activate failed: {e}");
                    }
                }
            } else {
                eprintln!("vardictrust: warning: profile heap controller uninitialized; jemalloc profiling unavailable");
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
                // SAFETY: Restoring pprof's original SIGPROF signal handler upon profiler session end.
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
                    if let Err(e) = std::fs::write(cpu.dir.join("cpu-meta.json"), meta.as_bytes()) {
                        eprintln!("vardictrust: warning: profile write cpu-meta.json failed: {e}");
                    }

                    match report.pprof() {
                        Ok(profile) => {
                            use pprof::protos::Message;
                            let mut buf = Vec::new();
                            match profile.encode(&mut buf) {
                                Ok(()) => {
                                    if let Err(e) =
                                        std::fs::write(cpu.dir.join("cpu-profile.pb"), buf)
                                    {
                                        eprintln!(
                                            "vardictrust: warning: profile write cpu-profile.pb failed: {e}"
                                        );
                                    }
                                }
                                Err(e) => {
                                    eprintln!(
                                        "vardictrust: warning: profile pprof encode failed: {e}"
                                    );
                                }
                            }
                        }
                        Err(e) => {
                            eprintln!(
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
                    if let Err(e) = std::fs::write(
                        cpu.dir.join("cpu-collapsed.txt"),
                        collapsed_content.as_bytes(),
                    ) {
                        eprintln!(
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
                            eprintln!(
                                "vardictrust: warning: profile render flamegraph SVG failed: {e}"
                            );
                        }
                    }
                    if let Err(e) = std::fs::write(cpu.dir.join("cpu-flamegraph.svg"), svg_bytes) {
                        eprintln!(
                            "vardictrust: warning: profile write cpu-flamegraph.svg failed: {e}"
                        );
                    }
                }
                Err(e) => {
                    eprintln!("vardictrust: warning: profile CPU report failed: {e}");
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
                        if let Err(e) = std::fs::write(dir.join(filename), &pprof_data) {
                            eprintln!("vardictrust: warning: profile write heap pprof failed: {e}");
                        }
                    }
                    Err(e) => {
                        eprintln!("vardictrust: warning: profile dump heap pprof failed: {e}");
                    }
                }
                match ctl.dump_flamegraph() {
                    Ok(flamegraph_svg) => {
                        if let Err(e) =
                            std::fs::write(dir.join("heap-flamegraph.svg"), &flamegraph_svg)
                        {
                            eprintln!(
                                "vardictrust: warning: profile write heap flamegraph failed: {e}"
                            );
                        }
                    }
                    Err(e) => {
                        eprintln!("vardictrust: warning: profile dump heap flamegraph failed: {e}");
                    }
                }
            } else {
                eprintln!("vardictrust: warning: profile heap controller uninitialized");
            }
        }
    }
}
