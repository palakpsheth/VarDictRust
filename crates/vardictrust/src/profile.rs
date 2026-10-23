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

#[cfg(feature = "profile-cpu")]
unsafe extern "C" fn sigprof_handler(
    sig: libc::c_int,
    info: *mut libc::siginfo_t,
    ucontext: *mut libc::c_void,
) {
    TICKS_TOTAL.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let old_ptr = OLD_ACT.load(std::sync::atomic::Ordering::Relaxed);
    if !old_ptr.is_null() {
        let old = &*old_ptr;
        if old.sa_flags & libc::SA_SIGINFO != 0 {
            let act: extern "C" fn(libc::c_int, *mut libc::siginfo_t, *mut libc::c_void) =
                std::mem::transmute(old.sa_sigaction);
            act(sig, info, ucontext);
        } else if old.sa_sigaction != 0
            && old.sa_sigaction != libc::SIG_IGN
            && old.sa_sigaction != libc::SIG_DFL
        {
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
fn probe_dir(dir: &std::path::Path) -> Result<()> {
    std::fs::create_dir_all(dir)?;
    let probe_file = dir.join(".vardict_profile_write_probe");
    std::fs::write(&probe_file, b"")?;
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
            eprintln!("vardictrust: warning: VARDICT_PROFILE_CPU is set, but the feature is not compiled in");
        }

        #[cfg(not(feature = "profile-heap"))]
        if heap_dir.is_some() {
            eprintln!("vardictrust: warning: VARDICT_PROFILE_HEAP is set, but the feature is not compiled in");
        }

        #[cfg(not(feature = "profile-dhat"))]
        if dhat_dir.is_some() {
            eprintln!("vardictrust: warning: VARDICT_PROFILE_DHAT is set, but the feature is not compiled in");
        }

        #[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
        if heap_dir.is_some() {
            eprintln!("vardictrust: warning: VARDICT_PROFILE_HEAP requested, but global allocator does not match; disabling heap profiling");
            heap_dir = None;
        }

        #[cfg(feature = "profile-cpu")]
        let cpu = if let Some(dir) = cpu_dir {
            probe_dir(&dir)?;
            let guard = pprof::ProfilerGuardBuilder::default()
                .frequency(499)
                .blocklist(&["libc", "libgcc", "pthread", "vdso"])
                .build()?;

            let mut old_act: libc::sigaction = unsafe { std::mem::zeroed() };
            let old_act_opt =
                if unsafe { libc::sigaction(libc::SIGPROF, std::ptr::null(), &mut old_act) } == 0 {
                    let old_box = Box::into_raw(Box::new(old_act));
                    OLD_ACT.store(old_box, std::sync::atomic::Ordering::Release);
                    TICKS_TOTAL.store(0, std::sync::atomic::Ordering::Relaxed);

                    let mut new_act: libc::sigaction = old_act;
                    new_act.sa_sigaction = sigprof_handler as *const () as usize;
                    new_act.sa_flags |= libc::SA_SIGINFO;
                    unsafe {
                        libc::sigaction(libc::SIGPROF, &new_act, std::ptr::null_mut());
                    }
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
        let dhat = if let Some(dir) = dhat_dir {
            probe_dir(&dir)?;
            Some(
                dhat::Profiler::builder()
                    .file_name(dir.join("dhat-heap.json"))
                    .trim_backtraces(Some(32))
                    .build(),
            )
        } else {
            None
        };

        #[cfg(feature = "profile-heap")]
        let heap_dir_active = if let Some(dir) = heap_dir {
            probe_dir(&dir)?;
            if let Some(ctl) = jemalloc_pprof::PROF_CTL.as_ref() {
                let mut ctl = ctl.blocking_lock();
                if !ctl.activated() {
                    let _ = ctl.activate();
                }
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
                unsafe {
                    libc::sigaction(libc::SIGPROF, &old_act, std::ptr::null_mut());
                }
                let old_box =
                    OLD_ACT.swap(std::ptr::null_mut(), std::sync::atomic::Ordering::AcqRel);
                if !old_box.is_null() {
                    drop(unsafe { Box::from_raw(old_box) });
                }
            }

            if let Ok(report) = cpu.guard.report().build() {
                let samples_kept: u64 =
                    report.data.values().copied().map(|c| c.max(0) as u64).sum();
                let ticks_total = TICKS_TOTAL
                    .load(std::sync::atomic::Ordering::Relaxed)
                    .max(samples_kept);

                let meta = format!(
                    "{{\"ticks_total\":{},\"samples_kept\":{}}}\n",
                    ticks_total, samples_kept
                );
                let _ = std::fs::write(cpu.dir.join("cpu-meta.json"), meta.as_bytes());

                if let Ok(profile) = report.pprof() {
                    use pprof::protos::Message;
                    let mut buf = Vec::new();
                    if profile.encode(&mut buf).is_ok() {
                        let _ = std::fs::write(cpu.dir.join("cpu-profile.pb"), buf);
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
                let _ = std::fs::write(
                    cpu.dir.join("cpu-collapsed.txt"),
                    collapsed_content.as_bytes(),
                );

                let mut svg_bytes = Vec::new();
                if !lines.is_empty() {
                    let mut opts = inferno::flamegraph::Options::default();
                    let _ = inferno::flamegraph::from_lines(
                        &mut opts,
                        lines.iter().map(|s| s.as_str()),
                        &mut svg_bytes,
                    );
                }
                let _ = std::fs::write(cpu.dir.join("cpu-flamegraph.svg"), svg_bytes);
            }
        }

        #[cfg(feature = "profile-heap")]
        if let Some(dir) = self.heap_dir.take() {
            if let Some(ctl) = jemalloc_pprof::PROF_CTL.as_ref() {
                let mut ctl = ctl.blocking_lock();
                if let Ok(pprof_data) = ctl.dump_pprof() {
                    let filename = if pprof_data.starts_with(&[0x1f, 0x8b]) {
                        "heap.pb.gz"
                    } else {
                        "heap.pb"
                    };
                    let _ = std::fs::write(dir.join(filename), &pprof_data);
                }
                if let Ok(flamegraph_svg) = ctl.dump_flamegraph() {
                    let _ = std::fs::write(dir.join("heap-flamegraph.svg"), &flamegraph_svg);
                }
            }
        }
    }
}
