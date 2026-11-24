use crate::types::VardictError;
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::sync::{Arc, Condvar, Mutex};

/// Executes `worker_fn` in parallel over `items` using at most `num_threads`,
/// while streaming the results sequentially into `consumer_fn` in the exact
/// original order of `items`.
#[deprecated(since = "1.0.10", note = "use runner::run_parallel_for_bam")]
#[allow(deprecated)]
pub fn run_parallel<T, R, F, W>(
    items: Vec<T>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    let batches: Vec<Vec<T>> = items.into_iter().map(|item| vec![item]).collect();
    run_parallel_batched(batches, num_threads, worker_fn, consumer_fn)
}

/// Executes `worker_fn` in parallel over `items` using the worker stack floor
/// appropriate for `bam_spec` (64 MiB when any path is `.cram`, 16 MiB otherwise).
pub fn run_parallel_for_bam<T, R, F, W>(
    bam_spec: &str,
    items: Vec<T>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    let batches: Vec<Vec<T>> = items.into_iter().map(|item| vec![item]).collect();
    run_parallel_batched_for_bam(bam_spec, batches, num_threads, worker_fn, consumer_fn)
}

/// Weighted crate-internal counterpart to [`run_parallel_for_bam`] that accounts
/// out-of-order pending result heap bytes via `weight_fn` (`AR-8`).
pub(crate) fn run_parallel_for_bam_weighted<T, R, F, W>(
    bam_spec: &str,
    items: Vec<T>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
    weight_fn: impl Fn(&R) -> usize + Send + Sync,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    let batches: Vec<Vec<T>> = items.into_iter().map(|item| vec![item]).collect();
    run_parallel_batched_for_bam_weighted(
        bam_spec,
        batches,
        num_threads,
        worker_fn,
        consumer_fn,
        weight_fn,
    )
}

/// Worker thread stack floor (16 MiB) for Rayon tile workers (`M-DOCUMENTED-MAGIC`).
///
/// Deep CRAM container decompression (`noodles-cram` RANS/bzip2/lzma codecs) combined
/// with recursive structural-variant and indel realignment frames on 50x–150x WGS tiles
/// exceeds Rayon's default 2 MiB worker stack, causing `SIGABRT` stack overflows if lower.
pub const WORKER_THREAD_STACK_FLOOR_BYTES: usize = 16 * 1024 * 1024;

/// Minimum stack size (64 MiB) for worker threads when reading CRAM inputs
/// (`noodles-cram` RANS/LZMA/bzip2 slice decoding + deep realignment recursion).
pub const CRAM_WORKER_THREAD_STACK_FLOOR_BYTES: usize = 64 * 1024 * 1024;

/// Returns `true` if any path in a VarDict `-b` specification (`tumor.cram|normal.bam`,
/// `a.cram:b.bam`, or single path) has a `.cram` extension (case-insensitive).
#[inline]
pub(crate) fn bam_spec_has_cram(bam_spec: &str) -> bool {
    bam_spec
        .split(['|', ':'])
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .any(|part| {
            std::path::Path::new(part)
                .extension()
                .is_some_and(|ext| ext.eq_ignore_ascii_case("cram"))
        })
}

fn worker_thread_stack_bytes_for_env_with_floor(
    rust_min_stack: Option<&str>,
    floor_bytes: usize,
) -> usize {
    let requested = rust_min_stack.and_then(parse_stack_bytes).unwrap_or(0);
    requested.max(floor_bytes)
}

/// Pool `stack_size` from `RUST_MIN_STACK` with a 16 MiB BAM floor (`.cram`
/// inputs automatically enforce the 64 MiB [`CRAM_WORKER_THREAD_STACK_FLOOR_BYTES`]
/// floor via [`worker_thread_stack_bytes_for_bam`]). Accepts plain byte counts
/// (`67108864`) and case-insensitive `K`/`M`/`G` suffixes (`64M` = 67_108_864);
/// absent, empty, unparsable, or overflowing values fall back to the floor.
fn worker_thread_stack_bytes_for_env(raw: Option<&str>) -> usize {
    worker_thread_stack_bytes_for_env_with_floor(raw, WORKER_THREAD_STACK_FLOOR_BYTES)
}

/// Parse a byte count with an optional case-insensitive `K`/`M`/`G` suffix.
/// Saturating via `checked_mul`; `None` on empty/garbage/overflow so the
/// caller falls back to the floor. No panics: suffix split and digit parse
/// are both `Option`-chained (§9).
fn parse_stack_bytes(raw: &str) -> Option<usize> {
    let s = raw.trim();
    let (digits, mult) = match s.as_bytes().last() {
        Some(b'k' | b'K') => (s.get(..s.len().saturating_sub(1))?, 1024usize),
        Some(b'm' | b'M') => (s.get(..s.len().saturating_sub(1))?, 1024 * 1024),
        Some(b'g' | b'G') => (s.get(..s.len().saturating_sub(1))?, 1024 * 1024 * 1024),
        _ => (s, 1),
    };
    digits.parse::<usize>().ok()?.checked_mul(mult)
}

pub(crate) fn worker_thread_stack_bytes() -> usize {
    worker_thread_stack_bytes_for_env(std::env::var("RUST_MIN_STACK").ok().as_deref())
}

/// Effective worker stack size for a given `-b` BAM/CRAM specification.
pub(crate) fn worker_thread_stack_bytes_for_bam(bam_spec: &str) -> usize {
    let floor = if bam_spec_has_cram(bam_spec) {
        CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
    } else {
        WORKER_THREAD_STACK_FLOOR_BYTES
    };
    let env_val = std::env::var("RUST_MIN_STACK").ok();
    worker_thread_stack_bytes_for_env_with_floor(env_val.as_deref(), floor)
}

struct CurrentWorkerRegion {
    active: bool,
    chr: String,
    start: u32,
    end: u32,
    #[cfg(test)]
    override_label: Option<String>,
}

impl CurrentWorkerRegion {
    const MAX_RETAINED_CHR_CAP: usize = 256;
    const SHRUNK_CHR_CAP: usize = 64;

    /// Marks the worker region inactive and reclaims oversized contig buffers across worker reuse.
    #[inline]
    fn deactivate(&mut self) {
        self.active = false;
        if self.chr.capacity() > Self::MAX_RETAINED_CHR_CAP {
            self.chr.clear();
            self.chr.shrink_to(Self::SHRUNK_CHR_CAP);
        }
    }
}

thread_local! {
    static CURRENT_WORKER_REGION: std::cell::RefCell<CurrentWorkerRegion> = const {
        std::cell::RefCell::new(CurrentWorkerRegion {
            active: false,
            chr: String::new(),
            start: 0,
            end: 0,
            #[cfg(test)]
            override_label: None,
        })
    };
    static CURRENT_WORKER_STACK_BYTES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// RAII guard tracking the active genomic region on a worker thread without
/// per-tile `format!` heap allocation; formats `{chr}:{start}-{end}` only if a
/// worker panic is caught (`AR-16`).
pub(crate) struct WorkerRegionGuard(());

impl Drop for WorkerRegionGuard {
    fn drop(&mut self) {
        if std::thread::panicking() {
            return;
        }
        let _ = CURRENT_WORKER_REGION.try_with(|slot| {
            if let Ok(mut reg) = slot.try_borrow_mut() {
                reg.deactivate();
            }
        });
    }
}

pub(crate) fn enter_worker_region(chr: &str, start: u32, end: u32) -> WorkerRegionGuard {
    let _ = CURRENT_WORKER_REGION.try_with(|slot| {
        if let Ok(mut reg) = slot.try_borrow_mut() {
            if reg.chr != chr {
                reg.chr.clear();
                reg.chr.push_str(chr);
            }
            reg.start = start;
            reg.end = end;
            reg.active = true;
            #[cfg(test)]
            {
                reg.override_label = None;
            }
        }
    });
    WorkerRegionGuard(())
}

/// Sets the region label currently being processed on the worker thread so
/// caught panics in [`run_parallel_batched`] include the active genomic region.
#[cfg(test)]
pub(crate) fn set_current_worker_region(region: Option<String>) {
    let _ = CURRENT_WORKER_REGION.try_with(|slot| {
        if let Ok(mut reg) = slot.try_borrow_mut() {
            reg.active = false;
            reg.override_label = region;
        }
    });
}

fn clear_current_worker_region() {
    let _ = CURRENT_WORKER_REGION.try_with(|slot| {
        if let Ok(mut reg) = slot.try_borrow_mut() {
            reg.deactivate();
            #[cfg(test)]
            {
                reg.override_label = None;
            }
        }
    });
}

fn take_current_worker_region() -> Option<String> {
    CURRENT_WORKER_REGION
        .try_with(|slot| {
            let mut reg = slot.try_borrow_mut().ok()?;
            #[cfg(test)]
            if let Some(lbl) = reg.override_label.take() {
                reg.active = false;
                return Some(lbl);
            }
            if !reg.active {
                return None;
            }
            let formatted = format!("{}:{}-{}", reg.chr, reg.start, reg.end);
            reg.deactivate();
            Some(formatted)
        })
        .ok()
        .flatten()
}

fn run_with_worker_stack_bytes<F, R>(stack_bytes: usize, f: F) -> Result<R, VardictError>
where
    F: FnOnce() -> Result<R, VardictError> + Send,
    R: Send,
{
    if CURRENT_WORKER_STACK_BYTES.with(|c| c.get()) >= stack_bytes {
        return f();
    }
    std::thread::scope(|s| {
        let handle = std::thread::Builder::new()
            .name("vardictrust-worker".to_string())
            .stack_size(stack_bytes)
            .spawn_scoped(s, move || {
                CURRENT_WORKER_STACK_BYTES.with(|c| c.set(stack_bytes));
                f()
            })
            .map_err(VardictError::Io)?;
        match handle.join() {
            Ok(res) => res,
            Err(payload) => std::panic::resume_unwind(payload),
        }
    })
}

/// Executes `f` on a scoped worker thread configured with `worker_thread_stack_bytes`
/// (16 MiB floor, respecting `RUST_MIN_STACK`) so single-threaded (`-th 1`) and
/// amplicon runs receive the same stack headroom as multi-threaded Rayon workers (`R5`).
#[deprecated(since = "1.0.10", note = "use runner::run_with_worker_stack_for_bam")]
pub fn run_with_worker_stack<F, R>(f: F) -> Result<R, VardictError>
where
    F: FnOnce() -> Result<R, VardictError> + Send,
    R: Send,
{
    run_with_worker_stack_bytes(worker_thread_stack_bytes(), f)
}

/// Executes `f` on a scoped worker thread configured with `worker_thread_stack_bytes_for_bam`
/// (64 MiB floor when `bam_spec` includes `.cram`, 16 MiB otherwise, respecting `RUST_MIN_STACK`).
pub fn run_with_worker_stack_for_bam<F, R>(bam_spec: &str, f: F) -> Result<R, VardictError>
where
    F: FnOnce() -> Result<R, VardictError> + Send,
    R: Send,
{
    run_with_worker_stack_bytes(worker_thread_stack_bytes_for_bam(bam_spec), f)
}

/// Per-worker in-flight batch multiplier (`max_in_flight = num_threads.max(1) * 2`,
/// `M-DOCUMENTED-MAGIC`): keeps one active batch plus one pre-queued batch per worker
/// thread so Rayon work-stealing stays saturated without unbounded task queueing (`AR-7`).
pub(crate) const MAX_IN_FLIGHT_FACTOR: usize = 2;

/// Dispatch-window multiplier (`4 * num_threads`, `M-DOCUMENTED-MAGIC`) for public unweighted
/// runner entrypoints where `size_of::<R>()` does not account owned `Vec`/`String` heap buffers;
/// crate-internal `_weighted` callers pass `None` and bound actual heap bytes via `weight_fn`.
const UNWEIGHTED_DISPATCH_WINDOW_FACTOR: usize = 4;

/// Maximum out-of-order pending output bytes (`64 MiB`, `M-DOCUMENTED-MAGIC`) buffered in
/// the ordered consumer `BTreeMap` while an earlier slow tile batch is still running (`AR-7`).
pub(crate) const MAX_PENDING_BYTES: usize = 64 * 1024 * 1024;

/// Maximum out-of-order pending batches (`262,144`, `M-DOCUMENTED-MAGIC`): covers `> 1,000,000`
/// BED rows at `DEFAULT_MAX_BATCH_TILES = 4` so `773k`-row `mkbed()` AML/exome panels with
/// empty/near-empty batches never stall workers on batch count while a hot exon runs (`AR-7`).
pub(crate) const MAX_PENDING_BATCHES: usize = 262_144;

/// Estimated `BTreeMap` node + `(Vec<R>, usize)` bookkeeping overhead (`64 B`, `M-DOCUMENTED-MAGIC`)
/// charged per out-of-order batch stored in `pending`.
pub(crate) const PENDING_BATCH_OVERHEAD_BYTES: usize = 64;

/// Executes `worker_fn` in parallel over pre-grouped `batches` (`Vec<Vec<T>>`),
/// streaming flattened results into `consumer_fn` in exact deterministic order (`EBO-3`).
#[deprecated(since = "1.0.10", note = "use runner::run_parallel_batched_for_bam")]
pub fn run_parallel_batched<T, R, F, W>(
    batches: Vec<Vec<T>>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    run_parallel_batched_with_stack_weighted(
        batches,
        num_threads,
        worker_thread_stack_bytes(),
        false,
        Some(
            num_threads
                .max(1)
                .saturating_mul(UNWEIGHTED_DISPATCH_WINDOW_FACTOR),
        ),
        worker_fn,
        consumer_fn,
        |_| std::mem::size_of::<R>(),
    )
}

/// Executes `worker_fn` in parallel over pre-grouped `batches` (`Vec<Vec<T>>`)
/// using the worker stack floor appropriate for `bam_spec` (64 MiB for CRAM, 16 MiB otherwise).
pub fn run_parallel_batched_for_bam<T, R, F, W>(
    bam_spec: &str,
    batches: Vec<Vec<T>>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    run_parallel_batched_with_stack_weighted(
        batches,
        num_threads,
        worker_thread_stack_bytes_for_bam(bam_spec),
        true,
        Some(
            num_threads
                .max(1)
                .saturating_mul(UNWEIGHTED_DISPATCH_WINDOW_FACTOR),
        ),
        worker_fn,
        consumer_fn,
        |_| std::mem::size_of::<R>(),
    )
}

/// Weighted crate-internal counterpart to [`run_parallel_batched_for_bam`] that bounds
/// out-of-order pending heap bytes via `weight_fn` (`AR-7`, `AR-8`).
pub(crate) fn run_parallel_batched_for_bam_weighted<T, R, F, W>(
    bam_spec: &str,
    batches: Vec<Vec<T>>,
    num_threads: usize,
    worker_fn: F,
    consumer_fn: W,
    weight_fn: impl Fn(&R) -> usize + Send + Sync,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    run_parallel_batched_with_stack_weighted(
        batches,
        num_threads,
        worker_thread_stack_bytes_for_bam(bam_spec),
        true,
        None,
        worker_fn,
        consumer_fn,
        weight_fn,
    )
}

struct FlowState {
    next_expected_idx: usize,
    in_flight: usize,
    pending_bytes: usize,
    pending_batches: usize,
    shutdown: bool,
}

#[allow(clippy::too_many_arguments)]
fn run_parallel_batched_with_stack_weighted<T, R, F, W, G>(
    batches: Vec<Vec<T>>,
    num_threads: usize,
    stack_bytes: usize,
    enforce_serial_stack: bool,
    max_window_batches: Option<usize>,
    worker_fn: F,
    mut consumer_fn: W,
    weight_fn: G,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
    G: Fn(&R) -> usize + Send + Sync,
{
    if batches.is_empty() {
        return Ok(());
    }
    if num_threads == 1 {
        if !enforce_serial_stack || CURRENT_WORKER_STACK_BYTES.with(|c| c.get()) >= stack_bytes {
            for batch in batches {
                for item in batch {
                    consumer_fn(worker_fn(item)?)?;
                }
            }
        } else {
            for batch in batches {
                for item in batch {
                    let val = run_with_worker_stack_bytes(stack_bytes, || worker_fn(item))?;
                    consumer_fn(val)?;
                }
            }
        }
        return Ok(());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .stack_size(stack_bytes)
        .build()
        .map_err(|e| {
            VardictError::Io(std::io::Error::other(format!(
                "Failed to build rayon thread pool: {}",
                e
            )))
        })?;

    let max_in_flight = num_threads.max(1).saturating_mul(MAX_IN_FLIGHT_FACTOR);
    let (tx, rx) = mpsc::channel::<(usize, Result<Vec<R>, VardictError>)>();

    type SharedState = (Mutex<FlowState>, Condvar, std::sync::atomic::AtomicBool);
    let shared_state: Arc<SharedState> = Arc::new((
        Mutex::new(FlowState {
            next_expected_idx: 0,
            in_flight: 0,
            pending_bytes: 0,
            pending_batches: 0,
            shutdown: false,
        }),
        Condvar::new(),
        std::sync::atomic::AtomicBool::new(false),
    ));

    struct ShutdownGuard {
        state: Arc<SharedState>,
    }
    impl Drop for ShutdownGuard {
        fn drop(&mut self) {
            self.state
                .2
                .store(true, std::sync::atomic::Ordering::Release);
            let mut st = self.state.0.lock().unwrap_or_else(|e| e.into_inner());
            st.shutdown = true;
            self.state.1.notify_all();
        }
    }

    let batches_len = batches.len();
    let worker_fn = Arc::new(worker_fn);

    std::thread::scope(|s| {
        let _guard = ShutdownGuard {
            state: Arc::clone(&shared_state),
        };

        let shared_clone = Arc::clone(&shared_state);
        let pool_clone = &pool;
        let worker_fn_clone = Arc::clone(&worker_fn);

        // Feeder thread: dispatches work whenever a worker slot is available, any optional
        // unweighted batch-window bound (`4 * num_threads`) holds, and either the batch is
        // next-expected or pending out-of-order output is below 64 MiB / 262k batches.
        s.spawn(move || {
            let (lock, cvar, _) = &*shared_clone;
            for (batch_idx, batch) in batches.into_iter().enumerate() {
                let mut st = match lock.lock() {
                    Ok(g) => g,
                    Err(p) => p.into_inner(),
                };
                while !st.shutdown
                    && !(st.in_flight < max_in_flight
                        && max_window_batches
                            .is_none_or(|win| batch_idx < st.next_expected_idx.saturating_add(win))
                        && (batch_idx == st.next_expected_idx
                            || (st.pending_bytes < MAX_PENDING_BYTES
                                && st.pending_batches < MAX_PENDING_BATCHES)))
                {
                    st = match cvar.wait(st) {
                        Ok(g) => g,
                        Err(p) => p.into_inner(),
                    };
                }
                if st.shutdown {
                    break;
                }
                st.in_flight = st.in_flight.saturating_add(1);
                drop(st);

                let tx_clone = tx.clone();
                let w_fn = Arc::clone(&worker_fn_clone);
                let state_for_worker = Arc::clone(&shared_clone);
                pool_clone.spawn(move || {
                    let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        let mut out = Vec::with_capacity(batch.len());
                        for (item_idx, item) in batch.into_iter().enumerate() {
                            if item_idx > 0
                                && state_for_worker
                                    .2
                                    .load(std::sync::atomic::Ordering::Relaxed)
                            {
                                break;
                            }
                            out.push(w_fn(item)?);
                        }
                        clear_current_worker_region();
                        Ok(out)
                    }))
                    .unwrap_or_else(|_| {
                        let msg = match take_current_worker_region() {
                            Some(reg) => format!(
                                "Worker thread panicked on batch {} (region {})",
                                batch_idx, reg
                            ),
                            None => format!("Worker thread panicked on batch {}", batch_idx),
                        };
                        Err(VardictError::Io(std::io::Error::other(msg)))
                    });
                    let _ = tx_clone.send((batch_idx, res));
                });
            }
        });

        // Consumer loop on main thread
        let mut next_expected_idx = 0usize;
        let mut pending: BTreeMap<usize, (Vec<R>, usize)> = BTreeMap::new();
        let (lock, cvar, _) = &*shared_state;

        for (batch_idx, res) in rx {
            let batch_vals = res?;
            if batch_idx == next_expected_idx {
                next_expected_idx = next_expected_idx.saturating_add(1);
                {
                    let mut st = match lock.lock() {
                        Ok(g) => g,
                        Err(p) => p.into_inner(),
                    };
                    st.in_flight = st.in_flight.saturating_sub(1);
                    st.next_expected_idx = next_expected_idx;
                }
                cvar.notify_all();

                for val in batch_vals {
                    consumer_fn(val)?;
                }

                while let Some((vals, stored_bytes)) = pending.remove(&next_expected_idx) {
                    next_expected_idx = next_expected_idx.saturating_add(1);
                    for val in vals {
                        consumer_fn(val)?;
                    }
                    {
                        let mut st = match lock.lock() {
                            Ok(g) => g,
                            Err(p) => p.into_inner(),
                        };
                        st.next_expected_idx = next_expected_idx;
                        st.pending_bytes = st.pending_bytes.saturating_sub(stored_bytes);
                        st.pending_batches = st.pending_batches.saturating_sub(1);
                    }
                    cvar.notify_all();
                }
            } else {
                let batch_bytes = PENDING_BATCH_OVERHEAD_BYTES
                    .saturating_add(
                        batch_vals
                            .capacity()
                            .saturating_mul(std::mem::size_of::<R>()),
                    )
                    .saturating_add(
                        batch_vals
                            .iter()
                            .map(&weight_fn)
                            .fold(0usize, usize::saturating_add),
                    );
                pending.insert(batch_idx, (batch_vals, batch_bytes));
                {
                    let mut st = match lock.lock() {
                        Ok(g) => g,
                        Err(p) => p.into_inner(),
                    };
                    st.in_flight = st.in_flight.saturating_sub(1);
                    st.pending_bytes = st.pending_bytes.saturating_add(batch_bytes);
                    st.pending_batches = st.pending_batches.saturating_add(1);
                }
                cvar.notify_all();
            }

            if next_expected_idx == batches_len {
                break;
            }
        }

        if next_expected_idx < batches_len {
            return Err(VardictError::Io(std::io::Error::other(format!(
                "Tile batch {} failed or panicked; pipeline aborted early",
                next_expected_idx
            ))));
        }

        Ok(())
    })
}

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;

    #[test]
    fn test_ordered_output() {
        let items: Vec<usize> = (0..100).collect();
        let mut output = Vec::new();

        run_parallel(
            items,
            4,
            |item| {
                if item == 0 {
                    std::thread::sleep(std::time::Duration::from_millis(50));
                } else if item % 2 == 0 {
                    std::thread::sleep(std::time::Duration::from_millis(1));
                }
                Ok(item * 2)
            },
            |val| {
                output.push(val);
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(output.len(), 100);
        for (i, &val) in output.iter().enumerate() {
            assert_eq!(val, i * 2, "Output at index {} is out of order!", i);
        }
    }

    #[test]
    fn test_early_abort_no_hang() {
        let (tx, rx) = mpsc::channel();

        std::thread::spawn(move || {
            let items: Vec<usize> = (0..100).collect();
            let res = run_parallel(
                items,
                4,
                |item| {
                    std::thread::sleep(std::time::Duration::from_millis(1));
                    Ok(item)
                },
                |val| {
                    if val == 3 {
                        return Err(VardictError::Io(std::io::Error::other("early abort")));
                    }
                    Ok(())
                },
            );
            tx.send(res.is_err()).unwrap();
        });

        // The thread should return promptly after hitting 3 and not hang waiting for 100
        let result = rx.recv_timeout(std::time::Duration::from_secs(10));
        assert!(result.is_ok(), "Test timed out, feeder thread likely hung");
        assert!(result.unwrap(), "Expected an error from early abort");
    }

    #[test]
    fn test_worker_panic_returns_err() {
        let (tx, rx) = mpsc::channel();

        std::thread::spawn(move || {
            let items: Vec<usize> = (0..20).collect();
            let res = run_parallel(
                items,
                4,
                |item| {
                    if item == 7 {
                        panic!("Simulated worker panic");
                    }
                    std::thread::sleep(std::time::Duration::from_millis(1));
                    Ok(item)
                },
                |_| Ok(()),
            );
            tx.send(res).unwrap();
        });

        let result = rx.recv_timeout(std::time::Duration::from_secs(10));
        assert!(result.is_ok(), "Test timed out, feeder thread likely hung");
        let run_result = result.unwrap();
        assert!(run_result.is_err(), "Expected an error from worker panic");
    }

    #[test]
    fn test_batched_ordered_output_and_early_abort() {
        let batches: Vec<Vec<usize>> = (0..20).map(|b| (b * 8..(b + 1) * 8).collect()).collect();
        let mut results = Vec::new();
        run_parallel_batched(
            batches.clone(),
            4,
            |item| {
                if item == 0 {
                    std::thread::sleep(std::time::Duration::from_millis(30));
                }
                Ok(item * 2)
            },
            |val| {
                results.push(val);
                Ok(())
            },
        )
        .unwrap();
        let expected: Vec<usize> = (0..160).map(|i| i * 2).collect();
        assert_eq!(results, expected);

        // Early abort mid-batch must return Err without hanging
        let mut seen = Vec::new();
        let res = run_parallel_batched(batches, 4, Ok, |val| {
            seen.push(val);
            if val == 5 {
                Err(VardictError::Other("abort mid-batch".to_string()))
            } else {
                Ok(())
            }
        });
        assert!(res.is_err());
        assert_eq!(seen, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn test_run_parallel_single_thread_synchronous_and_error() {
        let batches = vec![vec![1usize, 2, 3], vec![4, 5]];
        let caller_tid = std::thread::current().id();
        let mut results = Vec::new();
        run_parallel_batched(
            batches.clone(),
            1,
            move |item| {
                assert_eq!(
                    std::thread::current().id(),
                    caller_tid,
                    "num_threads == 1 must execute synchronously on the calling thread"
                );
                Ok(item * 10)
            },
            |val| {
                results.push(val);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(results, vec![10, 20, 30, 40, 50]);

        let err_res = run_parallel_batched(
            batches.clone(),
            1,
            |item| {
                if item == 3 {
                    Err(VardictError::Other("boom".to_string()))
                } else {
                    Ok(item)
                }
            },
            |_| Ok(()),
        );
        assert!(matches!(err_res, Err(VardictError::Other(ref m)) if m == "boom"));

        let panic_res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            run_parallel_batched(
                batches,
                1,
                |item| {
                    if item == 2 {
                        panic!("worker panic");
                    }
                    Ok(item)
                },
                |_| Ok(()),
            )
        }));
        assert!(
            panic_res.is_err(),
            "single-thread caller panics must unwind naturally"
        );
    }

    #[test]
    fn worker_stack_env_floor_pin() {
        // Absent or unparsable input falls back to the 16 MiB floor.
        assert_eq!(worker_thread_stack_bytes_for_env(None), 16 * 1024 * 1024);
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("bogus")),
            16 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("")),
            16 * 1024 * 1024
        );
        // Values below the floor clamp up; values above pass through.
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("4194304")),
            16 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("67108864")),
            67108864
        );
        // Suffixed forms (case-insensitive) parse to the same byte counts.
        assert_eq!(worker_thread_stack_bytes_for_env(Some("64M")), 67108864);
        assert_eq!(worker_thread_stack_bytes_for_env(Some("64m")), 67108864);
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("16384K")),
            16 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("1G")),
            1024 * 1024 * 1024
        );
        // Suffix-only, trailing garbage, and overflow all fall back to the floor.
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("M")),
            16 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("64MB")),
            16 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env(Some("999999999999999999999G")),
            16 * 1024 * 1024
        );
        // Live reading never drops below the floor whatever the environment holds.
        assert!(worker_thread_stack_bytes() >= 16 * 1024 * 1024);
    }

    #[test]
    fn worker_pool_supports_deep_stack_frames_for_cram_and_sv_recursion() {
        #[inline(never)]
        fn recurse_stack(depth: usize, seed: u8) -> u64 {
            let mut frame = [seed; 512 * 1024];
            std::hint::black_box(&mut frame);
            if depth == 0 {
                frame[0] as u64 + frame[frame.len() - 1] as u64
            } else {
                frame[0] as u64 + recurse_stack(depth - 1, seed.wrapping_add(1))
            }
        }

        // 8 frames * 512 KiB = 4 MiB stack (> default 2 MiB Rayon worker stack, <= 16 MiB)
        let mut out = Vec::new();
        run_parallel(
            vec![1u8, 2u8],
            2,
            |seed| Ok(recurse_stack(8, seed)),
            |val| {
                out.push(val);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(out.len(), 2);

        // Even when invoked from a tiny 512 KiB caller thread, `run_with_worker_stack`
        // provides the full >= 16 MiB worker stack so `-th 1` and amplicon mode never
        // overflow the caller's stack.
        let val = std::thread::Builder::new()
            .stack_size(512 * 1024)
            .spawn(|| run_with_worker_stack(|| Ok(recurse_stack(8, 3))))
            .unwrap()
            .join()
            .unwrap()
            .unwrap();
        assert!(val > 0);

        // Worker panic includes region context when set via `set_current_worker_region`.
        let err = run_parallel(
            vec![1u8, 2u8],
            2,
            |seed| {
                set_current_worker_region(Some(format!("chr1:{seed}-100")));
                if seed == 2 {
                    panic!("simulated region panic");
                }
                Ok(seed)
            },
            |_| Ok(()),
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("region chr1:2-100"),
            "unexpected error: {err}"
        );

        // RAII `enter_worker_region` clears active state on normal exit and formats lazily on panic.
        {
            let _g = enter_worker_region(&"chr_long_".repeat(40), 10, 20);
        }
        assert_eq!(take_current_worker_region(), None);
        let guard_err = run_parallel(
            vec![1u32, 2u32],
            2,
            |seed| {
                let _g = enter_worker_region("chr7", seed, 500);
                if seed == 2 {
                    panic!("simulated guard panic");
                }
                Ok(seed)
            },
            |_| Ok(()),
        )
        .unwrap_err();
        assert!(
            guard_err.to_string().contains("region chr7:2-500"),
            "unexpected guard panic error: {guard_err}"
        );
    }

    #[test]
    fn cram_worker_stack_floor_and_bam_spec_detection() {
        assert!(bam_spec_has_cram("sample.cram"));
        assert!(bam_spec_has_cram("tumor.CRAM|normal.bam"));
        assert!(bam_spec_has_cram("a.bam:b.cram"));
        assert!(!bam_spec_has_cram("sample.bam"));
        assert!(!bam_spec_has_cram(""));

        assert_eq!(
            worker_thread_stack_bytes_for_env_with_floor(
                None,
                CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
            ),
            64 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env_with_floor(
                Some("16777216"),
                CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
            ),
            64 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env_with_floor(
                Some("134217728"),
                CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
            ),
            128 * 1024 * 1024
        );
        assert_eq!(
            worker_thread_stack_bytes_for_env_with_floor(None, WORKER_THREAD_STACK_FLOOR_BYTES),
            16 * 1024 * 1024
        );

        assert!(
            worker_thread_stack_bytes_for_bam("sample.cram")
                >= CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
        );
        assert!(worker_thread_stack_bytes_for_bam("sample.bam") >= WORKER_THREAD_STACK_FLOOR_BYTES);

        // Serial (num_threads == 1) `run_parallel_batched_for_bam` without an outer wrapper
        // must still execute `worker_fn` with the 64 MiB CRAM stack floor while keeping
        // `consumer_fn` on the calling thread.
        let caller_tid = std::thread::current().id();
        let mut consumed = Vec::new();
        run_parallel_batched_for_bam(
            "sample.cram",
            vec![vec![10usize, 20]],
            1,
            move |x| {
                assert!(
                    CURRENT_WORKER_STACK_BYTES.with(|c| c.get())
                        >= CRAM_WORKER_THREAD_STACK_FLOOR_BYTES
                );
                assert_ne!(std::thread::current().id(), caller_tid);
                Ok(x + 1)
            },
            |val| {
                assert_eq!(std::thread::current().id(), caller_tid);
                consumed.push(val);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(consumed, vec![11, 21]);
    }

    #[test]
    fn test_head_of_line_slow_batch_does_not_starve_workers() {
        let batches: Vec<Vec<usize>> = (0..200).map(|i| vec![i]).collect();
        let progress = Arc::new((Mutex::new(0usize), Condvar::new()));
        let progress_for_workers = Arc::clone(&progress);
        let batch0_obs = Arc::new(Mutex::new((0usize, true)));
        let batch0_obs_worker = Arc::clone(&batch0_obs);
        let mut consumed = Vec::with_capacity(200);

        run_parallel_batched_for_bam_weighted(
            "panel.bam",
            batches,
            4,
            move |idx| {
                if idx == 0 {
                    let (lock, cvar) = &*progress_for_workers;
                    let guard = lock.lock().unwrap();
                    let (guard, timeout_res) = cvar
                        .wait_timeout_while(guard, std::time::Duration::from_secs(3), |done| {
                            *done < 100
                        })
                        .unwrap();
                    *batch0_obs_worker.lock().unwrap() = (*guard, timeout_res.timed_out());
                    Ok(0usize)
                } else {
                    let (lock, cvar) = &*progress_for_workers;
                    let mut done = lock.lock().unwrap();
                    *done += 1;
                    if *done >= 100 {
                        cvar.notify_all();
                    }
                    Ok(idx)
                }
            },
            |val| {
                consumed.push(val);
                Ok(())
            },
            |_| std::mem::size_of::<usize>(),
        )
        .unwrap();

        let (done_before_batch0_exit, timed_out) = *batch0_obs.lock().unwrap();
        assert!(
            !timed_out,
            "batch 0 timed out waiting for non-blocked workers (done = {done_before_batch0_exit})"
        );
        assert!(
            done_before_batch0_exit >= 100,
            "expected >= 100 batches completed before batch 0 exited, got {done_before_batch0_exit}"
        );
        let expected: Vec<usize> = (0..200).collect();
        assert_eq!(consumed, expected);
    }

    #[test]
    fn test_pending_byte_cap_throttles_feeder() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        const FOUR_MIB: usize = 4 * 1024 * 1024;
        const NUM_THREADS: usize = 4;
        const TOTAL_BATCHES: usize = 40;

        struct WeightedItem {
            idx: usize,
            bytes: usize,
            live_bytes: Arc<AtomicUsize>,
        }
        impl Drop for WeightedItem {
            fn drop(&mut self) {
                self.live_bytes.fetch_sub(self.bytes, Ordering::SeqCst);
            }
        }

        let live_bytes = Arc::new(AtomicUsize::new(0));
        let peak_bytes = Arc::new(AtomicUsize::new(0));
        let nonzero_done = Arc::new((Mutex::new(0usize), Condvar::new()));

        let batches: Vec<Vec<usize>> = (0..TOTAL_BATCHES).map(|i| vec![i]).collect();
        let live_w = Arc::clone(&live_bytes);
        let peak_w = Arc::clone(&peak_bytes);
        let done_w = Arc::clone(&nonzero_done);
        let mut drained_order = Vec::with_capacity(TOTAL_BATCHES);

        run_parallel_batched_for_bam_weighted(
            "panel.bam",
            batches,
            NUM_THREADS,
            move |idx| {
                if idx == 0 {
                    let (lock, cvar) = &*done_w;
                    let guard = lock.lock().unwrap();
                    let (_guard, timeout_res) = cvar
                        .wait_timeout_while(
                            guard,
                            std::time::Duration::from_secs(3),
                            |done| *done < 16,
                        )
                        .unwrap();
                    assert!(
                        !timeout_res.timed_out(),
                        "batch 0 timed out waiting for 16 x 4 MiB batches to reach MAX_PENDING_BYTES"
                    );
                    std::thread::sleep(std::time::Duration::from_millis(30));
                }

                let prev = live_w.fetch_add(FOUR_MIB, Ordering::SeqCst);
                peak_w.fetch_max(prev + FOUR_MIB, Ordering::SeqCst);
                let item = WeightedItem {
                    idx,
                    bytes: FOUR_MIB,
                    live_bytes: Arc::clone(&live_w),
                };

                if idx > 0 {
                    let (lock, cvar) = &*done_w;
                    let mut done = lock.lock().unwrap();
                    *done += 1;
                    if *done >= 16 {
                        cvar.notify_all();
                    }
                }
                Ok(item)
            },
            |item| {
                drained_order.push(item.idx);
                Ok(())
            },
            |item| item.bytes,
        )
        .unwrap();

        let max_in_flight = NUM_THREADS * MAX_IN_FLIGHT_FACTOR;
        let max_allowed_bytes = MAX_PENDING_BYTES + max_in_flight * FOUR_MIB;
        let observed_peak = peak_bytes.load(Ordering::SeqCst);
        assert!(
            observed_peak <= max_allowed_bytes,
            "peak in_flight + pending bytes ({observed_peak}) exceeded cap ({max_allowed_bytes})"
        );
        assert_eq!(live_bytes.load(Ordering::SeqCst), 0);
        let expected_order: Vec<usize> = (0..TOTAL_BATCHES).collect();
        assert_eq!(drained_order, expected_order);

        // Early consumer Err(BrokenPipe) mid-batch without hang or usize underflow.
        let (done_tx, done_rx) = mpsc::channel();
        let live_abort = Arc::new(AtomicUsize::new(0));
        let live_abort_check = Arc::clone(&live_abort);
        std::thread::spawn(move || {
            let batches: Vec<Vec<usize>> =
                (0..24).map(|b| (b * 4..(b + 1) * 4).collect()).collect();
            let live_w = Arc::clone(&live_abort);
            let mut seen = Vec::new();
            let res = run_parallel_batched_for_bam_weighted(
                "panel.bam",
                batches,
                NUM_THREADS,
                move |idx| {
                    if idx == 0 {
                        std::thread::sleep(std::time::Duration::from_millis(20));
                    }
                    live_w.fetch_add(FOUR_MIB, Ordering::SeqCst);
                    Ok(WeightedItem {
                        idx,
                        bytes: FOUR_MIB,
                        live_bytes: Arc::clone(&live_w),
                    })
                },
                |item| {
                    seen.push(item.idx);
                    if item.idx == 6 {
                        return Err(VardictError::Io(std::io::Error::from(
                            std::io::ErrorKind::BrokenPipe,
                        )));
                    }
                    Ok(())
                },
                |item| item.bytes,
            );
            done_tx.send((res, seen)).unwrap();
        });

        let (abort_res, seen) = done_rx
            .recv_timeout(std::time::Duration::from_secs(5))
            .expect("weighted BrokenPipe mid-batch abort timed out");
        let err = abort_res.expect_err("expected BrokenPipe error");
        assert!(err.is_broken_pipe(), "expected BrokenPipe, got {err:?}");
        assert_eq!(seen, vec![0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(live_abort_check.load(Ordering::SeqCst), 0);

        // Unbatched `run_parallel_for_bam_weighted` also streams in exact 0..N order.
        let mut unbatched = Vec::new();
        run_parallel_for_bam_weighted(
            "panel.bam",
            (0..12usize).collect(),
            NUM_THREADS,
            |i| Ok(format!("row_{i}")),
            |s| {
                unbatched.push(s);
                Ok(())
            },
            String::capacity,
        )
        .unwrap();
        let expected_rows: Vec<String> = (0..12).map(|i| format!("row_{i}")).collect();
        assert_eq!(unbatched, expected_rows);
    }

    #[test]
    fn test_unweighted_runner_preserves_four_x_threads_dispatch_window() {
        const NUM_THREADS: usize = 2;
        const TOTAL_BATCHES: usize = 40;
        const EXPECTED_WINDOW: usize = NUM_THREADS * 4;

        for use_for_bam in [true, false] {
            let batches: Vec<Vec<usize>> = (0..TOTAL_BATCHES).map(|i| vec![i]).collect();
            // (downstream_completed_while_batch0_held, batch0_holding)
            let state = Arc::new((Mutex::new((0usize, true)), Condvar::new()));
            let state_w = Arc::clone(&state);
            let observed_downstream = Arc::new(Mutex::new(0usize));
            let observed_w = Arc::clone(&observed_downstream);
            let mut drained = Vec::with_capacity(TOTAL_BATCHES);

            let worker = move |idx: usize| -> Result<Vec<u8>, VardictError> {
                if idx == 0 {
                    let (lock, cvar) = &*state_w;
                    let guard = lock.lock().unwrap();
                    let (mut guard, _timeout) = cvar
                        .wait_timeout_while(
                            guard,
                            std::time::Duration::from_millis(200),
                            |(done, _holding)| *done < EXPECTED_WINDOW,
                        )
                        .unwrap();
                    *observed_w.lock().unwrap() = guard.0;
                    guard.1 = false;
                    Ok(vec![0u8; 1024])
                } else {
                    let (lock, cvar) = &*state_w;
                    let mut guard = lock.lock().unwrap();
                    if guard.1 {
                        guard.0 += 1;
                        cvar.notify_all();
                    }
                    Ok(vec![(idx & 0xff) as u8; 1024])
                }
            };

            if use_for_bam {
                run_parallel_batched_for_bam("panel.bam", batches, NUM_THREADS, worker, |buf| {
                    drained.push(buf[0]);
                    Ok(())
                })
                .unwrap();
            } else {
                run_parallel_batched(batches, NUM_THREADS, worker, |buf| {
                    drained.push(buf[0]);
                    Ok(())
                })
                .unwrap();
            }

            let downstream_while_held = *observed_downstream.lock().unwrap();
            let dispatched_ahead = 1 + downstream_while_held;
            assert!(
                dispatched_ahead <= EXPECTED_WINDOW,
                "unweighted runner (for_bam={use_for_bam}) dispatched {dispatched_ahead} batches \
                 ({downstream_while_held} downstream while batch 0 held), exceeding 4 * num_threads = {EXPECTED_WINDOW}"
            );
            let expected_drained: Vec<u8> = (0..TOTAL_BATCHES).map(|i| (i & 0xff) as u8).collect();
            assert_eq!(drained, expected_drained);
        }
    }
}
