use crate::types::VardictError;
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::sync::{Arc, Condvar, Mutex};

/// Executes `worker_fn` in parallel over `items` using at most `num_threads`,
/// while streaming the results sequentially into `consumer_fn` in the exact
/// original order of `items`.
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

/// Worker thread stack size (16 MiB) for Rayon tile workers (`M-DOCUMENTED-MAGIC`).
///
/// Deep CRAM container decompression (`noodles-cram` RANS/bzip2/lzma codecs) combined
/// with recursive structural-variant and indel realignment frames on 50x–150x WGS tiles
/// exceeds Rayon's default 2 MiB worker stack, causing `SIGABRT` stack overflows if set lower.
const WORKER_THREAD_STACK_BYTES: usize = 16 * 1024 * 1024;

/// Executes `worker_fn` in parallel over pre-grouped `batches` (`Vec<Vec<T>>`),
/// streaming flattened results into `consumer_fn` in exact deterministic order (`EBO-3`).
pub fn run_parallel_batched<T, R, F, W>(
    batches: Vec<Vec<T>>,
    num_threads: usize,
    worker_fn: F,
    mut consumer_fn: W,
) -> Result<(), VardictError>
where
    T: Send + Sync + 'static,
    R: Send + 'static,
    F: Fn(T) -> Result<R, VardictError> + Send + Sync + 'static,
    W: FnMut(R) -> Result<(), VardictError>,
{
    if batches.is_empty() {
        return Ok(());
    }
    if num_threads == 1 {
        for batch in batches {
            for item in batch {
                consumer_fn(worker_fn(item)?)?;
            }
        }
        return Ok(());
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .stack_size(WORKER_THREAD_STACK_BYTES)
        .build()
        .map_err(|e| {
            VardictError::Io(std::io::Error::other(format!(
                "Failed to build rayon thread pool: {}",
                e
            )))
        })?;

    let window_size = if num_threads > 0 { num_threads * 4 } else { 4 };
    let (tx, rx) = mpsc::channel::<(usize, Result<Vec<R>, VardictError>)>();

    // Shared state for throttling: (next_expected_batch_idx, is_shutdown) + lock-free AtomicBool
    type SharedState = (Mutex<(usize, bool)>, Condvar, std::sync::atomic::AtomicBool);
    let shared_state: Arc<SharedState> = Arc::new((
        Mutex::new((0_usize, false)),
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
            let mut state = self.state.0.lock().unwrap_or_else(|e| e.into_inner());
            state.1 = true;
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

        // Feeder thread
        s.spawn(move || {
            let (lock, cvar, _) = &*shared_clone;
            for (batch_idx, batch) in batches.into_iter().enumerate() {
                let mut state = match lock.lock() {
                    Ok(g) => g,
                    Err(p) => p.into_inner(),
                };
                while batch_idx >= state.0 + window_size && !state.1 {
                    state = match cvar.wait(state) {
                        Ok(g) => g,
                        Err(p) => p.into_inner(),
                    };
                }
                if state.1 {
                    break;
                }
                drop(state);

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
                        Ok(out)
                    }))
                    .unwrap_or_else(|_| {
                        Err(VardictError::Io(std::io::Error::other(format!(
                            "Worker thread panicked on batch {}",
                            batch_idx
                        ))))
                    });
                    let _ = tx_clone.send((batch_idx, res));
                });
            }
        });

        // Consumer loop on main thread
        let mut next_expected_idx = 0;
        let mut pending = BTreeMap::new();
        let (lock, cvar, _) = &*shared_state;

        for (batch_idx, res) in rx {
            let batch_vals = res?;
            pending.insert(batch_idx, batch_vals);

            let mut advanced = false;
            while let Some(vals) = pending.remove(&next_expected_idx) {
                for val in vals {
                    consumer_fn(val)?;
                }
                next_expected_idx += 1;
                advanced = true;
            }

            if advanced {
                let mut state = match lock.lock() {
                    Ok(g) => g,
                    Err(p) => p.into_inner(),
                };
                state.0 = next_expected_idx;
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
        let result = rx.recv_timeout(std::time::Duration::from_secs(2));
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

        let result = rx.recv_timeout(std::time::Duration::from_secs(2));
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
    }
}
