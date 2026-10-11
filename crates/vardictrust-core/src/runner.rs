use crate::types::VardictError;
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::sync::{Arc, Condvar, Mutex};

/// Parallel task runner that processes items in parallel using a custom Rayon
/// thread pool, collects them, and streams them out in the original input order.
///
/// Implements bounded backpressure by using a sliding window. Workers are only
/// fed new tasks when the number of in-flight tasks is within the window limit.
pub fn run_parallel<T, R, F, W>(
    items: Vec<T>,
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
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| {
            VardictError::Io(std::io::Error::other(format!(
                "Failed to build thread pool: {}",
                e
            )))
        })?;

    let window_size = if num_threads > 0 { num_threads * 4 } else { 4 };
    let (tx, rx) = mpsc::channel::<(usize, Result<R, VardictError>)>();

    // Mutex state: (next_expected_idx, is_shutdown)
    let shared_state = Arc::new((Mutex::new((0_usize, false)), Condvar::new()));
    let worker_fn = Arc::new(worker_fn);
    let pool = Arc::new(pool);

    // Use a Drop guard to ensure shutdown flag is set and condvar is notified on early exit
    struct ShutdownGuard<'a> {
        state: &'a (Mutex<(usize, bool)>, Condvar),
    }
    impl<'a> Drop for ShutdownGuard<'a> {
        fn drop(&mut self) {
            let (lock, cvar) = self.state;
            let mut state = lock.lock().unwrap_or_else(|e| e.into_inner());
            state.1 = true;
            cvar.notify_all();
        }
    }

    std::thread::scope(|s| {
        let _guard = ShutdownGuard {
            state: &shared_state,
        };

        let feeder_state = shared_state.clone();
        let pool_clone = pool.clone();
        let items_len = items.len();

        s.spawn(move || {
            let (lock, cvar) = &*feeder_state;
            for (idx, item) in items.into_iter().enumerate() {
                let mut state = lock.lock().unwrap_or_else(|e| e.into_inner());
                while idx >= state.0 + window_size && !state.1 {
                    state = cvar.wait(state).unwrap_or_else(|e| e.into_inner());
                }
                if state.1 {
                    break;
                }
                drop(state);

                let tx_clone = tx.clone();
                let w_fn = worker_fn.clone();
                pool_clone.spawn(move || {
                    let res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| w_fn(item)))
                        .unwrap_or_else(|_| {
                            Err(VardictError::Io(std::io::Error::other(format!(
                                "Worker panicked processing tile {}",
                                idx
                            ))))
                        });
                    let _ = tx_clone.send((idx, res));
                });
            }
        });

        let mut next_expected_idx = 0;
        let mut pending = BTreeMap::new();
        let (lock, cvar) = &*shared_state;

        for (idx, res) in rx {
            let val = res?;
            pending.insert(idx, val);

            let mut advanced = false;
            while let Some(val) = pending.remove(&next_expected_idx) {
                consumer_fn(val)?;
                next_expected_idx += 1;
                advanced = true;
            }
            if advanced {
                let mut state = lock.lock().unwrap_or_else(|e| e.into_inner());
                state.0 = next_expected_idx;
                cvar.notify_all();
            }
            if next_expected_idx == items_len {
                break;
            }
        }

        if next_expected_idx < items_len {
            return Err(VardictError::Io(std::io::Error::other(format!(
                "Tile {} failed or panicked; pipeline aborted early",
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
}
