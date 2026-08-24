use crate::types::VardictError;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::sync::mpsc;

/// Parallel task runner that processes items in parallel using a custom Rayon
/// thread pool, collects them, and streams them out in the original input order.
///
/// Implements bounded backpressure by using a bounded sync channel (capacity 10).
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
    // Use Rayon global pool or build a custom one if threads > 0
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| {
            VardictError::Io(std::io::Error::other(format!(
                "Failed to build thread pool: {}",
                e
            )))
        })?;

    let (tx, rx) = mpsc::sync_channel::<(usize, Result<R, VardictError>)>(10);

    let tx_clone = tx.clone();
    pool.spawn(move || {
        let enumerated: Vec<(usize, T)> = items.into_iter().enumerate().collect();
        enumerated
            .into_par_iter()
            .for_each_with(tx_clone, |sender, (idx, item)| {
                let res = worker_fn(item);
                let _ = sender.send((idx, res));
            });
    });
    drop(tx);

    let mut next_expected_idx = 0;
    let mut pending = BTreeMap::new();

    for (idx, res) in rx {
        let val = res?;
        pending.insert(idx, val);

        while let Some(val) = pending.remove(&next_expected_idx) {
            consumer_fn(val)?;
            next_expected_idx += 1;
        }
    }

    Ok(())
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
                // Induce potential out-of-order execution
                if item % 2 == 0 {
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
}
