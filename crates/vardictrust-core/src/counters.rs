//! Opt-in per-thread work counters for profiling data-dependent complexity drivers.
//!
//! When compiled with `--features profile-counters`, each thread accumulates
//! counts into a thread-local `[Cell<u64>; NUM_COUNTERS]` array without atomic
//! cache-line contention on hot paths. When `profile-counters` is disabled (the
//! default for all release/dist builds), [`add`], [`inc`], and [`count!`]
//! compile down to zero-cost no-ops.

#[cfg(feature = "profile-counters")]
use std::cell::Cell;
#[cfg(feature = "profile-counters")]
use std::sync::Mutex;

/// Total number of registered work counters.
pub const NUM_COUNTERS: usize = 16;

/// Work counter identifiers keyed by engine mechanism.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum Ctr {
    /// Calls to `sv::partial_pipeline`.
    SvPartialPipeline = 0,
    /// Position keys snapshotted into `BTreeSet`s inside `sv::partial_pipeline`.
    SvSnapshotKeys = 1,
    /// Calls to `realign::partial_pipeline`.
    RealignPartialPipeline = 2,
    /// `loaded_regions` entries scanned in `Reference::get_base`.
    RefRegionsScanned = 3,
    /// `.crai` index entries examined during CRAM region queries.
    CraiRowsExamined = 4,
    /// `.crai` index entries matched during CRAM region queries.
    CraiRowsMatched = 5,
    /// CRAM container records visited during slice decode.
    CramRecordsVisited = 6,
    /// CRAM container records yielded after alignment-span filtering.
    CramRecordsYielded = 7,
    /// CRAM data containers decoded.
    CramContainerDecodes = 8,
    /// Calls to `sv::is_best_sv_used`.
    SvBestUsedCalls = 9,
    /// Soft-clip clusters visited inside `sv::is_best_sv_used`.
    SvSclipsVisited = 10,
    /// Candidate pairs visited in `realignlgins30`.
    Lgins30PairsVisited = 11,
    /// In-window candidate pairs evaluated in `realignlgins30`.
    Lgins30PairsInWindow = 12,
    /// Calls to `bam::query_region_iter`.
    QueryRegionCalls = 13,
    /// BAM/CRAM reader opens.
    ReaderOpens = 14,
    /// BAI/CRAI index bytes read.
    IndexBytesRead = 15,
}

impl Ctr {
    /// All counter variants in index order.
    pub const ALL: [Self; NUM_COUNTERS] = [
        Self::SvPartialPipeline,
        Self::SvSnapshotKeys,
        Self::RealignPartialPipeline,
        Self::RefRegionsScanned,
        Self::CraiRowsExamined,
        Self::CraiRowsMatched,
        Self::CramRecordsVisited,
        Self::CramRecordsYielded,
        Self::CramContainerDecodes,
        Self::SvBestUsedCalls,
        Self::SvSclipsVisited,
        Self::Lgins30PairsVisited,
        Self::Lgins30PairsInWindow,
        Self::QueryRegionCalls,
        Self::ReaderOpens,
        Self::IndexBytesRead,
    ];

    /// Stable snake_case key used in JSON reports.
    #[inline]
    pub const fn name(self) -> &'static str {
        match self {
            Self::SvPartialPipeline => "sv_partial_pipeline",
            Self::SvSnapshotKeys => "sv_snapshot_keys",
            Self::RealignPartialPipeline => "realign_partial_pipeline",
            Self::RefRegionsScanned => "ref_regions_scanned",
            Self::CraiRowsExamined => "crai_rows_examined",
            Self::CraiRowsMatched => "crai_rows_matched",
            Self::CramRecordsVisited => "cram_records_visited",
            Self::CramRecordsYielded => "cram_records_yielded",
            Self::CramContainerDecodes => "cram_container_decodes",
            Self::SvBestUsedCalls => "sv_best_used_calls",
            Self::SvSclipsVisited => "sv_sclips_visited",
            Self::Lgins30PairsVisited => "lgins30_pairs_visited",
            Self::Lgins30PairsInWindow => "lgins30_pairs_in_window",
            Self::QueryRegionCalls => "query_region_calls",
            Self::ReaderOpens => "reader_opens",
            Self::IndexBytesRead => "index_bytes_read",
        }
    }

    /// Numeric index into the counter array.
    #[inline(always)]
    pub const fn index(self) -> usize {
        self as usize
    }
}

/// Immutable snapshot of counter values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CounterSnapshot {
    values: [u64; NUM_COUNTERS],
}

impl Default for CounterSnapshot {
    fn default() -> Self {
        Self::new()
    }
}

impl CounterSnapshot {
    /// Creates a zero-initialized counter snapshot.
    #[inline]
    pub const fn new() -> Self {
        Self {
            values: [0u64; NUM_COUNTERS],
        }
    }

    /// Returns the recorded count for `counter`.
    #[inline]
    pub fn get(&self, counter: Ctr) -> u64 {
        self.values.get(counter.index()).copied().unwrap_or(0)
    }

    /// Returns the raw counter slice in [`Ctr::ALL`] order.
    #[inline]
    pub fn as_slice(&self) -> &[u64; NUM_COUNTERS] {
        &self.values
    }

    /// Serializes the snapshot as a deterministic JSON object string.
    pub fn to_json(&self) -> String {
        use std::fmt::Write;
        let mut out = String::with_capacity(512);
        out.push_str("{\n");
        for (i, ctr) in Ctr::ALL.iter().enumerate() {
            let val = self.get(*ctr);
            let _ = write!(out, "  \"{}\": {}", ctr.name(), val);
            if i + 1 < NUM_COUNTERS {
                out.push(',');
            }
            out.push('\n');
        }
        out.push('}');
        out
    }
}

#[cfg(feature = "profile-counters")]
thread_local! {
    static LOCAL: [Cell<u64>; NUM_COUNTERS] = const { [const { Cell::new(0) }; NUM_COUNTERS] };
}

#[cfg(feature = "profile-counters")]
static SHARDS: Mutex<[u64; NUM_COUNTERS]> = Mutex::new([0u64; NUM_COUNTERS]);

/// Returns `true` when `profile-counters` instrumentation is compiled in.
#[inline(always)]
pub const fn is_enabled() -> bool {
    cfg!(feature = "profile-counters")
}

/// Adds `delta` to `counter` on the current thread when `profile-counters` is enabled.
/// Compiles to a no-op when `profile-counters` is disabled.
#[cfg(feature = "profile-counters")]
#[inline(always)]
pub fn add(counter: Ctr, delta: u64) {
    LOCAL.with(|cells| {
        if let Some(cell) = cells.get(counter.index()) {
            cell.set(cell.get().saturating_add(delta));
        }
    });
}

/// Adds `delta` to `counter` on the current thread when `profile-counters` is enabled.
/// Compiles to a no-op when `profile-counters` is disabled.
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn add(_counter: Ctr, _delta: u64) {}

/// Increments `counter` by 1 on the current thread when `profile-counters` is enabled.
/// Compiles to a no-op when `profile-counters` is disabled.
#[inline(always)]
pub fn inc(counter: Ctr) {
    add(counter, 1);
}

/// Reads the current thread's local counter value (0 when `profile-counters` is disabled).
#[cfg(feature = "profile-counters")]
#[inline]
pub fn get_local(counter: Ctr) -> u64 {
    LOCAL.with(|cells| cells.get(counter.index()).map(Cell::get).unwrap_or(0))
}

/// Reads the current thread's local counter value (0 when `profile-counters` is disabled).
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn get_local(_counter: Ctr) -> u64 {
    0
}

/// Flushes the current thread's local counters into the process-wide shard accumulator
/// and resets the thread-local cells to zero.
#[cfg(feature = "profile-counters")]
pub fn flush_thread_local() {
    let mut local_vals = [0u64; NUM_COUNTERS];
    let mut any_nonzero = false;
    LOCAL.with(|cells| {
        for (dst, cell) in local_vals.iter_mut().zip(cells.iter()) {
            let v = cell.replace(0);
            *dst = v;
            if v != 0 {
                any_nonzero = true;
            }
        }
    });
    if any_nonzero {
        if let Ok(mut totals) = SHARDS.lock() {
            for (tot, v) in totals.iter_mut().zip(local_vals.iter()) {
                *tot = tot.saturating_add(*v);
            }
        }
    }
}

/// Flushes the current thread's local counters into the process-wide shard accumulator.
/// Compiles to a no-op when `profile-counters` is disabled.
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn flush_thread_local() {}

/// Returns a snapshot of the current thread's counters without touching process-wide totals.
#[cfg(feature = "profile-counters")]
pub fn snapshot_local() -> CounterSnapshot {
    let mut values = [0u64; NUM_COUNTERS];
    LOCAL.with(|cells| {
        for (dst, cell) in values.iter_mut().zip(cells.iter()) {
            *dst = cell.get();
        }
    });
    CounterSnapshot { values }
}

/// Returns a snapshot of the current thread's counters (all zeros when disabled).
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn snapshot_local() -> CounterSnapshot {
    CounterSnapshot::new()
}

/// Returns a combined snapshot of flushed worker totals plus the current thread's
/// unflushed counters.
#[cfg(feature = "profile-counters")]
pub fn snapshot() -> CounterSnapshot {
    let mut values = [0u64; NUM_COUNTERS];
    if let Ok(totals) = SHARDS.lock() {
        values = *totals;
    }
    LOCAL.with(|cells| {
        for (dst, cell) in values.iter_mut().zip(cells.iter()) {
            *dst = dst.saturating_add(cell.get());
        }
    });
    CounterSnapshot { values }
}

/// Returns a combined snapshot of counters (all zeros when `profile-counters` is disabled).
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn snapshot() -> CounterSnapshot {
    CounterSnapshot::new()
}

/// Resets the current thread's local counters to zero.
#[cfg(feature = "profile-counters")]
pub fn reset_local() {
    LOCAL.with(|cells| {
        for cell in cells.iter() {
            cell.set(0);
        }
    });
}

/// Resets the current thread's local counters to zero (no-op when disabled).
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn reset_local() {}

/// Resets both the current thread's local counters and the process-wide accumulator.
#[cfg(feature = "profile-counters")]
pub fn reset_all() {
    reset_local();
    if let Ok(mut totals) = SHARDS.lock() {
        *totals = [0u64; NUM_COUNTERS];
    }
}

/// Resets all counters (no-op when `profile-counters` is disabled).
#[cfg(not(feature = "profile-counters"))]
#[inline(always)]
pub fn reset_all() {}

/// Zero-default-cost work counter macro.
///
/// Expands to nothing when the `profile-counters` feature is disabled, so
/// arguments are not evaluated and no code is emitted in default builds.
#[macro_export]
macro_rules! count {
    ($c:ident, $n:expr) => {{
        #[cfg(feature = "profile-counters")]
        {
            $crate::counters::add($crate::counters::Ctr::$c, $n as u64);
        }
    }};
    ($c:ident) => {{
        #[cfg(feature = "profile-counters")]
        {
            $crate::counters::inc($crate::counters::Ctr::$c);
        }
    }};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counter_names_and_indices_are_bijective() {
        for (idx, ctr) in Ctr::ALL.iter().copied().enumerate() {
            assert_eq!(ctr.index(), idx);
            assert!(!ctr.name().is_empty());
        }
    }

    #[cfg(feature = "profile-counters")]
    #[test]
    fn counters_increment_and_snapshot_when_enabled() {
        assert!(is_enabled());
        reset_local();
        let before = snapshot_local();
        assert_eq!(before.get(Ctr::SvPartialPipeline), 0);

        inc(Ctr::SvPartialPipeline);
        add(Ctr::SvSnapshotKeys, 42);
        crate::count!(RealignPartialPipeline);
        crate::count!(RefRegionsScanned, 7usize);

        let after_local = snapshot_local();
        assert_eq!(after_local.get(Ctr::SvPartialPipeline), 1);
        assert_eq!(after_local.get(Ctr::SvSnapshotKeys), 42);
        assert_eq!(after_local.get(Ctr::RealignPartialPipeline), 1);
        assert_eq!(after_local.get(Ctr::RefRegionsScanned), 7);
        assert_eq!(get_local(Ctr::SvSnapshotKeys), 42);

        let before_global = snapshot();
        flush_thread_local();
        assert_eq!(get_local(Ctr::SvSnapshotKeys), 0);
        let after_flush = snapshot();
        assert!(
            after_flush.get(Ctr::SvSnapshotKeys) >= before_global.get(Ctr::SvSnapshotKeys),
            "flushed snapshot must preserve thread-local increments"
        );
        assert!(after_flush.to_json().contains("\"sv_snapshot_keys\": "));
    }

    #[cfg(not(feature = "profile-counters"))]
    #[test]
    fn counters_are_zero_cost_noop_when_disabled() {
        assert!(!is_enabled());
        reset_local();
        inc(Ctr::SvPartialPipeline);
        add(Ctr::SvSnapshotKeys, 999);
        let side_effect_evaluated = std::cell::Cell::new(false);
        crate::count!(RefRegionsScanned, {
            side_effect_evaluated.set(true);
            123usize
        });
        assert!(
            !side_effect_evaluated.get(),
            "count! macro must not evaluate expression when profile-counters is disabled"
        );
        assert_eq!(get_local(Ctr::SvPartialPipeline), 0);
        assert_eq!(get_local(Ctr::SvSnapshotKeys), 0);
        let snap = snapshot();
        for ctr in Ctr::ALL {
            assert_eq!(snap.get(ctr), 0);
        }
    }
}
