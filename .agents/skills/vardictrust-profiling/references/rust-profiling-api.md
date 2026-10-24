# VarDictRust Profiling Subsystem API & Allocator Reference

Reference companion for [`.agents/skills/vardictrust-profiling/SKILL.md`](../SKILL.md). Documents the exact Rust types, feature gates, allocator configurations, and JSON artifact schemas implemented in `crates/vardictrust/src/main.rs` and `crates/vardictrust/src/profile.rs`.

---

## 1. Global Allocator & `MALLOC_CONF` Configuration (`crates/vardictrust/src/main.rs`)

VarDictRust selects its global allocator at compile time based on Cargo features:

```rust
// crates/vardictrust/src/main.rs:11-26
#[cfg(all(
    not(target_env = "msvc"),
    not(target_os = "macos"),
    not(feature = "profile-dhat")
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

// Production default (profile-heap OFF, profile-dhat OFF):
// Rapid dirty-page decay (250 ms), muzzy decay disabled, 8 arenas to bound multi-thread RSS.
#[cfg(all(
    not(target_env = "msvc"),
    not(target_os = "macos"),
    not(feature = "profile-heap"),
    not(feature = "profile-dhat")
))]
#[allow(non_upper_case_globals)]
#[export_name = "_rjem_malloc_conf"]
pub static malloc_conf: &[u8; 69] =
    b"background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0";

// Heap profiling build (profile-heap ON):
// Unprefixed malloc symbol (_rjem_ -> malloc_conf) with sampling every 2^19 bytes (512 KiB).
#[cfg(all(
    not(target_env = "msvc"),
    not(target_os = "macos"),
    feature = "profile-heap"
))]
#[allow(non_upper_case_globals)]
#[export_name = "malloc_conf"]
pub static malloc_conf: &[u8; 105] =
    b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0";
```

When `feature = "profile-dhat"` is enabled (`crates/vardictrust/src/profile.rs:14-16`):
```rust
#[cfg(feature = "profile-dhat")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;
```

Enabling both `profile-heap` and `profile-dhat` is rejected at compile time (`crates/vardictrust/src/profile.rs:9-12`):
```rust
#[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
compile_error!(
    "features `profile-heap` (jemalloc) and `profile-dhat` (dhat::Alloc) are mutually exclusive global allocators"
);
```

---

## 2. `ProfilingSession` Lifecycle API (`crates/vardictrust/src/profile.rs`)

| Function / Method | Signature | Behavior |
| :--- | :--- | :--- |
| `ProfilingSession::start` | `pub fn start(cli: &Cli) -> Result<Self>` | Validates that requested `--profile-*` flags have their Cargo features compiled in (returns actionable `anyhow::bail!` otherwise). Starts `dhat::Profiler`, `pprof::ProfilerGuardBuilder::default().frequency(hz).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build()`, and snapshots initial jemalloc `stats::allocated`/`active`/`resident`/`mapped`/`retained`. |
| `ProfilingSession::finish` | `pub fn finish(mut self, cli: &Cli) -> Result<()>` | Stops active guards, writes `cpu-flamegraph.svg` + `cpu-profile.pb` (CPU), `heap-summary.json` + `heap-final.prof` via `prof_ctl::dump` (Heap), `dhat-heap.json` (DHAT), and `profile-meta.json` in each output directory. |
| `harness/engine_fingerprint.sh` | `GATE_FILES` exclusion | `crates/vardictrust/src/profile.rs` is explicitly excluded from `GATE_FILES` (`engine_fingerprint.sh:10-12`) so profiling instrumentation tweaks do not invalidate variant-calling parity evidence. |

---

## 3. Output Artifact Schemas

### `profile-meta.json`
Written to every `--profile-cpu`, `--profile-heap`, or `--profile-dhat` directory:
```json
{
  "schema_version": 1,
  "version": "1.0.3",
  "mode": "simple",
  "threads": 8,
  "wall_seconds": 1.4231,
  "peak_rss_kib": 184320,
  "cpu_profile": true,
  "cpu_hz": 999,
  "heap_profile": true,
  "dhat_profile": false
}
```

### `heap-summary.json`
Written when `--profile-heap <DIR>` is used:
```json
{
  "schema_version": 1,
  "version": "1.0.3",
  "mode": "simple",
  "threads": 8,
  "wall_seconds": 1.4231,
  "peak_rss_kib": 184320,
  "jemalloc_start": {
    "allocated_bytes": 40960,
    "active_bytes": 65536,
    "metadata_bytes": 4505600,
    "resident_bytes": 4653056,
    "mapped_bytes": 8388608,
    "retained_bytes": 0
  },
  "jemalloc_finish": {
    "allocated_bytes": 524288,
    "active_bytes": 786432,
    "metadata_bytes": 5120000,
    "resident_bytes": 12582912,
    "mapped_bytes": 25165824,
    "retained_bytes": 4194304
  },
  "heap_dump_file": "heap-final.prof",
  "heap_dump_status": "ok"
}
```

---

## 4. Automated Smoke Tests (`crates/vardictrust/tests/profile_smoke.rs`)

- **Default build (features disabled):** Verifies `--profile-cpu`, `--profile-heap`, and `--profile-dhat` fail fast with an explicit error message instructing the user which `--features` flag to compile with.
- **`--features "profile-cpu,profile-heap"`:** Runs `vardictrust` on `golden/ce11_chrM.fa` + `golden/test.bam` and asserts `cpu-flamegraph.svg`, `cpu-profile.pb`, `heap-summary.json`, and `profile-meta.json` are non-empty.
- **`--features "profile-dhat"`:** Runs `vardictrust` with `--profile-dhat` and asserts `dhat-heap.json` and `profile-meta.json` are generated and valid JSON.
