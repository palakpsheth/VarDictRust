# VarDictRust Profiling Subsystem API & Allocator Reference

Reference companion for [`.agents/skills/vardictrust-profiling/SKILL.md`](../SKILL.md). Documents the exact Rust types, feature gates, allocator configurations, fingerprint inclusion rules, and artifact schemas implemented in `crates/vardictrust/src/main.rs` and `crates/vardictrust/src/profile.rs`.

---

## 1. Global Allocator & `MALLOC_CONF` Configuration (`crates/vardictrust/src/main.rs`)

VarDictRust selects its global allocator at compile time using priority-based feature gates (`crates/vardictrust/src/main.rs:7-30`):

```rust
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
#[cfg(all(feature = "jemalloc", not(feature = "profile-heap")))]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 69] =
    b"background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8\0";
#[cfg(feature = "profile-heap")]
#[export_name = "_rjem_malloc_conf"]
pub static MALLOC_CONF: &[u8; 105] = b"prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0\0";
```

When both `profile-heap` and `profile-dhat` are compiled in simultaneously, `dhat::Alloc` owns `#[global_allocator]`, and `profile::Session::start()` emits a stderr warning and disables heap profiling (`crates/vardictrust/src/profile.rs:125-129`):
```rust
#[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
if heap_dir.is_some() {
    eprintln!("vardictrust: warning: VARDICT_PROFILE_HEAP requested, but global allocator does not match; disabling heap profiling");
    heap_dir = None;
}
```

---

## 2. `profile::Session` Lifecycle API (`crates/vardictrust/src/profile.rs`)

| Item / Method | Signature / Location | Behavior |
| :--- | :--- | :--- |
| `Session::start` | `pub fn start() -> Result<Self>` (`profile.rs:81-218`, called in `main.rs:43` before CLI parsing) | Reads `VARDICT_PROFILE_CPU`, `VARDICT_PROFILE_HEAP`, and `VARDICT_PROFILE_DHAT` from the process environment. Fails fast with `anyhow::bail!("{var_name} must be a non-empty directory path")` if any is set to `""`. If a variable is set when its Cargo feature is disabled, prints `vardictrust: warning: <VAR> is set, but the feature is not compiled in` to stderr and creates no files. When enabled, probes target directories for write access via `.vardict_profile_write_probe` (`probe_dir`), starts `pprof::ProfilerGuardBuilder::default().frequency(499).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build()` with chained `SIGPROF` tick counting, starts `dhat::Profiler::builder().file_name(dir.join("dhat-heap.json")).trim_backtraces(Some(32)).build()`, and activates `jemalloc_pprof::PROF_CTL`. |
| `<Session as Drop>::drop` | `fn drop(&mut self)` (`profile.rs:220-383`) | Restores the original `SIGPROF` handler and writes `cpu-meta.json`, `cpu-profile.pb`, `cpu-collapsed.txt`, and `cpu-flamegraph.svg` (CPU); dumps `heap.pb.gz` (or `heap.pb`) and `heap-flamegraph.svg` via `jemalloc_pprof::PROF_CTL` (Heap); and drops `_dhat: Option<dhat::Profiler>` to flush `dhat-heap.json` (DHAT). |
| `harness/engine_fingerprint.sh` | `GATE_FILES` inclusion (`engine_fingerprint.sh:10,30`) | `harness/engine_fingerprint.sh` hashes the entire `crates/vardictrust/src` tree in both default (`FP`) and `--perf` (`PERF_FP`) modes. Under the Two-Tier testing architecture, routine PRs touching `profile.rs` run `bash harness/pr_check.sh` (`<25s`), while release certification PRs (`release/vX.Y.Z`) run `bash harness/recertify.sh` (including `scale_bench.sh` when `PERF_FP` changes). |

---

## 3. Output Artifact Schemas

### `cpu-meta.json` (`VARDICT_PROFILE_CPU=<DIR>`)
Written on `Session` drop alongside `cpu-profile.pb`, `cpu-collapsed.txt`, and `cpu-flamegraph.svg`:
```json
{"ticks_total":142,"samples_kept":138}
```
When `harness/profile.sh` runs with `PROFILE_LOOPS=N` (`N > 1`), per-iteration collapsed stacks and `cpu-meta.json` counts are merged:
```json
{"ticks_total":426,"samples_kept":414,"loops":3}
```

### `heap.pb.gz` / `heap-flamegraph.svg` (`VARDICT_PROFILE_HEAP=<DIR>`)
Written on `Session` drop via `jemalloc_pprof::PROF_CTL`:
- `heap.pb.gz` (or `heap.pb` if uncompressed): pprof protobuf profile of live in-use heap bytes at exit, decoded by `harness/profile_catalog.py` into `heap-collapsed.txt` and Markdown hotspot tables.
- `heap-flamegraph.svg`: interactive flamegraph SVG rendered directly by `jemalloc_pprof`.

### `dhat-heap.json` (`VARDICT_PROFILE_DHAT=<DIR>`)
Written by `dhat::Profiler` on drop (`dhat-rs` JSON format with `ftbl` frame table and `pps` program points containing `tb` total bytes, `tbk` total blocks, `mb` max live bytes, and `fs` innermost-first frame indices). Consumed by `harness/profile_catalog.py` to emit `dhat-bytes-collapsed.txt`, `dhat-blocks-collapsed.txt`, `dhat-bytes.svg`, and `dhat-blocks.svg`.

---

## 4. Automated Smoke Tests (`crates/vardictrust/tests/profile_smoke.rs`)

- **Always-active tests (3 tests in default build):**
  - `compiled_out_warning_when_features_disabled_and_no_files_created`: verifies `VARDICT_PROFILE_{CPU,DHAT,HEAP}` emit the compiled-out warning to stderr and create zero files when features are disabled.
  - `version_and_help_with_no_env_vars_creates_no_files`: verifies `--version` and `--help` create no files when `VARDICT_PROFILE_*` env vars are unset.
  - `empty_profile_env_var_fails_fast_with_actionable_error`: verifies `VARDICT_PROFILE_{CPU,HEAP,DHAT}=""` exits non-zero with `<VAR> must be a non-empty directory path`.
- **Feature-gated tests:**
  - `cpu_profile_files_exist_and_dir_probe_fails_fast` (`#[cfg(feature = "profile-cpu")]`): runs `fixtures/micro/micro.bam` with `VARDICT_PROFILE_CPU`, asserts `cpu-collapsed.txt`, `cpu-profile.pb`, `cpu-meta.json`, and `cpu-flamegraph.svg` exist, and verifies an unwritable path fails fast before the run.
  - `heap_profile_files_exist_and_dir_probe_fails_fast` (`#[cfg(all(feature = "profile-heap", not(feature = "profile-dhat")))]`): asserts `heap.pb.gz` (or `heap.pb`) and `heap-flamegraph.svg` exist and unwritable paths fail fast.
  - `dhat_profile_files_exist_and_dir_probe_fails_fast` (`#[cfg(feature = "profile-dhat")]`): asserts `dhat-heap.json` exists and unwritable paths fail fast.
  - `heap_dhat_allocator_mismatch_warning` (`#[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]`): verifies the allocator mismatch warning when both features are compiled in.
