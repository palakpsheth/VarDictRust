---
name: vardictrust-profiling
description: Use when profiling VarDictRust CPU hotspots, jemalloc heap residency, or DHAT allocation churn, running harness/profile.sh or harness/profile_catalog.py, or validating performance optimizations against VarDictJava.
---

# VarDictRust CPU, Heap & Allocation Profiling (`vardictrust-profiling`)

## Overview

VarDictRust includes an opt-in, zero-default-overhead profiling subsystem in `crates/vardictrust/src/profile.rs` (`profile::Session::start()` + `Drop`) backed by three Cargo features (`profile-cpu`, `profile-heap`, `profile-dhat`), activated via `VARDICT_PROFILE_*` environment variables, and orchestrated by `harness/profile.sh` and `harness/profile_catalog.py`.

**Core principle:** Scale optimizations must never alter variant output or violate the no-panic rule (`AGENTS.md` §§8–10). Every optimization must be guided by empirical flamegraph/heap/DHAT data and validated against VarDictJava v1.8.3 (`bash harness/rgate.sh` and `bash harness/compare.sh`).

---

## Quick Reference: Cargo Feature & Environment Variable Matrix

| Cargo Feature | Backend Crate | Environment Variable | Output Artifacts in `<DIR>` | Allocator / Mutual Exclusion Rules |
| :--- | :--- | :--- | :--- | :--- |
| *(default release)* | `tikv-jemallocator` (`jemalloc` feature) | *(none)* | None (zero profiling overhead; compiled-out warning if `VARDICT_PROFILE_*` set) | Static `_rjem_malloc_conf`: `background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8` |
| `profile-cpu` | `pprof` (`0.15`, `flamegraph` + `protobuf-codec`) + `inferno` | `VARDICT_PROFILE_CPU=<DIR>` | `cpu-flamegraph.svg`, `cpu-profile.pb`, `cpu-collapsed.txt`, `cpu-meta.json` | 499 Hz `SIGPROF` sampling with chained tick counter (`ticks_total` vs `samples_kept`) |
| `profile-heap` | `tikv-jemallocator` (`profiling`) + `jemalloc_pprof` | `VARDICT_PROFILE_HEAP=<DIR>` | `heap.pb.gz` (or `heap.pb`), `heap-flamegraph.svg` | Static `_rjem_malloc_conf`: `prof:true,prof_active:true,lg_prof_sample:19,background_thread:true,dirty_decay_ms:5000,muzzy_decay_ms:0` |
| `profile-dhat` | `dhat::Alloc` + `dhat::Profiler` (`trim_backtraces(Some(32))`) | `VARDICT_PROFILE_DHAT=<DIR>` | `dhat-heap.json` | Replaces `#[global_allocator]` with `dhat::Alloc` (`main.rs:9-14`); if `profile-heap` is also compiled in, heap profiling is disabled with a stderr warning |

See [`references/rust-profiling-api.md`](references/rust-profiling-api.md) for the complete `profile::Session` lifecycle (`start()` + `Drop`), static `MALLOC_CONF` byte arrays, `engine_fingerprint.sh` inclusion rules, and artifact schemas.

---

## Profiling Workflows & Commands

### 1. Multi-Modal Profiling Runner (`harness/profile.sh`)
Use `harness/profile.sh` to build the three `--profile profiling` binaries (`target/profile-{cpu,heap,dhat}/profiling/vardictrust`), run named benchmark cells under `harness/courtesy.sh` (and `systemd-run --user --scope -p MemoryMax=14G` where available), and generate `scratch/profiles/<tag>/catalog.md`:
```bash
# Dry-run the profiling plan without building or executing
bash harness/profile.sh baseline --dry-run

# Profile a single slice cell (cpu + heap + dhat) into scratch/profiles/baseline/s1a/
bash harness/profile.sh baseline --only s1a

# Profile an opt-in bounded DHAT sub-slice or full SEQC2 exome
bash harness/profile.sh optimized --only hg008-sub
bash harness/profile.sh optimized --only hg008-cram-sub
bash harness/profile.sh optimized --only seqc2-slices-sub
PROFILE_THREADS=4 bash harness/profile.sh optimized --only seqc2-exome

# Run all 32 default README modality cells (germline, somatic, --fisher, amplicon, pileup -p,
# splicing -i, SV, CRISPR -J/-j, CRAM twin, extended CLI flags, hg002, hg008, wgs-scale) with 3 slice CPU loops
PROFILE_LOOPS=3 PROFILE_THREADS=8 bash harness/profile.sh optimized
```

### 2. Ranking Hotspots into a Markdown Catalog (`harness/profile_catalog.py`)
`harness/profile_catalog.py` consumes one or more existing profile cell directories, attributes samples/allocations to the innermost `vardictrust_core::` (fallback `vardictrust::`) frame, and writes `dhat-bytes-collapsed.txt`, `dhat-blocks-collapsed.txt`, `heap-collapsed.txt` (plus `dhat-bytes.svg` / `dhat-blocks.svg` if `inferno-flamegraph` is on `PATH`):
```bash
# Rank hotspots across existing cell directories into catalog.md (folding rows < 2% into 'other')
python3 harness/profile_catalog.py scratch/profiles/optimized/s1a scratch/profiles/optimized/hg008 -o scratch/profiles/optimized/catalog.md --min-pct 2 --strict
```

### 3. Mandatory Optimization Verification Loop
```bash
# 1. Verify smoke tests for all profiling feature combinations
cargo test -p vardictrust --test profile_smoke
cargo test -p vardictrust --test profile_smoke --features "profile-cpu,profile-heap"
cargo test -p vardictrust --test profile_smoke --features "profile-dhat"

# 2. Verify zero output divergence vs VarDictJava v1.8.3
bash harness/rgate.sh

# 3. Any change under crates/vardictrust-core/src or crates/vardictrust/src (including profile.rs)
#    changes FP and PERF_FP — re-certify live evidence before opening a PR:
# bash harness/courtesy.sh bash harness/recertify.sh
```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "I'll enable `profile-heap` and `profile-dhat` together to collect jemalloc and DHAT profiles in one run." | `dhat::Alloc` overrides `#[global_allocator]` when `profile-dhat` is enabled (`main.rs:9-14`), causing `Session::start()` to disable `VARDICT_PROFILE_HEAP` (`profile.rs:125-129`). Build and run `profile-heap` and `profile-dhat` separately (as `harness/profile.sh` does). |
| "Replacing `IndexMap` or `JavaHashMap` with `std::collections::HashMap` or `FxHashMap` is 10% faster." | **Forbidden (`AGENTS.md` §8).** Map iteration order reaches realignment and `.var`/VCF output; changing it silently breaks VarDictJava v1.8.3 parity. |
| "Using unchecked indexing `seq[i]` in `vardictrust-core` will speed up the inner loop." | **Forbidden (`AGENTS.md` §9).** `vardictrust-core` enforces `#![forbid(unsafe_code)]` and zero panics/unwrap/bare indexing on untrusted input. Use `.get(i)` or slice iterators. |

## Red Flags — STOP Immediately If You Catch Yourself

- Leaving `profile-cpu`, `profile-heap`, or `profile-dhat` in `default = [...]` in `crates/vardictrust/Cargo.toml`.
- Editing `crates/vardictrust/src/profile.rs` without running `bash harness/recertify.sh` (`crates/vardictrust/src` is hashed by `harness/engine_fingerprint.sh`).
- Committing `scratch/profiles/`, `/tmp/prof-*`, `*.svg`, `*.pb`, `*.pb.gz`, or `dhat-heap.json` artifacts into the Git repository.
