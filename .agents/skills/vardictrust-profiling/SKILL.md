---
name: vardictrust-profiling
description: Use when profiling VarDictRust CPU hotspots, jemalloc heap residency, or DHAT allocation churn, running harness/profile.sh or harness/profile_catalog.py, or validating performance optimizations against VarDictJava.
---

# VarDictRust CPU, Heap & Allocation Profiling (`vardictrust-profiling`)

## Overview

VarDictRust includes an opt-in, zero-default-overhead profiling subsystem in `crates/vardictrust/src/profile.rs` backed by three Cargo features (`profile-cpu`, `profile-heap`, `profile-dhat`) and two harness runners (`harness/profile.sh` and `harness/profile_catalog.py`).

**Core principle:** Scale optimizations must never alter variant output or violate the no-panic rule (`AGENTS.md` §§8–10). Every optimization must be guided by empirical flamegraph/heap/DHAT data and validated against VarDictJava v1.8.3 (`bash harness/rgate.sh` and `bash harness/compare.sh`).

---

## Quick Reference: Cargo Feature & CLI Flag Matrix

| Cargo Feature | Backend Crate | CLI Flag | Output Artifacts in `<DIR>` | Allocator / Mutual Exclusion Rules |
| :--- | :--- | :--- | :--- | :--- |
| *(default release)* | `tikv-jemallocator` | *(none)* | None (zero profiling overhead) | Static `MALLOC_CONF`: `background_thread:true,dirty_decay_ms:250,muzzy_decay_ms:0,narenas:8` |
| `profile-cpu` | `pprof` (`0.15`, `flamegraph` + `protobuf-codec`) | `--profile-cpu <DIR>` | `cpu-flamegraph.svg`, `cpu-profile.pb`, `profile-meta.json` | Compatible with `profile-heap` (`--mode both`); default `PROFILE_HZ=999` |
| `profile-heap` | `tikv_jemalloc_ctl` (`stats` + `prof_ctl`) | `--profile-heap <DIR>` | `heap-summary.json`, `heap-final.prof`, `profile-meta.json` | Requires `unprefixed_malloc_on_supported_platforms`; static `MALLOC_CONF` sets `prof:true,prof_active:true,lg_prof_sample:19` |
| `profile-dhat` | `dhat::Alloc` | `--profile-dhat <DIR>` | `dhat-heap.json`, `profile-meta.json` | **Mutually exclusive with `profile-heap`** (`compile_error!` at `profile.rs:10`); run with `-th 1` to avoid global allocator mutex contention |

See [`references/rust-profiling-api.md`](references/rust-profiling-api.md) for the complete Rust API (`ProfilingSession`), static `MALLOC_CONF` byte slices, and JSON schemas.

---

## Profiling Workflows & Commands

### 1. Single-Case Profiling (`harness/profile.sh`)
Use `harness/profile.sh` to build the required feature binary (`target/release-profile/vardictrust` or `target/release-dhat/vardictrust`) and profile a preset (`s1a`, `s1b`, `s2`, `s3`, `s4`, `s6`, `s7`) or custom BAM/BED slice:
```bash
# Profile both CPU flamegraph (999 Hz) and jemalloc heap on preset s1a (CEPH1463 chr20)
bash harness/profile.sh --mode both --preset s1a --threads 4 --out /tmp/prof-s1a

# Profile DHAT allocation churn (total bytes/blocks allocated, short-lived Vec/HashMap churn)
bash harness/profile.sh --mode dhat --preset s1a --threads 1 --out /tmp/prof-s1a-dhat

# Run all three modes sequentially (cpu+heap binary first, then separate dhat binary)
bash harness/profile.sh --mode all --preset s6 --threads 8 --out /tmp/prof-s6
```
- On shared or memory-constrained hosts, wrap long runs in `harness/courtesy.sh`:
  `bash harness/courtesy.sh bash harness/profile.sh --mode all --preset s1b --threads 8`

### 2. Full 15-Cell Multi-Modal Catalog (`harness/profile_catalog.py`)
To benchmark and profile across the 8 representative workloads (`s1a`, `s1b`, `s3`, `rna-i`, `hg002-sub`, `hg008-sub`, `hg008`, `wgs-scale`) at `th1` and `th8`:
```bash
# List all 15 cells (8 workloads x th1/th8; s3 amplicon is single-thread only)
python3 harness/profile_catalog.py --list

# Run a targeted subset of cells with CPU + jemalloc + DHAT profiling
python3 harness/profile_catalog.py --cells th1-s1a,th8-hg008-sub --out-dir /tmp/profile-run
```

### 3. Inspecting Hotspots & Allocation Churn
- **CPU (`cpu-flamegraph.svg` / `cpu-profile.pb`):** Look for wide stacks in `cigar::parse`, `realign`, `tovars`, `dup_var_Consensus`, or `noodles` BAM/CRAM record decoding.
- **Heap Residency (`heap-summary.json`):** Inspect `allocated_bytes`, `active_bytes`, `resident_bytes`, and `retained_bytes` across `start` vs `finish` snapshots to confirm tile state drops cleanly (RAII, `AGENTS.md` §10).
- **DHAT Churn (`dhat-heap.json`):** Sort call sites by `total_bytes` and `total_blocks` to spot per-read or per-position heap allocations that can use stack buffers, `SmallVec`, or cleared buffer reuse.

### 4. Mandatory Optimization Verification Loop
```bash
# 1. Verify smoke tests for all profiling feature combinations
cargo test -p vardictrust --test profile_smoke
cargo test -p vardictrust --test profile_smoke --features "profile-cpu,profile-heap"
cargo test -p vardictrust --test profile_smoke --features "profile-dhat"

# 2. Verify zero output divergence vs VarDictJava v1.8.3
bash harness/rgate.sh

# 3. If engine code changed, re-certify live evidence (including 3-run scale.json)
# bash harness/courtesy.sh bash harness/recertify.sh
```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "I'll enable `profile-cpu`, `profile-heap`, and `profile-dhat` in one `cargo build --all-features` call." | `profile-heap` (`tikv-jemallocator`) and `profile-dhat` (`dhat::Alloc`) both replace `#[global_allocator]` and trigger `compile_error!`. Use `--features "profile-cpu,profile-heap"` and `--features "profile-dhat"` separately (or `harness/profile.sh --mode all`). |
| "Replacing `IndexMap` or `JavaHashMap` with `std::collections::HashMap` or `FxHashMap` is 10% faster." | **Forbidden (`AGENTS.md` §8).** Map iteration order reaches realignment and `.var`/VCF output; changing it silently breaks VarDictJava v1.8.3 parity. |
| "Using unchecked indexing `seq[i]` in `vardictrust-core` will speed up the inner loop." | **Forbidden (`AGENTS.md` §9).** `vardictrust-core` enforces `#![forbid(unsafe_code)]` and zero panics/unwrap/bare indexing on untrusted input. Use `.get(i)` or slice iterators. |

## Red Flags — STOP Immediately If You Catch Yourself

- Leaving `profile-cpu`, `profile-heap`, or `profile-dhat` in `default = [...]` in `crates/vardictrust/Cargo.toml`.
- Caching whole-chromosome pileup or read state across tiles (`AGENTS.md` §10 requires bounded per-tile + halo memory).
- Committing `/tmp/prof-*`, `*.svg`, `*.pb`, `*.prof`, or `dhat-heap.json` artifacts into the Git repository.
