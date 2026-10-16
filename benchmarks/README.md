# VarDictRust Benchmark Suite & Evidence Ledger

This directory houses the performance evidence, execution logs, benchmark drivers, and environment specifications for **VarDictRust v1.0.1**.

Every benchmark claim is backed by a verifiable, machine-readable evidence artifact in [`benchmarks/evidence/`](evidence/), bound to the certified release engine fingerprint `c57e9661d77ef254`. For complete statistical analyses, biological cohort descriptions, and methodology, see the [Omnibus Empirical Validation Report](../docs/validation/omnibus-empirical-validation.md) and the [Master Test Matrix](../docs/testing/test-matrix-and-ci-audit.md).

---

## 1. Pinned Hardware & Host Environment

All release benchmarks were executed under standardized host constraints (below).
Every performance number in this directory is measured on the stated host
against VarDictJava as the baseline — paper figures are never quoted as ours.
- **CPU/host:** Intel i7-13700H, Linux x86_64, 8 worker threads (Ubuntu 24.04 LTS).
- **Concurrency:** Normalized to **8 CPU threads** (`-th 8`) as a standard multi-tenant cluster baseline (not an engine scaling ceiling).
- **Storage:** Local NVMe SSD.
- **Toolchain:** Rust 1.96.0 (`x86_64-unknown-linux-musl` and `x86_64-unknown-linux-gnu`), OpenJDK 25.0.2 / OpenJDK 8 for oracle (`vardict-oracle`).

---

## 2. Four-Target Performance Budget

Evaluated via `harness/perf_budget.py`, comparing `vardictrust` against the `VarDictJava v1.8.3` pipeline (Java + R + Perl). Requirement: Rust wall time $\le$ Java wall time AND Rust peak RSS $\le$ Java peak RSS.

| Workload / Target | Modality | Java Wall | Rust Wall | Speedup | Java Peak RSS | Rust Peak RSS | RSS Scaling | Cpp Wall / RSS (vardictcpp 3 hc52dbad_0, advisory¹) | Evidence File |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Targeted Amplicon Panel** | AmpliconMode (`-a 10:0.95`) | 0:02.03 | 0:00.07 | **29.0×** | 123.1 MB | 52.7 MB | **0.43× (2.3× lower)** | 0:03.55 / 89.2 MB, advisory [SKIP (CPP-3)](../docs/upstream-bugs.md#cpp-3) | [`amplicon-tier2.json`](evidence/amplicon-tier2.json) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | SimpleMode (`-f 0.01`) | 1:34.18 | 0:45.63 | **2.1×** | 4.1 GB | 1.1 GB | **0.26× (3.8× lower)** | 0:16.87 / 855 MB, advisory [MISMATCH (CPP-9)](../docs/upstream-bugs.md#cpp-9) + [CPP-1](../docs/upstream-bugs.md#cpp-1) | [`hg002-tier2.json`](evidence/hg002-tier2.json) |
| **SEQC2 HCC1395 Somatic Exome** | SomaticMode (`tumor\|normal`, `-th 4`²) | 52:09.24 | 17:13.73 | **3.0×** | 3.5 GB | 1.0 GB | **0.28× (3.6× lower)** | 7:05.12 / 2.1 GB, advisory [MISMATCH (CPP-9)](../docs/upstream-bugs.md#cpp-9) | [`seqc2-exome-tier2.json`](evidence/seqc2-exome-tier2.json) |
| **GIAB HG008 Somatic PDAC** | SomaticMode (`tumor\|normal`, `-th 4`²) | 1:15.71 | 0:48.72 | **1.6×** | 3.3 GB | 0.8 GB | **0.22× (4.5× lower)** | 0:20.54 / 311 MB, advisory [MISMATCH (CPP-9)](../docs/upstream-bugs.md#cpp-9) | [`hg008-tier2.json`](evidence/hg008-tier2.json) |

¹ Cpp (pinned **vardictcpp 3 `hc52dbad_0` bioconda**, same `vardict-oracle` env) runs beside each gate arm via `harness/perf_budget.py`. The budget verdict stays Rust-vs-Java only. Per-cell cpp wall/RSS above are transcribed from the linked evidence `extra` (`cpp_wall_s`, `cpp_maxrss_kb`, ratios). Known cpp divergences are normalized in [`docs/upstream-bugs.md`](../docs/upstream-bugs.md#cpp-divergences). A missing `cpp_time` records a loud `SKIP`, never a silent pass.
² hg008/seqc2-exome at `-th 4`: Java `-th 8` is OOM-killed on 31G hosts (bisected 1✓/4✓/8✗); exact cmds in evidence.

### Reproduction Command:
```bash
python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
```

---

## 3. Master Evidence Ledger (16 Release Artifacts)

| Workload / Target | Evidence Artifact | Verification Tier | Scope / Modality | Output Metric | Result |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **GIAB HG002 Germline** | [`hg002-tier2.json`](evidence/hg002-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | chr20 (10–30 Mbp, 30×) | 227,562 variants | **MATCH** (0 diffs) |
| **GIAB HG002 Truth** | [`hg002-truth.json`](evidence/hg002-truth.json) | [Tier-3](../docs/testing/test-matrix-and-ci-audit.md#tier-3) Concordance | NIST v4.2.1 High-Confidence | **$F_1 = 0.9699$** | **EQUAL** |
| **SEQC2 Somatic Exome** | [`seqc2-exome-tier2.json`](evidence/seqc2-exome-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | HCC1395 / HCC1395BL Somatic | 623,444 variants | **MATCH** (0 diffs) |
| **SEQC2 Somatic Truth** | [`seqc2-truth.json`](evidence/seqc2-truth.json) | [Tier-3](../docs/testing/test-matrix-and-ci-audit.md#tier-3) Concordance | SEQC2 v1.2.1 Somatic Truth | $F_1 = 0.1555$ | **EQUAL** (0 site diffs) |
| **SEQC2 Somatic Slices** | [`seqc2-slices-tier2.json`](evidence/seqc2-slices-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | Multi-Chr Slices (chr1–21) | 4,151 variants | **MATCH** (0 diffs) |
| **GIAB HG008 Somatic PDAC** | [`hg008-tier2.json`](evidence/hg008-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | Tumor/Normal WGS 60×/40× | 18,426 variants | **MATCH** (0 diffs) |
| **GIAB HG008 CRAM Twin** | [`hg008-cram-twin.json`](evidence/hg008-cram-twin.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | CRAM vs BAM Reference Decode | Bitwise identity | **MATCH** (0 diffs) |
| **GIAB HG008 Truth** | [`hg008-truth.json`](evidence/hg008-truth.json) | [Tier-3](../docs/testing/test-matrix-and-ci-audit.md#tier-3) Concordance | NIST Draft v0.3 Somatic Truth | $F_1 = 0.3856$ | **EQUAL** |
| **Targeted Amplicon Panel** | [`amplicon-tier2.json`](evidence/amplicon-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | Primer-trimmed panel (`-a`) | 2 variants | **MATCH** (0 diffs) |
| **GIAB HG005 Truth** | [`hg005-truth.json`](evidence/hg005-truth.json) | [Tier-3](../docs/testing/test-matrix-and-ci-audit.md#tier-3) Concordance | NIST v4.2.1 Han Chinese WGS | **$F_1 = 0.9667$** | **EQUAL** |
| **GIAB NA24385 b37 Truth** | [`na24385-b37-truth.json`](evidence/na24385-b37-truth.json) | [Tier-3](../docs/testing/test-matrix-and-ci-audit.md#tier-3) Concordance | NIST v4.2.1 Ashkenazi Son b37 | **$F_1 = 1.0000$** | **EQUAL** |
| **Extended Modalities** | [`extended-coverage-tier2.json`](evidence/extended-coverage-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | 17 Modes (RNA `-i`, SVs, flags) | 4,993 variants | **MATCH** (0 diffs) |
| **Multi-Modal Slices** | [`multimodal-slices-tier2.json`](evidence/multimodal-slices-tier2.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Differential | 18 Slices across 6 Cohorts | 1,129 variants | **MATCH** (0 diffs) |
| **Chromosome 1 Scaling** | [`scale.json`](evidence/scale.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Scaling | Full chr1 Span (249 Mbp) | $S_8 = 3.63\times$, RSS 1.40× | **PASS** |
| **Performance Budget** | [`perf-budget.json`](evidence/perf-budget.json) | Regression Gate | 4 Target Workload Regression | Wall/RAM gates | **PASS** |
| **Musl CI Packaging** | [`musl-ci.json`](evidence/musl-ci.json) | Release Gate | Static PIE Linux Packaging | Bitwise identity | **PASS** |

---

## 4. Whole Chromosome 1 Scaling & Flat Memory Footprint

Evaluated via `harness/scale_bench.sh` on HCC1395 across all 248,956,422 bp of chromosome 1. Values below are the single-run certified evidence recorded in `scale.json`.

### Multi-Thread Parallel Speedup (chr1 full span):
| Threads (`-th`) | Wall Clock Time | Speedup Ratio |
| :--- | :--- | :--- |
| 1 core | 26m 57s (1616.54 s) | 1.00× (baseline) |
| 2 cores | 13m 49s (828.93 s) | 1.95× |
| 4 cores | 8m 46s (525.72 s) | 3.07× |
| 8 threads | 7m 26s (445.91 s) | **3.63× parallel** |

### Bounded Resident Memory vs Genomic Span (at `-th 8`):
| Genomic Span | Peak RSS |
| :--- | :--- |
| 1 Mbp | 816.1 MB (835,652 KB) |
| 10 Mbp | 1.93 GB (2,023,444 KB) |
| 100 Mbp | 2.14 GB (2,242,148 KB) |
| 249 Mbp (Full chr1) | **2.70 GB (2,829,376 KB)** |

**Observed RSS Ratio (249 Mbp vs 10 Mbp): 1.40×** (Passing pre-registered threshold $\le 1.5\times$). Memory remains bounded across deep genomic spans via tiled execution with boundary halo deduplication.

Pinned **vardictcpp 3 `hc52dbad_0`** runs beside each loop above (`harness/scale_bench.sh`, advisory only): per-thread cpp wall (`cpp_wall_s`, exome BED workload) and per-span cpp RSS (`cpp_maxrss_kb`, `-R` full-span workload) land separately in `scale.json` `extra` — at `-th 8` BED wall 160.9 s, full-span (249 Mbp) RSS 22.5 GB ([CPP-5](../docs/upstream-bugs.md#cpp-5)); `--chunk` is deliberately unpinned (cpp default rewrites the `Seg` column, so pinning a value would fork output shape).

### Reproduction Command:
```bash
DATA_DIR="${DATA_DIR:-scratch/g4}"
bash harness/scale_bench.sh \
  "${DATA_DIR}/ref/d1vd1/GRCh38.d1.vd1.fa" \
  "${DATA_DIR}/seqc2/HCC1395.bam" \
  "${DATA_DIR}/seqc2/exome.engine.bed" \
  HCC1395 \
  "${DATA_DIR}/scale"
```

---

## 5. Multi-Thread Determinism

Guarantees bitwise-identical VCF records across varying worker thread allocations (1, 2, 4, 8 threads).

```bash
cargo test -p vardictrust --test thread_determinism
```

---

## 6. Micro-Benchmarks (Criterion)

In-memory pileup and realignment micro-benchmarks (`crates/vardictrust-core/benches/engine.rs`).

```bash
cargo bench            # Full run, report under target/criterion/
cargo bench --no-run   # Compile-only gate
```

| Benchmark | Target Function | Fixture Profile | Baseline Latency |
| :--- | :--- | :--- | :--- |
| `parse_cigar/100_reads_mixed_cigars` | `vardictrust_core::cigar::parse::pileup` | 100 sorted 150 bp reads over 20 kb synthetic ACGT ref; mixed M/S/D/I | $\approx 21.6\,\mu\text{s}$ / read |
| `find_bp/60bp_consensus_indel40` | `vardictrust_core::realign::find_bp` | 60 bp synthetic consensus (~10% mismatches), indel size 40 | $\approx 5.48\,\mu\text{s}$ |

---

## 7. Dataset Provenance & Asset Acquisition

All external genomic benchmark assets are tracked in `harness/assets.tsv` and recorded upon download in `benchmarks/MANIFEST`:

```bash
# Ingest public references and truth sets
bash harness/fetch_data.sh refs hg002 hg008 seqc2
```
