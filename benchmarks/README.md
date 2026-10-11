# VarDictRust Benchmark Suite & Evidence Ledger

This directory houses the performance evidence, execution logs, benchmark drivers, and environment specifications for **VarDictRust v1.0.0**.

---

## 1. Pinned Hardware & Host Environment

All release benchmarks were executed under strict host constraints:
- **CPU:** AMD Ryzen (Linux x86_64, Ubuntu 24.04 LTS).
- **Concurrency:** Pinned to **8 CPU threads** (`-th 8`).
- **Storage:** Local NVMe SSD.
- **Toolchain:** Rust 1.96.0 (`x86_64-unknown-linux-musl` and `x86_64-unknown-linux-gnu`), OpenJDK 25.0.2 / OpenJDK 8 for oracle (`vardict-oracle`).

---

## 2. Four-Target Performance Budget

Evaluated via `harness/perf_budget.py`, comparing `vardictrust` against the `VarDictJava v1.8.3` pipeline (Java + R + Perl). Requirement: Rust wall time $\le$ Java wall time AND Rust peak RSS $\le$ Java peak RSS.

| Workload / Target | Modality | Java Wall | Rust Wall | Speedup | Java Peak RSS | Rust Peak RSS | RSS Scaling | Evidence File |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Targeted Amplicon Panel** | AmpliconMode (`-a 10:0.95`) | 0:02.10 | 0:00.07 | **30.0×** | 120.8 MB | 50.8 MB | **0.42× (2.4× lower)** | [`amplicon-tier2.json`](evidence/amplicon-tier2.json) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | SimpleMode (`-f 0.01`) | 4:16.98 | 1:22.47 | **3.1×** | 3.6 GB | 1.1 GB | **0.30× (3.3× lower)** | [`hg002-tier2.json`](evidence/hg002-tier2.json) |
| **SEQC2 HCC1395 Somatic Exome** | SomaticMode (`tumor\|normal`) | 54:23.00 | 19:45.73 | **2.8×** | 4.1 GB | 1.4 GB | **0.36× (2.8× lower)** | [`seqc2-exome-tier2.json`](evidence/seqc2-exome-tier2.json) |
| **GIAB HG008 Somatic PDAC** | SomaticMode (`tumor\|normal`) | 1:10.99 | 0:49.48 | **1.4×** | 5.7 GB | 1.2 GB | **0.21× (4.8× lower)** | [`hg008-tier2.json`](evidence/hg008-tier2.json) |

### Reproduction Command:
```bash
python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2
```

---

## 3. Whole Chromosome 1 Scaling & Flat Memory Footprint

Evaluated via `harness/scale_bench.sh` on HCC1395 across all 248,956,422 bp of chromosome 1.

### Multi-Thread Parallel Speedup (chr1 full span):
| Threads (`-th`) | Wall Clock Time | Speedup Ratio |
| :--- | :--- | :--- |
| 1 core | 57m 16s (3436.81 s) | 1.00× (baseline) |
| 2 cores | 31m 00s (1859.64 s) | 1.85× |
| 4 cores | 14m 56s (895.79 s) | 3.84× |
| 8 cores | 11m 49s (709.22 s) | **4.85× parallel** |

### Bounded Resident Memory vs Genomic Span (at `-th 8`):
| Genomic Span | Peak RSS |
| :--- | :--- |
| 1 Mbp | 813.6 MB (833,156 KB) |
| 10 Mbp | 2.06 GB (2,114,212 KB) |
| 100 Mbp | 2.30 GB (2,355,920 KB) |
| 249 Mbp (Full chr1) | **2.87 GB (2,940,504 KB)** |

**Observed RSS Ratio (249 Mbp vs 10 Mbp): 1.39×** (Passing pre-registered threshold $\le 1.5\times$). Memory remains bounded across deep genomic spans via tiled execution with boundary halo deduplication.

### Reproduction Command:
```bash
bash harness/scale_bench.sh scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa scratch/g4/seqc2/HCC1395.bam scratch/g4/seqc2/exome.engine.bed HCC1395 scratch/g4/scale
```

---

## 4. Multi-Thread Determinism

Guarantees bitwise-identical VCF records across varying worker thread allocations (1, 2, 4, 8 threads).

```bash
cargo test -p vardictrust --test thread_determinism
```

---

## 5. Micro-Benchmarks (Criterion)

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

## 6. Dataset Provenance & Asset Acquisition

All external genomic benchmark assets are tracked in `harness/assets.tsv` and recorded upon download in `benchmarks/MANIFEST`:

```bash
# Ingest public references and truth sets
bash harness/fetch_data.sh refs hg002 hg008 seqc2
```
