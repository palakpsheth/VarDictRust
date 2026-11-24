# VarDictRust Benchmark Suite & Evidence Ledger

This directory houses the performance evidence, execution logs, benchmark drivers, and environment specifications for **VarDictRust v1.0.10**.

Every benchmark claim is backed by a verifiable, machine-readable evidence artifact in [`benchmarks/evidence/`](evidence/) (evaluating default `--compat=vardictjava` parity execution), bound to the certified release engine fingerprint `2d82664ff81000f7`. For complete statistical analyses, biological cohort descriptions, and methodology, see the [Omnibus Empirical Validation Report](../docs/validation/omnibus-empirical-validation.md) and the [Master Test Matrix](../docs/testing/test-matrix-and-ci-audit.md).

---

## 1. Pinned Hardware & Host Environment

All release benchmarks were executed under standardized host constraints (below).
Every performance number in this directory is measured on the stated host
against VarDictJava as the baseline — paper figures are never quoted as ours.
- **CPU/host:** Isolated GCP `n2-standard-32` (Intel Xeon @ 2.80 GHz, Linux x86_64, 8 worker threads / `$P = 4$` pinned physical cores `0-7` for `scale.json`, Ubuntu LTS).
- **Concurrency:** Normalized to **8 CPU threads** (`-th 8`) as a standard multi-tenant cluster baseline (not an engine scaling ceiling).
- **Storage:** Local NVMe SSD.
- **Toolchain:** Rust 1.96.0 (`x86_64-unknown-linux-musl` and `x86_64-unknown-linux-gnu`), OpenJDK 25.0.2 / OpenJDK 8 for oracle (`vardict-oracle`).

---

## 2. Five-Target Performance Budget

Evaluated via `harness/perf_budget.py`, comparing `vardictrust` against the `VarDictJava v1.8.3` pipeline (Java + R + Perl). Requirement: Rust wall time $\le$ Java wall time AND Rust peak RSS $\le$ Java peak RSS. For the 5th target (`hg008-cram-twin`), CRAM-to-BAM wall time ratio must remain $\le 2.50\times$ with 100% bitwise VCF identity.

| Workload / Target | Modality | Java Wall | Cpp Wall¹ | Rust Wall | Speedup (vs Java) | Java Peak RSS | Cpp Peak RSS¹ | Rust Peak RSS | Memory Scaling (vs Java) | Cpp vs Java Parity¹ | Evidence File |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Targeted Amplicon Panel** | AmpliconMode (`-a 10:0.95`) | `0:02.05` | `0:03.75` | **`0:00.04`** | **51.2× faster** | `118.9 MiB` | `88.9 MiB` | **`11.9 MiB`** | **10.0× less** | **[SKIP (`CPP-3`)](../docs/upstream-bugs.md#cpp-3)** | [`amplicon-tier2.json`](evidence/amplicon-tier2.json) |
| **GIAB HG002 (chr20:10M–30M, 30×)** | SimpleMode (`-f 0.01`) | `1:31.72` | `0:18.87` | **`0:08.20`** | **11.2× faster** | `4.02 GiB` | `853.1 MiB` | **`261.7 MiB`** | **15.7× less** | **[MISMATCH (`CPP-9`)](../docs/upstream-bugs.md#cpp-9)** · **[`CPP-1`](../docs/upstream-bugs.md#cpp-1)** | [`hg002-tier2.json`](evidence/hg002-tier2.json) |
| **SEQC2 HCC1395 Somatic Exome** | SomaticMode (`tumor\|normal`, `-th 4`²) | `51:10.05` | `7:11.99` | **`2:28.63`** | **20.7× faster** | `3.70 GiB` | `2.15 GiB` | **`459.0 MiB`** | **8.3× less** | **[MISMATCH (`CPP-9`)](../docs/upstream-bugs.md#cpp-9)** | [`seqc2-exome-tier2.json`](evidence/seqc2-exome-tier2.json) |
| **GIAB HG008 Somatic PDAC** | SomaticMode (`tumor\|normal`, `-th 4`²) | `1:19.70` | `0:21.96` | **`0:10.36`** | **7.7× faster** | `3.34 GiB` | `303.6 MiB` | **`226.8 MiB`** | **15.1× less** | **[MISMATCH (`CPP-9`)](../docs/upstream-bugs.md#cpp-9)** | [`hg008-tier2.json`](evidence/hg008-tier2.json) |
| **GIAB HG008 CRAM-vs-BAM Twin** | SomaticMode CRAM vs BAM (`-th 4`²) | `0:09.00` (Rust BAM) | — | **`0:19.58` (Rust CRAM)** | **2.18× CRAM/BAM** | `235.5 MiB` (BAM) | — | **`1.78 GiB` (CRAM)** | **Bounded `512 MiB` cache** | **MATCH (`18,426` Rust; `18,425` Cpp CRAM==BAM)** | [`hg008-cram-twin.json`](evidence/hg008-cram-twin.json) |

¹ Cpp (pinned **vardictcpp 3 `hc52dbad_0` bioconda**, same `vardict-oracle` env) runs beside each gate arm via `harness/perf_budget.py`. The budget verdict stays Rust-vs-Java only. Per-cell cpp wall/RSS above are transcribed from the linked evidence `extra` (`cpp_wall_s`, `cpp_maxrss_kb`, ratios). Known cpp divergences are normalized in [`docs/upstream-bugs.md`](../docs/upstream-bugs.md#cpp-divergences). A missing `cpp_time` records a loud `SKIP`, never a silent pass.
² hg008/seqc2-exome at `-th 4`: Java `-th 8` is OOM-killed on 31G hosts (bisected 1✓/4✓/8✗); exact cmds in evidence.

### Reproduction Command:
```bash
python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 hg002-tier2 amplicon-tier2 hg008-cram-twin
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
| **Chromosome 1 Scaling** | [`scale.json`](evidence/scale.json) | [Tier-2](../docs/testing/test-matrix-and-ci-audit.md#tier-2) Scaling | Full chr1 Span (249 Mbp) | $S_8 = 3.20\times$, RSS 1.34× | **PASS** |
| **Performance Budget** | [`perf-budget.json`](evidence/perf-budget.json) | Regression Gate | 5 Target Workload Regression | Wall/RAM gates | **PASS** |
| **Musl CI Packaging** | [`musl-ci.json`](evidence/musl-ci.json) | Release Gate | Static PIE Linux Packaging | Bitwise identity | **PASS** |

- **92-Chromosome Per-Contig Parity Suite (`v1.0.7`):** Full per-chromosome breakdowns across `HG008_WGS_Somatic` (`21` chrs), `SEQC2_Paired_Somatic` (`24` chrs), `SEQC2_Single_Sample_Scale` (`24` chrs), `HG002_chr20_Full_VCF` (`1` chr), and `HG002_Full_WGS_GIAB_v4.2.1` (`22` autosomes `chr1..chr22`, `8,111,171` variants, `5/5` `rtg vcfeval` files byte-identical to `VarDictJava v1.8.3`) are recorded in [`chr-parity-v107.json`](chr-parity-v107.json) and [`docs/validation/chr-parity-v107.md`](../docs/validation/chr-parity-v107.md) (**`92 / 92 MATCH`**).

---

## 4. Whole Chromosome 1 Scaling & Flat Memory Footprint

Evaluated via `harness/scale_bench.sh` on HCC1395 across all 248,956,422 bp of chromosome 1. Values below are the 3-run median certified evidence recorded in `scale.json` (`perf_fingerprint: 052f5458a6dd77bb`).

### Multi-Thread Parallel Speedup (chr1 full span, 3-run medians):
| Threads (`-th`) | Wall Clock Time | Speedup Ratio |
| :--- | :--- | :--- |
| 1 core | 3m 49s (228.78 s) | 1.00× (baseline) |
| 2 cores | 2m 03s (123.34 s) | 1.85× |
| 4 cores | 1m 17s (77.24 s) | 2.96× |
| 8 threads | 1m 11s (71.44 s) | **3.20× parallel** |

### Bounded Resident Memory vs Genomic Span (at `-th 8`, 3-run medians):
| Genomic Span | Peak RSS |
| :--- | :--- |
| 1 Mbp | 430.6 MiB (440,928 KB) |
| 10 Mbp | 997.3 MiB (1,021,188 KB) |
| 100 Mbp | 1.04 GiB (1,094,556 KB) |
| 249 Mbp (Full chr1) | **1.31 GiB (1,369,100 KB)** |

**Observed RSS Ratio (249 Mbp vs 10 Mbp): 1.34× (`1.341×`)** (Passing pre-registered threshold $\le 1.5\times$). Memory remains bounded across deep genomic spans via tiled execution with boundary halo deduplication (`SplitVarDeduper` `200 kbp` cross-tile `.var` ownership window and `SpillableVcfAccumulator` `8 MiB` disk-backed whole-stream VCF spill bound).

Pinned **vardictcpp 3 `hc52dbad_0`** runs beside each loop above (`harness/scale_bench.sh` via `harness/run_cpp.sh`, advisory only): per-thread cpp wall (`cpp_wall_s`, exome BED workload) and per-span cpp RSS (`cpp_maxrss_kb`, `-R` full-span workload at `-th 8`) land separately in `scale.json` `extra` — at `-th 8` BED wall 182.80 s, full-span (249 Mbp) unchunked `-R` RSS 24.34 GiB (25,521,616 KB per [CPP-5](../docs/upstream-bugs.md#cpp-5)).

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
cargo test -p vardictrust-core --test engine_behavior concurrency::
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

## 7. Dataset Provenance & End-to-End Cold-Start Bootstrap

All external genomic benchmark assets are tracked with cryptographic SHA-256 (and legacy MD5) checksums in [`harness/assets.tsv`](../harness/assets.tsv) and recorded upon download in [`benchmarks/MANIFEST`](MANIFEST). A bare repository clone can reproduce every reference genome, BAM/CRAM slice, truth VCF/BED, and all 16 release evidence artifacts end-to-end with the commands below:

```bash
# 1. Provision pinned oracle + C++ conda environments (vardict-oracle, vardict-cpp, g4-tools)
bash harness/setup_oracle.sh
conda env create -f harness/env/g4-tools.yml

# 2. Fetch, SHA-256-verify, unpack, index, and stage all references & biological cohorts into scratch/g4/
#    (or run individual targets: refs | amplicon | hg002 | hg008 | seqc2 | multimodal)
bash harness/stage_g4.sh all

# 3. Run full live certification (populates scratch/golden-store/ on first run and stamps benchmarks/evidence/*.json)
bash harness/recertify.sh
```

Individual cohort staging scripts invoked by [`harness/stage_g4.sh`](../harness/stage_g4.sh):
- **Public assets + SHA-256 check:** `bash harness/fetch_data.sh refs hg002 hg008 seqc2`
- **Targeted amplicon panel (`scratch/g4/amplicon/`):** `bash harness/get_amplicon_fixture.sh`
- **GIAB HG002 chr20 slice (`scratch/g4/hg002/`):** `bash harness/get_hg002_slice.sh`
- **GIAB HG008 paired PDAC slices + CRAM twins (`scratch/g4/hg008/`):** `bash harness/get_hg008_slice.sh`
- **SEQC2 WES alignment + BED derivation (`scratch/g4/seqc2/`):** `bash harness/stage_g4.sh seqc2` (calls [`harness/seqc2_align.sh`](../harness/seqc2_align.sh) and [`harness/select_slice_regions.py`](../harness/select_slice_regions.py))
- **18 Multi-modal slices (`scratch/g4/multimodal/`):** `bash harness/get_multimodal_slices.sh`
