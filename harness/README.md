# VarDictRust Differential Testing & Benchmarking Harness

This directory contains the testing, benchmarking, oracle comparison, and verification harness for **VarDictRust**, built under the [rewrites.bio](https://rewrites.bio) manifesto.

The primary mission of the harness is **uncompromising fidelity verification**: proving byte-for-byte variant equivalence and float concordance ($\varepsilon = 10^{-6}$) against **VarDictJava v1.8.3** on valid inputs, while certifying execution speedup and bounded memory usage across deep sequencing datasets.

---

## 1. Architecture Overview

The harness operates in two distinct tiers:
- **Intermediate Table Parity ([Tier 1](../docs/testing/test-matrix-and-ci-audit.md#tier-1)):** Evaluates intermediate `.var` tables generated via `--emit-intermediate` using `compare_var.py`. Validates column count integrity (Simple 36/38, Amplicon 38/40, Somatic 55/61), strand bias counters, Fisher p-values, odds ratios, and variant classifications prior to VCF serialization.
- **VCF Differential Parity ([Tier 2](../docs/testing/test-matrix-and-ci-audit.md#tier-2)):** Evaluates final VCFv4.2 output against VarDictJava v1.8.3 running in an isolated reference oracle environment (`vardict-oracle`) using `compare_vcf.py`. Aligns records by `CHROM/POS/REF/ALT`, checks genotypes, filters, and info tags, enforcing banker's half-even rounding tolerances ($\varepsilon = 10^{-6}$).

---

<a id="key-scripts--runners"></a>
## 2. Core Drivers & Comparators

| Script / Tool | Description & Usage |
| :--- | :--- |
| [`compare.sh`](compare.sh) | **Core Differential Driver:** Runs both VarDictJava and VarDictRust against a given reference FASTA, BAM, and BED, comparing final VCF outputs using `compare_vcf.py`. Supports both single BAM and somatic paired BAMs (`tumor.bam\|normal.bam`). |
| [`compare_vcf.py`](compare_vcf.py) | **VCF Differential Comparator:** Python comparator that parses and aligns records by coordinate and allele, comparing FILTER, INFO, and FORMAT fields with exact float tolerances ($\varepsilon = 10^{-6}$) and validating VCF header metadata. Supports `--strict-bytes`. |
| [`compare_var.py`](compare_var.py) | **Table Differential Comparator:** Compares raw intermediate `.var` output tables row-by-row, validating allele frequency ties and column alignment. |
| [`golden_store.sh`](golden_store.sh) | **Golden Baseline Store CLI:** Manages cached, hash-keyed oracle baseline outputs (`scratch/golden-store/`) across Java and C++ runs with atomic locking, verification, and telemetry. |
| [`run_oracle.sh`](run_oracle.sh) | **Oracle Execution Wrapper:** Executes upstream `VarDictJava v1.8.3` within the isolated `vardict-oracle` conda environment (OpenJDK 8 + Perl/R pipeline). |
| [`setup_oracle.sh`](setup_oracle.sh) | **Oracle Environment Setup:** One-time configuration script provisioning the `vardict-oracle` conda environment. |
| [`rgate.sh`](rgate.sh) | **Rapid Developer Smoke Gate:** Runs a fast $(<30\text{s})$ differential smoke-test against micro-fixtures, ideal for fast inner-loop development. |
| [`slice.sh`](slice.sh) | **BAM Slicer Utility:** Extracts local genomic intervals using `samtools` for rapid inner-loop fixture generation. |
| [`trace_compare.sh`](trace_compare.sh) | **Trace Debugger:** Fine-grained differential runner producing record-level diff traces on mismatched outputs. |

---

## 3. Release Acceptance, CI Gates & Self-Test Suite

| Script | Purpose & Gate Coverage | Runtime |
| :--- | :--- | :--- |
| [`acceptance.sh`](acceptance.sh) | **Master Acceptance Gate:** Executes the 8 mandatory release criteria defined in [`AGENTS.md` §15](../AGENTS.md#15-definition-of-done-every-change). Verifies workspace unit tests, 91/1 corpus parity, Tier-2 micro/macro parity, crash robustness, 4-target performance budget, truth concordance, hygiene, and distribution equivalence. | $<1\text{m}$ |
| [`pr_check.sh`](pr_check.sh) | **Standard PR Regression Suite:** Enforces a 6-stage verification gate on all PRs: Hygiene (`cargo fmt`, `clippy -D warnings`, `cargo audit`, `cargo deny`), Crash Robustness, 91/1 Corpus, Tier-2 Micro Parity, 4-Target Performance Budget, and Gitleaks secret scanning. | $<35\text{s}$ |
| [`run_self_tests.sh`](run_self_tests.sh) | **Diagnostic Self-Test Harness:** Discovers and executes all 24 diagnostic test scripts (16 shell tests and 8 python test suites in `harness/tests/`), reporting unified status and elapsed milliseconds. | $<10\text{s}$ |
| [`setup_hooks.sh`](setup_hooks.sh) | **Git Hooks Configuration:** Configures repository Git hook path to `.githooks/` (`git config core.hooksPath .githooks`), activating fast local pre-commit hygiene and secret checks. | $<1\text{s}$ |
| [`engine_fingerprint.sh`](engine_fingerprint.sh) | **Release Engine Fingerprinter:** Computes a 16-character SHA-256 hash across engine source code (`crates/*/src`), comparator drivers, and locked dependencies (`Cargo.lock`). Guarantees evidence records are tied to an exact, unmutated engine commit. | $<1\text{s}$ |
| [`bin_equiv.sh`](bin_equiv.sh) | **Binary Equivalence Gate:** Validates that statically-compiled `x86_64-unknown-linux-musl` release binaries produce identical variant calls to glibc binaries. | $<30\text{s}$ |

---

## 4. Benchmark Runners, Truth Evaluators & Asset Management

| Script / Asset | Description |
| :--- | :--- |
| [`assets.tsv`](assets.tsv) | **Asset Registry:** Declarative 5-column manifest (`group`, `name`, `url`, `md5`, `dest`) tracking all external public genomic references and truth sets (SEQC2, GIAB HG008, GIAB HG002). |
| [`fetch_data.sh`](fetch_data.sh) | **Safe Asset Downloader:** Reads `assets.tsv`, downloads external biological assets via `curl`, verifies MD5 checksums, and records provenance in `benchmarks/MANIFEST`. Supports `-` for dynamic calculation. |
| [`stage_g4.sh`](stage_g4.sh) | **Turnkey Cohort Stager:** Idempotent orchestrator (`all`, `refs`, `amplicon`, `hg002`, `hg008`, `seqc2`, `multimodal`) that unpacks reference FASTAs, builds `.fai` and RTG `.sdf` indices, aligns SEQC2 BAMs, derives engine/eval BEDs, and streams all benchmark slices. |
| [`convert_fastas.py`](convert_fastas.py) | **CI Reference Materializer:** Reconstructs dense reference `.fa` and `.fa.fai` files from upstream `VarDictJava/testdata/fastas/*.fa.csv` in $<10$ seconds using pure Python (zero external dependencies). |
| [`run_integration.sh`](run_integration.sh) | **Corpus Runner:** Executes all 92 upstream VarDictJava integration test cases (`VarDictJava/testdata/integrationtestcases/`). Auto-materializes FASTAs via `convert_fastas.py` if missing. Verifies actual output against pinned live-oracle SHA-256 for XFAIL fixtures. Asserts 91 PASS / 1 verified XFAIL (0 FAIL / 0 SKIP / 0 XPASS). |
| [`integration_xfail.tsv`](integration_xfail.tsv) | **Defect Registry:** Declarative 6-column TSV tracking authorized upstream fixture defects ([UB-4](../docs/upstream-bugs.md#ub-4-hard_clip_casefa-fixture-drift-between-csv-chunks-and-reconstructed-fasta)), rationale, and live-oracle SHA-256 signatures (`expected_sha256`) to prevent silent regressions. |
| [`get_hg002_slice.sh`](get_hg002_slice.sh) | **HG002 HTTP BAM Slicer:** Downloads the remote 20 Mbp euchromatic GIAB HG002 slice from public AWS S3 without fetching the 100+ GB whole-genome BAM. |
| [`get_hg008_slice.sh`](get_hg008_slice.sh) | **HG008 HTTP BAM Slicer & CRAM Twin Builder:** Streams paired Tumor (195×) and Normal (169×) GIAB HG008 slices over HTTPS, derives `slices.bed`/`engine.bed`/`eval.bed`, and encodes `HG008-{T,N}.slice.cram` twins. |
| [`run_hg002.sh`](run_hg002.sh) | **HG002 Benchmark Driver:** First-class automated runner executing HG002 Tier-2 dual-engine comparison and RTG vcfeval truth concordance. |
| [`scale_bench.sh`](scale_bench.sh) | **Scaling Benchmark Driver:** Runs whole-chromosome 1 scaling benchmarks across thread allocations (1, 2, 4, 8 threads), measuring speedup and memory linearity. |
| [`profile.sh`](profile.sh) | **Profiling Driver (issue #30):** `profile.sh {baseline,optimized} [--only <cell>] [--dry-run]` runs the Rust engine only on named cells (`tiny`, `s1a`, `s8`, `s9`, `s3`, `rna-i`, `hg002`, `hg008`, `wgs-scale`) with the `profile-cpu` / `profile-heap` / `profile-dhat` builds, each through `courtesy.sh` (and a `systemd-run` `MemoryMax=14G` scope for dhat and the large cells; dhat never runs on HG002/HG008/WGS), records `/usr/bin/time -v`, and writes `scratch/profiles/<tag>/<cell>/` plus `catalog.md`. Refuses to write under `benchmarks/evidence/`. |
| [`profile_catalog.py`](profile_catalog.py) | **Profile Catalog Ranker:** stdlib-only `profile_catalog.py <cell-dir>... -o catalog.md [--min-pct 2]`; ranks CPU self time, dhat bytes/blocks and jemalloc heap bytes by the innermost `vardictrust_core::` frame, decodes the pprof `heap.pb.gz` with a hand-written protobuf reader, writes full-stack `*-collapsed.txt` files, renders `dhat-*.svg` when `inferno-flamegraph` is on `PATH`, and links every flamegraph SVG from the catalog. |
| [`perf_budget.py`](perf_budget.py) | **Performance Budget Verifier:** Evaluates execution time and peak memory against VarDictJava baselines across 4 targets (HG008, SEQC2, HG002, Amplicon), asserting Rust $\le$ Java wall time and RSS. |
| [`run_multimodal_slices.sh`](run_multimodal_slices.sh) | **Multi-Modal Multi-Slice Parity Suite:** Executes dual-engine live oracle verification across 18 non-overlapping biological slices from 6 diverse public datasets spanning 9 distinct chromosomes. Asserts 18/18 MATCH and writes `benchmarks/evidence/multimodal-slices-tier2.json`. |
| [`get_multimodal_slices.sh`](get_multimodal_slices.sh) | **Multi-Modal Slice Streamer:** Ingests and streams 18 non-overlapping slices across 6 public datasets (COLO829, Platinum NA12878, SEQC2, GIAB HG005, Amplicon, RNA-seq) consuming $<150\text{ MB}$ disk space. |
| [`run_extended_coverage.sh`](run_extended_coverage.sh) | **Extended Modality Parity Suite:** Verifies dual-engine parity across 17 caller modalities including RNA splicing (`-i`), CRISPR cut sites (`-J/-j`), `<INV>`, `<DUP>`, and CLI flag permutations (`-u -UN`, `--adaptor`, `-k 0`, `-U`, `-K`, `-t`, `-3`, `-O`, `-x`, `-s -e`, `-M`, `--deldupvar`, `--DP ERR`). Asserts 17/17 MATCH. |
| [`concordance_multimodal.sh`](concordance_multimodal.sh) | **Multi-Modal Truth Concordance Suite:** Automated runner executing dual-engine calling and `rtg vcfeval` concordance for GIAB HG005 (GRCh38) and GIAB NA24385 (GRCh37/b37) against NIST v4.2.1 truth sets with offline caching in `scratch/g4/truth_cache/`. Asserts EQUAL (0 site diffs). |
| [`concordance.sh`](concordance.sh) | **Somatic Truth Evaluator:** Executes Real Time Genomics `rtg vcfeval` with `--squash-ploidy` for somatic variant truth comparisons. |
| [`concordance_germline.sh`](concordance_germline.sh) | **Germline Truth Evaluator:** Assesses germline variants against NIST truth, reporting dual thresholds (unfiltered discovery baseline and optimal `QUAL >= 128` operating point). |
| [`concordance_pair.sh`](concordance_pair.sh) | **Dual-Engine Truth Concordance:** Compares RTG truth concordance metrics between VarDictJava and VarDictRust, asserting 100% identical site lists. |
| [`courtesy.sh`](courtesy.sh) | **Host Etiquette Wrapper:** Cooperative multi-tenant cluster resource etiquette: limits concurrency to $\le 8$ threads, runs at `nice -n19 ionice -c3`, ensures disk/memory headroom (`MIN_FREE_GB=150`, `MIN_AVAIL_GB=12`), and guards against execution during scheduled host maintenance windows. |
| [`memwatch.sh`](memwatch.sh) | **Memory Watchdog:** Background daemon monitoring process tree RSS to prevent system OOM conditions. |

---

## 5. Slices Catalogue

The [`slices.tsv`](slices.tsv) index and [`slices/`](slices/) directory define targeted genomic intervals extracted from deep sequencing datasets to exercise specific caller behaviors:

- `s1a`, `s1b`: High-depth NA24385 exome intervals exercising single-nucleotide variant calling and HashMap iteration stability.
- `s2`: Complex indel intervals exercising local realignment and Smith-Waterman gap alignments.
- `s3`, `s4`: Paired tumor/normal exome intervals (Colo829) exercising SomaticMode classification and LOH detection.
- `s5`, `s6`: Targeted amplicon assay intervals exercising primer-trimming logic and `AMPBIAS` filtering.
- `s7`, `s8`, `s9`: Deep WGS high-depth regions exercising structural variant calling (`<DUP>`, `<INV>`, deletions).
- `tiny`: 4-variant micro interval for rapid inner-loop smoke testing.
- `data/grch38_exclude_ig_hla.bed`: Genomic mask excluding hypervariable immunoglobin and HLA loci on GRCh38.

---

## 6. Meta-Test Suite (`harness/tests/`)

The harness contains an internal self-test suite (24 test scripts) verifying that the harness tools, comparators, and zero-dependency environment contracts function correctly:

```bash
# Run all harness meta-tests using the diagnostic self-test harness
bash harness/run_self_tests.sh
```

The runner dynamically discovers and executes:
- 16 Shell Tests: `test_compare_sh.sh`, `test_concordance.sh`, `test_concordance_sv.sh`, `test_courtesy.sh`, `test_edge_cases.sh`, `test_env_pinned.sh`, `test_evidence.sh`, `test_fetch_data.sh`, `test_golden_wiring.sh`, `test_perf_budget.sh`, `test_profile_sh.sh`, `test_release_deps.sh`, `test_run_cpp.sh`, `test_run_oracle.sh`, `test_scale_bench.sh`, `test_seqc2_align.sh`.
- 8 Python Test Suites: `test_bioconda_recipe.py`, `test_compare_vcf.py` (10 tests), `test_docs_ui.py` (10 tests), `test_golden_store.py` (39 tests), `test_profile_catalog.py` (20 tests), `test_select_slice_regions.py` (4 tests), `test_stdlib_purity.py` (AST zero-PyPI import guard), `test_vcf_header.py` (1 test).

All 24 diagnostic test scripts pass 100% green.

---

## 7. Golden baseline store

The Golden Baseline Store (`harness/golden_store.sh`, implemented via `harness/golden_store.py` in pure Python stdlib) caches and replays pinned oracle baseline outputs across differential testing runs. Because VarDictJava v1.8.3 is archived and immutable (commit `0886261`) and `vardictcpp 3` is pinned to `hc52dbad_0` in `benchmarks/env/vardict-oracle.yml`, an oracle invocation on identical biological inputs and CLI arguments produces identical output. Caching these baselines eliminates redundant oracle compute, reducing local test and re-certification cycles from hours to seconds while preserving byte-for-byte and strict float differential rigor.

### 7.1 Store Layout & Storage Mechanics

By default, the store resides at `scratch/golden-store/` (configurable via `GOLDEN_STORE_DIR`). CI environments and fresh repository clones start with an empty store directory and therefore execute all oracle arms live:

```
scratch/golden-store/
├── cells/
│   └── <key>/
│       ├── manifest.json
│       ├── <n>/
│       │   ├── out[.gz]
│       │   ├── err
│       │   └── time[.gz]
│       └── .lock
├── stats.log
└── .shacache.json
```

- **`cells/<key>/`**: Unique directory per execution signature key (a 16-character hex SHA-256 digest).
- **`manifest.json`**: Cell metadata log recording the key, schema version, oracle kind (`java` or `cpp`), package pins, environment recipe SHA-256, wrapper script SHA-256, input file SHA-256 hashes, normalized argument list, environment variables, artifact metadata, host profile, UTC timestamp, and revision history.
- **`<n>/`**: Revision subdirectory (`0`, `1`, ...) containing run artifacts:
  - `out` or `out.gz`: Captured stdout / output VCF or table (automatically compressed with gzip if raw size $\ge 1\text{ MB}$).
  - `err`: Captured stderr diagnostic stream.
  - `time` or `time.gz`: Captured `/usr/bin/time -v` resource telemetry (wall clock time, peak RSS, exit code).
- **`.lock`**: Atomic, portable `mkdir`-based directory lock protecting concurrent access to the cell across worker processes. Stale locks older than 600 seconds (`LOCK_STALE_S`) are automatically reclaimed.
- **`stats.log`**: Append-only TSV event ledger tracking UTC timestamps, `GOLDEN_RUN_ID`, cache events (`HIT`, `MISS`, `REGEN`, `CORRUPT`, `FAIL`, `REFUSED`, `UNSTABLE`, `PINFAIL`, `OFF`), oracle kind, cell key, and invocation labels.
- **`.shacache.json`**: Mtime/size-keyed SHA-256 cache preventing repeated hashing of large FASTA reference genomes and BAM files (bypassed when `GOLDEN_REHASH=1`).

### 7.2 Key Derivation & Exclusions

The cell key is computed as a SHA-256 digest over everything that influences oracle execution:
1. Oracle kind (`java` or `cpp`).
2. Package pin (`vardict-java=1.8.3` or `vardictcpp=3`) and `benchmarks/env/vardict-oracle.yml` SHA-256.
3. Wrapper script SHA-256 (`harness/run_oracle.sh` or `harness/run_cpp.sh`, or `GOLDEN_ORACLE_CMD` / `GOLDEN_CPP_CMD` overrides).
4. Content SHA-256 hashes of all input files (reference FASTA, FASTA `.fai` index, BAM/CRAM files, BED intervals).
5. Normalized command-line arguments: input file paths in argv are normalized to `@IN0`, `@IN1`, ... so different worktrees or checkouts on the same machine share identical cache cells.
   - **Literal arguments:** The reference FASTA path argument (`-G <ref>`) is intentionally kept literal rather than normalized because VarDictJava embeds the reference path verbatim in the `##reference=` VCF header line, and VCF header comparators validate exact match.
   - **Thread count:** Thread arguments (`-th <N>`) remain in the key so that performance timing records (`time -v`) are never conflated across thread counts.
6. Execution environment variables: `AF_THR`, `LC_ALL`, `LANG` (and `CPP_CHUNK` for C++).

### 7.3 Operational Modes (`GOLDEN_MODE`)

The operating mode is configured via the `GOLDEN_MODE` environment variable:
- **`auto`** (default): Uses cached baseline cells when available; executes live and pins on miss.
- **`regen`**: Forces live execution of the oracle and appends a new revision `<n>/` to the cell manifest, re-pinning the baseline.
- **`off`**: Bypasses the golden store entirely. **`GOLDEN_MODE=off` is how to re-prove the oracle environment**, forcing live execution of VarDictJava and vardictcpp against all test cases. In default mode, the pre-push micro compare replays pinned Java once the store is populated, and evidence `time_ms` for test suites (such as `run_multimodal_slices.sh`) measures VarDictRust execution plus replay time.

### 7.4 Pin Rejection Conditions

An oracle execution is pinned into the store on exit code 0 *only* if all validity invariants hold. If an invariant is violated, the run exits normally with its native return code, but the output is discarded from the store under a logged rejection event:
- **`FAIL`**: Non-zero exit code or missing target output file.
- **`REFUSED`**: Implausible output, such as empty output and empty stderr, or a Java VCF output lacking the mandatory `#CHROM` header line.
- **`UNSTABLE`**: An input file's size or timestamp changed while the oracle was running.
- **`PINFAIL`**: Filesystem write failure or permissions error during artifact storage.

### 7.5 Kill-Switch, Corruption Recovery & Verification

- **Automatic Kill-Switch & Live Fallback:** When a cached cell is read during replay, `golden_store.py` verifies the SHA-256 checksum of every artifact (`out`, `err`, `time`). If any artifact is corrupted or truncated, replay is aborted, a `CORRUPT` event is logged, and the harness immediately falls back to a live oracle execution, storing a clean new revision `<n+1>/`.
- **Integrity Verification:**
  ```bash
  bash harness/golden_store.sh verify
  ```
  Scans all manifests across `scratch/golden-store/cells/` and verifies every artifact against its recorded SHA-256 hash. Exits 1 if any cell is corrupt. Because `verify` re-hashes every artifact in the store, runtime grows with total store size.
- **Remedy for Obsolete or Corrupted Cells:**
  To safely retire a cell without corrupting the append-only history or racing active workers:
  ```bash
  # 1. Locate the cell key
  bash harness/golden_store.sh list
  # 2. Tombstone the obsolete key with a documented rationale
  bash harness/golden_store.sh tombstone <key> "corrupt artifact recovered live"
  ```
  Tombstoning appends `{"tombstone": true, "reason": "..."}` to `manifest.json`. Subsequent invocations treat the cell as absent, executing live and re-pinning a fresh entry.

### 7.6 Host-Bound Performance Telemetry

Replay of `/usr/bin/time -v` resource metrics (`time` or `time.gz`) is strictly restricted to the **same host** where measurements were taken. `golden_store.py` records host hardware attributes (`cpu` model name and `nproc` core count from `/proc/cpuinfo`) in `manifest.json`. If a cell is accessed on a different machine, output artifacts (`out` and `err`) are replayed, but the time sidecar is omitted and timing falls back to live measurement. This preserves the invariant of AGENTS.md §14: deterministic variant calling outputs gate globally, but host-dependent runtime and memory metrics are never conflated across hardware architectures.

### 7.7 Telemetry, Stamping & Evidence Blocks

- Setting `GOLDEN_RUN_ID=<uuid>` groups all operations from a test run in `stats.log`.
- Benchmark runners query store statistics to embed cache provenance into committed evidence files:
  ```bash
  bash harness/golden_store.sh stats --run "$GOLDEN_RUN_ID" --json
  ```
  This emits a structured JSON block (`hit`, `live`, `events`) stored directly under `extra.golden` in `benchmarks/evidence/*.json`.
- Gate scripts verify cache hit rates using the `--min-hit` assertion:
  ```bash
  bash harness/golden_store.sh stats --run "$GOLDEN_RUN_ID" --min-hit 0.80
  ```

### 7.8 No-Waiver Gates & Scale Carry Rules

- **Zero Acceptance Waivers:** Every acceptance gate in `harness/acceptance.sh` is strictly PASS or FAIL. The acceptance gate enforces zero waivers; every gate must pass cleanly or fail loudly.
- **Scale Carry Rule:** Whole Chromosome 1 scale benchmarks (`benchmarks/evidence/scale.json`, Item 5) require substantial host compute time ($\approx 7\text{ hours}$ across 3 median runs). Benchmark evidence is carried forward between releases *only* while performance-relevant engine code remains unchanged, as validated by:
  ```bash
  bash harness/engine_fingerprint.sh --perf
  ```
  If `extra.perf_fingerprint` in `scale.json` differs from the current engine perf fingerprint, Item 5 fails and requires recertification.
- **Certifying Scale Invocation:**
  To generate certified scale benchmark evidence on qualifying hardware:
  ```bash
  EST_MIN=450 bash harness/courtesy.sh bash harness/scale_bench.sh --no-cpp <ref> <bam> <bed> <sample> scratch/g4/scale
  ```
  - Requires courtesy scheduling wrapper (`COURTESY_ACTIVE=1`, `nice -n 19`, disk headroom $\ge 150\text{ GB}$, memory $\ge 12\text{ GB}$).
  - Scheduled run start must fall between **05:15 and ~20:30** to avoid the nightly 03:40–05:15 host maintenance window.
  - Requires default core pin `SCALE_CPUS=0-7` with $P \ge 4$ distinct physical cores to satisfy the pre-registered speedup bar $\text{bar}(P) = 0.6 \cdot \min(8, P) = 2.4\times$.
- **C++ Region Arm Notice:** In `scale_bench.sh`, the C++ region scaling arm is executed as a BED run through `run_cpp.sh` (converting `-R chr:start-end` into an equivalent 0-based BED interval). Consequently, `cpp_*` numbers in `scale.json` are not directly comparable to pre-#31 measurements.
