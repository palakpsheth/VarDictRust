# VarDictRust Differential Testing & Benchmarking Harness

This directory contains the testing, benchmarking, oracle comparison, and verification harness for **VarDictRust**, built under the [rewrites.bio](https://rewrites.bio) manifesto.

The primary mission of the harness is **uncompromising fidelity verification**: proving byte-for-byte variant equivalence and float concordance ($\varepsilon = 10^{-6}$) against **VarDictJava v1.8.3** on valid inputs, while certifying execution speedup and bounded memory usage across deep sequencing datasets.

---

## 1. Architecture Overview

The harness operates in two distinct tiers:
- **Tier 1 (Intermediate Representation):** Evaluates intermediate `.var` tables generated via `--emit-intermediate`. Validates column count integrity (Simple 36/38, Amplicon 38/40, Somatic 55/61), strand bias counters, Fisher p-values, odds ratios, and variant classifications prior to VCF serialization.
- **Tier 2 (Full VCF Parity):** Evaluates final VCFv4.2 output against VarDictJava v1.8.3 running in an isolated reference oracle environment (`vardict-oracle`). Aligns records by `CHROM/POS/REF/ALT`, checks genotypes, filters, and info tags, enforcing banker's half-even rounding tolerances ($\varepsilon = 10^{-6}$).

---

## 2. Core Drivers & Comparators

| Script / Tool | Description & Usage |
| :--- | :--- |
| [`compare.sh`](compare.sh) | **Core Differential Driver:** Runs both VarDictJava and VarDictRust against a given reference FASTA, BAM, and BED, comparing final VCF outputs using `compare_vcf.py`. Supports both single BAM and somatic paired BAMs (`tumor.bam\|normal.bam`). |
| [`compare_vcf.py`](compare_vcf.py) | **VCF Differential Comparator:** Python comparator that parses and aligns records by coordinate and allele, comparing FILTER, INFO, and FORMAT fields with exact float tolerances ($\varepsilon = 10^{-6}$) and validating VCF header metadata. |
| [`compare_var.py`](compare_var.py) | **Table Differential Comparator:** Compares raw intermediate `.var` output tables row-by-row, validating allele frequency ties and column alignment. |
| [`run_oracle.sh`](run_oracle.sh) | **Oracle Execution Wrapper:** Executes upstream `VarDictJava v1.8.3` within the isolated `vardict-oracle` conda environment (OpenJDK 8 + Perl/R pipeline). |
| [`rgate.sh`](rgate.sh) | **Rapid Developer Smoke Gate:** Runs a fast $(<30\text{s})$ differential smoke-test against micro-fixtures, ideal for fast inner-loop development. |
| [`trace_compare.sh`](trace_compare.sh) | **Trace Debugger:** Fine-grained differential runner producing record-level diff traces on mismatched outputs. |

---

## 3. Release Acceptance & CI Gates

| Script | Purpose & Gate Coverage |
| :--- | :--- |
| [`acceptance.sh`](acceptance.sh) | **Master Acceptance Gate:** Executes the 8 mandatory release criteria defined in `PLAN.md` §15. Verifies workspace unit tests, 91/1 corpus parity, Tier-2 micro/macro parity, crash robustness, 4-target performance budget, truth concordance, hygiene, and distribution equivalence. |
| [`pr_check.sh`](pr_check.sh) | **Standard PR Regression Suite:** Enforces a 6-stage verification gate on all PRs: Hygiene (`cargo fmt`, `clippy -D warnings`, `cargo audit`, `cargo deny`), Crash Robustness, 91/1 Corpus, Tier-2 Micro Parity, 4-Target Performance Budget, and Gitleaks secret scanning. |
| [`setup_hooks.sh`](setup_hooks.sh) | **Git Hooks Configuration:** Configures repository Git hook path to `.githooks/` (`git config core.hooksPath .githooks`), activating fast local pre-commit hygiene and secret checks. |
| [`engine_fingerprint.sh`](engine_fingerprint.sh) | **Release Engine Fingerprinter:** Computes a 16-character SHA-256 hash across engine source code (`crates/*/src`), comparator drivers, and locked dependencies (`Cargo.lock`). Guarantees evidence records are tied to an exact, unmutated engine commit. |
| [`bin_equiv.sh`](bin_equiv.sh) | **Binary Equivalence Gate:** Validates that statically-compiled `x86_64-unknown-linux-musl` release binaries produce identical variant calls to glibc binaries. |

---

## 4. Benchmark Runners & Asset Management

| Script / Asset | Description |
| :--- | :--- |
| [`assets.tsv`](assets.tsv) | **Asset Registry:** Declarative 5-column manifest (`group`, `name`, `url`, `md5`, `dest`) tracking all external public genomic references and truth sets (SEQC2, GIAB HG008, GIAB HG002). |
| [`fetch_data.sh`](fetch_data.sh) | **Safe Asset Downloader:** Reads `assets.tsv`, downloads external biological assets via `curl`, verifies MD5 checksums, and records provenance in `benchmarks/MANIFEST`. |
| [`convert_fastas.py`](convert_fastas.py) | **CI Reference Materializer:** Reconstructs dense reference `.fa` and `.fa.fai` files from upstream `VarDictJava/testdata/fastas/*.fa.csv` in $<10$ seconds using pure Python (zero external dependencies). |
| [`run_integration.sh`](run_integration.sh) | **Corpus Runner:** Executes all 92 upstream VarDictJava integration test cases (`VarDictJava/testdata/integrationtestcases/`). Auto-materializes FASTAs via `convert_fastas.py` if missing. Asserts 91 PASS / 1 documented stale-fixture FAIL. |
| [`get_hg002_slice.sh`](get_hg002_slice.sh) | **HTTP BAM Slicer:** Downloads the remote 20 Mbp euchromatic GIAB HG002 slice from public AWS S3 without fetching the 100+ GB whole-genome BAM. |
| [`run_hg002.sh`](run_hg002.sh) | **HG002 Benchmark Driver:** First-class automated runner executing HG002 Tier-2 dual-engine comparison and RTG vcfeval truth concordance. |
| [`scale_bench.sh`](scale_bench.sh) | **Scaling Benchmark Driver:** Runs whole-chromosome 1 scaling benchmarks across thread allocations (1, 2, 4, 8 threads), measuring speedup and memory linearity. |
| [`perf_budget.py`](perf_budget.py) | **Performance Budget Verifier:** Evaluates execution time and peak memory against VarDictJava baselines across 4 targets (HG008, SEQC2, HG002, Amplicon), asserting Rust $\le$ Java wall time and RSS. |
| [`concordance.sh`](concordance.sh) | **Somatic Truth Evaluator:** Executes Real Time Genomics `rtg vcfeval` with `--squash-ploidy` for somatic variant truth comparisons. |
| [`concordance_germline.sh`](concordance_germline.sh) | **Germline Truth Evaluator:** Assesses germline variants against NIST truth, reporting dual thresholds (unfiltered discovery baseline and optimal `QUAL >= 128` operating point). |
| [`concordance_pair.sh`](concordance_pair.sh) | **Dual-Engine Truth Concordance:** Compares RTG truth concordance metrics between VarDictJava and VarDictRust, asserting 100% identical site lists. |
| [`courtesy.sh`](courtesy.sh) | **Host Etiquette Wrapper:** Enforces Host Policy A3: limits concurrency to $\le 8$ threads, runs at `nice -n19 ionice -c3`, and guards against execution during maintenance windows. |
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

The harness contains an internal self-test suite (13 test scripts) verifying that the harness tools and comparators themselves function correctly:

```bash
# Run all harness meta-tests
for t in harness/tests/test_*.sh; do bash "$t"; done
python3 -m unittest harness/tests/test_*.py
```

All 13 test scripts pass 100% green.
