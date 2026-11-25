# VarDictRust Documentation Portal

Welcome to the **VarDictRust** documentation portal.

VarDictRust is a fast, memory-safe Rust somatic and germline NGS variant caller that delivers 100% byte-identical output to **VarDictJava v1.8.3** by default (`--compat=vardictjava` / `--parity`)—with opt-in `v1.1+` `--compat=off` (`CompatProfile`, 21 correctness switches) and `--vcf-all-alleles` multi-allelic VCF emission—while consolidating the Java→R→Perl pipeline (plus native multithreaded BGZF and `.tbi`/`.csi` indexing) into a single statically-compiled binary that scales to deep whole-genome sequencing (WGS). Built under the **[rewrites.bio](https://rewrites.bio)** manifesto.

---

## Quick Navigation

### 1. Getting Started & User Guide
- **[README](https://github.com/palakpsheth/VarDictRust#readme)** — Executive overview, performance scoreboards, quickstart commands, and installation.
- **[Installation Guide](https://github.com/palakpsheth/VarDictRust#installation)** — Pre-built static/native binaries (`cargo-dist` SLSA Level 3), Cargo install, Seqera Wave containers (Docker/Podman/Apptainer), Conda/Bioconda (`conda`/`mamba`/`pixi`), and Nextflow (`nf-core/modules`).
- **[Usage Recipes & Modalities](https://github.com/palakpsheth/VarDictRust#usage-recipes)** — Diátaxis recipes for Germline WGS, Somatic Tumor-Normal, Targeted Amplicon (`-a` / `-sp`), RNA-seq Splicing (`-i`), UNIX/BGZF streaming pipelines (`--output`, `--tbi`/`--csi`), and `v1.1+` `--compat=off` + `--vcf-all-alleles` calling.
- **[Migration Guide](https://github.com/palakpsheth/VarDictRust/blob/main/docs/migration.md)** — Side-by-side Rosetta Stone from `VarDictJava` + R + Perl (`var2vcf_valid.pl` / `var2vcf_paired.pl`) pipelines to single-binary `vardictrust`.
- **[CLI Parameter Reference](https://github.com/palakpsheth/VarDictRust#cli-parameter-reference)** — Complete 91-flag reference across 7 functional groups, including `--compat=vardictjava|off`, `--amplicon-single-primer` (`-sp`), all 18/18 consolidated `var2vcf_*.pl` post-filters (`--vcf-all-alleles`), and native BGZF/`.tbi`/`.csi` output flags.
- **[Contributing Guidelines & Release Runbook](https://github.com/palakpsheth/VarDictRust/blob/main/CONTRIBUTING.md#release-process)** — Development prerequisites, code style, pre-commit hygiene, test harnesses, and release process.

### 2. Empirical Validation & Parity
- **[Omnibus Empirical Validation](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md)** — Comprehensive 9-pillar validation manuscript covering the Master 10-Cohort Roster, NIST v4.2.1 HG002 benchmark, 17 extended modalities, and 18 multi-modal slices.
- **[92-Chromosome Parity Suite (v1.0.7)](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/chr-parity-v107.md)** — Per-chromosome parity receipt (`92 / 92 MATCH`, 0 diffs) across 5 cohorts, including full 22-autosome `GIAB HG002` WGS (`8,111,171` variants, byte-identical `rtg vcfeval` vs `VarDictJava v1.8.3`).
- **[Ground Truth Concordance](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#pillar-1)** — Precision, recall, and F1 concordance benchmarks evaluated with RTG `vcfeval` (HG002, HG005, NA24385 b37, HG008, SEQC2).
- **[Extended Modalities & SVs](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#pillar-8)** — RNA-seq splicing (`-i`), CRISPR (`-J/-j`), `<INV>`, `<DUP>`, and operational flags.
- **[Multi-Modal Public Cohorts](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#pillar-9)** — 18 slices across 6 public sequencing cohorts (HG005, NA24385 b37, NA12878, COLO829, HCC2218, SEQC2).
- **[Benchmark Evidence Ledger](https://github.com/palakpsheth/VarDictRust/blob/main/benchmarks/README.md)** — 16 release evidence artifacts, hardware specifications, and scaling metrics.

### 3. Architecture & Methodology
- **[Divergences & Bug Tracking](https://github.com/palakpsheth/VarDictRust/blob/main/DIVERGENCES.md)** — Categorized audit of VarDictJava quirks, edge cases, and deterministic HashMap emulation.
- **[Upstream Bug Catalog](https://github.com/palakpsheth/VarDictRust/blob/main/docs/upstream-bugs.md)** — Consolidated catalog of upstream VarDictJava v1.8.3 defects preserved by default (`--compat=vardictjava`: `UB-1` through `UB-7` output-affecting, plus environmental `UB-8`) and their opt-in `--compat=off` / `--vcf-all-alleles` resolutions.
- **[Ecosystem Survey](https://github.com/palakpsheth/VarDictRust/blob/main/docs/ecosystem-survey.md)** — Comparative survey of modern NGS variant callers (Mutect2, Strelka2, Octopus, VarDictJava) and rewrites.bio principles.

### 4. Testing, Harness & Quality Assurance
- **[Comprehensive Test Matrix & CI Audit](https://github.com/palakpsheth/VarDictRust/blob/main/docs/testing/test-matrix-and-ci-audit.md)** — 7-tier test architecture, 773 cargo unit/integration tests, and CI workflows.
- **[Test Harness Documentation](https://github.com/palakpsheth/VarDictRust/blob/main/harness/README.md)** — Architecture of Tier-1/Tier-2 comparators, 37 self-tests, and cluster etiquette guidelines.
- **[v1.0.10 Release PR Specification](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.10-pr.md)** — Release manifest, cryptographic fingerprints, and certification (`v1.0.9` in [`docs/release/v1.0.9-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.9-pr.md), `v1.0.8` in [`docs/release/v1.0.8-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.8-pr.md), `v1.0.7` in [`docs/release/v1.0.7-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.7-pr.md), `v1.0.6` in [`docs/release/v1.0.6-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.6-pr.md), `v1.0.5` in [`docs/release/v1.0.5-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.5-pr.md), `v1.0.4` in [`docs/release/v1.0.4-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.4-pr.md), `v1.0.3` in [`docs/release/v1.0.3-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.3-pr.md), `v1.0.1` in [`docs/release/v1.0.1-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.1-pr.md), `v1.0.0` in [`docs/release/v1.0.0-pr.md`](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.0-pr.md)).

### 5. Project Roadmap & Governance
- **[ROADMAP](https://github.com/palakpsheth/VarDictRust/blob/main/ROADMAP.md)** — Governance for `v1.1+` opt-in `--compat=vardictjava|off` (`CompatProfile`: all 21 `v1.1.0` accuracy, SV/mate-pileup, somatic multi-allelic/amplicon, and SAM `=`/`X` switches), `--vcf-all-alleles`, `--amplicon-single-primer` (`-sp`), and planned `v2.0` UMI deduplication, orientation-bias modeling, and long-read support.
- **[Universal Agent Contract](https://github.com/palakpsheth/VarDictRust/blob/main/AGENTS.md)** — Engineering rules, coordinate indexing, no-panic mandate, definition of done, and repository Agent Skills ([`vardictrust-issue`](https://github.com/palakpsheth/VarDictRust/blob/main/.agents/skills/vardictrust-issue/SKILL.md), [`vardictrust-pr`](https://github.com/palakpsheth/VarDictRust/blob/main/.agents/skills/vardictrust-pr/SKILL.md), [`vardictrust-profiling`](https://github.com/palakpsheth/VarDictRust/blob/main/.agents/skills/vardictrust-profiling/SKILL.md), [`vardictrust-release`](https://github.com/palakpsheth/VarDictRust/blob/main/.agents/skills/vardictrust-release/SKILL.md)).
