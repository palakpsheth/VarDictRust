# VarDictRust Documentation Portal

Welcome to the **VarDictRust** documentation portal.

VarDictRust is an idiomatic Rust rewrite of **VarDictJava v1.8.3** (a somatic and germline variant caller for NGS), consolidating its fragmented Java→R→Perl pipeline into a single, statically-compiled, memory-safe binary that scales to deep whole-genome sequencing (WGS). Built under the **[rewrites.bio](https://rewrites.bio)** manifesto.

---

## Quick Navigation

### 1. Getting Started & User Guide
- **[README](../../README.md)** — Executive overview, performance scoreboards, quickstart commands, and installation.
- **[Installation Guide](../../README.md#installation)** — Pre-built static musl binaries, cargo install, containers (Docker/Podman/Singularity), and conda/bioconda instructions.
- **[Usage Recipes & Modalities](../../README.md#usage-recipes)** — Diátaxis recipes for Germline WGS, Somatic Tumor-Normal, Targeted Amplicon (`-a`), RNA-seq Splicing (`-i`), and UNIX streaming pipelines.
- **[CLI Parameter Reference](../../README.md#cli-parameter-reference)** — Complete 60-flag reference across 6 functional groups.
- **[Contributing Guidelines](../../CONTRIBUTING.md)** — Development prerequisites, code style, pre-commit hygiene, and test harnesses.

### 2. Empirical Validation & Parity
- **[Omnibus Empirical Validation](../validation/omnibus-empirical-validation.md)** — Comprehensive 9-pillar validation manuscript covering the Master 10-Cohort Roster, NIST v4.2.1 HG002 benchmark, 9 extended modalities, and 18 multi-modal slices.
- **[Ground Truth Concordance](../validation/omnibus-empirical-validation.md#pillar-1)** — Precision, recall, and F1 concordance benchmarks evaluated with RTG `vcfeval` (HG002, HG005, NA24385 b37, HG008, SEQC2).
- **[Extended Modalities & SVs](../validation/omnibus-empirical-validation.md#pillar-8)** — RNA-seq splicing (`-i`), CRISPR (`-J/-j`), `<INV>`, `<DUP>`, and operational flags.
- **[Multi-Modal Public Cohorts](../validation/omnibus-empirical-validation.md#pillar-9)** — 18 slices across 6 public sequencing cohorts (HG005, NA24385 b37, NA12878, COLO829, HCC2218, SEQC2).
- **[Benchmark Evidence Ledger](../../benchmarks/README.md)** — 16 release evidence artifacts, hardware specifications, and scaling metrics.

### 3. Architecture & Methodology
- **[Divergences & Bug Tracking](../../DIVERGENCES.md)** — Categorized audit of VarDictJava quirks, edge cases, and deterministic HashMap emulation.
- **[Upstream Bug Catalog](../upstream-bugs.md)** — Consolidated catalog of upstream VarDictJava v1.8.3 defects preserved for parity (UB-1 through UB-7).
- **[Ecosystem Survey](../ecosystem-survey.md)** — Comparative survey of modern NGS variant callers (Mutect2, Strelka2, Octopus, VarDictJava) and rewrites.bio principles.

### 4. Testing, Harness & Quality Assurance
- **[Comprehensive Test Matrix & CI Audit](../testing/test-matrix-and-ci-audit.md)** — 7-tier test architecture, 208 cargo unit/integration tests, and CI workflows.
- **[Test Harness Documentation](../../harness/README.md)** — Architecture of Tier-1/Tier-2 comparators, 15 self-tests, and cluster etiquette guidelines.
- **[Release Waivers & Parity Audit](../release/waivers.md)** — Formal justification of intentional behavioral divergences (0 active waivers).
- **[v1.0.0 Release PR Specification](../release/v1.0.0-pr.md)** — Release manifest, cryptographic fingerprints, and certification.

### 5. Project Roadmap & Governance
- **[ROADMAP](../../ROADMAP.md)** — Governance for post-1.0 opt-in `--fix` flags, UMI deduplication, and long-read support.
- **[Universal Agent Contract](../../AGENTS.md)** — Engineering rules, coordinate indexing, no-panic mandate, and definition of done.
