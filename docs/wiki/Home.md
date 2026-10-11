# VarDictRust Documentation Portal

Welcome to the **VarDictRust** documentation portal.

VarDictRust is an idiomatic Rust rewrite of **VarDictJava v1.8.3** (a somatic and germline variant caller for NGS), consolidating its fragmented Java→R→Perl pipeline into a single, statically-compiled, memory-safe binary that scales to deep whole-genome sequencing (WGS). Built under the **[rewrites.bio](https://rewrites.bio)** manifesto.

---

## Quick Navigation

### 1. Getting Started
- **[README](https://github.com/palakpsheth/VarDictRust#readme)** — Executive overview, performance scoreboards, quickstart commands, and installation.
- **[Installation Guide](https://github.com/palakpsheth/VarDictRust#installation)** — Pre-built static musl binaries, cargo install, and verification.
- **[Contributing Guidelines](https://github.com/palakpsheth/VarDictRust/blob/main/CONTRIBUTING.md)** — Development prerequisites, code style, and test harnesses.

### 2. Empirical Validation & Parity
- **[Omnibus Empirical Validation](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md)** — Comprehensive 9-pillar validation manuscript covering the Master 10-Cohort Roster, NIST v4.2.1 HG002 benchmark, 9 extended modalities, and 18 multi-modal slices.
- **[Ground Truth Concordance](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#2-pillar-1-full-nist-v421-hg002-benchmark-evaluation)** — Precision, recall, and F1 concordance benchmarks evaluated with RTG `vcfeval` (HG002, HG005, NA24385 b37, HG008, SEQC2).
- **[Extended Modalities & SVs](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#8-pillar-8-extended-caller-modality-parity-suite-9-modes)** — RNA-seq splicing, CRISPR, `<INV>`, `<DUP>`, and operational flags.
- **[Multi-Modal Public Cohorts](https://github.com/palakpsheth/VarDictRust/blob/main/docs/validation/omnibus-empirical-validation.md#9-pillar-9-multi-modal-multi-slice-parity-expansion-18-slices-across-6-cohorts)** — 18 slices across 6 public sequencing cohorts (HG005, NA24385 b37, NA12878, COLO829, HCC2218, SEQC2).

### 3. Architecture & Methodology
- **[Divergences & Bug Tracking](https://github.com/palakpsheth/VarDictRust/blob/main/DIVERGENCES.md)** — Categorized audit of VarDictJava quirks, edge cases, and deterministic HashMap emulation.
- **[Upstream Bugs](https://github.com/palakpsheth/VarDictRust/blob/main/docs/upstream-bugs.md)** — Consolidated catalog of upstream VarDictJava v1.8.3 defects preserved for parity.
- **[Ecosystem Survey](https://github.com/palakpsheth/VarDictRust/blob/main/docs/ecosystem-survey.md)** — Comparative survey of modern NGS variant callers.

### 4. Testing & Quality Assurance
- **[Comprehensive Test Matrix & CI Audit](https://github.com/palakpsheth/VarDictRust/blob/main/docs/testing/test-matrix-and-ci-audit.md)** — 7-tier test architecture, 208 cargo unit/integration tests, and CI workflows.
- **[Release Waivers & Parity Audit](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/waivers.md)** — Formal justification of intentional behavioral divergences.
- **[v1.0.0 Release PR Specification](https://github.com/palakpsheth/VarDictRust/blob/main/docs/release/v1.0.0-pr.md)** — Release manifest, cryptographic fingerprints, and certification.

### 5. Project Roadmap
- **[ROADMAP](https://github.com/palakpsheth/VarDictRust/blob/main/ROADMAP.md)** — Governance for post-1.0 opt-in `--fix` flags, UMI deduplication, and long-read support.
