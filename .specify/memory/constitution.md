<!--
[SYNC IMPACT REPORT]
Version change: 0.0.0-template -> 1.0.0
List of modified principles:
  - [PRINCIPLE_1_NAME] -> I. Exact Output Parity (The Prime Directive)
  - [PRINCIPLE_2_NAME] -> II. Strict Two-Tier Validation
  - [PRINCIPLE_3_NAME] -> III. Library-First & No-Panic Rust
  - [PRINCIPLE_4_NAME] -> IV. Scale-First Tiled Engine
  - [PRINCIPLE_5_NAME] -> V. Coordinate and Data Fidelity
Added sections:
  - Core Principles (populated from templates)
  - Approved Technology Stack and Dependencies (replaces Section 2)
  - Development Workflow and Git Discipline (replaces Section 3)
  - Governance (populated)
Removed sections: None
Templates requiring updates:
  - .specify/templates/plan-template.md (✅ updated/aligned)
  - .specify/templates/spec-template.md (✅ updated/aligned)
  - .specify/templates/tasks-template.md (✅ updated/aligned)
  - .specify/templates/checklist-template.md (✅ updated/aligned)
Follow-up TODOs: None (all placeholders populated).
-->

# VarDictRust Constitution

## Core Principles

### I. Exact Output Parity (The Prime Directive)
The VarDictRust implementation MUST produce output that is byte-equivalent and structurally identical
to VarDictJava v1.8.3 on valid inputs. This includes the exact variant set (chrom/pos/ref/alt),
every FILTER decision, genotype format/values, somatic STATUS, and intermediate column counts
(36/38 for Simple, 38/40 for Amplicon, 55/61 for Somatic). Any oracle output bugs MUST be replicated
in v1.0 and tracked as issues, while oracle crash behaviors MUST be handled gracefully to produce the
equivalent non-crashing output. Output-changing bug fixes or enhancements are strictly deferred
post-1.0 as opt-in configurations.

### II. Strict Two-Tier Validation
Never claim a task is completed without verifiable proof of output equivalence. Validation requires
both running the Cargo test suite and performing two-tier differential testing: Tier-1 intermediate
`.var` comparison (`harness/compare_var.py`) to verify engine output, and Tier-2 final VCF comparison
(`harness/compare.sh`) to verify full pipeline parity. Faking a pass by adjusting float tolerances,
narrowing compared columns, or editing oracle expected outputs is strictly forbidden.

### III. Library-First & No-Panic Rust
All core genomics, alignment, and variant calling logic MUST reside within the library crate
`crates/vardict-core`, which enforces `#![forbid(unsafe_code)]`. Library code MUST NOT panic: use of
`unwrap`, `expect`, `panic!`, `unreachable!`, and bare slice indexing `x[i]` is forbidden.
All errors MUST be handled using Rust's `Result` type paired with `thiserror` for clean propagation,
ensuring CLI robustness.

### IV. Scale-First Tiled Engine
To support deep whole-genome sequencing (WGS), the engine MUST process reads in region tiles using a
halo strategy. This ensures memory usage is bounded regardless of sequencing depth. Loading entire
chromosomes into memory is prohibited. rayon-based work-stealing MUST be implemented to parallelize
computation across tiles while ensuring that output variants are deterministically sorted and identical
to single-region runs.

### V. Coordinate and Data Fidelity
Data representation MUST match the oracle's coordinate systems and numeric formatting. Coordinate basis
MUST be explicitly typed: `Pos0` for 0-based noodles BAM positions, and `Pos1` for 1-based VCF outputs,
converting only via explicit conversion methods. To prevent non-deterministic output ordering,
order-preserving collections (`IndexMap` or `BTreeMap`) MUST be used in realignment and print loops.
Float formatting MUST implement R's 5-decimal `round_half_even` banker's rounding to maintain parity.

## Approved Technology Stack and Dependencies

- **Rust Crate Ecosystem**: Pure Rust libraries are preferred. The project depends on `noodles` (for
  BAM/CRAM/FASTA/VCF/BGZF I/O), `bio` (for alignment/sequence utilities), `block-aligner` (for AVX2-accelerated
  realignment), `statrs` and `fishers_exact` (for stats validations), `rayon` (for parallelism), `clap` (for CLI
  args), `indexmap` (for order preservation), and `thiserror`/`anyhow` (for error handling).
- **Prohibited Libraries**: Out-of-tree or hallucinated dependencies such as `genomicframe-core`,
  `rustalign-simd`, `ksw2rs`, `seqair`, and `cyanea-align` MUST NOT be used.
- **Custom Statistics**: Brent root-finding MLE for Fisher odds ratio calculation MUST be custom
  implemented in Rust (matching the Java oracle's `FisherExact.java` and `UnirootZeroIn.java`), rather
  than relying on crates that only supply p-values.

## Development Workflow and Git Discipline

- **Surgical Port Loop**: Developers MUST follow the port loop: pick exactly one VarDictJava function,
  read and document its implementation (citing the oracle source `file:lines` in doc-comments), write
  the Rust implementation + unit tests, run differential testing (`harness/compare.sh`), and commit
  only upon a verified `MATCH`.
- **Quality Gates**: All checks MUST pass prior to committing or raising PRs. The required execution order
  is: `cargo fmt --all` -> `cargo clippy --all-targets --all-features -- -D warnings` -> `cargo test`.
- **Commit Standards**: Commits MUST be atomic, referencing the relevant `PLAN.md` task ID (e.g. `Txxx: port
  CigarModifier`). Include `Co-Authored-By: Claude <noreply@anthropic.com>` in metadata. No large data,
  condenser environments, or target artifacts may be committed.

## Governance

- The principles outlined in this Constitution are non-negotiable and supersede all other developer habits.
- Verification of compliance with this Constitution and `CLAUDE.md` is mandatory for all Pull Request reviews.
- Changes to the principles or stack requirements require a formal amendment procedure, which includes bumping
  the Constitution version and updating the Sync Impact Report.
- Development must be guided by `CLAUDE.md` and `PLAN.md` at runtime.

**Version**: 1.0.0 | **Ratified**: 2026-06-19 | **Last Amended**: 2026-06-19
