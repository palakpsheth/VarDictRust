# Feature Specification: VarDictRust Port

**Feature Branch**: `001-vardict-rust-port`

**Created**: 2026-06-19

**Status**: Approved

**Input**: User description: "Rewrite VarDictJava v1.8.3 (a somatic/germline NGS variant caller) from scratch in idiomatic Rust as a single, statically-compiled binary that scales to deep WGS and produces output identical to VarDictJava v1.8.3 on valid inputs."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Single-Sample (Germline) Variant Calling (Priority: P1)
As a genomic analyst, I want to run variant calling on a single BAM/CRAM sample using a reference genome and a target BED file of interest, so that I can call germline variants with exact output parity to VarDictJava's Simple mode.

**Why this priority**: This represents the foundational MVP. Somatic and amplicon calling pipelines build upon the core single-sample calling infrastructure.

**Independent Test**: Fully tested by executing `vardictrust` in germline mode against the single-sample test cases, comparing the intermediate output with `harness/compare_var.py` and final VCF output with `harness/compare.sh`.

**Acceptance Scenarios**:
1. **Given** a valid single-sample BAM/CRAM file, a reference FASTA, and a target region BED file, **When** `vardictrust` is executed, **Then** it produces a VCF file matching the oracle's output exactly, and exits with code 0.
2. **Given** a missing reference FASTA index, **When** `vardictrust` is executed, **Then** it exits immediately with a descriptive error code instead of crashing.

---

### User Story 2 - Somatic (Paired Tumor/Normal) Variant Calling (Priority: P2)
As an oncology researcher, I want to perform somatic variant calling on paired tumor and normal BAM/CRAM samples, so that I can distinguish somatic mutations from germline variants and sequencing noise with exact output parity to VarDictJava's Somatic mode.

**Why this priority**: Crucial oncology variant calling capability required to support clinical research pipelines.

**Independent Test**: Validated by running paired tumor/normal sample integration test cases and verifying somatic STATUS and statistical outputs match the oracle.

**Acceptance Scenarios**:
1. **Given** paired tumor and normal BAM/CRAM files, a reference FASTA, and a target BED, **When** `vardictrust` is run in somatic mode (`-b "tumor.bam|normal.bam"`), **Then** it produces a somatic VCF matching the oracle's VCF and intermediate `.var` outputs (55/61 columns).
2. **Given** mismatched sample names or BAM counts, **When** `vardictrust` is executed, **Then** it errors loudly without panic.

---

### User Story 3 - Amplicon Mode and Splicing/SV Support (Priority: P3)
As a genomic researcher, I want to process PCR-based amplicon sequencing datasets or detect structural variants (SVs) and splicing variants, so that I can achieve full functional parity with all modes of VarDictJava v1.8.3.

**Why this priority**: Required for comprehensive feature completion of the v1.8.3 specification rewrite.

**Independent Test**: Validated by running amplicon and splicing/SV integration test cases, comparing outputs against golden oracle files.

**Acceptance Scenarios**:
1. **Given** an amplicon BAM file, a reference FASTA, and a target BED, **When** `vardictrust` is run in amplicon mode (`-a`), **Then** it produces variant calls matching the oracle's Amplicon mode outputs (38/40 columns).

---

### Edge Cases

- **Oracle Nondeterminism**: VarDictJava has known run-to-run output order nondeterminism (issue #282). The engine must maintain a deterministic ordering or align records by chromosome/position/ref/alt so that VCF comparative tests do not fail due to sorting alone.
- **Complex CIGAR Boundaries**: Reads containing complex insertions, deletions, or clipping at tile/region boundaries must be processed without causing out-of-bounds or null pointer exceptions (fixing the `modifyCigar` crash family).
- **Extreme Sequencing Depth**: High-depth genomic regions (e.g., >100,000x coverage) must be processed within a bounded memory tile, ensuring the system does not run out of memory (OOM).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST output VCF format matching VarDictJava v1.8.3 structurally and content-wise.
- **FR-002**: System MUST support all four VarDict modes: Simple, Somatic, Amplicon, and Splicing/SV.
- **FR-003**: System MUST execute local realignment matching the oracle's VariationRealigner heuristics.
- **FR-004**: System MUST perform two-sided Fisher's Exact test and Brent root-finding odds ratio MLE calculation to replicate `--fisher` outputs within $\epsilon = 10^{-5}$ tolerance.
- **FR-005**: System MUST compute somatic STATUS (somatic, germline, LOH, noise, etc.) matching the Java SomaticPostProcessModule logic.
- **FR-006**: System MUST utilize Noodles for pure Rust BAM/CRAM/FASTA/VCF file reading and writing.
- **FR-007**: System MUST support parallel region tiling with work-stealing while ensuring deterministic ordered output.
- **FR-008**: System MUST forbid unsafe code and panics in the core engine crate, returning proper `Result` types.

### Key Entities *(include if feature involves data)*

- **Alignment Read**: Represents a mapped sequencing read with fields for chromosome, start position, CIGAR, sequence, quality, and flags.
- **Region Tile**: Represents a genomic slice with a coordinate interval and a boundary halo for local sequencing coverage.
- **Variant Call**: A parsed genomic variant containing genomic position, reference allele, alternate allele, allele frequency, and calculated statistics (depth, qualities, SBF, SSF, Fisher test scores).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of integration test cases in the test corpus produce VCF and `.var` output matching the VarDictJava oracle (Tier-1 and Tier-2 tests report MATCH).
- **SC-002**: Float differences for Fisher p-values, odds ratios, and allele frequencies MUST be within $\epsilon = 10^{-5}$ after Banker's rounding.
- **SC-003**: Average execution speed on WGS datasets is at least 3x faster than VarDictJava v1.8.3 under equivalent thread counts.
- **SC-004**: Peak RAM usage remains below 4GB per thread regardless of sequencing depth or chromosome size.

## Assumptions

- Input BAM and CRAM files are well-formed and index files (`.bai` / `.crai`) exist.
- Reference genomes have matching `.fai` index files.
- The reference oracle VarDictJava v1.8.3 is available for differential testing.
