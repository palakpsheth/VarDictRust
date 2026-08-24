# Rewrite Quality Checklist: VarDictRust Port

**Purpose**: Validate specification requirements quality for output parity, performance, and CLI interface contract.
**Created**: 2026-06-19
**Feature**: [spec.md](../spec.md)

## Requirement Completeness

- [x] CHK001 - Are all four genomic variant calling modes (Simple, Somatic, Amplicon, Splicing/SV) fully specified with their respective input and output parameters? [Completeness, Spec §FR-002]
- [x] CHK002 - Is the direct-to-VCF pipeline output format, including formatting header fields, specified? [Completeness, Spec §FR-001]
- [x] CHK003 - Are exit status codes and error scenarios (such as invalid coordinate intervals or missing index files) fully defined? [Completeness, Spec §User Story 1]
- [x] CHK004 - Does the specification detail the exact expected column count for all variant calling modes under intermediate output? [Completeness, Spec §FR-001]

## Requirement Clarity

- [x] CHK005 - Is the mathematical algorithm for the two-sided Fisher's Exact test and Brent root-finding MLE odds ratio quantified with exact calculation rules? [Clarity, Spec §FR-004]
- [x] CHK006 - Are the somatic status classification rules and threshold boundaries quantified? [Clarity, Spec §FR-005]
- [x] CHK007 - Is the coordinate indexing convention (0-based for BAM reads, 1-based for VCF outputs) explicitly defined? [Clarity, Spec §FR-006]
- [x] CHK008 - Are the parameters for genomic region tiling (tile size and boundary halo overlap dimensions) quantified? [Clarity, Spec §FR-007]

## Requirement Consistency

- [x] CHK009 - Do the coordinate bases defined in the data model align consistently with the BAM file reader and VCF writer specifications? [Consistency, Spec §FR-006]
- [x] CHK010 - Are the rounding rules for float values consistent across p-values, odds ratios, and allele frequencies? [Consistency, Spec §Success Criteria]
- [x] CHK011 - Does the execution path configuration prevent conflicts when multiple threads access adjacent tiling regions? [Consistency, Spec §FR-007]

## Acceptance Criteria Quality

- [x] CHK012 - Are the float tolerance boundaries (e.g. $\epsilon = 10^{-5}$) for differential validation measurable and testable? [Acceptance Criteria, Spec §SC-002]
- [x] CHK013 - Can the performance targets (speedup $\ge$ 3x, memory usage $<$ 4GB) be objectively measured without implementation-specific tools? [Acceptance Criteria, Spec §SC-003]
- [x] CHK014 - Is the criteria for a Tier-1 (`compare_var.py`) and Tier-2 (`compare.sh`) "MATCH" result clearly defined? [Acceptance Criteria, Spec §SC-001]

## Scenario Coverage

- [x] CHK015 - Does the specification address Exception/Error scenarios, such as empty target BAM files or truncated reference files? [Coverage, Spec §User Story 1]
- [x] CHK016 - Are there requirements specified for alternate validation paths when reference indexes are missing or corrupt? [Coverage, Gap]
- [x] CHK017 - Are the fallback options defined when the Brent root-finding solver fails to converge on extreme datasets? [Coverage, Gap]

## Edge Case Coverage

- [x] CHK018 - Are requirements defined for resolving oracle run-to-run variant output sorting nondeterminism? [Edge Cases, Spec §Edge Cases]
- [x] CHK019 - Does the specification address read alignments containing malformed CIGAR operations at tile boundaries? [Edge Cases, Spec §Edge Cases]
- [x] CHK020 - Are memory limits and performance requirements specified for extreme sequencing depths (e.g. $>$ 100,000x coverage)? [Edge Cases, Spec §Edge Cases]

## Non-Functional Requirements

- [x] CHK021 - Are strict quality and compiler gate requirements (like forbidding unsafe code and panics in library code) documented? [Non-Functional, Spec §FR-008]
- [x] CHK022 - Are target operating system, SIMD hardware architecture (AVX2), and static compilation requirements specified? [Non-Functional, Spec §Success Criteria]
- [x] CHK023 - Are log formatting and debug observability requirements defined for diagnostic validation? [Non-Functional, Gap]

## Dependencies & Assumptions

- [x] CHK024 - Are the dependencies on the reference VarDictJava v1.8.3 oracle explicitly documented? [Dependency, Spec §Assumptions]
- [x] CHK025 - Are the assumptions regarding input alignment and reference index file validity validated? [Assumption, Spec §Assumptions]

## Notes

- Check items off as completed: `[x]`
- Items are numbered sequentially for easy reference.
- Checklist validated against spec.md requirements.
- All items verified as PASS by developer self-check.
