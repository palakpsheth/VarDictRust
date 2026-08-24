# Data Model: VarDictRust

This document outlines the core domain entities, data representations, and state transitions for the Rust rewrite of VarDict.

## 1. Domain Entities

### AlignmentRead
Represents a mapped sequencing read retrieved from BAM/CRAM alignment files via the `noodles` library.

| Attribute | Type | Description |
| :--- | :--- | :--- |
| `name` | `Vec<u8>` | Query template name. |
| `flags` | `u16` | SAM flags mapping read pairing, mappings, and strand. |
| `reference_sequence_id` | `usize` | Target reference sequence ID mapping to the BAM header. |
| `alignment_start` | `Pos0` | 0-based coordinate of the leftmost alignment start. |
| `mapping_quality` | `u8` | MAPQ score. |
| `cigar` | `Cigar` | Read mapping CIGAR representation. |
| `sequence` | `Vec<u8>` | Query sequence bases. |
| `quality_scores` | `Vec<u8>` | Quality scores in Phred+33 format. |

### RegionTile
Represents a partitioned slice of a chromosome mapped to a thread-level unit of work.

| Attribute | Type | Description |
| :--- | :--- | :--- |
| `chromosome` | `Vec<u8>` | Target chromosome/contig name. |
| `start` | `Pos1` | 1-based start coordinate of the active region. |
| `end` | `Pos1` | 1-based end coordinate of the active region. |
| `halo_start` | `Pos0` | 0-based start including boundary padding for read alignment overlap. |
| `halo_end` | `Pos0` | 0-based end including boundary padding for read alignment overlap. |

### VariantCall
Represents a candidate somatic or germline variant called by the engine before filtering or status classification.

| Attribute | Type | Description |
| :--- | :--- | :--- |
| `chromosome` | `Vec<u8>` | Contig name. |
| `position` | `Pos1` | 1-based coordinate of the variant. |
| `ref_allele` | `Vec<u8>` | Reference allele sequence. |
| `alt_allele` | `Vec<u8>` | Alternate allele sequence. |
| `depth` | `u32` | Total read depth at the site. |
| `vd` | `u32` | Variant-supporting read depth. |
| `allele_frequency` | `f64` | Frequency of the alternate allele (`vd / depth`). |
| `strand_bias` | `StrandBias` | Forward/reverse strand coverage metrics. |
| `qual_mean` | `f64` | Mean base quality of supporting reads. |
| `fisher_p_value` | `Option<f64>` | Fisher's Exact test p-value if requested (`--fisher`). |
| `odds_ratio` | `Option<f64>` | Odds ratio MLE calculation if requested (`--fisher`). |
| `somatic_status` | `SomaticStatus` | somatic categorization determined by postprocessing filters. |

## 2. Validation and Coordinate Safety

### coordinate Indexing
- BAM/CRAM file mappings are retrieved as **0-based half-open** coordinates, represented by the `Pos0` type.
- Internal pileup calculations and VCF outputs are calculated in **1-based** coordinates, represented by the `Pos1` type.
- Conversions MUST only occur through explicit methods:
  - `Pos0::to_one_based()` $\rightarrow$ `Pos1`
  - `Pos1::to_zero_based()` $\rightarrow$ `Pos0`
- Raw integer additions to coordinates are prohibited to prevent coordinate basis conversion bugs.

### Floating-Point Rounding
- All floating-point fields written to outputs must implement R's **HALF_EVEN** (Banker's) rounding to exactly 5 decimal places.
- Iterative math routines for Brent's Root solver MLE calculation must converge with a tolerance of `ulp(1)^0.25`.
