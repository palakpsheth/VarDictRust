### Table 1. Multi-Modal Algorithmic Equivalence & Intermediate (.var) / VCF Parity Summary (78/87 Sweep Sub-Arms Exact Match)

| Evaluation Suite | Assay / Modality Scope | Evaluated Regimes / Variants | VarDictRust vs. VarDictJava Verdict |
| :--- | :--- | :--- | :--- |
| Amplicon Baseline | EGFR/MRE11a amplicon (-a) calling & coordinate parity | MATCH (2 variants) | 100% Exact Match |
| Somatic Benchmark Slice | NIST GIAB HG008 tumor-normal somatic slice parity | MATCH (18426 variants) | 100% Exact Match |
| Multi-Region Reference Offset | HG002 Chr20 multi-region reference coordinate offset verification | MATCH (2572 variants) | 100% Exact Match |
| Tile-Size Invariance | Halo-tiled equivalence (0.5, 1, 2, 10 Mbp tiles) on SEQC2 HCC1395 | 4 tile sizes x 5,090 variants | 100% Exact Match across all tile sizes |
| 13-Mode Functional Suite | 13-mode 3-way functional CLI equivalence sweep | 13 CLI modes | 100% Exact Match (13/13 modes) |
| Transcriptome Splice-Aware RNA-Seq | 51.8 GiB HG002 GM24385 mRNA (-i splice-aware mode) | 12,668 splice rows | 100% Exact Match (diff=0 B) |
| All-Positions Pileup (-f 0) | Chromosome 1 10 Mbp all-positions (-f 0) pileup equivalence | MATCH (7032 variants) | 100% Exact Match |
| Lossless BAM-vs-CRAM Twin | HG008 tumor-normal BAM-vs-CRAM bitwise equivalence | Bitwise VCF & .var comparison | 100% Exact Match (0 byte differences) |
| Extended & Multimodal Suites | 17-modality extended coverage & 18-slice multimodal regression suites | 17/17 modes; 18/18 slices (1,129 vars) | 100% Exact Match |
| Full-Genome Germline WGS | Full-genome HG002 WGS on <=50 kbp GIAB callable BED intervals | 8,111,171 callable variants | 100% Exact Match on callable intervals |
| Multi-Regime Sweep: 4-Assay Regression | Targeted amplicon, HG002, HG008 & SEQC2 slices (.var & .vcf) | 4 sub-arms | 4/4 Sub-Arms Exact Match |
| Multi-Regime Sweep: Contiguous Chr1 Spans | Chromosome 1 contiguous spans (1 Mbp to 249 Mbp; .var & .vcf) | 1 sub-arms | 0/1 Sub-Arms Match (1 Unchunked Multi-Tile Boundary Diffs) |
| Multi-Regime Sweep: HG002 Chr20 Exome | HG002 Chromosome 20 whole-exome slice (.var & .vcf) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: Multimodal & CRAM Twin | 18-slice multimodal suite & CRAM-vs-BAM twin (.var & .vcf) | 4 sub-arms | 4/4 Sub-Arms Exact Match |
| Multi-Regime Sweep: 51.8 GiB RNA-Seq | HG002 GM24385 splice-aware RNA-seq (-i; .var & .vcf) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: Chr1 All-Positions | Chromosome 1 10 Mbp all-positions (-f 0; .var & .vcf) | 2 sub-arms | 0/2 Sub-Arms Match (2 Unchunked Multi-Tile Boundary Diffs) |
| Multi-Regime Sweep: 13-Mode CLI Matrix | 13 functional CLI modes across germline, somatic & amplicon regimes | 12 sub-arms | 12/12 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Panel SRR15006540 | Clinical AML germline panel SRR15006540 (.var & .vcf) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Panel SRR8657348 | Clinical AML germline panel SRR8657348 (.var & .vcf) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Somatic Pair | Clinical AML tumor-normal somatic panel (SRR15006376 / SRR15006386) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: Ultra-Deep MRD Stress | Singleton & all-positions MRD stress (-f 0 / -f 0.001 / -r 1) | 30 sub-arms | 28/30 Sub-Arms Match (2 Unchunked Multi-Tile Boundary Diffs) |
| Multi-Regime Sweep: WGS Low-VAF MRD | Low-VAF somatic MRD matrix (-f 0.002 / -f 0.0005 / -r 1) | 24 sub-arms | 20/24 Sub-Arms Match (4 Unchunked Multi-Tile Boundary Diffs) |

*All parity comparisons align records by CHROM/POS/REF/ALT with float tolerance eps=1e-6 after R banker's rounding.*
