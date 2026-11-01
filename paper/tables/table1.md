### Table 1. Multi-Modal Algorithmic Equivalence & Intermediate (.var) / VCF Parity Summary (87/87 Sweep Sub-Arms Exact Match)

| Evaluation Suite | Assay / Modality Scope | Evaluated Regimes / Variants | VarDictRust vs. VarDictJava [1] Verdict<sup>a</sup> |
| :--- | :--- | :--- | :--- |
| Amplicon Baseline | EGFR/MRE11A amplicon (-a) calling & coordinate parity | 2 variants | 100% Exact Match |
| Somatic Benchmark Slice | NIST GIAB HG008 [19] tumor-normal somatic slice parity | 18,426 variants | 100% Exact Match |
| Multi-Region Reference Offset | GIAB HG002 [13, 18] Chr20 multi-region reference coordinate offset verification | 2,572 variants | 100% Exact Match |
| Tile-Size Invariance<sup>b</sup> | Halo-tiled equivalence (0.5, 1, 2, 10 Mbp tiles) on SEQC2 HCC1395 [10] | 4 tile sizes x 5,090 variants | 100% Exact Match across all tile sizes<sup>b</sup> |
| 13-Mode Functional Suite | 13-mode 3-way functional CLI equivalence sweep [1, 11] | 13 CLI modes | 100% Exact Match (13/13 modes) |
| Transcriptome Splice-Aware RNA-Seq | 51.8 GiB GIAB HG002 GM24385 [13] mRNA (-i splice-aware mode) | 12,668 splice rows | 100% Exact Match (0 byte differences) |
| All-Positions Pileup (-f 0, Single-Window)<sup>c</sup> | Chromosome 1 10 Mbp all-positions pileup (-f 0; .var & VCF)<sup>c</sup> | 7,032 VCF variants (7,077 .var rows) | 100% Exact Match<sup>c</sup> |
| Lossless BAM-vs-CRAM Twin | GIAB HG008 [19] tumor-normal BAM-vs-CRAM [26] bitwise equivalence (-G) | Bitwise VCF & .var comparison | 100% Exact Match (0 byte differences) |
| Extended & Multimodal Suites | 17-modality extended coverage & 18-slice multimodal regression suites | 17/17 modes; 18/18 slices (1,129 variants) | 100% Exact Match |
| Full-Genome Germline WGS | Full-genome GIAB HG002 [13, 18] WGS on <=50 kbp callable BED intervals | 8,111,171 callable variants | 100% Exact Match on callable intervals |
| Multi-Regime Sweep: 4-Assay Regression | Targeted amplicon, HG002 [13, 18], HG008 [19] & SEQC2 [10] slices (.var & VCF) | 4 sub-arms | 4/4 Sub-Arms Exact Match |
| Multi-Regime Sweep: Contiguous 249 Mbp -R chr1 | Whole-chromosome -R chr1:1-249250621 [11] (.var)<sup>c</sup> | 1 sub-arm | 1/1 Sub-Arm Exact Match<sup>c</sup> |
| Multi-Regime Sweep: HG002 Chr20 Exome | GIAB HG002 Chromosome 20 whole-exome slice [11] (.var & VCF) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: Multimodal & CRAM Twin | 18-slice multimodal suite & CRAM-vs-BAM twin [26] (.var & VCF) | 4 sub-arms | 4/4 Sub-Arms Exact Match |
| Multi-Regime Sweep: 51.8 GiB RNA-Seq | GIAB HG002 GM24385 [13] splice-aware RNA-seq (-i; .splice.tsv & .var) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: 10 Mbp Chr1 (-f 0) | Chromosome 1 10 Mbp single BED interval (.var & VCF)<sup>c</sup> | 2 sub-arms | 2/2 Sub-Arms Exact Match<sup>c</sup> |
| Multi-Regime Sweep: 12-Sub-Arm CLI Matrix | 12 CLI mode sub-arms across germline, somatic, amplicon & CRAM regimes (.var) | 12 sub-arms | 12/12 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Panel SRR15006540 | Clinical AML germline panel SRR15006540 [11] (.var & VCF) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Panel SRR8657348 | Clinical AML germline panel SRR8657348 [11] (.var & VCF) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: AML Somatic Pair | Clinical AML tumor-normal somatic panel (SRR15006376 / SRR15006386) [11] (.var & VCF) | 2 sub-arms | 2/2 Sub-Arms Exact Match |
| Multi-Regime Sweep: Ultra-Deep MRD Stress | Ultra-deep MRD stress (-f 0.0005 / -f 0 + --fisher, up to 45kx depth)<sup>d</sup> | 30 sub-arms | 30/30 Sub-Arms Exact Match<sup>d</sup> |
| Multi-Regime Sweep: WGS Low-VAF MRD | Low-VAF WGS & WES MRD matrix (-f 0.001 / -f 0 + --fisher)<sup>d</sup> | 24 sub-arms | 24/24 Sub-Arms Exact Match<sup>d</sup> |

*<sup>a</sup>All parity comparisons align records by CHROM/POS/REF/ALT against the VarDictJava v1.8.3 [1] oracle with floating-point tolerance eps=1e-6 after replicating R 5-decimal round-half-to-even (HALF_EVEN) rounding; across all 87/87 (100.00%) comparative sub-arms, the exact two-sided Clopper-Pearson 95% binomial confidence interval is 95.85%–100.00%. <sup>b</sup>Halo-tiled streaming partitions intervals >1 Mbp into 100 kbp tiles with 500 bp overlap halos and deterministic boundary deduplication, yielding identical calls across 0.5, 1.0, 2.0, and 10.0 Mbp tile sizes. <sup>c</sup>Standard BED intervals (<=1 Mbp) execute untiled by default; for the two multi-megabase single-interval stress suites (10 Mbp chr1 at -f 0 and 249.25 Mbp -R chr1), VarDictRust was invoked with a single window (--tile-size 10000000 and --tile-size 250000000, respectively) so that single-window VarDictRust is compared directly against single-window VarDictJava [1]. <sup>d</sup>Evaluated both without --fisher (replicating R teststrandbias.R / testsomatic.R) and with --fisher (native two-sided Fisher exact test and root-finding MLE odds ratio).*
