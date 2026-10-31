### Supplementary Table S4 (Table S4). Verified Truth-Set Operating Points, Variant-Class Stratification & 18-Arm Somatic Noise-Floor Matrix

| Source Cell ID | Cohort / Platform | Stratum / Regime / Operating Point | Engine / Parity | Primary Metric 1 (TP / Candidates / SNV) | Primary Metric 2 (FP / PASS / InDel) | Summary Score / F1 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| s4_tier3_truth_and_qual_sweep | HG002 GRCh38 | All-PASS (Unfiltered) | VarDictRust == VarDictJava | TP=23,062 | FP=22,458, FN=813 | P=0.5038, R=0.9659, F1=0.6622 |
| s4_tier3_truth_and_qual_sweep | HG002 GRCh38 | Optimal QUAL (>=128.000) | VarDictRust == VarDictJava | TP=23,004 | FP=553, FN=871 | P=0.9763, R=0.9635, F1=0.9699 |
| s4_tier3_truth_and_qual_sweep | SEQC2 HCC1395 | All-PASS (Unfiltered) | VarDictRust == VarDictJava | TP=1,069 | FP=11,334, FN=140 | P=0.0852, R=0.8842, F1=0.1555 |
| s4_tier3_truth_and_qual_sweep | SEQC2 HCC1395 | Optimal QUAL (>=126.000) | VarDictRust == VarDictJava | TP=794 | FP=293, FN=415 | P=0.7279, R=0.6567, F1=0.6905 |
| s4_tier3_truth_and_qual_sweep | HG008 GIAB v0.3 | All-PASS (Unfiltered) | VarDictRust == VarDictJava | TP=619 | FP=1,565, FN=288 | P=0.2687, R=0.6825, F1=0.3856 |
| s4_tier3_truth_and_qual_sweep | HG008 GIAB v0.3 | Optimal QUAL (>=176.000) | VarDictRust == VarDictJava | TP=593 | FP=38, FN=314 | P=0.9354, R=0.6538, F1=0.7696 |
| evidence_certified | HG005 GRCh38 | All-PASS (Unfiltered) | VarDictRust == VarDictJava | TP=29 | FP=2, FN=1 | P=0.9355, R=0.9667, F1=0.9508 |
| evidence_certified | HG005 GRCh38 | Optimal QUAL (>=226.0) | VarDictRust == VarDictJava | TP=29 | FP=1, FN=1 | P=0.9667, R=0.9667, F1=0.9667 |
| evidence_certified | NA24385 GRCh37 | All-PASS (Unfiltered) | VarDictRust == VarDictJava | TP=5 | FP=81, FN=0 | P=0.0581, R=1.0000, F1=0.1099 |
| evidence_certified | NA24385 GRCh37 | Optimal QUAL (>=160.0) | VarDictRust == VarDictJava | TP=5 | FP=0, FN=0 | P=1.0000, R=1.0000, F1=1.0000 |
| s4_tier3_truth_and_qual_sweep | HG002 | Variant-Class & AF Breakdown | rust | SNV=209,030, Ins=2,305, Del=10,975 | MNP=2,065, SV=3,187 | PASS=48,447 / Total=227,562 |
| s4_tier3_truth_and_qual_sweep | HG002 | Variant-Class & AF Breakdown | java | SNV=209,030, Ins=2,305, Del=10,975 | MNP=2,065, SV=3,187 | PASS=48,447 / Total=227,562 |
| s4_tier3_truth_and_qual_sweep | HG002 | Variant-Class & AF Breakdown | cpp | SNV=209,066, Ins=2,309, Del=10,949 | MNP=2,049, SV=3,187 | PASS=48,431 / Total=227,560 |
| s4_tier3_truth_and_qual_sweep | SEQC2 | Variant-Class & AF Breakdown | rust | SNV=594,344, Ins=8,966, Del=15,716 | MNP=4,383, SV=35 | PASS=120,307 / Total=623,444 |
| s4_tier3_truth_and_qual_sweep | SEQC2 | Variant-Class & AF Breakdown | java | SNV=594,344, Ins=8,966, Del=15,716 | MNP=4,383, SV=35 | PASS=120,307 / Total=623,444 |
| s4_tier3_truth_and_qual_sweep | SEQC2 | Variant-Class & AF Breakdown | cpp | SNV=594,344, Ins=8,966, Del=15,715 | MNP=4,383, SV=34 | PASS=120,305 / Total=623,442 |
| s4_tier3_truth_and_qual_sweep | HG008 | Variant-Class & AF Breakdown | rust | SNV=15,213, Ins=1,161, Del=1,410 | MNP=614, SV=28 | PASS=11,358 / Total=18,426 |
| s4_tier3_truth_and_qual_sweep | HG008 | Variant-Class & AF Breakdown | java | SNV=15,213, Ins=1,161, Del=1,410 | MNP=614, SV=28 | PASS=11,358 / Total=18,426 |
| s4_tier3_truth_and_qual_sweep | HG008 | Variant-Class & AF Breakdown | cpp | SNV=15,214, Ins=1,161, Del=1,408 | MNP=614, SV=28 | PASS=11,356 / Total=18,425 |
| s5_cross_platform_short_read_chr20_panel | Element AVITI 111x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=61 | PASS Somatic=13 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Element AVITI 111x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=61 | PASS Somatic=13 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Element AVITI 111x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=186 | PASS Somatic=13 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | PacBio Onso 136x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=76 | PASS Somatic=19 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | PacBio Onso 136x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=78 | PASS Somatic=19 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | PacBio Onso 136x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=185 | PASS Somatic=19 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Illumina 195x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=115 | PASS Somatic=24 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Illumina 195x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=131 | PASS Somatic=24 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Illumina 195x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=1,199 | PASS Somatic=23 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima ppmSeq 157x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=574 | PASS Somatic=133 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima ppmSeq 157x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=581 | PASS Somatic=136 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima ppmSeq 157x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=2,791 | PASS Somatic=159 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Roche SBX 119x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=770 | PASS Somatic=293 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Roche SBX 119x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=788 | PASS Somatic=295 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Roche SBX 119x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=3,354 | PASS Somatic=296 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima UG100 306x | Min VAF 1%, >=2 Reads (-f 0.01 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=944 | PASS Somatic=376 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima UG100 306x | All VAFs, >=2 Reads (-f 0 -r 2) | 100% Exact Match (Java == Rust) | Total Candidates=1,519 | PASS Somatic=472 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s5_cross_platform_short_read_chr20_panel | Ultima UG100 306x | Singleton Noise Floor (-f 0 -r 1) | 100% Exact Match (Java == Rust) | Total Candidates=4,994 | PASS Somatic=538 | 0 diffs vs. VarDictJava (eps=1e-6) |
| s2_fig2_hg002_chr20_truth | HG002 Chr20 Exome (Daniels et al. 2024) | 7-Caller Truth Benchmark | VarDictRust == VarDictJava | Calls=7,422 | SNV F1=0.9858 | InDel F1=0.8956 (see Table 3) |
| s5_grch38_field_callers_truth | HG002, SEQC2 & HG008 GRCh38 | 17-Run Native GRCh38 Caller Matrix | VarDictRust == VarDictJava | HG002 SNV F1=0.6619 | HG002 InDel F1=0.7177 | 100% Exact Match (see Table 3) |
| s6_low_vaf_titration_truth_curve | HG008 Tumor-Purity Titration (20% to 1%) | 30-Run Low-VAF Somatic Curve (-f 0.002) | VarDictRust == VarDictJava | 20% Purity SNV F1=0.7512 | 1% Purity SNV F1=0.0569 | 100% Exact Match across all 5 purities (see Table 3) |
| s6_hg002_v50q_cmrg_mosaic_eval | HG002 v4.2.1, v5.0q T2T-Q100, CMRG & Mosaic | 12-Evaluation Next-Gen Germline & Mosaic Suite | VarDictRust == VarDictJava (<=50 kbp BED) | v5.0q SNV F1=0.6650 | CMRG SNV F1=0.6791 | Mosaic Recall=98.82% (see Table 7) |
| s6_hg008_modern_ml_and_commercial_callers | HG008 GIAB v0.4 Somatic Benchmark | 11-Caller Modern ML & Commercial Landscape | VarDictRust v1.0.6 Full WGS & Slice | Full WGS SNV F1=0.7886 | Full WGS InDel F1=0.4898 | 25 evaluations across 11 callers (see Table 7) |

*VarDictRust and VarDictJava produce identical rtg vcfeval accuracy at both All-PASS and Optimal QUAL operating points across all 5 truth cohorts and 18/18 exact-match noise-floor regimes.*
