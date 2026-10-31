### Table 3. Germline & Somatic Truth Accuracy (H3), Native GRCh38 Field Callers & Low-VAF Tumor-Purity Titration (20% to 1% Purity)

| Source Cell ID | Cohort / Purity Stratum | Caller / Engine | Wall Time / Peak RSS | SNV Precision / Recall / F1 | InDel Precision / Recall / F1 | Concordance / Parity Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | VarDictRust | 7,422 calls | 0.9981 / 0.9738 / 0.9858 | 0.9747 / 0.8283 / 0.8956 | VarDictRust == VarDictJava EQUAL |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | VarDictJava | 7,422 calls | 0.9981 / 0.9738 / 0.9858 | 0.9747 / 0.8283 / 0.8956 | VarDictRust == VarDictJava EQUAL |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | vardictcpp | 7,421 calls | 0.9981 / 0.9738 / 0.9858 | 0.9747 / 0.8283 / 0.8956 | Field Comparator |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | octopus | 2,436 calls | 0.9973 / 0.9878 / 0.9925 | 0.9336 / 0.9056 / 0.9194 | Field Comparator |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | FreeBayes | 4,292 calls | 0.9972 / 0.9765 / 0.9867 | 0.9439 / 0.8670 / 0.9038 | Field Comparator |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | LoFreq | 2,598 calls | 0.9908 / 0.9701 / 0.9803 | 0.9358 / 0.8755 / 0.9047 | Field Comparator |
| s2_fig2_hg002_chr20_truth | HG002 chr20 Exome (Daniels 2024) | Mutect2 | 2,246 calls | 0.9995 / 0.9330 / 0.9651 | 0.9939 / 0.7039 / 0.8241 | Field Comparator |
| s4_tier3_truth_and_qual_sweep | HG002 Slice (Optimal QUAL>=128.000) | VarDictRust == VarDictJava | TP=23,004, FP=553 | 0.9763 / 0.9635 / 0.9699 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| s4_tier3_truth_and_qual_sweep | HG002 Slice (All-PASS -f 0.01) | VarDictRust == VarDictJava | TP=23,062, FP=22,458 | 0.5038 / 0.9659 / 0.6622 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| s4_tier3_truth_and_qual_sweep | SEQC2 HCC1395 (Optimal QUAL>=126.000) | VarDictRust == VarDictJava | TP=794, FP=293 | 0.7279 / 0.6567 / 0.6905 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| s4_tier3_truth_and_qual_sweep | SEQC2 HCC1395 (All-PASS -f 0.01) | VarDictRust == VarDictJava | TP=1,069, FP=11,334 | 0.0852 / 0.8842 / 0.1555 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| s4_tier3_truth_and_qual_sweep | HG008 Slice v0.3 (Optimal QUAL>=176.000) | VarDictRust == VarDictJava | TP=593, FP=38 | 0.9354 / 0.6538 / 0.7696 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| s4_tier3_truth_and_qual_sweep | HG008 Slice v0.3 (All-PASS -f 0.01) | VarDictRust == VarDictJava | TP=619, FP=1,565 | 0.2687 / 0.6825 / 0.3856 | Combined All-Variant vcfeval | 0 site-diffs (Java == Rust) |
| evidence_certified | HG005 GIAB v4.2.1 (Optimal / All-PASS) | VarDictRust == VarDictJava | TP=29 | Opt: 0.9667/0.9667/0.9667 | All-PASS: 0.9355/0.9667/0.9508 | EQUAL (0 site-diffs) |
| evidence_certified | NA24385 GRCh37 v4.2.1 (Optimal / All-PASS) | VarDictRust == VarDictJava | TP=5 | Opt: 1.0000/1.0000/1.0000 | All-PASS: 0.0581/1.0000/0.1099 | EQUAL (0 site-diffs) |
| s5_grch38_field_callers_truth | HG002_GRCh38 | VarDictRust | 0:21.30 (21.30s) / 474.8 MB | 0.5041 / 0.9634 / 0.6619 | 0.5781 / 0.9463 / 0.7177 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | HG002_GRCh38 | VarDictJava | 2:39.59 (159.59s) / 10732.2 MB | 0.5041 / 0.9634 / 0.6619 | 0.5781 / 0.9463 / 0.7177 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | HG002_GRCh38 | vardictcpp | 0:21.14 (21.14s) / 859.1 MB | 0.5043 / 0.9634 / 0.6620 | 0.5783 / 0.9463 / 0.7179 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG002_GRCh38 | FreeBayes | 6:08.72 (368.72s) / 206.7 MB | 0.9686 / 0.9598 / 0.9642 | 0.9373 / 0.9538 / 0.9455 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG002_GRCh38 | LoFreq | 6:16.78 (376.78s) / 119.0 MB | 0.8924 / 0.9909 / 0.9391 | 0.5809 / 0.9601 / 0.7238 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG002_GRCh38 | Strelka2 | 8:22.19 (502.19s) / 130.1 MB | 0.9931 / 0.9909 / 0.9920 | 0.9438 / 0.9836 / 0.9633 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG002_GRCh38 | octopus | 3:11.50 (191.50s) / 3558.9 MB | 0.9911 / 0.9861 / 0.9886 | 0.9051 / 0.9810 / 0.9415 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG008_v0.3 | VarDictRust | 0:17.14 (17.14s) / 510.5 MB | 0.9757 / 0.6552 / 0.7840 | 0.9370 / 0.4034 / 0.5640 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | HG008_v0.3 | VarDictJava | 1:19.15 (79.15s) / 17448.2 MB | 0.9757 / 0.6552 / 0.7840 | 0.9370 / 0.4034 / 0.5640 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | HG008_v0.3 | vardictcpp | 0:17.30 (17.30s) / 565.1 MB | 0.9757 / 0.6552 / 0.7840 | 0.9444 / 0.4034 / 0.5653 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG008_v0.3 | Mutect2 | 7:27.31 (447.31s) / 1263.1 MB | 0.9825 / 0.6438 / 0.7779 | 0.9355 / 0.3932 / 0.5537 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG008_v0.3 | Strelka2 | 0:36.35 (36.35s) / 119.0 MB | 0.9896 / 0.7794 / 0.8720 | 0.9545 / 0.2847 / 0.4386 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | HG008_v0.3 | LoFreq | 1:18.71 (78.71s) / 251.5 MB | 1.0000 / 0.6814 / 0.8105 | 0.0000 / 0.0000 / 0.0000 (0 Indel calls) | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | VarDictRust | 4:39.29 (279.29s) / 514.3 MB | 0.7336 / 0.7532 / 0.7433 | 0.3077 / 0.4800 / 0.3750 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | VarDictJava | 1:06:07 (3967.00s) / 14275.5 MB | 0.7336 / 0.7532 / 0.7433 | 0.3077 / 0.4800 / 0.3750 | Java == Rust EQUAL |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | vardictcpp | 4:53.71 (293.71s) / 2259.8 MB | 0.7336 / 0.7532 / 0.7433 | 0.3077 / 0.4800 / 0.3750 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | Mutect2 | 44:01.05 (2641.05s) / 1580.9 MB | 0.6152 / 0.8361 / 0.7089 | 0.6087 / 0.8400 / 0.7059 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | Strelka2 | 2:52.88 (172.88s) / 118.8 MB | 0.3044 / 0.8689 / 0.4509 | 0.4324 / 0.6400 / 0.5161 | GRCh38 Field Caller |
| s5_grch38_field_callers_truth | SEQC2_v1.2.1 | LoFreq | 9:01.90 (541.90s) / 281.6 MB | 0.9526 / 0.7627 / 0.8471 | 0.0000 / 0.0000 / 0.0000 (0 Indel calls) | GRCh38 Field Caller |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | VarDictRust | 0:08.72 (8.72s) / 762.3 MB | 0.8790 / 0.6559 / 0.7512 | 0.6162 / 0.2075 / 0.3104 | MATCH |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | VarDictJava | 1:06.76 (66.76s) / 12027.7 MB | 0.8790 / 0.6559 / 0.7512 | 0.6162 / 0.2075 / 0.3104 | ORACLE |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | vardictcpp | 0:10.94 (10.94s) / 1048.3 MB | 0.8770 / 0.5707 / 0.6914 | 0.5823 / 0.1565 / 0.2466 | MISMATCH |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | Strelka2 | 0:36.62 (36.62s) / 118.8 MB | 1.0000 / 0.7070 / 0.8283 | 1.0000 / 0.1156 / 0.2073 | NA |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | LoFreq | 1:18.14 (78.14s) / 252.2 MB | 1.0000 / 0.4072 / 0.5787 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 20% Purity (10.0% VAF) | Mutect2 | 7:10.48 (430.48s) / 1130.0 MB | 1.0000 / 0.5826 / 0.7363 | 0.9032 / 0.0952 / 0.1723 | NA |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | VarDictRust | 0:08.76 (8.76s) / 839.4 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 | MATCH |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | VarDictJava | 1:06.58 (66.58s) / 13019.7 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 | ORACLE |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | vardictcpp | 0:10.91 (10.91s) / 1020.3 MB | 0.6726 / 0.5179 / 0.5852 | 0.4082 / 0.1361 / 0.2041 | MISMATCH |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | Strelka2 | 0:36.61 (36.61s) / 118.9 MB | 1.0000 / 0.4719 / 0.6412 | 1.0000 / 0.0476 / 0.0909 | NA |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | LoFreq | 1:14.64 (74.64s) / 251.8 MB | 1.0000 / 0.0716 / 0.1335 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 10% Purity (5.0% VAF) | Mutect2 | 7:16.35 (436.35s) / 1171.3 MB | 1.0000 / 0.3237 / 0.4891 | 0.9375 / 0.0510 / 0.0968 | NA |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | VarDictRust | 0:08.65 (8.65s) / 822.7 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 | MATCH |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | VarDictJava | 1:06.17 (66.17s) / 12117.8 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 | ORACLE |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | vardictcpp | 0:10.90 (10.90s) / 1077.7 MB | 0.4932 / 0.3680 / 0.4215 | 0.2183 / 0.1054 / 0.1422 | MISMATCH |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | Strelka2 | 0:36.56 (36.56s) / 118.8 MB | 1.0000 / 0.1891 / 0.3181 | 1.0000 / 0.0272 / 0.0530 | NA |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | LoFreq | 1:14.67 (74.67s) / 251.5 MB | 1.0000 / 0.0051 / 0.0102 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 5% Purity (2.5% VAF) | Mutect2 | 7:14.72 (434.72s) / 1131.1 MB | 1.0000 / 0.0988 / 0.1798 | 0.8889 / 0.0272 / 0.0528 | NA |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | VarDictRust | 0:08.60 (8.60s) / 813.6 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 | MATCH |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | VarDictJava | 1:05.16 (65.16s) / 13443.3 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 | ORACLE |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | vardictcpp | 0:11.00 (11.00s) / 1088.3 MB | 0.2430 / 0.1175 / 0.1584 | 0.1391 / 0.0544 / 0.0782 | MISMATCH |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | Strelka2 | 0:36.70 (36.70s) / 118.8 MB | 1.0000 / 0.0187 / 0.0368 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | LoFreq | 1:14.80 (74.80s) / 251.5 MB | 0.0000 / 0.0000 / 0.0000 (0 calls) | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 2% Purity (1.0% VAF) | Mutect2 | 6:59.88 (419.88s) / 1178.6 MB | 1.0000 / 0.0119 / 0.0236 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | VarDictRust | 0:08.53 (8.53s) / 803.0 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 | MATCH |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | VarDictJava | 1:05.41 (65.41s) / 12946.3 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 | ORACLE |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | vardictcpp | 0:10.95 (10.95s) / 1020.9 MB | 0.0856 / 0.0426 / 0.0569 | 0.0660 / 0.0238 / 0.0350 | MISMATCH |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | Strelka2 | 0:36.58 (36.58s) / 119.0 MB | 1.0000 / 0.0017 / 0.0034 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | LoFreq | 1:14.51 (74.51s) / 252.0 MB | 0.0000 / 0.0000 / 0.0000 (0 calls) | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |
| s6_low_vaf_titration_truth_curve | HG008 1% Purity (0.5% VAF) | Mutect2 | 6:47.62 (407.62s) / 1180.5 MB | 1.0000 / 0.0017 / 0.0034 | 0.0000 / 0.0000 / 0.0000 (0 calls) | NA |

*Undefined 0/0 precision/F1 at 0.0000 recall (e.g., LoFreq indels or 0.5-1.0% VAF) are normalized to 0.0000 with explicit 0-call annotations. Low-VAF titration uses in silico HG008 T/N read admixtures called with -f 0.002.*
