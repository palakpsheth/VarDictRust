### Table 7. Next-Generation Truth Benchmarks (GIAB v5.0q T2T-Q100, CMRG v1.00, Mosaic SNV v1.10) & HG008 v0.4 11-Caller Landscape

| Source Cell ID | Benchmark / Scope | Caller / Configuration | Sequencing Chemistry | SNV Precision / Recall / F1 | InDel Precision / Recall / F1 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v4.2.1_All | VarDictRust_FullContigs | Illumina 2x250 bp WGS | 0.5081 / 0.9577 / 0.6639 | 0.5698 / 0.8836 / 0.6928 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v4.2.1_All | VarDictRust_GIAB50kb | Illumina 2x250 bp WGS | 0.5081 / 0.9576 / 0.6639 | 0.5697 / 0.8835 / 0.6928 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v4.2.1_All | VarDictJava_GIAB50kb | Illumina 2x250 bp WGS | 0.5081 / 0.9576 / 0.6639 | 0.5697 / 0.8835 / 0.6928 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v5.0q_T2T_Q100 | VarDictRust_FullContigs | Illumina 2x250 bp WGS | 0.5155 / 0.9367 / 0.6650 | 0.6151 / 0.7007 / 0.6551 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v5.0q_T2T_Q100 | VarDictRust_GIAB50kb | Illumina 2x250 bp WGS | 0.5083 / 0.8866 / 0.6461 | 0.5708 / 0.5281 / 0.5487 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_v5.0q_T2T_Q100 | VarDictJava_GIAB50kb | Illumina 2x250 bp WGS | 0.5083 / 0.8866 / 0.6461 | 0.5708 / 0.5281 / 0.5487 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_CMRG_v1.00_273Genes | VarDictRust_FullContigs | Illumina 2x250 bp WGS | 0.5482 / 0.8920 / 0.6791 | 0.6310 / 0.7702 / 0.6937 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_CMRG_v1.00_273Genes | VarDictRust_GIAB50kb | Illumina 2x250 bp WGS | 0.5355 / 0.7173 / 0.6132 | 0.5803 / 0.4904 / 0.5316 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_CMRG_v1.00_273Genes | VarDictJava_GIAB50kb | Illumina 2x250 bp WGS | 0.5355 / 0.7173 / 0.6132 | 0.5803 / 0.4904 / 0.5316 |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_Mosaic_SNV_v1.10 | VarDictRust_FullContigs | Illumina 2x250 bp WGS | 0.0000 / 0.9882 / 0.0000 | SNV-Only Truth Set (84/85 Mosaic SNVs) |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_Mosaic_SNV_v1.10 | VarDictRust_GIAB50kb | Illumina 2x250 bp WGS | 0.0000 / 0.9882 / 0.0000 | SNV-Only Truth Set (84/85 Mosaic SNVs) |
| s6_hg002_v50q_cmrg_mosaic_eval | NIST_Mosaic_SNV_v1.10 | VarDictJava_GIAB50kb | Illumina 2x250 bp WGS | 0.0000 / 0.9882 / 0.0000 | SNV-Only Truth Set (84/85 Mosaic SNVs) |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | VarDictRust_v1.0.6 | Illumina_195x_169x | 0.8262 / 0.7542 / 0.7886 | 0.6288 / 0.4011 / 0.4898 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | VarDictRust_v1.0.6 | Illumina_195x_169x | 0.8258 / 0.7810 / 0.8028 | 0.6703 / 0.4832 / 0.5616 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | VarDictRust_v1.0.6 | Illumina_195x_169x | 0.9916 / 0.6048 / 0.7513 | 0.9479 / 0.3095 / 0.4667 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | DRAGEN_v4.2.4_Linear | Illumina_WGS | 0.9546 / 0.8643 / 0.9072 | 0.9140 / 0.4320 / 0.5867 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | DRAGEN_v4.2.4_Linear | Illumina_WGS | 0.9570 / 0.8848 / 0.9195 | 0.9307 / 0.7067 / 0.8034 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | DRAGEN_v4.2.4_Linear | Illumina_WGS | 0.9905 / 0.7104 / 0.8274 | 0.8585 / 0.3095 / 0.4550 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | DRAGEN_v4.2.4_Graph | Illumina_WGS | 0.9591 / 0.8633 / 0.9087 | 0.9170 / 0.4323 / 0.5876 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | DRAGEN_v4.2.4_Graph | Illumina_WGS | 0.9615 / 0.8836 / 0.9209 | 0.9335 / 0.7067 / 0.8044 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | DRAGEN_v4.2.4_Graph | Illumina_WGS | 0.9905 / 0.7104 / 0.8274 | 0.8585 / 0.3095 / 0.4550 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | Lancet2_NYGC | Illumina_WGS | 0.9933 / 0.7936 / 0.8823 | 0.9833 / 0.2179 / 0.3568 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | Lancet2_NYGC | Illumina_WGS | 0.9960 / 0.8231 / 0.9013 | 0.9867 / 0.3618 / 0.5295 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | Lancet2_NYGC | Illumina_WGS | 1.0000 / 0.5843 / 0.7376 | 1.0000 / 0.1599 / 0.2757 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | Hartwig_SAGE_2026 | Illumina_WGS | 0.9956 / 0.8161 / 0.8969 | 0.9710 / 0.6309 / 0.7649 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | Hartwig_SAGE_2026 | Illumina_WGS | 0.9958 / 0.8413 / 0.9121 | 0.9653 / 0.7883 / 0.8679 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | Hartwig_SAGE_2026 | Illumina_WGS | 1.0000 / 0.6048 / 0.7537 | 0.9858 / 0.4728 / 0.6391 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | DeepSomatic_v1.6.0 | PacBio_HiFi | 0.9886 / 0.8512 / 0.9148 | 0.9158 / 0.1845 / 0.3072 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | DeepSomatic_v1.6.0 | PacBio_HiFi | 0.9953 / 0.8763 / 0.9320 | 0.9140 / 0.3165 / 0.4702 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | DeepSomatic_v1.6.0 | PacBio_HiFi | 1.0000 / 0.6627 / 0.7971 | 1.0000 / 0.1395 / 0.2448 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | ClairS_v0.4.1 | PacBio_HiFi | 0.9793 / 0.8419 / 0.9054 | 0.9791 / 0.1302 / 0.2298 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | ClairS_v0.4.1 | PacBio_HiFi | 0.9793 / 0.8736 / 0.9234 | 0.9797 / 0.2264 / 0.3678 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | ClairS_v0.4.1 | PacBio_HiFi | 1.0000 / 0.6491 / 0.7872 | 1.0000 / 0.1122 / 0.2018 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_All | Ultima_DeepVariant_Somatic | Ultima_UG100 | 0.9749 / 0.9101 / 0.9414 | 0.9429 / 0.2459 / 0.3900 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_FullWGS_v0.4_NoGI | Ultima_DeepVariant_Somatic | Ultima_UG100 | 0.9751 / 0.9444 / 0.9595 | 0.9486 / 0.4218 / 0.5839 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | Ultima_DeepVariant_Somatic | Ultima_UG100 | 1.0000 / 0.7700 / 0.8701 | 0.9859 / 0.2381 / 0.3836 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | Mutect2 | Illumina_195x_169x | 0.9874 / 0.6661 / 0.7955 | 0.9431 / 0.3946 / 0.5564 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | Strelka2 | Illumina_195x_169x | 0.9896 / 0.8075 / 0.8893 | 0.9545 / 0.2857 / 0.4398 |
| s6_hg008_modern_ml_and_commercial_callers | HG008_Slice_v0.4_All | LoFreq | Illumina_195x_169x | 1.0000 / 0.7053 / 0.8272 | 0.0000 / 0.0000 / 0.0000 (0 Indel calls) |

*External HG008 v0.4 full-WGS callsets are sourced from NIST GIAB Liss_lab (DRAGEN, DeepSomatic, ClairS, Ultima DeepVariant, Hartwig SAGE, Lancet2); Mosaic SNV v1.10 reports unmasked recall (84/85 = 98.82%) as primary metric.*
