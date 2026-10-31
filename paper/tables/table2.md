### Table 2. Direct Replication of vardictcpp (Wolff 2026) Across Tier P AML Panels, Field Callers & Table 3 Whole-Chromosome -R chr1

| Cell ID | Sample / Cohort | Threads | BED Intervals | VarDictJava Wall (s) / Peak RSS (MB) | vardictcpp Wall (s) / Peak RSS (MB) | VarDictRust Wall (s) / Peak RSS (MB) | Rust vs Java Speedup / RSS Ratio | Parity Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| s2_tierp_v2_SRR15006375_th1 | SRR15006375 | 1 | 12,292 | 408.19 (400.74-418.94)s / 2537.2 MB | 109.66 (109.19-109.90)s / 63.8 MB | 121.93 (121.14-122.38)s / 210.3 MB | 3.35x / 12.1x leaner | MATCH (4056 rows); MATCH (3966 variants) |
| s2_tierp_v2_SRR15006375_th8 | SRR15006375 | 8 | 12,292 | 232.70 (229.14-235.27)s / 4410.5 MB | 15.35 (15.30-15.42)s / 328.3 MB | 49.47 (49.39-49.59)s / 270.1 MB | 4.70x / 16.3x leaner | MATCH (4056 rows); MATCH (3966 variants) |
| s2_tierp_v2_SRR15006376_th1 | SRR15006376 | 1 | 105,949 | 334.16 (326.10-340.26)s / 4244.4 MB | 83.71 (83.03-84.32)s / 74.0 MB | 80.41 (80.08-80.84)s / 225.7 MB | 4.16x / 18.8x leaner | MATCH (12739 rows); MATCH (12655 variants) |
| s2_tierp_v2_SRR15006376_th8 | SRR15006376 | 8 | 105,949 | 164.36 (150.76-174.29)s / 7442.2 MB | 12.85 (12.50-13.09)s / 354.0 MB | 48.74 (48.65-48.86)s / 281.0 MB | 3.37x / 26.5x leaner | MATCH (12739 rows); MATCH (12655 variants) |
| s2_tierp_v2_SRR15006386_th1 | SRR15006386 | 1 | 11,732 | 334.94 (331.87-338.80)s / 2423.4 MB | 98.20 (98.09-98.37)s / 56.2 MB | 113.52 (113.17-114.04)s / 187.5 MB | 2.95x / 12.9x leaner | MATCH (4014 rows); MATCH (3925 variants) |
| s2_tierp_v2_SRR15006386_th8 | SRR15006386 | 8 | 11,732 | 180.62 (178.99-182.17)s / 4434.7 MB | 15.55 (14.95-16.69)s / 298.6 MB | 43.78 (43.64-43.88)s / 305.4 MB | 4.13x / 14.5x leaner | MATCH (4014 rows); MATCH (3925 variants) |
| s2_tierp_v2_SRR15006540_th1 | SRR15006540 | 1 | 135,145 | 407.11 (402.18-412.42)s / 4572.5 MB | 111.51 (111.32-111.66)s / 76.3 MB | 126.00 (125.49-126.91)s / 194.6 MB | 3.23x / 23.5x leaner | MATCH (16031 rows); MATCH (15912 variants) |
| s2_tierp_v2_SRR15006540_th8 | SRR15006540 | 8 | 135,145 | 184.38 (179.54-187.21)s / 4087.8 MB | 16.69 (16.31-17.00)s / 339.0 MB | 78.04 (77.89-78.24)s / 224.8 MB | 2.36x / 18.2x leaner | MATCH (16031 rows); MATCH (15912 variants) |
| s2_tierp_v2_SRR8657348_th1 | SRR8657348 | 1 | 773,667 | 861.18 (857.96-863.07)s / 2028.7 MB | 174.48 (173.98-175.43)s / 143.5 MB | 154.98 (154.77-155.12)s / 474.3 MB | 5.56x / 4.3x leaner | MATCH (54448 rows); MATCH (53224 variants) |
| s2_tierp_v2_SRR8657348_th8 | SRR8657348 | 8 | 773,667 | 215.85 (212.70-217.72)s / 3697.7 MB | 24.68 (24.64-24.71)s / 251.8 MB | 39.63 (39.50-39.76)s / 493.2 MB | 5.45x / 7.5x leaner | MATCH (54448 rows); MATCH (53224 variants) |
| s2_tierp_v2_somatic_76_86_th1 | somatic_76_86 | 1 | 105,949 | 726.09 (720.47-729.27)s / 4307.4 MB | 215.16 (214.94-215.50)s / 144.9 MB | 230.53 (229.92-230.96)s / 247.8 MB | 3.15x / 17.4x leaner | MATCH (14126 rows); MATCH (14000 variants) |
| s2_tierp_v2_somatic_76_86_th8 | somatic_76_86 | 8 | 105,949 | 373.34 (371.55-374.68)s / 5679.9 MB | 31.66 (31.60-31.75)s / 786.0 MB | 135.43 (133.83-137.51)s / 382.1 MB | 2.76x / 14.9x leaner | MATCH (14126 rows); MATCH (14000 variants) |
| s2_table3_R_chr1_java_200g_n3 | SRR8657348 (-R chr1 249 Mbp) | 8 | 1 (249.25 Mbp) | 501.83 (492.93-509.65)s / 106.47 GiB (-Xmx4g OOM rc=3) | — | — | Rust 51.5x faster / 342.2x leaner | JAVA_XMX4G_OOM_RC3; JAVA_XMX200G_PASS_N3 |
| s2_table3_R_chr1_rust_and_cpp | SRR8657348 (-R chr1 249 Mbp) | 8 | 1 (249.25 Mbp) | — | 105.89 (105.07-106.58)s / 10.58 GiB | 9.75 (9.64-9.81)s / 0.31 GiB | Rust 10.9x faster / 34.0x leaner vs C++ | PASS_R_CHR1_RUST_AND_CPP_N3 |
| s2_fig1_fig3_field_callers | SRR8657348 & Somatic 76/86 | 1 & 8 | 773,667 / 105,949 | Somatic 8c: 373.34s | Somatic 8c: 31.66s | Somatic 8c: 135.43s | Somatic Venn VarDict=358 SNVs | PASS_16_FIELD_CALLER_RUNS (Rust==Java) |

*Timed n=3 runs report technical replicate mean (min-max) on identical GCP Intel Xeon Platinum 8373C Ice Lake silicon.*
