### Table 5. Full-Genome chr1..chr22 Germline (HG002 130.8 GB) & 364x Somatic (HG008 329.0 GB) WGS Scaling (Arm A vs Arm B)

| Cell ID | Cohort & BED Geometry (-th 32) | VarDictJava v1.8.3 (-Xmx200g) | vardictcpp (--chunk 50000) | VarDictRust v1.0.6 | Parity & Full-Genome rtg vcfeval Concordance |
| :--- | :--- | :--- | :--- | :--- | :--- |
| s5_wgs_hg002_full_contigs_3way_th32 | HG002 130.8 GB WGS Arm A (22 Un-Chunked Chromosomes) | OOM Crash (rc=3) at 213.70 GB RSS (10:26.73) | 18:50.80 / 31.91 GB RSS (9,953,338 vars; +309,569 boundary dupes) | 25:16.61 (1516.61s) / 4.93 GB (4.70 GiB) RSS (9,643,769 vars) | OOM-vs-completed (Rust 6.5x less RAM than C++, 0 boundary dupes) |
| s5_wgs_hg002_full_giab_3way_th32 | HG002 130.8 GB WGS Arm B (<=50 kbp GIAB Callable BED) | 3:08:15 / 34.14 GB RSS (8,111,171 vars) | 11:40.98 / 28.35 GB RSS (8,111,565 vars) | 20:21.21 (1221.21s) / 1.43 GB (1.37 GiB) RSS (8,111,171 vars) | 100% MATCH on callable intervals (9.25x faster, 23.8x leaner vs Java); GIAB v4.2.1 F1=0.9651 (P=0.9738, R=0.9566) EQUAL |
| s5_wgs_hg008_full_contigs_3way_th32 | HG008 329.0 GB 364x T/N Somatic WGS Arm A (22 Contigs, -U) | OOM Crash (rc=3) at 214.12 GB RSS (1:20:05) | rc=124 / 91.85 GB RSS | 2:49:52 (10192.00s) / 12.14 GB (11.57 GiB) RSS (13,067,287 vars) | OOM-vs-completed (Rust 7.6x less RAM than C++) |
| s5_wgs_hg008_full_giab_3way_th32 | HG008 329.0 GB 364x T/N Somatic WGS Arm B (<=50 kbp GIAB BED, -U) | Timeout (rc=124) at 62.51 GB RSS (4:00:00) | rc=124 / 21.16 GB RSS | 1:45:55 (6355.00s) / 4.88 GB (4.65 GiB) RSS (10,195,118 vars) | OOM/Timeout-vs-completed; GIAB v0.3 F1=0.6547, v0.4 F1=0.6585 |

*HG008 small-variant full-WGS runs use -U (--nosv) to evaluate small-variant somatic calling without structural-variant realignment overhead.*
