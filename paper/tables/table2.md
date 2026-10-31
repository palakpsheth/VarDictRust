### Table 2. Direct Replication of vardictcpp (Wolff et al. 2026) Across Clinical AML Panels, Field Callers & Unchunked Chromosome 1 (-R chr1)

| Sample / Cohort | Threads | BED Intervals | VarDictJava Wall (s) / Peak RSS | vardictcpp Wall (s) / Peak RSS | VarDictRust Wall (s) / Peak RSS | Rust vs. Java Speedup / Memory Reduction | Parity Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SRR15006375 | 1 | 12,292 | 408.19 (400.74-418.94)s / 2537.2 MB | 109.66 (109.19-109.90)s / 63.8 MB | 121.93 (121.14-122.38)s / 210.3 MB | 3.35x faster / 12.1x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006375 | 8 | 12,292 | 232.70 (229.14-235.27)s / 4410.5 MB | 15.35 (15.30-15.42)s / 328.3 MB | 49.47 (49.39-49.59)s / 270.1 MB | 4.70x faster / 16.3x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006376 | 1 | 105,949 | 334.16 (326.10-340.26)s / 4244.4 MB | 83.71 (83.03-84.32)s / 74.0 MB | 80.41 (80.08-80.84)s / 225.7 MB | 4.16x faster / 18.8x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006376 | 8 | 105,949 | 164.36 (150.76-174.29)s / 7442.2 MB | 12.85 (12.50-13.09)s / 354.0 MB | 48.74 (48.65-48.86)s / 281.0 MB | 3.37x faster / 26.5x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006386 | 1 | 11,732 | 334.94 (331.87-338.80)s / 2423.4 MB | 98.20 (98.09-98.37)s / 56.2 MB | 113.52 (113.17-114.04)s / 187.5 MB | 2.95x faster / 12.9x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006386 | 8 | 11,732 | 180.62 (178.99-182.17)s / 4434.7 MB | 15.55 (14.95-16.69)s / 298.6 MB | 43.78 (43.64-43.88)s / 305.4 MB | 4.13x faster / 14.5x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006540 | 1 | 135,145 | 407.11 (402.18-412.42)s / 4572.5 MB | 111.51 (111.32-111.66)s / 76.3 MB | 126.00 (125.49-126.91)s / 194.6 MB | 3.23x faster / 23.5x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006540 | 8 | 135,145 | 184.38 (179.54-187.21)s / 4087.8 MB | 16.69 (16.31-17.00)s / 339.0 MB | 78.04 (77.89-78.24)s / 224.8 MB | 2.36x faster / 18.2x leaner | 100% Exact Match (.var & .vcf) |
| SRR8657348 | 1 | 773,667 | 861.18 (857.96-863.07)s / 2028.7 MB | 174.48 (173.98-175.43)s / 143.5 MB | 154.98 (154.77-155.12)s / 474.3 MB | 5.56x faster / 4.3x leaner | 100% Exact Match (.var & .vcf) |
| SRR8657348 | 8 | 773,667 | 215.85 (212.70-217.72)s / 3697.7 MB | 24.68 (24.64-24.71)s / 251.8 MB | 39.63 (39.50-39.76)s / 493.2 MB | 5.45x faster / 7.5x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006376 / SRR15006386 (Somatic Pair) | 1 | 105,949 | 726.09 (720.47-729.27)s / 4307.4 MB | 215.16 (214.94-215.50)s / 144.9 MB | 230.53 (229.92-230.96)s / 247.8 MB | 3.15x faster / 17.4x leaner | 100% Exact Match (.var & .vcf) |
| SRR15006376 / SRR15006386 (Somatic Pair) | 8 | 105,949 | 373.34 (371.55-374.68)s / 5679.9 MB | 31.66 (31.60-31.75)s / 786.0 MB | 135.43 (133.83-137.51)s / 382.1 MB | 2.76x faster / 14.9x leaner | 100% Exact Match (.var & .vcf) |
| SRR8657348 (Unchunked -R chr1, 249 Mbp) | 8 | 1 (249.25 Mbp) | 501.83 (492.93-509.65)s / 106.47 GiB (4 GiB Heap = OOM) | — | — | Rust 51.5x faster / 342.2x leaner | 4 GiB Heap OOM Crash; Completes at 200 GiB Heap (n=3) |
| SRR8657348 (Unchunked -R chr1, 249 Mbp) | 8 | 1 (249.25 Mbp) | — | 105.89 (105.07-106.58)s / 10.58 GiB | 9.75 (9.64-9.81)s / 0.31 GiB | Rust 10.9x faster / 34.0x leaner vs. C++ | Both Complete in Bounded Memory (n=3) |
| SRR8657348 & SRR15006376/86 (7-Caller Panel) | 1 & 8 | 773,667 / 105,949 | Somatic 8c: 373.34s | Somatic 8c: 31.66s | Somatic 8c: 135.43s | Somatic Venn: 358 Concordant SNVs | 100% Exact Match Across All 16 Runs (Rust == Java) |

*Timed n=3 runs report technical replicate mean (min-max) on identical GCP Intel Xeon Platinum 8373C Ice Lake silicon.*
