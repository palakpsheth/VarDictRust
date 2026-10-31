### Table 4. Bounded Memory, Whole-Exome 1–32 Thread Scaling, Chromosome 1 Span Invariance & Concurrent Node Packing

| Workload & Configuration | VarDictJava Wall / Peak RSS | vardictcpp Wall / Peak RSS | VarDictRust Wall / Peak RSS | Speedup / Memory Reduction Verdict |
| :--- | :--- | :--- | :--- | :--- |
| Targeted Amplicon Panel (8 Threads, n=3) | 3.07s / 120.1 MB | 5.24s / 119.7 MB | 0.01s / 11.8 MB | 307.00x faster, 10.2x leaner vs. Java |
| HG002 Germline Benchmark Slice (8 Threads, n=3) | 153.25s / 10509.5 MB | 20.10s / 860.8 MB | 21.29s / 475.7 MB | 7.20x faster, 22.1x leaner vs. Java |
| HG008 Somatic Benchmark Slice (8 Threads, n=3) | 77.28s / 17318.4 MB | 16.91s / 569.0 MB | 17.11s / 515.9 MB | 4.52x faster, 33.6x leaner vs. Java |
| SEQC2 Somatic Exome Slice (8 Threads, n=3) | 4045.00s / 16085.0 MB | 287.70s / 2245.7 MB | 277.52s / 516.7 MB | 14.58x faster, 31.1x leaner vs. Java |
| 51.8 GiB HG002 RNA-seq (-i splice, 8 Threads, n=3) | 23.31 (23.03-23.70)s / 4.79 GiB | Unsupported (-i splice unimplemented) | 4.94 (4.90-5.01)s / 359.9 MB | 4.72x faster, 13.6x leaner vs. Java |
| SEQC2 Whole Exome (1 Threads, n=3) | 5757.67 (5717.00-5789.00)s / 1.79 GiB | 1005.71 (993.97-1012.37)s / 1.21 GiB | 812.75 (805.74-825.81)s / 343.5 MB | Self-speedup S_P=1.00x (7.1x faster, 5.3x leaner vs. Java) |
| SEQC2 Whole Exome (2 Threads, n=3) | 5802.33 (5762.00-5839.00)s / 4.74 GiB | 530.70 (529.51-531.93)s / 1.29 GiB | 409.14 (405.60-411.13)s / 415.1 MB | Self-speedup S_P=1.99x (14.2x faster, 11.7x leaner vs. Java) |
| SEQC2 Whole Exome (4 Threads, n=3) | 2270.51 (2254.41-2280.43)s / 6.79 GiB | 267.32 (267.11-267.51)s / 1.37 GiB | 212.30 (211.15-213.44)s / 443.1 MB | Self-speedup S_P=3.83x (10.7x faster, 15.7x leaner vs. Java) |
| SEQC2 Whole Exome (8 Threads, n=3) | 1814.63 (1752.96-1911.57)s / 9.05 GiB | 139.14 (138.32-140.40)s / 1.49 GiB | 118.90 (117.36-119.88)s / 536.9 MB | Self-speedup S_P=6.84x (15.3x faster, 17.3x leaner vs. Java) |
| SEQC2 Whole Exome (16 Threads, n=3) | 1893.66 (1866.57-1929.42)s / 17.44 GiB | 74.17 (73.64-75.12)s / 1.69 GiB | 75.77 (74.66-76.76)s / 607.5 MB | Self-speedup S_P=10.73x (25.0x faster, 29.4x leaner vs. Java) |
| SEQC2 Whole Exome (32 Threads, n=3) | 1991.64 (1969.99-2003.54)s / 21.77 GiB | 63.68 (63.41-64.15)s / 2.06 GiB | 67.86 (67.27-68.20)s / 747.2 MB | Self-speedup S_P=11.98x (29.4x faster, 29.8x leaner vs. Java) |
| Chromosome 1 Contiguous Span 1–249 Mbp (8 Threads, n=3) | 4 GiB Heap OOM (see Table 2 for 200 GiB Heap) | 249 Mbp: 30.33s / 0.79 GiB | 249 Mbp: 27.39s / 1.38 GiB (10 Mbp: 1.17 GiB) | 249M/10M RSS ratio = 1.176x (Bounded Memory Verified) |
| Chromosome 1 Contiguous Span 1–249 Mbp (16 Threads, n=3) | 4 GiB Heap OOM (see Table 2 for 200 GiB Heap) | 249 Mbp: 18.39s / 1.14 GiB | 249 Mbp: 14.85s / 2.31 GiB (10 Mbp: 1.79 GiB) | 249M/10M RSS ratio = 1.291x (Bounded Memory Verified) |
| Chromosome 1 Contiguous Span 1–249 Mbp (32 Threads, n=3) | 4 GiB Heap OOM (see Table 2 for 200 GiB Heap) | 249 Mbp: 16.23s / 1.70 GiB | 249 Mbp: 11.87s / 3.69 GiB (10 Mbp: 2.62 GiB) | 249M/10M RSS ratio = 1.405x (Bounded Memory Verified) |
| Chromosome 1 10 Mbp All-Positions (-f 0, 8 Threads, n=3) | 100% Exact Match via Oracle Baseline | 6.01 (5.99-6.04)s / 0.37 GiB | 19.58 (19.24-20.09)s / 5.90 GiB | 100% Exact Match (MATCH (7032 variants)) |
| Concurrent 2x4-Thread SEQC2 Exome Node Packing (n=3) | — | — | Job 1: 211.38s (437.5 MB); Job 2: 211.05s (442.0 MB) | Combined Peak RSS = 0.859 GiB (Zero OOMs) |

*All multi-core and span benchmarks executed on GCP n2-standard-32 (Intel Xeon Platinum 8373C Ice Lake) with n=3 technical replicates.*
