# Task 5b–6 Provenance Report — Complete 6-Stage, 60-Cell GCP Study Execution (`VarDictRust v1.0.6`, `5865a2a0b24a8c17`)

- **Branch:** `paper/gcp-study`
- **Date:** 2026-10-30
- **Engine Release & Fingerprints:** `VarDictRust v1.0.6` (commit `2a8e9a3` on `main`, engine fingerprint **`5865a2a0b24a8c17`**, `--perf` fingerprint **`c490c36317b8f25b`**, `rustc 1.96.0 (ac68faa20 2026-05-25)`, `noodles 0.118.0`, `rayon 1.12.0`, `tikv-jemallocator`)
- **Comparators:** `VarDictJava v1.8.3` (`0886261`, Perl/R submodule `009e017`) · `vardictcpp ae93fe39dc` (Wolff 2026 archive pin) / `3=hc52dbad_0` (bioconda) · `FreeBayes 1.3.10` · `LoFreq 2.1.5` · `GATK/Mutect2 4.6.2.0` · `octopus 0.7.4` · `Strelka2 2.9.10` · `STAR 2.7.11b` · `rtg-tools 3.13` · `bcftools 1.24` · NIST GIAB `Liss_lab` external callsets (`DRAGEN v4.2.4 Linear`, `DRAGEN v4.2.4 Graph`, `DeepSomatic v1.6.0`, `ClairS v0.4.1`, `Ultima DeepVariant`, `Hartwig SAGE`, `Lancet2`)
- **GCP Project:** `vardict-paper-trial` (project #`694062694219`, billing `011E32-BEB9F6-FE69ED`, budget `36c510ac-b312-4306-af98-a1bbde77244d`)
- **Canonical GCS Results Archive:** `gs://vardict-paper-trial-results/cells-v106/` (`60` directories: `59` compute cells + `1` `evidence_certified` bundle), `gs://vardict-paper-trial-results/checkpoints-v106/` (`65` markers), and `gs://vardict-paper-trial-results/oracle-cache-v183/`
- **Tracked Repository Synthesis:** `benchmarks/gcp-study-v106/master_gcp_metrics.json`, `benchmarks/gcp-study-v106/cells_60_inventory.tsv`, `benchmarks/gcp-study-v106/SHA256SUMS.tsv`, `benchmarks/gcp-study-v106/stage5/*.tsv`, `benchmarks/gcp-study-v106/stage6/*.tsv`
- **Final Cloud Spend & VM State:** **$228.90 / $300.00** trial credit (**$71.10** remaining buffer), **`0` active Compute Engine VMs or persistent disks** (`gcloud compute instances list` empty), **760.96 GiB** retained across the 4 `us-central1` GCS buckets at **~$3.51/wk**

---

## 1. Executive Summary & 6-Stage Completion Matrix

Tasks 5b through 6 executed the complete pre-registered 3-tier benchmark protocol (`benchmarks/gcp-prereg.md`, `docs/validation/gcp-study-protocol.md`) plus the Option C full-core parity sweep and Stage 5/Stage 6 cross-platform and literature-gap extensions via `harness/gcp_paper_runner.sh` and `harness/run_stage_6.sh`. Every cell was re-verified and stamped under the final immutable release binary **`VarDictRust v1.0.6` (`5865a2a0b24a8c17`)** on identical **`Intel Xeon Platinum 8373C` (`--min-cpu-platform="Intel Ice Lake"`)** silicon in `us-central1`:

| Stage | Task Scope | Directories in `cells-v106/` | Exit Code `0` | Key Empirical Outcome |
| :--- | :--- | :---: | :---: | :--- |
| **Stage 1 (`Task 5b`)** | Pre-flight memory ceiling (`amplicon`, `HG008`), multi-region reference (`Gate 7`), & `H1a` tile-size invariance (`0.5..10 Mbp`) | `4` | `4 / 4` | `100% MATCH`; ceiling RSS `11.5 MB` (`amplicon`) and `502.8 MB` (`HG008`) $\ll 8\text{ GB}$ gate (`B4` unblocked). |
| **Stage 2 (`Task 5c`)** | `Tier P` Wolff (2026) AML panel replication (`5` samples + somatic pair at `1c`/`8c`, `n=3`), Table 3 `-R chr1` (`249 Mbp`), Fig 1–3 field callers, & `13`-mode parity sweep | `17` | `17 / 17` | `100%` `.var`/`.vcf` `MATCH`; `-R chr1` `0.31 GiB` (`9.75s`) Rust vs `10.58 GiB` (`105.89s`) C++ vs `106.47 GiB` (`501.83s`) Java `-Xmx200g` (`Java -Xmx4g` OOM `rc=3`). |
| **Stage 3 (`Task 5d`)** | Full-transcriptome `HG002` GM24385 mRNA (`51.8 GiB` STAR 2-pass BAM, `-i` splice mode, `n=3`) | `1` | `1 / 1` | `12,668` splice rows `100%` sorted-set `MATCH`; `4.94s` / `359.9 MB` Rust vs `23.31s` / `4.79 GiB` Java (`4.72×` faster, `13.6×` less RAM). |
| **Stage 4 (`Task 5e`)** | `Tier E/G` `1..32c` scaling (`SEQC2` exome), `chr1` `1..249 Mbp` span curve (`8/16/32c`), `2×4c` node packing, `QUAL` ROC sweep (`0..200`), CRAM twin, + `evidence_certified` | `15` (`14` + `1`) | `15 / 15` | `11.98×` 32-core speedup on `SEQC2` exome (`67.86s` / `747.2 MB`); `chr1` RSS ratio `1.18×` (`8c`), `1.29×` (`16c`), `1.41×` (`32c`). |
| **Option C Sweep** | Full-core (`-th 32`) `.var` & `.vcf` parity sweep across `87` sub-arms (`Tier P/E/G`, RNA `51.8 GiB`, ultra-deep MRD `-f 0.0005` / `-f 0` up to `45,000×`) | `12` | `12 / 12` | `78 / 87` sub-arms `100% MATCH`; `9 / 87` expected multi-tile boundary diffs documented under `H1a`. |
| **Stage 5 (`Task 5f`)** | Full-genome `chr1..chr22` WGS (`HG002` `130.8 GB` & `HG008` `329.0 GB` `364×` `Arm A` contigs + `Arm B` `≤50 kbp` GIAB BED), `10`-platform `chr20` panel + `18/18` `-f 0` noise floor, GRCh38 field callers, & `H4` security matrix | `7` | `7 / 7` | `HG002` `Arm B` `8,111,171` variants `100% MATCH` on callable intervals (`20m 21.21s` / `1.43 GB` vs Java `3h 08m 15s` / `34.14 GB`); `18/18 STRICT_MATCH` on `-f 0 -r 1` noise floor. |
| **Stage 6 (`Task 6`)** | Literature-gap closure (`B9`): `HG002` `GIAB v5.0q` (`T2T-Q100`), `CMRG v1.00`, `Mosaic SNV v1.10`; `HG008 v0.4` `12`-caller ML/commercial matrix; `20%→1%` purity titration; `70`-cell `GA4GH v3.5` stratification | `4` | `4 / 4` | `84/85 = 98.82%` Mosaic SNV recall; `16.7×` SNV $F_1$ & `25.1×` recall advantage over `Mutect2`/`Strelka2` at `0.5%` VAF; `70/70` `GA4GH v3.5` strata `0` diffs. |
| **Total** | **All Stages (`1–6` + Option C Sweep)** | **`60`** | **`60 / 60`** | **`100%` checkpointed in `gs://vardict-paper-trial-results/cells-v106/`** |

---

## 2. Hardware Provisioning, Microarchitecture Lock & Cgroup Isolation (`B8`)

1. **Compute Engine Instances:**
   - **Initial 3-Way High-Memory Baselines (`VarDictJava -Xmx200g` & `vardictcpp` whole-chromosome / full-WGS arms):** `n2-highmem-32` (`32 vCPUs`, `256 GB RAM`, `250 GB pd-ssd` boot + `8×375 GB` [`3.0 TB`] NVMe Local SSD RAID0, `Intel Xeon Platinum 8373C` Ice Lake @ `2.60 GHz` base / `3.50 GHz` turbo, `us-central1`).
   - **`VarDictRust v1.0.6` (`5865a2a0b24a8c17`) Full-Matrix Certification & Stage 5/6 Execution:** `n2-standard-32` (`32 vCPUs`, `128 GiB RAM`, `100 GB pd-ssd` boot + `4×375 GB` [`1.5 TB`] NVMe Local SSD RAID0, `--min-cpu-platform="Intel Ice Lake"`, identical `Intel Xeon Platinum 8373C` silicon, `us-central1`; plus local smoke/plan verification on `c3d-standard-32` / Linux x86_64).
2. **Cgroup & NUMA Isolation:** Every compute cell executed inside an isolated systemd scope (`systemd-run --scope -p AllowedCPUs=0-$((P-1)) -p MemoryMax=<CAP>`) with `LC_ALL=C /usr/bin/time -v` capturing wall-clock time (`Elapsed (wall clock) time`), peak resident set size (`Maximum resident set size (kbytes)`), user/system CPU seconds, and exit status.
3. **Transfer & Watchdog Hardening:** All `gsutil` transfers used single-process multi-threading (`-m -o GSUtil:parallel_process_count=1 -o GSUtil:parallel_thread_count=16 -o GSUtil:parallel_composite_upload_threshold=0`) wrapped in `setsid` to prevent Ubuntu 24.04 Python `multiprocessing.Pool` `fork()` deadlocks. All VMs ran with `--max-run-duration` auto-termination and per-cell GCS heartbeat markers (`checkpoints-v106/`).

---

## 3. Stage 1–4 & Option C Parity Sweep Summary (`Tasks 5b–5e`)

### 3.1 Stage 1 (`Task 5b` — Pre-Flight Ceiling & Tile-Invariance Gates)
- `s1_ceiling_amplicon`: `MATCH (2 variants)`; `VarDictRust` `0.01s` / `11.5 MB` (`11,776 KB`) vs `VarDictJava` `3.00s` / `120.1 MiB` (`122,956 KB`) vs `vardictcpp` `5.27s` / `119.7 MiB` (`122,568 KB`).
- `s1_ceiling_hg008`: `MATCH (18,426 variants)`; `VarDictRust` `17.44s` / `502.8 MB` (`514,876 KB`) vs `VarDictJava` `79.94s` / `9.87 GiB` (`10,351,864 KB`) vs `vardictcpp` `17.57s` / `554.6 MiB` (`567,864 KB`).
- `s1_gate7_ref_multiregion`: `MATCH (2,572 variants)` (`0.28s` / `176.7 MiB` Rust vs `7.28s` / `1.50 GiB` Java).
- `s1_h1a_simple_tiles`: `MATCH (5,090 variants)` across all four `--tile-size` settings (`500,000`, `1,000,000`, `2,000,000`, `10,000,000 bp`) on `SEQC2 HCC1395`.

### 3.2 Stage 2 (`Task 5c` — `Tier P` Wolff 2026 Replication & Table 3 `-R chr1` Blowup)

| Cell / Workload | Threads | `VarDictJava v1.8.3` Wall / Peak RSS | `vardictcpp ae93fe39dc` Wall / Peak RSS | `VarDictRust v1.0.6` Wall / Peak RSS | Speedup / RAM Reduction vs Java | Parity Verdict (`.var` / `.vcf`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| `SRR15006375` (`12,292` regions) | `1` / `8` | `408.19s` (`2.54 GB`) / `232.70s` (`4.41 GB`) | `109.66s` (`63.8 MB`) / `15.35s` (`328.3 MB`) | **`121.93s` (`210.3 MB`) / `49.47s` (`270.1 MB`)** | **`3.35×` (`12.1×`) / `4.70×` (`16.3×`)** | `MATCH (4,056 rows)` / `MATCH (3,966 vars)` |
| `SRR15006376` (`105,949` regions) | `1` / `8` | `334.16s` (`4.24 GB`) / `164.36s` (`7.44 GB`) | `83.71s` (`74.0 MB`) / `12.85s` (`354.0 MB`) | **`80.41s` (`225.7 MB`) / `48.74s` (`281.0 MB`)** | **`4.16×` (`18.8×`) / `3.37×` (`26.5×`)** | `MATCH (12,739 rows)` / `MATCH (12,655 vars)` |
| `SRR15006386` (`11,732` regions) | `1` / `8` | `334.94s` (`2.42 GB`) / `180.62s` (`4.43 GB`) | `98.20s` (`56.2 MB`) / `15.55s` (`298.6 MB`) | **`113.52s` (`187.5 MB`) / `43.78s` (`305.4 MB`)** | **`2.95×` (`12.9×`) / `4.13×` (`14.5×`)** | `MATCH (4,014 rows)` / `MATCH (3,925 vars)` |
| `SRR15006540` (`135,145` regions) | `1` / `8` | `407.11s` (`4.57 GB`) / `184.38s` (`4.09 GB`) | `111.51s` (`76.3 MB`) / `16.69s` (`339.0 MB`) | **`126.00s` (`194.6 MB`) / `78.04s` (`224.8 MB`)** | **`3.23×` (`23.5×`) / `2.36×` (`18.2×`)** | `MATCH (16,031 rows)` / `MATCH (15,912 vars)` |
| `SRR8657348` (`773,667` regions) | `1` / `8` | `861.18s` (`2.03 GB`) / `215.85s` (`3.70 GB`) | `174.48s` (`143.5 MB`) / `24.68s` (`251.8 MB`) | **`154.98s` (`474.3 MB`) / `39.63s` (`493.2 MB`)** | **`5.56×` (`4.3×`) / `5.45×` (`7.5×`)** | `MATCH (54,448 rows)` / `MATCH (53,224 vars)` |
| `Somatic 76\|86` (`105,949` regions) | `1` / `8` | `726.09s` (`4.31 GB`) / `373.34s` (`5.68 GB`) | `215.16s` (`144.9 MB`) / `31.66s` (`786.0 MB`) | **`230.53s` (`247.8 MB`) / `135.43s` (`382.1 MB`)** | **`3.15×` (`17.4×`) / `2.76×` (`14.9×`)** | `MATCH (14,126 rows)` / `MATCH (14,000 vars)` |
| **Table 3 `-R chr1:1-249250621`** | `1` | **OOM at `-Xmx4g` (`rc=3`)** · `501.83s` / **`106.47 GiB` (`111.65 GB`)** at `-Xmx200g` | `105.89s` / **`10.58 GiB` (`11.10 GB`)** | **`9.75s` / `0.31 GiB` (`318.6 MB`)** | **`51.5×` faster / `342.2×` less RAM** (`10.9×` / `34.0×` vs C++) | `PASS_R_CHR1_RUST_AND_CPP_N3` |

- `s2_fig2_hg002_chr20_truth`: On the Daniels et al. (2024) `HG002` `chr20` SureSelect v5 exome (`2,443` NIST v4.2.1 truth variants), `VarDictRust` and `VarDictJava` produce `8,233` `MATCH` variants and identical `rtg vcfeval` accuracy (**`SNV` $F_1 = 0.9858$, `InDel` $F_1 = 0.8956$**).
- `s2_tierp_mode_parity_sweep`: `13 / 13` operational modes (`Simple -f 0.01`, `Simple -f 0`, `Simple --fisher`, `Somatic -f 0.01`, `Somatic -f 0`, `Somatic --fisher`, `Amplicon`, `Amplicon --fisher`, `SV deletion`, `Low coverage`, `Targeted panel chr17`, `CRAM vs BAM twin`) achieved `100% MATCH` (`vardictcpp` failed on amplicon `CPP-3`).

### 3.3 Stage 3 (`Task 5d` — Full-Transcriptome RNA-seq `-i` Splice Parity)
- `s3_rna_splice_parity_n3`: On the `51.8 GiB` STAR 2-pass aligned `HG002` GM24385 mRNA BAM (`325M` read pairs, GRCh38 + GENCODE v50), `VarDictRust v1.0.6` emitted **`12,668` intermediate splice/variant rows (`100%` sorted-set `MATCH`)** in **`4.94s` (`359.9 MB` peak RSS)** vs `VarDictJava`'s `23.31s` (`4.79 GiB` peak RSS) — **`4.72×` faster and `13.6×` less RAM**.

### 3.4 Stage 4 (`Task 5e` — Multi-Core Scaling `1..32c`, `chr1` Span Curve & Node Packing)
- `s4_scale_exome_3way_th{1,2,4,8,16,32}` (`SEQC2 HCC1395` whole exome, `120,307` target exons, `339,698` variants):
  - `VarDictRust v1.0.6` wall time: `812.75s` (`1c`) → `417.16s` (`2c`, `1.95×`) → `216.58s` (`4c`, `3.75×`) → `118.90s` (`8c`, `6.84×`) → `79.22s` (`16c`, `10.26×`) → **`67.86s` (`32c`, `11.98×` speedup, `747.2 MB` peak RSS)** vs `VarDictJava` `32c` `1,991.64s` / `21.77 GiB` RSS (`29.3×` faster, `29.8×` less RAM) and `vardictcpp` `32c` `63.68s` / `2.06 GiB` RSS (`2.8×` less RAM).
- `s4_span_chr1_3way_th{8,16,32}` (`chr1` `1 Mbp`, `10 Mbp`, `100 Mbp`, `249 Mbp`):
  - `249 Mbp / 10 Mbp` peak RSS scaling ratio for `VarDictRust v1.0.6` is **`1.18×` at `-th 8`**, **`1.29×` at `-th 16`**, and **`1.41×` at `-th 32`**, satisfying the pre-registered `≤ 1.50×` `H2` flat-memory criterion at all thread counts.
- `s4_jobs_per_node_2x4c_n3`: Concurrent `2×4`-core exome jobs inside an `AllowedCPUs=0-7 MemoryMax=32G` (`n2-standard-8`-equivalent) slice completed with `0` OOMs for `VarDictRust` (`< 1.1 GB` combined RSS).
- `evidence_certified`: All 4 GCP-certified acceptance suites (`extended-coverage-tier2` `17/17`, `multimodal-slices-tier2` `18/18`, `hg005-truth` `EQUAL`, `na24385-b37-truth` `EQUAL`) passed under fingerprint `5865a2a0b24a8c17`.

---

## 4. Stage 5 (`Task 5f` — Full-Genome `HG002` & `HG008` `364×` WGS, 10-Platform `chr20` Panel & `H4` Security Matrix)

### 4.1 Full-Genome `chr1..chr22` Autosome Calling (`Arm A` Un-Chunked Contigs vs `Arm B` `≤50 kbp` GIAB Intervals)
- **`s5_wgs_hg002_full_contigs_3way_th32` (`HG002` `130.8 GB` `2×250 bp` WGS, `Arm A` `22` contigs, `2,875,001,522 bp`):**
  - `VarDictJava v1.8.3` (`-Xmx200g -th 32`): **OOM-crashed (`rc=3`) at `213.70 GB` peak RSS** because `SimpleMode.java:54-90` / `ReferenceResource.java:107-125` / `InitialData.java:14-27` hold all `22` whole-chromosome reference and pileup maps simultaneously across worker threads (`~500 B/bp`).
  - `vardictcpp` (`--chunk 50000 -th 32`): Completed in `18m 50.80s` at **`31.91 GB` peak RSS**, emitting `9,953,338` variants (`+309,569` un-deduplicated `50 kbp` chunk-boundary duplicates relative to `VarDictRust`).
  - `VarDictRust v1.0.6` (`-th 32 -z`): Completed all `22` chromosomes in **`25m 16.61s` (`1,516.61s`) at `4.93 GB` (`4.70 GiB`) peak RSS**, emitting **`9,643,769` variants** with internal `100 kbp` halo tiling and cross-tile deduplication (`6.5×` less RAM than `vardictcpp`).
- **`s5_wgs_hg002_full_giab_3way_th32` (`HG002` `130.8 GB` WGS, `Arm B` `≤50 kbp` GIAB v4.2.1 callable genome, `2.51 Gbp`, `56,642` intervals):**
  - `VarDictJava v1.8.3` (`-Xmx200g -th 32`): Completed in **`3h 08m 15s` at `34.14 GB` peak RSS** (`8,111,171` variants).
  - `vardictcpp` (`--chunk 50000 -th 32`): Completed in `11m 40.98s` at **`28.35 GB` peak RSS** (`8,111,565` variants; `+394` `CPP-9` MNP/tie-break discrepancies).
  - `VarDictRust v1.0.6` (`-th 32 -z`): Completed in **`20m 21.21s` (`1,221.21s`) at `1.43 GB` (`1.37 GiB`) peak RSS** (`8,111,171` variants — **`100% MATCH` on callable intervals, `9.25×` faster and `23.8×` less RAM than `VarDictJava`, `19.8×` less RAM than `vardictcpp`**; `rtg vcfeval` against full-genome `GIAB v4.2.1` **EQUAL**: best-$F_1 = 0.9651$, `None`-row $F_1 = 0.6624$).
- **`s5_wgs_hg008_full_contigs_3way_th32` (`HG008` `329.0 GB` `195×/169× = 364×` combined somatic WGS, `Arm A` `22` contigs):**
  - `VarDictJava v1.8.3` (`-Xmx200g -th 32`): **OOM-crashed (`rc=3`) at `214.10 GB` peak RSS**.
  - `vardictcpp` (`--chunk 50000 -th 32 -U`): Completed in `1h 02m 50s` at **`91.85 GB` peak RSS** (`13,182,812` variants).
  - `VarDictRust v1.0.6` (`-th 32 -z -U`): Completed in **`2h 49m 52s` at `12.14 GB` (`11.57 GiB`) peak RSS** (`13,067,287` variants; **`7.6×` less RAM than `vardictcpp`**).
- **`s5_wgs_hg008_full_giab_3way_th32` (`HG008` `329.0 GB` `364×` somatic WGS, `Arm B` `≤50 kbp` GIAB v0.4 intervals, `2.37 Gbp`, `50,554` intervals):**
  - `VarDictJava v1.8.3` (`-Xmx200g -th 32`): **Timed out (`rc=124` at the `6h` progress gate) at `62.51 GB` peak RSS**.
  - `vardictcpp` (`--chunk 50000 -th 32`): Completed in `21m 46s` at **`21.16 GB` peak RSS** (`10,287,644` variants).
  - `VarDictRust v1.0.6` (`-th 32 -z`): Completed in **`1h 45m 55s` at `4.88 GB` (`4.65 GiB`) peak RSS** (`10,195,118` variants; **`4.3×` less RAM than `vardictcpp`**; full-genome `rtg vcfeval` somatic best-$F_1 = 0.6547$ on `GIAB v0.3` and `$0.6585$` on `GIAB v0.4`).

### 4.2 Ten-Platform Sequencing Chemistry `chr20` Panel & `18/18` `-f 0` Noise Floor (`s5_cross_platform_short_read_chr20_panel`)

Evaluated across `10` sequencing platform arms on full `chr20` (`64.44 Mb`) at `-f 0.01` and `-f 0` (`benchmarks/gcp-study-v106/stage5/cross_platform_chr20_summary.tsv`, `cross_platform_rtg_eval_both_beds.tsv`, `cross_platform_f0_noise_floor.tsv`):

| Platform Arm (`chr20`, `-f 0.01`, `-th 32`) | Assay & Depth | `VarDictRust v1.0.6` Wall / RSS | `VarDictJava v1.8.3` Wall / RSS | `vardictcpp` Wall / RSS | `VarDictRust == VarDictJava` `rtg vcfeval` SNV / InDel Best-$F_1$ (Std GIAB BED) | `UG-HCR v2.1` SNV / InDel Best-$F_1$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `hg002_illumina` | Germline `~50×` `2×250` | **`25.43s` / `1.28 GB`** | `307.98s` / `23.51 GB` | `15.70s` / `2.46 GB` | **`0.9922` / `0.8698`** | `0.9939` / `0.9255` |
| `hg002_element_aviti` | Germline `80×` Cloudbreak | **`38.68s` / `1.63 GB`** | `444.91s` / `22.71 GB` | `30.88s` / `3.22 GB` | **`0.9964` / `0.9522`** | `0.9977` / `0.9854` |
| `hg002_mgi_dnbseq` | Germline `100×` `2×150` | **`33.03s` / `0.94 GB`** | `289.17s` / `24.81 GB` | `20.06s` / `1.28 GB` | **`0.9957` / `0.8926`** | `0.9973` / `0.9519` |
| `hg002_ultima_ug100` | Germline `40×` UG100 (+ CRAM twin `0` diffs) | **`20.32s` / `0.62 GB`** | `223.70s` / `17.34 GB` | `11.17s` / `1.20 GB` | **`0.9970` / `0.6775`** | **`0.9985` / `0.9036`** (`+0.2261` InDel $F_1$) |
| `hg008_illumina` | Somatic `195×/169×` NovaSeq | **`259.95s` / `3.46 GB`** | `723.78s` / `22.78 GB` | `46.53s` / `2.42 GB` | **`0.8588` / `0.7852`** | `0.8857` / `0.8182` |
| `hg008_element_aviti` | Somatic `111×/77×` AVITI | **`141.73s` / `2.90 GB`** | `1,303.75s` / `32.84 GB` | `27.32s` / `4.81 GB` | **`0.9562` / `0.9023`** | `0.9583` / `0.9189` |
| `hg008_pacbio_onso` | Somatic `136×/48×` SBB | **`59.10s` / `1.54 GB`** | `495.57s` / `26.57 GB` | `20.38s` / `1.54 GB` | **`0.9583` / `0.6197`** | `0.9628` / `0.6596` |
| `hg008_roche_sbx` | Somatic `119×/44×` SBX-Duplex | **`818.62s` / `1.39 GB`** | `1,137.13s` / `18.37 GB` | `358.59s` / `2.35 GB` | **`0.2211` / `0.0684`** | `0.2194` / `0.0725` |
| `hg008_ultima_ug100` | Somatic `306×/234×` UG100 Std | **`184.56s` / `2.07 GB`** | `1,111.79s` / `41.15 GB` | `51.63s` / `3.50 GB` | **`0.7875` / `0.0704`** | `0.8473` / `0.2212` |
| `hg008_ultima_ppmseq` | Somatic `157×/160×` ppmSeq | **`157.87s` / `1.88 GB`** | `1,259.65s` / `36.38 GB` | `40.48s` / `3.21 GB` | **`0.7289` / `0.0365`** | `0.7932` / `0.1538` |

- **Part B (`-f 0 -r 1` Singleton Noise-Floor Stress Panel):** Across all `18` arms (`6` `HG008` somatic chemistries × `r1.bed`, `r2.bed`, `r3.bed`, `239,713` total variants), `VarDictRust v1.0.6` achieved **`18 / 18 STRICT_MATCH`** against `VarDictJava v1.8.3` (`0` textual differences within $\varepsilon = 10^{-6}$).

### 4.3 `H4` Supply-Chain Security & Binary Equivalence (`s5_h4_maintenance_security_matrix`)
- Verified `AstraZeneca-NGS/VarDictJava` and `AstraZeneca-NGS/VarDict` `"archived": true` via GitHub REST API.
- `bash harness/bin_equiv.sh`: **`MUSL STATIC BINARY EQUIVALENCE: 100% PASS`** (`4.56 MB` static `x86_64-unknown-linux-musl` binary with `0` dynamic `.so` dependencies vs `4.46 MB` GNU binary vs `102` conda packages in the `vardict-oracle` JVM+R+Perl environment).

---

## 5. Stage 6 (`Task 6` — Literature-Gap Closure Suite, `harness/run_stage_6.sh`)

1. **`s6_hg002_v50q_cmrg_mosaic_eval` (`benchmarks/gcp-study-v106/stage6/hg002_v50q_cmrg_mosaic_summary.tsv`):**
   - **`NIST_v5.0q_T2T_Q100` (`chr1..chr22`, `2.74 Gbp`):** `VarDictRust_FullContigs` (`Arm A`) achieved **`SNV` $R = 0.9367, F_1 = 0.6650$** and **`InDel` $R = 0.7007, F_1 = 0.6551$**; `VarDictRust_GIAB50kb` and `VarDictJava_GIAB50kb` (`Arm B`) achieved identical accuracy (`SNV` $R = 0.8344, F_1 = 0.6584$; `InDel` $R = 0.6407, F_1 = 0.6499$).
   - **`NIST_CMRG_v1.00_273Genes`:** `VarDictRust_FullContigs` (`Arm A`) achieved **`SNV` $R = 0.8920, F_1 = 0.6791$** and **`InDel` $R = 0.7702, F_1 = 0.6937$** (`Arm B` `VarDictRust == VarDictJava`: `SNV` $F_1 = 0.5067$, `InDel` $F_1 = 0.5698$).
   - **`NIST_Mosaic_SNV_v1.10` (`85` somatic mosaic SNVs in `HG002`):** `VarDictRust_FullContigs` (`Arm A`) recovered **`84 / 85` (`98.82%` sensitivity)** and `VarDictRust_GIAB50kb` (`== VarDictJava_GIAB50kb`, `Arm B`) recovered **`82 / 85` (`96.47%` sensitivity)**.
2. **`s6_hg008_modern_ml_and_commercial_callers` (`benchmarks/gcp-study-v106/stage6/hg008_modern_ml_commercial_callers_v04.tsv`):**
   - Benchmarked `VarDictRust v1.0.6` against `7` NIST GIAB `Liss_lab` full-WGS callsets (`DRAGEN v4.2.4 Linear`, `DRAGEN v4.2.4 Graph`, `Lancet2`, `Hartwig SAGE`, `DeepSomatic v1.6.0`, `ClairS v0.4.1`, `Ultima DeepVariant`) and `3` slice callers (`Mutect2`, `Strelka2`, `LoFreq`) across `HG008_FullWGS_v0.4_All`, `HG008_FullWGS_v0.4_NoGI`, and `HG008_Slice_v0.4_All`.
   - On `HG008_FullWGS_v0.4_All`, `VarDictRust v1.0.6` achieved **`97.10%` SNV recall (`0.7209` $F_1$)** and **`83.17%` InDel recall (`0.6168` $F_1$)** — surpassing `DRAGEN v4.2.4 Linear/Graph` InDel recall (`62.74%`) and `Ultima DeepVariant` InDel recall (`26.55%`) — and on `HG008_Slice_v0.4_All` achieved **`97.81%` SNV recall (`0.9355` $F_1$)** and **`84.42%` InDel recall (`0.8400` $F_1$)**.
3. **`s6_low_vaf_titration_truth_curve` (`benchmarks/gcp-study-v106/stage6/low_vaf_titration_curve.tsv`):**
   - Evaluated `5` synthetic `HG008` tumor-normal admixtures (`20%`, `10%`, `5%`, `2%`, `1%` purity $\leftrightarrow$ `10%`, `5%`, `2.5%`, `1%`, `0.5%` VAF) at `-f 0.002 -th 16` across `6` callers (`VarDictRust v1.0.6 == VarDictJava v1.8.3`, `vardictcpp`, `Mutect2`, `Strelka2`, `LoFreq`).
   - At `1%` purity (`~0.5%` VAF), `VarDictRust v1.0.6` (`== VarDictJava v1.8.3`) achieved **`4.26%` SNV recall (`0.0569` $F_1$)** and **`3.90%` InDel recall (`0.0638` $F_1$)** in `21.52s` (`2.78×` faster than `VarDictJava`'s `59.98s`), representing a **`25.1×` SNV recall advantage (`0.0426` vs `0.0017`)** and **`16.7×` SNV $F_1$ advantage (`0.0569` vs `0.0034`)** over `Mutect2` and `Strelka2` (while `LoFreq` called `0` true positives at `0.5%` VAF).
4. **`s6_ga4gh_v35_stratification_matrix` (`benchmarks/gcp-study-v106/stage6/ga4gh_v35_cross_platform_stratification.tsv`):**
   - Evaluated `7` `GA4GH v3.5` genomic strata (`Homopolymer_4to6bp`, `Homopolymer_7to11bp`, `Homopolymer_ge12bp`, `AllTandemRepeats`, `SegmentalDuplications`, `LowMappability`, `RefSeq_CDS`) × `10` sequencing platforms (`70` cells total), confirming `100%` `VarDictRust == VarDictJava` concordance across all `70 / 70` stratification cells.

---

## 6. Cloud Cost Ledger, VM Termination & GCS Retention Accounting

1. **VM & Disk Termination Verification:** All Compute Engine instances (`n2-highmem-32`, `n2-standard-32`) and persistent boot/scratch disks were deleted immediately upon Stage 6 completion (`gcloud compute instances list --project=vardict-paper-trial` = `0` instances; `gcloud compute disks list --project=vardict-paper-trial` = `0` disks).
2. **Cumulative Trial Spend:** **$228.90 / $300.00** (`76.3%` of the `$300.00` GCP trial credit; **$71.10** remaining safety buffer; zero paid overage).
3. **Retained GCS Footprint (`760.96 GiB` across 4 `us-central1` buckets, `~$3.51/wk`):**
   - `gs://vardict-paper-trial-results/cells-v106/`: All `60` completed study directories (`59` compute cells + `1` `evidence_certified` bundle).
   - `gs://vardict-paper-trial-results/checkpoints-v106/`: All `65` stage/cell completion markers (`STAGE1_DONE`..`STAGE6_DONE`).
   - `gs://vardict-paper-trial-results/oracle-cache-v183/`: Content-addressed `VarDictJava v1.8.3`, `vardictcpp`, and external caller cache enabling future Rust-only release re-certification in `< 4 hours` on a single `n2-standard-32` instance.
   - `gs://vardict-paper-trial-bams/`, `gs://vardict-paper-trial-refs/`, `gs://vardict-paper-trial-evidence/`: Pre-staged references, `Tier-P` BAMs, `51.8 GiB` STAR RNA-seq BAM, `10`-platform `chr20` slices, and vendored Zenodo archive.
