# HG008 V0.4 Draft Benchmark + Ladder Artifacts — Data Guide

Truth and timing artifacts behind the slice-to-WGS estimator and the three NIST
example use-cases (plan: `2026-10-15-hg008-wgs-estimate-ladder`). All runs are
compute-only `vardictrust` on a 20T/31GB/no-GPU box — no oracle involved.

## 1. V0.4 draft truth (`scratch/g4/hg008/truth/`, tracked in `harness/assets.tsv`)

Source: `.../HG008/Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.4-20261006/GRCh38/`
(Oct 2026; **draft** — assembly-based tumor-vs-normal vcfeval, GRCh38).

| File | Size | md5 | Purpose |
|---|---|---|---|
| `..._tumorvariants.vcf.gz` (+`.tbi`) | 280K | `d5b0321c…` | Somatic truth variants (18,871 genome-wide) |
| `..._all.bed` | 885K | `2e5ddf4b…` | Benchmark eval regions (37,937 intervals) |
| `..._nogermlineinterference.bed` | 973K | `20872f86…` | Stratum excluding germline-interference zones |

Left out deliberately: `..._nogermlinewithin50bp.bed` (77M) and
`..._full_tumornormal.vcf.gz` (26M) — pull only if a stratum needs them.

## 2. Ladder artifacts (`scratch/g4/hg008/ladder/`)

Reusable inputs and evidence (all timings `/usr/bin/time -v`, flags
`-f 0.01 -c 1 -S 2 -E 3 -g 4 -z`):

| Artifact | Notes |
|---|---|
| `full_th8.vcf` (7.7M, 18426 vars) | Paired 0823p23 slices over `engine.bed` (4.028 Mbp), 38.5 s @ `-th 8` |
| `tumoronly_chr1.vcf` (1335 vars) | Single-sample chr1 subset (525 kb), 2.45 s |
| `tumoronly_nognomad.vcf` (364 vars) | Above minus any gnomAD v4.1 presence (971 dropped) |
| `p2T.chr1.bam` (+`.bai`, 21M) | Remote slice of p2 CRAM over chr1 subset, depth 118.5x |
| `p2T_chr1.vcf` (1364 vars) | Single-sample p2 calls, 2.26 s |
| `eval_*.bed` | `engine.bed ∩ V0.4` domains (eval must be intersected — calls span only slices) |
| `extra.json`, `*.time`, `v04_*_rtg.json` | Fit values and rtg summaries |

## 3. Key results

* **Estimator:** 9.57 s/Mbp @ `-th 8` (dual-BAM ~170x), RSS bounded 1.5 GB →
  same-depth WGS ≈ 8 h (6–10 h band); 60x/40x ≈ 2.8 h (2–4 h); RSS 1.5–2.5 GB.
* **V0.4 rtg (`concordance.sh`):** `all` P=0.387 R=0.696 F1=0.497 (881 truth in
  domain); `nogermline` P=0.386 R=0.902 F1=0.541.
* **Tumor-only + gnomAD (chr1, 109 truth):** F1 0.148 → 0.340 (P 0.081→0.232,
  R 0.927→0.633); presence-filtering costs 32 TP.
* **Passage drift (p2 vs 0823p23, same window/flags):** 909 shared /
  426 batch-only / 455 p2-only; 40 shared sites |ΔAF|≥0.2 (top: 1.0→0.03).
* **Edge finding:** BED vs `-R` modes agree on keys/order; 5/97 records differ
  in DP by ≤3 reads in MSI≥13 poly-A — always use BED mode for truth work.

## 4. Reproduce

```bash
# V0.4 truth (see harness/assets.tsv for URLs/md5):
bash harness/fetch_data.sh hg008   # or curl the 4 V0.4 rows into scratch/g4/hg008/truth/
# Ladder + evals (g4-tools env for rtg/bcftools/bedtools):
export PATH=/home/palak/miniconda3/envs/g4-tools/bin:$PATH
./target/release/vardictrust -G scratch/g4/ref/giabv3/GRCh38_GIABv3.fa -f 0.01 \
  -N "HG008-T|HG008-N" -b "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam" \
  -c 1 -S 2 -E 3 -g 4 -z -th 8 scratch/g4/hg008/engine.bed > full_th8.vcf
bedtools intersect -a scratch/g4/hg008/engine.bed \
  -b scratch/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed > eval.bed
bash harness/concordance.sh full_th8.vcf <truth.vcf.gz> eval.bed \
  scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf outdir "HG008-T"
```

Caveats: V0.4 is draft; single-sample VCFs carry no `STATUS` (use PASS-only
pre-filter, not `concordance.sh` as-is); chr1 tumor-only numbers are
small-domain (109 truth); precision figures are slice-domain only.
