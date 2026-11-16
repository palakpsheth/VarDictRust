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
* **p13 trajectory (p2→p13→p23, same 525 kb chr1 window):** p13 slice 96x
  (~85X claim ✓), 1268 calls in 1.1 s. 893 shared all3; 500/893 (56%) monotonic;
  30 monotonic sites |ΔAF|≥0.2 with clean fixation curves (e.g. 0.03→1.0→1.0,
  0.015→0.51→0.60) vs 10 non-monotonic — evolution signal over noise.
* **AF-threshold sweep (tumor-only chr1, 109 truth):** any-presence F1=0.340 →
  af0.01 F1=0.357 (P 0.229 / R 0.807); F1 flat ~0.356 across cutoffs while recall
  climbs 0.63→0.86 — precision is capped by non-germline FPs (needs PoN, not AF tuning).
* **Full-slice sweep (all 21 chroms, 881 truth, 1.07M gnomAD keys):** F1 flat
  ~0.31 across cutoffs (any 0.311 / af0.001 0.314 / af0.01 0.310 / af0.05 0.304);
  recall 0.54→0.73 as the cutoff loosens, precision stuck ~0.20 — confirms the
  chr1 verdict at 8x the truth: AF tuning cannot fix precision; artifact/PoN
  filtering is the lever. Full-slice tumor-only ran 13.6 s (12,400 vars).
* **V0.6 SV/CNV (BLOCKED, filed as finding):** truth pulled (368 KB, tracked in
  `assets.tsv`). truvari 5.5.0 bench of PASS somatic SV calls: 0 TP / 210 FN —
  engine SV records (`<DEL>`/`<DUP>`/`<INV>`) carry no SVLEN/END, so truvari
  discards all 5 query SVs (`comp cnt: 0`). Only 5/28 SV-class calls are even
  PASS+somatic. V0.6 eval is blocked on SV output representation, not biology;
  no engine tuning attempted here. CNV benchmarking has no established method
  per the V0.6 README — not attempted.
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
bash harness/concordance.sh full_th8.vcf \
  scratch/g4/hg008/truth/HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz \
  eval.bed scratch/g4/ref/giabv3/GRCh38_GIABv3.sdf outdir "HG008-T"
# Tumor-only + gnomAD (single-sample VCFs lack STATUS: PASS-only pre-filter,
# manual vcfeval — concordance.sh does not apply):
G=https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz
tabix -h $G <first chr1.bed row> > gnomad_chr1slice.vcf   # + remaining rows; 745 MB
python3 -c "drop calls whose (chrom,pos,ref,alt) is in gnomAD (any-presence or AF>=cutoff)" \
  > tumoronly_nognomad.vcf
# Passage slices (remote CRAM, range requests — no bulk download):
samtools view -b --reference scratch/g4/ref/giabv3/GRCh38_GIABv3.fa -L chr1.bed \
  https://ftp-trace.ncbi.nlm.nih.gov/.../BCM_ILMN-HG008Tp2_20250521/BCM_HTB3066_MGH3066_p2_MGH3066_p2_1.cram \
  -o p2T.chr1.bam   # same for p13; ~20 MB each
# SV eval (needs truvari>=5 in an isolated env, NOT g4-tools):
truvari bench -b <V0.6_PASS.vcf.gz> -c query_sv.vcf.gz --reference <GRCh38_GIABv3.fa> \
  --includebed <V0.6_clonal_and_subclonal.bed> --sizemax -1 --passonly --pick multi -o out
```

## 5. Evidence receipts

All runs `/usr/bin/time -v`, flags `-f 0.01 -c 1 -S 2 -E 3 -g 4 -z`, `-th 8`
unless noted. Timing receipts (`*.time`) and these summaries live in
`scratch/g4/hg008/ladder/`; committed `benchmarks/evidence/hg008-truth.json`
covers V0.3 only — the V0.4/V0.6 numbers below are supplementary in-guide results.

Estimator fit (`extra.json`):
`sec_per_mbp_th8=9.57, speedup_8=3.79, rss_full_kb=1496928,
rss_ratio_full_vs_chr1=1.07, full_th4 51.90s vs tier2 48.68s (+6.6%)`

V0.4 rtg (`None` rows): `all`: TP-base 613 / TP-call 570 / FP 903 / FN 268,
P=0.387 R=0.696 F1=0.497. `nogermline`: TP-base 551 / TP-call 520 / FP 826 /
FN 60, P=0.386 R=0.902 F1=0.541.

Tumor-only chr1 (`None` rows, 109 truth): unfiltered TP-base 101 FP 1117
P=0.081 R=0.927 F1=0.148; nognomad TP-base 69 FP 218 P=0.232 R=0.633 F1=0.340;
af0.001 F1=0.356; af0.01 TP-base 88 FP 286 P=0.229 R=0.807 F1=0.357;
af0.05 F1=0.356.

Tumor-only full-slice (`None` rows, 881 truth): any TP-base 478 FP 1509
P=0.218 R=0.543 F1=0.311; af0.001 F1=0.314; af0.01 TP-base 598 FP 2158
P=0.200 R=0.679 F1=0.310; af0.05 F1=0.304.

V0.6 truvari (`summary.json`): TP-base 0 / TP-comp 0 / FP 0 / FN 210
(`comp cnt: 0` — sizeless query SVs discarded).

Caveats: V0.4 is draft; single-sample VCFs carry no `STATUS` (use PASS-only
pre-filter, not `concordance.sh` as-is); chr1 tumor-only numbers are
small-domain (109 truth); precision figures are slice-domain only.
