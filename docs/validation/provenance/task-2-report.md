# Task 2 report — vardictcpp Archive.zip vendoring (paper/gcp-study)

**Status: DONE_WITH_CONCERNS** (fix round 1/5: re-verified — see Evidence) · date 2026-10-16 · no commit (new files under gitignored `docs/superpowers/` only) · no MANIFEST/assets.tsv touched.

**Evidence (fix round 1/5):** text-only extraction re-saved to
`docs/superpowers/specs/paper-plan/task-2-text-extract/` (112 files, 2.0 MB,
same include/exclude patterns as §2; `diff -r` vs `/tmp` original = IDENTICAL)
+ full listing `task-2-archive-listing.txt`. Every §2/§4 enumerated detail
re-checked against these saved files — zero disagreements, no content changes.
Behavior note: per parent-session instruction, no long hardcoded sleeps were used
after the midpoint — short bounded wait-loops / blocking commands only.

## 0. Auth + preconditions

- Env override only: `CLOUDSDK_CORE_ACCOUNT=vardict-agent@…`
  `CLOUDSDK_CORE_PROJECT=vardict-paper-trial` (no activate/revoke run).
- `gcloud storage buckets list` shows all 4 buckets:
  `vardict-paper-trial-{bams,evidence,refs,results}` (us-central1, versioned). Proceeded.

## 1. Download + verification (step 1)

- `mkdir -p /tmp/opencode/vardictcpp-archive` (2.0 TB free confirmed).
- `curl -C - -L --retry 5 --retry-delay 10 -o Archive.zip
  https://zenodo.org/api/records/22955424/files/Archive.zip/content` (single stream).
- Result: **4,075,298,167 B, md5 `7bbcd523e3c3c38b8ce56f926ddf30c8`** —
  both match plan values exactly. No resume needed (first try completed).

## 2. Enumeration (step 2, no full extraction)

Full `unzip -l` (223 files, exit 0) saved as `task-2-archive-listing.txt`
(same directory as this report). Text-only extraction (~2 MB) via
`unzip '*.bed' '*.sh' '*.txt' '*.tsv' '*.md' '*.out' '*.py' '*.rc' '*.wall'
'*.peak' '*.log' '*.stats'` with `-x '*.bam*' '*.vcf.gz*' '*.tbi' '*.png'
'*.pdf' '*.tsv.gz' '*.fa*'` — no BAM/VCF/fastq touched.

**(a) Exact BED paths inside the archive (11 total):**
`beds/{simple,panel17,sv,lowcov,amplicon,c17,c17amp,big,bigamp}.bed`
(line counts: simple 2000, panel17 1500, sv **2**, lowcov 800, amplicon 800,
c17 6000, c17amp 6000, big 8000, bigamp 8000; all `chrN start end R` except
truth BEDs) + `truth/conf20.bed` (10,140 lines) + `truth/target20.bed`
(4,955 lines, `20`-style b37 naming, no 4th column).

**(b) Full BAM list — YES, both SRR15006375 and SRR15006540 exist:**
`bams/SRR8657348.bam` (949,697,894 B), `bams/SRR15006375.bam` (515,231,925),
`bams/SRR15006376.bam` (353,645,177), `bams/SRR15006386.bam` (498,686,670),
`bams/SRR15006540.bam` (499,864,311), `bams/hg002.exome.20.bam` (244,592,675)
— each with `.bai`. (Sizes are compressed-BAM bytes, cf. README's "BAM size".)
Per-archive README §2.1: all five WES BAMs mapped in-house
(bwa mem **0.7.19-r1273** → hg19, RG `ID:<acc> SM:<acc> PL:ILLUMINA`).

**(c) hg002.exome.20 provenance — NO SRA accession; pre-aligned NIST GIAB
download** (archive `README.md` §2.2; nothing about it in `01/02_*.sh` —
those only fetch/map the 5 WES SRRs):
- URL: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/
  AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/
  151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam`
- Original processing (BAM `@PG`): Agilent SureSelect v5, bwa mem **0.7.12** →
  `human_g1k_v37_decoy` (**b37**, no `chr` prefix), Picard MarkDuplicates 1.124.
- Vendored subset: chr`20` streamed via `samtools view` → `hg002.exome.20.bam`
  (3,632,411 reads, 244,592,675 B; md5 `fa8090887f81cf3563ef3fd5db810067`
  per `MD5SUMS.txt`). Truth: NIST GIAB **v4.2.1** HG002 chr20
  (`truth/truth20.norm.vcf.gz`, `conf20.bed`); scoring ref `ref20.fa`
  (chr20 renamed `20`). **There is no `giab_accuracy.sh`** — equivalents are
  `scripts/hg002_callers.sh` (caller driver) + `scripts/eval_caller.sh`
  (rtg vcfeval, best-F, PASS-only).

**(d) Field-caller versions — there is NO `benchmark.sh`.** Benchmark drivers are
`performance_benchmark/scripts/{single_sample_speed,somatic_speed}.sh`
(commands below verbatim from those files); versions pinned in
`performance_benchmark/PROVENANCE.md` (all bioconda):
| caller | version | build | env |
|---|---|---|---|
| FreeBayes | **1.3.10** | `hbefcdb2_0` | `cl_freebayes` |
| LoFreq | **2.1.5** | `py313ha7493ab_16` | `cl_lofreq` |
| Mutect2 (GATK4) | **4.6.2.0** | `py310hdfd78af_1` | `cl_gatk4` (run under Java 21) |
| octopus | **0.7.4** | `ha3c1580_2` | `cl_octopus` |
| Strelka2 | **absent** | — | zero hits for `strelka` in listing + all text |
Exact commands: FreeBayes `freebayes -f $REF -t $BED $BAM` (1c) /
`freebayes-parallel fb_regions.txt 8 -f $REF $BAM` (8c);
LoFreq single `lofreq call --call-indels -f $REF -l $BED -o out $IQ`
(`$IQ` = dindel-recalibrated `srr.iq.bam` via `lofreq indelqual --dindel`) /
`lofreq call-parallel --pp-threads 8 …`; LoFreq somatic
`lofreq somatic -n $NOR -t $TUM -f $REF -l $BED -o out_ --threads {1|8}`;
Mutect2 `gatk Mutect2 -R $REF -I $BAM -L $BED --native-pair-hmm-threads {1|8}
-O out.vcf.gz` (tumor-only in single-sample; `-I $TUM -I $NOR -normal $NN`
in somatic); octopus `octopus -R $REF -I $BAM -t $BED -o out.vcf --threads
{1|8}` (somatic: `-I $TUM $NOR --normal-sample $NN`); HG002 arm adds
`FilterMutectCalls` for Mutect2 (`hg002_callers.sh`).

## 3. GCS staging (step 3) — md5 table (local vs `ls -L`)

| Object | Size (B) | md5 local | md5 GCS | verdict |
|---|---|---|---|---|
| `gs://…-evidence/vendor/README.md` | 4,370 | `8bb795e3…542b53` | `i7eV4wc5dUqB7Ej0XlQrUw==` | MATCH |
| `gs://…-evidence/vendor/Figure_data.zip` | 13,637,860 | `8b797698…6343a896` (= plan md5) | `i3l2mInkQk4tDbkCY0Oolg==` | MATCH |
| `gs://…-evidence/vendor/biorxiv.pdf` | 396,827 | sha256 `d9505294…25e3dca2` (= plan sha) | `CY4QZXmS0gEI3seRNpTMWA==` | MATCH |
| `gs://…-refs/vardictcpp/Archive.zip` | 4,075,298,167 | `7bbcd523…ddf30c8` | `e7zVI+PDw4uM5W+Sbd8wyA==` | MATCH |

First upload of Archive.zip landed as a 32-component composite object
(parallel composite upload default; no server MD5), so it was re-uploaded
single-stream (`CLOUDSDK_STORAGE_PARALLEL_COMPOSITE_UPLOAD_ENABLED=False`),
MD5-verified, and the stale composite generation (`…#1792189655762589`)
deleted — one live version remains (gen `1792190839008539`).
**Spend note:** bytes stored = 14,039,057 (vendor/) + 4,075,298,167 (refs/) =
**4,089,337,224 B (~3.81 GiB)** Standard us-central1 → first-month ≈
**$0.08 storage + <$0.01 ops ≈ $0.08–0.10**, egress $0 (matches plan estimate).

## 4. Concerns (for Wk1 owner / prereg)

1. **Filenames ≠ plan shorthand:** `benchmark.sh` / `giab_accuracy.sh` do not
   exist — cite `performance_benchmark/scripts/{single_sample_speed,
   somatic_speed}.sh` + `scripts/{hg002_callers,eval_caller}.sh` in prereg.
2. **No Strelka2 anywhere** in the archive; Strelka2 pins/digests stay open
   (already in Appendix C AUTHOR_INPUT_NEEDED).
3. **hg002 arm is b37 chr20-only, pre-aligned** — CPP-7/8 surface; keep the
   b37-vs-hg19 split explicit in prereg (plan already does).
4. **Covered-target BEDs are derived** (`mkbed()` in `01_download_and_map.sh`,
   not capture-kit BEDs): reproducing Table-2 region counts (11,732 / 12,292 /
   105,949 / 135,145) requires regenerating per-sample covered BEDs — consistent
   with the plan's Tier-P BED-vs-`-R` gate, but the regeneration step must be in
   prereg, not assumed.
5. **Commit-hash inconsistency in archive docs:** top-level `README.md` says
   vardictcpp commit `5ab3759`, `PROVENANCE.md` says `ae93fe39dc` (same branch
   `fix/somatic-amplicon-parity`). Wk1 owner to confirm which binary the raw
   outputs came from before citing.
6. `MD5SUMS.txt` (209 lines) covers every archive file — reuse for input
   verification at rebuild; top-level `README.md` §7 gives the
   `samtools view … ##idx##… 20` rebuild one-liner for hg002 chr20.
