# GCP Study Pre-registration Freeze (branch `paper/gcp-study`)

**Status:** FROZEN sections A = spend-gating. Nothing in A changes without owner sign-off + new freeze entry. Section B = open items, each one line + recommendation + why.
**Frozen:** 2026-10-16 (updated 2026-10-25 for `v1.0.4` / `552856835a1eec01`) · **Owner sign-off:** B1–B3, B5–B6 FROZEN; B4 gated on Task 5b `Gate 7 ref-multiregion` · **Spend to date:** ≈ $0.71 one-time + ≈ $1.45/mo (GCS vendor + refs + Task 5a Tier-P BAMs & RNA FASTQs, `task-5a-report.md §5`) — « $5 task cap, « $100 checkpoint.
**Requirements parent:** `docs/validation/gcp-study-protocol.md` (§2 H1–H4 + stop rules, §4 tiers + utility table, §5 Tier-3 + field, §6 prereg gates 1–8, §10-Wk1, Appendix A vendoring bullet).

---

## A. FROZEN (verified — artifact + line per value)

### A1. GCP identity + region + budget stops

| Item | Value (verbatim) | Source |
|---|---|---|
| Project | `vardict-paper-trial`, project #`694062694219` | `gcp-study-protocol.md:4` |
| SA | `vardict-agent@vardict-paper-trial.iam.gserviceaccount.com` (key 600, outside repo, never `git add`) | `gcp-study-protocol.md:4` |
| Billing | `011E32-BEB9F6-FE69ED`, trial $300/90d LINKED | `gcp-study-protocol.md:4` |
| Budget | `36c510ac-b312-4306-af98-a1bbde77244d`, alerts $100/$200/$300 = HARD STOPS (halt + re-evaluate + explicit resume) | `gcp-study-protocol.md:4`, `:44` |
| Region | `us-central1` (4 buckets: refs/bams/results/evidence, versioned, uniform access; ≤2 versions, Nearline after 30d) | `gcp-study-protocol.md:43-44` |
| Trial cap | 64 concurrent vCPUs account-wide (lifted from 12 on 2026-10-17, `eedf8d3`; UPGRADED full account 2026-10-16; 16/32-thread matrix UNBLOCKED, staged after 8-thread green; one 8-vCPU VM at a time remains courtesy norm) | `gcp-study-protocol.md:4`, `:44` |
| Spend order | P → E-span → G-span → HG008-full last (a cap halt still leaves publishable tiers) | `gcp-study-protocol.md:44` |
| Kill-vs-drain | Batch `timeout` per cell; killed = `SKIP (budget-stop)`, spend logged, never partial VCFs; UB-8 rerun attestation syncs `KEEP_DIR` to GCS on every failure path, counts inside the retry budget | `gcp-study-protocol.md:44` |

### A2. Hypotheses H1–H4 + stop rules (protocol §2, `gcp-study-protocol.md:30-36`)

- **H1 (parity, scoped):** Tier-2 `compare.sh` (eps 1e-6, records-aligned CHROM/POS/REF/ALT) java==rust 0 diffs on ALL enrolled cohorts; Tier-3 `rtg vcfeval --squash-ploidy` 0 site-diffs ONLY on truth-eligible cohorts (HG002/SEQC2/HG008/HG005/b37). Amplicon/RNA/SV-only/COLO829-structural = Tier-2-only, labeled ineligible — except HG008-V0.6-SV/CNV via truvari (gated, §A9). Gate stays java==rust; cpp never gates (`CPP_SKIP=1` escapes).
- **H1a (tile==single):** tiled output == single-region output (`--tile-size` sweep vs `10,000,000` / oracle untiled, boundary-dedup proof) MATCH before any H2 scale claim. Simple-mode and Splice (`-i`) first (`splice-tile-ownership` + `mate-window-ownership` landed in `v1.0.3`, `PR #29`); **somatic multi-tile `H1a` remains gated on `audit-somatic-dedup` (`ROADMAP.md:85`) landing with its own `compare.sh` MATCH** (somatic tiles currently emit without cross-tile `pending` deduplication).
- **H2 (memory):** Rust peak RSS ≤ Java per workload AND flat span scaling (RSS_249Mbp/RSS_10Mbp ≤ 1.5×). Java/cpp OOMs (incl. full-genome 3-way) scored `OOM-vs-completed` (Rust wins iff Rust completes with MATCH), never a ratio. H1 MATCH applies where oracle completes.
- **H3 (field):** pre-registered operating points (§A6) + full QUAL sweep (§A6) + MRD `-f 0` arm (Tier-P SRR + Tier-G chr1; Table-3 blowup is region-size-driven — equal at both `-f` — so `-f 0` is a correctness/operating-point cell, not a second scale cell). Rust EQUALs Java F1; field callers reported without rank-picking; MNP/SV strata separate.
- **H4 (maintenance):** archive proof (`gh api` archived flag + timestamp + pin, not screenshots) + static-musl/SBOM/SLSA/CVE/determinism/tabix matrix = migration case.
- **Stop rules:** any java≠rust diff → BLOCK, single-function port loop, no relaxation; UB-8 flake = ONE attested rerun on quiesced host preserving original `oracle-fail-*.log` (no retry loops); systematic OOM → highmem cell + report-both, never heap-retune silently.
- **Stats honesty:** timed n=3 = technical runs (same image/inputs/region; NOT biological replicates) → mean±range/min-max, no 95% CI; truth n=1 → point + site-diff=0; CIs only n≥5/bootstrap. Independent unit = cohort/sample/interval (variant counts are observations, never n). QUAL sweep reports ALL thresholds (multiplicity by full disclosure); somatic filter rule PASS-only + AF/Fisher frozen.

### A3. Vendored paper artifacts (PDF sha + Zenodo md5s + GCS paths)

| Object | md5 / sha256 | Bytes | GCS path | Source |
|---|---|---|---|---|
| `README.md` | md5 `8bb795e30739754a81ec48f45e542b53` (Zenodo-API-declared; task-text value `8bb7954e…` was a transcription typo — local matches API, task-1-report §7-C1) | 4,370 | `gs://vardict-paper-trial-evidence/vendor/README.md` | task-1-report §1; task-2-report §3 |
| `Figure_data.zip` | md5 `8b79769889e4424e2d0db9026343a896` | 13,637,860 | `gs://vardict-paper-trial-evidence/vendor/Figure_data.zip` | task-1-report §1; task-2-report §3 |
| `biorxiv.pdf` (v1, 9 pp, PDF-1.5, DOI 10.64898/2026.09.23.753771) | sha256 `d9505294ec26f02810ce67e0928c9597f8ba8bedcb81b40139bc77cf25e3dca2` (full, task-1-report §5); md5 `098e10657992d20108dec7913694cc58` (MANIFEST computed-locally basis from the Task-1 vendored file — Wk1 re-verify against GCS object before citation) | 396,827 | `gs://vardict-paper-trial-evidence/vendor/biorxiv.pdf` | task-1-report §5; Appx A |
| `Archive.zip` | md5 `7bbcd523e3c3c38b8ce56f926ddf30c8` (single generation `1792190839008539`, composite deleted) | 4,075,298,167 | `gs://vardict-paper-trial-refs/vardictcpp/Archive.zip` | task-2-report §1, §3 |
| Zenodo record | concept DOI 10.5281/zenodo.21836658, version DOI 10.5281/zenodo.22955424 (2026-09-25, CC-0); v1 DOI 10.5281/zenodo.21836659 cited alongside latest | — | — | `gcp-study-protocol.md:127`,`:135` |
| MANIFEST biorxiv URL note | versionless DOI `https://doi.org/10.64898/2026.09.23.753771` resolves to latest (v1 suffix dropped); no MANIFEST change needed | — | — | fix-round-1 |

### A4. Tier-P cohort pins (SRR + BioProjects + somatic pair + build)

- Accessions: **SRR15006375, SRR15006376, SRR15006386, SRR15006540** (adult AML bone-marrow, FLT3-ITD, FiLT3r, Illumina MiSeq) + **SRR8657348** (MV4-11 AML cell line, CCLE) — all NATIVELY hg19 (no liftover). BioProjects **PRJNA742684** (panel) / **PRJNA523380** (MV4-11). Somatic pair = **SRR15006376 tumor + SRR15006386 normal**. (task-1-report §4; `gcp-study-protocol.md:70`)
- Table-1 parity region/row counts (PDF-verbatim; Table 1 = 5-sample parity, Table 2 = single-sample runtime/memory replicates, verified in `task-1-report.md §4` and `task-5a-report.md §3`): 15006386: 11,732/4,014; 15006375: 12,292/4,056; 15006376: 105,949/12,739; 15006540: 135,145/16,031; SRR8657348: 773,667/54,448 — all 100% byte-identical. SRR8657348 single-sample grid role = Fig-1 perf arm (Tier-2 parity-only, NOT Fig-2). (task-1-report §4; `gcp-study-protocol.md:16`)

### A5. Archive BAM list + mapping provenance (bwa 0.7.19-r1273, RG)

- `bams/SRR8657348.bam` (949,697,894 B), `SRR15006375.bam` (515,231,925), `SRR15006376.bam` (353,645,177), `SRR15006386.bam` (498,686,670), `SRR15006540.bam` (499,864,311), `hg002.exome.20.bam` (244,592,675) — each + `.bai`. (task-2-report §2b; Archive `MD5SUMS.txt` 209 lines reused for input verification at rebuild)
- Five WES BAMs mapped in-house: `bwa mem -t 16 -R "@RG\tID:<acc>\tSM:<acc>\tPL:ILLUMINA" hg19 <acc>_1.fq.gz <acc>_2.fq.gz | samtools sort`, **bwa 0.7.19-r1273**, ref hg19/GRCh37 UCSC `chr`-naming (`hg19.fa` md5 `530d89d3ef07fdb2a9b3c701fb4ca486`). (extract `PROVENANCE.md:24-27`; extract `README.md:37-41`; extract `scripts/01_download_and_map.sh:40-42`)
- FASTQ source: ENA `fastq_ftp` (paired `_1`/`_2.fastq.gz`) via filereport API (`01_download_and_map.sh` `dl()`); exact per-accession ENA URLs + md5s in `harness/assets.tsv` (group `tierp`, ENA API 2026-10-16).

### A6. hg002 b37 chr20-only provenance (CPP-7/8 surface — keep b37-vs-hg19 split explicit)

- NIST GIAB pre-aligned download (NOT mapped by us): `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam` — Agilent SureSelect v5, **bwa mem 0.7.12 → `human_g1k_v37_decoy` (b37, no `chr` prefix), Picard MarkDuplicates 1.124**; chr`20` streamed via `samtools view` → `hg002.exome.20.bam` (3,632,411 reads); scoring ref `ref20.fa` (chr20 renamed `20`); truth `truth20.norm.vcf.gz` + `conf20.bed`. (task-2-report §2c; extract `README.md:60-67`)

### A7. BED inventory (11 paths) + mkbed-regeneration requirement

- `beds/{simple,panel17,sv,lowcov,amplicon,c17,c17amp,big,bigamp}.bed` (line counts 2000/1500/2/800/800/6000/6000/8000/8000; all `chrN start end R`) + `truth/conf20.bed` (10,140 lines) + `truth/target20.bed` (4,955 lines, b37 `20`-naming, no 4th column). (task-2-report §2a)
- Covered-target BEDs are DERIVED via `mkbed()` (mapped, primary, MAPQ≥1 spans merged within 100 bp — extract `scripts/01_download_and_map.sh:47-57`), NOT capture-kit BEDs: regenerate per-sample covered BEDs in prereg to reproduce Table-1 region counts; do not assume capture-kit BEDs. (`gcp-study-protocol.md:88`; task-2-report §4-4)

### A8. Fig-2 truth + `eval_caller.sh` invocation (verbatim)

- Fig-2 = **GIAB HG002 chr20 exome** (Agilent SureSelect v5) vs **NIST GIAB v4.2.1** over high-confidence exome-callable regions (**depth≥20, 2.31 Mbp, 2,443 truth variants**), scored with **`rtg vcfeval`**, regions/method per Daniels et al. 2024 [16]; F1 denominator = the 2,443 v4.2.1 chr20 HC truth variants. Caller set: octopus / vardictcpp (=Java, identical) / FreeBayes / LoFreq / Mutect2-tumor-only. Vendored F1s (re-measure at Wk1): SNV — octopus 0.9925, vardictcpp=Java 0.9854, FreeBayes 0.9822, LoFreq 0.9771, Mutect2 0.6005; Indel — 0.9071 octopus / 0.8277 FreeBayes / 0.8009 vardict / 0.7761 LoFreq / 0.4874 Mutect2-tumor-only. MV4-11 is the Fig-1 arm, NOT Fig-2. (task-1-report §2; `gcp-study-protocol.md:81`,`:96`)
- Invocation verbatim (extract `scripts/eval_caller.sh:1-30`): `bash eval_caller.sh <label> <caller.vcf[.gz]>` → normalise (`bcftools reheader --fai $D/ref20.fa.fai | bcftools norm -f $D/ref20.fa -m-any | bcftools view -f PASS,. | bcftools sort`, bgzip+tabix), then per `snp`/`indel`: `rtg vcfeval -b /tmp/tr_$cls.vcf.gz -c /tmp/ev_${label}_$cls.vcf.gz -t $D/ref20.sdf --evaluation-regions $D/target20.bed -o /tmp/ev_${label}_$cls --vcf-score-field=QUAL`, best-F row from `summary.txt` (`P=$6 R=$7 F=$8`); emits `<label> calls=<n> SNV P/R/F1 … INDEL P/R/F1 …`.

### A9. Somatic / germline filter rules (verbatim from harness)

- Somatic (`harness/concordance.sh:7-8,12-13`): `bcftools view -f PASS -i 'INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"' "$CALLS" | bcftools sort -Oz …`; truth `bcftools annotate -x FILTER` (stripped, else rtg drops it); `rtg vcfeval -b truth.nofilter.vcf.gz -c calls.som.vcf.gz -t $SDF -e $BED --squash-ploidy --sample "ALT,$TS" --vcf-score-field QUAL -T 8`; gate on the `None` row (tp_baseline/tp_call/fp/fn/P/R/F1); fail on empty baseline.
- Germline (`harness/concordance_germline.sh:31-37`): `bcftools view -f PASS "$VCF" | bcftools sort -Oz …`; truth FILTER-stripped the same way; `rtg vcfeval -b $TRUTH_NOFILTER -c calls.pass.vcf.gz -t $SDF -e $BED --sample "$SAMPLE,$SAMPLE" --vcf-score-field QUAL -T 8`; gate = identical `None`-row counts AND identical tp/tp-baseline/fp/fn site lists java==rust; cpp third arm advisory only (`concordance_germline.sh:128-129`: `vardictcpp 3 hc52dbad_0 bioconda`), rank line advisory, never gated; `CPP_SKIP=1` escapes.

### A10. Operating points + sweep (frozen)

- Germline operating points: **QUAL ≥ 128 / 226 / 160** (frozen set; per-cohort assignment asserted in wrapper `extra`). Somatic operating point: **PASS + STATUS rule (§A9), no cut beyond PASS+STATUS** — no AF/Fisher numeric thresholds exist in `concordance.sh` (PASS + `StrongSomatic`/`LikelySomatic` only) or the protocol, so none are pinned. (`gcp-study-protocol.md:96`)
- Full ROC sweep table, title **"ROC/QUAL sweep"**: QUAL **0/30/50/75/100/128/150/200 + VD equivalents**, ALL thresholds reported (multiplicity by full disclosure); transition reported as result not premise. (`gcp-study-protocol.md:96`)
- MRD `-f 0` arm: Tier-P SRR + Tier-G chr1 (correctness/operating-point cell, not a second scale cell). (`gcp-study-protocol.md:33`)
- Truths: GIAB **v4.2.1** (HG002/HG005/b37), SEQC2 **v1.2.1**, HG008 **v0.3** (continuity) + **v0.4 smvar** in parallel via rtg (`assets.tsv:18-21`); neither v0.3/v0.4 is SV+CNV — that is V0.6. (`gcp-study-protocol.md:96`)
- Stratify SNP/InDel/MNP/SV via named tool (`bcftools` strata + masks + counting rule — `rtg` alone insufficient) + AF bins (≤1%/2-read sensitivity) + depth bins. (`gcp-study-protocol.md:96`)
- cpp zero-`##contig` truth-header completion (CPP-6): truth-first, `##fileformat` line 1, zero record-byte changes, asserted in wrapper. (`gcp-study-protocol.md:96`)

### A11. Tier × flag expected MATCH/MISMATCH/SKIP matrix (one row per tier×engine surface)

Gate stays java==rust everywhere; cpp advisory, binned CPP-1..9, never gates. (`gcp-study-protocol.md:30`,`:88`,`:95`)

| # | Tier × surface | java | rust | cpp (advisory) | Tier-3 | Pinned keys in `extra` |
|---|---|---|---|---|---|---|
| P1 | P amplicon EGFR/MRE11a full (8-col BED, primer-trim + AMPBIAS) | MATCH | MATCH | SKIP CPP-3 (TSV-only) | ineligible (assay) | budget + ceiling cells |
| P2 | P SRR panel 15006375/76/86/540 + SRR8657348 (hg19, covered BEDs) | MATCH (Table-1 counts, FP=FN=0) | MATCH | advisory (archive `parity.tsv` 100% on 5 samples); MISMATCH/SKIP binned CPP-1..10 | ineligible (SRR8657348 = Fig-1 perf arm, parity-only) | `--chunk=$CPP_CHUNK`, BED-vs-`-R` per engine, `interval_mechanism_{java,rust,cpp}` (CPP-5; cpp `-R`+BED never co-passed) |
| P3 | P somatic pair 76/86 | MATCH | MATCH | advisory | Tier-3 via SEQC2/HG008 only (SRR pair = Fig-3 Venn characterization, no truth rank) | same as P2 |
| P4 | P panel SV calls | MATCH | MATCH | advisory, CPP-9 drop rates (1–2 SV alleles) reported | Tier-2 MATCH required | same as P2 |
| P5 | P panel CRAM twin | N/A (Java 1.8.3 cannot read CRAM: `NoClassDefFoundError` BZip2 — extract `README.md:99-102`) | MATCH vs own BAM run (byte-identical) | MATCH vs own BAM run | — | bitwise-identity cell, separate from P2 |
| E1 | E SEQC2 HCC1395/BL full exome + exome-wide SVs | MATCH (623k) | MATCH | MISMATCH CPP-9 (1–2 SV: `<DEL>` + 495 bp DEL at `chr22:36511000`) | EQUAL SEQC2 v1.2.1 | budget + jobs/node + 16/32 staged |
| E2 | E exome CRAM twin | N/A (as P5) | MATCH vs own BAM | MATCH vs own BAM | — | as P5 |
| G1 | G HG002 full WGS germline | MATCH | MATCH | MISMATCH CPP-9 (346/227,562 key-level loci [15 MNP + 331 tie-breaks] + 2 dropped DELs [287 bp, 747 bp] = 0.15%) | EQUAL GIAB v4.2.1 (F1 0.9699) | budget + span + 16/32 |
| G2 | G HG008 T/N full WGS | MATCH | MATCH (somatic multi-tile `H1a` gated on `audit-somatic-dedup` + own `compare.sh` MATCH) | MISMATCH CPP-9 (SV/long-DEL drop: 744 bp DEL) | EQUAL v0.3 + v0.4 smvar parallel rtg; V0.6 SV+CNV via truvari GATED (BLOCKED, §A9/§A15) | budget + ceiling |
| G3 | G span curve chr1 1/10/100/249 Mbp + `-R` repro + ceiling cells | `-R` cell: `-R chr1:1-249250621` SRR8657348 hg19, JDK25 `-Xmx 200g` 107.6 GB/6:10 vs OOM at `-Xmx 4g` (reference arm) | flat RSS (ratio ≤1.5×; `v1.0.3` 1.223× at 1.49 GB / 1.42 GiB); ceiling ≤8 GB else BLOCK full-genome (§A12) | `-R` contained repro 10.6 GB/2:03 (reference arm) | — | cgroup/`ulimit`/`timeout`, `KEEP_DIR` preserves `cpp_R*.time` + exit code, no `rm -f` |
| G4 | G CRAM twins | N/A (as P5) | MATCH vs own BAM | MATCH vs own BAM | — | as P5 |
| R1 | RNA HG002 GM24385 rep-1 mRNA full transcriptome (NovaSeq 6000 S4 2×150, ~325M clusters; matched LCL to Tier-G DNA) | MATCH sorted-set | MATCH sorted-set | SKIP (`-i` unsupported) | ineligible (parity-only) | STAR build+align wall/RSS on `n2-highmem-8` (reuse staged `HG002_rep1_Aligned.sortedByCoord.out.bam`); `-i` order nondeterministic by design (plain-HashMap iteration, `DIVERGENCES.md`) → sorted line-SET identity (`--emit-intermediate`), never byte order |
| S1 | SV COLO829(-BL) full WGS (`<INV>`/`<DUP>`/CRISPR `-J`/`-j`) | MATCH | MATCH | CPP-2 (`-J`/`-j`) | ineligible (Valle-Inclán 2022 = parity); HG008 v0.4 gives SVs a truth leg | SV strata: Fig5 MNP = HG002 germline; SV = HG008-V0.6 Tier-3 where unblocked, else COLO829 Tier-2 (no double-count) |
| D1 | Diversity/b37 strata HG005 + NA24385-b37 full (tier of origin) | MATCH | MATCH | CPP-7/8 surfaces | EQUAL v4.2.1 | b37-vs-hg19 split explicit |
| T1 | 16/32-thread matrix | report-only | gates S_16>S_8 + RSS flat | report-only, contained | — | Tier-G span + Tier-G genome Tier-2 all 3 engines (Rust `n2-standard-16/32`; Java/cpp `n2-highmem-16/32`); Tier-E @16/32 only after Tier-G clean; Tier-P excluded; UNBLOCKED (64-vCPU quota confirmed 2026-10-17) |
| T2 | Determinism 1/2/4/8 threads (≥1 real somatic + ≥1 SV/MNP cohort) | records-exact + named header policy (NOT full-byte: harness strips `##`; RNA `-i` sorted-set only) | same | same harness | — | `gcp-study-protocol.md:97` |
| C1 | CLI concat-flags surface | — | — | CPP-4 (concat flags exit 120) expected | — | 60-flag CLI matrix (thread-default difference VERIFIED: Java `CmdParser.java:617-628` → 1 when `-th` absent; cpp §2 default 1; ours all-cores — always pass `-th` explicitly, record defaults) |

### A12. Compute gates 1–8 (all pinned before spend — `gcp-study-protocol.md:103`)

1. **Ceiling FIRST** — Rust peak RSS **≤8 GB on both ceiling cells** (1000× amplicon + HG008 WGS slice at `-th 8`) else **BLOCK** full-genome (span `1.223×` at `1.49 GB` / `1.42 GiB` in `v1.0.3` was exome-depth, does not transfer to 195x/169x unbounded vecs).
2. **Java sizing & Single-VM Cgroup Architecture** — All GCP cells (`Tasks 5b–5e`) execute on a single dedicated **`n2-highmem-32` (`32 vCPUs`, `256 GB RAM`, `2 TB pd-ssd` [`--boot-disk-size=2048GB --boot-disk-type=pd-ssd`], `us-central1-a`)** via `harness/gcp_paper_runner.sh`. Each cell is isolated via `systemd-run --scope -p AllowedCPUs=0-$((P-1)) -p MemoryMax=<CAP>` with `-Xmx` = **0.7× slice RAM** (`ORACLE_XMX_GB`) in `run_oracle.sh` (`200G` on `MemoryMax=240G` for the Table-3 `-R chr1` Java completion cell; `4G` for the Table-3 OOM cell); **137-path split**: exit-137-under-cgroup-cap = OOM path (score, no UB-8 rerun); 137-with-headroom = UB-8 path (one rerun, log preserved).
3. **cpp containment** — `MemoryMax=32G` cgroup slice on `n2-highmem-32`, `KEEP_DIR` synced to GCS on every exit path, expected `CPP: SKIP (cgroup OOM)` counts as containment PASS.
4. **RNA alignment on `n2-highmem-32` (`2 TB pd-ssd`)** — `STAR 2.7.11b --runThreadN 32` build+align on the `2 TB pd-ssd` workspace, FASTQ pull + index + aligned BAM (`gs://vardict-paper-trial-bams/rna-hg002/aligned/HG002_rep1_Aligned.sortedByCoord.out.bam`) synced to GCS with **≤$10 sub-cap**, index RSS in `extra`.
5. **SRR BAM provenance per accession (CLOSED in Task 5a, `task-5a-report.md §2`)** — All 6 pre-aligned BAMs + `.bai` extracted from `gs://vardict-paper-trial-refs/vardictcpp/Archive.zip`, MD5-verified against `MD5SUMS.txt` (`SRR15006375.bam 23cb57e7…`, `SRR15006376.bam c735d817…`, `SRR15006386.bam 98625e7f…`, `SRR15006540.bam 69828a79…`, `SRR8657348.bam f4b66ac8…`, `hg002.exome.20.bam fa809088…`), and staged to `gs://vardict-paper-trial-bams/tier-p/` (`12/12` objects server-MD5 verified; zero FASTQ realignment required).
6. **jobs/node method** — 2×4-thread jobs inside an `AllowedCPUs=0-7 MemoryMax=32G` (`n2-standard-8`-equivalent) cgroup slice on `n2-highmem-32` (no second VM needed), else defer with scope cut.
7. **ref-multiregion gate** — 2+-region oracle MATCH proving per-buffer offsets (`reference.rs:301-318`), else full-genome restricted to single-region streaming (logged).
8. **2-Strike Spot → On-Demand Escalation Policy (locked 2026-10-24)** — Launch `n2-highmem-32` (`32 vCPUs`, `256 GB RAM`, `250 GB pd-ssd` boot + `8×375 GB` [`3.0 TB`] NVMe Local SSD RAID0) as `SPOT` (`--provisioning-model=SPOT`) with per-cell GCS checkpointing; if `SPOT` is preempted/killed twice (`2` preemptions), escalate to `STANDARD` (On-Demand) on attempt 3 and resume from the last completed GCS cell checkpoint.

Plus: `export GOLDEN_MODE=off` on all GCP timed/OOM cells; same BAM+ref+BED+flags+threads per cell; pre-flight ≥2× 8-thread RSS free or SKIP loud; metrics wall/peak RSS/ratio/S_8/$/sample/jobs/node + `tabix -p vcf`/`bcftools index` pass per cohort + loud SKIP per field caller (java==rust untouched); publish unfavorable cells too (e.g., SEQC2 exome `431.99s` C++ vs `466.34s` Rust at `4.7×` higher RAM, HG008 slice `311 MB` C++ vs `321 MB` Rust, and `scale.json` `1 Mbp` small-window startup floor vs `1.49 GB` / `1.42 GiB` Rust at `249 Mbp`). Amplicon memory: non-`-p` VCF is variant-bounded (`Task 4b`) and `-p --emit-intermediate` is streamed (`Task 4d`, `ROADMAP.md:89`); `-p` VCF stays out of scope (`gcp-study-protocol.md:125`). `-R` full-chromosome relevance: Galaxy wrapper exposes region for full chromosomes only (tools-iuc `macros.xml` L27) — mechanism, not just our finding. (`gcp-study-protocol.md:88`,`:103`; `:125` for open `audit-somatic-dedup`/`audit-stream-pileup`/`audit-fisher-stream`/`audit-ref-multiregion` gates in `ROADMAP.md:85-88`)

### A13. Field callers — versions + exact commands (Task-2-verified, all bioconda)

| Caller | Version | Build | Env | Source |
|---|---|---|---|---|
| FreeBayes | **1.3.10** | `hbefcdb2_0` | `cl_freebayes` | extract `PROVENANCE.md:14` |
| LoFreq | **2.1.5** | `py313ha7493ab_16` | `cl_lofreq` | extract `PROVENANCE.md:15` |
| Mutect2 (GATK4) | **4.6.2.0** (run under Java 21) | `py310hdfd78af_1` | `cl_gatk4` | extract `PROVENANCE.md:16` |
| octopus | **0.7.4** | `ha3c1580_2` | `cl_octopus` | extract `PROVENANCE.md:17` |
| Strelka2 | **2.9.10** | `hdfd78af_2` | `g4-tools` / bioconda (ABSENT in vardictcpp archive; our extension arm, `B2` frozen) | `gcp-study-protocol.md:76` |
| vardictcpp | **archive commit `ae93fe39dc`** (`fix/somatic-amplicon-parity`, verified in Task 5a, `B6` frozen) · **local bioconda `3=hc52dbad_0`** | `ae93fe39dc` / `hc52dbad_0` | `paper-vardictcpp-build` / `vardict-oracle` | `task-5a-report.md §3`; `benchmarks/env/vardict-oracle.yml:97` |
| VarDictRust | **1.0.4** (`2503907`, engine fingerprint **`552856835a1eec01`**, `--perf` fingerprint **`9346a6e5419cc5e2`**) | `rustc 1.96.0 (ac68faa20 2026-05-25)`, `noodles 0.117.0`, `rayon 1.12.0`, `tikv-jemallocator` | release tag `v1.0.4` | `benchmarks/env/toolchain.txt:1`; `Cargo.lock` |
| VarDictJava | 1.8.3 (`hdfd78af_0`: `VarDict-1.8.3.jar` + commons-cli 1.2, commons-math3 3.6.1, htsjdk 2.21.1, jregex 1.2_01) | — | `vardict` / `vardict-oracle` | extract `README.md:9-10` |
| STAR | **2.7.11b** | `h5ca1c30_8` | bioconda (`B1` frozen; with GENCODE v50 basic-CHR GTF) | `gcp-study-protocol.md:76` |
| rtg-tools | **3.13** | — | `g4-tools` | `benchmarks/env/g4-tools.yml:81` |
| bcftools (stratifier) | **1.24** | — | `g4-tools` | `benchmarks/env/g4-tools.yml:10` |
| bwa | 0.7.19-r1273 | — | `g4-tools` | `benchmarks/env/g4-tools.yml` (`bwa=0.7.19`) |

Exact commands (verbatim; `$REF`=hg19, `$BED`=per-sample covered BED):

- Single-sample — FreeBayes 1c: `freebayes -f $REF -t $BED $BAM`; 8c: `freebayes-parallel fb_regions.txt 8 -f $REF $BAM` (`fb_regions.txt` = `awk '{print $1":"$2"-"$3}' $BED`). (`single_sample_speed.sh:47-48`)
- Single-sample — LoFreq: `lofreq indelqual --dindel -f $REF -o srr.iq.bam $BAM` (once) then `lofreq call --call-indels -f $REF -l $BED -o lf1.vcf $IQ` / `lofreq call-parallel --pp-threads 8 … -o lf8.vcf $IQ`. (`single_sample_speed.sh:50-51`; `PROVENANCE.md:39-41`)
- Single-sample — Mutect2 (tumor-only): `gatk Mutect2 -R $REF -I $BAM -L $BED --native-pair-hmm-threads {1|8} -O m2_{1,8}.vcf.gz`. (`single_sample_speed.sh:52-53`)
- Single-sample — octopus: `octopus -R $REF -I $BAM -t $BED -o oc{1,8}.vcf --threads {1|8}`. (`single_sample_speed.sh:54-55`)
- Single-sample — vardictcpp/Java reference: `vardictcpp -G $REF -f 0.01 -N s -b $BAM -c 1 -S 2 -E 3 -g 4 $BED -th {1|8}`; Java: `java -Xmx8g -cp $JCP com.astrazeneca.vardict.Main` + same args. (`single_sample_speed.sh:43-46`; `PROVENANCE.md:47`)
- Somatic 76/86 — vardictcpp/Java: `-N 'SRR15006376|SRR15006386' -b '$TUM|$NOR'` + same flags. (`somatic_speed.sh:39-42`; `PROVENANCE.md:48`)
- Somatic — Mutect2: `gatk Mutect2 -R $REF -I $TUM -I $NOR -normal $NN -L $BED --native-pair-hmm-threads {1|8}`. (`somatic_speed.sh:43-44`)
- Somatic — octopus: `octopus -R $REF -I $TUM $NOR --normal-sample $NN -t $BED --threads {1|8} -o …`. (`somatic_speed.sh:45-46`)
- Somatic — LoFreq: `lofreq somatic -n $NOR -t $TUM -f $REF -l $BED -o lf{1,8}_ --threads {1|8}`. (`somatic_speed.sh:47-48`)
- HG002 arm — FreeBayes/LoFreq/octopus as above with `$REF=ref20.fa $BED=target20.bed`; Mutect2 adds `FilterMutectCalls -R $REF -V m2.vcf.gz -O m2.filt.vcf.gz`; LoFreq needs `indelqual` BAM + `samtools index`. (`hg002_callers.sh:11-19`)
- Fig-2 germline field = octopus + FreeBayes; Fig-3 somatic field = Mutect2 (PASS) + octopus (somatic) + LoFreq (somatic); +Strelka2 extension labeled; pinned digests + PASS-vs-all filter + shared BED/SDF + `squash-ploidy` justification per cell. (`gcp-study-protocol.md:96`)

### A14. Fig-1 / Fig-3 replication targets (vendored `Figure_data.zip`, 46 files + 7 dirs)

- `caller_speed.tsv` (Fig-1 single-sample, SRR8657348 full 773,667-region exome): vardictcpp 1c 101.74 s/0.12 GB, 8c 13.94/0.25; Java 479.36/0.96, 139.68/2.02; FreeBayes 668.78/0.86, 77.16/4.62; LoFreq 747.22/2.95, 81.75/0.26; Mutect2 603.27/1.05, 408.39/1.64; octopus 891.94/0.84, 785.14/1.51. (extract `results/caller_speed.tsv`)
- `caller_speed_somatic.tsv` (Fig-3 somatic 76/86): vardictcpp 184.1/0.15, 34.8/0.75; Java 579.7/1.87, 302.6/2.83; Mutect2 1610.8/2.44, 679.7/2.53; octopus 1307.2/0.95, 401/2.15; LoFreq 1845.8/2.97, 371.6/1.93 (calls: vardict 14126, Mutect2 478, LoFreq 111). (extract `performance_benchmark/results/somatic_speed.tsv`)
- `caller_accuracy.tsv` + `caller_accuracy_rtg.out` (§A8 values). (extract `results/caller_accuracy.tsv`)
- Venn (PDF-verbatim): VarDict **737** SNVs; sole-called **606**; overlaps **21/1/0/1/9/86/1**; all-four **13**; Fig-3 somatic Venn has NO orthogonal truth set (concordance characterization, not accuracy rank — task-1-report §7-C6). Sets vendored as `fig3_somatic_venn/{vardict,mutect2,octopus,lofreq}.snv` (vardict 737 / mutect2 125 / octopus 17 / lofreq 111 lines). Reconciled note: PDF-verbatim `606` is **VarDict-sole SNVs** ($606 + 118 + 13 = 737$ VarDict SNVs), whereas `concordance/somatic_snv_overlap.txt` `exactly-1: 610` is the **sum of sole-called SNVs across all 4 callers** ($610 + 32 + 88 + 13 = 743$ 4-caller union; both numbers are self-consistent).

### A15. HG008-V0.6 truvari gating

- Enroll HG008 V0.6 stvar/CNV (PASS/ALL VCFs + clonal and clonal+subclonal beds + CNV bedpe/bed, `assets.tsv:26-33`) via **truvari (NOT rtg)** for the SV arm — **GATED on V0.6 representation fix** (sizeless `<DEL>/<DUP>/<INV>` discarded, CNV method none; eval currently BLOCKED per guide). COLO829 stays Tier-2 caller-parity; HG008 v0.4 gives SVs a truth leg. (`gcp-study-protocol.md:96`)

---

## B. RESOLVED / GATED SIGN-OFF ITEMS

- **B1 STAR + GENCODE GTF (FROZEN 2026-10-17, `53de256`):** **STAR 2.7.11b=h5ca1c30_8** (bioconda) + **GENCODE Release 50 (GRCh38.p14) basic-CHR GTF**. Contig-prefix check: both `GRCh38_GIABv3.fasta` and GENCODE v50 basic-CHR GTF use `chr1..chr22,chrX,chrY,chrM` naming (no prefix translation needed).
- **B2 Strelka2 (FROZEN 2026-10-17, `53de256`):** bioconda **`strelka=2.9.10=hdfd78af_2`** (our somatic extension arm; absent in vardictcpp archive).
- **B3 Venue (FROZEN 2026-10-17, `53de256` / `839da96`):** **Genome Biology Research (PRIMARY)** + **Bioinformatics Original Article (FALLBACK #1)** preceded by **bioRxiv** preprint; author guides checked `2026-10-16` (local HTML copies).
- **B4 Full-genome invocation (GATED on Task 5b `Gate 7 ref-multiregion`):** per-chromosome regions over GRCh38 GIABv3 contigs streamed per-region (no whole-chromosome `-R` except the contained Table-3 repro cell), BED-vs-`-R` tabulated per cell per engine, `--chunk` = `$CPP_CHUNK` (`50000` default in `run_cpp.sh`) recorded in `extra`.
- **B5 16/32-thread (FROZEN / UNBLOCKED 2026-10-17, `eedf8d3`):** 64-vCPU account quota confirmed; 16/32-thread matrix proceeds staged after 8-thread green.
- **B6 Commit hash (FROZEN / RESOLVED 2026-10-17, `f8afc46` → `ae93fe39dc`):** Task 5a fresh builds (`cmake 4.4.4` / `htslib 1.24` / `g++ 11.4`) with recorded somatic invocation reproduce archived raw-stream md5 `6674ec0b…` at `-th 1` AND `-th 8` (`14,126` rows), whereas `5ab3759` emits `14,127` rows. Cite `ae93fe39dc`.
