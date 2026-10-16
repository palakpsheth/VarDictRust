# Task 1 Report — Zenodo vardictcpp vendoring (local-network only, zero GCP spend)

**Date:** 2026-10-16 · **Branch:** `paper/gcp-study` · **Work dir:** `/tmp/opencode/zenodo-vardictcpp/`
**Status: DONE_WITH_CONCERNS** (all downloads verified; one md5-expectation typo in task text; PDF fetched try-1)

## 1. md5 verification

| File | Source | md5sum (observed) | Expected | Verdict |
|---|---|---|---|---|
| `README.md` (4,370 B) | `zenodo.org/api/records/22955424/files/README.md/content` | `8bb795e30739754a81ec48f45e542b53` | task text `8bb7954e27839754e2d0db9026343a896` | **MATCH vs Zenodo API-declared checksum; MISMATCH vs task text (see concern C1)** |
| `Figure_data.zip` (13,637,860 B) | `zenodo.org/api/records/22955424/files/Figure_data.zip/content` | `8b79769889e4424e2d0db9026343a896` | `8b79769889e4424e2d0db9026343a896` | **MATCH** |
| `Archive.zip` | — | — | `7bbcd523e3c3c38b8ce56f926ddf30c8` | **NOT downloaded** (4.07 GB; explicitly out of scope). Zenodo API record listing confirms the file exists at 4,075,298,167 B with declared `md5:7bbcd523e3c3c38b8ce56f926ddf30c8` — matches plan value. |

sha256 (for MANIFEST convenience, not required): `biorxiv.pdf d9505294…3dca2`, `README.md 702c2d77…cb70e`, `Figure_data.zip ed33aa22…a695b4811c` (full values in §5).

## 2. Finding (a) — Fig-2 germline field-F1 truth source

**Answered by README.md + vendored PDF (no longer TBD).** Fig-2 is **NOT** the MV4-11/SRR panel:

- Sample: **GIAB HG002 chromosome 20 exome** (Agilent SureSelect v5), BAM `bams/hg002.exome.20.bam` (in Archive.zip).
- Truth: **NIST GIAB v4.2.1** (`truth20.norm.vcf.gz`), scored with **`rtg vcfeval`** over high-confidence exome-callable regions (**depth ≥ 20; 2.31 Mbp; 2,443 truth variants** — PDF §3).
- So the F1 denominator = the 2,443 GIAB v4.2.1 chr20 high-confidence truth variants.
- Fig-2 caller set (README `tool_vcfs/` + PDF): **octopus, vardictcpp (= VarDictJava, identical), FreeBayes, LoFreq, Mutect2 (tumor-only)** — i.e. octopus/FreeBayes/vardict F1s all share this same truth denominator.
- Vendored summary (`caller_accuracy.tsv`): SNV F1 — octopus 0.9925, vardictcpp=Java 0.9854, FreeBayes 0.9822, LoFreq 0.9771, Mutect2 0.6005; Indel F1 — octopus 0.9071, FreeBayes 0.8277, vardict 0.8009, LoFreq 0.7761, Mutect2 0.4874. (Existing Figure_data values; our re-measure at Wk1 will confirm.)
- **Plan correction:** `gcp-study-protocol.md` §§4–5 currently say "Fig-2 truth TBD from PDF (MV4-11); else parity-only" and "require the paper's truth source — identify from GCS-vendored PDF before enrolling MV4-11 for Tier-3". That TBD is now resolved: **MV4-11 (SRR8657348) is the Fig-1 tumor-only performance arm, not the Fig-2 accuracy arm.** Tier-3 enrollment for the Fig-2 replicate should be **HG002 chr20 exome vs GIAB v4.2.1** (already in our Tier-G truth family), not MV4-11. Recommend a one-line prereg correction at Wk1 freeze.

## 3. Finding (b) — Archive.zip layout (README + record listing + PDF; zip NOT downloaded)

Per README.md §§"What is here vs in Archive.zip" + PDF §5 Data Availability:

- **Analysis-ready BAMs: YES.** `bams/`: `SRR8657348.bam`, `SRR15006376.bam`, `SRR15006386.bam`, `hg002.exome.20.bam` — **all hg19**.
- **Per-tool performance outputs: YES.** `performance_benchmark/raw/single_sample/` and `.../somatic/`: each caller's VCF/TSV plus `.wall`, `.peak`, `.rc` timing files (Figs 1 & 3 inputs).
- **Reference + truth sets: YES** (PDF: "the BAMs, the reference and truth sets, the caller output files, and the benchmark scripts and results, are available on Zenodo").
- **Target BEDs: YES per PDF** ("The exact caller commands, target BEDs, and the recorded wall-clock and peak-RSS values … are in the Zenodo archive") — but README names no BED paths, so exact BED filenames/locations inside Archive.zip are unconfirmed until Wk1 (see concern C3).
- **Wk1 action stands:** check Archive.zip for the analysis-ready BAMs/BEDs before rebuilding SRR BAMs from FASTQs — the BAMs at least are confirmed present, so SRR restaging may be reducible to download + digest-verify.

## 4. Finding (c) — SRR accessions + reference build (PDF-verbatim)

- **Tumor/normal AML panel:** `SRR15006375`, `SRR15006376` (tumor), `SRR15006386` (normal), `SRR15006540` — adult AML bone-marrow, FLT3-ITD, FiLT3r study, **BioProject PRJNA742684, Illumina MiSeq** (PDF §5; matches plan Tier-P pin).
- **Germline/performance arm:** `SRR8657348` — MV4-11 AML cell line, CCLE, **BioProject PRJNA523380** (PDF §5). Used in Fig-1 (tumor-only, 773,667 target regions) and Table-3 chr1 `-R` blowup (`chr1:1-249250621`).
- **Reference build: hg19 throughout** (PDF Tables 1–3, Figs 1–3 captions). Truth-side naming uses GRCh37 (`HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi` present in Figure_data `tool_vcfs/`).
- **Table-1 parity counts (PDF-verbatim, confirm plan §1.1 snippet quarantine lifted for these):** SRR15006386 11,732 regions / 4,014 rows; SRR15006375 12,292 / 4,056; SRR15006376 105,949 / 12,739; SRR15006540 135,145 / 16,031; SRR8657348 773,667 / 54,448 — all 100% byte-identical. SRR8657348 total row 91,288 / 91,288 (100%).
- **Measurement provenance (for prereg operating-point notes):** Table 2 = median of 5 replicates, `/usr/bin/time -v`, Java on **JDK 25 with `-Xmx 8g`**, AMD Ryzen 7950X; Fig-1 measured on **AMD Epyc 9654**; Table 3 Java at `-Xmx 200g` (107.6 GB RSS, 6:10) vs cpp 10.6 GB (2:03) on chr1 `-R`. Paper DOI as printed: **10.64898/2026.09.23.753771** (matches plan; PDF header confirms). Code: `github.com/MHH-Bioinformatics-Hematology/vardictcpp`, MIT, Bioconda + Galaxy.

## 5. PDF outcome

- **SUCCESS on try 1 (no backoff needed):** HTTP 200, saved to `/tmp/opencode/zenodo-vardictcpp/biorxiv.pdf`, **396,827 B, 9 pages, PDF-1.5**, sha256 `d9505294ec26f02810ce67e0928c9597f8ba8bedcb81b40139bc77cf25e3dca2`.
- Biorxiv rate-limit (5× 429 from author network) did not reproduce from this network. No hammering performed (single request).
- **Recommended MANIFEST line (Wk1):** `biorxiv-10.64898-2026.09.23.753771v1.pdf sha256:d9505294…3dca2` — full: `d9505294ec26f02810ce67e0928c9597f8ba8bedcb81b40139bc77cf25e3dca2`.

## 6. Figure_data.zip listing (46 files, 21 MB unzipped; names + sizes)

```
computed_output/fig1_performance/caller_speed.tsv (292 B)
computed_output/fig1_performance/caller_speed_somatic.tsv (250 B)
computed_output/fig2_accuracy/caller_accuracy.tsv (359 B)
computed_output/fig2_accuracy/tool_vcfs/:
  cpp.vcf 3,345,571 B | java.vcf 3,345,571 B  (equal size — consistent w/ byte-identical claim)
  cpp.germ.vcf.gz 209,194 B | java.germ.vcf.gz 209,196 B (+.tbi ~8 KB each)
  cpp.norm.vcf.gz 630,883 B | java.norm.vcf.gz 630,883 B (+.tbi)
  cpp.pass.vcf.gz 572,309 B (+.tbi); cpp.germ.snp/indelsplit VCFs (+.tbi)
  fb.vcf 1,888,985 B | lofreq.vcf 175,251 B | lofreq.gt.vcf 193,506 B
  m2.vcf.gz 147,478 B | m2.filt.vcf.gz 159,562 B (+.tbi)
  octopus.vcf 281,590 B | vd_small.vcf 51,972 B
  truth20.norm.vcf.gz 3,633,042 B | truth.snp.vcf.gz 2,912,735 B | truth.indel.vcf.gz 651,549 B (+.tbi each)
  HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi 1,645,774 B (index only — full benchmark VCF lives in Archive.zip)
computed_output/fig3_somatic_venn/: vardict.snv 14,108 B (737 lines) | mutect2.snv 2,325 B (125) | lofreq.snv 2,060 B (111) | octopus.snv 307 B (17)
scripts/: benchmark.sh (1,888 B) | giab_accuracy.sh (3,507 B) | make_caller_figs.py (9,102 B) | make_venn.py (891 B) | data/{caller_accuracy.tsv, caller_speed.tsv, caller_speed_somatic.tsv} (mirror copies)
```

No deep analysis performed beyond line counts and the three summary TSVs quoted in §2/§4.

## 7. Concerns for Wk1 owner

- **C1 (integrity-note, non-blocking):** Task text's README.md expected md5 (`8bb7954e…`) does not match the Zenodo API-declared checksum (`md5:8bb795e30739754a81ec48f45e542b53`); the download matches the API declaration byte-for-byte, so this looks like a transcription typo in the task, not a corrupt download. Recommend correcting the expected value wherever it is recorded.
- **C2 (sample coverage):** README's Archive `bams/` list names only 4 BAMs (SRR8657348, SRR15006376, SRR15006386, hg002.exome.20) but PDF Table 1 covers 5 samples incl. **SRR15006375 and SRR15006540**. Their BAMs may live under different names/paths in Archive.zip or only as per-tool outputs. Verify at Wk1 before assuming all 5 BAMs are vendored.
- **C3 (BED paths):** PDF asserts target BEDs are in the archive; README names none. Enumerate BED paths inside Archive.zip at Wk1 (needed for Tier-P `--chunk`/BED-vs-`-R` prereg row).
- **C4 (hg002.exome.20 provenance):** No source SRR/accession for the HG002 chr20 exome BAM is stated in README or the extracted PDF text; check `scripts/giab_accuracy.sh` + Archive.zip at Wk1 (matters for BAM-provenance prereg gate).
- **C5 (record DOI variance):** README cites `doi:10.5281/zenodo.21836659` while the plan cites concept DOI `10.5281/zenodo.21836658` / version DOI `10.5281/zenodo.22955424`. Likely v1-record vs concept-record numbering — confirm canonical citation DOI at manuscript stage.
- **C6 (scope reminder):** Fig-3 somatic Venn has **no orthogonal truth set** (PDF states this explicitly) — stays a concordance characterization, not an accuracy rank; already consistent with plan §2 (Tier-3 somatic = SEQC2/HG008, not SRR76/86).
- Skipped per instructions: Archive.zip download, deep Figure_data analysis, GCS/MANIFEST writes (Wk1 owner), commit (scratch only).
