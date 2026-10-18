# Task 5a report — downloads + verifications + commit-hash arbitration (paper/gcp-study)

**Status: DONE** · 2026-10-18 · branch `paper/gcp-study` · no commit (report only)
· Requirements: protocol §4 (Tier-P + Tier-E + RNA arm) + §10-Wk1 + `gcp-prereg.md` §A (verbatim accessions/flags/md5s)

## 0. Summary

| Item | Result |
|---|---|
| RNA FASTQs (R1+R2, rep-1 mRNA) | DOWNLOADED (VM) → md5-OK both → in GCS at exact sizes |
| 6 Tier-P BAMs+bai | Extracted from GCS Archive.zip → BAM md5s OK vs `MD5SUMS.txt` → in GCS, 15/15 objects md5 MATCH |
| Truth files (conf20/target20/truth20) | In GCS, md5 MATCH |
| Commit-hash arbitration | **VERDICT: `ae93fe39dc`** — fresh builds reproduce archived raw stream md5 exactly at -th 1 AND -th 8; `5ab3759` does not |
| VM | created spot, deleted after; `instances list` = 0 |
| /tmp | 69,577,785,522 B FASTQs deleted; task5a dir now 20 MB (build trees + logs + md5.in) |
| Spend | egress $0.46 + VM ~$0.25 + storage ~$1.35/mo → running total ≈ $1.45/mo, « $5 task cap and $100 checkpoint |

## 1. RNA FASTQs (Tier-G RNA arm + cross-tier matched LCL)

Source dir (protocol §4 verbatim):
`ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_RNAseq/AshkenazimTrio/HG002_NA24385_son/UNC_Illumina/reads/`
(+ `md5.in`, 4,500 B, kept at `/tmp/opencode/task5a/rna/md5.in`)

| File | Bytes (exact) | Expected md5 | VM-side `md5sum -c` | GCS size | GCS hash |
|---|---|---|---|---|---|
| `…AGGTTATA-CAGTTCCG_R1_001.fastq.gz` | 34,164,581,094 | `5a3f63aecf5ca8daebfaeaa005089348` | OK 2026-10-18 17:07 UTC | 34,164,581,094 | CRC32C `H31Swg==` (composite, no MD5 — see §5-C1) |
| `…AGGTTATA-CAGTTCCG_R2_001.fastq.gz` | 35,413,204,428 | `1d6a09975bea9065ed19888762e8255f` | OK 2026-10-18 ~20:3x UTC | 35,413,204,428 | CRC32C `A2zJqQ==` (composite, no MD5 — see §5-C1) |

GCS: `gs://vardict-paper-trial-bams/rna-hg002/` — 2 objects, 69,577,785,522 B (64.80 GiB).
Local /tmp copies md5-verified identically before deletion (R1 `5a3f63ae…`, R2 `1d6a0997…`).

### Download saga (stated plainly)
1. Home-link `curl -C -` R1 pull stalled twice (byte-identical 27,771,871,232 B, mtime 11+ min, no transfer process) — controller-ordered stops.
2. **Dual-writer corruption:** two local processes wrote the same R1 path concurrently; md5 mismatch caught it. Clean re-pull with single-writer `aria2c` → R1 md5 OK, R2 md5 OK (35,413,204,428 B). Lesson: one writer per path, always.
3. Home-link `gcloud storage cp` uploads SIGTERM-killed twice with 0 objects landed (70 GB at 5.7 MB/s home uplink ≈ 3.5 h exposure per attempt).
4. Pivot (controller order): ephemeral VM downloader — R1+R2 fresh-pulled on VM at ~75–90 MB/s, md5-gated, gsutil-uploaded intra-region. Full success.

## 2. SRR BAMs + truth (Tier-P, from GCS Archive.zip — no SRA re-download)

Pulled `gs://vardict-paper-trial-refs/vardictcpp/Archive.zip` (4,075,298,167 B, md5 `7bbcd523e3c3c38b8ce56f926ddf30c8` MATCH;
one-time GCS→internet egress 3.80 GiB × $0.12 = **$0.46**). Extracted ONLY `bams/*.bam`+`.bai`,
`truth/*.bed`+`*.vcf.gz`, `scripts/*` (+ `reference/hg19.fa[.fai]` and 2 archived vardictcpp somatic
`.tsv.gz` solely to run the arbitration — deviation D1, deleted after; local Archive copy deleted).

BAM md5s vs archive `MD5SUMS.txt`: **all 6 + bais OK**
(`SRR15006375.bam 23cb57e7…`, `SRR15006376.bam c735d817…`, `SRR15006386.bam 98625e7f…`,
`SRR15006540.bam 69828a79…`, `SRR8657348.bam f4b66ac8…`, `hg002.exome.20.bam fa809088…`;
hg19.fa `530d89d3…` matches PROVENANCE exactly).

GCS (all Standard us-central1, single-generation after composite re-upload fix — first upload landed
10-component composite with CRC32C-only, re-uploaded single-stream):
- `gs://vardict-paper-trial-bams/tier-p/` — 12 objects (6 BAMs + 6 bais, 2.82 GB), **GCS base64-md5 decode-match vs local/MD5SUMS for all 12**
- `gs://vardict-paper-trial-refs/tier-p-truth/` — `conf20.bed`, `target20.bed`, `truth20.norm.vcf.gz` (3.9 MB), **all 3 md5 MATCH**
- Total: **15/15 objects verified.**

## 3. Commit-hash arbitration — VERDICT `ae93fe39dc`

- Cloned `MHH-Bioinformatics-Hematology/vardictcpp` (full history); both SHAs exist as commits.
- Built EACH at its commit in new conda env `paper-vardictcpp-build` (cmake 4.4.4, htslib 1.24, system g++ 11.4.0
  — vs archive's htslib 1.22 + g++ 13.3.0, recorded as drift D3; pinned envs untouched).
- Invocation verbatim from `somatic_speed.sh` (Table-1 parity cell, SRR15006386.bam as normal):
  `vardictcpp -G hg19.fa -f 0.01 -N 'SRR15006376|SRR15006386' -b 'SRR15006376.bam|SRR15006386.bam' -c 1 -S 2 -E 3 -g 4 <BED> -th {1|8}`
  with covered BED regenerated via verbatim `mkbed()` → **105,949 regions / 21.36 Mb** (matches Table 1 exactly).
- Results (uncompressed stream md5 vs `zcat` of archived raw):
  - archived `…vardictcpp.1core.tsv.gz` stream: `6674ec0b40111c986e01f90bd00b65c1` (gz md5 `595d3138…` = MD5SUMS)
  - `ae93fe39dc` -th 1: `6674ec0b…` **MATCH** (14,126 rows = published call count)
  - `ae93fe39dc` -th 8: `6674ec0b…` **MATCH** (thread-deterministic)
  - `5ab3759` -th 1: `3c13cb23…` MISMATCH (14,127 rows, net +1 extra row vs archive)
- **Verdict: the published raw outputs were produced by `ae93fe39dc`** (branch `fix/somatic-amplicon-parity`).
  Top-level archive `README.md` claim (`5ab3759`) is stale; `PROVENANCE.md` (`ae93fe39dc`) is correct.
  Cite `ae93fe39dc`; prereg B6 closes on this digest match, not prose.

## 4. VM lifecycle

`paper-dl-tmp`, us-central1-a, **e2-standard-2, SPOT** (on-demand fallback not needed), 60 GB→200 GB
pd-standard boot disk (resized live: 60 GB pd-standard I/O-starved — md5 of 35 GB took 71 min, uploads
disk-bound at ~12 MB/s; 200 GB ≈ 3× IOPS baseline; +~$0.02). SA `vardict-agent@…` as attached SA
(default compute SA lacked `iam.serviceAccountUser`). Start 2026-10-18 16:18 UTC → deleted same day
after R2 verified in GCS; wall ≈ 4.4 h × ~$0.05/h + disk ≈ **~$0.25 total**. `instances list` = **0 VMs**
(paper-star-build already gone, nothing else appeared). No idle burn.

## 5. Spend accounting (running totals)

| Item | Amount |
|---|---|
| One-time egress (Archive.zip pull 3.80 GiB) | $0.46 |
| VM paper-dl-tmp (spot + disk) | ~$0.25 |
| New Standard storage | RNA 64.80 GiB ($1.30/mo) + tier-p 2.87 GiB ($0.06/mo) = **~$1.35/mo** |
| Ingress to GCS | $0 (free) |
| Prior spend (vendor + refs) | ~$0.08–0.10/mo |
| **Running total** | **≈ $0.71 one-time + ≈ $1.45/mo** — « $5 task cap, « $100 checkpoint |

## 6. Concerns / deviations

- **C1 (integrity record):** gsutil auto-parallel-composites files >150 MB → RNA GCS objects carry CRC32C
  only, no server MD5. Integrity rests on: VM-side `md5sum -c` OK for both files + exact byte sizes in GCS
  + local md5 matches + gsutil per-chunk CRC32C transfer checks. Future >150 MB uploads: pass
  `gsutil -o GSUtil:parallel_composite_upload_threshold=0` (or working single-stream `gcloud storage cp`)
  to retain server MD5s.
- **C2 (undiagnosed):** `gcloud storage cp` of the 34 GB file from the VM hung at ~0 B/s across two
  attempts (45+ min, process alive) while 10 MB uploads worked and `gsutil cp` worked. Uploader-specific
  stall, root cause unknown — avoid `gcloud storage cp` for multi-GB single objects until retested.
- **C3 (disk):** 60 GB pd-standard is I/O-inadequate for 35 GB file hashing/upload streaming (71-min md5);
  resized to 200 GB live. Future VM tasks: request ≥200 GB or balanced/SSD boot disks up front.
- **D1 (extract deviation):** pulled `reference/hg19.fa[.fai]` + 2 archived vardictcpp `.tsv.gz` beyond the
  task's ONLY list — required to execute the arbitration; all deleted after (Archive.zip removed, x/ removed).
- **D3 (toolchain drift):** arbitration built with htslib 1.24 + g++ 11.4 vs archive's 1.22 + 13.3.0 —
  `ae93fe39dc` still byte-identical, so immaterial to the verdict; recorded for the methods section.
- **C4:** /tmp FASTQs (69,577,785,522 B) deleted per close-out order after GCS verification; VM-side
  copies died with the VM by design. `gs://vardict-paper-trial-bams/rna-hg002/` (versioned bucket) is now
  the sole copy — both objects at exact byte sizes with md5-gated pull provenance (§1). Re-pull from NCBI
  FTP remains available if ever needed.
