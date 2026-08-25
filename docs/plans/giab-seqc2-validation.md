# External Validation Plan — GIAB / SEQC2 (deferred from v0.1.0; executes as T167–T169)

> **Superseded by the execution-grade specification:**
> [`docs/plans/2026-08-24-g4-execution.md`](2026-08-24-g4-execution.md) — subagent-ready
> tasks with problem/theory/implementation/validation per task, media-stack coexistence
> milestones (M0–M3), and dispatch guide. This file remains as background rationale.

> **Status: PLANNED, NOT STARTED.** Deferred 2026-08-24 (multi-GB downloads, hours of
> runtime). This file is written to be executable cold, by an agent or human, without
> session context. Read `CLAUDE.md` first; §0b integrity governs every claim below.

## 1. Goals & success criteria (maps to PLAN.md Phase 12)

| Task | Goal | Pass criterion |
|---|---|---|
| **T167** | Full SEQC2 HCC1395 WES, tumor\|normal somatic run, both engines | `compare.sh` **MATCH** (float-only diffs within ε=1e-5 after HALF_EVEN rounding) |
| **T168** | GIAB HG008 (PDAC) somatic slices, both engines | `compare.sh` MATCH per slice |
| **T169** | Truth-set concordance (secondary tier) | `rtg vcfeval` + `truvari bench`: vardictrust precision/recall ≈ VarDictJava 1.8.3 (delta ≈ 0; any delta → new DIVERGENCES entry + investigation before proceeding) |

Primary bar stays **equivalence to the oracle** (§5); truth-set tier only confirms we
didn't faithfully reproduce a *wrong* result. These three tasks are prerequisite
evidence for the T170 acceptance sweep and the M3 / `v1.0.0` gate definition
(PLAN.md: M3 requires all-mode MATCH on corpus **+ SEQC2/GIAB slices**).

## 2. Prerequisites & pre-flight checks

```bash
df -h .                      # need ≥80 GB free scratch
nproc                        # runtime scales with cores (dev box: 12c/24GB)
conda env list | grep vardict-oracle   # oracle env must exist (harness/setup_oracle.sh)
command -v samtools bcftools           # slicing/CRAM work
```

Tool installs (one-time, conda): `rtg-tools` (Java), `truvari`, `bcftools`.
Reference genome: GRCh38 FASTA+`.fai` (~3 GB download if absent locally; reuse
`refdata/` if present). Contig naming must match the BAMs (`chr1` vs `1`).

## 3. Data acquisition — `harness/fetch_data.sh` (new script to write)

Downloads with md5 verification, writing `benchmarks/MANIFEST` (URLs, md5s, exact
commands, dates) per PLAN.md §benchmarks:

| Dataset | Source pointer (verify at execution time — FTP layouts drift) |
|---|---|
| SEQC2 HCC1395 WES tumor/normal BAMs + truth VCF + high-conf BED | NIST/GIAB somatic hosting: `ftp.ncbi.nlm.nih.gov/giab/…/data_somatic/HCC1395…` (PLAN.md cites this tree) |
| GIAB HG008 PDAC somatic BAM(s) + truth | `…/data_somatic/HG008/Liss_lab/` (path already validated in PLAN.md research) |

Rules:
- Slice-first strategy: fetch full files once, then `bash harness/slice.sh <in.bam> <chr:start-end> <out.bam>` around benchmark loci (exome BED chunks on chr20/chr22 first) — validate the pipeline on slices BEFORE committing to full-WES runs.
- Never commit data. MANIFEST only.
- If a URL 404s: search the GIAB index page for the renamed asset, record the redirect in MANIFEST. No silent substitutions.

## 4. Execution steps

1. `bash harness/fetch_data.sh` → verify md5s against published values.
2. **Slice phase:** pick 3 exome regions with known variants (truth VCF intersect).
   For each: `bash harness/run_oracle.sh ref.fa "tumor.bam|normal.bam" bed out.vcf sample`
   then `bash harness/compare.sh ref.fa "tumor.bam|normal.bam" bed sample` → expect MATCH.
3. **Full-WES phase (T167):** same commands with whole-exome BED; expect long runtime
   (hours/engine on 12c). Record wall time + peak RSS into `benchmarks/README.md`
   (extends the P10 slice evidence to a production-scale dataset).
4. **HG008 slices (T168):** repeat step 2 pattern.
5. **Truth tier (T169):**
   ```bash
   rtg vcfeval -b truth.vcf -c ours.vcf -t SDF --squash-ploidy -o eval-ours
   rtg vcfeval -b truth.vcf -c java.vcf -t SDF --squash-ploidy -o eval-java
   truvari bench -b truth.vcf -c ours.vcf -o bench-ours   # SV incl.
   ```
   Compare PR metrics between engines (they should agree because outputs match);
   report absolute numbers vs truth honestly. Any engine-vs-truth surprise that is
   NOT shared by Java = our bug → fix loop under §7 port discipline.
6. **CRAM twin piggyback (feeds T156):** convert one sliced BAM to CRAM, rerun
   compare.sh, expect MATCH (#249/#380 evidence).

## 5. Failure handling

- MISMATCH → stop, minimize, follow the standard debug loop; fail-set documentation
  rules apply (a real mismatch reported honestly beats a fake pass; never weaken
  comparators).
- Float-only diffs beyond ε → check locale/rounding path first (P2/P8 history),
  then open DIVERGENCES entry.
- Runtime blowup → tile-size tuning is allowed (perf-only, no output change);
  re-verify determinism after.

## 6. Deliverables on completion

- `benchmarks/MANIFEST` + updated `benchmarks/README.md` (full-dataset rows)
- DIVERGENCES.md status updates (if anything new surfaced)
- README Validation-status section upgraded: "external truth sets" line
- Unblocks T170 acceptance sweep and the annotated **v1.0.0** tag (T173)
