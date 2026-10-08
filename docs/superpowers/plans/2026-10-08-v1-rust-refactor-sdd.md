# VarDictRust v1.0.0 — Pure-Parity Release Plan (SDD, revision 3)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`.
> - One fresh implementer per task, followed by spec review, then code-quality review.
> - Ledger: `.superpowers/sdd/2026-10-08-v1-rust-refactor-sdd/progress.md`. Steps use `- [ ]`.
> - **Controller-only** tasks (T0, T25) and every ⛔ step need the human's explicit approval in chat.

**Supersedes** revision 2 (2026-10-08 06:44Z). Incorporates all user decisions from 2026-10-08 (ε = 1e-6, SEQC2 gates v1.0.0, keep working until parity is achieved, cargo-dist adopted, CodSpeed dropped and replaced by a pre-registered performance budget, no Co-Authored-By trailers) plus external research (rtg `--squash-ploidy`, cargo-dist v0.33.0, rewrites.bio governance, VCF header parity, edge cases). Every finding is listed in **§2 Discrepancy Register** together with the task that resolves it.

**Goal.** Ship `vardictrust` **v1.0.0** with:
- output identical to VarDictJava **v1.8.3** (`0886261`, submodule `009e017`) on valid inputs;
- more robustness where Java crashes;
- parity evidence on the corpus plus GIAB HG008 and SEQC2 HCC1395;
- a GitHub musl release via cargo-dist and crates.io publication.

**Architecture.** Engine work is limited to four things:
1. Emulating Java's skip-read for the PR #402 family.
2. A deterministic `java.util.HashMap` iteration-order model (s1a).
3. Exact Fisher/SBF math matching R's `dhyper` / `relErr` summation (mandatory under ε = 1e-6).
4. Header and formatting residues.

The rest is harness hardening, external-truth evidence and release engineering. Each piece of evidence is a JSON record bound to an **engine fingerprint**, so it goes stale automatically if the engine changes.

**Tech stack.**
- Rust 2021 workspace (`vardict-core`, `vardict-cli`), with a bash/python3 harness.
- Conda env `vardict-oracle` (exists) and `g4-tools` (new): samtools, bcftools, htslib, bedtools, rtg-tools, bwa, gitleaks.
- `cargo-dist` (v0.33.0), `gh`.

---

## 0. Locked decisions (user, 2026-10-07 & 2026-10-08)

| ID | Decision |
|---|---|
| **A1** | **v1.0.0 = pure byte parity with v1.8.3.** No `--fix` flags and no output-changing corrections. Obvious upstream bugs are *replicated* and catalogued in `ROADMAP.md` for v1.1.0. There, each fix gets an opt-in `--fix <id>` (default off) plus a diff-confinement gate. |
| **A2** | Both datasets: HG008 first, then SEQC2 HCC1395, aligned by us from ENA FASTQ. |
| **A3** | **The host is not touched.** No crontab, systemd, docker or service changes, and no quiescing. Heavy jobs run via `harness/courtesy.sh`, which applies `nice -n19 ionice -c3` and refuses to start unless there is ≥150 GB free disk and ≥12 GB available RAM, and the run won't overlap 03:40–05:15. Threads are capped at 8: pass `-th 8` to **both** engines, because vardictrust defaults to all 20 cores. |
| **A4** | Full release: push `main`, open the PR, CI green, merge, tag `v1.0.0` (which triggers cargo-dist `release.yml`). The first crates.io publish is manual with the user's token; trusted publishing is set up after that. |
| **A5** | **ε = 1e-6 (raw float comparison).** String equality OR finite float diff `\|a - b\| <= 1e-6 (+ 1e-12)`. No 5-dp `HALF_EVEN` pre-rounding masking 5-decimal float differences. Applies to both `compare_vcf.py` and `compare_var.py`. Tighter than CLAUDE §5's 1e-5 (allowed by §0b). OQ-1 resolved. |
| **A6** | **SEQC2 gates v1.0.0.** External HCC1395 WES parity & concordance must pass before tagging v1.0.0. OQ-2 resolved. |
| **A7** | **Keep working until parity is achieved.** Mismatches, budget failures, and scale failures are not stop/waiver reasons; they trigger immediate fix tasks (`superpowers:systematic-debugging`). Timeboxes become progress reporting checkpoints. Stop only for controller actions (user approvals, credentials, host policies). OQ-10 resolved. |
| **A8** | **Adopt cargo-dist (`dist`).** Use cargo-dist for automated musl release artifacts, checksums, shell installers, and OIDC trusted publishing to crates.io. OQ-6 resolved. |
| **A9** | **Pre-registered performance budget replaces CodSpeed.** Rust wall <= Java wall and Rust peak RSS <= Java peak RSS on Tier-2 evidence (`hg008-tier2` and `seqc2-exome-tier2`). CodSpeed SaaS is dropped. If budget fails, perform profiling and optimization (keep working). OQ-7 resolved. |
| **A10** | **No `Co-Authored-By` trailer in new commits.** Trailers removed from commit templates. CLAUDE §18 amended. AI transparency maintained in README "AI Disclosure" and `CITATION.cff`. OQ-9 resolved. |
| **A11** | **Crate name = `vardictrust` (direct `cargo install vardictrust`).** Package `vardict-cli` renamed to `vardictrust` in `Cargo.toml`. Binary and package names are identical. Direct installation via `cargo install vardictrust` is supported on crates.io. OQ-3 resolved. |
| **A12** | **Release target `x86_64-unknown-linux-musl` with multi-arch CI matrix.** Parity-validated musl artifact packaged by cargo-dist for release. `.github/workflows/ci.yml` introduces a build/smoke matrix covering `x86_64-unknown-linux-gnu`, `x86_64-unknown-linux-musl`, and `aarch64-unknown-linux-musl` (cross-compilation check) to ensure portability across CPU architectures. OQ-11 resolved. |
| **A13** | **Performance baseline: Rust <= Java wall time & max RSS.** Confirmed baseline: CodSpeed dropped; pre-registered budget `perf_budget.py` requires Rust wall time <= Java wall time and Rust max RSS <= Java max RSS. OQ-13 resolved. |

## 1. Global constraints (every task)

- **G1 Parity gate.** `compare_vcf.py` must exit 0:
  - exact: CHROM/POS/REF/ALT multiset, FILTER, INFO/FORMAT key sets, non-numeric values, and VCF header structure;
  - numeric values: raw float comparison within ε = 1e-6 (A5);
  - `--strict-bytes` (T11) is **report-only** (surfaces literal formatting differences within ε).
- **G2 Harness only tightens (CLAUDE §0b).**
  - Each harness change gets its own commit, a red→green test, and a green **R-GATE**.
  - Never relax ε, skip records, or hand-edit `golden/`.
  - **Decide by exit codes**, never by grepping MATCH text.
- **G3 Host etiquette (A3).**
  - Anything expected to run >10 min goes through `EST_MIN=<n> bash harness/courtesy.sh …` with `-th 8`.
  - Java runs with `JAVA_OPTS=-Xmx8g` (31 GB host, ~15 GB available).
- **G4 Data.** All data lives under `scratch/` (gitignored). Only these are committed: `benchmarks/MANIFEST`, `benchmarks/evidence/*.json`, small BEDs, and KB-size fixtures.
- **G5 No secrets** in files, logs or commits. The crates.io token is entered by the user in their own terminal (`cargo login`).
- **G6 CLAUDE.md engine rules.**
  - No panics in `vardict-core`.
  - `Pos0`/`Pos1` for positions.
  - `IndexMap`/`BTreeMap`, or the T10 `JavaHashMap` where Java iterates a HashMap. Never std `HashMap` in output paths.
  - `/// ports X.java:a-b` comments.
- **G7 Commit template** (CLAUDE §15/§18; no trailer per A10). For docs-only commits, the validation line reads `Validation: N/A (docs only); fmt/clippy/test green`.
  ```
  T<nn>: <imperative summary>

  <why / what changed>

  Validated vs VarDictJava 1.8.3 (0886261/009e017), dataset <name>, cmd `<exact cmd>` -> <result>
  ```
- **G8 Evidence.** "Done" requires an observed exit code or diff, pasted into the task report. Never claim a run you did not do.
- **G9 Checkpoints (A7).** Timeboxes become progress reporting checkpoints, not stops. The implementer continues working, reporting status, diffs, and next steps at each checkpoint. Nobody quietly narrows scope or skips to a waiver.

### R-GATE (`harness/rgate.sh`, created in T2). Run after every harness or engine change.

All four must hold:
1. `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` and `cargo test --workspace` all pass.
2. `harness/run_integration.sh` prints exactly `PASS: 91  FAIL: 1  SKIP: 0`. The only `MISMATCH:` line is the documented stale fixture `Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt`.
3. `compare.sh` exits 0 on `fixtures/micro`.
4. `harness/run_slice.sh` exits 0 for every id in `RGATE_SLICES`. The default is `s1b s5 s7`; T10 adds `s1a`.

---

## 2. Discrepancy Register

Severity: **C** critical · **I** important · **M** minor.

| ID | Sev | Finding (evidence) | Resolution | Task |
|---|---|---|---|---|
| D-01 | C | `compare.sh` writes to `mktemp` and then deletes it. The persisted `ours.vcf`/`java.vcf` that later steps consumed never existed (`compare.sh:5,14`). | `KEEP_DIR` persists `java.vcf`, `rust.vcf`, `compare.txt`, `*.time`. | T5 |
| D-02 | C | On a mismatch, `set -e` aborts at `compare.sh:12` (`python3 …; rc=$?`), so there is no MISMATCH line, no diff and no cleanup. | Use `… \| tee … \|\| rc=$?` under pipefail; clean up with a trap. | T5 |
| D-03 | I | `grep -c MATCH` also counts `MISMATCH` lines, and each passing run prints 2 MATCH lines. | Decide by exit code only (G2). Evidence JSON records the result. | T5, T6 |
| D-04 | I | `\| tee` without pipefail masks exit codes (rev-1 plan steps). | pipefail everywhere; acceptance reads exit codes. | T5, T6 |
| D-05 | C | `MATCH (0 variants)` is a vacuous pass. | `--min-variants` / `MIN_VARIANTS` floors, plus a vacuous warning. | T3, T5 |
| D-06 | C | The inner `bash -c` pipelines in `run_oracle.sh` lack pipefail, so a Java crash gives an empty or partial VCF with exit 0. | `set -euo pipefail` inside both inner scripts. | T4 |
| D-07 | C | The `run_oracle.sh` splitter cuts multi-char options. `-th 8` becomes `-t h 8`, and Java `-t` means dedup. `-UN -VS -DP -mfreq -nmfreq -fisher -adaptor -chimeric -deldupvar` are cut the same way. | Whitelist the multi-char options taken from `CmdParser.java`. | T4 |
| D-08 | I | `compare_vcf.py` stores records in a dict, so duplicate keys collapse, and it never compares record counts. | Group by key, compare counts per key, then compare records pairwise after sorting. | T3 |
| D-09 | I | No harness timeouts. | `CMP_TIMEOUT` per engine (default 12h). | T5 |
| D-10 | I | `run_integration.sh` exits 1 on the expected 91/1 result. | `rgate.sh`/`acceptance.sh` parse the exact summary line and the exact fail-set. | T2, T6 |
| D-11 | I | `compare.sh --all` does not exist, and the rev-1 acceptance list was invented rather than taken from PLAN §15. | `acceptance.sh` maps PLAN §15 items 1–8 verbatim. | T6 |
| D-12 | C | `bam.rs:510,528,603` hardcode `/mnt/hbmon/...`, causing 3 failing tests. | `CARGO_MANIFEST_DIR`-relative paths. | T1 |
| D-13 | M | Stale paths and host model in `KICKOFF_PROMPT.md:16` and the G4 doc. | KICKOFF made repo-relative; G4 doc gets a SUPERSEDED banner. | T1 |
| D-14 | C | The `ci.yml` checkout has no submodules, and `.gitmodules` uses an **ssh** URL that CI cannot clone. | `submodules: recursive`; switch `.gitmodules` to https. | T1 |
| D-15 | C | PR #402 family. `modify.rs:1102-1105` emits a negative `rm`, and `parse.rs:208-214` keeps the stale cigar. Java throws on this read and skips it (`CigarParser.java:89`). | Emulate the skip at Java's exact throw point, proved by the oracle. | T9 |
| D-16 | C | UB-3, the s1a AF-tie. Java iterates a `HashMap<Integer,Vars>` (`ToVarsBuilder.java:93`, `SimplePostProcessModule.java:35`); Rust uses a BTreeMap. | `JavaHashMap` model validated against real Java; keep working until closed (A7). | T10 |
| D-17 | I | `--fisher` raw-cell forms (`NM=1.0`, `HIAF=0.0020`) are within ε but not byte-equal. | `--strict-bytes` audit, then port the Java formatters. | T11 |
| D-18 | C | Under ε = 1e-6 (A5), SBF 1-ulp residue at the 5th decimal is a fatal mismatch. | Mandatory port of R's exact hypergeometric summation (`dhyper`, `1+1e-7` relErr threshold). | T11 |
| D-19 | I | `README.md:26` says somatic mode writes `.var`, but it writes VCF (`lib.rs:585,638`). | Fix the README. | T8 |
| D-20 | C | `--version` is bare clap (`cli.rs:58`). It has no v1.8.3/Lai 2016 credit and no AI disclosure. | `long_version` plus an integration test. | T8 |
| D-21 | I | The README AI disclosure names Claude only, and the §18 trailer is hardcoded to Claude. | Build the tool inventory from repo evidence; the trailer names the actual tool. | T7, T22 |
| D-22 | C | The user's exception directive conflicts with CLAUDE §5/§11 and was not resolved in rev 1. | A1 (v1.0 pure), the `ROADMAP.md` v1.1.0 table, and policy amendments. | T7 |
| D-23 | C | No public SEQC2 WES BAMs exist. | Align ENA FASTQ SRR7890879/80 to GRCh38.d1.vd1 ourselves. | T13, T18 |
| D-24 | C | Rev 1 used the genome-wide HC BED (709,209 intervals, 2.49 Gb) as the WES region BED. | Engine BED = Agilent v6+UTR exome; eval BED = exome ∩ HC. | T19, T20 |
| D-25 | C | Rev 1 used one FASTA for both datasets (GIABv3: 218 contigs, N-masked; d1.vd1: 2,779 contigs). | Two references, two SDFs, and an @SQ-vs-`.fai` gate per dataset. | T13, T15, T18 |
| D-26 | I | The HG008 somatic truth is only DRAFT v0.3. | Label it DRAFT; use it for parity, not headline accuracy. | T17, T22 |
| D-27 | C | The rev-1 rtg gate used ΔF1 < 0.001 and parsed the wrong columns. It stripped no FILTERs and selected no sample. | Exact equality of the `None` row plus identical site lists; `--sample ALT,<tumor>`; strip FILTER; require a non-empty baseline. | T12 |
| D-28 | I | Oracle out-of-memory on large regions (#197/#384). | Cap engine BED windows at ≤20 kb; set `-Xmx8g`; add a RAM guard. | T12, T15, T19 |
| D-29 | C | Rev 1 quiesced the media stack. | Dropped per A3; use `courtesy.sh` instead. | T12 |
| D-30 | M | The host has 20 cores, but PLAN §15.5 says "12 cores". | Measure 1/2/4/8 threads; anything above 8 needs a waiver. | T7, T21 |
| D-31 | C | The remote is empty and there is no local `main`, but CI triggers only on `main`. Rev 1 tagged on `init-scaffold`. | T0 creates `main@4b295eb` and `release/v1.0.0`; T25 opens the PR. | T0, T25 |
| D-32 | C | The crates have no `description`, the path dependency has no `version`, and package is named `vardict-cli`. | Rename package to `vardictrust` in `Cargo.toml` (A11), add metadata, `version = "1.0.0"`, and run a `--dry-run`. | T23 |
| D-33 | I | `CHANGELOG.md` still says "[0.1.0] — unreleased", and its claims are unverified. | Release as `[1.0.0]` and verify every claim. | T23 |
| D-34 | C | Tracked meta files and no secret scan before a public push. | gitleaks plus a ⛔ review (OQ-8). | T25 |
| D-35 | I | ±2 kb windows yield 0 truth clusters. | Chained clustering plus an IG/HLA exclude. | T14 |
| D-36 | I | Java has no CRAM support, so the CRAM twin has no oracle. | Compare ours(CRAM) vs oracle(BAM), and ours(CRAM) vs ours(BAM). | T16 |
| D-37 | M | truvari and som.py: there is no SV truth. | Mark N/A with a waiver. | T6, T17, T22 |
| D-38 | I | CodSpeed is not configured, and user requested no external bench SaaS. | Pre-registered performance budget: Rust wall <= Java wall, Rust peak RSS <= Java peak RSS (A9). | T21 |
| D-39 | I | `cargo audit` is missing from CI and locally. | Run it locally and add it to `security.yml`. | T23 |
| D-40 | I | `AF_THR` cannot force a mismatch, because it applies to both engines. | A fake `cargo` on PATH; the oracle stays real. | T5 |
| D-41 | C | Rev 1 had placeholder steps. | Full code below; data-driven docs are cross-checked by a reviewer. | all |
| D-42 | C | Host and publish side effects were dispatched to subagents. | Controller-only ⛔ steps. | T0, T25 |
| D-43 | I | G4 R2–R6/F1–F2 requirements were silently changed. | Superseded explicitly and re-homed. | T1, T13–T24 |
| D-44 | I | `-th` was not passed to either engine; Rust defaults to all cores. | `-th 8` to both (after D-07 is fixed). | T4, T15, T19 |
| D-45 | I | Java exits after >10 exceptions (`Utils.java:261-277`). | Rust never aborts (CLAUDE §5); document in DIVERGENCES. | T9 |
| D-46 | I | The W4.6 slice set was session-local. | `harness/slices/` + `slices.tsv` + `run_slice.sh`. | T2 |
| D-47 | I | Evidence could predate a later engine change. | `engine_fingerprint.sh` value inside every evidence JSON. | T6 |
| D-48 | M | `--fisher` residues at VCF level are unmeasured. | Strict audit on micro and s1b with `--fisher`. | T11 |
| D-49 | I | PLAN specifies dist; hand-rolled `release.yml` lacked cross-compilation & OIDC publish. | Adopt cargo-dist v0.33.0 for musl builds, shell installer, and gated OIDC publishing (A8). | T23 |
| D-50 | C | VCF header parity: Rust emits extra `##reference` and `##contig` lines; INFO/FORMAT/FILTER descriptions differ from Java var2vcf. | `run_oracle.sh` passes `-G "$REF" -b <contigs.bed>` to var2vcf; header compare with allowlist (`##fileDate`); fix Rust header text to match Java exactly. | T5a |
| D-51 | I | Edge-case parity: empty BAM (0 reads), zero-coverage region, single-read BAM, contig boundary reads. | Automated test suite verifying identical exit 0 and empty/matching VCF records across all 4 edge cases. | T5b |
| D-52 | C | Float comparison ε = 1e-5 with 5-dp HALF_EVEN pre-rounding masked 5th decimal differences. | Tighten to ε = 1e-6 on raw floats with no r5 pre-rounding in `compare_vcf.py` and `compare_var.py` (A5). | T3 |
| D-53 | C | Under ε = 1e-6, SBF 1-ulp difference is a fatal mismatch. | Mandatory port of R's exact hypergeometric summation (`dhyper`, `1+1e-7` relErr threshold). | T11 |
| D-54 | C | `rtg vcfeval` on somatic calls requires `--squash-ploidy` for sites-only baseline vs tumor calls. | Add `--squash-ploidy` to `concordance.sh`. | T12 |
| D-55 | I | No continuous bench SaaS (CodSpeed dropped per user). | Pre-registered performance budget: Rust wall <= Java wall, Rust peak RSS <= Java peak RSS on Tier-2 runs. | T21 |
| D-56 | I | Automated release pipeline. | cargo-dist with x86_64-musl, shell installer, attestations, and OIDC publish. | T23 |
| D-57 | I | rewrites.bio governance gaps: no CITATION.cff, unpinned oracle env, no upstream policy in README. | Add CITATION.cff (Lai 2016 + VarDictJava), pin env exports in `benchmarks/env`, add upstream tracking policy in README. | T12, T23 |

---

## 3. Task graph & dispatch

```mermaid
flowchart LR
  T0["T0 branches (controller)"] --> T1 --> T2 --> T3 --> T4 --> T5 --> T5a --> T5b --> T6 --> T7 --> T8
  T8 --> T9 --> T10 --> T11
  T11 --> T12 --> T13 --> T14 --> T15 --> T16 --> T17
  T17 --> T18 --> T19 --> T20
  T20 --> T21 --> T22 --> T23 --> T24 --> T25["T25 publish (controller)"]
```

The chain is sequential because there is one shared working tree, and the engine fingerprint must be frozen before the data runs. T13's downloads may run in the background once T12 is done.

| Task | Label | Model | Checkpoint / Timebox |
|---|---|---|---|
| T0, T25 | controller ⛔ | — | — |
| T1–T5b, T6–T8, T12–T14, T22–T24 | dispatchable | inherit | 1 session |
| T9, T11 | dispatchable | pro | 1 session (T11 SBF math: checkpoint at 2h) |
| T10 | dispatchable | pro | checkpoint at 6h (A7: continue until closed) |
| T15–T21 | dispatchable, long-running via courtesy | inherit | per step |

---

## 4. Tasks

### T0 — Branches & ledger (controller)
- [ ] `git branch main 4b295eb && git switch -c release/v1.0.0 main`. The `[branch "main"]` stanza already exists in `.git/config`. Do not push.
- [ ] Create the ledger `.superpowers/sdd/2026-10-08-v1-rust-refactor-sdd/progress.md`: a header with A1–A10, then one line per task: `Txx | status | commit | evidence`.
- [ ] Commit the plan doc: `T00: plan rev 3 (resolves adversarial review D-01..D-57)`.

### T1 — Portable paths, CI submodules, stale docs (D-12, D-13, D-14, D-43)
- [ ] **Red:** `cargo test -p vardict-core bam::` → 3 failures.
- [ ] In the three tests at `bam.rs:510,528,603`, replace the literal path with:
  ```rust
  let bam_path = Path::new(concat!(
      env!("CARGO_MANIFEST_DIR"),
      "/../../VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam"
  ));
  ```
- [ ] `.gitmodules`: set `url = https://github.com/AstraZeneca-NGS/VarDictJava.git`, then run `git submodule sync`.
- [ ] `ci.yml` checkout step gets `with: { submodules: recursive }`.
- [ ] `KICKOFF_PROMPT.md:16`: replace the path with a repo-relative one.
- [ ] Prepend this banner to `docs/plans/2026-08-24-g4-execution.md`:
  > **SUPERSEDED (2026-10-08)** by `docs/superpowers/plans/2026-10-08-v1-rust-refactor-sdd.md`. The host model, paths and quiesce steps here are stale and MUST NOT be executed (decision A3). Re-homed: R2/R3 → T15/T19, R4 → T16, R5/R6 → T17/T20, F1 → T21, F2 → T22–T24.
- [ ] `git grep -nE '/mnt/hbmon|/Users/' -- ':!VarDictJava' ':!docs/plans/2026-08-24-g4-execution.md'` → 0 hits.
- [ ] **Green:** `cargo test --workspace` → 178/178 pass.
- [ ] Commit `T01: portable test paths, https submodule + CI checkout, supersede G4 doc`.

### T2 — Reproducible slice set + R-GATE (D-46, D-10)
- [ ] Copy the persisted BEDs from `.superpowers/sdd/2026-08-21-parity-to-v1.0/w4-sweep/` into `harness/slices/`:

  | Slice | Source BED |
  |---|---|
  | s1a | `na24385_20156208.bed` |
  | s1b | `na24385_38386265.bed` |
  | tiny | `tiny.bed` |
  | s3 | `colo829_chr1.bed` |
  | s4 | `hcc2218_chrX.bed` |
  | s5 | `zeroin_chr19.bed` |
  | s6 | `colo829_amp.bed` |
  | s7 | `mre11a_amp.bed` |
  | s8 | `invsv.bed` |
  | s9 | `dupdiscsv.bed` |

- [ ] Write `harness/slices.tsv`, tab-separated. Columns: `id ref bam bed sample min_variants extra`.
  - Take the reference, BAM(s) and flags from `task-w46-report.md:61-72`, `w4-sweep/NOTES.md` and the matching corpus case file.
  - Find the somatic pairs for s3/s4/s5 via `ls VarDictJava/testdata/integrationtestcases | grep '^Somatic'`.
  - Flags: s1a/s1b/tiny/s4/s8/s9 use `-f 0.001`; s5 uses `-f 0.001 -X 3`; s6 uses `-a 10:0.95` (the documented no-`-D` substitution); s7 uses `-a 10:0.95 -p`; s3 has none.
  - The sample is `SAMPLE`, or `T|N` names for somatic.
  - `min_variants` starts at 1.
- [ ] `harness/run_slice.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: run_slice.sh <id>  — one harness/slices.tsv row through compare.sh. Exit = compare.sh exit.
  set -euo pipefail
  ID=$1
  row=$(awk -F'\t' -v id="$ID" '$1==id' harness/slices.tsv)
  [ -n "$row" ] || { echo "run_slice: unknown id $ID" >&2; exit 2; }
  IFS=$'\t' read -r _ REF BAM BED SAMPLE MINV EXTRA <<< "$row"
  read -ra EX <<< "${EXTRA:-}"
  KEEP_DIR="${SLICE_KEEP:-scratch/slices}/$ID" MIN_VARIANTS="$MINV" \
    bash harness/compare.sh "$REF" "$BAM" "$BED" "$SAMPLE" "${EX[@]}"
  ```
- [ ] `harness/rgate.sh`:
  ```bash
  #!/usr/bin/env bash
  # Regression gate (plan rev2 §R-GATE). Exit 0 = PASS. Never edit this to make a change pass.
  set -uo pipefail
  OUT=${RGATE_DIR:-scratch/rgate}; mkdir -p "$OUT"; fail=0
  echo "== cargo fmt/clippy/test"
  { cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings \
      && cargo test --workspace; } > "$OUT/cargo.txt" 2>&1 || { echo "FAIL cargo ($OUT/cargo.txt)"; fail=1; }
  echo "== integration corpus"
  bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1
  if grep -qx 'PASS: 91  FAIL: 1  SKIP: 0' "$OUT/integration.txt" \
     && [ "$(grep -c '^MISMATCH: ' "$OUT/integration.txt")" -eq 1 ] \
     && grep -qxF 'MISMATCH: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt' "$OUT/integration.txt"
  then :; else echo "FAIL integration ($OUT/integration.txt)"; fail=1; fi
  echo "== Tier-2 micro"
  MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam \
    fixtures/micro/micro.bed > "$OUT/micro.txt" 2>&1 || { echo "FAIL micro ($OUT/micro.txt)"; fail=1; }
  for id in ${RGATE_SLICES:-s1b s5 s7}; do
    echo "== slice $id"
    bash harness/run_slice.sh "$id" > "$OUT/$id.txt" 2>&1 || { echo "FAIL $id ($OUT/$id.txt)"; fail=1; }
  done
  [ "$fail" -eq 0 ] && echo "RGATE PASS" || echo "RGATE FAIL"
  exit "$fail"
  ```
- [ ] **Baseline** (the current `compare.sh` already exits non-zero on mismatch):
  - Run `mkdir -p scratch/baseline; for id in s1a s1b tiny s3 s4 s5 s6 s7 s8 s9; do bash harness/run_slice.sh $id > scratch/baseline/$id.txt 2>&1; echo "$id rc=$?"; done`.
  - Expected: every slice rc=0 except **s1a**, whose only diff should be java `('20','20202875','TA','T')` vs rust `('20','20202875','T','A')`. Counts are expected to match DIVERGENCES W4.6b (312/4/16/14/61/4/96/20/71).
  - Write each observed record count into `min_variants`. For s1a, use the java count.
  - Any other deviation → STOP.
- [ ] `bash harness/rgate.sh` → `RGATE PASS`.
- [ ] Commit `T02: committed W4.6 slice set, run_slice.sh, rgate.sh (baseline 9/10 MATCH; s1a = documented AF-tie)`.

### T3 — Tighten `compare_vcf.py` and `compare_var.py` to ε = 1e-6 (D-05, D-08, D-52)
- [ ] **Red** `harness/tests/test_compare_vcf.py`:
  ```python
  import subprocess, sys, tempfile, os, unittest
  HDR = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
  def rec(pos, alt="T", qual="30", info="DP=10;AF=0.5"):
      return f"chr1\t{pos}\t.\tA\t{alt}\t{qual}\tPASS\t{info}\tGT:DP\t0/1:10\n"
  def run(java, rust, *args):
      d = tempfile.mkdtemp(); j, r = os.path.join(d, "j.vcf"), os.path.join(d, "r.vcf")
      open(j, "w").write(HDR + "".join(java)); open(r, "w").write(HDR + "".join(rust))
      p = subprocess.run([sys.executable, "harness/compare_vcf.py", *args, j, r], capture_output=True, text=True)
      return p.returncode, p.stdout + p.stderr
  class T(unittest.TestCase):
      def test_identical(self):
          rc, out = run([rec(1), rec(2)], [rec(2), rec(1)]); self.assertEqual(rc, 0); self.assertIn("MATCH (2 variants)", out)
      def test_duplicate_key_counts(self):
          rc, out = run([rec(1)], [rec(1), rec(1)]); self.assertEqual(rc, 1); self.assertIn("record counts differ", out)
      def test_min_variants_vacuous(self):
          self.assertEqual(run([], [], "--min-variants", "1")[0], 1)
      def test_vacuous_warns(self):
          rc, out = run([], []); self.assertEqual(rc, 0); self.assertIn("vacuous", out)
      def test_eps(self):
          # Under raw float comparison with eps=1e-6 (A5, D-52):
          self.assertEqual(run([rec(1, qual="1.0000004")], [rec(1, qual="1.0000000")])[0], 0)
          self.assertEqual(run([rec(1, qual="1.000001")], [rec(1, qual="1.0")])[0], 0)
          self.assertEqual(run([rec(1, qual="1.00001")], [rec(1, qual="1.0")])[0], 1)
      def test_set_difference(self):
          self.assertEqual(run([rec(1)], [rec(2)])[0], 1)
      def test_malformed_row(self):
          rc, out = run(["chr1\t1\t.\tA\n"], ["chr1\t1\t.\tA\n"]); self.assertEqual(rc, 1); self.assertIn("malformed", out)
  if __name__ == "__main__": unittest.main()
  ```
  Running `python3 -m unittest harness/tests/test_compare_vcf.py` → fails on duplicate, min_variants, vacuous, malformed, and test_eps.
- [ ] **Implement.** Drop `r5` pre-rounding. Replace `eq()` in `harness/compare_vcf.py` and `harness/compare_var.py` with raw float comparison within `(eps + 1e-12)`:
  ```python
  import math
  def eq(a, b, eps):
      if a == b: return True
      try:
          fa, fb = float(a), float(b)
          if math.isnan(fa) or math.isnan(fb): return False
          return abs(fa - fb) <= (eps + 1e-12)
      except (ValueError, TypeError):
          return False
  ```
  Change default `--eps` to `1e-6` in both `compare_vcf.py` and `compare_var.py`.
  Replace `load`/`main` in `compare_vcf.py` with:
  ```python
  from collections import defaultdict

  def load(path):
      out = []
      with open(path) as fh:
          for n, ln in enumerate(fh, 1):
              if ln.strip() and not ln.startswith("#"):
                  cols = ln.rstrip("\n").split("\t")
                  if len(cols) < 8:
                      print(f"MISMATCH: malformed VCF row {path}:{n}: {ln.strip()[:80]!r}"); sys.exit(1)
                  out.append(cols)
      return out

  def group(rows):
      g = defaultdict(list)
      for r in rows: g[keyof(r)].append(r)
      for k in g: g[k].sort(key="\t".join)
      return g

  def compare_pair(k, j, r, eps, errs):
      w = f"{k[0]}:{k[1]} {k[2]}>{k[3]}"
      if j[6] != r[6]: errs.append(f"{w} FILTER: {j[6]!r} != {r[6]!r}")
      if not eq(j[5], r[5], eps): errs.append(f"{w} QUAL: {j[5]!r} != {r[5]!r}")
      ij, ir = kv(j[7]), kv(r[7])
      if set(ij) != set(ir):
          errs.append(f"{w} INFO keys: only-java={sorted(set(ij)-set(ir))} only-rust={sorted(set(ir)-set(ij))}")
      else:
          for key in ij:
              if ij[key] is True or ir[key] is True:
                  if ij[key] != ir[key]: errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
              elif not eq(str(ij[key]), str(ir[key]), eps):
                  errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
      if len(j) > 9 and len(r) > 9:
          if j[8] != r[8]: errs.append(f"{w} FORMAT keys: {j[8]!r} != {r[8]!r}")
          else:
              fkeys = j[8].split(":")
              for si in range(9, min(len(j), len(r))):
                  for key, xj, xr in zip(fkeys, j[si].split(":"), r[si].split(":")):
                      if not eq(xj, xr, eps): errs.append(f"{w} sample{si-9} {key}: {xj!r} != {xr!r}")

  def main():
      ap = argparse.ArgumentParser()
      ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
      ap.add_argument("--min-variants", type=int, default=0)
      a = ap.parse_args()
      LJ, LR = load(a.java), load(a.rust)
      J, R = group(LJ), group(LR)
      if set(J) != set(R):
          oj = sorted(map(str, set(J) - set(R)))[:10]; orr = sorted(map(str, set(R) - set(J)))[:10]
          print(f"MISMATCH: variant sets differ. java-only(<=10)={oj} rust-only(<=10)={orr}"); sys.exit(1)
      bad = sorted(k for k in J if len(J[k]) != len(R[k]))
      if bad:
          print(f"MISMATCH: record counts differ for {len(bad)} key(s) (<=10): "
                f"{[(str(k), len(J[k]), len(R[k])) for k in bad[:10]]}"); sys.exit(1)
      if len(LJ) < a.min_variants:
          print(f"MISMATCH: only {len(LJ)} records compared (< --min-variants {a.min_variants})"); sys.exit(1)
      errs = []
      for k in sorted(J):
          for j, r in zip(J[k], R[k]): compare_pair(k, j, r, a.eps, errs)
          if len(errs) >= 20: break
      if errs:
          print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
      if not LJ: print("WARNING: vacuous comparison (0 records)", file=sys.stderr)
      print(f"MATCH ({len(LJ)} variants)"); sys.exit(0)
  ```
- [ ] **Green:** unit tests pass and `bash harness/rgate.sh` → PASS.
- [ ] Commit `T03: compare_vcf and compare_var tightened to eps=1e-6 raw floats, multiset check, min-variants`.

### T4 — Fix `run_oracle.sh` (D-06, D-07, D-44, D-50)
- [ ] **Red** `harness/tests/test_run_oracle.sh`. It uses PATH stubs only; the real conda is never touched:
  ```bash
  #!/usr/bin/env bash
  set -uo pipefail
  T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
  printf '#!/usr/bin/env bash\nprintf "%%s\\n" "$@" > "$STUB_LOG"; shift 3; exec "$@"\n' > "$T/conda"
  printf '#!/usr/bin/env bash\n[ -n "${FAKE_VJ_FAIL:-}" ] && exit 7; echo fake\n' > "$T/vardict-java"
  for s in teststrandbias.R var2vcf_valid.pl testsomatic.R var2vcf_paired.pl; do printf '#!/usr/bin/env bash\ncat\n' > "$T/$s"; done
  chmod +x "$T"/*; export PATH="$T:$PATH" STUB_LOG="$T/args"
  bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o.vcf" S -th 8 -UN -mfreq 0.2 -f0.001 >/dev/null
  [ "$(grep -A1 -x -- '-th' "$T/args" | tail -1)" = 8 ] || { echo "FAIL: -th 8 not preserved"; fail=1; }
  grep -qx -- '-UN' "$T/args"    || { echo "FAIL: -UN split"; fail=1; }
  grep -qx -- '-mfreq' "$T/args" || { echo "FAIL: -mfreq split"; fail=1; }
  grep -qx -- '0.001' "$T/args"  || { echo "FAIL: -f0.001 not split"; fail=1; }
  if FAKE_VJ_FAIL=1 bash harness/run_oracle.sh r.fa b.bam x.bed "$T/o.vcf" S >/dev/null 2>&1; then
    echo "FAIL: java crash masked (exit 0)"; fail=1; fi
  [ $fail -eq 0 ] && echo "PASS test_run_oracle"; exit $fail
  ```
  Run it → FAIL on `-th`, `-UN`, `-mfreq` and the masked crash.
- [ ] **Implement.** Make this the first statement of the splitter's `-[a-zA-Z]*)` branch:
  ```bash
        case "${arg:1}" in
          adaptor|chimeric|deldupvar|fisher|DP|mfreq|nmfreq|th|UN|VS)
            # Java multi-char single-dash options (CmdParser.java buildOptions) — never split.
            EXTRA_ARGS+=("$arg"); continue ;;
        esac
  ```
  Insert `set -euo pipefail` as the first line inside **both** `bash -c '…'` scripts.
  In addition, when `${1}.fai` exists, generate a temp contigs BED `awk '{print $1"\t0\t"$2}' "${1}.fai" > "$tmp_contigs"` and pass `-G "$1" -b "$tmp_contigs"` to `var2vcf_valid.pl` and `var2vcf_paired.pl`.
- [ ] **Green:** the test reports PASS; `bash harness/rgate.sh` → PASS.
- [ ] Commit `T04: run_oracle preserves Java multi-char options (-th/-UN/...), pipefail, var2vcf -G/-b`.

### T5 — Rewrite `compare.sh` (D-01, D-02, D-05, D-09, D-40)
- [ ] **Red** `harness/tests/test_compare_sh.sh`. The oracle is real; a fake `cargo` drops our last record:
  ```bash
  #!/usr/bin/env bash
  set -uo pipefail
  cargo build -q --release -p vardict-cli
  T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
  M="fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed"
  mkdir -p "$T/bin"; cat > "$T/bin/cargo" <<EOF
  #!/usr/bin/env bash
  while [ "\$1" != "--" ]; do shift; done; shift
  "$PWD/target/release/vardictrust" "\$@" | head -n -1
  EOF
  chmod +x "$T/bin/cargo"
  out=$(PATH="$T/bin:$PATH" KEEP_DIR="$T/k1" bash harness/compare.sh $M 2>&1); rc=$?
  [ $rc -ne 0 ] || { echo "FAIL1 rc=0"; fail=1; }
  grep -qx MISMATCH <<<"$out" || { echo "FAIL1 no MISMATCH line"; fail=1; }
  grep -q '^[<>]' <<<"$out" || { echo "FAIL1 no diff"; fail=1; }
  [ -s "$T/k1/java.vcf" ] && [ -f "$T/k1/rust.vcf" ] || { echo "FAIL1 files not kept"; fail=1; }
  [ "$(grep -vc '^#' "$T/k1/java.vcf")" -ge 1 ] || { echo "FAIL1 fixture vacuous"; fail=1; }
  KEEP_DIR="$T/k2" MIN_VARIANTS=1 bash harness/compare.sh $M > "$T/o2" 2>&1 || { echo "FAIL2 rc"; fail=1; }
  [ "$(tail -1 "$T/o2")" = MATCH ] && [ -f "$T/k2/java.time" ] || { echo "FAIL2 output"; fail=1; }
  MIN_VARIANTS=100000 bash harness/compare.sh $M >/dev/null 2>&1 && { echo "FAIL3 floor ignored"; fail=1; }
  mkdir -p "$T/tmp"; PATH="$T/bin:$PATH" TMPDIR="$T/tmp" bash harness/compare.sh $M >/dev/null 2>&1
  [ -z "$(ls -A "$T/tmp")" ] || { echo "FAIL4 temp dir leaked"; fail=1; }
  [ $fail -eq 0 ] && echo "PASS test_compare_sh"; exit $fail
  ```
  Run it → FAIL1, FAIL2, FAIL4.
- [ ] **Implement:** overwrite `harness/compare.sh` with:
  ```bash
  #!/usr/bin/env bash
  set -euo pipefail
  # Usage: compare.sh <ref> <bam|tumor|normal> <bed> [sample] [extra args...]
  # Env: AF_THR (default 0.01) | KEEP_DIR (persist java.vcf rust.vcf compare.txt *.time)
  #      MIN_VARIANTS (default 0) | CMP_TIMEOUT per engine (default 12h) | CMP_STRICT=1 (report-only, T11)
  # Exit 0 = MATCH; non-zero = MISMATCH or engine failure. Decide on the EXIT CODE, never on text.
  REF=$1; BAM=$2; BED=$3; SAMPLE=${4:-SAMPLE}; shift 4 || shift $#
  if [ -n "${KEEP_DIR:-}" ]; then mkdir -p "$KEEP_DIR"; OUT=$KEEP_DIR
  else OUT=$(mktemp -d); trap 'rm -rf "$OUT"' EXIT; fi
  TO=${CMP_TIMEOUT:-12h}
  # Mirror upstream IntegrationTest DEFAULT_ARGS ("-z -c 1 -S 2 -E 3 -g 4"): the
  # corpus and VarDictJava's own suite always run zero-based. run_oracle.sh
  # translates the flag-less -z into Java's explicit "-z 1" form.
  /usr/bin/time -v -o "$OUT/java.time" timeout "$TO" \
    bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$OUT/java.vcf" "$SAMPLE" -z "$@" \
    || { echo "ENGINE-FAIL: oracle"; exit 2; }
  /usr/bin/time -v -o "$OUT/rust.time" timeout "$TO" \
    cargo run -q --release -p vardict-cli -- -G "$REF" -f "${AF_THR:-0.01}" -N "$SAMPLE" \
    -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "$@" "$BED" > "$OUT/rust.vcf" \
    || { echo "ENGINE-FAIL: vardictrust"; exit 2; }
  rc=0
  python3 harness/compare_vcf.py --min-variants "${MIN_VARIANTS:-0}" "$OUT/java.vcf" "$OUT/rust.vcf" \
    | tee "$OUT/compare.txt" || rc=$?
  if [ "$rc" -eq 0 ]; then echo "MATCH"
  else echo "MISMATCH"; diff "$OUT/java.vcf" "$OUT/rust.vcf" | head -50 || true; fi
  exit "$rc"
  ```
- [ ] **Green:** the test reports PASS; `bash harness/rgate.sh` → PASS.
- [ ] Commit `T05: compare.sh KEEP_DIR/MIN_VARIANTS/timeouts, MISMATCH+diff on failure, trap cleanup`.

### T5a — VCF header parity (D-50)
- [ ] **Red** `harness/tests/test_vcf_header.py`:
  - Run `compare.sh` on `fixtures/micro`.
  - Parse headers (all `#` lines). Check `##fileformat`, `##source`, `##reference`, `##contig`, `##INFO`, `##FILTER`, `##FORMAT`, `#CHROM`.
  - Assert that all INFO/FORMAT/FILTER keys and their descriptions match Java's Perl script verbatim (`var2vcf_valid.pl:48-120`).
  - Allowed volatile line exceptions: `##fileDate` (generation timestamp) and `##source` (version tag).
  - In existing code, Rust `PMEAN` description is `"Mean position in reads"` whereas Java is `"The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call"`. Test FAILS.
- [ ] **Implement:**
  - In `crates/vardict-core/src/output/vcf.rs`: update all INFO/FILTER/FORMAT descriptions to match Java `var2vcf_valid.pl` and `var2vcf_paired.pl` exactly.
  - In `harness/compare_vcf.py`: add header comparison function `compare_headers(j_hdr, r_hdr, errs)` checking ordered header lines, ignoring only lines matching `^##(fileDate|source)=`.
- [ ] **Green:** `python3 -m unittest harness/tests/test_vcf_header.py` passes; `bash harness/rgate.sh` → PASS.
- [ ] Commit `T05a: VCF header parity: match Java INFO/FILTER/FORMAT descriptions, validate headers in compare_vcf`.

### T5b — Edge-case parity (D-51)
- [ ] **Red** `harness/tests/test_edge_cases.sh`:
  - Creates 4 edge-case fixtures under `scratch/test_edge_cases/`:
    1. `empty.bam`: valid SAM/BAM header, 0 alignments.
    2. `zero_cov.bam`: reads on chr1; query `-R chr2:1-100` (zero coverage in window).
    3. `single_read.bam`: exactly 1 read covering target.
    4. `contig_boundary.bam`: reads at contig start (`chr1:1-50`) and contig end.
  - Runs `compare.sh` for both engines on each case.
  - Asserts both exit 0 and produce matching output (empty variants table or matching headers).
  - Asserts no panics or crashes in `vardictrust`.
- [ ] **Implement:**
  - If any case causes a panic or diff, fix the edge condition in `vardict-core` (`cigar/parse.rs`, `bam.rs`, `output/vcf.rs`) with safe `Option`/bounds checks (`#![forbid(unsafe_code)]`).
- [ ] **Green:** `bash harness/tests/test_edge_cases.sh` PASS; `bash harness/rgate.sh` → PASS.
- [ ] Commit `T05b: edge-case parity: empty BAM, zero coverage, single read, contig boundary`.

### T6 — Evidence fingerprint, recorder, PLAN §15 acceptance (D-03, D-11, D-37, D-47)
- [ ] `harness/engine_fingerprint.sh`:
  ```bash
  #!/usr/bin/env bash
  # Content hash of everything that can change engine output or the Tier-2 gate.
  set -euo pipefail
  git diff --quiet HEAD -- crates harness Cargo.lock \
    || { echo "engine_fingerprint: uncommitted changes under crates/ harness/ Cargo.lock" >&2; exit 1; }
  { git rev-parse HEAD:crates/vardict-core/src HEAD:crates/vardict-cli/src \
      HEAD:harness/compare.sh HEAD:harness/run_oracle.sh HEAD:harness/compare_vcf.py
    cargo metadata --locked --format-version 1 | python3 -c 'import json,sys
  m=json.load(sys.stdin); print("\n".join(sorted(p["name"]+"@"+p["version"] for p in m["packages"] if p.get("source"))))'
  } | sha256sum | cut -c1-16
  ```
- [ ] `harness/record_evidence.py`:
  ```python
  #!/usr/bin/env python3
  """Usage: record_evidence.py <name> <result> <dir|-> -- <cmd...>
  Writes benchmarks/evidence/<name>.json; reads optional <dir>/compare.txt, java.time, rust.time, extra.json."""
  import datetime, json, pathlib, re, subprocess, sys
  name, result, d = sys.argv[1], sys.argv[2], sys.argv[3]
  cmd = sys.argv[sys.argv.index("--") + 1:]
  keep = None if d == "-" else pathlib.Path(d)
  def tv(p):
      if keep is None or not (keep / p).exists(): return None
      t = (keep / p).read_text()
      g = lambda k: (re.search(rf"{re.escape(k)}: (.+)", t) or [None, None])[1]
      return {"wall": g("Elapsed (wall clock) time (h:mm:ss or m:ss)"),
              "max_rss_kb": g("Maximum resident set size (kbytes)"), "exit": g("Exit status")}
  cmp_txt = (keep / "compare.txt").read_text() if keep and (keep / "compare.txt").exists() else ""
  m = re.search(r"^MATCH \((\d+) variants\)", cmp_txt, re.M)
  rec = {"name": name, "result": result, "n_variants": int(m.group(1)) if m else None,
         "fingerprint": subprocess.check_output(["bash", "harness/engine_fingerprint.sh"], text=True).strip(),
         "head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
         "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
         "cmd": cmd, "java_time": tv("java.time"), "rust_time": tv("rust.time"),
         "compare_head": cmp_txt.splitlines()[:25]}
  if keep and (keep / "extra.json").exists(): rec["extra"] = json.loads((keep / "extra.json").read_text())
  out = pathlib.Path("benchmarks/evidence") / f"{name}.json"
  out.parent.mkdir(parents=True, exist_ok=True)
  out.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n"); print(out)
  ```
- [ ] Create `docs/release/waivers.md` with the single header line `# v1.0.0 acceptance waivers (user-approved only)`. Lines take the form `- item <n>: <reason> — approved by <user> <date>`, and **only the controller writes them, after explicit user approval**.
- [ ] `harness/acceptance.sh`:
  ```bash
  #!/usr/bin/env bash
  # PLAN.md §15 acceptance, items 1-8 verbatim. Exit 0 only if every item is PASS or user-WAIVED.
  set -uo pipefail
  OUT=${ACCEPT_DIR:-scratch/acceptance}; mkdir -p "$OUT"; : > "$OUT/log.txt"; declare -A R
  FP=$(bash harness/engine_fingerprint.sh) || { echo "acceptance: commit first"; exit 2; }
  waived() { grep -qE "^- item $1:" docs/release/waivers.md 2>/dev/null; }
  mark() { if [ "$2" = 1 ]; then R[$1]=PASS; elif waived "$1"; then R[$1]=WAIVED; else R[$1]="FAIL: $3"; fi; }
  ev() { python3 - "benchmarks/evidence/$1.json" "$FP" <<'PY'
  import json, sys
  p, fp = sys.argv[1], sys.argv[2]
  try: d = json.load(open(p))
  except Exception as e: sys.exit(f"{p}: {e}")
  if d.get("fingerprint") != fp: sys.exit(f"{p}: stale fingerprint {d.get('fingerprint')} != {fp}")
  if d.get("result") not in ("MATCH", "EQUAL", "PASS"): sys.exit(f"{p}: result={d.get('result')}")
  PY
  }
  run() { "$@" >> "$OUT/log.txt" 2>&1; }
  # 1 Unit
  t=1; run cargo test --workspace || t=0; mark 1 $t "cargo test"
  # 2 Tier-1: corpus .var parity vs VarDictJava expected + micro .var golden
  bash harness/run_integration.sh > "$OUT/integration.txt" 2>&1
  t=0; grep -qx 'PASS: 91  FAIL: 1  SKIP: 0' "$OUT/integration.txt" \
    && [ "$(grep -c '^MISMATCH: ' "$OUT/integration.txt")" -eq 1 ] \
    && grep -qxF 'MISMATCH: Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt' "$OUT/integration.txt" \
    && run cargo build -q --release && run bash harness/smoke.sh && t=1
  mark 2 $t "integration/smoke"
  # 3 Tier-2: micro + every committed slice + external datasets
  t=1; MIN_VARIANTS=1 run bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed || t=0
  for id in $(awk -F'\t' '!/^#/{print $1}' harness/slices.tsv); do run bash harness/run_slice.sh "$id" || { t=0; echo "slice $id failed" >> "$OUT/log.txt"; }; done
  for e in hg008-tier2 seqc2-slices-tier2 seqc2-exome-tier2; do run ev "$e" || t=0; done
  mark 3 $t "Tier-2"
  # 4 Robustness: crash corpus (#341 family, PR #402, #187) + CRAM twin
  t=1; run cargo test -p vardict-cli --test robustness || t=0; run ev hg008-cram-twin || t=0; mark 4 $t "robustness/CRAM"
  # 5 Scale (thresholds pre-registered in T21)
  t=1; run ev scale || t=0; mark 5 $t "scale"
  # 6 Truth set: rtg exact java==rust equality on both datasets (som.py/truvari via waiver)
  t=1; for e in hg008-truth seqc2-truth; do run ev "$e" || t=0; done; mark 6 $t "truth"
  # 7 Hygiene
  t=1; run cargo fmt --all -- --check || t=0; run cargo clippy --all-targets --all-features -- -D warnings || t=0
  run cargo audit || t=0; run cargo deny check || t=0; mark 7 $t "hygiene"
  # 8 Distribution: musl CI evidence + --version content
  t=1; v=$(target/release/vardictrust --version 2>&1)
  for s in "VarDictJava v1.8.3" 0886261 009e017 Lai gkw227 AstraZeneca-NGS "AI disclosure"; do grep -qF "$s" <<<"$v" || t=0; done
  run ev musl-ci || t=0; mark 8 $t "distribution"
  rc=0; for i in 1 2 3 4 5 6 7 8; do echo "item $i: ${R[$i]}"; [[ "${R[$i]}" == FAIL* ]] && rc=1; done
  echo "fingerprint $FP (log: $OUT/log.txt)"; exit $rc
  ```
- [ ] **Red→green** `harness/tests/test_evidence.sh`:
  - A fake evidence file whose fingerprint is `deadbeef` must be rejected by the `ev` logic. Extract `ev` into a sourced snippet or run the python heredoc directly.
  - `record_evidence.py t-test MATCH - -- true` must write the current fingerprint.
  - Clean up the test JSON afterwards.
- [ ] Baseline run of `bash harness/acceptance.sh`:
  - Install `cargo-audit` and `cargo-deny` if missing (`cargo install --locked …`, user-level).
  - Expected: items 1, 2, 7 PASS; items 3–6 and 8 FAIL with missing evidence.
  - Record this honest baseline in the report.
- [ ] Commit `T06: engine fingerprint, evidence recorder, PLAN §15 acceptance.sh, waivers ledger`.

### T7 — Policy amendments + ROADMAP v1.1.0 (D-21, D-22, D-30, D-49, D-52, D-55, D-56)
`AGENTS.md` is a symlink to `CLAUDE.md`, so edit it once.
- [ ] CLAUDE §5:
  - Update tolerance from `ε = 1e-5` to `ε = 1e-6` (A5).
  - Append:
    > - **Release policy (user decision 2026-10-07 & 2026-10-08).** v1.0.0 = pure parity with v1.8.3: no output-changing fixes and no fix flags. Corrections ship from v1.1.0 per `ROADMAP.md`, each behind `--fix <id>` (default off) with a diff-confinement gate. The default output stays byte-identical to v1.0; with the flag, diffs are confined to the targeted records.
- [ ] CLAUDE §8, append:
  > - **Java `HashMap` iteration that reaches output** is emulated with `vardict_core::javamap::JavaHashMap`. It ports `java.util.HashMap` spread, bin, resize and treeify-trigger semantics, and is validated against real Java by the op-log tests. This is the only sanctioned non-sorted order. Never use `std::collections::HashMap` there.
- [ ] CLAUDE §11, append: `Post-1.0 = v1.1.0+, governed by ROADMAP.md.`
- [ ] CLAUDE §18: remove the `Co-Authored-By` trailer requirement from the commit template per user directive (A10). AI involvement is disclosed transparently in `README.md` and `CITATION.cff`.
- [ ] PLAN.md edits:
  - Update all `1e-5` references (lines 40, 173, 174, 201, 405, 910) to `1e-6` (A5).
  - §8: add "v1.0 = pure parity (A1); see ROADMAP.md".
  - §15.5: replace CodSpeed with pre-registered performance budget (Rust wall <= Java wall, Rust max RSS <= Java max RSS, A9); update threads to 8 (A3).
  - §15.8 and §17 Distribution: cargo-dist (`dist`) for static musl release + gated OIDC publish (A8).
- [ ] Create `ROADMAP.md`:
  ```markdown
  # VarDictRust roadmap

  v1.0.0 replicates VarDictJava v1.8.3 output exactly, including the upstream bugs below.
  From v1.1.0 each correction is opt-in via `--fix <id>` (default off) and lands only with
  (1) a manual reproduction on VarDictJava 1.8.3, (2) a unit test, (3) a `compare_fix.sh` run proving
  default output is unchanged and flagged output differs only on the targeted records.

  ## v1.1.0 candidates — unambiguous bugs with a clear correct behaviour
  | id | upstream | v1.8.3 behaviour (replicated in v1.0) | proposed `--fix` behaviour |
  |---|---|---|---|
  | 402 | PR #402 | modifyCigar yields `…D-1M…` for reads like `104M1I2M2D1M1D5M1D2M9D27M`; Java throws, read skipped | repair to `104M14D11I9D27M` per PR #402 |
  | 292 | #292 | SHIFT3 can give end < start with `-3` | per issue analysis |
  | 342 | #342 | `<dup-N>` in ALT (invalid VCF) | symbolic `<DUP>` |
  | 369 | #369 | calls emitted outside `-R` | restrict to region |
  | 378 | #378 | paired + amplicon column shift | correct columns (or loud Unsupported) |
  | ub1 | DIVERGENCES UB-1 | position-0 row in `.var` | drop row |

  ## Needs a specification first (debatable heuristics; no flag until specified)
  #92/#386/#323/#316/#400 (VD > DP depth) · #62 (over-filtering) · #363 family (realignment drops) ·
  #381 · #282 (soft-clip) · #264 · HICOV excludes insertions (ToVarsBuilder.java:523-538).

  ## Open upstream PRs reviewed (2026-10-08)
  #402 → row 402. #371, #372, #335 → no output impact on VarDictRust (N/A).

  ## Post-1.0 features
  UMI support, long reads — separate design docs.
  ```
- [ ] Commit `T07: policy — v1.0 pure parity, JavaHashMap sanction, ROADMAP v1.1.0, host-accurate PLAN §15`.

### T8 — `--version` credit + README output fix (D-19, D-20)
- [ ] **Red** `crates/vardict-cli/tests/version.rs`:
  ```rust
  #[test]
  fn version_credits_original_and_discloses_ai() {
      let out = std::process::Command::new(env!("CARGO_BIN_EXE_vardictrust"))
          .arg("--version").output().expect("run vardictrust");
      assert!(out.status.success());
      let s = String::from_utf8_lossy(&out.stdout);
      for needle in ["VarDictJava v1.8.3", "0886261", "009e017", "Lai", "gkw227", "AstraZeneca-NGS", "AI disclosure", "MIT"] {
          assert!(s.contains(needle), "--version lacks {needle:?}:\n{s}");
      }
  }
  ```
- [ ] Add to `#[command(...)]` in `cli.rs`:
  ```rust
      long_version = concat!(
          env!("CARGO_PKG_VERSION"), "\n",
          "Idiomatic Rust rewrite of VarDictJava v1.8.3 (commit 0886261; Perl/R submodule 009e017).\n",
          "Output validated against VarDictJava v1.8.3 - see README \"Validation status\".\n",
          "Original: VarDict/VarDictJava, AstraZeneca-NGS - https://github.com/AstraZeneca-NGS/VarDictJava\n",
          "Please cite: Lai Z, et al. VarDict: a novel and versatile variant caller for next-generation\n",
          "  sequencing in cancer research. Nucleic Acids Res. 2016;44(11):e108. doi:10.1093/nar/gkw227\n",
          "AI disclosure: written with AI coding assistants under human direction; see README \"AI Disclosure\".\n",
          "License: MIT"
      ),
  ```
- [ ] `README.md:26`: change `> output.var` to `> output.vcf`.
- [ ] **Green:** the test passes; R-GATE PASS.
- [ ] Commit `T08: --version credits VarDictJava 1.8.3/Lai 2016/AstraZeneca + AI disclosure; README somatic output is VCF`.

### T9 — PR #402 family: emulate Java's skip (D-15, D-45) — `pro`
- [ ] **Read & restate** in the report:
  - `CigarModifier.java` where `rm` is computed (Rust: `modify.rs:1102-1105`).
  - `CigarParser.java:80-95` (try/catch).
  - `Utils.java:261-277` and `AbstractMode.java:64-68`.
  - Answer: (a) which statement throws for a negative-length element; (b) which per-read side effects (coverage, `maxReadLength`, duplicate bookkeeping, …) have already happened before that throw.
- [ ] **Fixture.** Start with `gh pr view 402 -R AstraZeneca-NGS/VarDictJava --json title,body,files` and `gh pr diff 402 -R AstraZeneca-NGS/VarDictJava`. If the PR carries a read or reference, use it verbatim. Otherwise:
  - Build a 600 bp `fixtures/pr402/ref.fa`.
  - Build `pr402.sam` with 30 normal reads carrying one real SNV, plus the read with CIGAR `104M1I2M2D1M1D5M1D2M9D27M`. Place and sequence it so that the `rm<0` path actually runs; prove this in the next step.
  - `samtools view -b | sort | index` → `pr402.bam`, plus a 4-col `pr402.bed`.
  - Keep the fixture KB-size.
- [ ] **Red unit test** (`modify.rs` `#[cfg(test)]`):
  - Feeding the PR read to the modify entry point must return the new throw outcome (`Err(VardictError::…)` or `ModifyOutcome::JavaThrows`).
  - The outcome must contain no negative lengths and no partially-updated read.
  - Today it FAILS.
- [ ] **Oracle truth:**
  - Run `KEEP_DIR=scratch/pr402 MIN_VARIANTS=1 bash harness/compare.sh fixtures/pr402/ref.fa fixtures/pr402/pr402.bam fixtures/pr402/pr402.bed`.
  - Save the Java stderr exception.
  - Record the pre-fix result and diff.
- [ ] **Implement:**
  - When the modified CIGAR would have a non-positive length, return the throw outcome.
  - `parse.rs` skips the read at the point of Java's catch:
    - after replicating exactly the side effects found in step 1;
    - before overwriting `start_pos`, `query_seq` or `qual`.
  - Count the exception, but **never abort**. Cite `Utils.java:261-277` and CLAUDE §5.
  - Add `/// ports CigarModifier.java:<a-b>` and `/// ports CigarParser.java:80-95` comments.
- [ ] **Siblings:** add unit tests that drive the Rust counterparts of `twoDeletionsInsertionToComplex` and `threeIndels` with adversarial short-M cases. Apply the same outcome to any non-positive length; otherwise assert there is none.
- [ ] **Green:**
  - Unit tests pass.
  - The fixture `compare.sh` exits 0.
  - New `robustness.rs` case: 11 copies of the bad read → exit 0 and no panic. The output must equal the 1-copy output if step 1 found no pre-throw side effects; otherwise it must match the expectation derived from Java's semantics.
  - R-GATE PASS.
- [ ] DIVERGENCES.md: "PR #402 family — Java per-read skip replicated; >10 such reads: Java exits, VarDictRust continues (robustness, CLAUDE §5)".
- [ ] Commit `T09: emulate Java skip-on-throw for PR #402 modifyCigar family; never abort at 10`.

### T10 — `JavaHashMap` iteration order for s1a (D-16) — `pro`, ⏱ ≤6h, ≤2 fix rounds
- [ ] **Read & restate.** Read `task-p3-report.md`, `ToVarsBuilder.java:93`, `AlignedVarsData.java:14` and `SimplePostProcessModule.java:35`, plus every HashMap feeding `alignedVariants` (`nonInsertionVariants`, `insertionVariants`, …) and its somatic equivalent.
  - For each map, list its constructor and the **exact mutating methods** used. Note where its iteration reaches output.
  - State how Rust tiling relates to Java's per-region scope. The emulated insertion sequence must equal Java's per-region sequence.
  - If any method beyond put/putIfAbsent/computeIfAbsent/remove/get/containsKey/iteration is used, extend the model **and** the oracle program first.
- [ ] **Java oracle program** `crates/vardict-core/testdata/javamap/JavaMapOracle.java`:
  ```java
  import java.util.*;
  public class JavaMapOracle {
    public static void main(String[] a) {
      long seed = Long.parseLong(a[0]); int nOps = Integer.parseInt(a[1]);
      int base = Integer.parseInt(a[2]); int span = Integer.parseInt(a[3]); int stride = Integer.parseInt(a[4]);
      Random r = new Random(seed); HashMap<Integer, Integer> m = new HashMap<>(); StringBuilder o = new StringBuilder();
      for (int i = 0; i < nOps; i++) {
        int k = base + stride * r.nextInt(span); int op = r.nextInt(10);
        if (op < 5) { m.put(k, i); o.append("P ").append(k).append('\n'); }
        else if (op < 8) { m.computeIfAbsent(k, x -> 0); o.append("C ").append(k).append('\n'); }
        else { m.remove(k); o.append("R ").append(k).append('\n'); }
      }
      o.append("ORDER"); for (Integer k : m.keySet()) o.append(' ').append(k);
      System.out.println(o);
    }
  }
  ```
  Generate the op logs with:
  ```bash
  conda run -n vardict-oracle bash -c 'cd crates/vardict-core/testdata/javamap && javac JavaMapOracle.java \
    && java JavaMapOracle 1 3000 20202000 5000 1 > ops_1.txt \
    && java JavaMapOracle 2 20000 20150000 90000 1 > ops_2.txt \
    && java JavaMapOracle 3 400 0 300 1024 > ops_treeify.txt \
    && java -version 2>&1 | head -1 > JDK_VERSION && rm -f JavaMapOracle.class'
  ```
- [ ] **Model** `crates/vardict-core/src/javamap.rs`, registered with `pub mod javamap;` in `lib.rs`:
  ```rust
  //! Deterministic emulation of `java.util.HashMap<Integer, V>` iteration order (OpenJDK 8+).
  //! ports java.util.HashMap: hash() spread, putVal (tail append), computeIfAbsent (head insert),
  //! removeNode, resize() order-preserving lo/hi split, treeifyBin's resize-below-64 trigger.
  //! Treeified bins (>= 9 nodes with table >= 64) are not ordered by this model: `treeified()` reports them.
  const DEFAULT_CAPACITY: usize = 16;
  const TREEIFY_THRESHOLD: usize = 8;
  const MIN_TREEIFY_CAPACITY: usize = 64;

  #[derive(Debug, Clone)]
  pub struct JavaHashMap<V> { table: Vec<Vec<(i32, V)>>, len: usize, treeified: bool }

  impl<V> Default for JavaHashMap<V> { fn default() -> Self { Self::new() } }

  fn spread(key: i32) -> u32 { let h = key as u32; h ^ (h >> 16) }

  impl<V> JavaHashMap<V> {
      pub fn new() -> Self { Self { table: Vec::new(), len: 0, treeified: false } }
      pub fn len(&self) -> usize { self.len }
      pub fn is_empty(&self) -> bool { self.len == 0 }
      pub fn treeified(&self) -> bool { self.treeified }
      fn slot(&self, key: i32) -> Option<usize> {
          let n = self.table.len();
          if n == 0 { None } else { Some((spread(key) as usize) & (n - 1)) }
      }
      fn threshold(&self) -> usize { self.table.len() / 4 * 3 }
      /// `put` / `putIfAbsent` (putVal): replace in place, else append at bin tail.
      pub fn insert(&mut self, key: i32, value: V) -> Option<V> {
          if self.table.is_empty() { self.resize(); }
          let i = self.slot(key)?;
          let bin = self.table.get_mut(i)?;
          if let Some(e) = bin.iter_mut().find(|(k, _)| *k == key) { return Some(std::mem::replace(&mut e.1, value)); }
          bin.push((key, value));
          let bin_len = bin.len();
          if bin_len > TREEIFY_THRESHOLD { self.treeify_bin(); }
          self.len += 1;
          if self.len > self.threshold() { self.resize(); }
          None
      }
      /// `computeIfAbsent`: pre-resize check, new node at bin HEAD, no post-insert resize.
      pub fn compute_if_absent(&mut self, key: i32, f: impl FnOnce() -> V) -> Option<&mut V> {
          if self.table.is_empty() || self.len > self.threshold() { self.resize(); }
          let i = self.slot(key)?;
          if !self.table.get(i)?.iter().any(|(k, _)| *k == key) {
              let bin = self.table.get_mut(i)?;
              let had = bin.len();
              bin.insert(0, (key, f()));
              self.len += 1;
              if had >= TREEIFY_THRESHOLD - 1 { self.treeify_bin(); }
          }
          self.get_mut(key)
      }
      pub fn remove(&mut self, key: i32) -> Option<V> {
          let i = self.slot(key)?;
          let bin = self.table.get_mut(i)?;
          let p = bin.iter().position(|(k, _)| *k == key)?;
          self.len -= 1;
          Some(bin.remove(p).1)
      }
      pub fn get(&self, key: i32) -> Option<&V> {
          self.table.get(self.slot(key)?)?.iter().find(|(k, _)| *k == key).map(|e| &e.1)
      }
      pub fn get_mut(&mut self, key: i32) -> Option<&mut V> {
          let i = self.slot(key)?;
          self.table.get_mut(i)?.iter_mut().find(|(k, _)| *k == key).map(|e| &mut e.1)
      }
      pub fn contains_key(&self, key: i32) -> bool { self.get(key).is_some() }
      /// Java iteration order: bins ascending, each bin in list order.
      pub fn iter(&self) -> impl Iterator<Item = (i32, &V)> { self.table.iter().flatten().map(|(k, v)| (*k, v)) }
      pub fn keys(&self) -> impl Iterator<Item = i32> + '_ { self.iter().map(|(k, _)| k) }
      fn treeify_bin(&mut self) {
          if self.table.len() < MIN_TREEIFY_CAPACITY { self.resize(); } else { self.treeified = true; }
      }
      fn resize(&mut self) {
          let old = std::mem::take(&mut self.table);
          let n = old.len();
          let new_n = if n == 0 { DEFAULT_CAPACITY } else { n.saturating_mul(2) };
          let mut table: Vec<Vec<(i32, V)>> = (0..new_n).map(|_| Vec::new()).collect();
          for (j, bin) in old.into_iter().enumerate() {
              for (k, v) in bin {
                  let dst = if (spread(k) as usize) & n != 0 { j + n } else { j };
                  if let Some(b) = table.get_mut(dst) { b.push((k, v)); }
              }
          }
          self.table = table;
      }
  }
  ```
  Verify `computeIfAbsent`'s head-insert, its pre-resize and its treeify condition against the oracle JDK's `java/util/HashMap.java` (`lib/src.zip` in the conda env), and correct any line that differs. **The op-log test is the arbiter**, not this listing.
- [ ] **Red→green op-log test** in `#[cfg(test)]`:
  - Replay each `ops_*.txt` (`include_str!`) through P → `insert`, C → `compute_if_absent`, R → `remove`.
  - For `ops_1` and `ops_2`, assert `keys()` equals Java's `ORDER`.
  - For `ops_treeify`, assert `treeified()`.
- [ ] **Wire in.** Only at the step-1 sites whose **iteration** reaches output:
  - Replace the BTreeMap with `JavaHashMap`, filled with the same method Java uses and in the same per-region sequence.
  - If `treeified()` ever becomes true, `log::warn!` the region.
  - Leave lookup-only maps unchanged.
- [ ] **Green:**
  - `bash harness/run_slice.sh s1a` → **exit 0**.
  - Set rgate's default to `${RGATE_SLICES:-s1a s1b s5 s7}` (a tightening, same commit); `bash harness/rgate.sh` → PASS.
  - `for id in tiny s3 s4 s6 s8 s9; do bash harness/run_slice.sh $id || echo "FAIL $id"; done` → no FAIL.
  - `cargo test -p vardict-cli --test thread_determinism` passes.
- [ ] **Checkpoint report (A7).** At 6h: report s1a status, current diff, and next steps to the ledger and chat. Keep working until s1a exits 0. A divergence is not shipped in v1.0.0 (A1, A7).
- [ ] DIVERGENCES Cat 14 / s1a → RESOLVED with exact match.
- [ ] Commit `T10: emulate java.util.HashMap iteration order where it reaches output (closes s1a AF-tie)`.

### T11 — Strict-bytes audit + `--fisher` residues + mandatory SBF math (D-17, D-18, D-48, D-53) — `pro`
- [ ] **Red:** add `test_strict_bytes_reports` to `test_compare_vcf.py`. Java has `NM=1.0` and Rust has `NM=1`. Running with `--strict-bytes` must exit 0 and print `STRICT-BYTES: 1`.
- [ ] **Implement (report-only):**
  - Add `--strict-bytes` to `compare_vcf.py`. `compare_pair` takes an optional `strict` list: whenever raw strings differ but `eq()` passed, append `f"{w} {field}: {xj!r} vs {xr!r}"`.
  - Before the MATCH line, print `STRICT-BYTES: <n> textual difference(s) within eps` plus up to 20 examples. The exit code is unchanged.
  - In `compare.sh`, `CMP_STRICT=1` adds `--strict-bytes`.
  - Commit this as `T11.0`, with R-GATE.
- [ ] **Measure:**
  - `CMP_STRICT=1 KEEP_DIR=scratch/strict/micro-fisher MIN_VARIANTS=1 bash harness/compare.sh fixtures/micro/micro.ref.fa fixtures/micro/micro.bam fixtures/micro/micro.bed SAMPLE --fisher`.
  - The same on s1b with `--fisher`.
  - `CMP_STRICT=1` on every slice without `--fisher`.
  - Tabulate field → count → example.
- [ ] **Fix field by field:**
  - Locate the Java print site (`String.valueOf(double)`/`Double.toString`, `DecimalFormat`, `String.format`) and port the exact rule. Cite the oracle JDK from `JDK_VERSION`: JDK ≥19 `Double.toString` is shortest-repr.
  - Add a unit test with literal Java-produced strings.
  - Re-audit until that field's count is 0, with R-GATE.
  - Commit `T11.<k>: match Java formatting for <FIELD>; strict-bytes residue -> 0`.
- [ ] **SBF 1-ulp mandatory fix (D-18, D-53):**
  - Under ε = 1e-6 (A5), differences in the 5th decimal (1e-5) are fatal mismatches.
  - Port R's exact `dhyper` implementation and its `1 + 1e-7` `relErr` summation logic into `crates/vardict-core/src/stats/fisher.rs`.
  - Test against real R / `teststrandbias.R` output to ulp(1) precision.
  - There is no 'accept as residue' exit: keep working until all Fisher/SBF p-values and odds ratios match within ε = 1e-6 (A7). Checkpoint report at 2h if debugging.
- [ ] DIVERGENCES: SBF ulp residue → RESOLVED.
- [ ] Commit `T11.final: port R dhyper and relErr summation for exact SBF/Fisher p-values under eps=1e-6`.

### T12 — Tools env, courtesy wrapper, concordance (D-27, D-28, D-29)
- [ ] `harness/env/g4-tools.yml`:
  ```yaml
  name: g4-tools
  channels: [conda-forge, bioconda]
  dependencies: [samtools>=1.20, bcftools>=1.20, htslib>=1.20, bedtools>=2.31, rtg-tools>=3.12, bwa>=0.7.18, gitleaks>=8.18, python>=3.11]
  ```
  Run `conda env create -f harness/env/g4-tools.yml`. Then **pin versions** (rewrites.bio), with each file committed:
  - `conda env export -n vardict-oracle --no-builds | grep -v '^prefix:' > benchmarks/env/vardict-oracle.yml`
  - `conda env export -n g4-tools --no-builds | grep -v '^prefix:' > benchmarks/env/g4-tools.yml`
  - `{ rustc -Vv; cargo -V; conda run -n vardict-oracle java -version 2>&1; } > benchmarks/env/toolchain.txt`
  - The test `harness/tests/test_env_pinned.sh` asserts that `vardict-java` 1.8.3 appears in `benchmarks/env/vardict-oracle.yml`.
- [ ] `concordance.sh` uses `--squash-ploidy`. The rtg manual recommends it for somatic evaluation: it matches alleles regardless of genotype or zygosity, which is right for a sites-only baseline vs. tumor calls. vcfeval reads PASS records only by default, which is why the truth FILTERs are stripped.
- [ ] `harness/courtesy.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: EST_MIN=<minutes> courtesy.sh <cmd...>   — decision A3: never touch services.
  set -euo pipefail
  : "${EST_MIN:?courtesy: set EST_MIN (estimated runtime, minutes)}"
  free=$(df -BG --output=avail "${COURTESY_DIR:-.}" | tail -1 | tr -dc 0-9)
  [ "$free" -ge "${MIN_FREE_GB:-150}" ] || { echo "courtesy: ${free}G free < ${MIN_FREE_GB:-150}G" >&2; exit 3; }
  avail=$(( $(awk '/MemAvailable/{print $2}' /proc/meminfo) / 1048576 ))
  [ "$avail" -ge "${MIN_AVAIL_GB:-12}" ] || { echo "courtesy: ${avail}G RAM available < ${MIN_AVAIL_GB:-12}G" >&2; exit 5; }
  hhmm=${COURTESY_NOW:-$(date +%H%M)}
  now=$((10#${hhmm:0:2} * 60 + 10#${hhmm:2:2})); end=$((now + EST_MIN)); ws=220; we=315
  for off in 0 1440; do
    if [ "$now" -lt $((we + off)) ] && [ "$end" -gt $((ws + off)) ]; then
      echo "courtesy: [$hhmm +${EST_MIN}m] overlaps 03:40-05:15 maintenance; start after 05:15 or split" >&2; exit 4
    fi
  done
  exec nice -n 19 ionice -c 3 "$@"
  ```
- [ ] `harness/tests/test_courtesy.sh`:

  | Case | Expected exit |
  |---|---|
  | `COURTESY_NOW=0300 EST_MIN=60` | 4 |
  | `COURTESY_NOW=2300 EST_MIN=300` | 4 |
  | `COURTESY_NOW=0600 EST_MIN=60 … true` | 0 |
  | `MIN_FREE_GB=999999` | 3 |
  | `MIN_AVAIL_GB=999999` | 5 |

- [ ] `harness/concordance.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: concordance.sh <calls.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor-sample>
  # Calls -> FILTER=PASS & STATUS in {StrongSomatic,LikelySomatic}; truth FILTER stripped (else rtg drops it);
  # rtg vcfeval with sites-only baseline (sample ALT). Prints the 'None' row as JSON; fails on empty baseline.
  set -euo pipefail
  CALLS=$1; TRUTH=$2; BED=$3; SDF=$4; OUT=$5; TS=$6; mkdir -p "$OUT"
  bcftools view -f PASS -i 'INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"' "$CALLS" \
    | bcftools sort -Oz -o "$OUT/calls.som.vcf.gz"; tabix -f -p vcf "$OUT/calls.som.vcf.gz"
  bcftools annotate -x FILTER "$TRUTH" | bcftools sort -Oz -o "$OUT/truth.nofilter.vcf.gz"
  tabix -f -p vcf "$OUT/truth.nofilter.vcf.gz"
  rm -rf "$OUT/vcfeval"
  rtg vcfeval -b "$OUT/truth.nofilter.vcf.gz" -c "$OUT/calls.som.vcf.gz" -t "$SDF" -e "$BED" \
    --squash-ploidy --sample "ALT,$TS" --vcf-score-field QUAL -T 8 -o "$OUT/vcfeval"
  python3 - "$OUT/vcfeval/summary.txt" > "$OUT/summary.json" <<'PY'
  import json, sys
  rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith(("Threshold", "-"))]
  none = [r for r in rows if r[0] == "None"]
  if len(none) != 1: sys.exit("concordance: no unique 'None' row in summary.txt")
  r = none[0]
  d = dict(tp_baseline=int(r[1]), tp_call=int(r[2]), fp=int(r[3]), fn=int(r[4]),
           precision=r[5], sensitivity=r[6], f_measure=r[7])
  if d["tp_baseline"] + d["fn"] == 0: sys.exit("concordance: empty baseline in eval regions (vacuous)")
  print(json.dumps(d, sort_keys=True))
  PY
  for f in tp tp-baseline fp fn; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$OUT/vcfeval/$f.vcf.gz" > "$OUT/$f.sites"
  done
  cat "$OUT/summary.json"
  ```
- [ ] `harness/concordance_pair.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: concordance_pair.sh <name> <java.vcf> <rust.vcf> <truth.vcf.gz> <eval.bed> <ref.sdf> <outdir> <tumor>
  # Gate (D-27): identical 'None'-row counts AND identical tp/tp-baseline/fp/fn site lists. Not a dF1 tolerance.
  set -euo pipefail
  NAME=$1; J=$2; R=$3; TRUTH=$4; BED=$5; SDF=$6; OUT=$7; TS=$8
  bash harness/concordance.sh "$J" "$TRUTH" "$BED" "$SDF" "$OUT/java" "$TS" >/dev/null
  bash harness/concordance.sh "$R" "$TRUTH" "$BED" "$SDF" "$OUT/rust" "$TS" >/dev/null
  res=EQUAL; cmp -s "$OUT/java/summary.json" "$OUT/rust/summary.json" || res=DIFFERENT
  for f in tp tp-baseline fp fn; do cmp -s "$OUT/java/$f.sites" "$OUT/rust/$f.sites" || res=DIFFERENT; done
  python3 -c 'import json,sys; o=sys.argv[1]; json.dump({"java": json.load(open(o+"/java/summary.json")), "rust": json.load(open(o+"/rust/summary.json"))}, open(o+"/extra.json","w"), indent=2, sort_keys=True)' "$OUT"
  python3 harness/record_evidence.py "$NAME" "$res" "$OUT" -- "$0" "$@"
  [ "$res" = EQUAL ]
  ```
- [ ] **Red→green** `harness/tests/test_concordance.sh`, run under `conda run -n g4-tools`:
  - Build an SDF from `fixtures/micro/micro.ref.fa` with `rtg format`.
  - Take the first contig, length L. At positions L/4, L/2 and 3L/4, read REF via `samtools faidx`; ALT is `C` if REF is `A`, else `A`.
  - Write a sites-only truth VCF and a calls VCF. The calls VCF needs `##contig`, `##INFO=<ID=STATUS,Number=1,Type=String,…>`, `GT`, samples `TUMOR NORMAL`, PASS and `STATUS=StrongSomatic`.
  - The eval BED is the whole contig.
  - Assert `tp_baseline=3, fp=0, fn=0`. Dropping one call gives `fn=1`. An empty truth set must exit non-zero.
  - If rtg rejects `--sample ALT,<name>`, fix the invocation in this task and document why.
- [ ] Java heap: `conda run -n vardict-oracle bash -c 'grep -c JAVA_OPTS "$(which vardict-java)"'` must be ≥1. If it is 0, record whichever variable the launcher actually honours and use that in T15/T19.
- [ ] Commit `T12: g4-tools env, courtesy wrapper (A3), rtg concordance with exact-equality pair gate`.

### T13 — Data fetch + provenance (D-23, D-25)
- [ ] `harness/assets.tsv` (tab-separated; `-` = no upstream md5):
  ```
  # group	name	url	md5	dest
  refs	GRCh38.d1.vd1.fa	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa	-	ref/d1vd1/GRCh38.d1.vd1.fa
  refs	GRCh38.d1.vd1_BWA.tar.gz	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/technical/reference_genome/GRCh38/GRCh38.d1.vd1_BWA.tar.gz	015f5223bddd93b6e8f7a038c171f7be	ref/d1vd1/GRCh38.d1.vd1_BWA.tar.gz
  refs	GRCh38_GIABv3.fasta.gz	https://ftp.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz	-	ref/giabv3/GRCh38_GIABv3.fasta.gz
  seqc2	SRR7890879_1	https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/009/SRR7890879/SRR7890879_1.fastq.gz	52a511c2581212126c2e6f5a1743ec24	seqc2/fastq/SRR7890879_1.fastq.gz
  seqc2	SRR7890879_2	https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/009/SRR7890879/SRR7890879_2.fastq.gz	49469650695770b476851a61ed555113	seqc2/fastq/SRR7890879_2.fastq.gz
  seqc2	SRR7890880_1	https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/000/SRR7890880/SRR7890880_1.fastq.gz	355185fb242a4e672dffad2f34892094	seqc2/fastq/SRR7890880_1.fastq.gz
  seqc2	SRR7890880_2	https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/000/SRR7890880/SRR7890880_2.fastq.gz	6e03fe945d378a6ab1c42fcfa0e93999	seqc2/fastq/SRR7890880_2.fastq.gz
  seqc2	sSNV_v1.2.1	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/release/v1.2.1/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz	-	seqc2/truth/sSNV.vcf.gz
  seqc2	sINDEL_v1.2.1	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/release/v1.2.1/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz	-	seqc2/truth/sINDEL.vcf.gz
  seqc2	HC_regions_v1.2	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/release/v1.2.1/High-Confidence_Regions_v1.2.bed	-	seqc2/truth/HC.bed
  seqc2	exome_v6_UTR_hg38	https://ftp-trace.ncbi.nlm.nih.gov/seqc/ftp/Somatic_Mutation_WG/technical/reference_genome/Exome_Target_bed/S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz	-	seqc2/exome.bed6.gz
  ```
- [ ] **Verify every URL** with `curl -sfIL "$url" | head -1`, and probe `"$url.md5"`; if a published md5 exists, put it in the tsv. **Any 404 → STOP.** Do not substitute files.
- [ ] **Add the HG008 rows** (group `hg008`): the DRAFT v0.3 truth VCF (+`.tbi`) and the `_all`, `_nogermlineinterference` and `_nogermlinewithin50bp` BEDs.
  - Get the exact filenames from `curl -s https://ftp.ncbi.nlm.nih.gov/giab/ftp/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-smvar_DraftBenchmark_V0.3-20260425/ | grep -o 'href="[^"]*"'`.
  - The BAMs are **not** downloaded; T15 slices them remotely.
- [ ] `harness/fetch_data.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: fetch_data.sh <group>...  (refs seqc2 hg008). Resumable; md5-verified; appends benchmarks/MANIFEST.
  set -euo pipefail
  DATA=${DATA_DIR:-scratch/g4}; MAN=benchmarks/MANIFEST; mkdir -p "$DATA" benchmarks
  [ -f "$MAN" ] || printf '# name\turl\tmd5\tmd5_source\tbytes\tfetched\n' > "$MAN"
  want=" $* "
  while IFS=$'\t' read -r group name url md5 dest; do
    case "$group" in ''|\#*) continue ;; esac
    [[ "$want" == *" $group "* ]] || continue
    out="$DATA/$dest"; mkdir -p "$(dirname "$out")"
    if [ -f "$out.done" ]; then echo "have: $name"; continue; fi
    echo "fetch: $name"; curl -fL --retry 5 --retry-delay 30 -C - -o "$out" "$url"
    got=$(md5sum "$out" | cut -d' ' -f1); src=computed-locally
    if [ "$md5" != "-" ]; then
      [ "$got" = "$md5" ] || { echo "fetch: MD5 MISMATCH $name got=$got want=$md5" >&2; exit 1; }; src=upstream
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "$src" "$(stat -c %s "$out")" "$(date -Is)" >> "$MAN"
    touch "$out.done"
  done < harness/assets.tsv
  ```
- [ ] **Run** (≈25 GB; may run in the background): `EST_MIN=240 bash harness/courtesy.sh bash harness/fetch_data.sh refs hg008 seqc2`.
- [ ] Prepare the references:
  - `gunzip -k` the GIABv3 FASTA.
  - Untar the BWA index.
  - `samtools faidx` both FASTAs.
  - Build an SDF for each with `rtg format -o <ref>.sdf <ref>.fa`.
  - Append the md5s of the derived files to the MANIFEST.
- [ ] Commit `T13: assets.tsv (live-verified URLs/md5), fetch_data.sh, benchmarks/MANIFEST`.

### T14 — Slice-region selection (D-35)
- [ ] `harness/data/grch38_exclude_ig_hla.bed`:
  ```
  chr2	88857360	90235368	IGK
  chr6	28510119	33480577	MHC_HLA
  chr14	105586436	106879844	IGH
  chr22	22026075	22922913	IGL
  ```
- [ ] **Red first:** write `harness/tests/test_select_slice_regions.py` (unittest, temp files) and watch it fail (module absent). It must assert:
  - 6 chr1 variants within 1 kb plus 2 chr2 variants, with `--min-count 5` → exactly one chr1 region, padded.
  - An exclude BED covering the chr1 cluster → `SystemExit`.
  - An eval BED that drops 2 of the 6 chr1 variants → `SystemExit`.
  - Two equal-size clusters with `--max-regions 1` → the lexicographically first chromosome is chosen.
- [ ] `harness/select_slice_regions.py`:
  ```python
  #!/usr/bin/env python3
  """Pick dense truth-variant clusters as slice regions (deterministic).
  Output: 4-col BED (0-based, half-open) `chrom start end slice_<i>_<count>`, merged, sorted."""
  import argparse, bisect, gzip, sys

  def merge(iv):
      out = []
      for s, e in sorted(iv):
          if out and s <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], e))
          else: out.append((s, e))
      return out

  def read_bed(path):
      d = {}
      with open(path) as fh:
          for ln in fh:
              if not ln.strip() or ln.startswith(("#", "track", "browser")): continue
              c, s, e = ln.rstrip("\n").split("\t")[:3]; d.setdefault(c, []).append((int(s), int(e)))
      return {c: merge(v) for c, v in d.items()}

  def inside(iv, c, p0):
      lst = iv.get(c, []); i = bisect.bisect_right(lst, (p0, float("inf"))) - 1
      return i >= 0 and lst[i][0] <= p0 < lst[i][1]

  def vcf_sites(path):
      op = gzip.open if path.endswith(".gz") else open
      with op(path, "rt") as fh:
          for ln in fh:
              if ln.startswith("#"): continue
              f = ln.split("\t", 5); yield f[0], int(f[1]) - 1, len(f[3])

  def clusters(sites, chain, max_span):
      cur = []
      for c, p, rl in sites:
          if cur and (c != cur[0][0] or p - cur[-1][1] > chain or p + rl - cur[0][1] > max_span):
              yield cur; cur = []
          cur.append((c, p, rl))
      if cur: yield cur

  def main(argv=None):
      ap = argparse.ArgumentParser()
      ap.add_argument("--truth", action="append", required=True); ap.add_argument("--eval", required=True)
      ap.add_argument("--exclude"); ap.add_argument("--chain", type=int, required=True)
      ap.add_argument("--min-count", type=int, required=True); ap.add_argument("--pad", type=int, default=1000)
      ap.add_argument("--max-span", type=int, required=True); ap.add_argument("--max-regions", type=int, required=True)
      a = ap.parse_args(argv)
      ev = read_bed(a.eval); ex = read_bed(a.exclude) if a.exclude else {}
      sites = sorted({s for t in a.truth for s in vcf_sites(t) if inside(ev, s[0], s[1]) and not inside(ex, s[0], s[1])})
      cl = [k for k in clusters(sites, a.chain, a.max_span) if len(k) >= a.min_count]
      if not cl: sys.exit("select_slice_regions: no cluster meets --min-count; change parameters deliberately and record why")
      cl.sort(key=lambda k: (-len(k), k[0][0], k[0][1]))
      regs = {}
      for k in cl[:a.max_regions]:
          c = k[0][0]; regs.setdefault(c, []).append((max(0, k[0][1] - a.pad), max(p + rl for _, p, rl in k) + a.pad))
      n = 0
      for c in sorted(regs):
          for s, e in merge(regs[c]):
              n += 1; cnt = sum(1 for cc, p, _ in sites if cc == c and s <= p < e)
              print(f"{c}\t{s}\t{e}\tslice_{n}_{cnt}")

  if __name__ == "__main__": main()
  ```
- [ ] **Green:** the tests pass.
- [ ] Commit `T14: deterministic truth-cluster slice selection + GRCh38 IG/HLA exclude`.

### T15 — HG008 paired Tier-2 on remote slices (D-25, D-28, D-44)
- [ ] **Set up the base URL:** `U=https://ftp.ncbi.nlm.nih.gov/giab/ftp/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313`. Confirm the index suffix (`.bam.bai` or `.bai`) with `curl -sfI`.
- [ ] **@SQ gate:** the BAM header contigs must match the FASTA index exactly; otherwise STOP.
  - For T: `samtools view -H "$U/HG008-T_Illumina_195x_GRCh38-GIABv3.bam" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2);sub("LN:","",$3);print $2"\t"$3}' | sort > sq.T`.
  - Do the same for N, writing `sq.N`.
  - `cut -f1,2 GRCh38_GIABv3.fa.fai | sort` must be identical to both.
- [ ] **Regions:**
  - Run `python3 harness/select_slice_regions.py --truth <draft.vcf.gz> --eval <draft_all.bed> --exclude harness/data/grch38_exclude_ig_hla.bed --chain 20000 --min-count 3 --pad 500 --max-span 100000 --max-regions 200 > scratch/g4/hg008/slices.bed`.
  - Engine BED (D-28): `bedtools makewindows -b slices.bed -w 20000 -i src > engine.bed`.
  - Eval BED: `bedtools intersect -a <draft_all.bed> -b slices.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > eval.bed`.
  - Record the region count, total bp and truth count.
- [ ] **Remote slice.** Run in `scratch/g4/hg008/`, because htslib writes the `.bai` there. Wrap it as `EST_MIN=180 bash harness/courtesy.sh bash -c '…'`:
  ```bash
  for s in T:HG008-T_Illumina_195x_GRCh38-GIABv3 N:HG008-N-D_Illumina_169x_GRCh38-GIABv3; do
    t=${s%%:*}; f=${s#*:}
    samtools view -@4 -b -M -L slices.bed -X -o HG008-$t.unsorted.bam "$U/$f.bam" "$U/$f.bam.bai"
    samtools sort -@4 -o HG008-$t.slice.bam HG008-$t.unsorted.bam && samtools index HG008-$t.slice.bam
    rm HG008-$t.unsorted.bam
  done
  ```
  Record the slice md5s in the MANIFEST.
- [ ] **Tier-2:**
  ```bash
  CMD=(env KEEP_DIR=scratch/g4/hg008/cmp MIN_VARIANTS=1000 JAVA_OPTS=-Xmx8g bash harness/compare.sh
       scratch/g4/ref/giabv3/GRCh38_GIABv3.fa
       "scratch/g4/hg008/HG008-T.slice.bam|scratch/g4/hg008/HG008-N.slice.bam"
       scratch/g4/hg008/engine.bed "HG008-T|HG008-N" -th 8)
  rc=0; EST_MIN=300 bash harness/courtesy.sh "${CMD[@]}" || rc=$?
  python3 harness/record_evidence.py hg008-tier2 "$([ $rc -eq 0 ] && echo MATCH || echo MISMATCH)" scratch/g4/hg008/cmp -- "${CMD[@]}"
  ```
- [ ] Report a floor failure as the observed count; lowering the floor requires user approval (§6.6). Meanwhile, investigate why the count is low.
- [ ] A MISMATCH is dispatched **immediately** as an engine-debug task (A7): `superpowers:systematic-debugging`, one function, using a minimal fixture extracted from the mismatching region via `harness/slice.sh`. The task goes through R-GATE, then this step is re-run until it exits 0. A harness change is never the fix.
- [ ] Commit `T15: HG008 paired Tier-2 (remote slices, -th 8) -> <MATCH n>`.

### T16 — CRAM twin (D-36)
- [ ] M5 check:
  - If the BAM `@SQ` lines carry `M5:`, those values must equal the M5 values from `samtools dict GRCh38_GIABv3.fa`.
  - If they don't, compare names and lengths only, and record "no M5 in BAM header".
- [ ] Convert both slices: `samtools view -@4 -C -T <GIABv3.fa> -o HG008-$t.slice.cram HG008-$t.slice.bam && samtools index HG008-$t.slice.cram`.
- [ ] Ours on CRAM, using **exactly** compare.sh's Rust argument vector:
  ```bash
  D=scratch/g4/hg008; mkdir -p $D/cram
  EST_MIN=60 bash harness/courtesy.sh /usr/bin/time -v -o $D/cram/rust.time target/release/vardictrust \
    -G scratch/g4/ref/giabv3/GRCh38_GIABv3.fa -f 0.01 -N "HG008-T|HG008-N" \
    -b "$D/HG008-T.slice.cram|$D/HG008-N.slice.cram" -c 1 -S 2 -E 3 -g 4 -z -th 8 $D/engine.bed > $D/cram/rust.vcf
  ok=1
  cmp <(grep -v '^##' $D/cmp/rust.vcf) <(grep -v '^##' $D/cram/rust.vcf) || ok=0
  python3 harness/compare_vcf.py --min-variants 1000 $D/cmp/java.vcf $D/cram/rust.vcf | tee $D/cram/compare.txt || ok=0
  python3 harness/record_evidence.py hg008-cram-twin "$([ $ok = 1 ] && echo MATCH || echo MISMATCH)" $D/cram -- "see T16"
  ```
- [ ] Commit `T16: CRAM twin — ours(CRAM) == ours(BAM) and == oracle(BAM)`.

### T17 — HG008 draft-truth concordance (D-26, D-27, D-37)
- [ ] Run `EST_MIN=30 bash harness/courtesy.sh conda run -n g4-tools bash harness/concordance_pair.sh hg008-truth scratch/g4/hg008/cmp/java.vcf scratch/g4/hg008/cmp/rust.vcf <draft.vcf.gz> scratch/g4/hg008/eval.bed <GIABv3.sdf> scratch/g4/hg008/truth HG008-T` → exit 0 (EQUAL).
- [ ] Label the absolute numbers **"vs GIAB HG008 DRAFT v0.3 (not a released benchmark)"**.
- [ ] Truvari runs only if the truth contains SVs (`bcftools view -H -i 'INFO/SVTYPE!="."' <truth> | head -1`). Otherwise record "N/A — no SV truth" as input for the item-6 waiver.
- [ ] Commit `T17: HG008 draft-truth concordance — java == rust (exact counts + site lists)`.

### T18 — SEQC2 alignment (D-23)
- [ ] `harness/seqc2_align.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: seqc2_align.sh <SM> <R1.fq.gz> <R2.fq.gz> <bwa-index-prefix> <out.bam>  — <= 8 threads total (A3).
  # Deliberately simpler than SEQC2's trimmomatic+Picard pipeline: affects absolute recall, never parity.
  set -euo pipefail
  SM=$1; R1=$2; R2=$3; IDX=$4; OUT=$5; tmp="$OUT.tmp"; mkdir -p "$tmp"
  bwa mem -M -t 6 -R "@RG\tID:$SM\tSM:$SM\tPL:ILLUMINA\tLB:$SM" "$IDX" "$R1" "$R2" \
    | samtools fixmate -m - - | samtools sort -@1 -m 2G -T "$tmp/s" -o "$tmp/sorted.bam" -
  samtools markdup -@2 "$tmp/sorted.bam" "$OUT"; samtools index "$OUT"; rm -rf "$tmp"
  ```
- [ ] Align one sample at a time, starting after 05:15. Wrap each as `EST_MIN=420 bash harness/courtesy.sh conda run -n g4-tools bash harness/seqc2_align.sh …`:
  - `HCC1395` from `SRR7890879_{1,2}`;
  - then `HCC1395BL` from `SRR7890880_{1,2}`.
  - Output BAMs go to `scratch/g4/seqc2/`.
- [ ] @SQ gate against `GRCh38.d1.vd1.fa.fai`, as in T15.
- [ ] `samtools flagstat` must show ≥95% mapped; record the figure.
- [ ] Add the BAM md5s to the MANIFEST.
- [ ] Commit `T18: SEQC2 HCC1395/BL aligned from ENA FASTQ to GRCh38.d1.vd1 (script + MANIFEST)`.

### T19 — SEQC2 Tier-2: slices, then full exome (D-24, D-28, D-44)
- [ ] Engine BED: `zcat exome.bed6.gz | awk -F'\t' '$1 ~ /^chr([0-9]+|X|Y)$/' | cut -f1-4 | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o first > exome.engine.bed`. Record the interval count and bp (≈90 Mb).
- [ ] Eval BED: `bedtools intersect -a exome.engine.bed -b HC.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > eval.bed`. Count the truth variants inside it. The review computed ≈1,159 SNVs and 48 indels. Record the actual numbers; explain any difference rather than forcing a match.
- [ ] Slices:
  - `select_slice_regions.py --truth sSNV.vcf.gz --truth sINDEL.vcf.gz --eval eval.bed --exclude harness/data/grch38_exclude_ig_hla.bed --chain 250000 --min-count 5 --pad 1000 --max-span 2000000 --max-regions 27 > slices.bed`
  - `bedtools intersect -a exome.engine.bed -b slices.bed | sort -k1,1 -k2,2n | bedtools merge -c 4 -o first > slices.engine.bed`
- [ ] Slice Tier-2 → evidence `seqc2-slices-tier2`:
  - Settings: `KEEP_DIR=scratch/g4/seqc2/cmp-slices MIN_VARIANTS=200 JAVA_OPTS=-Xmx8g`, sample `"HCC1395|HCC1395BL"`, `-th 8`, `EST_MIN=120`.
  - Use T15's `CMD`/`record_evidence.py` pattern.
- [ ] **Only after the slices MATCH**, run full-exome Tier-2 → evidence `seqc2-exome-tier2`:
  - Settings: `KEEP_DIR=…/cmp-exome MIN_VARIANTS=2000`, `EST_MIN=480`.
  - Per-engine wall time and RSS land in the JSON automatically.
- [ ] A floor that isn't met → report it; changing a floor needs user approval.
- [ ] Commit `T19: SEQC2 HCC1395 Tier-2 — slices <MATCH n>, full exome <MATCH n> (-th 8)`.

### T20 — SEQC2 truth concordance
- [ ] `tabix` both truth files, then `bcftools concat -a -Oz -o truth.vcf.gz sSNV.vcf.gz sINDEL.vcf.gz`. FILTERs are stripped inside `concordance.sh`.
- [ ] `EST_MIN=30 bash harness/courtesy.sh conda run -n g4-tools bash harness/concordance_pair.sh seqc2-truth cmp-exome/java.vcf cmp-exome/rust.vcf truth.vcf.gz eval.bed <d1vd1.sdf> scratch/g4/seqc2/truth HCC1395` → EQUAL.
- [ ] Commit `T20: SEQC2 v1.2.1 concordance (exome∩HC) — java == rust`.

### T21 — Scale evidence (D-30, D-38). Thresholds are pre-registered here and must not move.
- [ ] `harness/scale_bench.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: scale_bench.sh <ref> <bam> <bed> <sample> <outdir>
  # Pre-registered (plan rev2): PASS iff wall(th1)/wall(th8) >= 4.8 AND maxRSS(chr1 full) <= 1.5 x maxRSS(10 Mb)
  # AND outputs identical across thread counts. Thresholds must not change after seeing data.
  set -euo pipefail
  REF=$1; BAM=$2; BED=$3; SM=$4; OUT=$5; BIN=${VARDICT_BIN:-target/release/vardictrust}; mkdir -p "$OUT"
  for th in 1 2 4 8; do
    /usr/bin/time -v -o "$OUT/th$th.time" "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z -th "$th" "$BED" > "$OUT/th$th.vcf"
  done
  for th in 2 4 8; do cmp <(grep -v '^##' "$OUT/th1.vcf") <(grep -v '^##' "$OUT/th$th.vcf"); done
  for sz in 1000000 10000000 100000000 248956422; do
    /usr/bin/time -v -o "$OUT/R$sz.time" "$BIN" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -th 8 -R "chr1:1-$sz" > "$OUT/R$sz.vcf"
  done
  python3 - "$OUT" > "$OUT/extra.json" <<'PY'
  import json, re, sys, pathlib
  o = pathlib.Path(sys.argv[1])
  def get(f, k): return re.search(rf"{re.escape(k)}: (.+)", (o / f).read_text()).group(1)
  def secs(w):
      p = [float(x) for x in w.split(":")]
      return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
  W = {t: secs(get(f"th{t}.time", "Elapsed (wall clock) time (h:mm:ss or m:ss)")) for t in (1, 2, 4, 8)}
  M = {s: int(get(f"R{s}.time", "Maximum resident set size (kbytes)")) for s in (1000000, 10000000, 100000000, 248956422)}
  sp = W[1] / W[8]; rr = M[248956422] / M[10000000]
  print(json.dumps({"wall_s": W, "maxrss_kb": M, "speedup_8": sp, "rss_ratio_chr1_vs_10mb": rr,
                    "pass": sp >= 4.8 and rr <= 1.5,
                    "note": ">8 threads not measured (host policy A3)"}, indent=2))
  PY
  ```
- [ ] Run it on the SEQC2 tumor BAM (single-sample `HCC1395`, bed = `exome.engine.bed`). Set `EST_MIN` to 4 × the Rust wall time from `seqc2-exome-tier2`, plus 60.
- [ ] Record `record_evidence.py scale <PASS|FAIL per extra.json "pass"> scratch/g4/scale -- …`.
- [ ] **Performance budget (A9, replaces CodSpeed; pre-registered, overridable only by you via OQ-13).** `harness/perf_budget.py`:
  ```python
  #!/usr/bin/env python3
  """Usage: perf_budget.py <evidence-name>...   Pre-registered (plan rev3, A9/OQ-13):
  PASS iff for every named evidence: rust wall <= java wall AND rust peak RSS <= java peak RSS.
  Inputs are the /usr/bin/time -v records compare.sh captured (both engines -th 8, same host, same courtesy).
  Java side = vardict-java | R | var2vcf pipeline (the thing VarDictRust consolidates); RSS = largest
  single child (getrusage RUSAGE_CHILDREN), i.e. the JVM vs vardictrust. Exit 0 = PASS."""
  import json, sys
  def secs(w):
      p = [float(x) for x in w.split(":")]
      return p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) == 3 else 0)
  rows, ok = {}, True
  for n in sys.argv[1:]:
      d = json.load(open(f"benchmarks/evidence/{n}.json"))
      j, r = d.get("java_time"), d.get("rust_time")
      if not j or not r or None in (j["wall"], r["wall"], j["max_rss_kb"], r["max_rss_kb"]):
          sys.exit(f"perf_budget: {n} lacks java/rust timing")
      jw, rw, jm, rm = secs(j["wall"]), secs(r["wall"]), int(j["max_rss_kb"]), int(r["max_rss_kb"])
      p = rw <= jw and rm <= jm; ok = ok and p
      rows[n] = {"java_wall_s": jw, "rust_wall_s": rw, "java_rss_kb": jm, "rust_rss_kb": rm,
                 "wall_ratio": round(rw / jw, 4), "rss_ratio": round(rm / jm, 4), "pass": p,
                 "evidence_fingerprint": d["fingerprint"]}
  print(json.dumps({"budget": "rust<=java wall and peak RSS (A9)", "runs": rows, "pass": ok}, indent=2, sort_keys=True))
  sys.exit(0 if ok else 1)
  ```
  - **Red→green** `harness/tests/test_perf_budget.sh`. Using two temp evidence JSONs:
    - rust faster and smaller → exit 0;
    - rust wall greater by 1 s → exit 1;
    - rust RSS greater by 1 kB → exit 1;
    - a missing `rust_time` → non-zero.
    - Clean up afterwards.
  - Run: `mkdir -p scratch/g4/perf && rc=0; python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2 > scratch/g4/perf/extra.json || rc=$?`, then `record_evidence.py perf-budget "$([ $rc -eq 0 ] && echo PASS || echo FAIL)" scratch/g4/perf -- python3 harness/perf_budget.py hg008-tier2 seqc2-exome-tier2`.
  - Add `ev perf-budget` to acceptance item 5 in the same commit (a tightening).
- [ ] **On any FAIL (scale or budget): keep working (A7).** Thresholds stay as they are.
  - Profile with `perf record -g` / `cargo flamegraph` on the failing input (user-level tools only) and fix the hot spot. Per CLAUDE §10, this may include `block-aligner` SIMD on the realignment path.
  - Each fix must keep R-GATE green. It changes the fingerprint, so re-run T15–T21 evidence.
  - Report the numbers at each checkpoint.
- [ ] Keep the existing criterion benches. `cargo bench --no-run` must compile in CI; the numbers are informational.
- [ ] Commit `T21: scale evidence + Rust<=Java wall/RSS budget -> <PASS>`.

### T22 — Evidence-driven docs (D-21, D-26, D-37)
- [ ] **AI-tool inventory from repo evidence:**
  - Historical trailers: `git log --format='%(trailers:key=Co-Authored-By)' | sort -u`. Old commits keep their trailers; history is not rewritten. New commits carry none (A10).
  - Meta dirs: `.claude` → Claude Code; `.qoder` → Qoder; `.specify` → spec-kit; the `copilot-instructions` reference → GitHub Copilot; this plan → Gemini/Antigravity CLI.
  - Put the list in the report so the user can confirm it at T25.
- [ ] README "AI Disclosure". Under A10 this README section is where AI involvement is disclosed:
  - every tool with its role (implementation, review, SDD orchestration);
  - what was validated, against which version and dataset: corpus 91/1, 10 slices, HG008 DRAFT, SEQC2 v1.2.1. Include the gate: ε = 1e-6 on raw values, header lines exact;
  - what was **not** validated: any untested mode such as Splicing, more than 8 threads, som.py, and non-x86_64 targets.
- [ ] README "Validation status": a table whose fields are copied from `benchmarks/evidence/*.json` (name, result, n_variants, fingerprint, date). The reviewer cross-checks every number against the JSON.
- [ ] `docs/validation/giab-seqc2-validation.md` covers:
  - datasets, with URLs and md5s from the MANIFEST;
  - references;
  - exact commands, from each evidence file's `cmd`;
  - results;
  - the DRAFT caveat;
  - how our SEQC2 preprocessing differs from theirs.
- [ ] DIVERGENCES.md: final statuses for s1a, PR #402, the `--fisher` residues, the SBF ulp, and >10-exception robustness. `docs/upstream-bugs.md`: add a ROADMAP id for each entry.
- [ ] Commit `T22: README/validation/DIVERGENCES from evidence JSON; full AI-tool disclosure`.

### T23 — Release engineering with cargo-dist (D-32, D-33, D-39, D-49, D-56, D-57)
- [ ] Root `[workspace.package]`: `version = "1.0.0"`, `readme = "README.md"`, `keywords = ["bioinformatics", "variant-calling", "genomics", "vcf", "ngs"]`, `categories = ["science", "command-line-utilities"]`. `repository` must be `https://github.com/palakpsheth/VarDictRust`; dist requires it.
- [ ] `vardict-core`:
  - `description = "Variant-calling engine of VarDictRust, an exact-output Rust port of VarDictJava 1.8.3"`
  - `readme`, `keywords`, `categories` and `repository` all use `.workspace = true`.
- [ ] Rename `crates/vardict-cli` package to `vardictrust` in `crates/vardict-cli/Cargo.toml` (`name = "vardictrust"` per A11):
  - `description = "vardictrust: single-binary VarDict variant caller, output-identical to VarDictJava 1.8.3"`
  - the same workspace keys;
  - `vardict-core = { path = "../vardict-core", version = "1.0.0" }`.
  - Enables direct `cargo install vardictrust` from crates.io.
- [ ] Run `cargo publish --workspace --dry-run`.
  - If `--workspace` is unsupported: run `cargo publish -p vardict-core --dry-run` plus `cargo package -p vardictrust --no-verify --list`, and record why.
  - Each package must be <10 MB, and no `VarDictJava/` or `scratch/` paths may be included.
- [ ] `cargo audit` and `cargo deny check` must report 0 errors.
  - If a dependency has to change, do it in a separate commit and rerun R-GATE.
  - This changes the fingerprint, so T15–T21 are re-run (A7; no pause).
- [ ] `security.yml`: add a `cargo install --locked cargo-audit && cargo audit` step.
- [ ] **cargo-dist (A8, A12).** Research basis: dist v0.33.0 (2026-09-11), docs `axodotdev/cargo-dist` `book/src/reference/config.md`.
  - Install it user-level: `cargo install cargo-dist --locked`. Then `dist --version` must print `0.33.x`.
  - `dist init --yes`, then make `dist-workspace.toml` exactly:
    ```toml
    [workspace]
    members = ["cargo:."]

    [dist]
    cargo-dist-version = "0.33.0"
    ci = "github"
    installers = ["shell"]
    targets = ["x86_64-unknown-linux-musl"]   # official parity-validated release artifact (A12)
    install-path = "CARGO_HOME"
    pr-run-mode = "upload"                    # PRs build + upload artifacts -> musl evidence (T25)
    cargo-auditable = true                    # embed dep tree for `cargo audit bin`
    github-attestations = true                # build provenance (experimental in dist)
    checksum = "sha256"
    publish-jobs = ["./publish-crates"]       # custom crates.io job, gated below

    [dist.dependencies.apt]
    musl-tools = "*"
    ```
  - Make the `[profile.dist]` that `dist init` adds to the root `Cargo.toml` exactly `inherits = "release"`, with no `lto` line, so the shipped code generation equals the validated release profile.
  - Run `dist generate`. It **replaces** the hand-rolled `.github/workflows/release.yml`. Diff the old file, and re-home anything worth keeping: the smoke gate moves to the `ci.yml` musl job and the T25 steps.
  - `.github/workflows/publish-crates.yml` (reusable; dist gives custom publish jobs `id-token: write`):
    ```yaml
    name: publish-crates
    on:
      workflow_call:
        inputs:
          plan: { required: true, type: string }
    jobs:
      publish:
        # Off until the first manual publish + trusted-publishing setup (T25 checkpoint 3).
        if: ${{ vars.CRATES_TRUSTED_PUBLISHING == 'true' && !fromJson(inputs.plan).announcement_is_prerelease }}
        runs-on: ubuntu-latest
        permissions: { id-token: write, contents: read }
        steps:
          - uses: actions/checkout@v4
            with: { submodules: recursive }
          - uses: dtolnay/rust-toolchain@stable
          - uses: rust-lang/crates-io-auth-action@v1
            id: auth
          - run: cargo publish -p vardict-core && cargo publish -p vardictrust
            env: { CARGO_REGISTRY_TOKEN: "${{ steps.auth.outputs.token }}" }
    ```
  - Verify locally:
    - `dist plan` lists exactly one artifact set, for `x86_64-unknown-linux-musl`.
    - `dist build --artifacts=local --target x86_64-unknown-linux-musl` succeeds (it uses the host `musl-gcc` for `bzip2-sys`/`lzma-sys`).
    - Extract `target/distrib/vardictrust-x86_64-unknown-linux-musl.tar.xz`.
    - `file <bin>` shows `statically linked`.
    - `VARDICT_BIN=<bin> bash harness/smoke.sh` → exit 0.
- [ ] **Binary equivalence.** This proves the shipped musl/dist binary produces the same output as the binary that produced the evidence. `harness/bin_equiv.sh`:
  ```bash
  #!/usr/bin/env bash
  # Usage: bin_equiv.sh <binA> <binB> — compare.sh's exact Rust argv on micro + every slices.tsv row;
  # outputs must be byte-identical. Exit 0 = identical everywhere.
  set -euo pipefail
  A=$1; B=$2; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fail=0
  rows() { printf 'micro\tfixtures/micro/micro.ref.fa\tfixtures/micro/micro.bam\tfixtures/micro/micro.bed\tSAMPLE\t1\t\n'
           grep -v '^#' harness/slices.tsv; }
  while IFS=$'\t' read -r id REF BAM BED SM _ EXTRA; do
    read -ra EX <<< "${EXTRA:-}"
    "$A" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "${EX[@]}" "$BED" > "$T/$id.A.vcf"
    "$B" -G "$REF" -f 0.01 -N "$SM" -b "$BAM" -c 1 -S 2 -E 3 -g 4 -z "${EX[@]}" "$BED" > "$T/$id.B.vcf"
    if cmp -s "$T/$id.A.vcf" "$T/$id.B.vcf"; then echo "SAME $id"; else echo "DIFF $id"; fail=1; fi
  done < <(rows)
  exit "$fail"
  ```
  - `bash harness/bin_equiv.sh target/release/vardictrust <dist-musl-bin>` → exit 0.
  - Record it with `record_evidence.py dist-bin-equiv PASS - -- …`.
  - Add `dist-bin-equiv` to acceptance item 8 in the same commit (a tightening).
- [ ] `ci.yml`: multi-arch matrix job with checkout (`submodules: recursive`) (A12):
  - Matrix targets:
    1. `x86_64-unknown-linux-gnu` (standard Linux build & test)
    2. `x86_64-unknown-linux-musl` (`sudo apt-get install -y musl-tools`, static build + `VARDICT_BIN=target/x86_64-unknown-linux-musl/release/vardictrust bash harness/smoke.sh`)
    3. `aarch64-unknown-linux-musl` (cross-compilation check)
    4. `x86_64-apple-darwin` or `aarch64-apple-darwin` (macOS runner build check)
  - Asserts that all crates compile cleanly and smoke tests pass on multiple target architectures.
- [ ] **rewrites.bio governance:**
  - `CITATION.cff`:
    - Use `cff-version: 1.2.0`, `title: VarDictRust`, `version: 1.0.0`, `license: MIT`, and `repository-code`.
    - Take `authors` from `[workspace.package].authors` / `git log --format=%an | sort -u`, confirmed by the user at T25 ⛔1.
    - `preferred-citation`: the article "VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research", Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, Dry JR. *Nucleic Acids Research* 44(11):e108, 2016, doi `10.1093/nar/gkw227`. Verify the author list against the DOI landing page.
    - `references`: VarDictJava 1.8.3 (AstraZeneca-NGS, `https://github.com/AstraZeneca-NGS/VarDictJava`, commit `0886261`).
    - Validate with `uvx --from cffconvert cffconvert --validate` → exit 0.
  - README "Upstream version policy":
    - VarDictRust v1.x tracks VarDictJava **v1.8.3** exactly. `--version` prints the oracle commit.
    - If upstream releases a new version, parity is re-established against it in a new minor release, with fresh evidence.
    - Corrections ship only as opt-in `--fix` flags per ROADMAP.md (A1).
    - Add the install lines generated by dist (shell installer, `cargo install vardictrust`).
- [ ] CHANGELOG:
  - Rename the section to `## [1.0.0] — <date>`.
  - Check every claim against a test or evidence file; rewrite any claim nothing supports (e.g. "Splicing/SV").
  - Add "Validated against VarDictJava 1.8.3" with the evidence names.
- [ ] R-GATE PASS, and `bash harness/engine_fingerprint.sh` must still equal the T15–T21 evidence fingerprint. If it doesn't, re-run the affected evidence (A7).
- [ ] Commit `T23: package vardictrust, cargo-dist musl release + multi-arch CI, bin-equivalence, CITATION.cff, upstream policy`.

### T24 — Final acceptance + whole-branch review
- [ ] Run `bash harness/acceptance.sh | tee scratch/acceptance/final.txt`.
  - Expected: items 1–7 PASS or WAIVED; item 8 FAILs **only** because `musl-ci` evidence is missing (T25 produces it).
  - Anything else becomes a new fix task (A7), and acceptance is re-run afterwards.
- [ ] Draft the proposed waivers in the report, but do **not** write them to `waivers.md`. Only items with no possible engineering fix qualify:
  - item 5: >8 threads not measured (host policy A3);
  - item 6: som.py and truvari N/A (no SV truth; rtg exact equality is the gate).
  - A parity, budget or threshold failure is **never** waivable (A7).
- [ ] Whole-branch review with `superpowers:requesting-code-review` over `main..release/v1.0.0`. It must cover:
  - §0b integrity: `git diff main -- harness/` contains tightenings only;
  - no panics and no std `HashMap` in core;
  - every D-xx closed with evidence.
- [ ] Write `docs/release/v1.0.0-pr.md`: the D-register closure table, the acceptance output, evidence links, and the provenance line.
- [ ] Commit `T24: final acceptance report + PR body`.

### T25 — Publish (controller only, ⛔)
- [ ] **⛔ Checkpoint 1.** Present the following to the user and get explicit approval for each:
  - the acceptance output;
  - the proposed waivers;
  - the AI-tool list;
  - the meta-file list (`.qoder .claude .specify specs *.docx KICKOFF_PROMPT.md docs/plans`; see OQ-8);
  - the secret scan: `conda run -n g4-tools gitleaks detect --source . --log-opts "main release/v1.0.0" --redact` → 0 leaks, plus a manual review of `git log -p main release/v1.0.0 | grep -niE 'sabnzbd|deluge|api_?key|ghp_|github_pat_'`.
  - Then write the approved waivers into `docs/release/waivers.md` and commit.
- [ ] `git push -u origin main && git push -u origin release/v1.0.0`. The first CI run on `main@4b295eb` is the pre-fix code and is expected to be red; say so in the PR.
- [ ] `gh pr create --base main --head release/v1.0.0 --title "v1.0.0: pure-parity release vs VarDictJava 1.8.3" --body-file docs/release/v1.0.0-pr.md`, then `gh pr checks --watch`. Everything must be green, including the `ci.yml` `musl` smoke job and dist's `plan`/`build-local-artifacts` jobs, which run because `pr-run-mode = "upload"`. A red check is engine/CI work (A7): fix, push, re-watch.
- [ ] Record the evidence:
  - `gh run download <dist-run-id> -D scratch/dist-pr` to fetch the dist PR artifact.
  - Extract it and run `VARDICT_BIN=<extracted>/vardictrust bash harness/smoke.sh` → exit 0.
  - Run `python3 harness/record_evidence.py musl-ci PASS - -- gh run view <ci-run-id> + gh run view <dist-run-id> + smoke`.
  - Commit and push, then `bash harness/acceptance.sh` must **exit 0**.
- [ ] **⛔ Checkpoint 2:** merge and tag.
  - `gh pr merge --merge`
  - `git switch main && git pull --ff-only`
  - `git tag -a v1.0.0 -m "vardictrust 1.0.0 — validated vs VarDictJava 1.8.3 (0886261/009e017); datasets: corpus 91/1, W4.6 slices 10/10, HG008 DRAFT v0.3 slices, SEQC2 HCC1395 WES v1.2.1; cmd harness/acceptance.sh"`
  - `git push origin v1.0.0`
  - `gh run watch` on dist's `release.yml`. Expected outputs:
    - the GitHub Release contains `vardict-cli-x86_64-unknown-linux-musl.tar.xz` (the name depends on OQ-3), its `.sha256`, `vardict-cli-installer.sh`, and `source.tar.gz`;
    - the attestation is present (`gh attestation verify <archive> -R palakpsheth/VarDictRust`);
    - `publish-crates` is skipped, because the variable is not set yet.
  - Download the released archive, check its sha256, and run `smoke.sh` on it.
- [ ] **⛔ Checkpoint 3:** crates.io.
  - The user runs `cargo login` in their own terminal.
  - Controller: `cargo publish -p vardict-core`, wait for the index, then `cargo publish -p vardict-cli`.
  - The user configures trusted publishing for both crates: GitHub, repo `palakpsheth/VarDictRust`, workflow **`release.yml`**. That is the caller workflow; the reusable `publish-crates.yml` runs inside it.
  - Then run `gh variable set CRATES_TRUSTED_PUBLISHING --body true`. From v1.0.1 on, dist publishes via OIDC.
- [ ] Update the ledger and add a CHANGELOG `[Unreleased]` section. The archive branch stays local (OQ-5).

---

## 5. PLAN §15 → evidence map

| §15 item | Gate in `acceptance.sh` | Produced by |
|---|---|---|
| 1 Unit | `cargo test --workspace` | T1, T3–T11 |
| 2 Tier-1 | Corpus `PASS: 91  FAIL: 1  SKIP: 0` with the exact fail-set; `smoke.sh` | existing + R-GATE |
| 3 Tier-2 | micro + 10 slices; `hg008-tier2`, `seqc2-slices-tier2`, `seqc2-exome-tier2` | T2, T10, T15, T19 |
| 4 Robustness | `robustness` tests; `hg008-cram-twin` | T9, T16 |
| 5 Scale | `scale`: pre-registered thread/RSS thresholds **and** the Rust ≤ Java wall/RSS budget (replaces CodSpeed, A9); waiver only for >8 threads | T21, T25 |
| 6 Truth | `hg008-truth` and `seqc2-truth` with exact equality; som.py/truvari waiver | T17, T20 |
| 7 Hygiene | fmt, clippy, `cargo audit`, `cargo deny check` | T6, T23 |
| 8 Distribution | `--version` content; `musl-ci` (dist PR "upload" run **and** the `ci.yml` musl smoke job) | T8, T23, T25 |

## 6. Pause conditions (A7: keep working until parity)

**Engine and harness problems are never stop conditions.** A MISMATCH, a missed variant floor, a scale/budget failure or fingerprint drift is turned by the controller into the next task. The task is dispatched right away (`superpowers:systematic-debugging`, one function per change), then the affected evidence is re-run. Progress is reported at checkpoints. Nothing is narrowed, waived or shipped divergent.

The controller **pauses and asks you** only for things an agent cannot or must not do:
1. Approvals marked ⛔: push, PR merge, tag, crates.io publish, and writing waivers.
2. Secrets: `cargo login`, and setting up trusted publishing in the crates.io UI.
3. Host changes. Anything needing sudo, or touching services, crontab, docker or systemd (A3).
4. Unavailable data. The controller first tries every documented mirror (NCBI FTP, ENA, the AWS/GCS GIAB buckets). It pauses only if a file still 404s, an md5 still fails, or an @SQ gate fails on the *published* file.
5. A real policy question. Example: a v1.8.3 behaviour that can't be emulated without a crash, where CLAUDE §5's "don't crash" rule and byte parity conflict. The controller proposes a resolution with evidence.
6. Changing a pre-registered threshold or variant floor. The controller never does this on its own; it reports the observed value and keeps working on the engine meanwhile.

**Checkpoint reports** (replacing rev-2 timeboxes) go to the ledger and to chat:
- at T10 +6h and T11-SBF +2h, then every 6h of engine work;
- after each external dataset;
- before T25.

Each report covers: what closed, what is still open (with its diff), what's next, and the ETA.

## 7. Open questions

All engineering and architecture choices are resolved:
- **OQ-1:** RESOLVED (ε = 1e-6 raw floats, A5)
- **OQ-2:** RESOLVED (SEQC2 gates release, A6)
- **OQ-3:** RESOLVED (Crate named `vardictrust` for direct `cargo install vardictrust`, A11)
- **OQ-6:** RESOLVED (Adopted cargo-dist v0.33.0, A8)
- **OQ-7:** RESOLVED (CodSpeed dropped, A9)
- **OQ-9:** RESOLVED (No `Co-Authored-By` commit trailers, A10)
- **OQ-10:** RESOLVED (Keep working until parity, A7)
- **OQ-11:** RESOLVED (Official release `x86_64-unknown-linux-musl` + multi-arch CI matrix in `ci.yml`, A12)
- **OQ-13:** RESOLVED (Baseline: Rust <= Java wall time & max RSS, A13)

The remaining open questions are strictly controller-side procedural checkpoints:
- **OQ-4** Crates.io authentication. The first crates.io publish requires you to run `cargo login` in your own shell (T25 ⛔3).
- **OQ-5** Git branch cleanup. Default: do not push `archive/pre-squash-20260823` to origin.
- **OQ-8** Repository hygiene. Default: keep `.qoder/.claude/.specify/specs`, the docx, and KICKOFF in the repo for AI transparency, after verifying clean gitleaks scan.
- **OQ-12** VCF header contig policy (#320). Default (T5a): the parity reference is the oracle invoked with var2vcf's own header options `-G <ref> -b <contigs BED built from .fai>` (`var2vcf_valid.pl:16,51-52`, `var2vcf_paired.pl:16,56-57`). Under that invocation Java prints `##contig=<ID=name,length=end>` for every `.fai` contig, which is what VarDictRust prints. So the #320 behaviour stays in v1.0. Alternative: Rust prints no `##reference`/`##contig` lines unless asked, which matches the bare `var2vcf` invocation, but VCF 4.2 recommends contig lines.
