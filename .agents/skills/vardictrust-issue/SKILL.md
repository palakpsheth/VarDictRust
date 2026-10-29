---
name: vardictrust-issue
description: Use when filing, triaging, reproducing, or auditing GitHub issues in VarDictRust, classifying VarDictJava v1.8.3 parity mismatches versus crash bugs, or proposing post-1.0 --fix flags.
---

# VarDictRust Issue Authoring & Triage (`vardictrust-issue`)

## Overview

VarDictRust is an idiomatic-Rust rewrite of VarDictJava v1.8.3 (`0886261`) built under the [rewrites.bio](https://rewrites.bio) contract (`AGENTS.md` §§5, 11). Every GitHub issue must be classified against the **v1.0 Parity vs v1.1+ Corrections** release policy before code is touched.

**Core principle:** If VarDictJava v1.8.3 has a valid-input **output quirk**, `v1.0.x` replicates it and tracks the correction for `v1.1.0+` behind an opt-in `--fix <id>` flag in `ROADMAP.md`. If VarDictJava **crashes, panics, OOMs, or exhibits non-deterministic map iteration**, `v1.0.x` fixes it immediately.

---

## Quick Reference: Issue Classification Matrix

| Symptom / Report Type | Classification | Release Target | Template & Governance | Required Evidence |
| :--- | :--- | :--- | :--- | :--- |
| `harness/compare.sh` or `compare_var.py` shows `MISMATCH` vs VarDictJava v1.8.3 | **Parity Bug** (Rust divergence) | `v1.0.x` (fix immediately) | `.github/ISSUE_TEMPLATE/bug_report.yml` | Minimal `harness/slice.sh` BAM + `compare.sh` diff |
| VarDictRust panics, crashes, OOMs, or hangs on valid or malformed BAM/CRAM/FASTA/BED | **Robustness / Scale Bug** | `v1.0.x` (fix immediately) | `.github/ISSUE_TEMPLATE/bug_report.yml` (`AGENTS.md` §§9–11) | Stack trace / RSS log + minimal reproduction command |
| Both VarDictRust and VarDictJava v1.8.3 emit biologically questionable output (`VD > DP`, over-filtering, `SHIFT3`) | **Upstream Output Quirk** (`docs/upstream-bugs.md`) | `v1.1.0+` behind `--fix <id>` (`ROADMAP.md`) | `.github/ISSUE_TEMPLATE/feature_request.yml` (`v1.1+ --fix <id>`) | Upstream Java source citation (`Class.java:L1-L2`) + oracle output |
| Candidate `--fix <id>` flag-on matches v1.8.3 oracle everywhere measured & prior default was a Rust-only artifact | **Parity-Restoring Fix** (`AGENTS.md` §5) | `v1.0.x` (ship **ungated**, remove `--fix` id) | `.github/ISSUE_TEMPLATE/bug_report.yml` | 6-stage `harness/pr_check.sh` `MATCH` receipt per PR + 8-suite `recertify.sh`/`acceptance.sh` at release cut |
| Wall-clock or peak RSS regression vs certified `benchmarks/evidence/scale.json` | **Performance Regression** | `v1.0.x` | `.github/ISSUE_TEMPLATE/bug_report.yml` | `harness/profile.sh` flamegraph / heap summary + `scale.json` delta |

---

## Step-by-Step Triage & Filing Workflow

### 1. Verify Environment & Version
Always capture the exact binary version, Rust toolchain, and engine fingerprint:
```bash
./target/release/vardictrust --version
rustc --version
bash harness/engine_fingerprint.sh
bash harness/engine_fingerprint.sh --perf
```

### 2. Isolate a Minimal Reproducible Slice
Never attach full patient BAMs or commit genomic fixtures outside `golden/`. Slice the Failing region to a minimal window with `harness/slice.sh`:
```bash
# Extract a minimal region (e.g. 500 bp around the discordant locus)
bash harness/slice.sh <input.bam> <chr:start-end> /tmp/repro_slice.bam

# Run Tier-1 (.var intermediate columns) and Tier-2 (full VCF) differential oracle checks
bash harness/compare.sh <ref.fa> /tmp/repro_slice.bam <region.bed> [sample]
python3 harness/compare_var.py <oracle.var> <ours.var>
```

### 3. Check Existing Upstream Catalog Before Filing
1. Search `docs/upstream-bugs.md`, `DIVERGENCES.md`, and `ROADMAP.md` for known VarDictJava issues (`#341`, `#243`, `#253`, `#283`, `#282`, `#92`, `#62`, `#363`, `#292`, `#342`).
2. Search open and closed repository issues:
   ```bash
   gh issue list --repo palakpsheth/VarDictRust --state all --search "<keyword>"
   ```

### 4. Author the Issue via `gh issue create`
Include all required fields from `.github/ISSUE_TEMPLATE/bug_report.yml` or `.github/ISSUE_TEMPLATE/feature_request.yml`:
```bash
gh issue create --repo palakpsheth/VarDictRust \
  --title "<type>(<scope>): <concise summary>" \
  --label "<bug|enhancement>" \
  --body "$(cat <<'EOF'
### VarDictRust Version & Environment
- **Version:** `vardictrust 1.0.6` (commit `<sha>`, engine fingerprint `<fp>`)
- **OS / Arch:** `<uname -srm>`
- **Execution Mode:** `<Simple | Somatic | Amplicon | Splicing>` (`--fisher` / default)

### Exact Command & Minimal Reproduction
```bash
bash harness/compare.sh <ref.fa> <slice.bam> <region.bed>
```

### Observed vs Expected (VarDictJava v1.8.3 Oracle)
- **VarDictJava v1.8.3 (`0886261`):** `<oracle record or behavior>`
- **VarDictRust:** `<observed record or behavior>`
- **Upstream Source Citation:** `VarDictJava/src/main/java/com/astrazeneca/vardict/<Class>.java:<lines>`
EOF
)"
```

---

## Rationalization Counter-Table

| Excuse / Rationalization | Reality / Mandatory Action |
| :--- | :--- |
| "VarDictJava's output is clearly wrong here (`VD > DP`), so I'll file a v1.0 bug to change the default output." | Violates `AGENTS.md` §§5, 11. `v1.0.x` must match VarDictJava v1.8.3 on valid inputs. File under `v1.1+ --fix <id>` and record in `ROADMAP.md`. |
| "It's a small mismatch; I don't need to run `harness/compare.sh` or cite the Java source lines." | Every parity issue requires observed `compare.sh` / `compare_var.py` output and the exact `VarDictJava` class/line citation. |
| "I'll upload the 2 GB clinical BAM so the issue is easy to reproduce." | Never upload PHI/patient data or large BAMs. Use `bash harness/slice.sh` or synthetic tests in `crates/vardictrust-core`. |

## Red Flags — STOP Immediately If You Catch Yourself

- Proposing an ungated default output change on valid inputs that diverges from VarDictJava v1.8.3.
- Filing a parity issue without checking whether the difference is `HashMap` ordering nondeterminism (`#282` — `compare.sh` aligns by `chrom/pos/ref/alt`).
- Omitting `vardictrust --version`, command flags, or the execution mode (`Simple` 36/38, `Amplicon` 38/40, `Somatic` 55/61 columns).
