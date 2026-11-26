---
name: vardictrust-issue
description: Use when filing, triaging, reproducing, or auditing GitHub issues in VarDictRust, classifying VarDictJava v1.8.3 parity mismatches versus crash bugs, or proposing v1.1+ --compat=off CompatProfile fixes.
---

# VarDictRust Issue Authoring & Triage (`vardictrust-issue`)

## Overview

VarDictRust is an idiomatic-Rust rewrite of VarDictJava v1.8.3 (`0886261`) built under the [rewrites.bio](https://rewrites.bio) contract (`AGENTS.md` §§5, 11). Every GitHub issue must be classified against the **`--compat=vardictjava` Parity vs `--compat=off` Corrections** release policy before code is touched.

**Core principle:** If VarDictJava v1.8.3 has a valid-input **output quirk**, default mode (`--compat=vardictjava` / `--parity`) replicates it and tracks any bonafide correctness fix as a `CompatProfile` switch under opt-in `--compat=off` (`v1.1+`, stamping `##VarDictRustCompat=off` in VCF headers — never per-bug `--fix` CLI flags) in `ROADMAP.md`. If VarDictJava **crashes, panics, OOMs, or exhibits non-deterministic map iteration**, VarDictRust fixes it immediately in default mode.

---

## Quick Reference: Issue Classification Matrix

| Symptom / Report Type | Classification | Release Target | Template & Governance | Required Evidence |
| :--- | :--- | :--- | :--- | :--- |
| `harness/compare.sh` or `compare_var.py` shows `MISMATCH` vs VarDictJava v1.8.3 | **Parity Bug** (Rust divergence) | Default `--compat=vardictjava` (fix immediately) | `.github/ISSUE_TEMPLATE/bug_report.yml` | Minimal `harness/slice.sh` BAM + `compare.sh` diff |
| VarDictRust panics, crashes, OOMs, or hangs on valid or malformed BAM/CRAM/FASTA/BED | **Robustness / Scale Bug** | Default mode (fix immediately) | `.github/ISSUE_TEMPLATE/bug_report.yml` (`AGENTS.md` §§9–11) | Stack trace / RSS log + minimal reproduction command |
| Both VarDictRust (`--compat=vardictjava`) and VarDictJava v1.8.3 emit biologically questionable output (`VD > DP`, over-filtering, `SHIFT3`) | **Upstream Output Quirk** (`docs/upstream-bugs.md`) | `v1.1+` `CompatProfile` switch under `--compat=off` (`ROADMAP.md`) | `.github/ISSUE_TEMPLATE/feature_request.yml` (`v1.1+ --compat=off`) | Upstream Java source citation (`Class.java:L1-L2`) + oracle output |
| Candidate fix matches v1.8.3 oracle everywhere measured & prior default was a Rust-only artifact | **Parity-Restoring Fix** (`AGENTS.md` §5) | Ship **ungated** in default `--compat=vardictjava` mode (not gated under `--compat=off`) | `.github/ISSUE_TEMPLATE/bug_report.yml` | 6-stage `harness/pr_check.sh` `MATCH` receipt per PR + 8-suite `recertify.sh`/`acceptance.sh` at release cut |
| Wall-clock or peak RSS regression vs certified `benchmarks/evidence/scale.json` | **Performance Regression** | Default mode (fix immediately) | `.github/ISSUE_TEMPLATE/bug_report.yml` | `harness/profile.sh` flamegraph / heap summary + `scale.json` delta |

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
1. Search `docs/upstream-bugs.md`, `DIVERGENCES.md`, and `ROADMAP.md` for known VarDictJava issues and shipped `v1.1+` `--compat=off` / `--vcf-all-alleles` items (`#341`, `#243`, `#253`, `#283`, `#282`, `#116`, `#342`, `#292`, `#369`, `PR #402`, `--eqx`, `UB-5`, `UB-6`, `#370`, `#92`, `#62`, `#363`).
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
- **Version:** `vardictrust 1.1.0` (commit `<sha>`, engine fingerprint `<fp>`)
- **OS / Arch:** `<uname -srm>`
- **Execution Mode:** `<Simple | Somatic | Amplicon | Splicing>` (`--fisher` / default; `--compat=vardictjava` default vs `--compat=off` / `--vcf-all-alleles`)

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
| "VarDictJava's output is clearly wrong here (`VD > DP`), so I'll file a bug to change the default `--compat=vardictjava` output." | Violates `AGENTS.md` §§5, 11. Default `--compat=vardictjava` (`--parity`) must match VarDictJava v1.8.3 on valid inputs. Propose a `CompatProfile` switch under `--compat=off` (`v1.1+`) via `feature_request.yml` and record in `ROADMAP.md`. |
| "It's a small mismatch; I don't need to run `harness/compare.sh` or cite the Java source lines." | Every parity issue requires observed `compare.sh` / `compare_var.py` output and the exact `VarDictJava` class/line citation. |
| "I'll upload the 2 GB clinical BAM so the issue is easy to reproduce." | Never upload PHI/patient data or large BAMs. Use `bash harness/slice.sh` or synthetic tests in `crates/vardictrust-core`. |

## Red Flags — STOP Immediately If You Catch Yourself

- Proposing an ungated default output change on valid inputs that diverges from VarDictJava v1.8.3.
- Filing a parity issue without checking whether the difference is `HashMap` ordering nondeterminism (`#282` — `compare.sh` aligns by `chrom/pos/ref/alt`).
- Omitting `vardictrust --version`, command flags, or the execution mode (`Simple` 36/38, `Amplicon` 38/40, `Somatic` 55/61 columns).
