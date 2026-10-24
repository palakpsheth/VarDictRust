# AGENTS.md — VarDictRust

> **Universal Coding Agent Contract** (Open standard for Claude, Cursor, Copilot, Codex, Gemini CLI, Antigravity, and AI pair-programmers; formerly `CLAUDE.md`).
> In-code citations citing `CLAUDE.md §X` refer directly to section §X of this file (`AGENTS.md`).

VarDictRust is an **idiomatic-Rust rewrite of VarDictJava v1.8.3** (a somatic/germline variant caller for
NGS), consolidating its fragmented Java→R→Perl pipeline into a single, statically-compiled, memory-safe
binary that **scales to deep WGS**. It is built under the [rewrites.bio](https://rewrites.bio) manifesto.

- **Reference implementation:** `./VarDictJava` (v1.8.3, commit `0886261`; Perl/R submodule `./VarDictJava/VarDict`, commit `009e017`).
- **Prime directive:** a faster, more robust tool that produces output **identical to VarDictJava v1.8.3 on valid inputs**. We are *not* writing a new caller and *not* transliterating Java line-by-line.
- **Architecture & Governance:** see [ROADMAP.md](ROADMAP.md) for candidate fixes, opt-in flags, and post-1.0 governance.
- **Testing & Verification Matrix:** see [docs/testing/test-matrix-and-ci-audit.md](docs/testing/test-matrix-and-ci-audit.md) and [harness/README.md](harness/README.md). This file is the always-loaded agent contract.

---

<a id="0-meta-rule--evidence-before-completion"></a>
## 0. Meta-rule — evidence before completion
**Never claim a task is done without evidence.** Evidence = a passing [`harness/compare.sh`](harness/README.md#key-scripts--runners) **MATCH**, a passing
`cargo test`, or a diff that shows a match. "It should work" / "this looks right" is **not** evidence. If you
cannot produce evidence, say so explicitly and keep working.

---

## 0b. Integrity — never fake a pass
The success criterion exists to catch you being wrong. **Do not defeat it.** You must NEVER:
- relax `--eps`, narrow a `--cols` range, compare fewer records, or edit/skip/`#[ignore]`/comment-out/delete a failing test or anything in `harness/` to turn a MISMATCH into a MATCH;
- hand-edit `golden/`, `expected/`, or oracle output to match your code (golden files change ONLY after a verified upstream correction, with explicit provenance documented in `golden/PROVENANCE`);
- manually edit engine fingerprints, variant counts, or metrics in `benchmarks/evidence/*.json` or fake runner receipts to satisfy acceptance gates. Modifying evidence JSON files by hand without executing the underlying live validation runners ([`harness/run_extended_coverage.sh`](harness/run_extended_coverage.sh), [`harness/run_multimodal_slices.sh`](harness/run_multimodal_slices.sh), or [`harness/recertify.sh`](harness/recertify.sh)) is a **Class-0 integrity violation**. Success criteria must come from verifiable runner execution per [docs/testing/test-matrix-and-ci-audit.md](docs/testing/test-matrix-and-ci-audit.md). If a fingerprint changes due to a dependency or code update, re-certify honestly via `bash harness/recertify.sh`;
- report a MATCH or a passing test you did not actually run, or paste expected output you did not observe.

If a check fails, **your code (or your understanding) is wrong — fix the port, not the test.** If you can't make it pass honestly, **stop and report the MISMATCH with the diff.** A real mismatch reported honestly is worth far more than a fake match.

**Oracle nondeterminism:** VarDictJava has known run-to-run output-*ordering* nondeterminism (issue #282; see `VarDictJava/tests/integration/test_non_deterministic_behavior`). The comparators align records by `chrom/pos/ref/alt`, so ordering alone never causes a false mismatch — if order is the only difference it is already a MATCH; never change the engine to chase it.

---

*The four behavioral principles below are Andrej Karpathy's ([multica-ai/andrej-karpathy-skills](https://github.com/multica-ai/andrej-karpathy-skills)), reproduced and merged with VarDictRust-specific notes (marked **VR:**) as that file invites.*

**Tradeoff:** these guidelines bias toward caution over speed. For trivial tasks, use judgment. They override habit, not the project rules in §5–§18.

## 1. Think Before Coding
**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them — don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.
- **VR:** before porting any function, read its VarDictJava source and restate in one sentence what it does; flag any ambiguity in the Java before writing Rust.

## 2. Simplicity First
**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.
- **VR:** the one exception is **fidelity** — replicating VarDict's exact heuristics and edge cases *is* the task, so it is never "speculative." But never add capability beyond what the port requires.

## 3. Surgical Changes
**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it — don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: every changed line should trace directly to the task.

## 4. Goal-Driven Execution
**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Port `isGoodVar`" → "write the unit test encoding VarDictJava's accept/reject table, then make it pass"
- "Fix the mismatch" → "make `bash harness/compare.sh <ref> <bam> <bed>` print MATCH"
- "Refactor the pileup" → "`bash harness/run_integration.sh` and `bash harness/rgate.sh` are green before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification. Here the criterion is almost always **`harness/compare.sh` exits 0 (MATCH)** or the named `cargo test` passes.

---

**These principles are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

---

## 5. Emulate exactly (the rewrites.bio prime directive)
Output must match **VarDictJava v1.8.3** on valid inputs:
- **Exact:** the variant set (chrom/pos/ref/alt), every FILTER decision, the INFO/FORMAT keys present, genotypes, TYPE, somatic STATUS, and the **intermediate column counts & order** — Simple **36/38**, Amplicon **38/40**, Somatic **55/61** (no-`--fisher` / with-`--fisher`).
- **Float tolerance:** Fisher p-values/SBF/SSF, odds ratios, AF/AdjAF/HIAF, QUAL means within **ε = 1e-6**, after replicating R's 5-decimal **HALF_EVEN** rounding (Fisher MLE tol `ulp(1)^0.25`).
- **Bug policy (v1.0):** if v1.8.3 has an **output** bug, **replicate it** and open a tracked issue — do **not** silently "fix" it. Where v1.8.3 **crashes** (e.g. `modifyCigar`), do **not** crash: produce the same output it would on valid inputs. That is allowed and is the whole point.
- **Three commitments:** credit the original ([README.md](README.md)/`--version`/[Lai et al. 2016](https://doi.org/10.1093/nar/gkw227)/AstraZeneca-NGS); emulate exactly; be transparent about AI (document tools, role, and exactly what was/wasn't validated, against which version + dataset).
- **Release policy (v1.0 Parity vs v1.1+ Corrections):** v1.0.x = parity with v1.8.3: no output-changing fixes and no fix flags — except *parity-restoring* fixes, which ship UNGATED as default. A fix qualifies iff (a) flag-on MATCHes the v1.8.3 oracle everywhere measured, (b) the prior default never matched the oracle on that path (Rust-only artifact, no upstream behavior replicated), and (c) its `--fix` id is unreleased (plumbing removed, not kept as no-op). Precedent: `splice-tile-ownership`, `mate-window-ownership` (ROADMAP "Landed as Default"). True output-changing corrections still ship from v1.1.0 per [ROADMAP.md](ROADMAP.md), each behind `--fix <id>` (default off) with a diff-confinement gate.

---

<a id="6-coordinate-indexing--the-1-silent-footgun"></a>
## 6. Coordinate indexing — the #1 silent footgun
- BAM/SAM positions from noodles are **0-based, half-open** → type `Pos0`.
- VCF output and VarDict-internal logic are **1-based** → type `Pos1`.
- Convert **only** via `.to_one_based()` / `.to_zero_based()`. Never add a raw integer to a position.
- If you are unsure whether a value is 0- or 1-based: **STOP and re-read the Java.** Do not guess. Coordinate errors are silent and clinically catastrophic.

---

## 7. The port loop (rewrites.bio "work small")
```
Pick ONE VarDictJava function
  → read its source (cite it: /// ports CigarModifier.java:61-71)
  → write the Rust port + its unit test
  → run harness/compare.sh `<case>`  (or `cargo test <name>`)
        MATCH   → commit, pick the next function
        MISMATCH→ debug THIS function; change nothing else
```
Never port multiple functions in one change without per-function evidence. "Output doesn't match on paired-end
data" is actionable; "port the whole pipeline" is not.

---

## 8. Fidelity hazards (these break output equivalence silently)
- **Iteration order is correctness.** Use `indexmap::IndexMap` where Java uses `LinkedHashMap`, `BTreeMap` where Java uses `TreeMap`. **Never** iterate a `HashMap` in an output/realignment path.
- **Rounding.** Every float written to a `.var`/VCF field goes through [`round_half_even`](DIVERGENCES.md#statistical-rounding--floating-point-half-even) (banker's rounding), never Rust's default. Replicate the formatter split: no-`--fisher` keeps trailing zeros (`x==0?0:fmt`); `--fisher` strips trailing zeros.
- **Sequences are bytes.** Use `&[u8]`/`Vec<u8>`, never `String`. Replace Java's negative-index `charAt`/`substr` with bounds-checked slice helpers returning `Option`.
- **Quality is Phred+33.** `Phred(byte.saturating_sub(33))`.
- **`maxReadLength` propagation** matters: SomaticMode seeds BAM2's pipeline with BAM1's value.
- **Java `HashMap` iteration that reaches output** is emulated with [`JavaHashMap`](DIVERGENCES.md#ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family). It ports `java.util.HashMap` spread, bin, resize and treeify-trigger semantics, and is validated against real Java by the op-log tests. This is the only sanctioned non-sorted order. Never use `std::collections::HashMap` there.

---

<a id="9-idiomatic-rust--the-no-panic-rule"></a>
## 9. Idiomatic Rust + the no-panic rule
- **No panics in library code.** In `crates/vardictrust-core/src/**` (outside `#[cfg(test)]`) you may **not** use `unwrap`, `expect`, `panic!`, `unreachable!`, bare slice indexing `x[i]` without a prior bounds check, or arithmetic that can overflow on input. Use `Result<_, VardictError>` + `?`, `.get(i)?`, checked/saturating math. **This is what fixes the VarDictJava crash-issue family by design.**
- **`#![forbid(unsafe_code)]`** in `vardictrust-core`. `unsafe` only via vetted SIMD crates (e.g. `block-aligner`) behind safe APIs.
- Model the data contract as **types**, not strings: `enum VariantDesc/CigarOp/VarType/SomaticStatus` each with a `Display` that reproduces the exact VarDict tokens (`+seq`, `-len`, `&`, `^`, `#`, `<DUP>`, `StrongSomatic`, …).
- `thiserror` for the library error enum, `anyhow` in the CLI. Errors must be actionable (no cryptic stacktraces).
- Pass `Arc<Config>` explicitly — no mutable global singleton.
- Prefer iterators/pattern-matching over index loops and string-switches — **while preserving iteration order**.
- Decompose large Java god-classes (e.g. the 2376-line `CigarParser`) into cohesive functions. Add traits only where there are ≥2 real implementations.

---

## 10. Scale rules (optimize for large WGS by default)
- **Bounded memory at any depth.** The unit of work is a region **tile** with a **halo**; never hold whole-chromosome state. Stream reads through the per-tile pipeline and drop them promptly (RAII).
- Tiled output **must equal** a single-region run (overlap + boundary-variant dedup). Scale optimizations must never change output.
- `rayon` work-stealing across tiles; **deterministic** ordered output (stable sort + ordered collect; bounded backpressure).
- Sort only within a tile (no global sort). Multi-thread BGZF. AVX2 (`block-aligner`) for the realignment hot path — but only after correctness is established, and re-validate via `bash harness/acceptance.sh` after.
- **Cluster Resource Etiquette & Concurrency Guidelines:**
  - Multi-threaded testing and variant calling should be matched to available cores without exhausting host memory (`-th <N>`, default all cores; for tests, `RUST_TEST_THREADS=4` or `--test-threads=8`).
  - On shared multi-tenant clusters or memory-constrained machines, wrap long-running benchmarks or differential tests with [`harness/courtesy.sh`](harness/README.md) (`nice -n 19 ionice -c 3`, with disk headroom `MIN_FREE_GB=150` and memory headroom `MIN_AVAIL_GB=12`).

---

## 11. Which VarDictJava issues we fix (v1.0 = robustness only)
v1.0 fixes **only** issues that change **no valid-input output** (so v1.0 stays byte-comparable to v1.8.3 *and*
is strictly more robust): the `modifyCigar` crash family (#341/#243/#253/#283/#7/#38/…), NPE/AIOOBE/div-by-zero,
OOM (#197/#384 — fixed by the tiling in §10), CRAM (#249/#380 — noodles + `-G` decode), locale-independent
formatting, fail-fast on missing index (#187), lenient SAM validation, nondeterminism (#282/#181/#356), and
correct `##contig` lengths from the `.fai` (#320). **Output-changing** corrections (#92 VD>DP, #62 over-filtering,
#363 realignment drops, #292 SHIFT3, #342 `<dup-N>` ALT, …) and **new features** (UMI, long reads) are deferred
to post-1.0 behind individual opt-in flags. See [ROADMAP.md](ROADMAP.md) and [docs/upstream-bugs.md](docs/upstream-bugs.md). Post-1.0 = v1.1.0+, governed by [ROADMAP.md](ROADMAP.md).

---

## 12. Unsupported features error loudly
If a flag/mode isn't implemented yet: `return Err(VardictError::Unsupported("--foo not implemented in v1.0; use VarDictJava"))`.
**Never** silently ignore a flag or emit wrong output.

---

## 13. Do NOT use these crates (hallucinated or unfit)
`genomicframe-core`, `rustalign-simd`, `ksw2rs`, `seqair`, `cyanea-align`.
Use: **noodles** (BAM/CRAM/FASTA/VCF/BGZF), **bio** + (later) **block-aligner** (alignment), **statrs** +
**fishers_exact** (stats), **rayon**, **clap**, **anyhow**/**thiserror**, **indexmap**.

---

## 14. Commands
**Setup (run once):** see [CONTRIBUTING.md](CONTRIBUTING.md) (Prerequisites) — rustup, miniconda, samtools, JDK, and reference-genome setup. Confirm every version prints before starting development.  
**Testing Hierarchy Details:** see [docs/testing/test-matrix-and-ci-audit.md](docs/testing/test-matrix-and-ci-audit.md) and [harness/README.md](harness/README.md).
```bash
# Build / test / lint / format — ALL must pass before any commit (zero warnings):
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo build --release

# Reference oracle (VarDictJava 1.8.3 via bioconda, env "vardict-oracle"):
bash harness/setup_oracle.sh
bash harness/run_oracle.sh <ref.fa> <bam|tumor.bam|normal.bam> <bed> <out.vcf> [sample]

# Differential testing (Acceptance gates — see harness/README.md):
# Tier 1 (Intermediate Table Parity): validates .var schema (36/38, 38/40, 55/61 cols), ratios, p-values
python3 harness/compare_var.py golden/<case>.var <ours.var>

# Tier 2 (Full VCF Parity): validates VCF coordinates, alleles, INFO/FORMAT tags, half-even float rounding
bash harness/compare.sh <ref> <bam|tumor.bam|normal.bam> <bed> [sample]

# Fast developer inner-loop gate (<30s):
bash harness/rgate.sh

# Slice big BAMs for the dev loop:
bash harness/slice.sh <in.bam> <chr:start-end> <out.bam>
```
**Pinned baseline reuse (never re-prove the same oracle twice):** VarDictJava is EOL/archived — its pin is immutable, so recorded java outputs are permanent goldens. vardictcpp stays on its pinned version; bump it only by deliberate decision (then regenerate cpp goldens). Reuse a baseline IFF inputs + flags + pin are identical, verified by sha manifest; otherwise rerun the oracle and store the new outputs + manifest entry (new tests run only their own cells, then append). Deterministic outputs gate; host-dependent measurements (wall/RSS) are recorded with host context, never gated across hosts.

## 15. Definition of done (every change)
1. `cargo fmt --all`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test` all pass.
2. The relevant `harness/compare.sh` prints **MATCH** (or the named test passes). Pre-PR compound gate passes: `bash harness/pr_check.sh`.
3. `golden/` updated if output legitimately changed (with provenance).
4. For any dependency (`Cargo.lock`) or engine code changes merged to `main`, run `bash harness/recertify.sh` (re-running extended coverage and multimodal slice suites) to re-certify live parity and update evidence in `benchmarks/evidence/` under [docs/testing/test-matrix-and-ci-audit.md](docs/testing/test-matrix-and-ci-audit.md).
5. The commit/PR states: **"validated vs VarDictJava 1.8.3, dataset `<name>`, cmd `<…>`"**.
6. `gitleaks protect --staged` (commit) / `gitleaks detect --source .` (push) clean.
7. The PR body completes `.github/pull_request_template.md`, which is authoritative on what may be skipped: engine-gate boxes are either run or marked **N/A with a one-line justification** under the template's docs/infra-only note (the same code-change filter the pre-push hook enforces).

## 16. Repo layout
```
crates/vardictrust-core/   engine library (all logic; #![forbid(unsafe_code)])
crates/vardictrust/    bin "vardictrust" (clap front-end → core)
harness/               compare.sh, compare_var.py, compare_vcf.py, slice.sh, run_oracle.sh, setup_oracle.sh
golden/                versioned oracle outputs (.vcf, .var) + PROVENANCE
benchmarks/            hardware, datasets, exact commands
VarDictJava/           the v1.8.3 reference (read-only; the source of truth for ports)
ROADMAP.md             post-1.0 candidate flags and architecture roadmap
```

## 17. Quick map: where things live in VarDictJava (for porting)
`src/main/java/com/astrazeneca/vardict/` — `modes/{Simple,Somatic,Amplicon,Splicing}Mode` (orchestration);
`modules/CigarModifier` (read-level CIGAR rewrite → `src/cigar/modify.rs`); `modules/CigarParser` (the pileup
engine → `src/cigar/parse.rs`); `modules/VariationRealigner` (→ `src/realign.rs`); `modules/ToVarsBuilder`
(variant finalization → `src/tovars.rs`); `modules/StructuralVariantsProcessor` (→ `src/sv.rs`);
`data/fishertest/{FisherExact,UnirootZeroIn}` (→ `src/stats/fisher.rs`); `postprocessmodules/*` (filter +
somatic STATUS → `src/postprocess/`); `printers/*` (→ `src/output/`); `variations/*` (→ `src/variations/`).
Somatic STATUS is computed in **Java** (`SomaticPostProcessModule.determinateType`), not Perl.

## 18. Commits
Atomic commits per task: `git add -A && git commit -m "<type>: <short desc>"`. Branch off `origin/main` (`git fetch origin && git checkout -b <type>/<slug> origin/main`); never push `main` directly (the pre-push hook blocks it) — merge via a PR opened with `gh`, with the PR template filled. Run `bash harness/setup_hooks.sh` once per clone so the repo-managed hooks are active. Full gate semantics live in [CONTRIBUTING.md](CONTRIBUTING.md) §§5–6 and `harness/pr_check.sh` (6 stages).
Never commit secrets, **data** (BAM/CRAM/FASTA/VCF beyond the tiny `golden/` fixtures), `target/`, conda envs, or ephemeral agent scratch/planning trees (e.g. `docs/superpowers/` — agent plans belong in the external `<appDataDir>/brain/<conversation-id>` artifact store; see `.gitignore`). **Do** commit `Cargo.lock` (this is an application workspace). MIT-licensed, open-source project. Do not include `Co-Authored-By` commit trailers; machine assistance and AI tool participation are disclosed transparently in [README.md](README.md) and [CITATION.cff](CITATION.cff).

<a id="19-repository-specific-agent-skills"></a>
## 19. Repository-Specific Agent Skills (`.agents/skills/`)
Specialized workflow skills live under `.agents/skills/<name>/SKILL.md` (validated by `harness/tests/test_bioconda_recipe.py`):
- **[`vardictrust-issue`](.agents/skills/vardictrust-issue/SKILL.md):** Filing, triaging, reproducing, or auditing GitHub issues (`v1.0.x` parity/robustness vs `v1.1+ --fix <id>` classification, `harness/slice.sh` minimal slices, and `compare.sh` oracle receipts).
- **[`vardictrust-pr`](.agents/skills/vardictrust-pr/SKILL.md):** Opening, updating, reviewing, recertifying (`harness/recertify.sh`), or squash-merging pull requests, completing `.github/pull_request_template.md`, and replying inline to review comments.
- **[`vardictrust-profiling`](.agents/skills/vardictrust-profiling/SKILL.md):** Profiling CPU (`profile-cpu`), jemalloc heap (`profile-heap`), or DHAT allocation churn (`profile-dhat`) via `harness/profile.sh` and `harness/profile_catalog.py` ([`references/rust-profiling-api.md`](.agents/skills/vardictrust-profiling/references/rust-profiling-api.md)).
- **[`vardictrust-release`](.agents/skills/vardictrust-release/SKILL.md):** Multi-channel release engineering across GitHub Releases (`cargo-dist` SLSA Level 3), Crates.io, Bioconda, Nextflow `nf-core/modules`, and ReadTheDocs ([`scripts/check_release_consistency.py`](.agents/skills/vardictrust-release/scripts/check_release_consistency.py) and [`references/version-surface-map.md`](.agents/skills/vardictrust-release/references/version-surface-map.md)).

**Agent Discovery:** `.agents/skills/` is committed in-tree and auto-discovered natively by Gemini CLI, OpenAI Codex, GitHub Copilot, Cursor, and OpenCode.
