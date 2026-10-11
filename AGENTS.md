# AGENTS.md — VarDictRust

> **Universal Coding Agent Contract** (Open standard for Claude, Cursor, Copilot, Codex, Gemini CLI, Antigravity, and AI pair-programmers; formerly `CLAUDE.md`).
> In-code citations citing `CLAUDE.md §X` refer directly to section §X of this file (`AGENTS.md`).

VarDictRust is an **idiomatic-Rust rewrite of VarDictJava v1.8.3** (a somatic/germline variant caller for
NGS), consolidating its fragmented Java→R→Perl pipeline into a single, statically-compiled, memory-safe
binary that **scales to deep WGS**. It is built under the **rewrites.bio** manifesto.

- **Reference implementation:** `./VarDictJava` (v1.8.3, commit `0886261`; Perl/R submodule `./VarDictJava/VarDict`, commit `009e017`).
- **Prime directive:** a faster, more robust tool that produces output **identical to VarDictJava v1.8.3 on valid inputs**. We are *not* writing a new caller and *not* transliterating Java line-by-line.
- **Full plan:** see `./ROADMAP.md` for candidate fixes, opt-in flags, and post-1.0 governance. This file is the always-loaded agent contract.

---

## 0. Meta-rule — evidence before completion
**Never claim a task is done without evidence.** Evidence = a passing `harness/compare.sh` **MATCH**, a passing
`cargo test`, or a diff that shows a match. "It should work" / "this looks right" is **not** evidence. If you
cannot produce evidence, say so explicitly and keep working.

---

## 0b. Integrity — never fake a pass
The success criterion exists to catch you being wrong. **Do not defeat it.** You must NEVER:
- relax `--eps`, narrow a `--cols` range, compare fewer records, or edit/skip/`#[ignore]`/comment-out/delete a failing test or anything in `harness/` to turn a MISMATCH into a MATCH;
- hand-edit `golden/`, `expected/`, or oracle output to match your code (golden files change ONLY via `/snapshot-update`, after a real, justified, documented output change);
- report a MATCH or a passing test you did not actually run, or paste expected output you did not observe.

If a check fails, **your code (or your understanding) is wrong — fix the port, not the test.** If you can't make it pass honestly, **stop and report the MISMATCH with the diff.** A real mismatch reported honestly is worth far more than a fake match.

**Oracle nondeterminism:** VarDictJava has known run-to-run output-*ordering* nondeterminism (issue #282; see `VarDictJava/tests/integration/test_non_deterministic_behavior`). The comparators align records by `chrom/pos/ref/alt`, so ordering alone never causes a false mismatch — if order is the only difference it is already a MATCH; never change the engine to chase it.

---

*The four behavioral principles below are Andrej Karpathy's (multica-ai/andrej-karpathy-skills), reproduced and merged with VarDictRust-specific notes (marked **VR:**) as that file invites.*

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
- "Fix the mismatch" → "make `harness/compare.sh <case>` print MATCH"
- "Refactor the pileup" → "`compare.sh --all` is green before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification. Here the criterion is almost always **`harness/compare.sh <case>` exits 0 (MATCH)** or the named `cargo test` passes.

---

**These principles are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

---

## 5. Emulate exactly (the rewrites.bio prime directive)
Output must match **VarDictJava v1.8.3** on valid inputs:
- **Exact:** the variant set (chrom/pos/ref/alt), every FILTER decision, the INFO/FORMAT keys present, genotypes, TYPE, somatic STATUS, and the **intermediate column counts & order** — Simple **36/38**, Amplicon **38/40**, Somatic **55/61** (no-`--fisher` / with-`--fisher`).
- **Float tolerance:** Fisher p-values/SBF/SSF, odds ratios, AF/AdjAF/HIAF, QUAL means within **ε = 1e-6**, after replicating R's 5-decimal **HALF_EVEN** rounding (Fisher MLE tol `ulp(1)^0.25`).
- **Bug policy (v1.0):** if v1.8.3 has an **output** bug, **replicate it** and open a tracked issue — do **not** silently "fix" it. Where v1.8.3 **crashes** (e.g. `modifyCigar`), do **not** crash: produce the same output it would on valid inputs. That is allowed and is the whole point.
- **Three commitments:** credit the original (README/`--version`/Lai et al. 2016/AstraZeneca-NGS); emulate exactly; be transparent about AI (document tools, role, and exactly what was/wasn't validated, against which version + dataset).
- **Release policy (user decision 2026-10-07 & 2026-10-08).** v1.0.0 = pure parity with v1.8.3: no output-changing fixes and no fix flags. Corrections ship from v1.1.0 per `ROADMAP.md`, each behind `--fix <id>` (default off) with a diff-confinement gate. The default output stays byte-identical to v1.0; with the flag, diffs are confined to the targeted records.

---

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
- **Rounding.** Every float written to a `.var`/VCF field goes through `round_half_even` (banker's rounding), never Rust's default. Replicate the formatter split: no-`--fisher` keeps trailing zeros (`x==0?0:fmt`); `--fisher` strips trailing zeros.
- **Sequences are bytes.** Use `&[u8]`/`Vec<u8>`, never `String`. Replace Java's negative-index `charAt`/`substr` with bounds-checked slice helpers returning `Option`.
- **Quality is Phred+33.** `Phred(byte.saturating_sub(33))`.
- **`maxReadLength` propagation** matters: SomaticMode seeds BAM2's pipeline with BAM1's value.
- **Java `HashMap` iteration that reaches output** is emulated with `vardict_core::javamap::JavaHashMap`. It ports `java.util.HashMap` spread, bin, resize and treeify-trigger semantics, and is validated against real Java by the op-log tests. This is the only sanctioned non-sorted order. Never use `std::collections::HashMap` there.

---

## 9. Idiomatic Rust + the no-panic rule
- **No panics in library code.** In `crates/vardict-core/src/**` (outside `#[cfg(test)]`) you may **not** use `unwrap`, `expect`, `panic!`, `unreachable!`, bare slice indexing `x[i]` without a prior bounds check, or arithmetic that can overflow on input. Use `Result<_, VardictError>` + `?`, `.get(i)?`, checked/saturating math. **This is what fixes the VarDictJava crash-issue family by design.**
- **`#![forbid(unsafe_code)]`** in `vardict-core`. `unsafe` only via vetted SIMD crates (e.g. `block-aligner`) behind safe APIs.
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
- Sort only within a tile (no global sort). Multi-thread BGZF. AVX2 (`block-aligner`) for the realignment hot path — but only after correctness is established, and re-validate `compare.sh --all` after.

---

## 11. Which VarDictJava issues we fix (v1.0 = robustness only)
v1.0 fixes **only** issues that change **no valid-input output** (so v1.0 stays byte-comparable to v1.8.3 *and*
is strictly more robust): the `modifyCigar` crash family (#341/#243/#253/#283/#7/#38/…), NPE/AIOOBE/div-by-zero,
OOM (#197/#384 — fixed by the tiling in §10), CRAM (#249/#380 — noodles + `-G` decode), locale-independent
formatting, fail-fast on missing index (#187), lenient SAM validation, nondeterminism (#282/#181/#356), and
correct `##contig` lengths from the `.fai` (#320). **Output-changing** corrections (#92 VD>DP, #62 over-filtering,
#363 realignment drops, #292 SHIFT3, #342 `<dup-N>` ALT, …) and **new features** (UMI, long reads) are deferred
to post-1.0 behind individual opt-in flags. See ROADMAP.md and docs/upstream-bugs.md. Post-1.0 = v1.1.0+, governed by ROADMAP.md.

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
**Setup (run once):** see `CONTRIBUTING.md` (Prerequisites) — rustup, miniconda, samtools, JDK, and reference-genome setup. Confirm every version prints before starting development.
```bash
# Build / test / lint / format — ALL must pass before any commit (zero warnings):
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo build --release

# Reference oracle (VarDictJava 1.8.3 via bioconda, env "vardict-oracle"):
bash harness/setup_oracle.sh
bash harness/run_oracle.sh <ref.fa> <bam|tumor.bam|normal.bam> <bed> <out.vcf> [sample]

# Differential testing (the acceptance gate):
bash harness/compare.sh <ref> <bam> <bed> [sample]   # Tier-2: full VCF vs oracle  → prints MATCH/MISMATCH
python3 harness/compare_var.py golden/<case>.var <ours.var>   # Tier-1: intermediate table

# Slice big BAMs for the dev loop:
bash harness/slice.sh <in.bam> <chr:start-end> <out.bam>
```

## 15. Definition of done (every change)
1. `cargo fmt --all`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test` all pass.
2. The relevant `harness/compare.sh` prints **MATCH** (or the named test passes).
3. `golden/` updated if output legitimately changed (with provenance).
4. The commit/PR states: **"validated vs VarDictJava 1.8.3, dataset `<name>`, cmd `<…>`"**.

## 16. Repo layout
```
crates/vardict-core/   engine library (all logic; #![forbid(unsafe_code)])
crates/vardict-cli/    bin "vardictrust" (clap front-end → core)
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
Atomic commits per task: `git add -A && git commit -m "Txxx: <short desc>"`. Branch off `main`; open PRs via `gh`.
Never commit secrets, **data** (BAM/CRAM/FASTA/VCF beyond the tiny `golden/` fixtures), `target/`, or conda envs — see `.gitignore` (a single WGS BAM can exceed 150 GB). **Do** commit `Cargo.lock` (this is an application workspace). MIT-licensed, open-source project. No `Co-Authored-By` commit trailers per user directive (A10); AI involvement is disclosed transparently in `README.md` and `CITATION.cff`.
