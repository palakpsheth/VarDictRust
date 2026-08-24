# VarDictRust — Engineering Kickoff Prompt (hand this to a coding agent)

You are an autonomous coding agent picking up the **VarDictRust** project. Your job is to build it by
executing an existing, fully-specified plan. **Before writing any code, read `CLAUDE.md` and `PLAN.md` (in the
repo root) end to end** — they are the source of truth; this prompt only orients you.

---

## Mission
Rewrite **VarDictJava v1.8.3** (a somatic/germline NGS variant caller) from scratch in **idiomatic Rust** as a
**single, statically-compiled binary** that **scales to deep WGS** and produces output **identical to VarDictJava
v1.8.3 on valid inputs**. This follows the **rewrites.bio** manifesto: *a faster, more robust tool that produces
the same results* — **not** a new caller, **not** a line-by-line transliteration of the Java.

## Environment (already set up)
- Repository root (this machine): `/Users/palak.sheth/PycharmProjects/VarDictRust`; remote `git@github.com:palakpsheth/VarDictRust.git`.
- Dev/CI target: **Ubuntu 12c / 24 GB / amd64** → SIMD is **AVX2** (not Neon). Static target `x86_64-unknown-linux-musl`.
- Reference implementation is cloned at **`./VarDictJava`** (v1.8.3, commit `0886261`); the Perl/R submodule is checked out at **`./VarDictJava/VarDict`** (commit `009e017`). This is your **source of truth for every port** — read the Java, cite `file:line` in doc-comments.
- Test fixtures already on disk: **`./VarDictJava/testdata/integrationtestcases`** (94 cases; filenames encode `mode;ref;bam;region;flags`) and **`./VarDictJava/tests/integration`** (`expected/vardict.{simple,paired}.vcf` + `raw_input/raw.vardict.{simple,paired}.var` intermediate tables).

## How to work (this is the whole method)
0. **Set up first:** run `PLAN.md` Appendix **B0 (Prerequisites)** — confirm `rustc`/`cargo`/`conda`/`samtools` all print versions, and obtain any reference genome a fixture needs (B0 says which cases use the bundled tiny fastas vs. need `hg19`/`GRCh37`).
1. Read `CLAUDE.md` (the agent contract) and `PLAN.md` (strategic reference + **Appendix C: a 174-task execution manual**). New to the genomics terms or the `+ATG`/`-3` variant-description grammar? Skim `PLAN.md` Appendix **D** first.
2. **Execute `PLAN.md` Appendix C tasks T001 → T174, strictly in order.** Each task is self-contained:
   **Read** (which Java `file:lines` to study) → **Do** (exact Rust file + signature to create) → **Cmd** (how to validate) → **Expect** (what proves success). `PLAN.md` **B8** shows one task (T016) worked end-to-end — imitate its shape exactly.
3. **Never skip a validation.** A task is done only when its `Expect` passes (a `harness/compare.sh` **MATCH** or the named `cargo test`). Then commit: `git add -A && git commit -m "Txxx: <desc>"`.
4. After **every** code change run, in order, all passing with zero warnings:
   `cargo fmt --all` → `cargo clippy --all-targets --all-features -- -D warnings` → `cargo test`.
5. **The port loop (rewrites.bio "work small"):** pick ONE function → read its Java → port it + a unit test → run `compare.sh` → **MATCH = commit & next**; **MISMATCH = debug THIS function only, change nothing else.** Use the Tier-1 `.var` diff to localize a mismatch to the engine before the VCF stage.

## Locked decisions — do NOT re-litigate
- **Scope:** full parity (Simple, Somatic, Amplicon, Splicing/SV), equal depth; build the shared engine first.
- **Crate layout:** Cargo workspace — `crates/vardict-core` (engine lib, `#![forbid(unsafe_code)]`) + `crates/vardict-cli` (bin `vardictrust`).
- **I/O:** **noodles** (pure Rust; no C toolchain). Pin exact versions.
- **Reference oracle:** **bioconda** `vardict-java=1.8.3` (env `vardict-oracle`); `bash harness/setup_oracle.sh`.
- **CLI:** drop-in + ergonomic extras — all VarDictJava short flags/defaults **plus** long aliases/subcommands; **direct-to-VCF** (R+Perl folded in); `--emit-intermediate` reproduces the `.var`.
- **Distribution:** open-source **MIT**, GitHub release binaries via `dist`, publish to **crates.io via OIDC trusted publishing**.
- **Equivalence bar:** structural-exact + float tolerance **ε = 1e-5** after 5-dp **HALF_EVEN** rounding.
- **Issues:** v1.0 fixes **only robustness/crash/input/IO/memory/systems** issues (no valid-input output change). Output-changing fixes + new features are **post-1.0, opt-in**. (`PLAN.md` §8.)

## Non-negotiable rules (full detail in `CLAUDE.md`)
- **Integrity — never fake a pass (`CLAUDE.md` §0b).** Never relax `--eps`, narrow `--cols`, edit/skip/`#[ignore]`/delete a failing test or anything in `harness/`, hand-edit `golden/`/`expected/`, or report a result you didn't actually run. If a check fails, the **code** (or your understanding) is wrong — fix the port, not the test. Can't pass it honestly? **Stop and report the MISMATCH with the diff.**
- **Emulate exactly.** Match v1.8.3 output on valid inputs: variant set, FILTER, INFO/FORMAT keys, GT, TYPE, STATUS, and column counts/order. Where v1.8.3 **crashes**, do **not** crash — produce the output it would on valid inputs. If v1.8.3 has an **output** bug, replicate it and open a tracked issue (don't silently "fix").
- **Coordinates are typed.** BAM = 0-based → `Pos0`; VCF/VarDict = 1-based → `Pos1`; convert only via `.to_one_based()`/`.to_zero_based()`. Unsure of a value's basis? **STOP and re-read the Java.** Never guess.
- **No panics in `vardict-core`** (outside `#[cfg(test)]`): no `unwrap`/`expect`/`panic!`/unchecked indexing/overflowing math. Use `Result` + `?`, `.get(i)?`, checked math. *(This is what fixes the VarDictJava crash-issue family by design.)*
- **Fidelity hazards:** use `IndexMap`/`BTreeMap` (never `HashMap`) in output/realignment paths; `round_half_even` (never Rust default); sequences are `&[u8]`/`Vec<u8>` (never `String`); quality is Phred+33.
- **Scale-first:** process region **tiles with a halo**; bounded memory at any depth; stream + drop reads; `rayon` work-stealing; deterministic ordered output; tiled output must equal a single-region run. Don't optimize before correctness; re-run `compare.sh --all` after any perf change.
- **Evidence before completion.** "It should work" is not evidence. Only a MATCH / passing test counts.
- **Unsupported flags error loudly** (`VardictError::Unsupported`), never silently ignored.

## Critical gotchas (the source blueprint and a naive reading get these WRONG)
1. **Intermediate column schemas are the Java printer emissions: Simple 36/38, Amplicon 38/40, Somatic 55/61** (no-`--fisher`/with-`--fisher`). The R-script gate numbers (34/36/38, ≥48) are **input validators only** — do **not** use them as emission counts. Exact field orders: `PLAN.md` §4.3.
2. **Port the Fisher test from Java**, not from R: `data/fishertest/FisherExact.java` + `UnirootZeroIn.java`. Two-sided p-value = sum of null probs ≤ observed·(1+1e-7); odds ratio = **conditional MLE** via Brent root-find (tol `ulp(1)^0.25`). The `fishers_exact` crate gives **p-value only** (use it as a cross-check) — you must compute the odds ratio yourself.
3. **Somatic STATUS is computed in Java** (`SomaticPostProcessModule.determinateType`, with `isNoise`/`combineAnalysis`), **not** in Perl. Thresholds in `PLAN.md` §4.5.
4. The **highest-fidelity-risk** module is `VariationRealigner` (mutation- and iteration-order-sensitive); the most **crash-prone** is `CigarModifier.modifyCigar`. Port these one sub-function at a time, validating each.
5. **Do NOT add these crates** (hallucinated/unfit): `genomicframe-core`, `rustalign-simd`, `ksw2rs`, `seqair`, `cyanea-align`. Use noodles / bio / (later) block-aligner / statrs / fishers_exact.
6. **VarDictJava output ordering is nondeterministic on some inputs** (issue #282; `VarDictJava/tests/integration/test_non_deterministic_behavior`). `compare_vcf.py` aligns records by key, so a pure-ordering difference is **not** a mismatch — never change the engine to chase one.

## Validation (harness scripts + comparators are provided complete in `PLAN.md` B4)
- Oracle: `bash harness/setup_oracle.sh` then `bash harness/run_oracle.sh ...`.
- **Tier-1** (engine → intermediate `.var`): `python3 harness/compare_var.py golden/<case>.var <ours.var>` — localizes a mismatch to the engine, before R/Perl.
- **Tier-2** (full → VCF, the acceptance gate): `bash harness/compare.sh <ref> <bam|tumor.bam|normal.bam> <bed>` → prints `MATCH`/`MISMATCH`.
- Slice big BAMs for the dev loop: `bash harness/slice.sh <in.bam> <chr:start-end> <out.bam>`.

## Milestones (see `PLAN.md` §14 + Appendix C)
M0 scaffolding (T001–T014) → **M1 single-sample parity** (T097–T108) → **M2 somatic** (T109–T124) → amplicon (T125–T134) → **M3 full parity/SV** (T135–T146) → perf/scale (T147–T154) → robustness corpus + CI/release (T155–T166) → validation + **v1.0** (T167–T174).
**v1.0 complete** = M3 reached (all four modes MATCH on the testdata corpus + SEQC2/GIAB slices) **and** all §8 FIX-BY-DESIGN robustness items covered (crash corpus green, CRAM works, bounded memory verified) **and** the §15 acceptance sweep green **and** released under MIT to GitHub + crates.io.

## When to stop and ask the human
- Java is ambiguous about coordinate basis or an algorithm's intent.
- An output mismatch you cannot localize even with the Tier-1 `.var` diff (report it with the diff).
- A decision arises that `PLAN.md`/`CLAUDE.md` don't cover (don't assume — ask).

## What NOT to do
Don't transliterate Java idioms (write idiomatic Rust). Don't add hallucinated crates. Don't fix output-changing
VarDictJava bugs in v1.0 (replicate + file an issue; fixes are post-1.0 opt-in). Don't claim done without
evidence. Don't optimize before correctness.

---

**Begin now:** read `CLAUDE.md` and `PLAN.md` in full, run `PLAN.md` Appendix **B0** (confirm `rustc`/`cargo`/`conda`/`samtools` print; fetch any reference genome a fixture needs), then start at **task T001** in `PLAN.md` Appendix C and proceed in order, committing after each green task. **An honest MISMATCH beats a fake MATCH.**
