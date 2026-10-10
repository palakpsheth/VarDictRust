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
| ub2 | DIVERGENCES UB-2 | `-p` null-record row empty GENO and "0;0" bias | emit standard null placeholders or skip row |
| ub5 | DIVERGENCES UB-5 | `var2vcf_valid.pl` drops `<DUP>`/`<INV>` when `splitreads < 1` | preserve discordant-pair SVs lacking split reads |
| ub6 | DIVERGENCES UB-6 | `markDUPSV` pre-decrements `pe`, truncating boundary by 1bp | pass pre-decrement `pe + 1` window boundary |
| ub7 | DIVERGENCES UB-7 | single-winner multi-allelic collapsing drops secondary alleles | retain secondary qualifying alleles in VCF |

## Needs a specification first (debatable heuristics; no flag until specified)
#92/#386/#323/#316/#400 (VD > DP depth) · #62 (over-filtering) · #363 family (realignment drops) ·
#381 · #282 (soft-clip) · #264 · HICOV excludes insertions (ToVarsBuilder.java:523-538).

## Open upstream PRs reviewed (2026-10-08)
#402 → row 402. #371, #372, #335 → no output impact on VarDictRust (N/A).

## Post-1.0 features
UMI support, long reads — separate design docs.
