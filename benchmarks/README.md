# Benchmarks

Engine hot-path benches (criterion, `crates/vardict-core/benches/engine.rs`).
Fixtures are constructed in-code — no committed data.

## Run

```bash
cargo bench            # full run, report under target/criterion/
cargo bench --no-run   # compile-only gate
```

## What is measured

| Bench | Target | Fixture |
|---|---|---|
| `parse_cigar/100_reads_mixed_cigars` | `vardict_core::cigar::parse::pileup` (the `CigarParser::parse_cigar` loop incl. per-read `CigarModifier`, realign=1) | 100 coordinate-sorted 150 bp reads over a 20 kb synthetic ACGT reference; CIGAR mix M / S / D / I |
| `find_bp/60bp_consensus_indel40` | `vardict_core::realign::find_bp` (breakpoint scan inner loop of `realign_indels`) | 60 bp synthetic consensus (~10% mismatches), indel size 40, backward direction |

## Baseline (one local run, 2026-08-23, criterion 0.5, `cargo bench`)

```text
parse_cigar/100_reads_mixed_cigars
                        time:   [2.1352 ms 2.1621 ms 2.1932 ms]
                        (≈ 21.6 µs/read)

find_bp/60bp_consensus_indel40
                        time:   [5.4194 µs 5.4864 µs 5.5924 µs]
```

Environment caveat: single developer machine, no pinned hardware — treat as
order-of-magnitude baseline only; re-measure on target hardware before
accepting/rejecting optimization PRs.

## Scaling gates (deferred post-baseline)

Per the plan, AVX2 / `block-aligner` realignment acceleration and full tiling-scale
sweeps are **explicitly deferred** until after this baseline exists. These two
benches are the before/after gates for that work: any SIMD or tiling change
must keep `harness/run_integration.sh` at 91 PASS / 1 known-fail AND show its
effect here. Lite scale evidence below.

## Memory (T150-lite, 2026-08-23)

Peak RSS via `/usr/bin/time -v`, release binary, one NA24385 chr20 slice:
**chr20:20,525,032–20,604,921** (~80 kb BED region, **mean depth 96x, max 214x**),
`-c1 -S2 -E3 -g4 -f 0.001`. Raw logs + sample VCF:
`.superpowers/sdd/2026-08-21-parity-to-v1.0/p10-scale/`.

| `-th` | Peak RSS | Wall |
|---|---|---|
| 1 | 434,944 KB (~425 MB) | 8.93 s |
| 4 | 435,532 KB (~425 MB) | 3.33 s |
| 20 (all cores) | 435,908 KB (~425 MB) | 3.39 s |

Plan target (`PLAN.md`): peak memory **< 4 GB per
thread** — observed ~425 MB total regardless of thread count, i.e. per-thread
footprint stays far under target. RSS is flat from 1→20 threads, consistent
with §10 bounded-memory tiling (per-tile state only) at this depth; VCF output
was byte-identical across all three configs.

Not yet claimed: this is one region at moderate depth on one machine — **no
WGS-depth soak, no high-depth (>1000x) amplicon, no pinned hardware**, so
"bounded at any depth" remains an architectural claim, not a measured one.

## Scaling determinism (T151/T152-lite)

Guard test: `cargo test -p vardict-cli --test thread_determinism`
(`crates/vardict-cli/tests/thread_determinism.rs`). Runs vardictrust on a
synthetic multi-variant BAM over a 1 kb region forced into several tiles via
`--tile-size=300`, at `-th` 1 / 2 / 8 (capped by host cores), and asserts the
VCF bytes are identical across runs (§10 deterministic ordered output).
Observed runtime ≈ 0.2 s — far inside the <60 s budget. A regression in
ordered output or tile-boundary dedup fails the byte comparison.
