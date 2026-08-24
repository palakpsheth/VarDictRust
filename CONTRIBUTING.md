# Contributing to VarDictRust

Read [CLAUDE.md](CLAUDE.md) first — it is the project contract (prime directive,
fidelity hazards, integrity rules). This file is the quick-start.

## Build / lint / test

All of these must pass (zero warnings) before any commit:

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
cargo build --release
```

Library code in `crates/vardict-core` forbids `unsafe`, `unwrap`/`expect`,
panics, and unchecked indexing outside tests — the no-panic rule is what fixes
the VarDictJava crash family by design. Errors flow through
`Result<_, VardictError>` (`thiserror`) in core, `anyhow` in the CLI.

## Validation harness

The acceptance gate is **differential testing against VarDictJava v1.8.3** —
not "it looks right". One-time setup (rustup, miniconda, samtools, JDK,
reference genomes): see PLAN.md Appendix B0; oracle conda env:

```bash
bash harness/setup_oracle.sh        # builds env "vardict-oracle" (VarDictJava 1.8.3 via bioconda)
bash harness/run_oracle.sh <ref.fa> <bam> <bed> <out.vcf> [sample]   # reference output
```

Development loop:

```bash
bash harness/run_integration.sh     # Tier-0: corpus integration suite
bash harness/compare_var.py golden/<case>.var ours.var   # Tier-1: intermediate .var table
bash harness/compare.sh <ref> <bam> <bed> [sample]       # Tier-2: full VCF vs live oracle → MATCH/MISMATCH
bash harness/slice.sh <in.bam> <chr:start-end> <out.bam>  # slice big BAMs for fast iteration
```

A change is done when the relevant `compare.sh` prints **MATCH** (or the named
test passes), not before. Evidence before completion: never claim a pass you
did not observe.

## Integrity rules

Never weaken a gate to make it pass. Do not edit or skip failing tests, relax
comparators, narrow ranges, or hand-edit `golden/` outputs. If a check fails,
your code (or your understanding) is wrong — fix the port, not the test.
Full rules: [CLAUDE.md](CLAUDE.md) §0b. A real MISMATCH reported honestly is
worth more than a fake MATCH.

## Porting workflow

Work small: pick ONE VarDictJava function → read its Java source and cite it
(`/// ports CigarParser.java:2295-2341`) → port + unit test → run the relevant
harness check → commit. Never batch multiple functions without per-function
evidence. Coordinate types matter: BAM positions are 0-based half-open (`Pos0`),
VCF/VarDict logic 1-based (`Pos1`); convert only via `.to_one_based()` /
`.to_zero_based()`.

## Commits

Atomic commits per task, message prefixed with the task id:
`T123: short description`. Branch off `main`; open PRs via `gh`. Never commit
secrets, sequencing data (BAM/CRAM/FASTA beyond tiny `golden/` fixtures), or
build artifacts; do commit `Cargo.lock`.

## License

By contributing you agree your contributions are licensed under the MIT
License ([LICENSE](LICENSE)).
