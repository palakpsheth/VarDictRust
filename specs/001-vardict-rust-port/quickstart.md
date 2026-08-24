# Quickstart & Verification Guide: VarDictRust Port

This document details the runnable validation scenarios to verify output parity between VarDictRust and the reference VarDictJava oracle.

## Prerequisites

Ensure all dependencies, test datasets, and reference genomes are set up as per the project requirements:
1. Active conda environment `vardict-oracle` with `vardict-java=1.8.3` installed.
2. Compiled `vardictrust` binary under target folder.

## Parity Verification Scenarios

### Scenario 1: Single-Sample (Germline) Verification
Validate simple/germline mode output parity:

```bash
# 1. Run the VarDictJava oracle
bash harness/run_oracle.sh testdata/ref.fa testdata/sample.bam testdata/regions.bed expected_simple.vcf

# 2. Run VarDictRust
cargo run --release -- -G testdata/ref.fa -b testdata/sample.bam -N SAMPLE testdata/regions.bed > ours_simple.vcf

# 3. Perform two-tier validation
# Tier-1: Verify intermediate table parity
cargo run --release -- --emit-intermediate -G testdata/ref.fa -b testdata/sample.bam testdata/regions.bed > ours_simple.var
python3 harness/compare_var.py golden/expected_simple.var ours_simple.var

# Tier-2: Verify VCF equivalence
bash harness/compare.sh testdata/ref.fa testdata/sample.bam testdata/regions.bed
```

**Expected Outcome**: Comparative scripts output `MATCH`.

---

### Scenario 2: Somatic (Paired Tumor/Normal) Verification
Validate somatic mode output parity:

```bash
# 1. Run the VarDictJava oracle in somatic mode
bash harness/run_oracle.sh testdata/ref.fa "tumor.bam|normal.bam" testdata/regions.bed expected_somatic.vcf

# 2. Run VarDictRust
cargo run --release -- -G testdata/ref.fa -b "tumor.bam|normal.bam" -N "TUMOR|NORMAL" testdata/regions.bed > ours_somatic.vcf

# 3. Perform validation
bash harness/compare.sh testdata/ref.fa "tumor.bam|normal.bam" testdata/regions.bed
```

**Expected Outcome**: Output reports `MATCH`.

---

### Scenario 3: Quality Gates Check
Before committing code, verify formatting, linting, and internal tests pass:

```bash
cargo fmt --all
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

**Expected Outcome**: Exit code `0` for all validation commands, with no warnings.
