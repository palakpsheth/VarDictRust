# VarDictJava PR #402 / Issue #341 Crash Test Fixture

This directory contains the synthetic regression fixture derived from upstream **VarDictJava PR #402** and **Issue #341** ("Skip read on CigarModifier / modifyCigar exception").

---

## 1. Background & Upstream Failure Mode

In VarDictJava v1.8.3, certain complex CIGAR strings (specifically those containing sequences of insertions, deletions, and clips adjacent to matches) triggered `ArrayIndexOutOfBoundsException` or `StringIndexOutOfBoundsException` inside `CigarModifier.modifyCigar`.

Historically, when such an exception occurred during multi-threaded processing, unhandled exceptions caused the entire pipeline process to terminate abnormally (`SIGABRT` or JVM crash), discarding partially called variants. Upstream PR #402 addressed this by catching exceptions in `modifyCigar` and dropping the individual offending read rather than crashing the caller.

---

## 2. Fixture Structure

The fixture represents a 600 bp region on `chr7` designed to test caller resilience against malformed CIGAR inputs:

- [`ref.fa`](ref.fa), [`ref.fa.fai`](ref.fa.fai): 600 bp synthetic FASTA reference on `chr7` (file MD5: `0ec0bcf8abce297393e3c0cce5acabcb`).
- [`pr402.bed`](pr402.bed): Single-target interval (`chr7:1-600` targeting `GeneA`).
- [`pr402.sam`](pr402.sam), [`pr402.bam`](pr402.bam), [`pr402.bam.bai`](pr402.bam.bai):
  - 30 standard reads aligned to `chr7:1` with simple `142M` CIGARs.
  - 1 adversarial read (`bad_read`) with complex CIGAR:
    ```
    104M1I2M2D1M1D5M1D2M9D27M
    ```

---

## 3. Active Usage & Verification

This fixture is actively exercised by the integration test suite in [`crates/vardict-cli/tests/robustness.rs`](../../crates/vardict-cli/tests/robustness.rs):

```rust
// crates/vardict-cli/tests/robustness.rs:203
#[test]
fn test_pr402_11_bad_reads_does_not_abort() { ... }
```

### Invariants Verified:
1. **No-Panic Guarantee:** The caller never panics, aborts, or exits with a fatal signal when encountering malformed CIGAR strings.
2. **Selective Read Skipping:** The corrupted read is logged and skipped without mutating pileup state for surrounding reads.
3. **Output Invariance:** The caller emits identical variant calls when processing 1 corrupted read versus 11 duplicated corrupted reads.
