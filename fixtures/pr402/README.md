# VarDictJava PR #402 / Issue #341 Crash Test Fixture

This directory contains the synthetic regression fixture derived from upstream **VarDictJava PR #402** and **Issue #341** ("Skip read on CigarModifier / modifyCigar exception").

---

## 1. Background & Upstream Failure Mode

In VarDictJava v1.8.3, certain complex CIGAR strings containing closely spaced indels, deletions, and homopolymer tracts trigger unhandled exceptions inside `CigarModifier.modifyCigar`.

### Homopolymer CIGAR Corruption Math

When `modifyCigar` iteratively parses and adjusts alignment boundaries around complex indels, it calculates adjusted match lengths by subtracting indel spans from adjacent match blocks. For complex read alignments such as:

```text
104M1I2M2D1M1D5M1D2M9D27M
```

intermediate match segments (e.g. `2M`, `1M`) separated by successive deletions (`2D`, `1D`, `1D`, `9D`) cause arithmetic underflow when subtracting realignment offsets. This produces invalid negative match lengths (e.g. `-1M`). When VarDictJava constructs the modified CIGAR string with negative operators or slices sequence buffers using negative lengths, it triggers `StringIndexOutOfBoundsException` or `NegativeArraySizeException`.

### The 10-Exception Fatal Abort Ceiling

Upstream PR #402 caught exceptions in `modifyCigar` and skipped the offending read. However, upstream Java instituted a hardcoded abort threshold: if **more than 10** such exceptions occur during a run, VarDictJava re-throws the exception and fatally terminates the entire pipeline (`SIGABRT` / JVM crash), discarding all downstream variant calls. In deep target sequencing or clinical amplicon panels with high coverage across repetitive regions, encountering >10 noisy reads is routine, causing reproducible pipeline crashes.

### VarDictRust Crash-Safety Architecture

Under the **rewrites.bio** contract and the zero-panic rule ([`AGENTS.md §9`](../../AGENTS.md#9-idiomatic-rust--the-no-panic-rule)), VarDictRust replaces unchecked Java arithmetic with safe, saturating math and explicit `Result<ModifiedCigar, VardictError>` error handling:
1. **Crash Prevention:** Malformed CIGAR strings are safely intercepted before pileup generation.
2. **No Arbitrary Abort Limit:** VarDictRust does not abort after 10 corrupted reads; it reliably skips any number of corrupted reads while maintaining full caller throughput and variant calling fidelity for all valid alignments.

For architectural context and upstream bug taxonomy, see:
- [`DIVERGENCES.md § UB-1`](../../DIVERGENCES.md)
- [`docs/upstream-bugs.md § UB-1`](../../docs/upstream-bugs.md)
- [`ROADMAP.md`](../../ROADMAP.md)

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

This fixture is actively exercised by the integration test suite in [`crates/vardictrust/tests/robustness.rs`](../../crates/vardictrust/tests/robustness.rs):

```rust
// crates/vardictrust/tests/robustness.rs:203
#[test]
fn test_pr402_11_bad_reads_does_not_abort() { ... }
```

### Invariants Verified:
1. **No-Panic Guarantee:** The caller never panics, aborts, or exits with a fatal signal when encountering malformed CIGAR strings.
2. **Selective Read Skipping:** Corrupted reads are logged and skipped without mutating pileup state for surrounding reads.
3. **Output Invariance & Scalability:** The caller emits identical variant calls when processing 1 corrupted read versus 11 duplicated corrupted reads (surpassing upstream Java's 10-read crash ceiling).
