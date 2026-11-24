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

### VarDictRust Dual-Mode Architecture (`--compat=vardictjava` vs `--compat=off`)

Under the **rewrites.bio** contract and the zero-panic rule ([`AGENTS.md §9`](../../AGENTS.md#9-idiomatic-rust--the-no-panic-rule)), VarDictRust replaces unchecked Java arithmetic with safe, saturating math and explicit `Result<ModifiedCigar, VardictError>` error handling across both compatibility modes:
1. **Default Parity Mode (`--compat=vardictjava` / `--parity`):** Malformed `JavaThrows` CIGAR strings are safely intercepted and skipped before pileup generation without panicking and **without** upstream Java's arbitrary 10-exception abort ceiling (`validation_error_rate_threshold_enforced`), preserving 100% valid-read calling fidelity.
2. **Opt-In Correctness Mode (`--compat=off`, `CompatProfile::cigar_homopolymer_merge`):** Repairs upstream `PR #402` homopolymer multi-indel CIGARs (merging `104M1I2M2D1M1D5M1D2M9D27M` into `104M14D11I9D27M` with strict positive operator lengths and exact `(query, ref)` length conservation) so matched-sequence reads rescue their underlying complex/deletion variants (`vcf_all_alleles_and_compat_off`), while still safely skipping unrecoverable corrupt records without panicking.

For architectural context and upstream bug taxonomy, see:
- [`DIVERGENCES.md` (Category 1: CIGAR & Realignment Resilience)](../../DIVERGENCES.md#category-1-cigar-and-realignment-resilience)
- [`docs/upstream-bugs.md`](../../docs/upstream-bugs.md)
- [`ROADMAP.md` (Open Upstream PRs Reviewed)](../../ROADMAP.md#open-upstream-prs-reviewed-v100-baseline)

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

This fixture is actively exercised by the engine and CIGAR robustness test suites in [`crates/vardictrust-core/tests/engine_behavior.rs`](../../crates/vardictrust-core/tests/engine_behavior.rs) and [`crates/vardictrust-core/src/cigar/modify.rs`](../../crates/vardictrust-core/src/cigar/modify.rs):

```rust
// crates/vardictrust-core/tests/engine_behavior.rs
#[test]
fn malformed_cigar_returns_clean_result() { ... }
```

### Invariants Verified:
1. **No-Panic Guarantee:** The caller never panics, aborts, or exits with a fatal signal when encountering malformed CIGAR strings under either `--compat=vardictjava` or `--compat=off`.
2. **Default Parity Read Skipping (`--compat=vardictjava`):** Corrupted `JavaThrows` reads are counted and skipped without mutating pileup state for surrounding reads, emitting identical variant calls on 1 versus 11 duplicated corrupted reads (surpassing upstream Java's 10-read crash ceiling).
3. **Opt-In Homopolymer Multi-Indel Repair (`--compat=off`):** Matched-query `PR #402` alignments (`104M1I2M2D1M1D5M1D2M9D27M`) coalesce cleanly into `104M14D11I9D27M` (`cigar_homopolymer_merge`) with exact `(query, ref)` span conservation, rescuing true complex and deletion variants.
