# JavaMap JVM Oracle Operation Logs & Provenance

This directory contains the JVM test oracle generator and operation trace logs used to validate [`crates/vardictrust-core/src/javamap.rs`](../../src/javamap.rs).

---

## 1. Background: Eliminating Cross-Position Iteration Divergence ([UB-3](../../../../docs/upstream-bugs.md#ub-3))

In VarDictJava, internal variant collections across position keys (`nonInsertionVariants`, `insertionVariants`, `alignedVariants`) are stored in standard `java.util.HashMap` instances.

### Biological and Algorithmic Motivation

In deep NGS sequencing datasets or multiallelic somatic sites, multiple sub-clonal variants or indel candidates can occur at identical genomic coordinates with indistinguishable empirical allele frequencies (e.g. at `chr20:20202875` on dataset `s1a`, where candidate alleles exhibit identical $AF = 0.5000$). When the downstream pipeline script `var2vcf_valid.pl` filters and reformats the intermediate `.var` table into standard VCF 4.2 format, it encounters candidate lines sequentially. If two variants tie on allele frequency, `var2vcf_valid.pl` selects the first record encountered in the stream, silently discarding the latter.

Consequently, **which somatic allele survives in the final VCF depends strictly on `java.util.HashMap`'s internal iteration order**.

Standard Rust hash maps (`std::collections::HashMap`, `hashbrown`) use SipHash or SwissTable algorithms with randomized seeds, producing iteration orders different from Java. Using `BTreeMap` or sorted order resolves determinism, but produces output diverging from VarDictJava v1.8.3 on valid clinical inputs.

To achieve 100% parity under the [rewrites.bio](https://rewrites.bio) manifesto, VarDictRust implements `vardictrust_core::javamap::JavaHashMap`, faithfully replicating OpenJDK `java.util.HashMap` mechanics:
- Hash spreading: `(h = key.hashCode()) ^ (h >>> 16)`.
- Power-of-two bucket table indexing: `(capacity - 1) & hash`.
- Standard load factor (0.75) and capacity doubling.
- Node collision chaining and treeification threshold ($\ge 8$ items per bin).

For architectural context and upstream bug taxonomy, see:
- [`DIVERGENCES.md § UB-3`](../../../../DIVERGENCES.md#ub-3-javahashmap-iteration-and-allele-frequency-tie-breaking)
- [`docs/upstream-bugs.md § UB-3`](../../../../docs/upstream-bugs.md#ub-3)
- [`ROADMAP.md`](../../../../ROADMAP.md) (post-1.0 opt-in deterministic tie-breaking flag `--fix javamap-deterministic`)

---

## 2. Trace Files & Provenance

The operation logs in this directory were generated using the standalone Java harness [`JavaMapOracle.java`](JavaMapOracle.java):

```bash
# Compilation and execution with host OpenJDK
javac JavaMapOracle.java
java JavaMapOracle <seed> <numOps> <base> <span> <stride>
```

| File | Size / Operations | Description |
| :--- | :--- | :--- |
| [`JDK_VERSION`](JDK_VERSION) | `openjdk version "25.0.2-internal" 2026-01-20` | Reference JVM environment record. Note that OpenJDK's `java.util.HashMap` core hashing, bucket indexing, and treeification algorithms have remained invariant from Java 8 through Java 25. |
| [`JavaMapOracle.java`](JavaMapOracle.java) | 17 lines | Pure-Java harness that applies pseudo-random `put`, `computeIfAbsent`, and `remove` operations to a `HashMap<Integer, Integer>`, dumping the final `keySet()` iteration order. |
| [`ops_1.txt`](ops_1.txt) | 3,001 lines | Replay trace testing sequential inserts, lookups, and deletions. |
| [`ops_2.txt`](ops_2.txt) | 20,001 lines | High-volume replay trace exercising multiple table resizes and capacity growth phases. |
| [`ops_treeify.txt`](ops_treeify.txt) | 401 lines | Targeted collision sequence asserting correct bin collision threshold detection and treeify flag activation when bin count exceeds 8. |

---

## 3. Test Verification in Rust

These logs are compiled into the test binary via `include_str!` and replayed by the unit test suite in [`crates/vardictrust-core/src/javamap.rs`](../../src/javamap.rs):

```rust
// crates/vardictrust-core/src/javamap.rs:300
#[test]
fn test_javamap_oracle() {
    // Replays ops_1.txt, ops_2.txt: asserts exact keySet() sequence matching Java's output.
    // Replays ops_treeify.txt: asserts detection of collision threshold and telemetry flag activation.
}
```

This verification ensures that `JavaHashMap` produces iteration orders 100% bitwise-identical to the JVM reference implementation across all normal operational regimes, with active telemetry guards for treeified regimes (tracked post-1.0 in [`ROADMAP.md`](../../../../ROADMAP.md)).
