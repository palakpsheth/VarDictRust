# Upstream bugs — VarDictJava v1.8.3

Consolidated, actionable list of VarDictJava v1.8.3 behaviors that are part of the
specification under rewrites.bio §5 (replicate output bugs byte-for-byte, track
upstream, never silently fix), plus one upstream-repo defect that affects validation
infrastructure. Supersedes the drafts formerly embedded in `DIVERGENCES.md`.

**Version pin:** VarDictJava **v1.8.3 @ `0886261`**, Perl/R submodule @ `009e017`.
All line references below are to that pin.

| # | Item | Our stance |
|---|------|------------|
| [UB-1](#ub-1--p-spurious-all-zero-record-at-position-00) | `-p` spurious all-zero record at position 0/0 | **Replicated** (P8a) |
| [UB-2](#ub-2--p-null-record-row-cell-formatting-empty-geno-00-bias-empty-flank-seqs) | `-p` null-record row cell formatting | **Replicated** (P8a) |
| [UB-3](#ub-3--var-row-order-across-positions-follows-hashmap-iteration-order-282-family) | Cross-position `.var` emission order follows HashMap layout | **Not replicated** (ordered maps; ≤1-record AF-tie residue) |
| [UB-4](#ub-4--integration-fixture-inconsistency-overlapping-csv-reference-chunks) | Integration fixtures encode mutually inconsistent references | N/A (test-data defect; validated vs live oracle) |

---

## UB-1 — `-p`: spurious all-zero record at position 0/0

**Java behavior:** `SimplePostProcessModule.java:73-82` — for an `alignedVariants`
entry whose single variant has `startPosition != position`, with pileup on and
`referenceVariant == null`, Java prints a null-wrapper row immediately (:76-77)
and then continues processing with `refVar = new Variant()` (:78). That default
`Variant` keeps `startPosition = endPosition = 0`; it is added to `vrefs` and
printed at :100-101 via `SimpleOutputVariant`'s non-null path, whose ternaries map
every zero/default field to `"0"` (`SimpleOutputVariant.java:64/:81-82`) — leaking
the constructor sentinel into output coordinates. One bogus all-zero 36-col row at
position 0/0 per occurrence (position 0 is outside any requested region; the
downstream Perl VCF conversion receives a REF/ALT-less row at pos 0).

**Repro:** `vardict-java -G hard_clip_case.fa -f 0.0 -b hard_clip_next_to_del_test2.bam -c 1 -S 2 -E 3 -g 4 -z -p -r 1 <(printf 'test\t6433\t6564\ttestbed')`
→ stdout contains `hard_clip_next_to_del_test2	testbed	test	0	0	…` between the
6507-6533 and 6449-6468 rows.

**Our replication stance:** replicated byte-for-byte (P8a, commit `4d113b2`;
`postprocess/simple.rs` fallback branch keeping `Variant::new()` defaults +
`output/table.rs`). Verified sorted-diff-empty vs live oracle on the hard_clip case.

**Draft issue text:**

> **Title:** `-p` mode prints an all-zero record at position 0/0 when a position's only variant's raw start differs from its map key and no reference variant exists
>
> In v1.8.3, `SimplePostProcessModule.java:73-82`: when a position's single variant
> has `startPosition != position` and `referenceVariant == null`, a default
> `new Variant()` (:78) enters `vrefs` and is printed at :100-101 with
> `startPosition = endPosition = 0`, emitting an all-zero `.var` row at position 0
> outside the requested region. Fix direction: seed :78 from the map key (or skip
> printing the sentinel Variant). Output-changing, so it needs a release note.

---

## UB-2 — `-p` null-record row cell formatting (empty GENO, "0;0" bias, empty flank seqs)

**Java behavior:** same trigger branch as UB-1's first half —
`SimplePostProcessModule.java:75-77` wraps a null reference in
`new SimpleOutputVariant(null, …)` and prints it. The null path renders printer
defaults rather than Variant fields: genotype cell `""`
(`SimpleOutputVariant.java:22`), strand-bias cell `"0;0"` (:24), flanking-sequence
cells `""` (`OutputVariant.java:22-23`). The resulting `.var` row carries an empty
genotype token; downstream R/Perl consumers must special-case it.

**Our replication stance:** replicated byte-for-byte (P8a; `variations/mod.rs`
`null_record` flag + `output/table.rs` renderer). Cell-exact unit tests pin cols
14/16/31/32 to the oracle shape (`simple_36_null_record_renders_printer_defaults`).

**Draft issue text:**

> **Title:** `-p` mode emits a `.var` row with an empty GENO column when a position has no reference variant
>
> In v1.8.3, `SimplePostProcessModule.java:75-77` prints a
> `SimpleOutputVariant(null, …)` row whose genotype cell is the empty string,
> strand-bias is the literal `0;0`, and flank-sequence cells are empty
> (`SimpleOutputVariant.java:22,:24`, `OutputVariant.java:22-23`). Empty tokens in
> a tab-delimited table are silently mis-parsed by column-aware consumers; fix
> direction: skip the row or emit explicit placeholders. Output-changing, release
> note required.

---

## UB-3 — `.var` row order across positions follows HashMap iteration order (#282 family)

**Java behavior:** `SimplePostProcessModule.accept` iterates
`AlignedVarsData.alignedVariants`, a plain `HashMap<Integer, Vars>`
(`ToVarsBuilder.java:93`, `AlignedVarsData.java:14`), so `.var` row order across
genomic positions follows hash-table slot order rather than coordinate order
(`SimplePostProcessModule.java:35`). Deterministic for identical inputs/JVM, but
changes when upstream insertion history changes; downstream, `var2vcf_valid.pl`
keeps only the highest-AF row per raw-start bucket and breaks exact AF ties
first-in-stream — so WHICH allele survives a tie depends on hash layout.

**Our replication stance:** **NOT replicated** — Rust emits ascending position
order (BTreeMap), per CLAUDE.md §8 (hash-order iteration in output paths is
forbidden, #282 family). Record SETS are equal everywhere; the only observable
effect is which of two exact-AF-tie alleles sharing a raw-start bucket survives,
at most one VCF record per such locus. Known instance: NA24385 chr20:20202875
(`TA>T` java-first vs `T>A` rust-first, both af=0.0179; Tier-2 slice s1a).
Closed as permanent known limitation (P3 timebox; matching would require
emulating Java HashMap lifecycle wholesale).

**Draft issue text:**

> **Title:** .var row order across positions follows HashMap iteration order (#282 family)
>
> In v1.8.3, `SimplePostProcessModule.accept` iterates `AlignedVarsData.alignedVariants`,
> a `HashMap<Integer, Vars>` (ToVarsBuilder.java:93), so `.var` row order across
> genomic positions follows hash-table slot order rather than coordinate order.
> The order is stable for identical inputs/JVM but changes when upstream insertion
> history changes (region set, read tail), which manifests downstream:
> var2vcf_valid.pl keeps only the highest-AF row per raw-start bucket and breaks
> AF ties first-in-stream, so WHICH allele survives an exact AF tie depends on
> hash layout. Sorting `alignedVariants` by key (e.g. TreeMap) before iteration
> would make tie survival deterministic; note this is output-order-changing and
> should carry a release-note. Repro: NA24385 chr20:20156208-20234556 `-f 0.001`,
> rows at printed start 20202875 (map keys 20202875/20202876) emit del-before-snv.

---

## UB-4 — Integration-fixture inconsistency: overlapping CSV reference chunks

**Behavior (upstream test data):** the fasta-source CSV behind
`testdata/fastas/hard_clip_case.fa` declares two overlapping chunks for contig
`test` — `test,5234..7764` and `test,5475..8024` — which conflict at 1734 of the
2290 shared positions. The checked-in expected files for
`hard_clip_next_to_del_test2.bam` match a FIRST-chunk-wins reconstruction of that
overlap (proven: a live `vardict-java` run under a first-wins fasta reproduces the
fixture exactly), while the repository's converter script resolves overlaps
last-write-wins — and the sibling `test1` fixture matches last-wins. The suite
therefore contains mutually inconsistent (reference ↔ expected-file) pairs; no
single fasta satisfies both.

**Our stance:** nothing to replicate in an engine — this is an upstream test-data /
generator defect. We keep the corpus case as a documented permanent known non-bug
(DIVERGENCES.md "Known non-bugs") and validate the underlying behavior against the
live oracle instead (P7 dossier + P8a replication evidence). Provenance note: the
"repository's converter script" above was our unbundled copy
`harness/convert_fastas.py` (removed in audit task A3), whose reconstruction of the
CSV chunks resolved overlaps last-write-wins.

**Draft issue text:**

> **Title:** Integration fixtures encode mutually inconsistent references (overlapping CSV chunks, first-vs-last-write wins)
>
> The fasta-generator CSV for the integration suite contains overlapping chunks
> for the same contig (`test,5234..7764` vs `test,5475..8024` in
> `hard_clip_case.fa`) that disagree at 1734/2290 shared positions. Expected
> outputs were captured against different conflict resolutions: the
> `hard_clip_next_to_del_test2` fixture matches a first-chunk-wins reconstruction
> (reproducible exactly by running v1.8.3 under such a fasta), while the bundled
> converter applies last-write-wins and the `test1` fixture matches that. As
> shipped, `vardict-java` cannot reproduce one of its own expected files. Fix
> direction: reject or deterministically resolve overlapping chunks in the
> converter, then regenerate the affected fastas and expected files together.
