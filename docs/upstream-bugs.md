# Upstream bugs — VarDictJava v1.8.3

Consolidated, actionable list of VarDictJava v1.8.3 behaviors that are part of the
specification under rewrites.bio §5 (replicate output bugs byte-for-byte, track
upstream, never silently fix), plus one upstream-repo defect that affects validation
infrastructure. Supersedes the drafts formerly embedded in `DIVERGENCES.md`.

**Version pin:** VarDictJava **v1.8.3 @ `0886261`**, Perl/R submodule @ `009e017`.
All line references below are to that pin.

| # | Item | Our stance |
|---|------|------------|
| [UB-1](#ub-1--p-spurious-all-zero-record-at-position-00) | `-p` spurious all-zero record at position 0/0 | **Replicated** (commit `4d113b2`) |
| [UB-2](#ub-2--p-null-record-row-cell-formatting-empty-geno-00-bias-empty-flank-seqs) | `-p` null-record row cell formatting | **Replicated** (`variations/mod.rs`) |
| [UB-3](#ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family) | Cross-position `.var` emission order follows HashMap layout | **Replicated** (via JavaHashMap Determinism Engine) |
| [UB-4](#ub-4-integration-fixture-inconsistency-overlapping-csv-reference-chunks) | Integration fixtures encode mutually inconsistent references | N/A (test-data defect; validated vs live oracle) |
| [UB-5](#ub-5-simplemode-structural-variant-split-read-requirement-var2vcf_validpl197-205) | SimpleMode SV split-read requirement (`splitreads < 1` dropped) | **Replicated** (`lib.rs:403-414, 470, 499`) |
| [UB-6](#ub-6-mate-pair-boundary-parameter-in-mark_dup_sv-structuralvariantsprocessorjava1468) | Mate-pair boundary parameter in `mark_dup_sv` (`pe` decrement) | **Replicated** (`sv.rs:3039, 3105-3111`) |
| [UB-7](#ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205) | Multi-allelic site single-winner collapsing without `-A` | **Replicated** (`lib.rs:467-505`) |

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

**Our replication stance:** replicated byte-for-byte (commit `4d113b2`;
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

**Our replication stance:** replicated byte-for-byte (`variations/mod.rs`
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

**Our replication stance:** **Replicated** via JavaHashMap Determinism Engine (commit `d41a6f1`) via
[`vardictrust_core::javamap::JavaHashMap`](../crates/vardictrust-core/src/javamap.rs).
It faithfully emulates OpenJDK 8+ `java.util.HashMap` hash-spreading (`h ^ (h >> 16)`),
power-of-two table indexing, load factor 0.75, resize threshold doubling, and bucket
iteration order across variant maps (`non_insertion_variants`, `insertion_variants`,
and `aligned_variants`). This completely eliminated the single residual AF-tie at
`chr20:20202875` on dual-engine test slice `s1a`, emitting `TA>T` before `T>A` and achieving
100% byte-for-byte parity (1048/1048 variants, 0 diffs). Per AGENTS.md §8, `JavaHashMap`
is the only sanctioned non-sorted map in output-reaching paths.

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
live oracle instead. In `harness/integration_xfail.tsv`, this case is tracked as `XFAIL`
pinned to live-oracle SHA-256 `96abe5c4bb6e3ff4572e361fddae1cd61354f102551f494488456a96d17829b2`.
Both live VarDictJava 1.8.3 and VarDictRust emit identical 7 rows under the canonical FASTA
(4 complex indels, 1 reference pileup call, and 2 upstream quirk rows UB-1 and UB-2).
`harness/run_integration.sh` verifies this exact SHA-256 signature, treating any deviation
as `FAIL/REGRESSION` and unexpected fixture pass as `XPASS`. Provenance note: the
"repository's converter script" above was our unbundled copy
`harness/convert_fastas.py` (removed during codebase audit), whose reconstruction of the
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

---

## UB-5 — SimpleMode structural variant split-read requirement (`var2vcf_valid.pl:197-205`)

**Java / Perl behavior:** In `var2vcf_valid.pl:197-205`, when `-U` (`--nosv`) is not passed, SV calling is active. The Perl converter parses the SV field `$tamp` to extract split reads:
```perl
unless( $isamp ) {
    ($splitreads, $spanpairs, $cluster) = (split(/-/, $tamp)) if ( defined($tamp) && $tamp =~ /-/ ); 
    if ( $alt =~ /</ ) {
        # Ignore SV's without split read support for now until a better criteria
        unless ( $splitreads >= $opt_T ) {
            $pinfo1 = "";
            next;
        }
```
Any structural variant allele containing `<` (e.g., `<DUP>`, `<INV>`) that has only discordant mate-pair cluster evidence but zero split-read evidence (`splitreads < $opt_T`, where default `$opt_T = 1`) is unconditionally discarded during VCF formatting. Notably, these variants **are emitted** in the intermediate `.var` table, creating a discrepancy between the table and VCF output.

**Our replication stance:** **Replicated** in `crates/vardictrust-core/src/lib.rs:403-414, 470, 499` (`is_sv && splitreads < 1` suppresses VCF line generation for SimpleMode), ensuring exact bitwise parity with the full VarDictJava pipeline output.

**Draft issue text:**

> **Title:** SimpleMode VCF conversion silently drops discordant-pair SVs lacking split-read support
>
> In `var2vcf_valid.pl:197-205`, SV alleles containing `<` (such as `<DUP>` or `<INV>`) are dropped when `splitreads < $opt_T` (default 1). These candidate structural variants are detected by `StructuralVariantsProcessor` and printed to the `.var` table with `SV` tags (e.g., `0-12-1`), but never reach the final VCF output. If discordant pairs provide sufficient statistical evidence, dropping them without split reads should be controlled by a dedicated CLI flag rather than an unconditional filter.

---

## UB-6 — Mate-pair boundary parameter in `mark_dup_sv` (`StructuralVariantsProcessor.java:1468`)

**Java behavior:** In `StructuralVariantsProcessor.java:1413-1468`, the cluster end coordinate `pe` is decremented prior to marking duplication SVs:
```java
mlen = pe - tbp;
bp = tbp;
pe--;
end = pe;
...
markDUPSV(bp, pe, Collections.singletonList(svStructures.svrdup), maxReadLength);
```
Inside `markDUPSV`, the overlap predicate `isOverlap(start1, end1, start2, end2, rlen)` checks `if (start1 >= end2 || start2 >= end1) return false;`. Because Java passes the post-decrement `pe` as `end1`, any reverse-strand mate cluster whose reads begin exactly at the original un-decremented `pe` fails the overlap condition.

**Our replication stance:** **Replicated** in `crates/vardictrust-core/src/sv.rs:3039, 3105-3111`, passing `pe` (post-decrement) rather than `pe + 1` to match Java's exact discordant pair clustering and breakpoint boundaries.

**Draft issue text:**

> **Title:** Off-by-one window truncation in `markDUPSV` due to pre-decremented `pe` parameter
>
> In `StructuralVariantsProcessor.java:1416`, `pe` is decremented (`pe--`) before being passed to `markDUPSV` at line 1468. This shifts the search window right-boundary by 1 base to the left, causing discordant mate pairs that start at the exact cluster boundary to be excluded from duplication clustering. Passing `pe + 1` or the pre-decrement boundary would correctly include boundary-aligned discordant mates.

---

## UB-7 — Multi-allelic site single-winner collapsing (`var2vcf_valid.pl:135-205`)

**Java / Perl behavior:** In `var2vcf_valid.pl:135-139`, unless `-A` is explicitly passed, the multi-allelic candidate count is clamped to `$ALL = 1`:
```perl
foreach my $p (@pos) {
    my @tmp = sort { $b->[14] <=> $a->[14] } @{ $hash{ $chr }->{ $p } };
    my $ALL = $opt_A ? @tmp + 0 : 1;
    for(my $i = 0; $i < $ALL; $i++) {
        my ($sample, ..., $ref, $alt, ..., $af) = @{ $tmp[$i] };
        next unless ( $ref );
```
Only the single candidate allele with the highest allele frequency (`$tmp[0]`) is considered at each genomic locus. If that winner candidate is skipped (for example, if `$ref` is empty from a `-p` null-record, or if an SV candidate has `splitreads < 1` per UB-5), the loop executes `next` and terminates immediately. As a result, legitimate alternative variants at the same locus with lower allele frequencies are completely excluded from the VCF.

**Our replication stance:** **Replicated** in `crates/vardictrust-core/src/lib.rs:467-505`, maintaining a `pending: BTreeMap<i32, SimpleVcfCandidate>` that retains only the top-AF candidate per locus in single-sample VCF mode.

**Draft issue text:**

> **Title:** Single-winner multi-allelic collapsing drops secondary variants when primary candidate is filtered
>
> In `var2vcf_valid.pl:135-205`, when `-A` is omitted, only the highest-AF candidate is inspected. If that candidate triggers any `next` condition (e.g. empty REF or zero split reads), the entire position is abandoned, dropping any secondary variants that would have passed filters. The loop should evaluate candidate records until an eligible variant is found or emit all qualifying alleles in standard multi-allelic VCF notation.

