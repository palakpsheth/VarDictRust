# Upstream bugs — VarDictJava v1.8.3

Consolidated, actionable list of VarDictJava v1.8.3 behaviors that are part of the
specification under [rewrites.bio](https://rewrites.bio) §5 (replicate output bugs byte-for-byte, track
upstream, never silently fix), plus one upstream-repo defect that affects validation
infrastructure. Supersedes the drafts formerly embedded in [`DIVERGENCES.md`](../DIVERGENCES.md).

**Version pin:** VarDictJava **v1.8.3 @ [`0886261`](https://github.com/AstraZeneca-NGS/VarDictJava/commit/0886261)**, Perl/R submodule @ [`009e017`](https://github.com/AstraZeneca-NGS/VarDict/commit/009e017).
All line references below are to that pin.

| # | Class | Item | Our stance |
|---|-------|------|------------|
| [UB-1](#ub-1) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | `-p` spurious all-zero record at position 0/0 | **Replicated** ([`postprocess/simple.rs`](../crates/vardictrust-core/src/postprocess/simple.rs)) |
| [UB-2](#ub-2) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | `-p` null-record row cell formatting | **Replicated** ([`variations/mod.rs`](../crates/vardictrust-core/src/variations/mod.rs)) |
| [UB-3](#ub-3) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | Cross-position `.var` emission order follows HashMap layout ([#282](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282) family) | **Replicated** (via [JavaHashMap Engine](../crates/vardictrust-core/src/javamap.rs)) |
| [UB-4](#ub-4) | [Class D](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | Integration fixtures encode mutually inconsistent references | N/A (test-data defect; validated vs live oracle) |
| [UB-5](#ub-5) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | SimpleMode SV split-read requirement (`splitreads < 1` dropped) | **Replicated** ([`lib.rs:403-414`](../crates/vardictrust-core/src/lib.rs)) |
| [UB-6](#ub-6) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | Mate-pair boundary parameter in `mark_dup_sv` (`pe` decrement) | **Replicated** ([`sv.rs:3039`](../crates/vardictrust-core/src/sv.rs)) |
| [UB-7](#ub-7) | [Class A](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | Multi-allelic site single-winner collapsing without `-A` | **Replicated** ([`lib.rs:467-505`](../crates/vardictrust-core/src/lib.rs)) |
| [UB-8](#ub-8) | [Class D](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11) | Oracle JVM SIGKILL flake under shared-host memory pressure (exit 137, no output change) | **Tracked** (environmental; diagnostics in [`harness/run_oracle.sh`](../harness/run_oracle.sh)) |

**Third arm:** pinned **vardictcpp 3 `hc52dbad_0` bioconda** runs beside every gate
arm (advisory only; gate stays java==rust; `CPP_SKIP=1` escapes). Its known
divergences below are **normalized, not engine bugs** — see
[Known vardictcpp third-arm divergences](#cpp-divergences).

---

<a id="ub-1"></a>
<a id="ub-1--p-spurious-all-zero-record-at-position-00"></a>
## UB-1 — `-p`: spurious all-zero record at position 0/0

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

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

**Our replication stance:** replicated byte-for-byte in [`crates/vardictrust-core/src/postprocess/simple.rs`](../crates/vardictrust-core/src/postprocess/simple.rs) fallback branch keeping `Variant::new()` defaults and [`crates/vardictrust-core/src/output/table.rs`](../crates/vardictrust-core/src/output/table.rs). Verified sorted-diff-empty vs live oracle on the hard_clip case under [Tier-1](testing/test-matrix-and-ci-audit.md#tier-1) table differential testing via [`harness/compare_var.py`](../harness/compare_var.py).

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

<a id="ub-2"></a>
<a id="ub-2--p-null-record-row-cell-formatting-empty-geno-00-bias-empty-flank-seqs"></a>
## UB-2 — `-p` null-record row cell formatting (empty GENO, "0;0" bias, empty flank seqs)

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

**Java behavior:** same trigger branch as [UB-1](#ub-1)'s first half —
`SimplePostProcessModule.java:75-77` wraps a null reference in
`new SimpleOutputVariant(null, …)` and prints it. The null path renders printer
defaults rather than Variant fields: genotype cell `""`
(`SimpleOutputVariant.java:22`), strand-bias cell `"0;0"` (:24), flanking-sequence
cells `""` (`OutputVariant.java:22-23`). The resulting `.var` row carries an empty
genotype token; downstream R/Perl consumers must special-case it.

**Our replication stance:** replicated byte-for-byte in [`crates/vardictrust-core/src/variations/mod.rs`](../crates/vardictrust-core/src/variations/mod.rs) (`null_record` flag) and [`crates/vardictrust-core/src/output/table.rs`](../crates/vardictrust-core/src/output/table.rs) renderer. Cell-exact unit tests pin cols 14/16/31/32 to the oracle shape (`simple_36_null_record_renders_printer_defaults`).

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

<a id="ub-3"></a>
<a id="ub-3-var-row-order-across-positions-follows-hashmap-iteration-order-282-family"></a>
## UB-3 — `.var` row order across positions follows HashMap iteration order ([#282](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282) family)

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

**Java behavior:** `SimplePostProcessModule.accept` iterates
`AlignedVarsData.alignedVariants`, a plain `HashMap<Integer, Vars>`
(`ToVarsBuilder.java:93`, `AlignedVarsData.java:14`), so `.var` row order across
genomic positions follows hash-table slot order rather than coordinate order
(`SimplePostProcessModule.java:35`). Deterministic for identical inputs/JVM, but
changes when upstream insertion history changes; downstream, `var2vcf_valid.pl`
keeps only the highest-AF row per raw-start bucket and breaks exact AF ties
first-in-stream — so WHICH allele survives a tie depends on hash layout.

**Our replication stance:** **Replicated** via JavaHashMap Determinism Engine in
[`vardictrust_core::javamap::JavaHashMap`](../crates/vardictrust-core/src/javamap.rs).
It faithfully emulates OpenJDK 8+ `java.util.HashMap` hash-spreading (`h ^ (h >> 16)`),
power-of-two table indexing, load factor 0.75, resize threshold doubling, and bucket
iteration order across variant maps (`non_insertion_variants`, `insertion_variants`,
and `aligned_variants`). This completely eliminated the single residual AF-tie at
`chr20:20202875` on [Tier-2](testing/test-matrix-and-ci-audit.md#tier-2) dual-engine test slice `s1a`, emitting `TA>T` before `T>A` and achieving
100% byte-for-byte parity (1048/1048 variants, 0 diffs). Per the VarDictRust parity architecture (see [`CONTRIBUTING.md`](../CONTRIBUTING.md)), `JavaHashMap`
is the only sanctioned non-sorted map in output-reaching paths.

**Draft issue text:**

> **Title:** .var row order across positions follows HashMap iteration order ([#282](https://github.com/AstraZeneca-NGS/VarDictJava/issues/282) family)
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

<a id="ub-4"></a>
<a id="ub-4-hard_clip_casefa-fixture-drift-between-csv-chunks-and-reconstructed-fasta"></a>
<a id="ub-4-hard_clip_next_to_del_test2-chunk-overlap-disagreement"></a>
<a id="ub-4-integration-fixture-inconsistency-overlapping-csv-reference-chunks"></a>
## UB-4 — Integration-fixture inconsistency: overlapping CSV reference chunks

**Parity Class:** [Class D (Known Upstream Defects & Fixture Flaws)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

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
([`DIVERGENCES.md` ("Class D Known Non-Bugs")](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)) and validate the underlying behavior against the
live oracle instead. In [`harness/integration_xfail.tsv`](../harness/integration_xfail.tsv) across the [Tier-1](testing/test-matrix-and-ci-audit.md#tier-1) integration test suite, this case is tracked as `XFAIL`
pinned to live-oracle SHA-256 `96abe5c4bb6e3ff4572e361fddae1cd61354f102551f494488456a96d17829b2`.
Both live VarDictJava 1.8.3 and VarDictRust emit identical 7 rows under the canonical FASTA
(4 complex indels, 1 reference pileup call, and 2 upstream quirk rows [UB-1](#ub-1) and [UB-2](#ub-2)).
[`harness/run_integration.sh`](../harness/run_integration.sh) verifies this exact SHA-256 signature, treating any deviation
as `FAIL/REGRESSION` and unexpected fixture pass as `XPASS`. Provenance note: the test fixture reference materializer [`harness/convert_fastas.py`](../harness/convert_fastas.py) reconstructs FASTA files from upstream CSV chunks resolving overlaps last-write-wins.

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

<a id="ub-5"></a>
<a id="ub-5-simplemode-structural-variant-split-read-requirement-var2vcf_validpl197-205"></a>
## UB-5 — SimpleMode structural variant split-read requirement (`var2vcf_valid.pl:197-205`)

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

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

**Our replication stance:** **Replicated** in [`crates/vardictrust-core/src/lib.rs:403-414, 470, 499`](../crates/vardictrust-core/src/lib.rs) (`is_sv && splitreads < 1` suppresses VCF line generation for SimpleMode), ensuring exact bitwise parity with the full VarDictJava pipeline output.

**Draft issue text:**

> **Title:** SimpleMode VCF conversion silently drops discordant-pair SVs lacking split-read support
>
> In `var2vcf_valid.pl:197-205`, SV alleles containing `<` (such as `<DUP>` or `<INV>`) are dropped when `splitreads < $opt_T` (default 1). These candidate structural variants are detected by `StructuralVariantsProcessor` and printed to the `.var` table with `SV` tags (e.g., `0-12-1`), but never reach the final VCF output. If discordant pairs provide sufficient statistical evidence, dropping them without split reads should be controlled by a dedicated CLI flag rather than an unconditional filter.

---

<a id="ub-6"></a>
<a id="ub-6-mate-pair-boundary-parameter-in-mark_dup_sv-structuralvariantsprocessorjava1468"></a>
## UB-6 — Mate-pair boundary parameter in `mark_dup_sv` (`StructuralVariantsProcessor.java:1468`)

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

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

**Our replication stance:** **Replicated** in [`crates/vardictrust-core/src/sv.rs:3039, 3105-3111`](../crates/vardictrust-core/src/sv.rs), passing `pe` (post-decrement) rather than `pe + 1` to match Java's exact discordant pair clustering and breakpoint boundaries.

**Draft issue text:**

> **Title:** Off-by-one window truncation in `markDUPSV` due to pre-decremented `pe` parameter
>
> In `StructuralVariantsProcessor.java:1416`, `pe` is decremented (`pe--`) before being passed to `markDUPSV` at line 1468. This shifts the search window right-boundary by 1 base to the left, causing discordant mate pairs that start at the exact cluster boundary to be excluded from duplication clustering. Passing `pe + 1` or the pre-decrement boundary would correctly include boundary-aligned discordant mates.

---

<a id="ub-7"></a>
<a id="ub-7-multi-allelic-site-single-winner-collapsing-var2vcf_validpl135-205"></a>
## UB-7 — Multi-allelic site single-winner collapsing (`var2vcf_valid.pl:135-205`)

**Parity Class:** [Class A (Faithfully Emulated Upstream Quirks)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

**Java / Perl behavior:** In `var2vcf_valid.pl:135-139`, unless `-A` is explicitly passed, the multi-allelic candidate count is clamped to `$ALL = 1`:
```perl
foreach my $p (@pos) {
    my @tmp = sort { $b->[14] <=> $a->[14] } @{ $hash{ $chr }->{ $p } };
    my $ALL = $opt_A ? @tmp + 0 : 1;
    for(my $i = 0; $i < $ALL; $i++) {
        my ($sample, ..., $ref, $alt, ..., $af) = @{ $tmp[$i] };
        next unless ( $ref );
```
Only the single candidate allele with the highest allele frequency (`$tmp[0]`) is considered at each genomic locus. If that winner candidate is skipped (for example, if `$ref` is empty from a `-p` null-record, or if an SV candidate has `splitreads < 1` per [UB-5](#ub-5)), the loop executes `next` and terminates immediately. As a result, legitimate alternative variants at the same locus with lower allele frequencies are completely excluded from the VCF.

**Our replication stance:** **Replicated** in [`crates/vardictrust-core/src/lib.rs:467-505`](../crates/vardictrust-core/src/lib.rs), maintaining a `pending: BTreeMap<i32, SimpleVcfCandidate>` that retains only the top-AF candidate per locus in single-sample VCF mode.

**Draft issue text:**

> **Title:** Single-winner multi-allelic collapsing drops secondary variants when primary candidate is filtered
>
> In `var2vcf_valid.pl:135-205`, when `-A` is omitted, only the highest-AF candidate is inspected. If that candidate triggers any `next` condition (e.g. empty REF or zero split reads), the entire position is abandoned, dropping any secondary variants that would have passed filters. The loop should evaluate candidate records until an eligible variant is found or emit all qualifying alleles in standard multi-allelic VCF notation.

---

<a id="ub-8"></a>
<a id="ub-8-oracle-jvm-sigkill-flake-under-shared-host-memory-pressure"></a>
## UB-8 — Oracle JVM SIGKILL flake under shared-host memory pressure (exit 137)

**Parity Class:** [Class D (Known Upstream Defects & Fixture Flaws)](../DIVERGENCES.md#the-four-canonical-parity-classes-rewritesbio-5-8-11)

**Behavior (environmental, not a VarDictJava logic bug):** on the full-chr20 hg002
oracle arm (`-f 0.01 -th 8`, ~90 s JVM runtime, MaxRSS ≈ 4.1 GB under `-Xmx8g`),
the `vardict-java` process is intermittently SIGKILLed ~85–90 s in — near
completion — while the downstream R/Perl stages exit cleanly. `compare.sh`
reports `ENGINE-FAIL: oracle`. An identical-input retry passes byte-identically
(same 227562-variant MATCH), so this is a flake, not input-dependent behavior.

**Evidence (2026-10-13, 5 consecutive `compare.sh` hg002 runs):** 1 MATCH, then 4
consecutive oracle `ENGINE-FAIL`s, every failure with the same signature —
`java.time` `Exit status: 137`, wrapper stderr:

```
_: line 18: 1178932 Killed                  vardict-java -G "$1" -f "$2" ...
     1178933 Done                    | teststrandbias.R
     1178934 Done                    | var2vcf_valid.pl -N "$3" -E -f "$2" "${V2V_ARGS[@]}" > "$6"
oracle-stage-fail: vardict-java exit 137
```

Prime suspect: `systemd-oomd` — controller ground truth confirms it is active on
this host, with journal policy kills on `user.slice` pressure >50% for >20s
(Oct 11 18:19, 22:32, whole scopes). oomd kills are journal-only by design, so a
silent kernel ring through this run's kill window matches that signature exactly
rather than ruling oomd out. The in-run sender stays unattributed (no journal
entries that day) — prime suspect, not a confirmed attribution.

Ruled out: our wrapper (pipe/trap/timeout/env verified; exit code faithfully
propagated), JVM heap exhaustion (4.1 GB < 8 GB cap, no Java exception), and
system-wide exhaustion (`free -h` ~23 GB available at kill time). Context is
shared-host userspace pressure (parallel agent sessions; this box has a
documented history of host-pressure kills and manual `kill -9` cleanups).

**Our stance:** nothing to replicate (no output effect — killed runs produce no
VCF) and nothing to fix in our wrapper; retry loops are forbidden as a fix.
[`harness/run_oracle.sh`](../harness/run_oracle.sh) now persists
`$KEEP_DIR/oracle-fail-<stage>.log` (stage stderr + exit code + `free -h`) on any
engine-stage failure so future flakes arrive with their signature attached.
Guidance: re-run `compare.sh`; keep ambient load down
([`harness/courtesy.sh`](../harness/README.md), `-th 8`); do not chase this by
retuning oracle heap/threads.

---

<a id="cpp-divergences"></a>
## Known vardictcpp third-arm divergences (advisory, non-gating)

Applies to pinned **vardictcpp 3 `hc52dbad_0` bioconda** (`benchmarks/env/vardict-oracle.yml`).
Each item was observed live on real fixtures; none touches the java==rust gate.
Cpp numbers live in evidence `extra` (`cpp_version`, `cpp_head`, `cpp_f1`,
`cpp_wall_s`, `cpp_maxrss_kb`, `cpp_skip_reason`) and are never duplicated here.

| # | Divergence | Surface |
|---|-----------|---------|
| [CPP-1](#cpp-1) | VCF header provenance + `END` INFO key deltas | `harness/compare.sh` third arm (`compare_cpp.txt`) |
| [CPP-2](#cpp-2) | `-J`/`-j` parsed but ignored | Extended-coverage M4 cell |
| [CPP-3](#cpp-3) | `-a` amplicon mode is TSV-only (`-vcf` ignored) | Tier-2 amplicon arm (honest SKIP) |
| [CPP-4](#cpp-4) | Concatenated single-dash flags rejected | Tier-1 concatenated-flags fixtures |
| [CPP-5](#cpp-5) | Boolean `-z`, positional BED, unpinned `--chunk` | `harness/run_cpp.sh` contract |
| [CPP-6](#cpp-6) | Zero-`##contig` VCF emitter (tabix needs completion) | Truth-concordance scaffolding |
| [CPP-7](#cpp-7) | CHROM resolved from reference, not BAM, on duplicate-name FASTAs | b37 extended-coverage cells |
| [CPP-8](#cpp-8) | N-REF null rows emitted in native `-vcf` | b37 extended-coverage cells |
| [CPP-9](#cpp-9) | Caller deltas: MNP joining + large SV alleles | hg002/exome/hg008 Tier-2 arms |
| [CPP-10](#cpp-10) | `-DP`/`--DP` parsed but ignored (printer stays stdout) | Extended-coverage M17 cell |

<a id="cpp-1"></a>
### CPP-1 — VCF header provenance + `END` INFO key

cpp emits its own header lineage (no `##reference`, its own `##INFO` set incl.
`END` on every germline SNV/indel) where VarDictJava emits `##reference`/`##contig`
and restricts `INFO END` to structural variants. Variant *calls* and numeric/string
values agree 100% on `micro` (`2/2`), 4 of 10 dev slices (`s3`, `s4`, `s8`, `s9`),
2 of 17 extended-coverage modes (`ext-m2` `<INV>`, `ext-m3` `<DUP>`), and 12 of 15
non-RNA multi-modal slices (`HG005_1A/1B/1C`, `SEQC2_2A/2B/2C`, `HG002_3C`,
`HG008_4A/4C`, `NA24385_5A/5B/5C`); [`compare_vcf.py`](../harness/compare_vcf.py)
still reports `CPP: MISMATCH` on header/`END`-key deltas. Known-normalized:
never an engine bug, never gating.

<a id="cpp-2"></a>
### CPP-2 — `-J`/`-j` parsed but ignored

cpp accepts CRISPR cut-site flags (`-J 50454941 -j 25`) but its own banner warns
"not all affect output yet": on the M4 cell it emits all regional rows where
Java target-windows to 3, so an advisory `CPP-MISMATCH` on row count is expected.

<a id="cpp-3"></a>
### CPP-3 — `-a` amplicon mode is TSV-only (`-vcf` ignored)

Amplicon mode without a BED errors (exit 1 → loud `CPP: SKIP`), and with a
BED cpp emits its 38/40-col TSV table even when `-vcf` is passed (verified:
`-vcf` before and after the BED both yield TSV). When `-a` and `-vcf` are passed,
`vardictcpp` runs Fisher exact strand bias and outputs a 40-column TSV table
matching Java's `teststrandbias.R` output (`compare_var.py` MATCH on `amplicon`
and `s6`). A VCF-vs-TSV comparison would MISMATCH meaninglessly, so `compare.sh`
(the Tier-2 VCF comparator) honestly records
`CPP: SKIP: cpp emitted non-VCF surface` naming the cause (locked by CPP11).
TSV-level amplicon parity belongs to a Tier-1 `compare_var.py` cell, not the
VCF gate.

<a id="cpp-4"></a>
### CPP-4 — Concatenated single-dash flags rejected

cpp rejects concatenated forms such as `-f0.0025-F0x700-Q10` (`unknown option`,
exit 120) that VarDictJava accepts. Affected fixtures record loud
`CPP: SKIP (cpp engine exit 120)`; no flag-splitting is added (would break the
identical-flags contract).

<a id="cpp-5"></a>
### CPP-5 — Boolean `-z`, positional BED, unpinned `--chunk` (retired `-R` noise)

`run_cpp.sh` mirrors the oracle exactly: BED travels positionally and bare
`-z` stays bare — cpp's `-z` takes no value, and a mirrored `-z 1` leaks a
stray positional that cpp opens as a BED file (`cannot open BED file '1'`;
`-R` runs only survived because the stray was ignored with a warning).
`--chunk` is deliberately unpinned (cpp default rewrites the `Seg` column);
the peak jobs-per-node selector stays out of the gate by design.

<a id="cpp-6"></a>
### CPP-6 — Zero-`##contig` VCF emitter

cpp's native `-vcf` emitter carries zero `##contig` lines, so `tabix -p vcf`
refuses to index its calls file. Truth concordance completes the missing lines
from the truth header (truth first, awk-ID dedup; `##fileformat` stays line 1)
as evaluation scaffolding only — no record, FILTER, QUAL, or INFO byte is
touched, so the reported cpp F1 is unaffected. Handled forward-compatibly if a
future pin emits contigs.

<a id="cpp-7"></a>
### CPP-7 — CHROM resolved from reference, not BAM, on ambiguous FASTAs

Root cause verified in audit: `VarDictJava/testdata/fastas/GRCh37.fa` (materialized
from sparse `GRCh37.fa.csv`) contains *both* `>20` (length 62,376,064; real biological sequence)
and `>chr20` (length 60,470,068; sparse stub that is all `N`s across `20156208-20234556`).
When a BED specifies `chr20` against a BAM with `@SQ SN:20`, Java and Rust resolve
the contig through the BAM dictionary to `20`, whereas `vardictcpp` queries `>chr20`
directly in the FASTA and emits `REF=N` calls across the window (CPP-8) or exits 1
when the `chr`-prefixed contig is absent (`chr11` in `GRCh37.fa`). When the BED uses
the BAM's `20` contig, `vardictcpp` emits 0 `REF=N` rows and matches 4/4 (100%) on
`tiny` and 1,042/1,048 (99.4%) on `s1a`. In truth concordance, evidence records both
`cpp_f1` (unnormalized) and `cpp_f1_normalized`. Recorded advisory, never gating.

<a id="cpp-8"></a>
### CPP-8 — N-REF null rows emitted in native `-vcf`

Direct consequence of CPP-7: when `vardictcpp` queries the decoy `>chr20` N-stub in
`GRCh37.fa`, it reads all `N`s and emits a variant call with `REF=N` at every covered
position across the window (up to 78,304 false-positive rows on b37 slices), whereas
Java's `var2vcf` drops null rows. When contig names match the BAM (`20`), `vardictcpp`
emits zero `REF=N` rows. Truth concordance evidence records both `cpp_f1` and
`cpp_f1_normalized` (filtering out `REF=^N+$` calls and aligning contig naming).
Advisory only, never gating.

<a id="cpp-9"></a>
### CPP-9 — Caller deltas: MNP joining, tie-breaking, large SV deletions, and somatic post-processing (measured, advisory)

Beyond serialization ([CPP-1](#cpp-1)), record-level `CHROM/POS/REF/ALT` and
field-level comparison shows genuine caller-behavior deltas in cpp-vs-java
while `rust == java` `MATCH` holds on the same inputs:

- **Dropped multi-hundred-bp deletions (2 per deep WGS/WES dataset):**
  - `hg002-tier2`: drops `chr20:17402396` (287 bp Alu deletion, `SPLITREAD=2;SPANPAIR=2`) and `chr20:18661577` (747 bp deletion, `SPLITREAD=2;SPANPAIR=2`).
  - `hg008-tier2` (and `HG008_4B`): drops `chr1:6173467` (228 bp deletion, emitting an un-realigned `chr1:6173467 A→C` SNV instead) and `chr3:127563710` (301 bp deletion).
  - `seqc2-exome-tier2`: drops `chr2:178443383 T→<DEL>` (symbolic deletion SV) and `chr22:36511000` (495 bp deletion).
- **MNP decomposition & equal-AF single-winner tie-breaking (`hg002` germline, 346/227,562 loci = 0.15%):**
  All 346 Cpp-only `(CHROM, POS, REF, ALT)` keys share the exact same `(CHROM, POS)` coordinate with 346 Java-only keys:
  - **15 loci** are joint MNPs in Java decomposed into single-base SNVs at the first base in cpp (e.g. `chr20:10016944 TGA→CGC` in Java/Rust vs `T→C` in cpp, also seen in `HG002_3A`; `chr20:10268361 CT→AA` vs `C→A`).
  - **331 loci** are equal-AF (`VD=1` or `2`) single-winner tie-breaks where Java's `java.util.HashMap` bucket order ([UB-3](#ub-3) / [UB-7](#ub-7)) selects one low-frequency allele (e.g. `chr20:10093752 C→T`, `AF=0.013`; `chr20:20008963 G→T` and `20012549 G→T` in `HG002_3B`) while cpp's deterministic `A < C < G < T` map order selects another (`C→A`, `G→C`, `G→A`) at identical `AF`.
- **Unrounded `nmfreq` comparison in somatic `STATUS` classification (`seqc2-exome`):**
  At `chr3:195784364 T→C`, every numeric INFO and FORMAT field is identical between Java/Rust and cpp (`Tumor: DP=975, VD=49, AF=0.0503`; `Normal: DP=1021, VD=51, AF=0.05`, `SSF=0.52814`). Java rounds normal `v2.freq` to 4 decimals (`0.0500 >= nmfreq 0.05`) before `SomaticPostProcessModule.determinateType` and emits `STATUS=Germline`, whereas cpp compares unrounded `51 / 1021 = 0.0499510284 < 0.05` and misclassifies the variant as `STATUS=LikelySomatic`.
- **Indel/STR/`<DUP>` pileup & realignment count shifts:**
  - `hg008-tier2` (5 loci): `chr10:18036774 ATATATATATATATATATA→T` flips `STATUS=StrongLOH` (`VD=0` in tumor) to `STATUS=Germline` (`VD=15` in tumor) with an adjacent `chr10:18036792 ATA→TTT` depth shift (`DP=58` vs `73`); `chr1:72290009 C→<DEL>` emits `sample1 VD=254 > DP=155` in cpp (vs `DP=180, VD=178` in Java/Rust); `chr1:237117671` and `chr13:31976202` differ in normal `DP`/`VD`.
  - `seqc2-exome-tier2` (3 loci): `chr3:195779856 G→<DUP>` (`sample1 DP: 29` vs `422`), `chr3:195783263 T→<DUP>` (`StrongSomatic` `VD=0` vs `LikelySomatic` `VD=12`), and `chr8:8027922` (367 bp deletion where cpp emits a negative mismatch average `sample1 NM=-0.1`).
  - `hg002-tier2` (7 `<DUP>` SV loci): `QUAL`/`VD`/`AD`/`AF`/`ALD` count differences.
- **1-ULP (`1e-5`) Fisher exact rounding (`ODDRATIO`, `SSF`, `SOR`, `SBF`):**
  Native C++ `double` Fisher exact computation vs R's extended-precision `f80` `fisher.test()` accounts for **5/5** field diffs in `seqc2-slices-tier2` (`4,151 / 4,151` keys match), **941/945** diff records in `seqc2-exome-tier2`, and **45/50** diff records in `hg008-tier2`.
- **Counts:** `hg002` `227,562` vs `227,560`; `seqc2-exome` `623,444` vs `623,442`; `hg008` `18,426` vs `18,425`; `seqc2-slices` `4,151` vs `4,151` (0 key diffs).

Directionally this favors vardictrust's Java fidelity on these datasets,
but the absolute rates are small, mode-dependent (germline MNP vs somatic
SV), and paper-protocol data (hg19 WES, operating-point filtered) may
differ — no accuracy-rank claim is made from inter-caller comparison alone
(the Fig-3 no-truth caveat applies). Recorded advisory, never gating.

<a id="cpp-10"></a>
### CPP-10 — -DP parsed but ignored

Pinned vardictcpp accepts `-DP`/`--DP ERR` (listed under "Accepted (parsed,
VarDict-compatible; not all affect output yet)" in `vardictcpp --help`) but still
prints variants to stdout, so the M17 stderr comparison reports an advisory
`CPP-MISMATCH`; Java and Rust both honour it (M17 gate arm MATCH); advisory only,
no engine impact; reproduced with `-N` (without `-N` cpp exits 1 — that was a
harness omission fixed here).

