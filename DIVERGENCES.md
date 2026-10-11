# DIVERGENCES.md — VarDictRust Parity Tracking

Tracks every known deviation from VarDictJava v1.8.3 output.
Format: affected test(s), root cause, status, remediation.

> **History note (2026-08-23):** git history was squashed to a single initial
> commit before first deploy. Commit hashes cited below refer to the archived
> pre-squash line (`archive/pre-squash-20260823` branch); milestone evidence
> tags were re-pointed with identical messages.

**Principle:** Every bug, quirk, and accident in VarDictJava v1.8.3 is part of the
specification. We either match it exactly or prove the difference is unreachable.

---

## Master Parity Status & Taxonomy Summary

```
Verified Parity Baseline (release/v1.0.0, commit cb39a64):
- Integration Tests:  91 PASS / 1 FAIL / 0 SKIP (100% of valid inputs pass; sole FAIL is stale fixture UB-4)
- Core Library Tests: 188 PASS / 0 FAIL (cargo test -p vardict-core --lib)
- Total Workspace:    208 PASS / 0 FAIL (cargo test --workspace across core, cli, robustness, thread-determinism)
- Tier-2 Sweeps:      15/15 strict MATCH (s1a-s9 all 10 slices, HG008, SEQC2 exome/slices, HG002, Amplicon)
- Truth Concordance:  3/3 EQUAL (GIAB HG002 NIST v4.2.1, SEQC2 v1.2.1, GIAB HG008 NIST v3.3.2)
- Scaling Budget:     PASS (4.85x speedup, 1.39x flat RSS scaling across 249 Mbp chromosome 1)
- Final Acceptance:   8/8 PASS under release engine fingerprint 7e9b1be766f73b17
```

### The Four Canonical Parity Classes (rewrites.bio §5, §8, §11)

| Class | Designation | Tracked Items | v1.0.0 Parity Contract | Roadmap Stance (`ROADMAP.md`) |
| :--- | :--- | :--- | :--- | :--- |
| **Class A** | **Faithfully Emulated Upstream Quirks** | **UB-1** (`-p` pos-0 row), **UB-2** (`-p` null formatting), **UB-3** (HashMap iteration order via `JavaHashMap`), **UB-5** (SimpleMode SV split-read filter), **UB-6** (`markDUPSV` boundary), **UB-7** (single-winner collapsing) | **Byte-for-byte replication**. VarDictRust output matches VarDictJava v1.8.3 exactly on valid inputs (§5). | Candidate opt-in flags for v1.1.0+: `--fix ub1`, `--fix ub2`, `--fix ub5`, `--fix ub6`, `--fix ub7`. |
| **Class B** | **Accepted Divergences by Design** | **`##source`** header (`VarDictRust_v1.0.0`), **`.fai` contig lengths** (fixes #320), **CLI arg splitting** (`-f 0.001`), **`--fisher` raw floats**, **Splice-mode stdout order** | **Functional variant record parity**. Deviations confined strictly to VCF `#`-headers or CLI parsing ergonomics. | Permanent design standards; no flags required. |
| **Class C** | **Crash-Resilience & Robustness** | **PR #402** `modifyCigar` negative match length, **Issue #187** missing index fail-fast, **Issues #197 / #384** bounded halo tiling (2.9GB peak RSS) | **Zero crash, zero panic** (`#![forbid(unsafe_code)]`). Produces identical output on valid inputs; safely skips/recovers on malformed inputs. | Permanent engine baseline; `--fix 402` in v1.1.0+ will repair CIGAR topology. |
| **Class D** | **Known Non-Bugs / Upstream Flaws** | **UB-4** (`hard_clip_case.fa` conflicting CSV reference chunks), **GRCh37.fa** Colo829 chr20 all-N fixture drift | **Proven parity vs live Java oracle**. Confirmed that live Java under same inputs produces 100% byte-identical output (0 diffs). | Documented permanent fixture artifacts; not counted as regressions in the 91/1 corpus. |

*Upstream Bug Specification & Empirical Validation:*
- Detailed filing-ready draft specifications for UB-1 through UB-7 are consolidated in [`docs/upstream-bugs.md`](docs/upstream-bugs.md).
- Full empirical validation reports, biological provenance tables, algorithmic patterns, coordinate contracts, and filtering threshold catalogs are documented in [`docs/validation/omnibus-empirical-validation.md`](docs/validation/omnibus-empirical-validation.md).

---

## Category 1: CigarParser Mismatch Extension (process_deletion)

### Root Cause

When `process_deletion` handles a deletion followed by an insertion (`is_next_insertion`),
it extends the variant description by scanning the subsequent match region for mismatches.
The extension uses `vext` (default 2) to stop after N consecutive matches.

**Divergence:** Rust extends too far into the match region compared to Java. The mismatch
pattern between read and reference differs, causing different extension lengths.

### Affected Tests

- `Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1`
  **STATUS: Quirks Replicated (Live Oracle Parity: 100% MATCH).** Output is byte-identical to the live
  VarDictJava v1.8.3 oracle under canonical FASTA input (7/7 table rows), reproducing both `-p` null-record quirks:
  1. **Field-defaults row at 6447/6447:** In `SimplePostProcessModule.java:73-82`, a `Vars` entry at key 6447
     holds a single realigned variant with raw start 6449 != key and a null `referenceVariant`. Lines 76-77 print
     `SimpleOutputVariant(null, ...)` yielding empty genotype (col 14), bias "0;0" (col 16), and empty flanking
     sequences (cols 31-32; upstream bug UB-2).
  2. **Spurious all-zero row at 0/0:** Line 78 creates fallback default `new Variant()` with start=end=0, printed
     at lines 100-101 (upstream bug UB-1).
  Replicated in `crates/vardict-core/src/postprocess/simple.rs` (`make_null_record`) and `crates/vardict-core/src/output/table.rs`.
  **VCF Path Filtering:** In Java's pipeline, `var2vcf_valid.pl:141` (`next unless ($ref)`) drops records with
  empty/falsy REF, discarding both null rows so only 5 variants reach final VCF. VarDictRust ports this filter
  via `output::vcf::perl_row_skipped`, achieving byte-identical 5/5 VCF variant MATCH.
  **Fixture Inconsistency:** The upstream static fixture file was generated from an older test CSV containing
  two conflicting, overlapping chunks (`test,5234..7764` vs `test,5475..8024` with 1734/2290 mismatching bases)
  under first-chunk-wins precedence. `convert_fastas.py` uses last-write-wins. Against live VarDictJava with the
  same reference, output is 100% identical.
- `Simple;hg19.fa;Colo829.chr10.LGDEL3.bam;chr10;127187158-127201055;-f 0.001;`
  **RESOLVED (W1.2)** — root causes were NOT the process_deletion extension scan itself but
  two find_inv_sub divergences (see Resolved table): wrong INV desc tail (`revcomp(consensus)`
  instead of findMatchRev's `extra` suffix) and marking a CLONED Sclip used instead of the
  live map entry.
- `Simple;hg19.fa;Colo829.chr20.DUPDiscSV.bam;chr20;47124929-47131802;-f 0.001;`
  **RESOLVED (W1.2)** — partial_pipeline merged soft clips with `or_insert` (discarding
  counts) where Java's shared maps accumulate across overlapping partial pipelines; the
  resulting cnt=1 (< minr=2) cluster at 47131875 suppressed the +342 lgins insertion and the
  SV-info fields. See Resolved table.
- `Simple;hg19.fa;Colo829.chr20.INVSV.bam;chr20;46841555-46845755;-f 0.001;`
  **RESOLVED (W1.2)** — wrong INV desc tail made realignIndels' findMM5 miss the T@47214669
  absorption (Java cnt 2→7, cov 12→17 via Realigndel Adj; Rust stayed 2/12).

### Investigation Notes

- CigarModifier correctly transforms `1D13M1D` → `15D13I` (verified)
- The `is_insertion_or_deletion_with_next_matched` function matches Java (vext threshold)
- The mismatch extension loop logic matches Java (`vsn > vext` break condition)
- **Key finding:** Reference at positions 6484-6506 is `TCACCAGATGGCTACCATTCCAG`
  but Java's expected output implies different mismatch pattern — needs instrumented
  Java oracle trace to compare internal state

**threeIndels hi_freq @ 62360386 (Java 0.8936 vs Rust 0.9333) — root cause pinned:**
hi_freq = hicnt/hicov for the `+T` insertion. hicnt=42 in both; hicov differs
(Java 47, Rust 45). hicov = sum of non-insertion variants' hicnt at the position.
Instrumented Java ToVarsBuilder shows Java has a **T SNV at 62360386 with
varsCount=2, hq=2** that Rust does NOT register (Rust non-insertion at 62360386:
C cnt=50/hicnt=45, A cnt=1/hicnt=0; no T). The missing T's hicnt=2 is the entire
hicov gap. The 2 T reads (tumor BAM, chrom `20`, both qual 37 / MAPQ 60,
non-duplicate, surviving `-F 0x700`): `A00130:81:HFFWLDSXX:2:1307:7563:3662`
(pos 62360262, `140M11S`, T at qpos124) and `A00130:81:HFFWLDSXX:3:2520:15374:21136`
(pos 62360380, `150M`, NM=3, T at qpos6). Rust's pileup/MNV path does not register
their T mismatch at 62360386 (the surrounding positions 62360380-85 are T-reference
with cnt~100; 62360386 is C-reference).

**UPDATE — localized to `realign_indels`:** a READTRACE in Rust
`add_variation_for_matching_part` shows BOTH reads DO register `s=T` at pos 62360386
(q=37, qbases=1) during pileup. A STAGE_TRACE across the BAM1 pipeline confirms:
T@62360386 cnt=**2 after pileup**, cnt=**2 after adjust_mnp**, cnt=**None after
realign_indels** (stays None through find_all_sv and adj_snv). So **`realign_indels`
is the exact stage that removes the 2-read T SNV** (Java keeps it). The T sits
immediately 5' of the `-1`/`-2`/`-3` deletions at 62360387. Next step: instrument
`realign_indels` (realign.rs) for position 62360386 to find which indel-realignment
path absorbs/deletes the T variant — likely the deletion realignment at 62360387
consumes the adjacent T mismatch reads in Rust but not in Java.

**Mechanism identified:** the removal is `realigndel`'s misp5 branch
(realign.rs:1014-1031): `if misp5 != 0 && mm5.len() == 1 && variant.cnt < dcnt → remove`,
**byte-identical** to Java `VariationRealigner.java:550-553`.

**CORRECTION — deletion counts match:** instrumented Java ToVarsBuilder shows the
deletion counts at 62360387 are IDENTICAL in Java and Rust (`-1`=6, `-2`=1, `-3`=3,
`T`=101, `C`=1). So the divergence is NOT the deletion count. With T.cnt=2 < dcnt=6
the misp5 branch would fire in BOTH implementations — therefore Java must compute a
different `misp5`/`mm5` for the `-1` deletion at 62360387 so the branch does not
target the T there. **The real divergence is realigndel's 5'-mismatch finding
(`findMM5` / `getMisp` / `mm5`) for the deletion at 62360387:** Rust resolves
`misp5=62360386, misnt5=T, mm5.size()==1` (removing the T), Java resolves it
differently (keeping the T). Next step: instrument realigndel's mismatch-finding for
the `-1` deletion at 62360387 in both implementations and diff the misp5/mm5/misnt5
values — the mismatch scan boundary or the reads feeding it differ.

**Bottom of chain — `find_mm5`:** `misp5`/`misnt5` come from `find_mm5`
(realign.rs:417). Its `misp`/`misnt` are set (realign.rs:465-525) only when the
single-base boundary mismatch (`str_bytes.len() == 1`, here the T) is followed by
`mn > 1` matching bases and then `n2 > 2` further matches, scanned over the filtered
consensus `seq` (wupseq minus `#`/`^`). The Java/Rust divergence is inside this scan:
either the consensus `wupseq` feeding `find_mm5` differs for the `-1` deletion at
62360387, or the `mn`/`n2` counting/boundaries differ. Next step: trace `find_mm5`
inputs (`wupseq`) and intermediates (`n`, `mn`, `n2`, `misp`) for the `-1` deletion
at 62360387 in both implementations and diff them.

### Remediation

Run VarDictJava with `-y` verbose on the hard_clip test case. Compare the exact
mismatch positions and extension lengths. The difference is likely in how
`read_pos_inc` / `multoffp` offsets are calculated after the insertion.

---

## Category 2: realignlgdel Soft Clip Consensus

### Status (updated W1.3, HEAD post-W1.2)

RESOLVED for 4 of 5 listed cases; remainder re-classified as Category 11:

- The three Colo829 cases (chr1 755917-756517, chr10 1368753-1369353,
  chr10 684892-685492) were fixed by W1.2's live `merge_sclip` partial-pipeline
  soft-clip sharing (a7870b4+c0f36a6), NOT by realignlgdel consensus changes.
  Re-verified MATCH at c0f36a6.
- `Zero_in_DUPSV chr19 290742-291636` was NOT a realignlgdel issue at all. Trace
  at HEAD: tumor lacks the variant entirely (StrongLOH path); Rust's existing
  second-sample `combine_analysis` correctly returns "" (keep) — combined run
  finds `-13^GCGGGGAGGGGAGAGG` @291563 with vd=4 vs normal vd=4. The sole
  divergence was output anchoring: `print_variations_from_second_sample` built
  the placeholder variant1 by copying alleles/positions from v2var BEFORE
  `adj_complex()`, while Java applies `adjComplex()` before constructing
  `SomaticOutputVariant(v2var, v2var, varForPrint, v2var, …)`
  (SomaticPostProcessModule.java:329-333). Fixed in W1.3 by moving the
  adj_complex call ahead of the fallback construction; regression test
  `test_strong_loh_complex_adj_before_placeholder`.
- `HCC2218_BL chrX 7811670-7811868`: spurious StrongSomatic +30 insertion
  @7811769. Java `-y` trace at HEAD shows exactly one
  `Start Combine 7811769 +GAACCACTGAGACAGGAGAGCGAGGTGGAA` with NO following
  `Combine:` line → combined-pipeline lookup miss → `vref==null` → FALSE →
  variant skipped (SomaticPostProcessModule.java:408-409,202-204). This is the
  reverted combineAnalysis path in `printVariationsFromFirstSample` =
  **Category 11**, not a bounded fix. See Category 11.

### Root Cause

The 5' section of `realignlgdel` creates deletions from soft clips. It:
1. Builds consensus via `find_conseq` from the `Sclip.nt` map
2. Searches for breakpoint via `find_bp` (within indel_size=50) or `find_match` (SEED=17)
3. Creates deletion variant and absorbs matching 3' soft clip (no min_reads gate on 3' clip)

**Divergence:** Our `find_conseq` produces shorter consensus (11bp) than Java needs (≥17bp
for `find_match`). The 5' soft clip at position 756260 in the Colo829 tumor BAM has only
2 reads (7S and 11S), producing an 11-base consensus that fails both `find_bp` (wrong
match position) and `find_match` (11 < SEED=17).

### Affected Tests

- `Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr1;755917-756517;`
- `Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr10;1368753-1369353;`
- `Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr10;684892-685492;`
- `Somatic;hg19.fa;HCC2218_BL_realigins.bam-HCC2218_C_realigins.bam;chrX;7811670-7811868;-f0.001;`
- `Somatic;hg19.fa;Zero_in_DUPSV_TUMOUR.bam-Zero_in_DUPSV_NORMAL.bam;chr19;290742-291636;-f0.001-X3`

### Investigation Notes

- Normal BAM creates "-40^CACTA" at 756220 from read with CIGAR `68M9D5M26D77M` (pileup)
- Tumor BAM has NO such read; Java trace confirms tumor's realigndel only processes 35D at 756286/756278
- **CRITICAL FINDING (Java trace):** The "-40^CACTA" at 756220 appears ONLY in the normal BAM's
  initial realigndel. The tumor BAM never creates this deletion in pileup or realignment.
  Yet Java's final output shows tumor VD=2, STATUS=Germline. The mechanism is unknown —
  likely in ToVarsBuilder or SomaticPostProcessModule (possibly `combineAnalysis` or
  variant lookup across samples).
- 5' soft clip at 756260: count=2, consensus="ATCACACCTAG" (11bp)
- `find_bp` fails: consensus doesn't match reference within 50bp upstream
- `find_match` fails: 11 < SEED=17
- `realignlgins` skips: 11 < 12 (minimum length)
- **Key mechanism (confirmed):** 5' section of realignlgdel creates deletion AND absorbs
  3' clip at sc3p via `adjCnt` with NO min_reads check → but this path is NOT taken here
  because find_bp/find_match both fail

### Remediation

Investigate Java's `ToVarsBuilder.process()` and `SomaticPostProcessModule` to find how
the tumor gets VD=2 for a variant that doesn't exist in its `non_insertion_variants`.
Run Java with `-y` and trace the ToVarsBuilder output for position 756220 in BAM1.

---

## Category 3: NA24385 Low-AF Filtering — RESOLVED (W2.2)

### Status

CLOSED. All six listed cases now MATCH end-to-end.
Corpus 88/4 → **90 PASS / 2 FAIL** (remaining FAILs are the out-of-scope
Colo829 amplicon `-D` case and the stale hard_clip fixture, see Cat 1/Cat 9).

### Root causes found and fixed (five distinct defects, all in the SV/soft-clip
collection path that low `-f` thresholds expose)

1. **Chimeric-clip seed check skipped for ≥20bp clips without an SA tag**
   (`process_soft_clip`, both 5' and 3' branches). Java:
   `if (len >= 20 && saTag != null) { SA check } else if (len >= SEED_1) { seed check }`
   — the seed check runs for ANY clip ≥17 without an SA tag. Rust's nested
   `if len >= 20 { if let Some(sa) … } else if len >= 17` structure skipped the
   seed check whenever len ≥20 and no SA tag existed. Effect: chimeric
   soft clips were registered, their clusters later produced spurious INV/DEL
   SVs (the C1 Complex `-62i` @38396474 record; C2's three extra INVs and
   Complex). Evidence: instrumented-Java verbose prints
   `AGTGTTTTTTCCCATATTAGGTA at 3' is a chimeric at 38396367 by SEED 17` for
   reads `20S107M23S`; Rust registered those clips instead.

2. **SV structures received raw NM instead of NM − (I+D)** (`prepare_sv_structures_for_analysis`
   call site). Java: `numberOfMismatches = NM − insertionDeletionLength`
   (CigarParser.java:227-236). Rust passed `record.nm.unwrap_or(0)`.

3. **Inter-chromosomal fusion branch filtered on the read's OWN MAPQ** and
   returned early, skipping Java's post-branch disc bumps. Java filters on the
   MC tag (mate cigar with two-sided clips) and the mate MQ TAG when present
   (CigarParser.java:2208-2215); the read's own MAPQ only gates via the outer
   `> 10` check. Same family as commit 676d6c6 which fixed the same-chromosome
   branch but missed this one.

4. **Inversion-tail discordant bumps used `.abs()`** where Java compares
   signed `(start - lastEnd) <= MIN_D` (CigarParser.java:2190-2202). The
   signed form fires far more often (any negative distance passes), so Rust
   systematically under-counted `disc`, flipping `filterSV`'s false-positive
   gate (`cnt/disc < 0.35 && !(…>=0.35 && cnt>=5)`) and leaving clusters
   unmarked-used — directly producing the spurious `<DUP919>` @44940571
   (`svrdup[55]`: Java disc=61→used, Rust disc=48→unused). This applies ONLY
   to the `mlen != 0` tail of the FR-inversion branch (:2190-2202); the
   sibling DEL-branch tail bumps all use `abs(start - lastEnd)` in Java
   (:2058-2077, e.g. the dup bump at :2064) and must stay absolute — a W2.2
   fix-round initially over-applied the signed form there and reverted.

5. **Integer truncation in realignins' adjRefFactor bound**
   (`realign.rs`). Java: `first3 - first5 < maxReadLength * 0.75` compares
   int against a double; Rust truncated the RHS to i32 first
   (`((R as f64) * 0.75) as i32`), so fractional bounds (e.g. 112 < 112.5)
   failed. This suppressed the ±factor reference-count adjustment at the +110
   insertion @44953220 (Java vd=8/hicov=19 vs Rust vd=5/hicov=72).

### Method / evidence

Fresh oracle runs (`vardict-java -z -y`) diffed at the `.var` level per case;
instrumented classpath-shadow builds of VarDictJava 1.8.3
(`javac` of a patched source copy shadowing `VarDict-1.8.3.jar`) dumped every
`addSV`/`adddisccnt` call and the svfdup/svrdup cluster lists at
`findDUPdisc` entry; equivalent Rust traces diffed line-by-line until the ADDSV
streams matched byte-for-byte and only condition-level divergences remained.
All five fixes verified against live Java output; unit tests added for fixes
1-4 (mutation-checked), fix 5 is covered by the corpus flip itself.

### Affected Tests (all now MATCH)

- `Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;20156208-20234556;-f0.001-fisher`
- `Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;20156208-20234556;-f0.001`
- `Simple;GRCh37.fa;NA24385_chr_20_36018823-36039513.bam;chr20;36018823-36039513;-f0.001`
- `Simple;GRCh37.fa;NA24385_chr_20_38386265-38409801.bam;chr20;38386265-38409801;-f0.001` (fixed by W2.2)
- `Simple;GRCh37.fa;NA24385_chr_20_44885701_44969115.bam;chr20;44885701;44969115;-f0.001` (fixed by W2.2)
- `Simple;GRCh37.fa;NA24631_sorted.chr1.1000k.bam;chr1;1;499702;-f0.001`

### Historical notes (pre-W2.2, superseded)

- Oracle script fixed: `-f0.001` now correctly split to `-f 0.001` for VarDictJava CLI
- Earlier "4 java-only + 10 rust-only" diffs predate W1.2/W1.4/W2.1 upstream
  changes (partial_pipeline merge, find_inv_sub early-return, combineAnalysis);
  fresh re-diff at W2.2 start showed 1 rust-only (C1) and 1 java-only + 6
  rust-only (C2) before the fixes above.

---

## Category 4: Somatic Post-Processing — RESOLVED

### Status

CLOSED. All six listed cases now MATCH end-to-end against the oracle:
- `threeIndels` fixed by `b3e2c38` (endVariant source) + `8fda935` (findMM5 dead loop replication).
- `ArgumentNA` (fisher and non-fisher) fixed by `bc93c0c` (-M minmatch separation).
- `TumorA.StrongLOH` fixed by `fe5962d` (DELETION zeroed tumor placeholder).
- `NA12878.part1` fixed by `c6269e5` (somatic_61 formatter for --fisher).
- `NA12889.StrSomSamSpec` fixed by `f7af5b7` (is_noise mutation in determinate_type).

### Affected Tests (all now MATCH)

- `Somatic;GRCh37.fa;threeIndels_zero_deletion_tumor.bam-threeIndels_zero_deletion_normal.bam;chr20;62354864-62374864;-f0.1-Q10-F0x700;`
- `Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;`
- `Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p;`
- `Somatic;hg1k.fa;TumorA.StrongLOH.bam-NormalA.StrongLOH.bam;chr21;10996701-11186054;-Q10-f0.01-r2-V0.05;`
- `Somatic;hs37d5.fa;NA12878.part1.1956957-1957177.bam-NA12889.part1.1956957-1957177.bam;chr1;1956957-1957177;-f0.001;`
- `Somatic;hs38DH.fa;NA12889.StrSomSamSpec.bam-NA12878.StrSomSamSpec.bam;chr20;43246926-43280248;-f0.001`

### Historical Investigation Notes (Resolved)

The notes below are preserved for provenance:
- Somatic pipeline structure fixed (added find_all_sv + adj_snv for both BAMs).
- Differences traced to upstream CigarParser/realigner and post-processing mutation propagation.
- STATUS determination logic (determinate_type) validated and verified.

---

## Category 5: CigarParser 3' Soft Clip Creation (Root Cause of VD Differences)

> **RESOLVED (W1.1).** The investigation below is retained for provenance. At HEAD,
> `adjust_mnp` DOES mark the 62358591 3' soft clip `used=true` (trace-proven), so
> realigndel's 3' loop skips it exactly as Java's does — the "issue in realigndel's
> 3' loop" hypothesized below was that its precondition (`!used`) now correctly holds.
> Regression test: `realign.rs::test_adjust_mnp_marks_3p_sclip_used`.
> The VD differences originally catalogued here were closed by W1.1/W1.2 (see Resolved
> table); remaining VCF-level residues are tracked in Category 14.

### Root Cause

The 1-2 read VD count differences across all failing tests trace to a single mechanism:
Rust's CigarParser creates 3' soft clips that Java's realigndel doesn't process.

**Specific example (threeIndels, position 62358581):**
- Read `146M2S` at position 62358445 creates a 3' soft clip at 62358591 (2 bases, "CT")
- Rust's realigndel processes this soft clip (adds 2 reads to the -1 deletion → VD=32)
- Java's realigndel does NOT process it (VD stays at 30)
- The "Mismatches" line shows "sclip3: 62358591" in both Java and Rust
- But Java's 3' loop skips it (either `softClips3End` doesn't have it, or it's marked `used`)

**Investigation results (updated with Java instrumentation):**
- CigarModifier's `capture_mis_softly_ms` correctly does NOT convert 2S→M (ref='T' ≠ query='C')
- Both Java and Rust have `while (flag && indel > 0)` loop; `captureMisSoftlyMS` is OUTSIDE the loop
- Quality is 37 (> LOWQUAL=10), so quality filtering doesn't explain the difference
- **CONFIRMED via Java instrumentation:** Java DOES store the soft clip at 62358591
  (SCLIP3_TRACE: start=62358591 hq=2 lq=0 region=62354864-62374864 cigarLen=2 readPos=146)
- The CigarModifier returns "146M2S" unchanged (CIGMOD_2S trace confirmed)
- The read passes all filters (MAPQ=60, NM=0, flags=0x93)
- **The issue is in realigndel's 3' loop:** despite the soft clip existing in softClips3End,
  the loop doesn't process it. Likely `findconseq` returns empty or `ismatch` fails
  for a subtle reason (possibly the `nt` map iteration order or consensus quality check)

### Remediation

**ROOT CAUSE FULLY IDENTIFIED:**

Java pipeline order: `filterAllSVStructures()` → `adjustMNP()` → `realignIndels()`

`adjustMNP()` (VariationRealigner.java line 346-357) iterates over `nonInsertionVariants`.
For each position with a variant, it checks `softClips3End.containsKey(position)`.
If the soft clip's consensus starts with the variant's mismatch string (`mnt`), and the
remaining consensus matches the reference (`ismatchref`), it marks `sc3v.used = true`.

This prevents `realigndel`'s 3' loop from processing the soft clip (since it checks
`!softClips3End.get(sc3pp).used`).

**In Rust:** `adjust_mnp` (realign.rs:4212-4240) has the same logic but doesn't
mark the soft clip at 62358591. Debug why:
1. Check if `non_insertion_variants` has an entry at position 62358591
2. If yes, check what `mnt` (variant description) is at that position
3. Check if `find_conseq` returns a consensus that starts with `mnt`
4. Check if `is_match_ref` succeeds for the remaining consensus

The fix is likely a subtle difference in how `adjust_mnp` iterates or how the
variant description (`mnt`) is extracted at position 62358591.

---

## Category 6: SV Info Field Count Difference

### Status (re-checked W3.3)

CLOSED. All three listed residues are resolved on the NA24385
chr20:20156208-20234556 case, re-diffed at W3.3 against the checked-in
expected table (corpus PASS): pairs at 20156480 now `3-18-2` (= Java, was
`3-17-2`, fixed by 676d6c6); the INV at 20158414 is a single `<INV>` line
with cov 73 (= Java's DP=73, was duplicated VD=6/DP=88) and SV field
`3-3-1` (= Java, was `6-3-1`, fixed by 3c5a6ce + W1.4). No remaining
divergence to trace.

### Root Cause

The SV info string format is `splits-pairs-clusters`. Multiple differences observed:

1. **Pairs count (middle value):** `3-18-2` vs `3-17-2` at position 20156480.
   1-count difference in discordant pair clustering.

2. **Splits count (first value):** `3-3-1` vs `6-3-1` at position 20158414 (INV).
   Rust detects 6 split reads vs Java's 3. This causes Rust to output 2 variant
   lines for the same INV (VD=3 and VD=6) while Java outputs 1 line (VD=3).

3. **DP difference at 20158414:** Java DP=73 vs Rust DP=88 (15 reads difference).
   Related to the SV splits count difference - Rust counts more SV-supporting reads.

### Soft Clip Filter Fix (dd43d8a)

Java's regex `^\d\dS.*\d\dS$` only drops reads where:
- First soft clip is exactly 2 digits (10-99bp)
- Last soft clip is at least 2 digits (>= 10bp)

Previous Rust code incorrectly dropped reads with first soft clip >= 100bp.
Fixed to match Java's behavior. This resolved the DP difference at position 20178200
(Java DP=387 vs Rust DP=268 → now matches).

### Remediation

Investigate SV detection logic (find_inv) to understand why Rust detects more
split reads than Java for the INV at 20158414. The 2x difference (6 vs 3) suggests
a systematic issue in how split reads are counted or filtered.

---

## Category 7: Complex Variant Position/Left_seq Shift

### Root Cause

Complex variants (e.g., `G/-13iGCGGGGAGGGGAGAGG`) show a 1-base position difference
between Java and Rust. This causes:
- left_seq to be shifted by 1 base (Java: `AAGAGAGCAAGACTCTACTG`, Rust: `AGAGAGCAAGACTCTACTGC`)
- start_position to differ by 1 (Java: 291564, Rust: 291563)

### Affected Tests

- `Somatic;GRCh37.fa;threeIndels_zero_deletion_tumor.bam-threeIndels_zero_deletion_normal.bam;chr20;62354864-62374864;-f0.1-Q10-F0x700;`
- `Somatic;hg19.fa;Zero_in_DUPSV_TUMOUR.bam-Zero_in_DUPSV_NORMAL.bam;chr19;290742-291636;-f0.001-X3`
- `Somatic;hg19.fa;HCC2218_BL_realigins.bam-HCC2218_C_realigins.bam;chrX;7811670-7811868;-f0.001;`

### Investigation Notes

**ROOT CAUSE FOUND (traced).** The shift is NOT in `adj_complex` (its left_seq
logic matches Java exactly). It is in which variant supplies the somatic output's
left/right seq + shift3/msi/msint.

Concrete evidence (threeIndels Complex variant, map pos 62361370, desc `+CACACA&CA`):
- tovars (pre-adj_complex): start=62361371, ref=CG, var=CACACACA,
  left_seq=`AAGAGAGCAAGACTCTACTG` (= Java's final output).
- `adj_complex` prefix-trims n=1: start→62361372, ref→G, var→ACACACA,
  left_seq→`AGAGAGCAAGACTCTACTGC` (shifted right 1). This is the Rust output.
- Java outputs the UN-shifted `AAGAGAGCAAGACTCTACTG`.

Why: Java `SomaticOutputVariant` takes startPosition/endPosition/refAllele/varAllele/
varType from `beginVariant` (=vref, the tumor, adj_complex'd) but takes
shift3/msi/msint/leftSequence/rightSequence from `endVariant`. In
`printVariationsFromFirstSample`'s `if (v2nt != null)` branch the constructor is
`SomaticOutputVariant(vref, v2nt, vref, v2nt, ...)`, so **endVariant = v2nt (the
normal variant), which is NOT adj_complex'd** — hence the un-shifted left_seq.
Rust `somatic_55`/`somatic_61` read left_seq/right_seq/shift3/msi/msint from
`variant1` (tumor, adj_complex'd) → shifted. For non-Complex variants tumor and
normal left_seq coincide, so only Complex variants diverge.

endVariant source by branch (Java constructor 2nd arg):
- both-samples `if v2nt`: v2nt = **variant2**
- both-samples else (StrongSomatic, varForPrint placeholder): vref = **variant1**
- one-sample DELETION `(variant,variant,null,variant)`: variant = **variant2**
- one-sample SAMPLE_SPECIFIC `(variant,variant,variant,null)`: variant = **variant1**
- LOH `(v1nt,v2var,v1nt,v2var)`: v2var = **variant2**
- StrongLOH `(v2var,v2var,varForPrint,v2var)`: v2var = **variant2**

### Remediation

Add an `end_is_variant2: bool` flag to `SomaticVariant` (true for the four
variant2 branches above, false for StrongSomatic-else and SAMPLE_SPECIFIC), and
have `somatic_55`/`somatic_61` take left_seq/right_seq/shift3/msi/msint from
variant2 when the flag is set, else variant1. Note: fixing this alone does not
pass threeIndels (a separate hi_freq difference at 62360386 remains — upstream
highQualityReadsCount/coverage divergence).

---

## Category 8: Duplicate INV with Different Insertion Sequences (20201754)

### Root Cause

Rust outputs TWO INV variants at position 20201754 (both display as `T/<INV303>`),
while Java outputs ONE. The two Rust variants have different insertion sequences:
- Variant 1: `-303^CTACTTGTTTGGGGGCTATTATTTAATTGGAAACAAGTG...<inv203>...` (extra 39-base prefix)
- Variant 2: `-303^TCACTTGTTTCC...<inv203>...`

Key difference: Java's SV info is `3-0-0` (pairs=0, clusters=0), meaning Java's
find_inv_sub does NOT process this INV. Rust's SV info is `6-2-1` (pairs=2, clusters=1),
meaning Rust's find_inv_sub DOES process it. The detection mechanism differs — Java
detects this INV only via a splits-only path (likely realignlgdel), while Rust also
detects it via find_inv_sub.

### Investigation update (this session)

Tested and RULED OUT the used-marking path: marking the soft clip map entry used
(Java line 663 — Rust sv.rs:1135 only sets a clone) AND marking the INV cluster used
(Java line 644 `inv.used = true`, which Rust never does) does NOT remove the duplicate.
The two variants arise from two different `softp` clusters that converge on the same
breakpoint via the `_5`/`_3` side-swap (bp/softp exchange, Java line 605-608) and build
two different insertion consensus sequences — it is genuinely two consensus builds, not
a re-processing of one cluster. Also noted: Rust's find_inv_sub gate (sv.rs:850) uses
`is_best_sv_used` whereas Java checks `scv.used` (line 534). The real fix must address
why Rust builds two consensus sequences (two soft-clip clusters survive) where Java
builds one — likely a difference in which soft-clip clusters qualify (cnt/used/consensus
length ≥17) or how `find_conseq`/`find_match_rev` resolves the breakpoint.

**Java `-y` pinpoints it:** Java's `findINVsub` produces exactly ONE INV at 20201754
— consensus `-303^TCACTTGTTTCCAATTAAATAATAGCCCCCAAACAAGTAGCTCATATGAT<inv203>...`,
Cnt:3 (the "Found INV: 20201754 BP: 20202056 Cov: 109 109 ... EXTRA:" trace). Rust
produces that SAME variant PLUS an extra `-303^CTACTTGTTTGGGGGCTATT...<inv203>`
variant. So the fix is narrow: find why Rust's `find_inv_sub` processes one additional
soft-clip cluster that Java's `findINVsub` skips — a difference in the consensus-length,
`find_match_rev` bp==0, or the `bp > softp && bp-softp > 150 && (bp-softp)/|mlen| < 1.5`
gate (Java line 627) for that one cluster. Trace Rust `find_inv_sub` for the cluster
yielding the `CTACTTGTTTGGGGGCTATT` consensus and compare its seq/bp/mlen against the
Java gate.

**Rust `INVSUB_TRACE` result:** the extra `CTACTTGTTTGGGGGCTATT` variant is created by
`find_inv_sub` from svfinv5 (svref_idx=0, dir=1, softp=20201754, bp=20202056, len=303,
mlen=520, cnt=2). It PASSES Rust's gate (303/520 = 0.58 < 1.5), and Java's identical gate
(line 627) would also pass it — so Java skips this cluster EARLIER, via the `used` flag
set by a prior INV detection across the four ordered `findINVsub` calls (svfinv5, svrinv5,
svfinv3, svrinv3), or a clustering count difference (Java may cluster cnt<2 here). The
Java-matching `TCACTTGTTTCC` variant is produced by a SEPARATE Rust path (not find_inv_sub).
Fix direction: reproduce Java's cross-call used-marking order so the svfinv5 cluster at
20201754 is already `used` when `find_inv_sub` reaches it (or match the cluster cnt that
Java computes for it).

**Ordering narrowed further:** Java `findINV()` calls `findINVsub`×4 FIRST (lines
481-484) and the 5'×5' nested loop AFTER (1138+) — the SAME order as Rust
(`find_inv_sub` at sv.rs:1237-1265 before the 5'×5' loop at 2631+). So the nested
loop's `invf5.used = true` (sv.rs:2714, matching Java 1196) runs AFTER `find_inv_sub`
and cannot block it. The 5'×5' pairing conditions match Java exactly (DISCPAIRQUAL>35,
cnt+rcnt>minr+5, is_overlap). Therefore the svfinv5[20201754] cluster must already be
marked `used` by **`findDEL`** (Java line 111, which runs before `findINV` at 115) via
`markSV`/`is_best_sv_used`. The fix is to find why Rust's `find_del` (or the
`mark_sv`/`is_best_sv_used` flow) does not mark the svfinv5 cluster at 20201754 used
before `find_inv_sub` reaches it — trace `mark_sv` calls in `find_del` for positions
overlapping 20201754 in both implementations.

**ROOT CAUSE (final):** `find_del` only marks deletion structures (svfdel/svrdel — Java
markSV at 228/299/391/464), NOT INV clusters, so the used-marking hypothesis is wrong.
The actual cause is a **soft-clip consensus difference**: the SAME `findINVsub(svfinv5)`
at softp=20201754 builds Java consensus `TCACTTGTTTCCAATTAAATAATAGCCCCCAAACAAGT...` but
Rust builds `CTACTTGTTTGGGGGCTATTATTTAATTGGAAACAAGTG` + that same Java consensus (i.e.
the Java consensus with an extra ~39-base prefix). So Rust's soft-clip cluster at
20201754 yields a LONGER `find_conseq` consensus than Java's — Rust includes extra
reads in the cluster or extends the consensus further. The fix is in the soft-clip
clustering / `find_conseq` for the 20201754 cluster (compare which reads populate
`soft_clips_5_end[20201754]` / `svfinv5` and how `find_conseq` builds the consensus in
both implementations).

**Bottom layer:** `Sclip.nt` is `BTreeMap<i32, BTreeMap<u8,i32>>` — sorted, matching
Java's `TreeMap<Character,Integer>` (same A<C<G<T iteration order), so the consensus
difference is NOT from map iteration/tie-breaking. The longer Rust consensus therefore
comes from **different reads populating the soft-clip cluster at 20201754** — the
divergence is one level deeper, in pileup soft-clip clustering (which reads are grouped
into `soft_clips_5_end[20201754]` and how their clipped sequences populate `nt`). Next
step: dump the reads/`nt` counts in `soft_clips_5_end[20201754]` in both implementations
and find the read(s) Rust includes that Java excludes (or vice versa) — likely a
difference in how soft clips are assigned to a position during `add_variation`/pileup.

### RESOLVED (W1.4): Java findINVsub early-return

After the W1.1/W1.2 fixes the case had collapsed to a single differing record
(`T/<INV303>` @20201754: Rust DP=92/SV=`3-2-1` vs Java DP=109/SV=`3-0-0`; desc and
splits already matched). Read-level dumps of the four INV cluster lists via an
instrumented-Java classpath override (`javac` of a /tmp source copy shadowing
`VarDict-1.8.3.jar`) showed:

- At `findINV()` time Java's svfinv5 contains BOTH `20157026` (soft={20157096=2}) and
  `20201601 end=20201730 mstart=20202120 mlen=520 varsCount=2 soft={}`. But
  `findINVsub(svfinv5)` iterates only the FIRST entry — it builds the INV@20157096 and
  hits **`return vref;`** (StructuralVariantsProcessor.java:675). The `20201601`
  cluster is NEVER considered by any findINVsub call.
- Rust's `find_inv_sub` looped through ALL clusters. It reached `20201601`, took the
  softp==0 scan branch (SC3@20201755, varsCount=3), built the SAME `-303^...<inv303>`
  INV at softp=20201754 that Java's *nested* loop later builds, marked the soft clips
  used, and stamped `sv.pairs += 2 / sv.clusters += 1` plus
  `refCoverage[20201754] = refCoverage[20201753] (=92)`.
- Java's record instead comes from the nested 3' loop ("Found INV:" trace, Cov: 109 109,
  Cnt: 3) which sets `pairs+=0/clusters+=0` and `incCnt(cov)+max(cov(bp))=109`.

Fix: `find_inv_sub` now returns after its first successful detection (sv.rs, ports the
Java `return vref;`). Regression test `test_find_inv_sub_returns_after_first_success`
constructs two viable clusters and asserts the second is untouched (verified to fail
without the fix). Case flipped PASS (both -f0.001 variants); corpus 85→87 PASS / 5 FAIL,
no new FAILs. Cat 13's splits over-count + duplicate INV symptoms on this same case are
resolved by the same fix.

### Latent divergences (RESOLVED in W3.3)

Rust-only behaviors in the `find_inv_sub` path with no counterpart in Java's
findINVsub were investigated and eliminated in milestone W3.3:

1. **mark_sv×4 sweep + shifted-softp used-marking in the find_inv_sub success
   path** — Rust marked the soft clip at the (side-swapped/shifted) `softp`
   used and ran `mark_sv` over all four INV cluster lists (sv.rs:1177-1199).
   **RESOLVED (W3.3a): genuinely extra — REMOVED.** Java's findINVsub success
   block performs exactly three mutations before `return vref;`:
   `inv.used=true` (:644), `scv.used=true` on the LIVE fetched entry (:663,
   already ported by the W1.2 scv_used_pos fix), and the embedded
   `realigndel(bams, dels5)` (:665-667). Traced every candidate Java path:
   `VariationRealigner.realigndel` (VR:407-968) contains **no markSV call**
   (its only `used` writes are soft-clip map entries at VR:573/:604, ported in
   Rust realign.rs); the markSV sites live in realignlgdel (VR:1032/:1244,
   DEL lists, runs earlier inside realignIndels) and in
   findDEL/findDELdisc/findINVdisc (SVP:228/299/391/464/1052/1123/1222/1313 —
   before findINV or after it with their own pair-derived coordinates).
   findAllSVs order is findDEL→findINV→findsv→findDELdisc→findINVdisc→
   findDUPdisc (SVP:107-132), so no Java path reproduces the sweep. The extra
   marking could only suppress detections Java performs (missing-record
   class). Regression test
   `test_find_inv_sub_success_marks_only_java_equivalents`
   (side=3/dir=1 swap scenario; mutation-checked to fail on the old code)
   pins the Java behavior.
2. **is_best_sv_used gate at sv.rs:850** — Rust skipped an INV cluster whose
   best soft-clip SV was already used. **RESOLVED (W3.3b): genuinely extra —
   REMOVED.** Java findINVsub has no SOFTP2SV check: its only skip is
   `scv.used` on the fetched clip entry (:534; :561 in the scan branch). The
   `SOFTP2SV.get(p).get(0).used` best-cluster gate exists solely in findsv
   (:705/:835, ported in find_sv where it belongs). No invariant ties cluster
   used-flags to clip-entry used-flags (filterSV VR:188/:193 and markSV mark
   clusters only; adjustMNP VR:355/:371 and realignlgdel/realignins mark clip
   entries only), so the gate dropped reachable records Java produces.
   Regression test `test_find_inv_sub_ignores_used_anchor_cluster`
   (used anchor + unused clip → still detected; mutation-checked to fail on
   the old code) pins the Java behavior.
3. **Success-gate off-by-one** — Rust computed `len = bp - softp + 1` first,
   then gated `len > 150 && len/|mlen| < 1.5`; Java gates on
   `bp-softp > 150 && (bp-softp)/|mlen| < 1.5`
   (StructuralVariantsProcessor.java:627-628), computing len after.
   **RESOLVED (W3.3c): real semantic divergence — FIXED to Java semantics**
   (gate on `bp - softp` before computing len). Divergence classes: gap==150
   accepted by old Rust, rejected by Java; ratio boundary gap=299/mlen=200
   rejected by old Rust ((gap+1)/200=1.5), accepted by Java (299/200≈1.495).
   Regression test `test_find_inv_sub_success_gate_uses_java_gap_semantics`
   covers gap∈{150,151} × mlen boundaries (mutation-checked to fail on the
   old gate). Oracle spot-checks: Colo829.chr20.INVSV MATCH vs live
   vardict-oracle VCF (20 variants, -f 0.001); NA24385 chr20 20156208 slice
   Rust VCF byte-identical pre/post fix (diff vs e79c54e build). Full corpus
   unchanged at 91 PASS / 1 FAIL.

### Affected Tests

- `Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;20156208;20234556;-f0.001`

### Remediation

Superseded by **RESOLVED (W1.4)** and **RESOLVED (W3.3)** above; all three latent items (W3.3a extra mark_sv sweep removed, W3.3b extra is_best_sv_used gate removed, W3.3c gap semantics fixed) are fully resolved and guarded by regression tests `test_find_inv_sub_success_marks_only_java_equivalents`, `test_find_inv_sub_ignores_used_anchor_cluster`, and `test_find_inv_sub_success_gate_uses_java_gap_semantics`. Full corpus remains 91 PASS / 1 FAIL.

---

## Category 9: Amplicon (MRE11a FIXED; Colo829 needs -D debug output)

### MRE11a — FIXED (9088d1a)

The +1 coverage (DP 847 vs 846) was the missing `parseCigarWithAmpCase` read filter:
Rust counted reads Java skips (edges > the -a distance from the amplicon or overlap
fraction below the -a threshold). Implementing the filter fixed MRE11a (MATCH).

### Colo829-18_S3 — FIXED (W3.1: `-D` debug output implemented)

With the read filter in place, Rust produces all 5 Colo829 amplicon variants correctly
(934099 G/C vd4, 934106 A/G vd3, 934123 T/+20 vd2, 934144 C/G vd6, 934144 C/+38 vd3 —
matching Java exactly). The `-D` debug column is now implemented and the case MATCHes:

- CLI: `-D`/`--debug` → `config.debug` (previously rejected as unsupported).
- `variant.DEBUG` accumulation: per-position `debug_lines` in `tovars::build_aligned_variants`
  — `create_variant`/`create_insertion` append one `Variant::debug_variants_content()`
  string each (insertions prefixed `I`; only when `config.debug`), then every finalized
  variant at the position gets `debug = join(" & ", debug_lines)` via
  `construct_debug_lines` (ports ToVarsBuilder.java:1002-1014; call sites :947/:996).
  Content format ports Variant.java:268-282
  (`n:cnt:F-fwd:R-rev:freq:bias:meanPos:pstd:meanQual:qstd:hifreq:mapq:qratio`);
  note mapq is printed with raw Java `Double.toString` semantics ("60.0"), the other
  floats via DecimalFormat-equivalent fixed precision.
- Amplicon output appends `"\t" + variant.debug + goodOutput + badOutput`
  (ports AmpliconOutputVariant.java:113-138). Good/bad entries:
  `\tGood<i> `/`\tBad<i> ` + 17 space-joined fields (debugAmpVariant, :258-278) + region.
  good/bad lists are tracked in `postprocess::amplicon` only when debug is on.
- Simple/Somatic printers also append `\t + variant.DEBUG` under debug
  (SimpleOutputVariant.java:116-117, SomaticOutputVariant.java:193-194); ported at the
  lib.rs emit sites. Coverage note: only the amplicon `-D` column has Tier-1 fixture
  coverage (Colo829-18_S3); the simple/somatic DEBUG appends have no corpus fixture.

**Known deviation (crash-family fix):** when a bad amplicon has no data at the position,
Java adds a *null* variant to badVariants (AmpliconPostProcessModule.java:178) which NPEs
in debugAmpVariant under `-D`. Rust renders an all-zero record
(`0 0 0 0 0 0 0 0 0 0.0 0 0.0 0 0.0 0.000 0 0`) instead of crashing.

**Full -D debug format spec (for implementation):** output line = 38 amplicon columns
+ `"\t"` + `DEBUG`, where `DEBUG = variant.DEBUG + goodOutput + badOutput`
(AmpliconOutputVariant.java:113-138).
- `variant.DEBUG` (ToVarsBuilder.constructDebugLines:1002-1014, only when config.debug):
  per-variant `debugVariantsContent` strings joined with `" & "`. `debugVariantsContent`
  (Variant.java:268-282) = `n:cnt:F-<fwd>:R-<rev>:<freq 0.0000>:<bias>:<meanPos 0.0>:<pstd>:<meanQual 0.0>:<qstd>:<hiFreq 0.0000>:<meanMapq>:<qratio 0.000>`,
  prefixed with `I` for insertions (debugVariantsContentInsertion).
- `goodOutput` = for each good variant index gvi: `"\tGood" + gvi + " " + join(" ", debugAmpVariant(" ", variant), region)`.
  `badOutput` likewise with `"Bad" + bvi`.
- `debugAmpVariant(delm, variant)` (AmpliconOutputVariant.java) = delm-joined:
  `totalPosCov posCov refFwd refRev varFwd varRev genotype("0" if null) freq(0.0000, "0" if 0) bias meanPos(0.0) pstd meanQual(0.0) qstd meanMapq(0.0) qratio(0.000) hiFreq(0.0000, "0" if 0) extraFreq(0.0000, "0" if 0)`.
- `goodVariants`/`badVariants` are `List<(Variant, regionString)>` produced by the
  amplicon post-processing (the good/bad amplicon-bias classification).
Implementation needs: a `-D` CLI flag → config.debug; accumulate variant.DEBUG in
ToVarsBuilder; track good/bad variants + region strings in amplicon post-process;
format and append DEBUG in amplicon_38 when debug is set.

### Historical notes (superseded by the fix above)
- Pileup/halo region, dedup, and supplementary filter all matched Java; the off-by-1 was
  purely the missing amplicon read filter.

---

## Category 10: Normal-sample Low-Quality Variant Registration — RESOLVED

### Status

CLOSED. Both affected cases MATCH end-to-end against the oracle.
- **Root Cause & Fix (commit `f7af5b7`):** `determinate_type` cloned `variant_to_compare` and called `is_noise` on the clone, discarding Java's coverage-zeroing mutation. Taking `&mut` and mutating in place fixed `NA12889.StrSomSamSpec` and `threeIndels`.
- In somatic mode, Rust now correctly reflects Java's in-place noise zeroing, matching normal block coverage and variant calls exactly.

### Affected Tests

- `Somatic;hs38DH.fa;NA12889.StrSomSamSpec.bam-NA12878.StrSomSamSpec.bam;chr20;43246926-43280248;-f0.001`
  (positions 43246972, 43246990 — last 2 remaining diffs after DELETION fix)
- `Somatic;GRCh37.fa;threeIndels_zero_deletion_tumor.bam-threeIndels_zero_deletion_normal.bam;chr20;62354864-62374864;-f0.1-Q10-F0x700;`
  (position 62356639)

### Investigation Notes

- NOT a base-quality gate: a thorough Java-vs-Rust audit of every `conf.goodq`
  (default 22.5) gate in CigarParser found faithful Rust equivalents for all. The
  matching-part SNV is added regardless of base quality in BOTH implementations
  (quality only routes to hicnt vs locnt). varCov = hicnt + locnt, no gate.
- Therefore the difference is a WHOLE READ that Rust processes but Java drops (or
  vice-versa), consistent with the +1 total coverage. Prime suspects: CigarModifier
  / realignment producing a different CIGAR (most common source of "1 extra read +
  1 low-qual variant"), the NM-absent read-drop branch, or the both-ends soft-clip
  drop regex semantics.

### Remediation

Enable the position-specific verbose trace in parse.rs for the offending position
and dump the same read in Java (`-y`) to see whether Java drops the read upstream
and at which filter, rather than at a quality gate.

### Deep findings (HCC2218 normal @ chrX:7811737, ref C → A)

- Offending read identified: `C0001:109:H03F6ADXX:1:2108:9475:83528`, flag 0x93
  (paired/proper/reverse/read2), CIGAR `1S98M`, **NM=7**, MAPQ 60. Base 'A' at
  query pos 49 (1-based 50), base quality **7**.
- Rust CigarModifier converts `1S98M` → `99M` (absorbs the 1S) and counts the
  qual-7 'A' as a variant: normal varCov=1, totalCov=212, meanQual=7.0, nm=7.0,
  meanPos=50.0.
- Java oracle intermediate output at 7811737 (normal block): totalCov=**211**,
  varCov=**0**, refFwd=99, refRev=112, BUT genotype=**C/A**, meanPos=**50.0**,
  meanQual=**7.0**, meanMapq=60.0, nm=**7.0**.
- **Key fingerprint:** Java's normal 'A' Variant carries the quality/position/nm
  stats of the single read yet has positionCoverage=varsCountOnForward=varsCountOnReverse=0.
  In `addVariationForMatchingPart` (CigarParser.java:1693) `varsCount++` is
  unconditional and precedes `meanQuality += q`, so quality-with-zero-count cannot
  come from a single normal pileup pass — the count was incremented then removed/moved
  (ToVarsBuilder count-adjustment ~357-373, adjustMNP, or realignment) while the
  quality accumulators remained. Java totalCov 211 (vs Rust 212) shows Java does not
  let this read contribute to refCoverage at 7811737.
- Ruled out: base-quality gate (none exists for matching-part SNV registration in
  either impl); NM gate (NM=7 ≤ mismatch=10); both-ends soft-clip regex (only a 5' 1S);
  `trim`/trimBasesAfter (default 0); chimeric detection (`-y` logs nothing for this read).
- `-y` verbose does NOT emit per-read pileup detail for this read, so the exact
  Java step that zeroes the count needs targeted Java instrumentation of
  nonInsertionVariants["A"] at 7811737 across the pipeline stages.

---

## Category 11: combineAnalysis (first-sample) — RESOLVED (W2.1)

### Status

CLOSED at W2.1a (`c1d50c8`) + W2.1b (`406db7f`). HCC2218 chrX flips to MATCH;
TumorA.StrongLOH and Colo829.MissedStrongLoh do NOT regress; corpus 87/5 → 88/4.

### What it fixes

- `Somatic;hg19.fa;HCC2218_BL_realigins.bam-HCC2218_C_realigins.bam;chrX;7811670-7811868;-f0.001;`
  Rust emitted a spurious StrongSomatic `+GAACCACTGAGACAGGAGAGCGAGGTGGAA` (+30)
  insertion at 7811769 (VD=3, meanMapq=0). Java creates it too (realignlgins30)
  but combineAnalysis re-runs the combined pipeline, does NOT find the insertion
  (vref==null → FALSE), and skips it. Re-CONFIRMED by W1.3 via live `-y` oracle
  trace. FIXED by W2.1b.

### The real root cause of the old regression (supersedes the realignlgdel theory)

The hypothesized "combined-pipeline realignment fidelity gap" **does not exist**.
Instrumented runs + oracle `-y` traces prove Rust's combined pipeline reproduces
every traced variant, at exactly the keys Java looks up:

| Variant | Java iteration key (`Start Combine`) | Rust combined-pipeline key | Verdict both |
|---|---|---|---|
| TumorA `-12^AAAAAATA` | 11035374 | 11035374 (pc=2) | `Combine: 1: 2 comb: 2` → "" keep |
| Colo829 `-11` | 59849477 | 59849477 (pc=3) | `Combine: 1: 3 comb: 3` → "" keep |
| Colo829 `-22` (second-sample path) | 59849679 | 59849679 (pc=2) | `Combine: 1: 2 comb: 2` → "" keep |

The actual bug: `combine_analysis` looked up the combined result at
`variant.start_pos` — which differs from the iteration position key once
adjComplex shifts anchoring (TumorA 11035374→11035376) or for raw deletion
anchoring (Colo829 -11 start_pos 59849476 vs key 59849477). Every lookup
missed ⇒ FALSE ⇒ record dropped. Java passes the pileup ITERATION position
key down through accept() → callingForBothSamples → printVariationsFrom*Sample
→ combineAnalysis and calls `getVarMaybe(vars, position, varn,
descriptionString)` (SomaticPostProcessModule.java:408). W2.1a ports that:
lookup by iteration key + description map.

### What landed

- W2.1a: faithful `combineAnalysis` port (:385-476) as tri-state
  Skip("FALSE")/Keep("")/Germline — the old port conflated FALSE with "";
  second-sample call site now skips on FALSE (previously fell back to
  StrongLOH); gate regex fixed from `starts_with('-')` to MINUS_NUM_NUM
  (`-\d\d`, Patterns.java:18); GERMLINE output header fields sourced from the
  ADJUSTED v2var (:333), closing the W1.3 review minor; verbose `Start Combine`
  prints the iteration key like Java (:391).
- W2.1b: first-sample call site (:177-217): gate `!SNV && (len>10 ||
  -\d\d) && pos_coverage < minr+3 && !contains('<')`; FALSE→skip / Germline→
  print adjusted placeholder as variant2 (:215) / ""→keep StrongSomatic (:212).

### Known simplifications (documented, not observed in corpus)

- Java updates its `maxReadLength` field from each combined run's result
  (:200/:311); Rust keeps the single seeded value. Only observable if read
  lengths differ between the main and combined runs; corpus shows no case.
- A pipeline callback error is treated as FALSE-skip. In Java the exception
  escapes combineAnalysis to accept()'s per-position try/catch
  (SomaticPostProcessModule.java:89-91): printExceptionAndContinue prints it
  and processing CONTINUES from the next position, aborting only after >10
  continued exceptions (Utils.java:261-278; MAX_EXCEPTION_COUNT,
  Configuration.java:300). The port keeps the skip control flow but emits the
  error under --verbose (loud-but-continue, not abort).
- Java inserts a fresh empty Variant into `v2.varDescriptionStringToVariants`
  before calling combineAnalysis (:195-196); inert in the port because nothing
  re-reads it downstream.

---

## Category 12: ArgumentNA systematic forward-read over-counting — RESOLVED

### Status

CLOSED. Both `ArgumentNA` cases (`-p-fisher` and `-p`) now MATCH end-to-end against the oracle.
- **Root Cause & Fix (commit `bc93c0c`):** Rust previously aliased `-M` to `max_mismatches` (`-m`), so `-M20` erroneously raised the mismatch threshold instead of setting the minimum matching read length (`minmatch`). Java defines `-m` as mismatch filter (default 8) and `-M` as a separate minmatch read-length filter (default 0, `CigarParser.java:303`). Adding `Config.minmatch` and the separate `-M` filter fixed both `ArgumentNA` cases with exact coverage parity.

### Historical Notes (Superseded)

The notes below are preserved for provenance:
- The initial hypothesis that this was not a flag-parsing gap was disproven by tracing `CmdParser.java:303` vs `Configuration.java`. Once `-M` was isolated as `minmatch`, forward reference read counts matched Java exactly.

---

## Category 13: NA24385 SV pairs under-count + splits over-count

### Symptom (after findMM5/-M fixes)

`Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;...;-f0.001` SV info
field (splits-pairs-clusters) differences:
- pairs systematically -1 in Rust: 20156480 `3-18-2`→`3-17-2`, 20161771 `2-3-1`→`2-2-1`,
  20179301 `3-3-1`→`3-2-1`, 20179472 `6-3-1`→`6-2-1`, 20224383 `2-10-1`→`2-9-1`
- splits over-count + duplicate INV: 20202056 Java `3-0-0` (one record) vs Rust
  `6-2-1` emitted TWICE; 20229671 Java `2-2-1` vs Rust `4-2-1` twice.

### Status

**RESOLVED (W1.4 / 676d6c6):** both symptoms on this case are fixed and the case
matches end-to-end (both -f0.001 variants). The pairs -1 was the spurious own-MAPQ<15
filter (commit 676d6c6). The splits over-count + duplicate INV at 20201754/20229671
was Rust's `find_inv_sub` continuing past its first success where Java's findINVsub
returns (`StructuralVariantsProcessor.java:675`) — see Category 8 "RESOLVED (W1.4)".
The remaining FAILs on other NA24385 chr20 cases (38386265, 44885701) were extra
records from other paths; both were fixed by W2.2 (see Category 3) and now MATCH.

**Historical analysis:** `check_pairs` (sv.rs:716) and `is_overlap` (sv.rs:594)
are byte-identical to Java `checkPairs`/`isOverlap` (StructuralVariantsProcessor.java
1742 / 1711). The pairs -1 was upstream in **discordant-pair clustering during
pileup** — how mates are grouped into `Sclip` clusters and their `varsCount`
(`prepare_sv_structures_for_analysis` / `add_sv`), which sets the cluster count that
`check_pairs` later reports as the max. Two distinct issues remained: (1) one discordant
pair is clustered/filtered differently in Rust (consistent -1); (2) split-read
over-counting at some INV/DEL breakpoints doubles the splits value and emits a
duplicate SV record (Category 8). sv.pairs is accumulated across many detection paths
(sv.rs:1100, 1398, 1537, 1707, 2021, 2441, 2735, 2860, 3043, 3249, ...).

### Remediation

Instrument check_pairs / the discordant-pair Sclip accumulation for 20156480 in both
implementations (Java toolchain) to identify the single pair Java counts but Rust
does not; separately trace find_inv_sub for 20202056 to find why Rust doubles the
splits and emits the INV twice.

---

## Category 14: Tier-2 VCF-Level Differential Verification vs Live Oracle

Swept 9 representative slices of all four modes against live VarDictJava 1.8.3 at
the VCF level (`harness/compare.sh`, oracle env `vardict-oracle`).

Result: **1 MATCH** (Colo829 INVSV chr20, 20 variants), **8 MISMATCH**, zero of
them engine-level: on every simple-mode slice tested the `.var` intermediates are
set-identical to the oracle's, so all diffs sit in the VCF emission path. Flag
parity was proven first (the W3.3 `-z` suspect turned out to be an invocation bug
in the harness itself, fixed gate-strengthening in W4a).

- **P1a — top-allele grouping key differs.** var2vcf_valid.pl:41 buckets rows
  by raw `.var` start column (`$hash{$chr}->{$a[3]}`), walks buckets in ascending
  numeric position (:129), sorts each by column 14 (af) descending (:135) and
  emits only the top-AF row per bucket (`$ALL = 1` without -A, :137). Our port of
  this rule (`postprocess/simple.rs`) truncated per ToVarsBuilder position group
  instead, whose key differs from raw start when indels/SNVs coexist nearby →
  records lost AND extra records, both directions (e.g. DEL `chr20:20179027 CT>C`
  java-only; DEL `chr20:20159222 AT>A` rust-only). **RESOLVED (W4.5a):** selection
  now buckets by the raw start column across the whole postprocessed list, stable
  af-desc sort per bucket, top row only; ties keep `.var` stream order like Perl's
  default mergesort.
- **P1b — extra AF<0.01 hard skip.** `lib.rs` dropped variants with
  frequency <0.01 && var_type non-empty in VCF mode; var2vcf_valid.pl never drops
  on AF without -S, it only annotates `f$Freq` (:159, `$Freq` = -f value, default
  0.02, :21). Examples: INS A/+333 @20178200 (af=0.0052), INS G/+440 @47131469
  (af=0.0044). The same `$Freq` drives the GT 0/1-vs-0/0 threshold (:220).
  **RESOLVED (W4.5a):** skip removed; the f-tag and GT threshold use config.freq
  (= the -f value both harness sides pass), matching perl exactly.
- **P1c — SVLEN off-by-one on non-INV SVs (found during W4.5a gating).**
  vcf.rs computed SVLEN = end-start+1 for all `<...>` ALTs; perl uses end-start
  and increments ONLY for INV (var2vcf_valid.pl:206-208). Masked at W4 because
  set-diffs short-circuited value comparison. **RESOLVED (W4.5a).**
- **P2 — Somatic/paired VCF writer.** Replicates `var2vcf_paired.pl` formatting in
  `output/vcf.rs` (`print_paired_vcf_header`, `format_paired_vcf_record`,
  `finish_paired_records`), wired into `lib.rs` run_somatic behind
  `!emit_intermediate` with raw-start top-AF bucketing. All
  three somatic slices now MATCH at the record level (comparator compares
  records only; headers are excluded per "Accepted divergences"): Colo829
  chr1 (16 variants), HCC2218 chrX -f0.001 (14), Zero_in_DUPSV chr19
  -f0.001 -X3 (61); corpus unchanged 91/1.
- **P3 — amplicon VCF wiring missing.** Same shape as P2; additionally Java's
  amplicon `-p` emits REF records (ALT='.') we don't produce.
  **RESOLVED (W4.5d).** Java routing truth: there is NO amplicon-specific
  writer — vardict pipelines single-sample output (simple OR amplicon) through
  `teststrandbias.R | var2vcf_valid.pl` (VarDictJava/tests/integration/
  tests_integration.sh:242; teststrandbias.R accepts 34/36/**38**-col input,
  "38 for amplicon mode"), and var2vcf_valid.pl:140 destructures 40 columns
  with amplicon-vs-simple decided purely by column presence:
  `:145 my $isamp = 1 if defined($ampflag)`. Amplicon deltas ported into the
  shared writer (`output/vcf.rs`): AMPBIAS filter (:185, pushed last after
  Cluster), pSTD suppression (:161), GDAMP/TLAMP/NCAMP/AMPFLAG INFO instead of
  DUPRATE/SV-block (:197,:233), ref==alt → ALT='.' GT=0/0 (:216-218, the `-p`
  REF records — already in the W4.5a port). Wiring mirrors W4.5a/W4.5c:
  `postprocess/amplicon.rs` raw-start top-AF bucketing behind
  `!emit_intermediate`; `lib.rs run_amplicon` prints the valid-writer header
  once + formats per record. **New fix surfaced by this task (writer-level):**
  with `-p` the engine sets freq=-1 (CmdParser.java:134-136), but perl's `$Freq`
  comes from the pipeline's OWN -f, not java's config; our single binary now
  carries both (`Config.af_threshold` = literal -f for all writer FILTER/GT/
  header sites, engine `freq` untouched). Gates: s7 MRE11a -a 10:0.95 -p MATCH
  (96 variants); s6 Colo829-18_S3 -a 10:0.95 MATCH (4 variants); tiny.bed
  simple regression MATCH (4); corpus unchanged 91/1. NOTE: s6 WITH -D cannot be a Tier-2 gate at all — under -D
  java appends tab-separated DEBUG fields to every .var row
  (AmpliconOutputVariant.java:113-138), line 1 exceeds teststrandbias.R's
  38-col check and the oracle pipeline dies ("Incorrect input detected"),
  emitting a headers-only VCF; upstream debug tables are simply not consumable
  by the R script. Our engine-level -D parity is covered by the corpus Tier-1
  fixture (Colo829-18_S3 PASSes with -D).

Tier-2 Verification Evidence:
- tiny.bed P1a repro: MATCH (4 variants); was java-only DEL CT>C@20179027.
- NA24385 chr20:38386265 (-f 0.001): MATCH (312 variants) after P1c fix.
- NA24385 chr20:20156208 (-f 0.001): record sets equal except ONE AF-tie bucket
  (chr20:20202875 T>A vs TA>T, both af=0.0179): stable-but-different engine
  emission order at an AF tie inside one bucket; matching requires engine-order
  replication — deferred (engine work, out of writer scope). Both engines'
  UNSORTED `.var` streams were inspected;
  SIX live oracle runs (3 plain + `-th 1`/`-th 4` + JVM SerialGC/-Xmx512m) are
  byte-identical (md5 `4d01a068…`) and all emit `TA>T` before `T>A`, while our
  stream emits `T>A` first — no oracle self-flip was captured, so this is NOT
  demonstrated nondeterminism (earlier "proven oracle-side" wording withdrawn).
  At chr20:20178765 both engines emit the same order and outputs match.
- Record-set delta on chr20:20156208 pre→post fix: restored exactly
  {20158483 AG>A (P1a), 20179027 CT>C (P1a), 20178200 A>+333 (P1b)}; removed 21
  lower-AF alleles sharing raw starts with higher-AF winners (Perl suppression).
- Full corpus unchanged: 91 PASS / 1 FAIL (the documented hard_clip stale
  fixture).

**Harness `-f` duplication (W4.5a, gate-strengthening).** run_oracle.sh put its
base `-f $AF` (AF_THR, default 0.01) BEFORE caller flags; commons-cli keeps the
FIRST repeated option value (CmdParser.java:605-608 getParsedOptionValue →
values[0]; config.freq at :114), so a case flag like `-f 0.001` never reached the
oracle engine/perl while compare.sh's rust side (clap, last-wins) honored it →
silent config asymmetry. run_oracle.sh now prefers an explicit case `-f` for
java+perl alike.

**GRCh37.fa fixture drift (environment, not a port bug).** The Colo829 SV slices
(DUPDiscSV chr20:47.12M, INVSV chr20:46.84M; BAM headers use `chr20`) currently
resolve to an ALL-N `chr20` block in VarDictJava/testdata/fastas/GRCh37.fa, so
BOTH vardict-java and vardictrust correctly emit zero records there (vacuous
MATCH). 12bp reference
flanks exist nowhere in the current
1GB file (mtimes preserved at Jun 2, content is not). The NA24385 slices are
unaffected (their BAMs use unprefixed names → correctChromosome strips to contig
"20", which holds real sequence at those loci). Regenerate the fasta from a real
GRCh37 before using the Colo829 slices as Tier-2 signal again.

STATUS: **RESOLVED (closed W4.6, 2026-08-23).** All writer items P1a/P1b/P1c/P2/P3
carry Tier-2 MATCH evidence. Engine (`.var`) parity unaffected: corpus remains
91 PASS / 1 FAIL.

### Tier-2 Benchmark Sweep

| Slice | Result |
|---|---|
| s1a NA24385 chr20:20156208 `-f 0.001` | MATCH (1048 variants) — AF-tie resolved via JavaHashMap (T10, 2026-10-08) |
| s1b NA24385 chr20:38386265 `-f 0.001` | MATCH (312 variants) |
| tiny.bed P1a repro `-f 0.001` | MATCH (4 variants) |
| s3 Colo829 chr1 somatic (defaults) | MATCH (16 variants) |
| s4 HCC2218 chrX somatic `-f 0.001` | MATCH (14 variants) |
| s5 Zero_in_DUPSV chr19 `-f 0.001 -X3` | MATCH (61 variants) |
| s6 Colo829 amplicon `-a 10:0.95` (no -D, documented substitution) | MATCH (4 variants) |
| s7 MRE11a amplicon `-a 10:0.95 -p` | MATCH (96 variants) |
| s8 INVSV `-f 0.001` vs hg19.fa (non-vacuous post-regen) | MATCH (20 variants) |
| s9 DUPDiscSV `-f 0.001` vs hg19.fa (non-vacuous post-regen) | MATCH (71 variants) |

Corpus rerun unchanged: 91 PASS / 1 FAIL. Milestone tags
`v0.2.0-somatic` and `v0.3.0-full-parity` cut at this state.

**Residual — RESOLVED (closed T10, 2026-10-08 via `JavaHashMap`).** Formerly permanent
limitation (P3 root-cause, 2026-08-23). The s1a AF-tie emission order at
chr20:20202875 was a **Java `java.util.HashMap` iteration-order artifact across
position keys**. In T10, OpenJDK 8+ `java.util.HashMap` spread, bin, resize, and
iteration semantics were faithfully implemented in `vardict_core::javamap::JavaHashMap`
(validated by Java oracle op-log replay tests). Wire-in at output-reaching variant
maps (`non_insertion_variants`, `insertion_variants`, and `aligned_variants`) replicates
Java's cross-position emission sequence, achieving 100% variant parity (1048/1048 variants,
calling `('20','20202875','TA','T')` exactly matching Java).

### Technical Mechanism: Hash Iteration Order in Allele Tie-Breaking

1. Both engines register the deletion `TA>T` at pileup map key **20202876**
   (first deleted base; CigarParser.java:1748 `getVariation(nonInsertionVariants,
   start, "-<dellen>")`) and derive the printed `.var` start as key−1 via the
   anchor-base prepend (ToVarsBuilder.java:735-739: `refallele =
   ref[position-1]; startPosition--`). The SNV `T>A` sits at map key 20202875.
   Registration is IDENTICAL in both engines — this is not a registration bug.
2. Within one position key the engines agree exactly: Java `sortVariants`
   (ToVarsBuilder.java:298-310: meanQuality×positionCoverage desc, then
   descriptionString asc) is ported verbatim (`tovars.rs sort_variants`). The
   same-key insertion tie at chr20:20178765 orders identically in both engines.
3. The divergence was **cross-position emission order**. Java iterates
   `alignedVariants` — a plain `HashMap<Integer, Vars>` (ToVarsBuilder.java:93,
   AlignedVarsData.java:14) — by hash-slot order in SimplePostProcessModule.java:35.
   With Integer keys (hashCode = value) that order is deterministic for a given
   insertion/rehash history — hence six oracle captures are byte-identical
   (md5 `4d01a068…`) — but it is NOT position-sorted: around
   the locus the java stream runs key sequence …20202802 → 20202876(del) →
   20202875(snv) → 20202850 → 20202978 → 20202742… (java-unsorted-run1.var
   lines 680-692), i.e. non-monotonic. By using `JavaHashMap` in VarDictRust,
   the identical sequence is produced, resolving the AF tie in exact agreement with Java.

Impact: record SETS are equal everywhere; the only observable effect is which of
two exact-AF-tie alleles sharing a raw-start bucket survives var2vcf_valid.pl's
first-in-stream rule when the tied rows originate from different map keys — at
most a one-record VCF difference per such locus (known instance: this one).

Upstream issue draft: consolidated as UB-3 in [docs/upstream-bugs.md](docs/upstream-bugs.md).

Somatic Differential Parity Evidence:
- Colo829-18/19 chr1:755917-756517 (defaults): MATCH (16 variants).
- HCC2218 BL/C chrX:7811670-7811868 (-f 0.001): MATCH (14 variants).
- Zero_in_DUPSV chr19:290742-291636 (-f 0.001 -X3): MATCH (61 variants).
- Full corpus rerun unchanged: 91 PASS / 1 FAIL.

### Reference Sequence Provenance & Verification

`harness/convert_fastas.py` was re-run over all bundled fastas; an in-session
md5 comparison observed every file regenerated **byte-identical** to its
pre-regen state — however the md5 artifacts were not retained (pre-state
unrecoverable, so no post-hoc re-verification is possible). The determinism
conclusion remains behaviorally corroborated: post-regen corpus 91 PASS / 1 FAIL
on the hg19/hg1k cases, GRCh37 slices still MATCH, `regen.log` rc=0.
The investigation also established the GRCh37.fa drift's true scope: its CSV
chunk source contains **no chr20 sequence near the Colo829 SV slices at all**
(its only `chr20`-named chunk sits at 60464465), so a real-sequence GRCh37.fa
at chr20:46.8M/47.1M was never generator-reproducible — the W4-era file was
not CSV-derived, and no fasta was hand-edited.

Resolution: the INVSV/DUPDiscSV Tier-2 slices run against **hg19.fa**, which is
what their upstream integration cases declare
(`Simple,hg19.fa,Colo829.chr20.INVSV.bam,…`) and what the corpus uses; hg19.fa's
CSV contains both regions under contig `chr20`. Recovery evidence
(empirical verification): oracle DUPDiscSV vs GRCh37.fa = 0 records
(headers-only VCF) → vs hg19.fa = **71 records**.

---

## Known non-bugs

| Item | Explanation |
|------|-------------|
| Integration case `Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1` fails | STALE FIXTURE ONLY: the fixture encodes a FIRST-chunk-wins reconstruction of the reference CSV, whose two overlapping `test` chunks conflict at 1734/2290 positions; today's converter (and test1's fixture) use last-write-wins. Live oracle under the first-wins fasta reproduces the fixture exactly — proven, not assumed. The `-p` null-record quirks are replicated in `postprocess/simple.rs` and `output/table.rs`; output is byte-identical to live oracle (7/7 rows sorted). Residual fixture diff is staleness-only. This fixture flaw is permanent until regenerated against a deterministic converter. |

| `addSV` stores `mend = mateStart_ms + rlen`, NOT the passed `mateEnd_me` (CigarParser.java:2325-2327) | Replicated verbatim in `add_sv` (cigar/parse.rs:2824-2826). Java's stored mend silently depends on every call site satisfying the identity `me == ms + rlen`; today they all do (`me` is built as `mate_start + total_length_including_soft_clipped`, parse.rs:2090), so outputs agree. If a future call site ever passes `me ≠ ms + rlen`, keep storing `ms + rlen` — do not "fix" toward `me` (§5 bug policy). |

### Upstream issue drafts — consolidated in docs/upstream-bugs.md

The v1.8.3 output bugs faithfully replicated per §5 (UB-1 spurious pos-0 record on `-p`,
UB-2 null-record cell formatting, UB-3 OpenJDK HashMap iteration order via `JavaHashMap`,
UB-5 SimpleMode SV split-read filter, UB-6 markDUPSV window truncation, UB-7 single-winner
multi-allelic collapsing), and the upstream integration-fixture inconsistency (UB-4)
are consolidated with filing-ready draft issue texts in
[docs/upstream-bugs.md](docs/upstream-bugs.md) (version pin 1.8.3 @ 0886261).

---

## Accepted Divergences (by design)

| Item | Reason |
|------|--------|
| VCF header `##source=VarDict_v1.8.2` vs our source string | Provenance honesty (rewrites.bio) |
| Oracle script splits `-f0.001` → `-f 0.001` | VarDictJava CLI parser limitation |
| Oracle script rewrites bare `-z` → `-z 1` (W4a) | Java's `-z` has an optional numeric arg (commons-cli); a bare trailing `-z` swallowed the BED path and crashed vardict-java silently. Explicit `-z 1` == upstream DEFAULT_ARGS semantics; use attached `-z0` to disable zero-based |
| compare.sh passes `-z` to both sides by default (W4a) | Mirrors corpus and upstream IntegrationTest DEFAULT_ARGS (`IntegrationTest.java:26`); amplicon mode lacks the 4-column-BED auto-default, so omitting it is a real flag divergence |
| Rust VCF emits `##reference` + `##contig` from the `.fai`; oracle VCF headers lack them (W4.5c fix 1) | Rust implements the issue-#320 fix (`##contig` lengths from the `.fai`) and emits both headers unconditionally; the oracle invocation passes no `-G`/`-b`, so vardict-java/var2vcf_*.pl never emit them. `compare_vcf.py` skips `#`-lines, so Tier-2 MATCH is record-level and does not cover headers. Header parity would require passing `-G`/`-b` to `harness/run_oracle.sh` — a gate change, deferred. |
| `--fisher` single-sample runs may diverge byte-wise in raw-cell fields (NM/HIAF/etc.) | With `--fisher` run_oracle.sh skips R entirely (:55/:71) and var2vcf_*.pl consume Java's RAW 38-col cells (unstripped DecimalFormat forms, e.g. golden fixture `NM=1.0` ×1410, `HIAF=0.0020` ×2). Our writers emit the R-echoed (stripped) forms everywhere except DUPRATE/SBF, which are already pipeline-split. No Tier-2 gate covers a `--fisher` single-sample slice; closing this needs a pipeline-aware formatter family keyed on `config.fisher` plus its own gate slice (post-1.0 candidate). Paired+fisher shares the same shape. |
| Splice-mode stdout intron-count lines have no stable order (both engines, by design) | Java's `spliceCount` is a plain `HashMap<String,int[]>` (CigarParser.java:65) printed via unordered iteration (:102); we mirror it (`HashMap`, iterated at lib.rs:233). The line SET is identical run-to-run; the ORDER is not — product-free, matching the oracle's own nondeterminism (#282 family). Would flake if ever snapshotted: keep splice-mode stdout out of byte-comparators or sort before comparing. |

---

## Resolved

| Item | Fix | Commit |
|------|-----|--------|
| SBF 1-ulp difference | Replaced dnhyper(1.0) with direct statrs::distribution::Hypergeometric pmf summation to match R's dhyper precision | T11.final |
| Soft clip filter regex mismatch | Fixed to match Java `^\d\dS.*\d\dS$` | dd43d8a |
| Somatic post-process coverage detection | Removed extra v1_empty/v2_empty check | a9579bb |
| DISCPAIRQUAL threshold 15→35 | Match Java Configuration.DISCPAIRQUAL=35 (fixed NA24631 extra DUP) | 8d0f12c |
| MINSVCDIST threshold 5→1.5 | Match Java Configuration.MINSVCDIST=1.5 | 72fb0ff |
| SVMAXLEN 2000000→150000 | Match Java Configuration.SVMAXLEN=150000 | 85e11cd |
| INV splits/VD double-count | Added SOFTP2SV check (is_best_sv_used) to find_sv tmp5/tmp3 loops; removed extra tv adj_cnt block from find_inv_sub. Fixed NA24385 INV at 20158414 (splits 6→3, VD 6→3, DP 88→73) | 3c5a6ce |
| Somatic --fisher ignored (always 55 cols) | Added somatic_61 formatter (per-sample + total pvalue/oddsratio, trailing-zero-strip) routed on config.fisher. Fixed NA12878 part1 fisher test | c6269e5 |
| Somatic DELETION emitted non-zero tumor cols | Build a zeroed tumor Variant carrying only header fields (Java passes tumorVariant=null). Fixed Colo829 somatic + TumorA StrongLOH tests | fe5962d |
| adjComplex suffix-trim loop bound | Changed `len >= n` → `len > n` to match Java `len-n > 0` (latent fidelity fix, no output change) | 3f36e6e |
| is_noise side effect discarded | determinate_type cloned variant_to_compare and called is_noise on the clone, discarding Java's coverage-zeroing mutation. Take &mut and mutate in place. Cracked via instrumented-Java classpath override. Fixed StrSomSamSpec | f7af5b7 |
| Complex left_seq endVariant source | somatic_55/61 read left/right seq + shift3/msi/msint from the adj_complex'd tumor; Java reads them from endVariant (the normal variant2, not adj_complex'd) in both-samples/LOH branches. Added end_is_variant2 flag; fixes Complex left_seq shift (threeIndels 62361372) | b3e2c38 |
| findMM5 dead n2-loop (dropped seq.length()) | Java findMM5 n2-loop guard `-1-n-1-n2 >= 0` is never true (seq.length() dropped), so misp is never set; Rust "fixed" the bound, over-triggering realigndel misp5 removal of adjacent SNVs (T@62360386, hicov 47→45, +T insertion hi_freq). Replicated Java's dead loop. Fixed threeIndels | 8fda935 |
| P2 writer byte-exactness family (4 items) | QUAL: both VCF writers now multiply log2(vd) by the HALF_EVEN 1dp table-rounded qual (`var2vcf_valid.pl:190` / `var2vcf_paired.pl:228`; cells from SimpleOutputVariant.java:155/:206). ODDRATIO (simple): reciprocal of a (0,1) value stringified with Perl's default NV form (%.15g) and non-reciprocal strings passed through verbatim (`var2vcf_valid.pl:147-151,:236`) — live oracle forms `1.84525676747919` / `10.0806451612903`. HIAF (simple): Java's RAW DecimalFormat("0.0000") cell, zeros kept (`SimpleOutputVariant.java:127-129/:211`) — teststrandbias.R's recycled colClasses force input col 23 to character so R echoes it verbatim; paired writer stays stripped (testsomatic.R leaves hiaf numeric). DUPRATE: 1dp without `--fisher` (36-col cell, :223), 2dp with (:175). Gates: tiny.bed / Zero_in_DUPSV chr19 / MRE11a all MATCH with oracle-vs-rust byte diffs eliminated except the SBF 1-ulp residue above; corpus 91/1 | P2 |
| INV desc tail built from revcomp(consensus) | find_inv_sub ignored findMatchRev's `matchedSequence` and appended revcomp(full consensus) to the `-N^seq<invM>seq` desc; Java appends the matched suffix as-is for dir==-1 / revcomp+prepend for dir==1 (StructuralVariantsProcessor.java:625-632). The wrong tail broke realignIndels' findMM5-based absorption of adjacent variants. Fixed Colo829 INVSV | W1.2 |
| Sclip used-flag set on a clone | find_inv_sub marked its cloned Sclip used instead of the live soft_clips_3_end/5_end entry, so realigndel re-absorbed the same clip (+7 reads). Java mutates the map value itself (StructuralVariantsProcessor.java:663), which survives later softp/bp shifts. Fixed Colo829 LGDEL3 | W1.2 |
| partial_pipeline discarded overlapping soft clips | Rust merged partial-pipeline soft_clips_3_end/5_end with entry(pos).or_insert (fresh counts dropped when position existed); Java's partialPipeline shares the SAME live maps so overlapping pipelines accumulate counts (e.g. cnt=2 clusters that realignlgins needs). Sum all Variation/nt/seq counts on merge. Fixed Colo829 DUPDiscSV (+342 insertion, SV info) | W1.2 |
| merge_sclip dropped SV state on re-merge | The count-summing merge kept Variation/nt/seq but discarded mates, soft-position counts, disc and the start/end/mstart/mend min/maxes; Java's shared-entry accumulation (CigarParser.addSV/addisccnt, CigarParser.java:2295-2341) keeps all of them. merge_sclip now unions mates, sums soft/disc and applies the addSV min/max conditionals (zero source fields skipped so SV-less entries don't clobber) | W1.2 fix1 |
| -M flag conflated with -m (minmatch vs mismatch) | Java -m=mismatch filter (default 8), -M=separate minmatch read-length filter (default 0, CigarParser.java:303). Rust aliased -M to max_mismatches, so -M20 raised the mismatch threshold instead of setting minmatch. Added minmatch Config field + -M option + read-length filter. Fixed ArgumentNA x2 (Category 12) | bc93c0c |
| Spurious own-MAPQ<15 filter in discordant-pair collection | Java filters discordant pairs only on the MC soft-clip pattern and the mate's MAPQ (MQ tag < 15); it does NOT filter on the read's own MAPQ. Rust had an extra `if mapping_quality < 15 { return; }` dropping ~1 pair per SV cluster → systematic pairs -1 (e.g. 3-18-2 vs 3-17-2). Removed the filter (kept the mq binding passed to add_sv). Fixed the NA24385 pairs under-count | 676d6c6 |
| find_inv_sub missing Java early-return | Java findINVsub returns after its FIRST successful INV (StructuralVariantsProcessor.java:675 `return vref;`); Rust's loop kept scanning and built an extra INV from a later svfinv5 cluster (20201601) Java never reached — stamping pairs/clusters into the record at 20201754 (`3-2-1`, DP=92 via refCoverage[softp-1]) where Java's nested loop yields `3-0-0`, DP=109. Added `return Ok(())` after the success block + regression test. Fixed NA24385 chr20:20201754 INV (both -f0.001 variants) and Cat 13 splits over-count | W1.4 |
| Missing amplicon read filter (parseCigarWithAmpCase) | Rust never ported CigarParser.parseCigarWithAmpCase: in amplicon mode, reads whose edges are >the -a distance from the amplicon or whose overlap fraction is below the -a threshold must be skipped. Without it Rust counted extra reads (MRE11a totalCov 847 vs 846). Implemented the -a "distance:fraction" parse (fallback 10:0.95) and the soft-clip/no-soft-clip overlap+distance checks. Fixed MRE11a amplicon | 9088d1a |
| Somatic pipeline missing find_all_sv + adj_snv | Added both calls for BAM1 and BAM2 | (prior session) |
| Oracle script crashes on concatenated options | Split `-fX` → `-f X` in run_oracle.sh | (prior session) |
| HashMap non-determinism in tovars | Changed to BTreeMap | e6c344b |
| Missing realign_indels in somatic pipeline | Added for both BAMs | (prior session) |
| SV distance guard missing | Added 1Mbp limit in find_match | (prior session) |
| Tiled streaming memory bounding | Halo intervals and tile streaming scaling whole chr1 to 2.9GB peak RSS (4.85x speedup) | 6ee4368, a4be62e |
| amplicon_40 table format | Added 40-column amplicon intermediate serialization under --fisher (AmpliconOutputVariant.java:145) | 3076620 |
| mark_dup_sv pe coordinate alignment | Fixed 1bp mate-pair window shift (pe instead of pe + 1) matching StructuralVariantsProcessor.java:1468 | a2d6c81 |
| SimpleMode SV splitread filter | Emulate var2vcf_valid.pl:197-205 dropping `<DUP>`/`<INV>` when splitreads < 1 in simple mode | a2d6c81 |

* **PR #402 family** — Java per-read skip replicated; >10 such reads: Java exits, VarDictRust continues (robustness, AGENTS.md §5)
