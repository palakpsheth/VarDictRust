# VarDictRust — Strategic Implementation Plan

> A from-scratch, **idiomatic-Rust** rewrite of **VarDictJava v1.8.3** into a single, statically-compiled,
> memory-safe binary that is **architected for large-scale WGS somatic data**, built under the
> **rewrites.bio** manifesto and **Karpathy AI-agent** discipline, and validated by **two-tier differential
> testing** against the original tool. All specs are ground-truthed against the **local v1.8.3 clone**
> (16,399 LOC Java + the checked-out Perl/R submodule + 94 integration fixtures), with `file:line` citations.

---

## 1. Context — Why this project exists

**The problem.** VarDict is a highly sensitive somatic/germline variant caller (Lai et al., *NAR*, 2016), prized in clinical oncology for amplicon-bias awareness, local realignment from soft-clipped reads, and complex-variant detection. Production VarDictJava v1.8.3 runs as a **fragmented 3-language pipeline**: Java core (BAM → tab-delimited intermediate) → R (`teststrandbias.R`/`testsomatic.R`, Fisher) → Perl (`var2vcf_valid.pl`/`var2vcf_paired.pl`, filter + VCF). It has verified, *outstanding* failure modes (the project is now **unmaintained** — confirmed on issue #396): a recurring **`StringIndexOutOfBoundsException`** crash family in `CigarModifier` (#341 still open, +#243/#253/#283/#7/#38/… ), **`OutOfMemoryError: GC overhead`** on large regions (#197/#384, closed as "split your BED"), **CRAM unsupported** (#249/#380), plus throughput loss from triple text serialization and R's `textConnection`.

**The outcome we want.** One Rust binary that does everything the 3 stages do in a single in-memory pass, **scales to modern high-depth WGS**, never falls over, and produces **output equivalent to VarDictJava v1.8.3 on valid inputs** (the prime directive; §3). v1.0 is "VarDictJava that doesn't crash, reads CRAM, and scales" — a pure, fully-validatable win.

**What prompted it.** A strategic blueprint (`Strategic Blueprint…docx`) proposed this citing *rewrites.bio*. This plan operationalizes it **after independently verifying every claim against primary sources, the local v1.8.3 source, and all 289 GitHub issues** — correcting the places where the blueprint *and* our own first-pass reading were wrong (§7).

**Local assets (verified).** `./VarDictJava/` (v1.8.3, `0886261`, 60 files/16,399 LOC); `./VarDictJava/VarDict/` (Perl/R submodule, **checked out**, pinned `009e017`); `./VarDictJava/testdata/integrationtestcases/` (**94** cases incl. the `TwoDelsInsStringBoundError` crash case); `./VarDictJava/tests/integration/` (`expected/vardict.{simple,paired}.vcf` + `raw_input/raw.vardict.{simple,paired}.var`). Remote: `git@github.com:palakpsheth/VarDictRust.git` (empty).

---

## 2. Governing principles

**rewrites.bio (3 commitments):** Credit the original (README/`--version`/Lai 2016/AstraZeneca-NGS); **Emulate exactly** ("a faster tool that produces the same results" — output equivalence to v1.8.3 is the definition of done); Be transparent about AI (document tools, role, validation coverage). *Think big* (fold 3 stages into one in-memory pass) / *work small* (the loop: pick ONE function → port → compare → match=next/mismatch=debug-this-one). "AI agents are unreliable for correctness — validate every step."

**Karpathy (4 principles):** Think before coding; Simplicity first; Surgical changes; Goal-driven execution (success criterion = `compare.sh <fn>` exits 0). Encoded in the project `CLAUDE.md` (§10).

**Scale-first (standing requirement).** Modern somatic NGS has exploded in depth and breadth (deep WGS). **Optimize for large data by default, everywhere** (§6). Peak memory must be bounded *regardless of region size* (this is what retires the OOM issue class by design); throughput must scale ~linearly with cores. Scale is a permanent acceptance criterion, not a late optimization.

**Idiomatic Rust (not transliteration).** Even though we port, we write **idiomatic Rust** (§9), not Java-in-Rust. rewrites.bio agrees: rethink, don't transliterate.

---

## 3. Locked decisions

1. **Scope: full parity, equal depth.** All four modes (Simple/Somatic/Amplicon/Splicing+SV) documented and built to comparable depth; construction follows the dependency DAG (shared engine first).
2. **Dev/CI: Ubuntu 12c / 24 GB / amd64 + GitHub Actions.** Dev loop on sliced/downsampled fixtures; CI on standard GH runners. amd64 ⇒ SIMD = **AVX2/SSE**; Linux ⇒ `apt`/conda tooling and a clean **x86_64-unknown-linux-musl** static binary.
3. **I/O foundation: noodles (pure Rust).** Single static binary, no C toolchain. Pin **exact** versions (API "experimental"); golden tests guard upgrades. Escape hatch: route one hot path through `rust-htslib` only if benchmarking proves a real gap.
4. **Equivalence bar: structural-exact + float tolerance.** Exact: identical variant set, FILTER decisions, INFO/FORMAT keys, GT, TYPE, STATUS, and **column counts/order** (36/38 single, 38/40 amplicon, 55/61 somatic). Tolerant: float fields within **ε=1e-6** after replicating R's 5-dp **HALF_EVEN** rounding (Fisher MLE tol `ulp(1)^0.25`).
5. **Issue-remediation scope: v1.0 = robustness only.** The initial release fixes **only** crash/robustness/input/IO/memory/systems/determinism issues (all **FIX-BY-DESIGN**, changing no valid-input output → fully validatable). Output-changing corrections and new features are deferred to a post-1.0 roadmap behind opt-in flags (§8).
6. **Scale-first by default** (see §2/§6) — bounded memory at any depth, linear core scaling.

---

## 4. Verified ground truth (`file:line` cited)
Paths under `VarDictJava/src/main/java/com/astrazeneca/vardict/`.

### 4.1 Execution flow & threading
`Main.main` (`Main.java:12`) → `CmdParser` → `Configuration` → `VarDictLauncher.start`. Launcher resolves chr lengths, sample names, regions (`RegionBuilder`), **amplicon detection** when BED has 8 numeric cols 7&8 (`VarDictLauncher.readBedFile:133-178`), adaptor seeds, the `GlobalReadOnlyScope` singleton. **Mode selection** (`VarDictLauncher.java:54-67`): splicing→`SplicingMode`; else region/non-amplicon → `hasBam2()? SomaticMode : SimpleMode`; else `AmpliconMode`.

**Per-region pipeline** is a `CompletableFuture` chain on a **`DirectThreadExecutor`** (inline/synchronous): `SAMFileParser → CigarParser(false) → VariationRealigner → StructuralVariantsProcessor → ToVarsBuilder` (`AbstractMode.pipeline:44-56`). **5 stages run sequentially on one thread per region.** `partialPipeline` is re-entered by `realignlgdel` to pull reads from extended regions.

**Region parallelism** (`AbstractMode.java:112-143`): `newFixedThreadPool(threads)` + a **bounded `LinkedBlockingQueue<Future>` (capacity 10)**; producer submits one `Callable`/region, consumer `take()`s **in submission order** → **deterministic output**. Per-thread `ThreadLocal` `SamReader`/`IndexedFastaSequenceFile`. Shared `splice` = `ConcurrentHashSet`. Abort after `MAX_EXCEPTION_COUNT=10`; fatal → `System.exit(1)`.
→ **Rust:** `rayon` over regions/tiles collecting `(index, bytes)`, emit in index order (bounded backpressure). Per-region pipeline = **plain synchronous calls**. `thread_local!` readers.

**SomaticMode** (`processBothBamsInPipeline:103-131`): BAM1 pipeline to completion, then BAM2 **seeded with BAM1's `maxReadLength`**, then `SomaticPostProcessModule` compares; shared `splice`.

### 4.2 CLI (defaults shown)
`-b` BAM(s) (`tumor|normal`), `-G` ref (or `hg19/hg38/mm10`), `-R chr:start-end`, `-N` sample, `-c/-S/-E/-s/-e/-g` BED cols (1-based→0-based), `-f` **0.01**, `-r` **2**, `-q` **22.5**, `-Q` (unset), `-O` **0**, `-m` **8**, `-P` **5**, `-o` **1.5**, `-V` **0.05**, `-B` **2**, `-X` **2**, `-I` **50**, `-th` (bare→all cores else 1), `-k` **1** (realign), `-a` `"10:0.95"`, `-z` (BED 0-based default), `-3`, `-F` **"0x504"**, `-U/--nosv`, `-w/-W/-A/-L` SV (**300/100/4/1000**), `-Y` **1200**, `--fisher`, `-p`, `-i/--splice`, `--chimeric`, `--deldupvar`, `-t/--dedup`, `-Z`, `-DP OUT|ERR`, `-J/-j` CRISPR, `-mfreq` **0.25**/`-nmfreq` **0.1**. Constants: `SEED_1=17, SEED_2=12, ADSEED=6, LOWQUAL=10, bias=0.05, SVMAXLEN=150000, SVFLANK=50, DISCPAIRQUAL=35, EXTENSION=5000`.

### 4.3 Intermediate column schemas — CORRECTED (from the printers)
> R-gates (34/36/38; ≥48) are *input validators*; the Java printers emit the counts below — *these* are what the engine must reproduce. Formatter rule: no-`--fisher` uses `x==0?0:DecimalFormat` (keeps trailing zeros); `--fisher` uses `getRoundedValueToPrint` (strips trailing zeros).

- **Single (`SimpleOutputVariant.java`): 36 / 38.** `1 sample·2 gene·3 chr·4 start·5 end·6 ref·7 alt·8 DP·9 VD·10 RefFwd·11 RefRev·12 AltFwd·13 AltRev·14 GT·15 AF·16 BIAS·17 PMEAN·18 PSTD·19 QUAL·20 QSTD·21 MQ·22 SN·23 HIAF·24 ADJAF·25 SHIFT3·26 MSI·27 MSILEN·28 NM·29 HICNT·30 HICOV·31 LSEQ·32 RSEQ·33 seg·34 TYPE·35 DUPRATE·36 SV`. `--fisher`: PVALUE,ODDRATIO inserted at 21–22 (rest +2).
- **Amplicon (`AmpliconOutputVariant.java`): 38 / 40.** Cols 1–34 as Simple, then **35 GDAMP·36 TLAMP·37 NCAMP·38 AMPFLAG** (no DUPRATE/SV).
- **Somatic (`SomaticOutputVariant.java`): 55 / 61.** 7 shared + 18 tumor (8–25) + 18 normal (26–43) + 12 tail (44–55; **50=varLabel(STATUS)**, var1/var2 duprate+sv). `--fisher`: per-sample pval/OR inside each block (tumor 8–27, normal 28–47), **total pval/OR at 60–61**.

### 4.4 Fisher — port `FisherExact.java`+`UnirootZeroIn.java`
`(refFwd,refRev,altFwd,altRev)` → `m=refFwd+refRev, n=altFwd+altRev, k=refFwd+altFwd, x=refFwd`, support `[max(0,k-n)…min(k,m)]`; density via Hypergeometric `logProbability` (→ statrs). **Two-sided p** = Σ{ null-probs ≤ observed·(1+1e-7) } (R's method). **OR = conditional MLE** via Brent (`zeroinC`) solving `E[X|ncp]=x`, tol `ulp(1)^0.25≈1.22e-4` (not naïve ad/bc). Final p/OR rounded **5-dp HALF_EVEN** (`round_as_r`; JDK-8227248 vs-R caveat → ε covers it). Cross-check p-value vs `fishers_exact` crate (p-value only). Columns: strand bias **10-13**; somatic tumor **10-13**, normal **28-31**, T-vs-N `c(9, dp1-vd1, 27, dp2-vd2)` two-sided-by-min.

### 4.5 Somatic STATUS — in Java (`SomaticPostProcessModule.java`)
Col 50/54. `lofreq(-V)=0.05`. No-cov-sample1→**Deletion**; no-cov-sample2→**SampleSpecific**. Both present (`determinateType:346-369`): if sample2 var good — `vref.freq>(1-lofreq) && 0.2<v2.freq<0.8`→**LikelyLOH**; `v2.freq<lofreq||cov≤1`→**LikelySomatic**; else **Germline**. If not good — `<lofreq||cov≤1`→**LikelySomatic** else **AFDiff**. Override `isNoise && SNV`→**StrongSomatic**. Sample1-only→**StrongSomatic** (short indels run `combineAnalysis` re-pileup → maybe **Germline**/drop). Sample2-only→**StrongLOH**. Port `determinateType`/`isNoise`/`combineAnalysis` verbatim.

### 4.6 `isGoodVar` (`Variant.java:321-370`) & var2vcf
`isGoodVar` rejects unless: refallele non-empty; `freq≥0.01`; `hicnt≥2`; `meanPosition≥5`; `meanQuality≥22.5`; mapq/length sanity when ref hicnt>minr & freq<0.25; not deletion at spliced junction; `qratio≥1.5`; **freq>0.30⇒accept**; `meanMappingQuality≥-O`; MSI cuts; strand-bias `"2;1"`+`freq<0.20`. `var2vcf_valid.pl` defaults `-d3 -v2 -f0.02 -p8 -q22.5 -Q10 -F0.2 -o1.5`; FILTER IDs `PASS,d3,v2,f0.02,p8,pSTD,q22.5,Q10,SN1.5,NM5.25,Bias,MSI12,InGap,InIns,Cluster0bp,LongMSI,AMPBIAS`; INFO/FORMAT tags incl. `REFBIAS=rfwd:rrev,VARBIAS=vfwd:vrev,AD(R),RD,ALD`; paired adds `STATUS,SSF,SOR`. (Confirm against local pinned `VarDict/var2vcf_*.pl` at impl time.)

### 4.7 Core algorithms
- **`CigarModifier.modifyCigar`**: ordered fixpoint sequence (strip leading/trailing D; I→S; chimeric softclip removal via seed; `while(flag)` indel-cleanup folds; mis-softclip ±3-match/mismatch rules; leading-clip merge). **Unguarded `substring`/`charAt` crash sites** (`:91,:294,:311,:778`) → Rust bounds-checks all (no panic; output unchanged on valid inputs). Preserve order & group semantics.
- **`CigarParser`** (2376 LOC, the **pileup engine**): mismatch gate `NM-(I+D)>8→drop`; `direction=getReadNegativeStrandFlag()` drives strand counts (`Variation.incDir`); N/S/I/D/M dispatch; MNV growth; indel-adjacent folding; `processSoftClip` walks matched-clipped bases into alignment, stores remainder in `Sclip` (`TreeMap`). `tp=min(posExclSC+1, rlen-posExclSC)`.
- **`VariationRealigner`** (largest fidelity risk; `realigndel→realignins→realignlgdel→realignlgins30→realignlgins`): `findconseq` (accept iff `match/total>0.9`+heuristics), `ismatch`/`ismatchref` (≤3 mm AND <15%), `findbp`/`findbi`/`findMatch` (SEED_1, ≥3 distinct bases), `adjCnt` transfers counts + sets `extracnt`(→ADJAF). **Iteration order load-bearing** (`fillAndSortTmp`: count desc then position) — `adjCnt` mutates shared maps. Constants `0.9/0.8/0.15/0.12/SEED_1/SEED_2/longmm=3/indelsize=50` exact.
- **`ToVarsBuilder`** (variant **finalization**): Variation→Variant; AF/means HALF_EVEN; MSI/shift3 windowed scan; genotype; IUPAC reduction; `<DEL>/<DUP>/<INV>` spelling. Desc grammar: `+seq`/`-len`/`&seq`/`^seq`/`#…^…`. Replicate the hicov-excludes-insertions bug (`:534`) in compat mode.

### 4.8 Java→Rust module map (REFINED)

| VarDictJava | Responsibility | Rust target |
|---|---|---|
| `Main`,`VarDictLauncher`,`CmdParser`,`Configuration`,`GlobalReadOnlyScope` | entry/CLI/config | `src/main.rs`,`cli.rs`,`config.rs`,`global.rs`(`Arc<Config>`) |
| `RegionBuilder`,`data/Region` | BED/region + 8-col amplicon (`AMP_BED_ROW_FORMAT=(0,1,2,6,7,3)`) | `regions.rs` (noodles-bed) |
| `data/Reference`,`ReferenceResource` | indexed FASTA + SEED k-mer map | `reference.rs` (noodles-fasta) |
| `data/SamView`,`SAMFileParser`,`RecordPreprocessor` | BAM/CRAM read, query, dedup/`-F`/downsample | `bam.rs` (noodles-bam/cram) |
| `modules/CigarModifier`,`data/ModifiedCigar` | read-level CIGAR rewrite | `cigar/modify.rs` |
| `modules/CigarParser` | **pileup engine** | `cigar/parse.rs` |
| `modules/VariationRealigner` | realignment | `realign.rs` |
| `modules/StructuralVariantsProcessor`,`data/SVStructures` | SV stage + shared `findMatch`/`markSV` | `sv.rs` |
| `modules/ToVarsBuilder` | **variant finalization** | `tovars.rs` |
| `variations/*`,`collection/VariationMap` | variant state | `variations/` (`IndexMap`+`BTreeMap`) |
| `data/fishertest/*` | Fisher | `stats/fisher.rs` (+statrs) |
| `modes/*` | orchestration + runner | `modes/*.rs`,`runner.rs` |
| `postprocessmodules/*` | filter + STATUS + amplicon-bias | `postprocess/{simple,somatic,amplicon}.rs` |
| `printers/*` | intermediate table + VCF | `output/{table.rs,vcf.rs}` (noodles-vcf) |
| `data/scopedata/*` | per-stage payload | `scope.rs` |
| `Utils`,`Patterns` | safe substr/charAt, HALF_EVEN, regexes | `util.rs`,`patterns.rs` |

---

## 5. Target architecture
**One statically-linked binary**, single in-memory pass, deterministic-ordered parallelism (see §6 for scale). No intermediate text; the 36/38/55/61 schema lives as typed structs, serialized to text only under `--emit-intermediate` (dev/test, tier-1 diffing). Per-worker independent readers ⇒ memory scales with *local depth*, not genome size.

---

## 6. Scale-first design (large WGS somatic, by default)
This is how decision #6 is realized concretely, and how the OOM issue class (#197/#384/#373) is fixed *by design*:
- **Bounded peak memory at any depth.** The unit of work is a region **tile** sized to cap RAM; never hold whole-chromosome state. Large regions (e.g. the chr1 mega-region in #384) are split into **overlapping windows with a halo** so realignment context is preserved — and boundary variants are de-duplicated so **tiled output is identical to a single-region run** (scale must not change output; relates to region-boundary issue #369).
- **Streaming, ownership-driven memory.** Reads flow through the per-tile pipeline and are dropped promptly (RAII), eliminating the JVM's billions-of-temp-objects GC collapse.
- **Work-stealing parallelism** (rayon) across tiles for load balance between deep amplicons and sparse WGS; ordered collect for deterministic output with bounded backpressure (the Java capacity-10 pattern).
- **Throughput**: multi-threaded BGZF decode/encode; **AVX2** alignment hot loops (block-aligner); compact, cache-friendly, data-oriented structures (byte slices, reused buffers, `SmallVec`/arena in hot paths); sort only within a tile (no global sort).
- **Predictability**: document memory as a function of (tile size × depth × threads); optionally expose `--tile-size`. Scale acceptance criteria in §15: peak RSS bounded regardless of `-R` size; ~linear scaling to 12 cores; deep-WGS slice benchmarks gated in CI (CodSpeed + peak-RSS tracking).

---

## 7. Corrections (AI "confidently wrong" — incl. our own first draft)
Single-sample is **36/38** (not 34/36/38 — those are R input validators); paired is **55/61** (not 48); somatic **STATUS is computed in Java** (not Perl); `fishers_exact` gives **p-value only** (port the cond-MLE OR); `genomicframe-core`/`rustalign-simd` **don't exist**, `ksw2rs` is misnamed/stale, `seqair`/`cyanea-align` are unproven v0.1 → use noodles + block-aligner; **Snyk lacks Cargo support** → `cargo-audit`+`cargo-deny`; GIAB HG008 is **PDAC** (not PNET); blueprint benchmark table is "synthesized" → measure our own; threading is `DirectThreadExecutor`(sync)+bounded-ordered-queue (not CompletableFuture parallelism); dev env is **Ubuntu 12c/24GB/amd64 (AVX2)**. Kept-correct: noodles, rust-htslib, rust-bio, block-aligner, rayon, statrs, clap, anyhow, thiserror, criterion/divan/insta/proptest, cargo-dist(`dist`)/release-plz/CodSpeed/cargo-audit/deny/crates-io-auth-action, truvari/rtg-vcfeval/nf-test/som.py, SEQC2 HCC1395 (+WES), GIAB HG008.

---

## 8. Issue-remediation strategy (289 issues triaged: 64 open / 225 closed)

v1.0 = pure parity (A1); see ROADMAP.md.

**v1.0 — FIX-BY-DESIGN only (on by default; changes no valid-input output; fully validatable):**
- **CIGAR-crash family** → bounds-checked slicing in `cigar/`: #341 (open, still unfixed), #243, #253, #283, #7, #6, #14, #22, #34, #37, #38, #48, #135, #189, #192, #193. One correct impl retires the whole class.
- **NPE / AIOOBE / divide-by-zero** → `Option`/ownership/guarded math: #100, #96, #1, #153, #166, #272, #281, #274, #89, #375.
- **NumberFormat / locale** → typed numerics + locale-independent formatting: #195, #167, #175, #206, #338, #384(locale).
- **OOM / memory / scale** → §6 tiled streaming: #197, #384, #373, #222, #80, #12, #64, #32, #216.
- **CRAM / IO robustness** → noodles + `-G`-driven decode, synthesize NM/MD, lenient-validate-and-warn, fail-fast on missing index, clear bgzf errors: #249, #380, #364, #379, #187, #13, #28, #61.
- **Nondeterminism** → stable sort + full tie-break keys: #282(ordering), #181, #390, #356.
- **VCF metadata correctness (not variant output)** → `##contig` lengths from the FASTA `.fai`: #320.
- **Packaging / build / CLI / systems** → single static musl binary removes JVM/Gradle/classpath/Perl-R packaging pain: #401, #385, #24, #233, #205, #97, #182, #111 (and #396 motivates the whole rewrite).

**Post-1.0 — FIX-OUTPUT (opt-in behind flags / an `--improved` mode; each validated individually so v1.0 byte-compat is never lost):** depth accounting #92/#386/#323/#316 (VD>DP, AF>1, AD≠DP); over-filtering #62 (emit + FILTER-tag instead of silent drop); realignment drops/mis-calls #363/#306/#55/#198/#194/#26; SHIFT3 coords #292 (`ToVarsBuilder.java:723-725`); non-spec `<dup-N>` ALT #342; SN definition #381; soft-clip depth #282; `--nosv` honoring #264 (verify).

**Later — features:** UMI #374, long-read/ONT #368, read-name output #391, `-L` for `<DUP>` #328, haploid #353, progress log #324, PMEAN #263.

**Explicitly excluded (not-legit/stale/user-error/off-topic):** #336 (off-topic), #331/#396 (meta/infra), #382/#398/#225 (mismatched-contig config), #43 (self-inflicted), #27 (external tool), and the large set of usage/question threads where no defect was isolated.

---

## 9. Rust idioms & best practices (port ≠ transliteration)
- **Coordinate safety as types.** Newtypes `Pos0(u32)` (0-based half-open, from BAM) and `Pos1(u32)` (1-based, VCF/VarDict internal) with explicit `to_1based()`; the #1 silent footgun becomes a **compile error**. Likewise `ReadOffset`, `RefPos`, `Phred(u8)`.
- **Sequences as `&[u8]`/`Vec<u8>`**, never `String` — no UTF-8 overhead; matches the blueprint's "safe byte-slice manipulation". Replace Java's negative-index `charAt`/`substr` with bounds-checked slice helpers returning `Option` (fixes the crash family idiomatically).
- **Model the data contract as types, not magic.** `enum VariantDesc { Snv(u8), Mnv(Vec<u8>), Ins(Vec<u8>), Del(u32), DelIns{del:u32, ins:Vec<u8>}, Complex{…} }` with a `Display` that reproduces the exact `+seq`/`-len`/`&`/`^`/`#` tokens. `enum CigarOp`, `enum VarType`, `enum SomaticStatus` with `Display` matching VarDict strings.
- **Errors:** `thiserror` enums in the library, `anyhow` in the binary; **no `unwrap`/`expect`/`panic!`** on paths reachable from input; `?` propagation; actionable messages (this *is* the fix for the cryptic-stacktrace issues). `#![forbid(unsafe_code)]` in core; `unsafe` only via vetted SIMD crates behind safe APIs.
- **Concurrency:** `Arc<Config>` passed explicitly (no mutable global singleton); `Send+Sync` data; rayon work-stealing; `thread_local!` readers; **deterministic output** via stable sort + ordered collect (retires the nondeterminism issues).
- **Numeric fidelity:** one `round_half_even` util; explicit per-mode formatter rules; be deliberate about f64 accumulation order where it affects the 5-dp result.
- **Data-oriented for scale:** compact structs, reused buffers, `SmallVec`/arena in hot loops, struct-of-arrays for pileup columns where it helps; stream, don't materialize.
- **Decompose, don't god-object:** split the 2376-line `CigarParser` into cohesive functions; prefer iterators/pattern-matching over index loops + `instanceof`/string-switch — **while preserving iteration order** (IndexMap/BTreeMap).
- **Traits only where ≥2 real impls** (e.g. a `Mode`/`PostProcess` trait, a `VariantSink`); otherwise keep it concrete (Karpathy simplicity).
- **Provenance comments:** each ported fn cites its source (`// ports CigarModifier.java:132-228`). `#[must_use]`, doc comments, clippy-pedantic where reasonable.
- **Testing:** colocated `#[cfg(test)]` units; `insta` snapshots for `.var`/VCF; `proptest` for BED/CIGAR parsers; criterion/divan benches; golden differential tests.

---

## 10. Phase 0 — Scaffolding & AI-agent operating model
**First commands (on approval):** verify submodule (already checked out); Ubuntu setup `apt-get install samtools openjdk-17-jdk` (+ conda `vardict-java=1.8.3` option) and `rustup` (clippy/rustfmt + `x86_64-unknown-linux-musl`).
**Cargo workspace:** `crates/vardict-core` (lib) + `crates/vardict-cli` (bin); `harness/` (compare.sh, compare_var.py, compare_vcf.py, slice.sh, run_oracle.sh); `golden/` (+PROVENANCE); `benchmarks/`; `.github/workflows/`; `rust-toolchain.toml`.
**`CLAUDE.md`+`AGENTS.md`:** Karpathy 4 + meta-rule (evidence before completion) + emulate-exactly + **coordinate newtypes** + the port loop + definition-of-done (`fmt --check`+`clippy -D warnings`+`cargo test`+`compare.sh` MATCH+golden updated+PR states version/dataset/cmd) + unsupported-flags-error-loudly + do-not-use crate list (§7) + fidelity hazards (iteration order, HALF_EVEN, Phred+33) + scale rules.
**Subagents/slash-commands:** `/port-function <Class.method>`, `/compare-outputs [case|--all]`, `/check-coordinates`, `/snapshot-update`; reuse the user's existing `oh-my-claudecode`/agentic-skills conventions (parent dir); stop-gate hook blocks evidence-free "done".
**Reference oracle (pinned):** local gradle (JDK 8/11/17) or bioconda or Docker; `run_oracle.sh` records exact version into golden PROVENANCE.
**Seed corpus (zero new data):** convert the 94 `testdata/integrationtestcases/` + `tests/integration/{expected/*.vcf, raw_input/*.var}` into `golden/`.

---

## 11. Two-tier differential testing
- **Tier 1 — engine → `.var`** (`--emit-intermediate`): vs `raw_input/raw.vardict.{simple,paired}.var` + oracle tables; `compare_var.py` (exact ints/strings; ε=1e-6 floats; HALF_EVEN; the §4.3 formatter rules). Localizes mismatches to the engine, before R/Perl.
- **Tier 2 — full → VCF**: vs `expected/vardict.{simple,paired}.vcf`; `compare_vcf.py` (exact chrom/pos/ref/alt/FILTER/GT/TYPE/STATUS/INFO-keys; ε=1e-6 floats).
- **Slicing (24 GB):** `samtools view -b … chr:start-end` around known variants + downsampling; full files occasional/manual.
- **Datasets:** SEQC2 HCC1395 WES (primary somatic), GIAB HG008 (PDAC). **Secondary** truth-set tier: `som.py`/`rtg vcfeval` (small variants), `truvari bench` (SV) — confirms we didn't faithfully reproduce a *wrong* result; primary bar stays equivalence to v1.8.3.

---

## 12. Implementation phases (Red→Green[=compare MATCH]→Refactor; all four modes at equal depth)
1. **I/O + CLI + scale skeleton:** clap (full §4.2), `RegionBuilder` incl 8-col amplicon + **tiling/halo (§6)**, indexed FASTA+seed map, indexed BAM+filters, rayon ordered collect.
2. **CIGAR + pileup:** `cigar/modify.rs` (ordered fixpoint, **bounds-checked — fixes crash family**) + `cigar/parse.rs` (strand counts, softclip→Sclip, MNV/indel folding).
3. **Realignment:** `realign.rs` + `sv.rs` shared helpers; `findconseq`/`ismatch`(15%)/`findbp`/`findMatch`/`adjCnt`; preserve sort order+constants. Red: EGFR ex19 15-bp del recovered with matching VD/AF.
4. **Finalization + schema:** `tovars.rs` + `isGoodVar`; emit exact **36-col** `.var` (HALF_EVEN, formatter rules). Red: Tier-1 byte-match vs `raw.vardict.simple.var`.
5. **Statistics:** port `FisherExact`+`UnirootZeroIn` (two-sided p + cond-MLE OR, `ulp^0.25`, 5-dp HALF_EVEN); cross-check vs `fishers_exact`; `--fisher` 38-col.
6. **VCF + var2vcf_valid (single-sample complete):** all FILTER IDs + INFO/FORMAT via noodles-vcf (BGZF; **`##contig` from `.fai` — fixes #320**); confirm vs local pinned Perl. **M1.**
7. **Somatic/paired:** two-BAM flow (BAM2 seeded w/ BAM1 maxReadLength) + `SomaticPostProcessModule` STATUS (§4.5) + **55/61** schema + `var2vcf_paired.pl`. Red: MATCH on SEQC2 WES slices + `expected/vardict.paired.vcf`. **M2.**
8. **Amplicon:** 8-col primer math, `-a 10:0.95`, masking, AMPBIAS, **38/40** tail, multiallelic iterate-all (#68). Red: amplicon fixtures MATCH.
9. **SV / Splicing:** `StructuralVariantsProcessor`+`SVStructures`, `SplicingMode`; `SVTYPE/SVLEN/SPLITREAD/SPANPAIR/DUPRATE`. Red: SV fixtures MATCH; truvari concordant. **M3: full functional parity.**
10. **Performance/scale hardening:** AVX2 realignment (block-aligner), tile-size tuning, arenas, BGZF threads; `compare.sh --all` green after each; CodSpeed + peak-RSS gates on deep-WGS slices.
11. **Robustness fixtures + CI/CD/security/release:** crash-corpus tests (the `TwoDelsInsStringBoundError`/hard-clip/CRAM/missing-index cases — **fixes #341 etc.**); `ci.yml`(fmt+clippy-Dwarnings+test+compare --all); `differential.yml`; `bench.yml`(CodSpeed); `security.yml`(cargo-audit+deny); `release.yml`(`dist` x86_64-musl + OIDC trusted publishing + release-plz). MIT.
12. **Validation tier & docs:** full SEQC2/GIAB (manual) w/ som.py/rtg/truvari; README provenance + **AI disclosure** + validated-version/dataset/commands; CONTRIBUTING/CHANGELOG, supported-flag matrix, **fixed-issues list (§8)**.

---

## 13. Risk register
| Risk | L | Mitigation |
|---|---|---|
| `VariationRealigner` fidelity (order/mutation) | High | Match comparators+`.used`+constants; Tier-1 `.var` diff isolates; port one realign-fn at a time. |
| Iteration-order drift | High | IndexMap/BTreeMap; assert order in tests. |
| HALF_EVEN/formatter rules | Med-High | Central util; replicate `--fisher` trailing-zero strip; ε=1e-6. |
| Tiling changes output near boundaries | Med | Halo + boundary-variant dedup; assert tiled==single-region on fixtures (§6). |
| `modifyCigar` crash inputs | Med | Bounds-check all; crash-corpus tests; graceful no-panic. |
| Coordinate 0/1 errors | High | Newtypes (compile error) + `/check-coordinates` + harness. |
| Cond-MLE OR mismatch | Med | Port Brent exactly; ε covers JDK-8227248. |
| noodles API churn | Med | Pin exact; golden tests; rust-htslib escape hatch. |
| Oracle build (JDK 8 vs host) | Med | Pinned conda/Docker; `run_oracle.sh` records version. |

## 14. Milestones & DAG
`P0→P1→P2→P3→P4→P5→P6(M1 single-sample)→{P7(M2 somatic),P8 amplicon,P9 SV(M3 full parity)}→P10 perf/scale→P11 robustness+CI/release→P12 validation.` Phases 7–9 share the P1–P6 engine; all built to equal depth. **v1.0 release = M3 + all §8 FIX-BY-DESIGN robustness + scale gates green.**

## 15. End-to-end verification (acceptance)
1. **Unit** (`cargo test`): Fisher matrices, CIGAR fixpoint/crash cases, 8-col BED math, HALF_EVEN.
2. **Tier-1** (`compare_var.py` vs `raw_input/*.var` + oracle) → engine parity.
3. **Tier-2 (primary gate)** (`compare.sh --all` vs `expected/*.vcf` + SEQC2/GIAB slices) — all modes (M1–M3).
4. **Robustness:** every crash-corpus case (#341 family, CRAM #249/#380, missing index #187) runs to completion with correct output (no panic).
5. **Scale:** peak RSS bounded regardless of `-R` size; ~linear scaling to 8 cores (A3); pre-registered performance budget: Rust wall <= Java wall, Rust max RSS <= Java max RSS (A9).
6. **Truth-set (secondary):** som.py/rtg vcfeval/truvari precision-recall ≈ v1.8.3.
7. **Hygiene:** `cargo fmt --check`, `clippy -D warnings`, `cargo audit`, `cargo deny check`.
8. **Distribution:** cargo-dist (`dist`) for static musl release + gated OIDC publish (A8); `--version` prints validated-against version + AI disclosure.

## 16. Provenance, licensing, attribution
MIT. README cites **Lai et al., NAR 2016**, links AstraZeneca-NGS/VarDict[Java], states validated-against **v1.8.3 (0886261)** + submodule `009e017`, lists **fixed issues** (§8) and **AI disclosure** (tools/role/validation coverage). Upstream bug reports only after manual reproduction with the original tool (never AI-generated cases).

## 17. Resolved decisions
- **Crate layout:** Cargo **workspace** — `crates/vardict-core` (engine lib) + `crates/vardict-cli` (bin `vardictrust`).
- **Oracle:** **bioconda** `vardict-java=1.8.3` in conda env `vardict-oracle`.
- **Distribution:** **open-source MIT**, cargo-dist (`dist`) for static musl release + gated OIDC publish (A8).
- **CLID:** **drop-in + ergonomic extras** — all VarDictJava short flags/defaults **plus** long aliases & subcommands, **direct-to-VCF**; `--emit-intermediate` reproduces the `.var` for legacy R/Perl piping.
- Locked defaults: name `vardictrust`; edition 2021; MSRV = current stable (`rust-toolchain.toml`); VCF 4.2; threads default = all cores (output deterministic regardless); tiling auto with `--tile-size` advanced; post-1.0 output fixes as individual opt-in flags.

---
---

# PART II — EXECUTION MANUAL (paint-by-numbers)

> **Audience:** an autonomous coding agent of *any* capability, including small/old models. Execute the tasks in **Appendix C in order**. Each task is self-contained: do exactly what it says, run its command, confirm the **Expect** line, then go to the next task. Never skip a validation. If a validation fails, fix **that task only** before continuing (the rewrites.bio loop). All "Java" paths are under `VarDictJava/src/main/java/com/astrazeneca/vardict/`; Perl/R under `VarDictJava/VarDict/`. The reference oracle is the conda env **`vardict-oracle`**.

## A. Global rules for the executor (read once; obey on every task)
1. **One task at a time.** A task = one focused change + its validation. Do not start the next task until the current one's **Expect** passes.
2. **After every code change, run in this order** (all must pass, zero warnings):
   `cargo fmt --all` → `cargo clippy --all-targets --all-features -- -D warnings` → `cargo test`.
3. **No panics in library code.** In `crates/vardict-core/src/**` (outside `#[cfg(test)]`) you may **not** use `unwrap()`, `expect()`, `panic!`, `unreachable!`, slice indexing `x[i]` without a prior bounds check, or integer ops that can overflow on input. Use `Result<_, VardictError>` + `?`, `.get(i)?`, and checked/saturating arithmetic. (This is what fixes the crash-issue family by design.)
4. **Coordinates are typed.** BAM/0-based values become `Pos0`; VCF/VarDict/1-based values become `Pos1`. Convert only via `.to_one_based()` / `.to_zero_based()`. Never add a raw integer to a position without going through these types. If a value's basis is unclear, **stop and re-read the Java** before coding.
5. **Numbers are rounded the VarDict way.** Any float written to a `.var`/VCF field goes through `round_half_even` (Appendix B7), never Rust's default formatting. Replicate the no-`--fisher` vs `--fisher` formatter rule (§4.3).
6. **Iteration order is part of correctness.** Use `indexmap::IndexMap` where Java uses `LinkedHashMap` and `BTreeMap` where Java uses `TreeMap`. Never iterate a `HashMap` in an output path.
7. **Definition of done (per task):** the task's stated `compare_*.py` prints `MATCH`, or its `cargo test`/command passes. Then commit: `git add -A && git commit -m "Txxx: <short desc>"`. Branch off `main`; never commit secrets.
8. **When porting a function:** (a) read the cited Java lines; (b) write one sentence in the Rust doc-comment describing what it does and cite the source (`/// ports CigarModifier.java:61-71`); (c) implement the given signature; (d) add the unit test named in the task; (e) validate. If the Java is ambiguous, port it literally first, make the test match the oracle, then simplify.
9. **Unsupported flag/mode:** `return Err(VardictError::Unsupported("--foo not implemented in v1.0; use VarDictJava"))`. Never silently ignore a flag.
10. **Do NOT add these crates** (they are hallucinated or unfit): `genomicframe-core`, `rustalign-simd`, `ksw2rs`, `seqair`, `cyanea-align`. Alignment = `bio` (+ `block-aligner` for SIMD later). Stats = `statrs` + `fishers_exact`.

## B. Scaffolding — exact contents (Phase 0 creates these verbatim)

### B0. Prerequisites — run once, confirm each prints (Ubuntu 12c/24GB/amd64)
```bash
# 1) Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
. "$HOME/.cargo/env"
rustup component add rustfmt clippy && rustup target add x86_64-unknown-linux-musl
rustc --version            # expect: rustc 1.8x.x
# 2) Conda (Miniforge) for the reference oracle
test -d "$HOME/miniforge3" || { wget -qO /tmp/mf.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && bash /tmp/mf.sh -b; }
source "$HOME/miniforge3/etc/profile.d/conda.sh"
# 3) Oracle + samtools (also installs the R/Perl reference scripts)
bash harness/setup_oracle.sh
conda run -n vardict-oracle samtools --version | head -1    # expect: samtools 1.x
```
**Do not start T001 until `rustc`, `cargo`, `conda`, and `samtools` all print versions.** If one is missing, install it — never work around a missing prerequisite.

**Reference genomes & test data.** Each of the 94 `VarDictJava/testdata/integrationtestcases/` filenames encodes its reference (`GRCh37.fa`, `hg19.fa`, or a tiny bundled fasta):
- **No download (start here):** cases naming `hard_clip_case.fa`, `example_genome`, or `test` — fastas are in `VarDictJava/testdata/fastas/`.
- **Needs a reference:** cases naming `GRCh37.fa` / `hg19.fa`. Fetch once, build the `.fai`, point `-G` at it, and **never commit it**:
  ```bash
  mkdir -p refdata && cd refdata
  wget -O hg19.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz && gunzip hg19.fa.gz
  conda run -n vardict-oracle samtools faidx hg19.fa
  ```
  Ensure contig naming matches the BAM (`chr1` vs `1`). `harness/fetch_data.sh` records exact URLs + md5 for the SEQC2/GIAB sets.

**`.gitignore`** (create in Phase 0, task T002; **commit `Cargo.lock`**, never commit data):
```gitignore
/target
/refdata/
/benchmarks/data/
*.bam
*.bai
*.cram
*.crai
*.pyc
.venv/
# keep the small golden fixtures even though they are .vcf/.var
!golden/
!golden/**
```

### B1. `rust-toolchain.toml`
```toml
[toolchain]
channel = "stable"
components = ["rustfmt", "clippy"]
targets = ["x86_64-unknown-linux-musl"]
```

### B2. Cargo manifests
`Cargo.toml` (workspace root):
```toml
[workspace]
members = ["crates/vardict-core", "crates/vardict-cli"]
resolver = "2"

[workspace.package]
version = "0.1.0"
edition = "2021"
license = "MIT"
repository = "https://github.com/palakpsheth/VarDictRust"
authors = ["Palak Sheth"]

[workspace.dependencies]
anyhow = "1"
thiserror = "2"
clap = { version = "4", features = ["derive"] }
rayon = "1"
indexmap = "2"
noodles = { version = "0.*", default-features = false, features = ["bam","cram","sam","fasta","vcf","bgzf","csi","tabix","core","bed"] }
statrs = "0.18"
fishers_exact = "1"
bio = "2"
log = "0.4"
env_logger = "0.11"

[profile.release]
lto = "thin"
codegen-units = 1
panic = "abort"
```
> In Phase 0 you will run `cargo add` (T004) to resolve and **pin exact** versions; replace each `0.*`/major-only spec with the resolved `=x.y.z` after `cargo add`.

`crates/vardict-core/Cargo.toml`:
```toml
[package]
name = "vardict-core"
version.workspace = true
edition.workspace = true
license.workspace = true

[dependencies]
thiserror.workspace = true
rayon.workspace = true
indexmap.workspace = true
noodles.workspace = true
statrs.workspace = true
fishers_exact.workspace = true
bio.workspace = true
log.workspace = true

[dev-dependencies]
insta = "1"
proptest = "1"
criterion = "0.5"

[[bench]]
name = "engine"
harness = false

[lints.rust]
unsafe_code = "forbid"
```
`crates/vardict-cli/Cargo.toml`:
```toml
[package]
name = "vardict-cli"
version.workspace = true
edition.workspace = true
license.workspace = true

[[bin]]
name = "vardictrust"
path = "src/main.rs"

[dependencies]
vardict-core = { path = "../vardict-core" }
clap.workspace = true
anyhow.workspace = true
env_logger.workspace = true
log.workspace = true
```

### B3. `CLAUDE.md` (full — place at repo root; symlink `AGENTS.md` → `CLAUDE.md`)
```markdown
# CLAUDE.md — VarDictRust
VarDictRust is an idiomatic-Rust rewrite of VarDictJava v1.8.3 (somatic/germline variant caller),
following the rewrites.bio manifesto. GOAL: a faster, scalable tool that produces output IDENTICAL to
VarDictJava v1.8.3 on valid inputs. Reference source: ./VarDictJava (v1.8.3, commit 0886261).

## 0. Evidence before completion (meta-rule)
Never say "done" without evidence: a passing `harness/compare.sh` MATCH, a passing `cargo test`, or a
diff that shows match. "Should work" is not evidence.

## 1. Think before coding   ## 2. Simplicity first   ## 3. Surgical changes   ## 4. Goal-driven
(Karpathy 4 — full text mirrored from the repo's CLAUDE.md; ask when unclear; minimum non-speculative
code; every changed line traces to the task; turn each task into a verifiable success criterion and loop.)

## 5. Prime directive — Emulate exactly (v1.0)
Output must match VarDictJava v1.8.3 on valid inputs: identical variants, FILTER, INFO/FORMAT keys, GT,
TYPE, STATUS, and column counts/order (Simple 36/38, Amplicon 38/40, Somatic 55/61). Floats within 1e-6
after 5-dp HALF_EVEN rounding. If v1.8.3 has an output bug, replicate it in v1.0 and open a tracked issue
(do not "fix" silently). Where v1.8.3 CRASHES, do not crash — produce the same output it would on valid
inputs (this is allowed and is the point).

## 6. Coordinate indexing — #1 footgun
BAM positions = 0-based half-open → type `Pos0`. VCF/VarDict-internal = 1-based → type `Pos1`. Convert
ONLY via `.to_one_based()`/`.to_zero_based()`. Never mix raw ints. Unsure? STOP and re-read the Java.

## 7. Fidelity hazards
- Iteration order: `IndexMap` (was LinkedHashMap), `BTreeMap` (was TreeMap). No HashMap in output paths.
- Rounding: use `round_half_even`; replicate the no-fisher vs --fisher formatter rule.
- Sequences are `&[u8]`/`Vec<u8>`, never String. Replace Java negative-index substr/charAt with
  bounds-checked slice helpers returning Option.
- Quality is Phred+33: `Phred(byte - 33)`.

## 8. The port loop
Pick ONE function -> read its Java source -> port to the given Rust signature + unit test ->
run harness/compare.sh (or cargo test) -> MATCH=commit+next ; MISMATCH=debug THIS one, change nothing else.

## 9. Definition of done (every change)
cargo fmt --all && cargo clippy --all-targets --all-features -- -D warnings && cargo test all pass;
relevant harness/compare.sh prints MATCH; golden updated if output changed; PR/commit states
"validated vs VarDictJava 1.8.3, dataset <name>, cmd <...>".

## 10. Commands
Build: cargo build           Test: cargo test [name]      Lint: cargo clippy --all-targets -- -D warnings
Fmt: cargo fmt --all         Oracle: bash harness/run_oracle.sh ...   Diff: bash harness/compare.sh <case|--all>

## 11. Do NOT use (hallucinated/unfit) crates
genomicframe-core, rustalign-simd, ksw2rs, seqair, cyanea-align. Use noodles / bio / block-aligner / statrs.

## 12. Scale rules
Bounded memory at any depth: process region TILES with a halo; stream reads (drop promptly); sort only
within a tile; rayon work-stealing with ordered output. Tiled output MUST equal single-region output.

## 13. Unsupported features error loudly
Return VardictError::Unsupported(...) — never silently drop a flag or emit wrong output.
```

### B4. Harness (place in `harness/`, `chmod +x` the .sh)
`setup_oracle.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
# Provision the VarDictJava 1.8.3 reference oracle via bioconda.
conda env list | grep -q vardict-oracle || \
  conda create -y -n vardict-oracle -c bioconda -c conda-forge \
    vardict-java=1.8.3 vardict=2019.06.04 samtools bcftools
echo "oracle ready: conda run -n vardict-oracle vardict-java --version"
```
`run_oracle.sh` (single-sample or paired; writes VCF):
```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: run_oracle.sh <ref.fa> <bam | tumor.bam|normal.bam> <bed> <out.vcf> [sample] [extra VarDict args...]
REF=$1; BAM=$2; BED=$3; OUT=$4; SAMPLE=${5:-SAMPLE}; shift 5 || shift $#
AF=${AF_THR:-0.01}
RUN(){ conda run -n vardict-oracle "$@"; }
if [[ "$BAM" == *"|"* ]]; then  # paired/somatic
  RUN vardict-java -G "$REF" -f "$AF" -N "$SAMPLE" -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$@" "$BED" \
    | RUN testsomatic.R | RUN var2vcf_paired.pl -N "$SAMPLE" -f "$AF" > "$OUT"
else
  RUN vardict-java -G "$REF" -f "$AF" -N "$SAMPLE" -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$@" "$BED" \
    | RUN teststrandbias.R | RUN var2vcf_valid.pl -N "$SAMPLE" -E -f "$AF" > "$OUT"
fi
RUN vardict-java --version 2>/dev/null || true
```
`compare.sh` (Tier-2 VCF differential gate):
```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: compare.sh <ref> <bam|tumor|normal> <bed> [sample] [extra args...]
REF=$1; BAM=$2; BED=$3; SAMPLE=${4:-SAMPLE}; shift 4 || shift $#
TMP=$(mktemp -d)
bash harness/run_oracle.sh "$REF" "$BAM" "$BED" "$TMP/java.vcf" "$SAMPLE" "$@"
cargo run -q --release -p vardict-cli -- -G "$REF" -f "${AF_THR:-0.01}" -N "$SAMPLE" \
  -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$@" "$BED" > "$TMP/rust.vcf"
python3 harness/compare_vcf.py "$TMP/java.vcf" "$TMP/rust.vcf"; rc=$?
[ $rc -eq 0 ] && echo "MATCH" || { echo "MISMATCH"; diff "$TMP/java.vcf" "$TMP/rust.vcf" | head -50; }
rm -rf "$TMP"; exit $rc
```
`compare_var.py` — **Tier-1** intermediate `.var` comparator (row order preserved):
```python
#!/usr/bin/env python3
"""Tier-1: compare two VarDict intermediate .var tables (TSV, no header).
Field-by-field, row order preserved: numeric within --eps after 5-dp HALF_EVEN, else exact.
--cols A-B restricts to the 1-based inclusive column range A..B. Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse
from decimal import Decimal, ROUND_HALF_EVEN

def r5(s):
    try: return float(Decimal(s).quantize(Decimal("0.00001"), rounding=ROUND_HALF_EVEN))
    except Exception: return None

def eq(a, b, eps):
    if a == b: return True
    ra, rb = r5(a), r5(b)
    return ra is not None and rb is not None and abs(ra - rb) <= eps

def load(p):
    with open(p) as f:
        return [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
    ap.add_argument("--cols", default=None, help="1-based inclusive range, e.g. 8-13")
    a = ap.parse_args()
    J, R = load(a.java), load(a.rust)
    if len(J) != len(R):
        print(f"MISMATCH: row count java={len(J)} rust={len(R)}"); sys.exit(1)
    lo, hi = 0, None
    if a.cols:
        s, e = a.cols.split("-"); lo, hi = int(s) - 1, int(e)
    errs = []
    for i, (j, r) in enumerate(zip(J, R)):
        if len(j) != len(r):
            errs.append(f"row {i}: field count java={len(j)} rust={len(r)}")
        else:
            top = min(hi, len(j)) if hi else len(j)
            for c in range(lo, top):
                if not eq(j[c], r[c], a.eps):
                    errs.append(f"row {i} col {c+1}: {j[c]!r} != {r[c]!r}")
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    print(f"MATCH ({len(J)} rows)"); sys.exit(0)

if __name__ == "__main__":
    main()
```

`compare_vcf.py` — **Tier-2** VCF comparator (records aligned by `CHROM/POS/REF/ALT`, so output order is ignored):
```python
#!/usr/bin/env python3
"""Tier-2: compare two VCFs for VarDictRust <-> VarDictJava equivalence.
Exact: CHROM/POS/REF/ALT, FILTER, INFO key-set, FORMAT key-set, per-sample non-numeric values (incl GT).
Tolerant: numeric values within --eps after 5-dp HALF_EVEN. Records aligned by key -> order ignored.
Exit 0=MATCH, 1=MISMATCH."""
import sys, argparse
from decimal import Decimal, ROUND_HALF_EVEN

def r5(s):
    try: return float(Decimal(s).quantize(Decimal("0.00001"), rounding=ROUND_HALF_EVEN))
    except Exception: return None

def eq(a, b, eps):
    if a == b: return True
    ra, rb = r5(a), r5(b)
    return ra is not None and rb is not None and abs(ra - rb) <= eps

def kv(field):
    d = {}
    if field in (".", ""): return d
    for p in field.split(";"):
        k, sep, v = p.partition("="); d[k] = v if sep else True
    return d

def load(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith("#"):
                out.append(ln.rstrip("\n").split("\t"))
    return out

def keyof(r): return (r[0], r[1], r[3], r[4])  # CHROM, POS, REF, ALT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("java"); ap.add_argument("rust"); ap.add_argument("--eps", type=float, default=1e-6)
    a = ap.parse_args()
    J = {keyof(r): r for r in load(a.java)}
    R = {keyof(r): r for r in load(a.rust)}
    if set(J) != set(R):
        oj = sorted(map(str, set(J) - set(R)))[:10]; orr = sorted(map(str, set(R) - set(J)))[:10]
        print(f"MISMATCH: variant sets differ. java-only(<=10)={oj} rust-only(<=10)={orr}"); sys.exit(1)
    errs = []
    for k in sorted(J):
        j, r = J[k], R[k]; w = f"{k[0]}:{k[1]} {k[2]}>{k[3]}"
        if j[6] != r[6]: errs.append(f"{w} FILTER: {j[6]!r} != {r[6]!r}")
        if not eq(j[5], r[5], a.eps): errs.append(f"{w} QUAL: {j[5]!r} != {r[5]!r}")
        ij, ir = kv(j[7]), kv(r[7])
        if set(ij) != set(ir):
            errs.append(f"{w} INFO keys: only-java={sorted(set(ij)-set(ir))} only-rust={sorted(set(ir)-set(ij))}")
        else:
            for key in ij:
                if ij[key] is True or ir[key] is True:
                    if ij[key] != ir[key]: errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
                elif not eq(str(ij[key]), str(ir[key]), a.eps):
                    errs.append(f"{w} INFO {key}: {ij[key]!r} != {ir[key]!r}")
        if len(j) > 9 and len(r) > 9:
            if j[8] != r[8]:
                errs.append(f"{w} FORMAT keys: {j[8]!r} != {r[8]!r}")
            else:
                fkeys = j[8].split(":")
                for si in range(9, min(len(j), len(r))):
                    for key, xj, xr in zip(fkeys, j[si].split(":"), r[si].split(":")):
                        if not eq(xj, xr, a.eps): errs.append(f"{w} sample{si-9} {key}: {xj!r} != {xr!r}")
        if len(errs) >= 20: break
    if errs:
        print("MISMATCH:"); [print("  " + e) for e in errs[:20]]; sys.exit(1)
    print(f"MATCH ({len(J)} variants)"); sys.exit(0)

if __name__ == "__main__":
    main()
```

> **Oracle nondeterminism:** VarDictJava can vary output **ordering** run-to-run (issue #282; see `VarDictJava/tests/integration/test_non_deterministic_behavior`). `compare_vcf.py` aligns by key, so ordering never causes a false mismatch. For Tier-1 `.var` (compared by row order), `sort` both files on the leading columns before diffing on any case known to be order-nondeterministic. **Never change the engine to chase a pure ordering difference.**
`slice.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: slice.sh <in.bam> <chr:start-end> <out.bam>   (requires samtools; subsets without loading whole file)
conda run -n vardict-oracle samtools view -b "$1" "$2" > "$3" && conda run -n vardict-oracle samtools index "$3"
```
`fetch_data.sh`: documents + downloads SEQC2 HCC1395 WES and GIAB HG008 slices (NCBI FTP `…/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/`), then calls `slice.sh` around benchmark loci; writes a `benchmarks/MANIFEST` with URLs, md5, and exact commands.

### B5. CI workflows (`.github/workflows/`, concise but complete)
- `ci.yml` (on push/PR): checkout → `dtolnay/rust-toolchain@stable` → `cargo fmt --all -- --check` → `cargo clippy --all-targets --all-features -- -D warnings` → `cargo test --all` → `cargo build --release`.
- `differential.yml` (on PR + nightly): setup miniconda → `bash harness/setup_oracle.sh` → build release → run `harness/compare.sh --all` over the `golden/` cases (sliced) → upload diffs as artifacts on failure.
- `bench.yml` (on PR): `CodSpeedHQ/action` running `cargo codspeed` (criterion harness) on deep-WGS slices; also record peak RSS; fail on >X% regression.
- `security.yml` (on push/PR + schedule): `cargo audit` and `cargo deny check`; (optional) Snyk Code SAST with SARIF upload.
- `release.yml` (on tag `v*`): `dist` builds `x86_64-unknown-linux-musl` (+ mac/win) archives → GitHub Release; then crates.io publish via `rust-lang/crates-io-auth-action` (OIDC trusted publishing); `release-plz` PR for version/changelog.

### B6. CLI (`crates/vardict-cli/src/cli.rs`) — clap derive, drop-in + ergonomic extras, direct-to-VCF
```rust
use clap::Parser;
/// VarDictRust: single-binary VarDict variant caller (VarDictJava 1.8.3-compatible).
#[derive(Parser, Debug)]
#[command(name = "vardictrust", version, about)]
pub struct Cli {
    /// Reference FASTA (.fai required). Aliases: --reference. Accepts hg19/hg38/mm10 shortcuts.
    #[arg(short = 'G', long = "reference")] pub reference: String,
    /// Indexed BAM/CRAM. Somatic: "tumor.bam|normal.bam". Alias: --bam.
    #[arg(short = 'b', long = "bam")] pub bam: String,
    /// Allele-frequency threshold. Alias: --af-threshold. [default: 0.01]
    #[arg(short = 'f', long = "af-threshold", default_value_t = 0.01)] pub freq: f64,
    /// Sample name (somatic: "tumor|normal"). Alias: --name.
    #[arg(short = 'N', long = "name")] pub sample: Option<String>,
    /// Region chr:start-end[:gene] (instead of BED).
    #[arg(short = 'R', long = "region")] pub region: Option<String>,
    /// BED columns (1-based): chr/start/end/gene.
    #[arg(short = 'c')] pub col_chr: Option<usize>,
    #[arg(short = 'S')] pub col_start: Option<usize>,
    #[arg(short = 'E')] pub col_end: Option<usize>,
    #[arg(short = 'g')] pub col_gene: Option<usize>,
    /// Threads (absent = all cores). Alias: --threads.
    #[arg(long = "threads", short = 't' )] pub threads: Option<usize>,
    /// Local realignment 0/1. [default: 1]
    #[arg(short = 'k', default_value_t = 1)] pub realign: u8,
    /// Amplicon mode "int:float". [default: 10:0.95 when 8-col BED]
    #[arg(short = 'a', long = "amplicon")] pub amplicon: Option<String>,
    /// Min variant reads [2]; base-qual [22.5]; min mapq filter; read-position [5]; max mismatches [8].
    #[arg(short = 'r', default_value_t = 2)] pub min_reads: u32,
    #[arg(short = 'q', default_value_t = 22.5)] pub base_qual: f64,
    #[arg(short = 'Q')] pub map_qual: Option<u32>,
    #[arg(short = 'P', default_value_t = 5.0)] pub read_pos: f64,
    #[arg(short = 'm', default_value_t = 8)] pub max_mismatches: u32,
    /// Min normal AF for somatic [0.05]; SV off; SV min len [1000].
    #[arg(short = 'V', default_value_t = 0.05)] pub min_normal_af: f64,
    #[arg(short = 'U', long = "nosv", default_value_t = false)] pub no_sv: bool,
    #[arg(short = 'L', default_value_t = 1000)] pub sv_min_len: u32,
    /// Emit the intermediate .var table instead of VCF (legacy R/Perl compatibility / Tier-1 testing).
    #[arg(long = "emit-intermediate", default_value_t = false)] pub emit_intermediate: bool,
    /// Use the in-engine Fisher test (always true here; flag kept for compatibility).
    #[arg(long = "fisher", default_value_t = false)] pub fisher: bool,
    /// Advanced: region tile size (bp) for bounded-memory WGS. [auto]
    #[arg(long = "tile-size")] pub tile_size: Option<u32>,
    /// Positional BED file (when -R absent).
    pub bed: Option<String>,
    // ... remaining VarDictJava flags (-x,-X,-3,-B,-o,-O,-I,-z,-F,-w,-W,-A,-Y,-Z,--chimeric,--deldupvar,
    //     -mfreq,-nmfreq,-J,-j,-i/--splice,-p,-DP) added one-per-task in Phase 1, each with the §4.2 default.
}
```

### B7. Core types (`crates/vardict-core/src/types.rs`) — newtypes, enums, rounding
```rust
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)] pub struct Pos0(pub u32); // 0-based (BAM)
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)] pub struct Pos1(pub u32); // 1-based (VCF/VarDict)
impl Pos0 { pub fn to_one_based(self) -> Pos1 { Pos1(self.0 + 1) } }
impl Pos1 { pub fn to_zero_based(self) -> Pos0 { Pos0(self.0 - 1) } }
#[derive(Clone, Copy, PartialEq, Eq, Debug)] pub struct Phred(pub u8); // base-33 already removed
impl Phred { pub fn from_char(b: u8) -> Self { Phred(b.saturating_sub(33)) } }

/// Variant description grammar (matches VarDict's +seq / -len / &seq / ^seq / #..^.. tokens).
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum VariantDesc { Snv(u8), Mnv(Vec<u8>), Ins(Vec<u8>), Del(u32),
                       DelIns { del: u32, ins: Vec<u8> }, Complex(Vec<u8>) }
impl std::fmt::Display for VariantDesc { /* reproduce exact token strings; cite ToVarsBuilder grammar */ }

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum VarType { Snv, Mnv, Insertion, Deletion, Complex, Dup, Inv, Del } // Display -> "SNV","MNV",...
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SomaticStatus { Germline, StrongSomatic, LikelySomatic, StrongLOH, LikelyLOH,
                         AFDiff, SampleSpecific, Deletion } // Display -> exact VarDict strings

/// Banker's rounding to `digits` decimals (matches Java DecimalFormat HALF_EVEN).
pub fn round_half_even(value: f64, digits: i32) -> f64 {
    if !value.is_finite() { return value; }
    let f = 10f64.powi(digits); let x = value * f; let r = x.round();
    let r = if (x - x.floor() - 0.5).abs() < f64::EPSILON {           // exactly .5 -> nearest even
        let fl = x.floor(); if (fl as i64) % 2 == 0 { fl } else { fl + 1.0 } } else { r };
    r / f
}

#[derive(thiserror::Error, Debug)]
pub enum VardictError {
    #[error("I/O: {0}")] Io(#[from] std::io::Error),
    #[error("unsupported: {0}")] Unsupported(String),
    #[error("parse: {0}")] Parse(String),
    #[error("{0}")] Other(String),
}
```

### B8. Worked example — the shape of EVERY task in Appendix C (imitate this)
Task **T016** (BED 3/4-col parsing), done end to end:
- **Read** `VarDictJava/src/main/java/com/astrazeneca/vardict/RegionBuilder.java:15-16` + `data/Region.java` — how a BED row maps to chr/start/end/gene, and that VarDict stores positions **1-based** internally.
- **Do** — create `crates/vardict-core/src/regions.rs` (`BedCols` = the small struct of 0-based column indices created in T015):
```rust
use crate::types::{Pos1, VardictError};

/// ports RegionBuilder.java:15-16 (BED row -> region)
#[derive(Debug, PartialEq, Eq)]
pub struct Region { pub chr: String, pub start: Pos1, pub end: Pos1, pub gene: String }

/// Parse a BED file. `cols` gives the 0-based column indices for chr/start/end/(optional)gene.
pub fn parse_bed(path: &str, cols: BedCols) -> Result<Vec<Region>, VardictError> {
    let text = std::fs::read_to_string(path)?;
    let mut out = Vec::new();
    for line in text.lines().filter(|l| !l.is_empty() && !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let u32_at = |i: usize| f.get(i).and_then(|s| s.parse::<u32>().ok())
            .ok_or_else(|| VardictError::Parse(line.to_string()));
        let start0 = u32_at(cols.start)?;             // BED start is 0-based, half-open
        let end = u32_at(cols.end)?;
        out.push(Region {
            chr: f.get(cols.chr).ok_or_else(|| VardictError::Parse(line.to_string()))?.to_string(),
            start: Pos1(start0 + 1),                   // -> VarDict 1-based
            end: Pos1(end),
            gene: cols.gene.and_then(|g| f.get(g)).unwrap_or(&"").to_string(),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bed3() {
        std::fs::write("/tmp/three.bed", "chr7\t55259400\t55259600\n").unwrap();
        let r = parse_bed("/tmp/three.bed", BedCols::default()).unwrap();
        assert_eq!(r, vec![Region { chr: "chr7".into(), start: Pos1(55259401), end: Pos1(55259600), gene: String::new() }]);
    }
}
```
- **Cmd** `cargo test -p vardict-core regions::tests::bed3`
- **Expect** `test regions::tests::bed3 ... ok`
- **Then** `cargo fmt --all && cargo clippy --all-targets -- -D warnings` (clean) → `git add -A && git commit -m "T016: BED 3/4-col parsing"`.

Every §C task has this shape: study the cited Java → write Rust + a test encoding the expected behavior → run the exact command → confirm the exact **Expect** → commit. **If you cannot state the Expect concretely, you don't understand the task yet — re-read the Java before coding.**

## C. Ordered task list (execute top-to-bottom; each: **Read** Java → **Do** Rust → **Cmd** → **Expect**)

> Format: `Read` = source to study (skip if "—"); `Do` = exact Rust file + item to create/port; `Cmd` = how to validate; `Expect` = what proves success. After each green task: `git add -A && git commit -m "Txxx: ..."`.

### Phase 0 — Scaffolding (T001–T014)
- **T001** Branch. Do: `git checkout -b init-scaffold`. Cmd: `git branch --show-current`. Expect: `init-scaffold`.
- **T002** Dirs. Do: `mkdir -p crates/vardict-core/src crates/vardict-cli/src harness golden benchmarks tests/integration .github/workflows`. Cmd: `ls crates`. Expect: `vardict-cli vardict-core`.
- **T003** Manifests. Do: write B2 files (`Cargo.toml` ×3) + `rust-toolchain.toml` (B1). Cmd: `cargo metadata --no-deps -q | head -c1`. Expect: prints (workspace parses).
- **T004** Pin deps. Do: in each crate run `cargo add <dep>@<resolved>` for every B2 dep (e.g. `cargo add noodles --features bam,cram,sam,fasta,vcf,bgzf,csi,tabix,core,bed`); then replace `0.*`/major specs with the exact resolved `=x.y.z`. Cmd: `cargo build`. Expect: builds (empty crates).
- **T005** Core stubs. Do: `crates/vardict-core/src/lib.rs` with `pub mod types;` + the other module decls (`cli? no`, `config, regions, reference, bam, cigar, realign, sv, tovars, variations, stats, modes, postprocess, output, util, patterns, scope, global, error`), each as an empty `pub mod` file. Cmd: `cargo build`. Expect: builds.
- **T006** Types. Do: write B7 into `types.rs`; add `#[test]` for `round_half_even(2.5,0)==2.0 && round_half_even(3.5,0)==4.0` and `Pos0(0).to_one_based()==Pos1(1)`. Cmd: `cargo test -p vardict-core types`. Expect: pass.
- **T007** CLI. Do: write B6 `cli.rs` + `crates/vardict-cli/src/main.rs` (parse `Cli`, `env_logger::init()`, call a stub `vardict_core::run(cfg)` returning `Ok(())`). Cmd: `cargo run -p vardict-cli -- --help`. Expect: help lists `-G -b -f -N -R ... --emit-intermediate`.
- **T008** CLAUDE.md. Do: write B3 to `CLAUDE.md`; `ln -s CLAUDE.md AGENTS.md`. Cmd: `test -f AGENTS.md`. Expect: exists.
- **T009** LICENSE + README. Do: MIT `LICENSE`; `README.md` crediting Lai et al. 2016 + AstraZeneca-NGS, stating "validated vs VarDictJava 1.8.3", AI-use disclosure section, supported-flag matrix stub. Cmd: `grep -c "Lai" README.md`. Expect: ≥1.
- **T010** Harness. Do: write all B4 scripts; `chmod +x harness/*.sh`. Cmd: `bash -n harness/compare.sh`. Expect: no syntax error.
- **T011** Oracle. Do: `bash harness/setup_oracle.sh`. Cmd: `conda run -n vardict-oracle vardict-java 2>&1 | head -1`. Expect: VarDict usage/version text.
- **T012** Seed golden. Do: copy `VarDictJava/tests/integration/expected/vardict.simple.vcf|vardict.paired.vcf` and `raw_input/raw.vardict.{simple,paired}.var` into `golden/`; write `golden/PROVENANCE` = "VarDictJava 1.8.3 / commit 0886261 / submodule 009e017". Cmd: `ls golden`. Expect: 4 files + PROVENANCE.
- **T013** CI. Do: write B5 workflows. Cmd: `for f in .github/workflows/*.yml; do python3 -c "import yaml,sys;yaml.safe_load(open(sys.argv[1]))" "$f"; done`. Expect: no error.
- **T014** Subagents/commands. Do: add `.claude/commands/{port-function,compare-outputs,check-coordinates,snapshot-update}.md` per §10; stop-gate hook in `.claude/settings.json` that blocks "done" without harness evidence. Cmd: `ls .claude/commands`. Expect: 4 files. **Commit; open PR; M0 done.**

### Phase 1 — I/O, regions, tiling, parallelism (T015–T032)
- **T015** Error wiring. Read: `Configuration.java`. Do: `config.rs` `Config` struct mirroring §4.2 fields+defaults; `From<Cli>`. Cmd: `cargo test config`. Expect: defaults match §4.2 (`freq==0.01`, `min_reads==2`, `amplicon=="10:0.95"` default, `flag_filter=="0x504"`).
- **T016** BED 3/4-col. Read: `RegionBuilder.java:15-16`, `data/Region.java`. Do: `regions.rs` `Region{chr,start:Pos1,end:Pos1,gene}` + `parse_bed(path,cols)->Vec<Region>`. Cmd: `cargo test regions::tests::bed3`. Expect: pass (intervals match).
- **T017** BED 8-col amplicon. Read: `RegionBuilder.java:17,104-146`, `VarDictLauncher.java:133-178`. Do: extend parser to detect 8 numeric cols 7&8 → `Region{amplicon:Some{insert_start,insert_end}}` using `AMP_BED_ROW_FORMAT=(0,1,2,6,7,3)`. Cmd: `cargo test regions::tests::bed8_amplicon`. Expect: insert interval computed.
- **T018** Region (`-R`). Read: `CmdParser.java:82`. Do: parse `chr:start-end[:gene]`. Cmd: `cargo test regions::tests::region_flag`. Expect: pass.
- **T019** Tiling. Read: §6. Do: `regions.rs::tile(region, tile_size, halo)->Vec<Tile>` (overlapping windows + halo). Cmd: `cargo test regions::tests::tiling_covers_and_overlaps`. Expect: union == region, neighbors overlap by halo.
- **T020** Reference. Read: `data/Reference.java`, `ReferenceResource.java`. Do: `reference.rs` indexed FASTA fetch via `noodles_fasta` (.fai); `fetch(chr,Pos1,Pos1)->Vec<u8>`. Cmd: `cargo test reference -- --ignored` (uses `VarDictJava/testdata/fastas`). Expect: known bases returned.
- **T021** Seed map. Read: `ReferenceResource.java:107-125`. Do: build SEED_1=17 & SEED_2=12 k-mer→positions map over the fetched window. Cmd: `cargo test reference::tests::seedmap`. Expect: unique seed resolves to one position.
- **T022** BAM open + region query. Read: `data/SamView.java`, `SAMFileParser.java`. Do: `bam.rs` indexed reader (`noodles_bam` + .bai/.csi), `query(chr,Pos0,Pos0)->impl Iterator<Record>`; `thread_local!` reader cache. Cmd: `cargo test bam -- --ignored`. Expect: record count matches `samtools view -c`.
- **T023** Record preprocess/filters. Read: `RecordPreprocessor.java`, `SamView.java:32`. Do: flag filter `(flags & 0x504)!=0 → skip`, MapQ/dedup/downsample, supplementary handling. Cmd: `cargo test bam::tests::flag_filter`. Expect: filtered counts match oracle on a slice.
- **T024** CRAM. Read: issue #249/#380 notes. Do: enable `noodles_cram` with `-G` reference; synthesize NM/MD if absent. Cmd: `cargo test cram -- --ignored`. Expect: a CRAM slice yields same records as its BAM twin (**fixes #249/#380**).
- **T025** Parallel runner. Read: `AbstractMode.java:112-143`. Do: `runner.rs` rayon over tiles, collect `(idx,Vec<u8>)`, emit in idx order, bounded backpressure. Cmd: `cargo test runner::tests::ordered_output`. Expect: output order == input order regardless of thread count.
- **T026** Determinism guard. Read: issues #282/#181/#356. Do: ensure all output sorts are stable with full tie-break keys. Cmd: run a fixture twice, `diff`. Expect: identical (**fixes nondeterminism**).
- **T027–T032** Add remaining CLI flags (`-x,-X,-3,-B,-o,-O,-I,-z,-F,-w,-W,-A,-Y,-Z,--chimeric,--deldupvar,-mfreq,-nmfreq,-J,-j,-i,-p,-DP`) one per task, each wired to `Config` with its §4.2 default + a parse test. Expect: each `cargo test config::flag_*` passes. **Commit; M-IO done.**

### Phase 2 — CIGAR modify + pileup (T033–T058)
- **T033** Cigar model. Read: `data/ModifiedCigar.java`, `Patterns.java`. Do: `cigar/mod.rs` `Cigar(Vec<(u32,CigarOp)>)` + parse/Display; `util.rs` bounds-checked `subseq(&[u8],i,len)->Option<&[u8]>`, `char_at(&[u8],isize)->Option<u8>` (negative-index-safe). Cmd: `cargo test cigar::tests::roundtrip`. Expect: parse∘display identity.
- **T034** modify: strip D. Read: `CigarModifier.java:61-71`. Do: `cigar/modify.rs::strip_leading_trailing_del`. Cmd: `cargo test cigar::modify::strip_del`. Expect: pass.
- **T035** modify: I→S. Read: `:73-82`. Do: `lead_trail_ins_to_softclip`. Cmd: `cargo test cigar::modify::ins_to_sclip`. Expect: pass.
- **T036** modify: chimeric 5'/3' softclip removal. Read: `:86-130`. Do: `remove_chimeric_softclip` (seed-unique, SEED_2). Cmd: `cargo test cigar::modify::chimeric`. Expect: pass.
- **T037** modify: indel-cleanup fixpoint loop. Read: `:132-228`. Do: `cleanup_indels` `while flag` over the sub-patterns. Cmd: `cargo test cigar::modify::cleanup`. Expect: pass.
- **T038** modify: complex folds. Read: `:486-786` (combineToCloseToOne/Correct, two/three-deletions, threeIndels, beginDigitM...). Do: port each helper; all walk via the bounds-checked `char_at`. Cmd: `cargo test cigar::modify::complex`. Expect: pass.
- **T039** modify: mis-softclip recovery. Read: `:358-478`. Do: `capture_mis_softly_ms` + `_3mismatches`. Cmd: `cargo test cigar::modify::missoft`. Expect: pass.
- **T040** modify: leading-clip merge. Read: `:254-352`. Do: `combine_dig_s_dig_m` + `combine_begin_dig_m`. Cmd: `cargo test cigar::modify::leadmerge`. Expect: pass.
- **T041** modify: crash fixtures. Read: the `TwoDelsInsStringBoundError` + `hard_clip_next_to_del*` testdata cases. Do: add tests feeding those exact CIGAR/seq inputs. Cmd: `cargo test cigar::modify::crash_corpus`. Expect: **no panic**, output equals oracle (**fixes #341/#243/#253/#283 family**).
- **T042** Variation struct. Read: `variations/Variation.java`. Do: `variations/mod.rs` `Variation{ cnt, fwd, rev, pmean, qmean, mqmean, nm, hicnt, locnt, pstd_flag, qstd_flag, ... }` + `inc_dir`, `add_cnt` (`:1407-1424`), `sub_cnt`. Cmd: `cargo test variations::tests::counts`. Expect: counters update like Java.
- **T043** Sclip/Cluster/Mate. Read: `variations/{Sclip,Cluster,Mate}.java`. Do: structs; `Sclip.nt: BTreeMap<Pos1,BTreeMap<u8,i32>>`, `seq: BTreeMap<...>`. Cmd: `cargo build`. Expect: builds.
- **T044** CigarParser scaffold. Read: `CigarParser.java:216-264`. Do: `cigar/parse.rs::pileup(records, ref, region)->VariationData`; `direction = record.flags().is_reverse_complemented()`. Cmd: `cargo build`. Expect: builds.
- **T045** parse: mismatch gate. Read: `:223-239`. Do: drop read if `NM-(I+D) > max_mismatches`. Cmd: `cargo test parse::tests::mismatch_gate`. Expect: pass.
- **T046** parse: cleanupCigar + run modify. Read: `:265-310,1530-1573`. Do: optional `CigarModifier`, then decode; drop reads soft-clipped >10bp both ends; track `max_read_length`. Cmd: `cargo test parse::tests::cleanup`. Expect: pass.
- **T047** parse: M base loop + addVariation. Read: `:374-610,1674-1740`. Do: per-base match/mismatch counting into `Variation` (+pstd/qstd flags, refCoverage). Cmd: `cargo test parse::tests::match_loop`. Expect: pass.
- **T048** parse: MNV growth. Read: `:411-473`. Do: extend variant over consecutive mismatches (qual≥goodq; vext look-ahead). Cmd: `cargo test parse::tests::mnv`. Expect: pass.
- **T049** parse: deletion (D). Read: `:651-851,1745-1793`. Do: `process_deletion` (`-len`, `-len&seq`). Cmd: `cargo test parse::tests::del`. Expect: pass.
- **T050** parse: insertion (I). Read: `:898-1102`. Do: `process_insertion` (`+seq`, adjInsPos homopolymer shift, ref-count subtraction so AF≤1). Cmd: `cargo test parse::tests::ins`. Expect: pass.
- **T051** parse: splice (N). Read: `:1289-1302`. Do: add to `splice` set + spliceCount. Cmd: `cargo test parse::tests::splice`. Expect: pass.
- **T052** parse: softclip. Read: `:1108-1283,1916-1976`. Do: `process_soft_clip` (walk matched-clipped bases back into alignment; store remainder into `Sclip`). Cmd: `cargo test parse::tests::softclip`. Expect: pass.
- **T053–T058** Tier-1 engine check (pre-realign): run `cargo run -p vardict-cli -- ... --emit-intermediate` on the smallest `Simple` testdata slice; compare the *match/indel counts* columns against the oracle's `.var`. Cmd: `python3 harness/compare_var.py --cols 8-13 golden/raw.vardict.simple.var <ours>`. Expect: depth/strand columns MATCH on no-realign cases. **Commit; M-pileup done.**

### Phase 3 — Realignment (T059–T075)
- **T059** VariationUtils. Read: `VariationUtils.java:69-154,262-334`. Do: `variations/util.rs` `findconseq`, `adjCnt`, `strandBias`, `correctCnt`. Cmd: `cargo test util::tests::findconseq`. Expect: consensus accept iff match/total>0.9 + heuristics.
- **T060** ismatch / ismatchref. Read: `VariationRealigner.java:2093-2125,2688-2725`. Do: `realign.rs::ismatch` (≤3 mm AND <15%). Cmd: `cargo test realign::tests::ismatch`. Expect: 15% gate enforced.
- **T061** findbp / findbi. Read: `:2337-2387,2212-2326`. Do: breakpoint scans (indelsize=50, maxmm=3−n/100, ≥3 distinct). Cmd: `cargo test realign::tests::findbp`. Expect: pass.
- **T062** findMM3 / findMM5. Read: `:2506-2622`. Do: flanking ≤3-consecutive-mismatch walkers. Cmd: `cargo test realign::tests::findmm`. Expect: pass.
- **T063** realigndel. Read: `:407-654`. Do: port (sort by count desc then pos; adjCnt transfers; `-N&seq` merge). Cmd: `cargo test realign::tests::realigndel`. Expect: pass.
- **T064** realignins. Read: `:661-962`. Do: port (wupseq/sanpseq flanks; long-insertion correction `NEWINS`). Cmd: `cargo test realign::tests::realignins`. Expect: pass.
- **T065** realignlgdel. Read: `:970-1363`. Do: large-deletion softclip realign (findconseq len≥7, findbp/findMatch, partialPipeline pull, re-run realigndel). Cmd: `cargo test realign::tests::lgdel`. Expect: pass.
- **T066** realignlgins30 + realignlgins. Read: `:1368-,1592-`. Do: paired 5'+3' clip reconstruction; `<dupN>`. Cmd: `cargo test realign::tests::lgins`. Expect: pass.
- **T067** SV helpers shared. Read: `StructuralVariantsProcessor.java:1884-1972`. Do: `sv.rs::find_match`, `markSV`, `isOverlap`. Cmd: `cargo test sv::tests::find_match`. Expect: pass.
- **T068–T075** Realignment integration: run the `Colo829.chr10.*` (CandDel/LGDEL/LGINS) + EGFR-type slices with realign on; Tier-1 compare full `.var`. Cmd: `python3 harness/compare_var.py golden/<case>.var <ours>`. Expect: MATCH (esp. the EGFR ex19 15-bp del row). **Commit; M-realign done.**

### Phase 4 — Finalization, isGoodVar, 36-col schema (T076–T088)
- **T076** Variant struct. Read: `variations/Variant.java`. Do: `variations/variant.rs` with all output fields (DP,VD,refF/R,altF/R,gt,af,bias,pmean,pstd,qual,qstd,mq,sn,hiaf,adjaf,shift3,msi,msilen,nm,hicnt,hicov,lseq,rseq,vartype,duprate,sv). Cmd: `cargo build`. Expect: builds.
- **T077** createVariant / createInsertion. Read: `ToVarsBuilder.java:322-488`. Do: `tovars.rs` (AF=`varsCount/ttcov` via `round_half_even`; means; qratio=`hicnt/(locnt|0.5)`; replicate hicov-excludes-insertions `:534`). Cmd: `cargo test tovars::tests::create`. Expect: fields match oracle row.
- **T078** sortVariants / collect. Read: `:275-310`. Do: sort by `meanQuality*posCoverage` desc, tie-break by desc string (stable). Cmd: `cargo test tovars::tests::sort`. Expect: deterministic order.
- **T079** MSI / shift3. Read: `:549-589`. Do: `find_msi` windowed motif scan; `shift3`. Cmd: `cargo test tovars::tests::msi`. Expect: pass.
- **T080** collectReferenceVariants + grammar. Read: `:598-967`. Do: build refallele/varallele/genotype; parse desc grammar (`+ - & ^ #`) into `VariantDesc`; `<DEL>/<DUP>/<INV>` spelling; IUPAC reduce `:1023-1031`. Cmd: `cargo test tovars::tests::desc_grammar`. Expect: `Display` reproduces VarDict tokens.
- **T081** isGoodVar. Read: `Variant.java:321-370`. Do: exact filter (§4.6) incl `freq>0.30` early-accept, MSI cuts, strand `"2;1"`. Cmd: `cargo test variant::tests::isgoodvar`. Expect: matches truth table.
- **T082** SimplePostProcess. Read: `SimplePostProcessModule.java:32-107`. Do: `postprocess/simple.rs`. Cmd: `cargo test postprocess::simple`. Expect: pass.
- **T083** 36-col printer. Read: `SimpleOutputVariant.java:184-227`. Do: `output/table.rs::simple_36` with the exact field order (§4.3) + formatter rule (no-fisher: `x==0?0:fmt`). Cmd: `python3 harness/compare_var.py golden/raw.vardict.simple.var <ours>`. Expect: **byte-MATCH** (mod Fisher cols).
- **T084–T088** Tighten until the entire single-sample `.var` matches across all `Simple` testdata cases (`--emit-intermediate`, no `--fisher`). Cmd: loop `compare_var.py` over Simple cases. Expect: all MATCH. **Commit; M-schema done.**

### Phase 5 — Fisher (T089–T096)
- **T089** Hypergeometric density. Read: `FisherExact.java:57-75`. Do: `stats/fisher.rs::logdc_dhyper` via `statrs Hypergeometric.ln_pmf`, round 7dp HALF_EVEN. Cmd: `cargo test fisher::tests::logdc`. Expect: matches a hand-checked value.
- **T090** dnhyper / mnhyper. Read: `:103-134`. Do: log-sum-exp normalized pmf; mean. Cmd: `cargo test fisher::tests::dn_mn`. Expect: pass.
- **T091** Two-sided p-value. Read: `:171-184`. Do: Σ{probs ≤ observed·(1+1e-7)}; `round_as_r` 5dp. Cmd: `cargo test fisher::tests::pvalue`. Expect: equals R `fisher.test` on fixed 2×2s (±1e-5).
- **T092** One-sided p. Read: `:186-197`. Do: pnhyper greater/less. Cmd: `cargo test fisher::tests::onesided`. Expect: pass.
- **T093** zeroinC (Brent). Read: `UnirootZeroIn.java:22-101`. Do: `stats/zeroin.rs` port; tol `f64::EPSILON.powf(0.25)`. Cmd: `cargo test zeroin::tests::roots`. Expect: finds roots of test fns.
- **T094** Conditional-MLE odds ratio. Read: `FisherExact.java:85-145`. Do: `mle(x)` via zeroin; Display `Inf`/integer/5dp. Cmd: `cargo test fisher::tests::oddratio`. Expect: equals R `estimate` (±1e-5).
- **T095** Cross-check. Do: assert our two-sided p == `fishers_exact` crate p (±1e-9) on random tables (proptest). Cmd: `cargo test fisher::tests::crosscheck`. Expect: pass.
- **T096** 38-col (`--fisher`). Read: `SimpleOutputVariant.java:135-179`. Do: insert PVALUE/ODDRATIO at 21–22 with `getRoundedValueToPrint` (trailing-zero strip). Cmd: `python3 harness/compare_var.py golden/<fisher case>.var <ours --fisher>`. Expect: MATCH (see testdata case `...-fisher.txt`). **Commit; M-stats done.**

### Phase 6 — VCF + var2vcf_valid → **M1** (T097–T108)
- **T097** VCF header. Read: `var2vcf_valid.pl` headers. Do: `output/vcf.rs` write VCF 4.2 header with all INFO/FORMAT (§4.6) and `##contig` lengths from the FASTA `.fai` (**fixes #320**). Cmd: `cargo test vcf::tests::header`. Expect: header keys == oracle header keys.
- **T098–T104** Port `var2vcf_valid.pl` filtering: one task per FILTER ID (`d3,v2,f0.02,p8,pSTD,q22.5,Q10,SN1.5,NM5.25,Bias,MSI12,...`) with its exact threshold + conditional (mapq only if af<0.8; pSTD only if pstd==0&&af<0.35&&!amp). Cmd: `cargo test vcf::filters::<id>`. Expect: each matches oracle FILTER assignment on a crafted row.
- **T105** Record emit. Do: serialize INFO (REFBIAS=`rfwd:rrev`, VARBIAS=`vfwd:vrev`, …) + FORMAT (`GT:DP:VD:AD:AF:RD:ALD`) via `noodles_vcf`; BGZF if `.gz`. Cmd: `cargo test vcf::tests::record`. Expect: pass.
- **T106** Wire `run()`. Do: `vardict-core::run(cfg)`: regions→tiles→rayon pipeline→VCF (or `.var` if `--emit-intermediate`). Cmd: `cargo run -p vardict-cli -- -G ref -f 0.01 -N S -b slice.bam -R chr:..-.. > out.vcf`. Expect: valid VCF.
- **T107** **M1 gate.** Cmd: `bash harness/compare.sh <ref> <Simple slice.bam> <bed>` for every `Simple` testdata case + 2 SEQC2/GIAB slices. Expect: **all MATCH** (Tier-2 VCF). 
- **T108** Snapshot. Do: `cargo insta test --accept` to store golden VCFs; commit. **Tag `v0.1.0-single-sample`. M1 done.**

### Phase 7 — Somatic/paired → **M2** (T109–T124)
- **T109** Two-BAM flow. Read: `SomaticMode.java:103-131`. Do: `modes/somatic.rs` (BAM1 to completion; BAM2 seeded with BAM1 `max_read_length`; shared `splice`). Cmd: `cargo build`. Expect: builds.
- **T110** 55-col printer. Read: `SomaticOutputVariant.java:279-341`. Do: `output/table.rs::somatic_55` (7+18+18+12, varLabel at 50). Cmd: `cargo test table::somatic55`. Expect: order matches §4.3.
- **T111** 61-col (`--fisher`). Read: `:208-274`. Do: per-sample pval/OR + total at 60–61. Cmd: `cargo test table::somatic61`. Expect: pass.
- **T112** isNoise. Read: `Variant.java:198-216`. Do: port (zeroing side-effect included). Cmd: `cargo test variant::isnoise`. Expect: pass.
- **T113** determinateType. Read: `SomaticPostProcessModule.java:346-369`. Do: LikelyLOH/LikelySomatic/Germline/AFDiff/StrongSomatic per §4.5 thresholds. Cmd: `cargo test somatic::determinate`. Expect: matches truth table.
- **T114** combineAnalysis. Read: `:385-476`. Do: joint re-pileup → Germline/drop/keep; SV guard `end-start>SVMINLEN`. Cmd: `cargo test somatic::combine`. Expect: pass.
- **T115** accept dispatch + one/both-sample. Read: `:55-128,156-336`. Do: Deletion/SampleSpecific/StrongLOH branches. Cmd: `cargo test somatic::accept`. Expect: pass.
- **T116** testsomatic.R parity. Do: somatic 3-Fisher (cols 10-13, 28-31, 8/9–26/27 two-sided-by-min). Cmd: `cargo test somatic::fisher3`. Expect: ±1e-5 vs R.
- **T117** var2vcf_paired. Read: `var2vcf_paired.pl`. Do: INFO `STATUS,SSF,SOR` + per-sample FORMAT; stringency gates. Cmd: `cargo test vcf::paired`. Expect: pass.
- **T118–T124** **M2 gate.** Cmd: `bash harness/compare.sh <ref> "tumor|normal" <bed>` vs `golden/vardict.paired.vcf` + SEQC2 HCC1395 WES slices. Expect: all MATCH. **Tag `v0.2.0-somatic`. M2 done.**

### Phase 8 — Amplicon (T125–T134)
- **T125** Amplicon parse path. Read: `CigarParser.java:1310-1367`. Do: `parse_cigar_with_amp` distance/overlap filter (`10:0.95`). Cmd: `cargo test parse::amp`. Expect: pass.
- **T126** AmpliconMode. Read: `modes/AmpliconMode.java`. Do: per-amplicon pileup + position→amplicon map. Cmd: `cargo build`. Expect: builds.
- **T127** AmpliconPostProcess + AMPBIAS. Read: `AmpliconPostProcessModule.java:32-272`. Do: `isAmpBiasFlag`, gvscnt/tvscnt/nocov; iterate-all-alleles (**#68**). Cmd: `cargo test postprocess::amplicon`. Expect: pass.
- **T128** 38/40-col printer. Read: `AmpliconOutputVariant.java:206-249,145-201`. Do: tail GDAMP/TLAMP/NCAMP/AMPFLAG. Cmd: `cargo test table::amplicon`. Expect: order matches §4.3.
- **T129–T134** **gate.** Cmd: `compare.sh` on the two `Amplicon;...` testdata cases. Expect: MATCH. **Commit.**

### Phase 9 — SV / Splicing → **M3** (T135–T146)
- **T135** SVStructures. Read: `data/SVStructures.java`. Do: structs. Cmd: build. Expect: builds.
- **T136** StructuralVariantsProcessor. Read: `StructuralVariantsProcessor.java` (`findDEL/findINV/findsv`). Do: `sv.rs` stage. Cmd: `cargo test sv::find`. Expect: pass.
- **T137** SplicingMode. Read: `modes/SplicingMode.java`, `CigarParser.java:101-116`. Do: `modes/splicing.rs` (emit splice counts). Cmd: `cargo test modes::splicing`. Expect: pass.
- **T138** SV VCF tags. Do: `SVTYPE/SVLEN/SPLITREAD/SPANPAIR/DUPRATE`; `-U/-L/-w/-W/-A/-Y`. Cmd: `cargo test vcf::sv`. Expect: pass.
- **T139–T146** **M3 gate.** Cmd: `compare.sh` on `NA24385.REF.SVDEL`, `Colo829.chr10.DEL*/LGINS*` cases; `truvari bench` concordant. Expect: MATCH. **Tag `v0.3.0-full-parity`. M3 done.**

### Phase 10 — Perf/scale (T147–T154)
- **T147** Bench harness. Do: `benches/engine.rs` (criterion) on a deep-WGS slice. Cmd: `cargo bench`. Expect: baseline recorded.
- **T148** AVX2 realign. Read: block-aligner docs. Do: feature-gate `block-aligner` for the SW hot path; keep `bio` fallback. Cmd: `compare.sh --all`. Expect: MATCH (no drift) + faster.
- **T149** Arenas/buffers. Do: reuse per-tile buffers; `SmallVec` in hot loops. Cmd: `compare.sh --all`. Expect: MATCH.
- **T150** Tile-size tuning + peak-RSS test. Do: `--tile-size` auto heuristic; test peak RSS bounded vs `-R` size. Cmd: a 10-Mb-region run with `/usr/bin/time -v`. Expect: peak RSS within budget (**fixes #197/#384 by design**).
- **T151–T154** Scaling test to 12 cores; CodSpeed gate; document throughput in `benchmarks/README`. Expect: ~linear scaling; CI bench passes. **Commit.**

### Phase 11 — Robustness corpus + CI/release (T155–T166)
- **T155** Crash corpus. Do: tests for every issue-#341 family input + missing-index (#187) + truncated bgzf (#61) + lenient SAM (#13/#28). Cmd: `cargo test robustness`. Expect: all run, no panic, correct/clear-error.
- **T156** CRAM e2e. Cmd: `compare.sh` on a CRAM slice vs BAM twin. Expect: MATCH (**#249/#380**).
- **T157–T161** Turn on all CI workflows; ensure green; add `cargo deny`/`audit` config; locale-independent formatting test (#384/#195). Expect: CI green.
- **T162** `dist` setup. Do: `dist init --ci=github`; target `x86_64-unknown-linux-musl`. Cmd: `dist plan`. Expect: plan prints musl target.
- **T163** Trusted publishing. Do: `release.yml` using `rust-lang/crates-io-auth-action`; configure crates.io trusted publisher. Cmd: dry-run. Expect: workflow validates.
- **T164–T166** `release-plz` PR automation; CHANGELOG/CONTRIBUTING; supported-flag matrix + **fixed-issues list (§8)** in README; AI-disclosure section. Expect: docs complete. **Commit.**

### Phase 12 — Validation tier + v1.0 (T167–T174)
- **T167** Full SEQC2 HCC1395 WES (manual). Cmd: `compare.sh` full WES tumor|normal. Expect: MATCH (or documented float-only diffs ≤1e-6).
- **T168** GIAB HG008 (PDAC) somatic. Cmd: `compare.sh` on HG008 slices. Expect: MATCH.
- **T169** Truth-set (secondary). Cmd: `som.py`/`rtg vcfeval`/`truvari bench` vs SEQC2/GIAB truth. Expect: precision/recall ≈ VarDictJava 1.8.3 (within documented delta).
- **T170** Acceptance §15 sweep. Cmd: run the 8 §15 checks. Expect: all pass.
- **T171** README provenance final pass (Lai 2016, AstraZeneca-NGS, validated version/dataset/commands). Expect: complete.
- **T172** Upstream issues (optional): for any defect surfaced, reproduce manually with the oracle, file with non-AI test case. Expect: documented.
- **T173** Tag **`v1.0.0`**; `dist`/crates.io publish via CI. Expect: GitHub Release + crates.io published.
- **T174** Post-1.0 backlog file: list the FIX-OUTPUT issues (§8) as individual opt-in-flag tasks for v1.1+. Expect: `ROADMAP.md` created.

### Definition of "v1.0 complete"
M3 reached (all four modes MATCH on the testdata corpus + SEQC2/GIAB slices) **and** all §8 FIX-BY-DESIGN robustness items covered (crash corpus green, CRAM works, bounded memory verified) **and** §15 acceptance sweep green **and** released under MIT to GitHub + crates.io.

---

## D. Glossary, grammar, and annotated samples (quick reference for any model)

### D1. Genomics terms
- **BAM/CRAM** — compressed aligned reads (CRAM is reference-compressed). **`.bai`/`.crai`** index for region queries; **`.fai`** indexes a FASTA.
- **CIGAR** — per-read alignment ops: `M` match/mismatch, `I` insertion, `D` deletion, `N` skip/intron, `S` soft-clip (bases present, not aligned), `H` hard-clip (bases absent). e.g. `30M2I68M`.
- **Pileup** — the per-reference-position tally of bases/indels across reads.
- **DP** total depth · **VD** variant depth · **AF/VAF** = VD/DP · **strand** fwd/rev (ref/alt counts = `.var` cols 10–13).
- **SNV / MNV / indel** — single / multi-nucleotide / insertion-or-deletion variant. **Soft-clip realignment** = recovering variant evidence from clipped read ends (VarDict's signature sensitivity).
- **Fisher's exact test** — 2×2 test → strand-bias / somatic p-value + odds ratio.
- **VCF** — variant output: **INFO** per-site (`KEY=VALUE;…`), **FORMAT**+sample cols per-genotype (`:`-joined), **FILTER** = PASS or a failure tag.
- **Somatic vs germline** — tumor-specific vs inherited; **STATUS** = VarDict's call (StrongSomatic, Germline, LikelyLOH, …).
- **0- vs 1-based** — BAM & BED are 0-based half-open; VCF & VarDict-internal are 1-based (hence `Pos0`/`Pos1`).

### D2. VarDict variant description-string grammar (`Variant.descriptionString`)
| Token | Meaning | Example |
|---|---|---|
| single base | SNV | `T` |
| `base&extra` | MNV | `T&CG` |
| `+SEQ` | insertion | `+ATG` |
| `-N` | deletion of N bases | `-3` |
| `-N&SEQ` | deletion + matched tail | `-3&AC` |
| `^SEQ` | deletion-then-insertion | `-2^TT` |
| `#SEQ…^…` | complex | — |
| `<DUP>`/`<INV>`/`<DEL>` | SV ≥ `-L` | `<DEL>` |

Model as `enum VariantDesc` whose `Display` reproduces these tokens **exactly** (Appendix B7).

### D3. Annotated intermediate `.var` row — single-sample, 36 cols, **tab-separated** (illustrative values; shape only)
```
SAMPLE  EGFR  chr7  55259515  55259515  T  G  812  9  400  403  4  5  G/G  0.0111  2;1  31.2  1  34.5  1  60.0  1.800  0.0125  0.0  0  0  0  0.5  9  812  ATGC..  ..GCTA  chr7:55259400-55259600  SNV  0.0  0
# 1 sample · 2 gene · 3 chr · 4 start · 5 end · 6 ref · 7 alt · 8 DP · 9 VD · 10 RefFwd · 11 RefRev · 12 AltFwd · 13 AltRev · 14 GT · 15 AF · 16 BIAS · 17 PMEAN · 18 PSTD · 19 QUAL · 20 QSTD · 21 MQ · 22 SN · 23 HIAF · 24 ADJAF · 25 SHIFT3 · 26 MSI · 27 MSILEN · 28 NM · 29 HICNT · 30 HICOV · 31 LSEQ · 32 RSEQ · 33 seg · 34 TYPE · 35 DUPRATE · 36 SV
```
(`--fisher` inserts SBF, ODDRATIO at cols 21–22 → 38 cols; somatic = 55/61; amplicon = 38/40 — see §4.3.)

### D4. Annotated VCF record — same variant after `var2vcf_valid` (illustrative; shape only)
```
chr7  55259515  .  T  G  .  PASS  SAMPLE=SAMPLE;TYPE=SNV;DP=812;VD=9;AF=0.0111;BIAS=2:1;REFBIAS=400:403;VARBIAS=4:5;PMEAN=31.2;PSTD=1;QUAL=34.5;SBF=0.9..;ODDRATIO=1.2..;MQ=60;SN=1.8;HIAF=0.0125;ADJAF=0;SHIFT3=0;MSI=0;MSILEN=0;NM=0.5;HICNT=9;HICOV=812;LSEQ=ATGC..;RSEQ=..GCTA  GT:DP:VD:AD:AF:RD:ALD  G/G:812:9:803,9:0.0111:400,403:4,5
```
Check equivalence with `compare_vcf.py` (Appendix B4) — never eyeball it.

