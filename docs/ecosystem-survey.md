# Ecosystem Survey — AI-assisted Rust rewrites of bioinformatics tools

**Author / Maintainer:** Palak Sheth <palakpsheth@gmail.com>  
**Repository:** [https://github.com/palakpsheth/VarDictRust](https://github.com/palakpsheth/VarDictRust)  
**Governing Manifesto:** [rewrites.bio](https://rewrites.bio)  

> **Point-in-time snapshot: 2026-08-24.** Compiled from web/GitHub sweeps (sources
> listed at the end; claims reflect pages as fetched that day). This document is
> provenance for decisions consolidated into project specifications and external validation reports in
> [`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md).
> Weaknesses are given equal standing to strengths — see §5.

## 1. The wave

AI coding assistants crossed a threshold where a domain expert can reimplement a
mature tool in days. **rewrites.bio** (Seqera, Apr 2026) is the movement's anchor
manifesto. Its principles:

1. **Credit the original authors** — README, papers, downstream reports; citations flow to the original.
2. **Emulate exactly** — byte-for-byte for deterministic tools; scientist-defined float tolerance otherwise; "everything counts": headers, column order, file naming.
3. **Be transparent about AI** — which tools, what role, how correctness was verified, where validation is thin.
4. **Think big / work small** — rethink architecture (fold pipeline steps, kill intermediate I/O), but implement in tight loops with output comparison at every step ("AI agents are unreliable for correctness — validate every step").
5. **Test and benchmark with real data**, documented hardware/datasets/commands.
6. **Build only what you need** — unsupported features must fail loudly pointing at the original; never silently ignore flags.
7. **Pin versions and document** — equivalence claims are about an exact version; document commands, data, dates, upstream-update policy.
8. **Stewardship** — maintain what you ship; preserve hot-swappable compatibility; open source; contribute upstream responsibly (verify manually; *never* AI-generated bug-report evidence).

Origin story: RustQC began as dupRadar rewrite over lunch, became a 15-tool QC suite at >60× speedup; weeks of parity iteration followed the weekend beta ("Stringency in validation will be key to community adoption and trust").

## 2. Exemplar rewrites and their evidence models

| Project | Replaces | Evidence model |
|---|---|---|
| [bwa-mem4](https://github.com/IPNP-BIPN/bwa-mem4) | bwa-mem2 2.3 | Byte-identity acceptance gate from first commit; whole-genome GIAB HG002 byte-identical SAM (353M SE + 707M PE records); 2000-job property test on SW kernel; release binaries must rebuild committed fixtures byte-identically per platform |
| [Bismark Rust suite](https://felixkrueger.github.io/Bismark/rust/overview/) | Perl Bismark (12 tools) | Byte-identical contract v0.25.1; worker-invariant `--multicore`; deliberate documented divergences; opt-in non-identical modes concordance-gated |
| [rammap](https://github.com/jwanglab/rammap) | minimap2 | Identical output; found 2 long-standing minimap2 bugs → fixed upstream v2.31 |
| [rustar-aligner](https://scverse.org/rustar-aligner/) | STAR | Public compatibility scoreboard; 99.88% mate faithfulness labeled "tie-adjusted"; 396 tests |
| [tama-rs](https://github.com/sguizard/tama-rs) | TAMA (Python 2) | Byte-identical goldens; dict-iteration-order exceptions compared as sets and documented |
| [RustQC](https://github.com/seqeralabs/RustQC) | 15 RNA-seq QC tools | nf-test snapshot harness as agent feedback loop; public AI & Provenance page incl. known gaps |
| [fastqc-compliant-rs](https://github.com/henriksson-lab/fastqc-compliant-rs) | Java FastQC | Golden char-for-char files; published parity matrix; refuses replacement claims until matrix green |
| [FastQC-Rust](https://github.com/ewels/FastQC-Rust) | Java FastQC | `UPSTREAM.toml` version pin; nightly CI bot files issues on upstream releases; equivalence report CI |
| [fgumi](https://github.com/fulcrumgenomics/fgumi) | fgbio (Scala) | 25× speedup, 100% output-reproduction claim |
| [fastars](https://github.com/necoli1822/fastars) | fastp | 100% identical sequences verified byte-for-byte; honest "limited samples" warning |

Also notable: sshash-rs/piscem (Patro, C++→Rust in <4 days, 100% concordance);
henriksson-lab's port program (minimap2, BLAST, HMMER, MUSCLE, DIAMOND, Kraken2…);
umi-collapse-rs, dist2mst-rs (~900×), EVidenceModeler_rust, bowtie2-rs ("parity
gates before core aligner work").

## 3. Community feedback

**Pro** (Heng Li, ["The AI Rewrite Dilemma"](https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma)): order-of-magnitude speedups for Python/R tools; static Rust binaries end dependency hell; welcomes rewrites of his own tools. Reversibility makes adoption risk-free — swap in, confirm on one run, roll back if needed.

**Con/skepticism** (Kai Blin thread): *"Who will maintain vibe-rewritten tools? Original maintainers won't."* Unmergeable AI PRs burden maintainers; fragmentation risk — subtle rounding/multimap differences poison cross-study comparability.

**Deeper critique** ([Clay McLeod](https://claymcleod.dev/blog/2026-07-13-beyond-bioinformatics-rewrites)): emulation preserves historical constraints without separating durable science from accidents; rewrites are seductively easy AI projects; noodles (spec-first foundation) may compound more than any single rewrite; **"a tool author can delegate implementation, but cannot delegate understanding"**; wants provenance-as-runtime-artifact.

**Counter-example teardown** ([polars-bio on RastQC](https://biodatageeks.org/polars-bio/blog/2026/07/04/cheap-to-rewrite-expensive-to-get-right-the-fastqc-compatible-trap/)): claimed FastQC-compatible with wrong core outputs, non-thread-invariant sampling, "speedup" that was parallelism-not-port (slower than the JVM single-threaded, silently gated above 50 MB). Thesis: *a rewrite is a validation contract; cheap generation needs disciplined validation* — spec-driven + test-driven development, golden assertions not header greps, partition-invariance checks, like-for-like benchmarks.

**Trust levers that convinced skeptics:** byte-diff against a pinned oracle on real data > any claim; published parity matrices; thread/partition-invariance proofs; manually-verified upstream bug reports; gating artifacts (not just source) before release.

## 4. How VarDictRust compares (strengths)

| Field practice | VarDictRust |
|---|---|
| Emulate exactly, bugs included | [`AGENTS.md`](../AGENTS.md) §5 bug-replication policy + [`DIVERGENCES.md`](../DIVERGENCES.md) ledger (rarer than byte-parity itself) |
| Pinned oracle + differential testing | VarDictJava v1.8.3 @`0886261` conda env; 92-case corpus; live-oracle VCF sweep |
| Work small, per-function evidence | §7 port loop with source citations — McLeod's knowledge-extraction ideal |
| Fail loudly on unsupported | §12 `Unsupported()` errors ([`cli.rs`](../crates/vardict-cli/src/cli.rs) `UNSUPPORTED_FLAGS`) |
| Real-data, like-for-like benchmarks | [`benchmarks/README.md`](../benchmarks/README.md); explicitly avoided the RastQC parallelism trap |
| Thread/partition invariance | Tested ([`thread_determinism.rs`](../crates/vardict-cli/tests/thread_determinism.rs), tile-span determinism) |
| Upstream bugs, manually verified | [`docs/upstream-bugs.md`](upstream-bugs.md) UB-1..4 drafted from live-oracle runs |
| Credit + AI transparency | [`README.md`](../README.md) Attribution/AI Disclosure; git release tags and evidence JSONs carry cryptographic provenance |

## 5. Gaps this survey exposed (and dispositions)

1. ~~No fuzzing~~ → Implemented `crates/vardict-cli/tests/robustness.rs` covering PR #402 hostile CIGAR sequences, missing indices, and zero-tile boundaries.
2. ~~Flag statuses asserted, not graded~~ → Added explicit validation table in `README.md` grading every CLI flag against VarDictJava v1.8.3.
3. ~~Release artifacts ungated~~ → Added Level 3 micro-fixture verification (`fixtures/micro/`) and static binary smoke test (`harness/smoke.sh`) gating the release artifact.
4. **External truth-set validation** → Completed Level 4 biological ground truth concordance suite ([`docs/validation/omnibus-empirical-validation.md`](validation/omnibus-empirical-validation.md)) evaluating GIAB HG002, SEQC2 HCC1395, GIAB HG008, GIAB HG005, and NA24385 b37 via `rtg vcfeval`.
5. Upstream-drift watcher: moot — VarDictJava v1.8.3 is frozen/archived; our pin (`0886261`) represents the fixed reference oracle contract.

## 6. Sources (as fetched 2026-08-24)

- https://rewrites.bio · https://github.com/seqeralabs/rewrites.bio · https://seqera.io/blog/rewrites-bio/ · https://seqera.io/blog/rustqc/
- https://seqeralabs.github.io/RustQC/about/ai-statement/
- https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma · https://lh3.github.io/2026/07/04/minibwa-is-the-new-bwa
- https://claymcleod.dev/blog/2026-07-13-beyond-bioinformatics-rewrites
- https://biodatageeks.org/polars-bio/blog/2026/07/04/cheap-to-rewrite-expensive-to-get-right-the-fastqc-compatible-trap/
- https://scverse.org/rustar-aligner/ · https://github.com/IPNP-BIPN/bwa-mem4/releases/tag/v4.0.0
- https://felixkrueger.github.io/Bismark/rust/overview/ · https://github.com/henriksson-lab/fastqc-compliant-rs
- https://github.com/ewels/FastQC-Rust · https://github.com/sguizard/tama-rs · https://github.com/necoli1822/fastars · https://github.com/jwanglab/rammap
- Phil Ewels, Nextflow Summit 2026 talk + launch threads (Kai Blin rebuttals on Bluesky/mu.social)
