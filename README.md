# VarDictRust

An idiomatic **Rust** rewrite of [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) v1.8.3 — a sensitive somatic and germline variant caller for next-generation sequencing (NGS) data.

VarDictRust consolidates the original 3-language pipeline (Java → R → Perl) into a **single, statically-compiled binary** that targets output **identical to VarDictJava v1.8.3** on valid inputs, while fixing crash/robustness issues by design and scaling to modern deep whole-genome sequencing.

Built under the [rewrites.bio](https://rewrites.bio) manifesto: *a faster, more robust tool that produces the same results*.

## Features

- **Single binary** — no Java, R, or Perl required
- **Engine-level output equivalent** to VarDictJava v1.8.3 (Tier-1 `.var` parity; see Validation status)
- **All four modes implemented**: Simple (germline), Somatic (paired tumor/normal), Amplicon, Splicing/SV
- **Crash-proof**: bounds-checked CIGAR processing fixes the entire #341 crash family
- **Scales to deep WGS**: bounded memory via region tiling, rayon parallelism
- **CRAM support** (fixes #249/#380)
- **Deterministic output** regardless of thread count (fixes #282)

## Quick Start

```bash
# Single-sample (germline)
vardictrust -G ref.fa -b sample.bam -f 0.01 -N SAMPLE regions.bed > output.vcf

# Somatic (paired tumor/normal)
vardictrust -G ref.fa -b "tumor.bam|normal.bam" -f 0.01 -N "TUMOR|NORMAL" regions.bed > output.var
```

## Validation status

Validated against VarDictJava v1.8.3, dataset = the VarDictJava/testdata
integration suite (92 cases):

- **Integration corpus: 91 PASS / 1 FAIL** — the sole FAIL is a stale expected
  fixture inside the read-only VarDictJava submodule (different reference
  version); live-oracle comparison proved our output byte-identical modulo two
  documented Java quirks (see DIVERGENCES.md "Known non-bugs"). This FAIL is
  permanent until upstream updates the fixture.
- **Upstream v1.8.3 bugs we replicate per §5** (with filing-ready issue texts):
  [docs/upstream-bugs.md](docs/upstream-bugs.md).
- **Tier-1 (intermediate `.var` tables): parity established** across the corpus.
- **Tier-2 (final VCF vs live oracle): parity established** — final W4.6 sweep
  (2026-08-23) of 9 representative slices across all four modes: 8 strict
  MATCH + 1 single-record AF-tie emission-order residue (NA24385 chr20:20202875;
  root-caused in P3 to Java HashMap iteration order across position keys and
  closed as a permanent known limitation). The two Colo829 SV slices are validated against
  `hg19.fa` — the reference their upstream test cases declare — after a fixture
  regen proved the bundled GRCh37.fa cannot carry those regions. Details:
  DIVERGENCES.md Category 14.
- **Milestone tags**: `milestone/somatic-parity`, `milestone/full-parity` (cut at W4.6
  closure; annotated with dataset + exceptions provenance).

## Supported Flags

Flag parity with VarDictJava v1.8.3 is defined against its option parser
(`CmdParser.java`) and verified by the differential-testing harness. Every
implemented flag, its long form, default, and **validation status**:

Status vocabulary (derived only from corpus case definitions + `DIVERGENCES.md`;
no flag is claimed validated beyond this evidence):

- `corpus` — active by default in every run, or explicitly varied in ≥1 of the
  92 integration-corpus cases / Tier-2 oracle sweep
- `unexercised` — parsed and wired, but never active in any corpus case;
  behavior follows the Java source but lacks direct differential evidence
- `no-output` — performance/infrastructure flag; cannot affect output bytes
- `parity-exempt` — accepted for CLI compatibility with no engine effect,
  matching upstream semantics

| Short | Long | Description | Default | Status |
|---|---|---|---|---|
| `-G` | `--reference` | Reference FASTA (`.fai` required; hg19/hg38/mm10 shortcuts accepted) | required | corpus |
| `-b` | `--bam` | Indexed BAM/CRAM; somatic: `"tumor.bam\|normal.bam"` | required | corpus |
| `-f` | `--af-threshold` | Allele-frequency threshold | `0.01` | corpus |
| `-N` | `--name` | Sample name (somatic: `"TUMOR\|NORMAL"`) | derived from BAM filename | corpus |
| `-R` | `--region` | Region `chr:start-end[:gene]` instead of BED | — | corpus |
| `-c` | — | BED chromosome column (1-based) | 3 | corpus |
| `-S` | — | BED start column | 7 | corpus |
| `-E` | — | BED end column | 8 | corpus |
| `-s` | — | BED segment-start column | 10 (= `-S` if set) | unexercised |
| `-e` | — | BED segment-end column | 11 (= `-E` if set) | unexercised |
| `-g` | — | BED gene column | 13 | corpus |
| — | `--threads` | Thread count (`-th`/`--th` aliases; bare `-th` = all cores) | all cores | no-output |
| `-k` | — | Local realignment 0/1 | `1` | corpus |
| `-a` | `--amplicon` | Amplicon mode `"int:float"` (distance:fraction) | auto `10:0.95` on 8-col BED | corpus |
| `-r` | — | Minimum variant reads | `2` | corpus |
| `-q` | — | Minimum base quality | `22.5` | corpus |
| `-Q` | — | Minimum mapping quality | — | corpus |
| `-P` | — | Mean read position filter | `5` | corpus |
| `-m` | — | Maximum mismatches | `8` | corpus |
| `-M` | — | Min match+insertion read length (shorter reads ignored) | `0` | unexercised |
| `-V` | — | Min normal AF (somatic) | `0.05` | corpus |
| `-U` | `--nosv` | Disable SV calling | off | unexercised |
| `-L` | — | SV minimum length | `1000` | corpus |
| — | `--emit-intermediate` | Emit legacy `.var` table instead of VCF | off | corpus |
| — | `--fisher` | In-engine Fisher exact test (p-value / odds ratio) | off | corpus |
| — | `--tile-size` | Region tile size in bp (bounded-memory WGS) | auto | no-output |
| `-x` | — | BED region extension (bp) | `0` | unexercised |
| `-X` | — | Nucleotide extension | `2` | corpus |
| `-3` | — | Move indels to 3′ (Java `moveIndelsTo3`) | off | unexercised |
| `-B` | — | Minimum reads for bias | `2` | corpus |
| `-o` | — | qratio filter | `1.5` | corpus |
| `-O` | — | Output MAPQ threshold | `0` | unexercised |
| `-I` | — | Indel size | `50` | corpus |
| `-z` | — | BED coordinates are 0-based | off (1-based) | corpus |
| `-F` | — | SAM flag filter | `0x504` | corpus |
| `-w` | — | SV read length | `300` | corpus |
| `-W` | — | SV flanking | `100` | corpus |
| `-A` | — | SV minimum AF | `4` | corpus |
| `-Y` | — | Max region length | `1200` | corpus |
| `-h` | `--header` | Print `.var` column header before results | off | corpus |
| `-v` | — | Accepted and ignored (never read by VarDictJava v1.8.3) | — | parity-exempt |
| — | `--chimeric` | Chimeric read handling | off | corpus |
| — | `--deldupvar` | Delete/dup variant handling | off | unexercised |
| — | `--mfreq` | Mono frequency threshold | `0.25` | corpus |
| — | `--nmfreq` | Normal mono frequency threshold | `0.1` | corpus |
| `-J` | — | CRISPR cutting site | — | unexercised |
| `-j` | — | CRISPR flanking size | — | unexercised |
| `-i` | `--splice` | Splice mode | off | unexercised |
| `-p` | — | Pileup-only mode | off | corpus |
| `-D` | `--debug` | Debug mode (append per-read debug columns to `.var`) | off | corpus |
| — | `--DP` | Debug print destination (`OUT` or `ERR`) | — | unexercised |
| `-y` | `--verbose` | Verbose mode | off | no-output |
| `-t` | `--dedup` | Remove duplicated reads (one pair per start position) | off | unexercised |
| — | `--adaptor` | Adaptor sequence(s), comma-separated, excluded from realignment | — | unexercised |
| `-d` | — | Delimiter for splitting region info | `\t` | corpus |
| `-u` | — | Unique mode: count overlapping mates once (forward read) | off | unexercised |
| — | `--UN` | Unique mode: count once using first read only | off | unexercised |
| `-K` | — | Include Ns in depth calculation | off | corpus |
| positional | — | BED file (used when `-R` absent) | — | corpus |

BED-column defaults follow VarDictJava's `RegionBuilder.DEFAULT_BED_ROW_FORMAT`
(0-based `2,6,7,9,10,12`), i.e. the multi-column VarDict pipeline BED format,
not plain 3-column BED — pass `-c 1 -S 2 -E 3` for simple BEDs.

**P1 flag-semantics remappings** (aligned with `CmdParser.java`, diverging from common Unix convention):

- `-t` is **dedup** (`--dedup`), not threads — threads is `-th`/`--threads` (`CmdParser.java:70` vs `617-628`)
- `-h` prints the **column header** (`--header`); help moved to `-H`/`--help` (`CmdParser.java:67` vs `222`)
- `-v` is accepted and ignored (defined but never read by VarDictJava v1.8.3, `CmdParser.java:224`)
- `-th`/`--th` normalize to `--threads`; absent `-th` → all cores (perf-only divergence from Java's 1-thread default; no output effect)

**Not implemented in v1.0** (rejected loudly with an actionable error, per CLAUDE.md §12):

`-n` (sample-name regexp), `-C` (numeric chromosomes, deprecated upstream),
`-T`, `--trim`, `--VS` (SAM validation stringency), `-Z`/`--downsample`.

See `vardictrust --help` for the generated reference.

## Attribution

This project is a faithful rewrite of **VarDictJava** by the [AstraZeneca-NGS](https://github.com/AstraZeneca-NGS) team.

> Lai Z, Markovets A, Ahdesmaki M, et al. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. *Nucleic Acids Research*, 2016;44(11):e108. doi:[10.1093/nar/gkw227](https://doi.org/10.1093/nar/gkw227)

Validated against VarDictJava v1.8.3 (commit `0886261`) with the Perl/R submodule at commit `009e017`.

## AI Disclosure

This project was developed with assistance from AI coding tools (Claude by Anthropic). Code was validated against the VarDictJava v1.8.3 reference implementation via two-tier differential testing:
- **Tier 1**: Intermediate `.var` table comparison (engine parity — established)
- **Tier 2**: Final VCF comparison against the live oracle (complete for the
  validated corpus and Tier-2 sweep; see Validation status)

The AI's role was implementation under human direction; all algorithmic decisions and validation criteria come from the original VarDictJava source code.

## License

[MIT](LICENSE)
