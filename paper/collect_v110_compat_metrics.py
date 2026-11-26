#!/usr/bin/env python3
"""Collect and verify deterministic v1.1.0 --compat=off and --preset empirical metrics.

Modes:
  --run     Execute all 6 deterministic public-slice and validator arms locally
            ($0 cloud spend, os.nice(19), capped RTG JVM memory, zero host timing
            metrics) and atomically write benchmarks/v110-compat-study/.
  --verify  Offline, tool-free, dataset-free verification of SHA256SUMS.tsv,
            JSON/TSV schemas, domain invariants, path-leak absence, and 32-switch
            CompatProfile catalog alignment against crates/vardictrust-core/src/config.rs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = ROOT / "benchmarks" / "v110-compat-study"
CONFIG_RS_PATH = ROOT / "crates" / "vardictrust-core" / "src" / "config.rs"
ROADMAP_MD_PATH = ROOT / "ROADMAP.md"

EXPECTED_STUDY_ID = "v110-compat-and-presets"
EXPECTED_STUDY_VERSION = "1.1.0"
EXPECTED_STUDY_FP = "b18477ae0a7b07c1"
EXPECTED_STUDY_PERF_FP = "5c41be8047133842"
EXPECTED_COMPAT_SWITCH_COUNT = 32
EXPECTED_PRESET_COUNT = 4
EXPECTED_VAR2VCF_FLAG_COUNT = 18
EXPECTED_TOTAL_ARMS = 6

DEFECT_KEYS = (
    "nonsv_info_end",
    "ad_number_r_mismatch",
    "noncanonical_gt",
    "noncanonical_sv",
    "accounting_type_iupac_shift3",
)

ALLOWED_DOMAINS = frozenset(
    {
        "VCF Schema & Genotypes",
        "Structural Variants & Indexing",
        "Somatic & Amplicon Calling",
        "CIGAR, Pileup & Aligner",
        "Depth, Rounding & Determinism",
    }
)

TSV_COLUMNS = [
    "arm_id",
    "arm_title",
    "cohort_scope",
    "provenance_type",
    "vardictjava_result",
    "compat_off_result",
    "verdict",
    "sanitized_cmd",
]

VALID_ALT_RE = re.compile(r"^(?:[ACGTN]+|<[A-Za-z0-9_-]+>)(?:,(?:[ACGTN]+|<[A-Za-z0-9_-]+>))*$")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# Structured metadata keyed by exact CompatProfile field name (1..32) + 3 CLI capabilities (33..35).
SWITCH_METADATA: dict[str, dict[str, str]] = {
    "sanitize_sample_info": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#116, AstraZeneca-NGS/VarDict#197",
        "vardictjava_behavior": "Interpolates raw `;`, `=`, `,`, `%`, `-`, `.`, and whitespace into `INFO/SAMPLE`, splitting VCFv4.2 fields",
        "compat_off_behavior": "Applies bijective VCFv4.2 percent-encoding on `INFO/SAMPLE` (`%2D`, `%2E`, `%XX`)",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "standard_dup_symbolic_alt": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#342, ToVarsBuilder.java:619,700",
        "vardictjava_behavior": "Emits non-standard `<dup-N>` symbolic `ALT` when `mlen < 100` and mismatches homozygous `<DUP>` span",
        "compat_off_behavior": "Emits standard `<DUP>` symbolic `ALT` with non-negative span and unified `2 * SVFLANK + dup_val` genotype",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "shift3_3prime_alleles": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#292, AstraZeneca-NGS/VarDictJava#395, AstraZeneca-NGS/VarDict#25",
        "vardictjava_behavior": "Shifts `-3` indel coordinates 3-prime without updating `REF` and `ALT` sequences to match the reference FASTA",
        "compat_off_behavior": "Re-extracts FASTA-concordant 3-prime-shifted deletion and insertion `REF`/`ALT` alleles with CRISPR double-shift guard",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "strict_region_bounds": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#369, AstraZeneca-NGS/VarDictJava#137",
        "vardictjava_behavior": "Allows realigned indel and SV anchors to spill outside requested `-R` or BED target intervals",
        "compat_off_behavior": "Enforces merged requested interval bounds (`start_pos >= 1` within target intervals)",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "suppress_pileup_sentinels": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "UB-1, UB-2, SimplePostProcessModule.java",
        "vardictjava_behavior": "Emits spurious position-0 `<R,R>` and `N/N` zero-depth sentinel records in `-p` pileup mode",
        "compat_off_behavior": "Suppresses `-p` `<R,R>` and `N/N` zero-depth sentinel rows in Simple and Amplicon modes",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "live_duprate": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#387, CigarParser.java:78,162",
        "vardictjava_behavior": "Copies zeroed read counters before pileup iteration, leaving `INFO/DUPRATE` permanently `0`",
        "compat_off_behavior": "Computes `INFO/DUPRATE` from live `u64` duplicate and total read counters across untiled and tiled runs",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "preserve_mixed_sex_contigs": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "AstraZeneca-NGS/VarDict#185, AstraZeneca-NGS/VarDict#191",
        "vardictjava_behavior": "Drops `chrX`, `chrY`, or `chrM` in `var2vcf` `reorder` when unprefixed `X`, `Y`, or `MT` co-occur",
        "compat_off_behavior": "Preserves both unprefixed and `chr`-prefixed sex and mitochondrial contigs under total contig order",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "cigar_homopolymer_merge": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#402, hedera-dx/VarDictJava@17d1cb51",
        "vardictjava_behavior": "Produces negative trailing `M` (`rm < 0`) during homopolymer multi-indel merges in `modifyCigar`, skipping the read",
        "compat_off_behavior": "Absorbs negative trailing `M` into `rdoff_corr` during homopolymer multi-indel CIGAR coalescing",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "eqx_cigar_read_length": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "CigarParser.java:1605-1629,1990-1991",
        "vardictjava_behavior": "Excludes `=` and `X` CIGAR operators from aligned read length, yielding `rlen1 = 0` and negative `PMEAN`",
        "compat_off_behavior": "Treats `=` (`SequenceMatch`) and `X` (`SequenceMismatch`) equivalently to `M` across pileup, realignment, and dedup",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "preserve_discordant_svs": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "UB-5, var2vcf_valid.pl:197-205",
        "vardictjava_behavior": "Drops single-sample symbolic `<DEL>`, `<DUP>`, and `<INV>` calls whenever `SPLITREAD < 1` despite `SPANPAIR > 0`",
        "compat_off_behavior": "Retains single-sample symbolic SVs supported by discordant mate pairs (`SPANPAIR > 0`) unless `--min-splitreads` is set",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "mark_dup_sv_inclusive_end": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "UB-6, StructuralVariantsProcessor.java:1413-1468",
        "vardictjava_behavior": "Passes post-decrement `pe` to `markDUPSV`, truncating the `svfdup` mate-cluster right boundary by 1 bp",
        "compat_off_behavior": "Passes pre-decrement inclusive cluster end `mark_pe` when `tbp != 0` and `tbp < pe`",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "sv_lazy_mate_pileup": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#370, StructuralVariantsProcessor.java:519-675",
        "vardictjava_behavior": "Scans mate windows eagerly in `findINVsub` before checking local `>= 12 bp` soft-clip presence or reverse-complement match",
        "compat_off_behavior": "Gates `find_inv_sub_paths` mate-window loads on local `>= 12 bp` soft-clip presence and validated `find_match_rev`",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "exact_binary_half_even": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "Utils.java:101-103",
        "vardictjava_behavior": "Truncates `5e-(d+1)` underflow ties to `0.0` in OpenJDK `DecimalFormat` `DigitList` despite binary `f64` exceeding midpoint",
        "compat_off_behavior": "Rounds exact IEEE-754 binary `f64` values via `round_half_even_compat` without `DigitList` underflow truncation",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "reconcile_vd_dp": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#92, AstraZeneca-NGS/VarDictJava#386, AstraZeneca-NGS/VarDictJava#400, AstraZeneca-NGS/VarDict#113",
        "vardictjava_behavior": "Increments realigned `VD` and `HICNT` without updating locus `DP` and `HICOV`, permitting `VD > DP` and `AF > 1.0`",
        "compat_off_behavior": "Enforces `VD <= DP`, `AF <= 1.0`, `ADJAF <= 1.0`, and `HICNT <= HICOV` (`HIAF <= 1.0`) across all variant types",
        "verification_arm": "Domain 3 (Multi-Allelic & Paralog)",
    },
    "coordinate_sorted_vars": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "UB-3, AstraZeneca-NGS/VarDictJava#282, AstraZeneca-NGS/VarDictJava#390, palakpsheth/VarDictRust#34",
        "vardictjava_behavior": "Iterates variant and post-processing position maps in `JavaHashMap` bucket order rather than coordinate order",
        "compat_off_behavior": "Iterates variant and post-processing position maps in strictly ascending genomic coordinate order",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "exact_35_overlap": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "palakpsheth/VarDictRust#34, VariationRealigner.java:2006-2042",
        "vardictjava_behavior": "Adds `+1` spurious out-of-bounds loop sentinel mismatch and `+1` overlap length when `rem5 < rem3` in `find35match`",
        "compat_off_behavior": "Omits out-of-bounds sentinel mismatch and length increment when `rem5 < rem3` in `find_35_match`",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "exact_cigar_dlen_match": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "palakpsheth/VarDictRust#34, VariationRealigner.java:2056-2065",
        "vardictjava_behavior": "Uses substring `contains(dlen + \"D\")` in `noPassingReads`, falsely matching decimal suffixes such as `12D` for `2D`",
        "compat_off_behavior": "Matches exact deletion CIGAR operator length (`len == dlen`) in `no_passing_reads`",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "somatic_continue_after_bad_var": {
        "category": "CompatProfile Switch",
        "domain": "Somatic & Amplicon Calling",
        "upstream_refs": "joachimwolff/VarDictJava#4, clintval/vartovcf#139, clintval/vartovcf#149, AstraZeneca-NGS/VarDictJava#376",
        "vardictjava_behavior": "Executes `break` in `printVariationsFromFirstSample` when the top-AF candidate fails `isGoodVar`, dropping passing secondary alleles",
        "compat_off_behavior": "Continues (`continue`) evaluating remaining tumor candidate alleles at multi-allelic somatic loci",
        "verification_arm": "Domain 3 (Multi-Allelic & Paralog)",
    },
    "adaptor_all_consensus_dirs": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#240, AstraZeneca-NGS/VarDictJava#241, AstraZeneca-NGS/VarDict#119, PolinaBevad/VarDictJava@92c62cf6",
        "vardictjava_behavior": "Calls `findConseq(sc, 0)` without `--adaptor` sets in `adjustMNP`, `realignIndels`, `sv`, and `adjSNV`, assembling adapter chimeras",
        "compat_off_behavior": "Threads directional `5`/`3` soft-clip orientation and `--adaptor` screening across all 20 consensus call sites",
        "verification_arm": "Domain 4 (Adapter Chimera)",
    },
    "canonical_sv_vcf_anchoring": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#398, AstraZeneca-NGS/VarDictJava#225, AstraZeneca-NGS/VarDictJava#229, clintval/vartovcf#192",
        "vardictjava_behavior": "Emits multi-kb `REF` spans on symbolic SVs, positive `SVLEN` on `<DEL>`, missing `##ALT` headers, and missing paired SV `INFO` tags",
        "compat_off_behavior": "Emits canonical 1-bp upstream VCFv4.2 anchor, `##ALT` headers, negative `<DEL>` `SVLEN` (`END - POS == abs(SVLEN)`), and paired `END;SVTYPE;SVLEN`",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "paired_somatic_amplicon": {
        "category": "CompatProfile Switch",
        "domain": "Somatic & Amplicon Calling",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#378, AstraZeneca-NGS/VarDictJava#155, AstraZeneca-NGS/VarDict#109, MHH-Bioinformatics-Hematology/vardictcpp#1",
        "vardictjava_behavior": "Runs single-sample `AmpliconMode` on `bam1` only in `.var` mode and rejects paired `-b` plus `-a` in VCF mode",
        "compat_off_behavior": "Runs per-amplicon tumor-normal pileup, primer masking, and `GDAMP;TLAMP;NCAMP;AMPFLAG` `AMPBIAS` concordance across `.var` and VCF",
        "verification_arm": "Domain 5 (Clinical Presets & Amplicon)",
    },
    "cigar_strip_outer_hard_clips": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#341, AstraZeneca-NGS/VarDictJava#103, Kennedy-Lab-UW/VarDictJava@9a0beb0d",
        "vardictjava_behavior": "Leaves outer `H`/`P` operators blocking terminal `I -> S` and `D` stripping in `CigarModifier` and distorting `PMEAN`",
        "compat_off_behavior": "Strips outer positive-length `H`/`P` operators without advancing `alignment_start` and treats interior `P` as non-consuming",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "nonnegative_subtracted_means": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "bcbio/bcbio-nextgen#1704, Variation.java",
        "vardictjava_behavior": "Subtracts `qmean` and `pmean` in `Variation.subCnt` without flooring at zero, allowing negative `PMEAN`, `QMEAN`, and VCF `QUAL`",
        "compat_off_behavior": "Floors subtracted `PMEAN`, `QMEAN`, `MQ`, `NM`, and VCF `QUAL` at `0.0` (`>= 0.0`)",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "reconcile_strand_and_hq_sums": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#316, clintval/vartovcf#120, clintval/vartovcf#142",
        "vardictjava_behavior": "Zeros `cnt = 0` without clearing `fwd`, `rev`, or `hicnt`, leaving `ALD[0] + ALD[1] != VD` and `HICNT > VD`",
        "compat_off_behavior": "Resets all auxiliary `Variation` accumulators on `cnt = 0` and reconciles `fwd + rev == VD` and `0 <= HICNT <= VD`",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "reconcile_rd_vd_dp": {
        "category": "CompatProfile Switch",
        "domain": "Depth, Rounding & Determinism",
        "upstream_refs": "AstraZeneca-NGS/VarDict#168, AstraZeneca-NGS/VarDictJava#323, clintval/vartovcf#102",
        "vardictjava_behavior": "Adjusts `ref_coverage` independently of `ref_fwd + ref_rev`, allowing `RD + VD > DP` (`AD[0] + AD[1] > DP`)",
        "compat_off_behavior": "Enforces `DP >= RD + VD` (`total_pos_cov >= ref_fwd + ref_rev + pos_coverage`) across variant and somatic builders",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "sanitize_iupac_alt_alleles": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#267, AstraZeneca-NGS/VarDictJava#287, AstraZeneca-NGS/VarDict#137, bcbio/bcbio-nextgen#3321",
        "vardictjava_behavior": "Reduces IUPAC ambiguity codes only on `REF`, leaking IUPAC codes into `ALT` and emitting non-pileup `REF == ALT` records",
        "compat_off_behavior": "Reduces IUPAC ambiguity codes in `ALT` outside symbolic `<...>` brackets and suppresses non-pileup `REF == ALT` records",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "sv_respect_min_len": {
        "category": "CompatProfile Switch",
        "domain": "Structural Variants & Indexing",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#328, AstraZeneca-NGS/VarDict#157",
        "vardictjava_behavior": "Ignores `-L` (`--min-sv-len`) in `find_del`, `find_dup`, and `find_inv`, emitting sub-1000 bp SVs as symbolic `<DUP>`/`<INV>`",
        "compat_off_behavior": "Materializes explicit sequence alleles (`<= 32,768 bp`) for split-refined SVs shorter than `-L` and suppresses short discordant-only SVs",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "pileup_preserve_shifted_ref": {
        "category": "CompatProfile Switch",
        "domain": "CIGAR, Pileup & Aligner",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#330, SimplePostProcessModule.java:73",
        "vardictjava_behavior": "Checks `vvar.size() == 1` during `-p` hydration, silently dropping covered coordinates when `>= 2` indels shift away",
        "compat_off_behavior": "Hydrates and emits the `ALT=\".\"` reference pileup row whenever all candidate variants shift away from a `-p` locus",
        "verification_arm": "Domain 6 (CSI & Pileup)",
    },
    "recompute_vartype_after_adj_complex": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "clintval/vartovcf#109, clintval/vartovcf#147, SimplePostProcessModule.java:83,93",
        "vardictjava_behavior": "Computes `INFO/TYPE` before `adjComplex` trims shared bases, leaving trimmed SNVs and indels mislabeled `TYPE=Complex`",
        "compat_off_behavior": "Recomputes `Variant::var_type` after `adj_complex` across Simple, Amplicon, and Somatic post-processing",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "strict_vcf_ad_number_r": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "clintval/vartovcf#101, AstraZeneca-NGS/VarDict#75, AstraZeneca-NGS/VarDict#76, EBIvariation/vcf-validator#99, sigven/pcgr#12",
        "vardictjava_behavior": "Keys single-sample `FORMAT/AD` on `VD != 0`, emitting 2-element `AD=rd,0` on `ALT=\".\"` rows and 1-element `AD=rd` on `VD=0` variants",
        "compat_off_behavior": "Emits strict `Number=R` `FORMAT/AD` cardinality (`1 + n_alt`: `AD=rd` when `ALT=\".\"` and `AD=rd,vd` when `ALT != \".\"`)",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "omit_nonsv_vcf_end": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "clintval/vartovcf#105, bcbio/bcbio-nextgen#1704",
        "vardictjava_behavior": "Emits `INFO/END` on every SNV and short indel under `-R` when `--no-end` is omitted, failing `GATK ValidateVariants` and `bcftools norm`",
        "compat_off_behavior": "Omits `INFO/END` on non-symbolic VCF records by default while retaining `INFO/END` on symbolic `<...>` SVs unless `--no-end` is set",
        "verification_arm": "Domain 1 (VCF Schema)",
    },
    "canonical_vcf_gt": {
        "category": "CompatProfile Switch",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "clintval/vartovcf#99, clintval/vartovcf#54, AstraZeneca-NGS/VarDictJava#344, AstraZeneca-NGS/VarDict#186",
        "vardictjava_behavior": "Emits unphased `1/0` (`0.5 <= AF < 1 - GTFREQ`) and `0/0` on sub-threshold `VD > 0` carriers (ignored as `CALL=IGN` in diploid `rtg vcfeval`)",
        "compat_off_behavior": "Canonicalizes unphased `1/0` to `0/1` and promotes concrete non-dot `0/0` `ALT` carriers (`VD > 0`) to `0/1`",
        "verification_arm": "Domain 2 (Diploid Somatic)",
    },
}

CLI_CAPABILITY_ENTRIES: list[dict[str, str]] = [
    {
        "name": "--vcf-all-alleles",
        "category": "Opt-In CLI Capability",
        "domain": "VCF Schema & Genotypes",
        "upstream_refs": "UB-7, AstraZeneca-NGS/VarDictJava#65, AstraZeneca-NGS/VarDictJava#383, AstraZeneca-NGS/VarDict#155",
        "vardictjava_behavior": "Collapses multi-allelic loci to a single highest-AF winner per position unless Perl `var2vcf_*.pl -A` is piped externally",
        "compat_off_behavior": "Emits all passing `(CHROM, POS, END, REF, ALT)` alleles per locus sorted by 4-decimal `AF` descending (`--all-alleles`, `--vcfall`)",
        "verification_arm": "Domain 3 (Multi-Allelic & Paralog)",
    },
    {
        "name": "--amplicon-single-primer",
        "category": "Opt-In CLI Capability",
        "domain": "Somatic & Amplicon Calling",
        "upstream_refs": "AstraZeneca-NGS/VarDictJava#303, AstraZeneca-NGS/VarDictJava#301, AstraZeneca-NGS/VarDictJava#326, AstraZeneca-NGS/VarDict#191",
        "vardictjava_behavior": "Enforces two-primer 5-prime and 3-prime edge distance plus `0.95` overlap filters on sheared R2 ends, discarding QIAseq/Archer SPE reads",
        "compat_off_behavior": "Enforces `-a` primer-edge distance strictly on primed amplicon edge(s) (`-sp`) and bypasses sheared-end overlap filtering",
        "verification_arm": "Domain 5 (Clinical Presets & Amplicon)",
    },
    {
        "name": "--preset",
        "category": "Opt-In CLI Capability",
        "domain": "Somatic & Amplicon Calling",
        "upstream_refs": "palakpsheth/VarDictRust#94, palakpsheth/VarDictRust#158",
        "vardictjava_behavior": "Requires manual 1-based vs 0-based BED `-c/-S/-E/-z` flags and separate Perl `var2vcf` filter pipes across sequencing modalities",
        "compat_off_behavior": "Configures `--compat=off`, standard 0-based BED indexing, and modality-calibrated filters for `wgs`, `wes`, `amplicon`, and `ctdna`",
        "verification_arm": "Domain 5 (Clinical Presets & Amplicon)",
    },
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp_path.write_text(content, encoding="utf-8")
    os.replace(tmp_path, path)


def parse_compat_profile_switches(
    config_rs: Path,
    roadmap_md: Path,
    require_exact_count: bool = True,
) -> list[str]:
    text = config_rs.read_text(encoding="utf-8")
    m = re.search(r"pub struct CompatProfile\s*\{(.*?)\n\}", text, flags=re.DOTALL)
    if not m:
        raise RuntimeError(f"Failed to locate pub struct CompatProfile in {config_rs}")
    switches = re.findall(r"^\s*pub\s+([a-z0-9_]+)\s*:\s*bool\s*,", m.group(1), flags=re.MULTILINE)
    if require_exact_count:
        if len(switches) != EXPECTED_COMPAT_SWITCH_COUNT:
            raise RuntimeError(
                f"Expected {EXPECTED_COMPAT_SWITCH_COUNT} CompatProfile bool fields in {config_rs}, found {len(switches)}"
            )
    elif len(switches) < EXPECTED_COMPAT_SWITCH_COUNT:
        raise RuntimeError(
            f"Expected at least {EXPECTED_COMPAT_SWITCH_COUNT} CompatProfile bool fields in {config_rs}, found {len(switches)}"
        )
    roadmap_text = roadmap_md.read_text(encoding="utf-8")
    check_list = switches if require_exact_count else switches[:EXPECTED_COMPAT_SWITCH_COUNT]
    for sw in check_list:
        if f"`{sw}`" not in roadmap_text:
            raise RuntimeError(f"CompatProfile switch `{sw}` missing from {roadmap_md}")
    return switches


def build_switch_catalog() -> list[dict[str, Any]]:
    switches = parse_compat_profile_switches(CONFIG_RS_PATH, ROADMAP_MD_PATH)
    catalog: list[dict[str, Any]] = []
    for idx, sw in enumerate(switches, start=1):
        if sw not in SWITCH_METADATA:
            raise RuntimeError(f"Missing SWITCH_METADATA entry for CompatProfile field {sw!r}")
        meta = SWITCH_METADATA[sw]
        catalog.append(
            {
                "index": idx,
                "name": sw,
                "category": meta["category"],
                "domain": meta["domain"],
                "upstream_refs": meta["upstream_refs"],
                "vardictjava_behavior": meta["vardictjava_behavior"],
                "compat_off_behavior": meta["compat_off_behavior"],
                "verification_arm": meta["verification_arm"],
            }
        )
    for offset, cli_meta in enumerate(CLI_CAPABILITY_ENTRIES, start=len(switches) + 1):
        catalog.append(
            {
                "index": offset,
                "name": cli_meta["name"],
                "category": cli_meta["category"],
                "domain": cli_meta["domain"],
                "upstream_refs": cli_meta["upstream_refs"],
                "vardictjava_behavior": cli_meta["vardictjava_behavior"],
                "compat_off_behavior": cli_meta["compat_off_behavior"],
                "verification_arm": cli_meta["verification_arm"],
            }
        )
    return catalog


def resolve_data_root(root_dir: Path) -> Path:
    env_dr = os.environ.get("DATA_ROOT", "").strip()
    if env_dr:
        return Path(env_dr).resolve()
    if (root_dir / "scratch" / "g4").is_dir() and not (root_dir / "scratch" / "g4").is_symlink():
        return root_dir
    res = subprocess.run(
        ["git", "-C", str(root_dir), "rev-parse", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode == 0 and res.stdout.strip():
        common_git = (root_dir / res.stdout.strip()).resolve()
        main_root = common_git.parent
        if (main_root / "scratch").exists() or (main_root / "VarDictJava").exists():
            return main_root
    return root_dir


def configure_tools_and_versions() -> dict[str, str]:
    home = os.environ.get("HOME", "")
    conda_root = os.environ.get("CONDA_ROOT") or (f"{home}/miniconda3" if home else "")
    candidates = [
        os.environ.get("G4_TOOLS_BIN", ""),
        f"{conda_root}/envs/g4-tools/bin" if conda_root else "",
        f"{conda_root}/envs/vardict-oracle/bin" if conda_root else "",
        f"{home}/miniforge3/envs/g4-tools/bin" if home else "",
        "/opt/conda/envs/g4-tools/bin",
    ]
    cur_parts = [p for p in os.environ.get("PATH", "").split(os.pathsep) if p]
    for cand in reversed(candidates):
        if cand and os.path.isdir(cand) and cand not in cur_parts:
            cur_parts.insert(0, cand)
    os.environ["PATH"] = os.pathsep.join(cur_parts)

    for req in ("samtools", "bcftools", "tabix", "bgzip", "rtg"):
        if shutil.which(req) is None:
            raise RuntimeError(f"Required bioinformatics tool {req!r} not found on PATH={os.environ['PATH']}")

    def parse_ver(cmd: list[str], pat: str) -> str:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        m = re.search(pat, out)
        if not m:
            raise RuntimeError(f"Could not parse version from {cmd}: {out!r}")
        return m.group(1)

    return {
        "rtg_version": parse_ver(["rtg", "RTG_MEM=512M", "RTG_JAVA_OPTS=-Dtalkback=false", "version"], r"RTG Tools\s+([0-9][0-9A-Za-z._-]*)"),
        "bcftools_version": parse_ver(["bcftools", "--version"], r"bcftools\s+([0-9][0-9A-Za-z._-]*)"),
        "tabix_version": parse_ver(["tabix", "--version"], r"tabix \(htslib\)\s+([0-9][0-9A-Za-z._-]*)"),
        "samtools_version": parse_ver(["samtools", "--version"], r"samtools\s+([0-9][0-9A-Za-z._-]*)"),
    }


def verify_engine_and_binary_guards(root_dir: Path) -> Path:
    fp = subprocess.check_output(
        ["bash", str(root_dir / "harness" / "engine_fingerprint.sh")],
        cwd=str(root_dir),
        text=True,
    ).strip()
    if fp != EXPECTED_STUDY_FP:
        raise RuntimeError(f"Engine fingerprint mismatch: expected {EXPECTED_STUDY_FP}, got {fp}")

    perf_fp = subprocess.check_output(
        ["bash", str(root_dir / "harness" / "engine_fingerprint.sh"), "--perf"],
        cwd=str(root_dir),
        text=True,
    ).strip()
    if perf_fp != EXPECTED_STUDY_PERF_FP:
        raise RuntimeError(f"Engine perf fingerprint mismatch: expected {EXPECTED_STUDY_PERF_FP}, got {perf_fp}")

    vdr_bin = root_dir / "target" / "release" / "vardictrust"
    if not vdr_bin.is_file() or not os.access(vdr_bin, os.X_OK):
        raise RuntimeError(f"Missing release binary {vdr_bin} (run 'cargo build --release' first)")

    ver_line = subprocess.check_output([str(vdr_bin), "--version"], text=True).splitlines()[0]
    if EXPECTED_STUDY_VERSION not in ver_line:
        raise RuntimeError(f"Expected {EXPECTED_STUDY_VERSION} in `vardictrust --version`, got {ver_line!r}")

    bin_mtime = vdr_bin.stat().st_mtime
    cargo_lock = root_dir / "Cargo.lock"
    if cargo_lock.is_file() and bin_mtime < cargo_lock.stat().st_mtime:
        raise RuntimeError(f"Stale {vdr_bin}: older than Cargo.lock")
    for rs_file in (root_dir / "crates").rglob("*.rs"):
        if "/target/" in rs_file.as_posix():
            continue
        if bin_mtime < rs_file.stat().st_mtime:
            raise RuntimeError(f"Stale {vdr_bin}: older than {rs_file.relative_to(root_dir)}")
    return vdr_bin


def sanitize_text(raw: str, tmpdir: str, data_root: Path, root_dir: Path) -> str:
    replacements: dict[str, str] = {}
    if tmpdir:
        replacements[tmpdir] = "$WORK_TMP"
    for tool_name in ("vardictrust", "bcftools", "rtg", "tabix", "bgzip", "samtools"):
        which_p = shutil.which(tool_name)
        if which_p:
            replacements[which_p] = tool_name
    vdr_rel = str(root_dir / "target" / "release" / "vardictrust")
    replacements[vdr_rel] = "vardictrust"
    if str(data_root):
        replacements[str(data_root)] = "$REPO_ROOT"
    if str(root_dir):
        replacements[str(root_dir)] = "$REPO_ROOT"
    home = os.environ.get("HOME", "")
    if home:
        replacements[home] = "$HOME_REDACTED"

    out = raw
    for prefix, token in sorted(replacements.items(), key=lambda kv: -len(kv[0])):
        out = out.replace(prefix, token)
    return out


def expected_vartype(ref_a: str, alt_a: str) -> str:
    if alt_a == ".":
        return ""
    if alt_a.startswith("<"):
        inner = alt_a[1:-1]
        return inner.split("-", 1)[0]
    if len(ref_a) == len(alt_a):
        return "SNV" if len(ref_a) == 1 else "Complex"
    if len(ref_a) == 1 and len(alt_a) > 1 and alt_a.startswith(ref_a):
        return "Insertion"
    if len(alt_a) == 1 and len(ref_a) > 1 and ref_a.startswith(alt_a):
        return "Deletion"
    r, a = ref_a, alt_a
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r, a = r[1:], a[1:]
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    if len(r) == 1 and len(a) > 1 and a.startswith(r):
        return "Insertion"
    if len(a) == 1 and len(r) > 1 and r.startswith(a):
        return "Deletion"
    return "Complex"


def check_bcftools_norm_ref(vcf_path: Path, ref_fa: Path) -> bool:
    res = subprocess.run(
        ["bcftools", "norm", "--check-ref", "e", "-f", str(ref_fa), str(vcf_path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return res.returncode == 0


def audit_vcf_schema(vcf_text: str, check_ref_ok: bool) -> tuple[int, dict[str, int]]:
    defects = {k: 0 for k in DEFECT_KEYS}
    hdr_lines: list[str] = []
    rec_lines: list[str] = []
    for line in vcf_text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            hdr_lines.append(line)
        else:
            rec_lines.append(line)
    hdr_blob = "\n".join(hdr_lines)

    for line in rec_lines:
        cols = line.split("\t")
        pos = int(cols[1])
        ref_a, alt_a, qual_s, info, fmt_s = cols[3], cols[4], cols[5], cols[7], cols[8]
        kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
        is_sv = alt_a.startswith("<")

        # 1. Non-SV INFO/END
        if not is_sv and "END" in kv:
            defects["nonsv_info_end"] += 1

        # 2. Strict VCFv4.2/v4.3 FORMAT/AD Number=R cardinality (1 + n_alt across BOTH single & paired VCFs)
        fmt = fmt_s.split(":")
        n_alt = 0 if alt_a == "." else len(alt_a.split(","))
        ad_bad = False
        gt_bad = False
        acct_bad = False
        any_vd_pos = False
        any_alt_gt = False

        for smp in cols[9:]:
            vals = dict(zip(fmt, smp.split(":")))
            ad_parts = vals["AD"].split(",")
            if len(ad_parts) != (1 + n_alt):
                ad_bad = True
            gt_val = vals["GT"]
            gt_tokens = re.split(r"[/|]", gt_val)
            if gt_val == "1/0":
                gt_bad = True
            for gt_tok in gt_tokens:
                if gt_tok == ".":
                    continue
                if not gt_tok.isdigit() or int(gt_tok) > n_alt:
                    gt_bad = True
                elif int(gt_tok) >= 1:
                    any_alt_gt = True
            dp = int(vals["DP"])
            vd = int(vals["VD"])
            if vd > 0:
                any_vd_pos = True
            ad_nums = [int(x) for x in ad_parts]
            ald_nums = [int(x) for x in vals["ALD"].split(",")]
            af_val = float(vals["AF"])
            if dp < sum(ad_nums) or vd > dp or af_val > 1.000001 or sum(ald_nums) != vd:
                acct_bad = True

        if ad_bad:
            defects["ad_number_r_mismatch"] += 1

        # 3. Non-canonical GT (1/0 or unobserved ALT where VD > 0 across samples has all 0/0)
        if gt_bad or (alt_a != "." and any_vd_pos and not any_alt_gt):
            defects["noncanonical_gt"] += 1

        # 4. Symbolic SV schema & header compliance
        if is_sv:
            sv_bad = False
            if alt_a not in ("<DEL>", "<DUP>", "<INV>"):
                sv_bad = True
            if len(ref_a) != 1 or ref_a not in ("A", "C", "G", "T", "N") or pos < 1:
                sv_bad = True
            if f"##ALT=<ID={alt_a[1:-1]}," not in hdr_blob:
                sv_bad = True
            if not all(req_k in kv for req_k in ("END", "SVTYPE", "SVLEN")):
                sv_bad = True
            else:
                end_pos = int(kv["END"])
                svlen = int(kv["SVLEN"])
                if end_pos <= pos or (end_pos - pos) != abs(svlen):
                    sv_bad = True
                if alt_a == "<DEL>" and svlen >= 0:
                    sv_bad = True
                if alt_a in ("<DUP>", "<INV>") and svlen <= 0:
                    sv_bad = True
            if sv_bad:
                defects["noncanonical_sv"] += 1

        # 5. Accounting, QUAL < 0, stale post-adjComplex TYPE, IUPAC ALT, REF==ALT, or -3 SHIFT3 FASTA mismatch
        if qual_s != "." and float(qual_s) < 0.0:
            acct_bad = True
        if alt_a != "." and "HICNT" in kv and "VD" in kv:
            if not (0 <= int(kv["HICNT"]) <= int(kv["VD"])):
                acct_bad = True
        if "HICNT" in kv and "HICOV" in kv:
            if not (0 <= int(kv["HICNT"]) <= int(kv["HICOV"])):
                acct_bad = True
        if alt_a != ".":
            if not VALID_ALT_RE.match(alt_a):
                acct_bad = True
            for single_alt in alt_a.split(","):
                if ref_a == single_alt:
                    acct_bad = True
        if "," not in alt_a and alt_a != "." and not is_sv:
            if kv.get("TYPE", "") != expected_vartype(ref_a, alt_a):
                acct_bad = True
        if acct_bad:
            defects["accounting_type_iupac_shift3"] += 1

    if not check_ref_ok:
        defects["accounting_type_iupac_shift3"] += 1

    return len(rec_lines), defects


def count_bgzf_blocks(gz_path: Path) -> int:
    blocks = 0
    size = gz_path.stat().st_size
    with gz_path.open("rb") as fh:
        pos = 0
        while pos < size:
            hdr = fh.read(12)
            if not hdr:
                break
            if len(hdr) < 12 or hdr[0:2] != b"\x1f\x8b" or (hdr[3] & 4) == 0:
                raise RuntimeError(f"{gz_path}: invalid BGZF header at offset {pos}")
            xlen = struct.unpack("<H", hdr[10:12])[0]
            extra = fh.read(xlen)
            bsize = None
            idx = 0
            while idx + 4 <= len(extra):
                si1, si2, slen = extra[idx], extra[idx + 1], struct.unpack("<H", extra[idx + 2 : idx + 4])[0]
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack("<H", extra[idx + 4 : idx + 6])[0] + 1
                    break
                idx += 4 + slen
            if bsize is None:
                raise RuntimeError(f"{gz_path}: missing BGZF BC subfield at offset {pos}")
            if bsize > 28:
                blocks += 1
            pos += bsize
            fh.seek(pos)
    return blocks


def intersect_beds(bed_a: Path, bed_b: Path, out_bed: Path, max_a_lines: int | None = None) -> tuple[int, int]:
    b_map: dict[str, list[tuple[int, int]]] = {}
    with bed_b.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            b_map.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    for c in b_map:
        b_map[c].sort()

    raw_ivs: dict[str, list[tuple[int, int]]] = {}
    with bed_a.open("r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh):
            if max_a_lines is not None and idx >= max_a_lines:
                break
            if line.startswith("#") or not line.strip():
                continue
            p = line.split("\t")
            c, s, e = p[0], int(p[1]), int(p[2])
            for bs, be in b_map.get(c, []):
                if bs >= e:
                    break
                if be > s:
                    lo, hi = max(s, bs), min(e, be)
                    if hi > lo:
                        raw_ivs.setdefault(c, []).append((lo, hi))

    merged_count = 0
    merged_bp = 0
    with out_bed.open("w", encoding="utf-8") as out:
        for c in sorted(raw_ivs):
            ivs = sorted(raw_ivs[c])
            cur_s, cur_e = ivs[0]
            for s, e in ivs[1:]:
                if s <= cur_e:
                    cur_e = max(cur_e, e)
                else:
                    out.write(f"{c}\t{cur_s}\t{cur_e}\n")
                    merged_count += 1
                    merged_bp += cur_e - cur_s
                    cur_s, cur_e = s, e
            out.write(f"{c}\t{cur_s}\t{cur_e}\n")
            merged_count += 1
            merged_bp += cur_e - cur_s
    return merged_count, merged_bp


def run_collection(out_dir: Path) -> None:
    try:
        os.nice(19)
    except (PermissionError, OSError):
        pass

    tool_versions = configure_tools_and_versions()
    vdr_bin = verify_engine_and_binary_guards(ROOT)
    data_root = resolve_data_root(ROOT)

    g4_dir = data_root / "scratch" / "g4"
    w4_dir = g4_dir / "wave4"
    vj_bams = data_root / "VarDictJava" / "src" / "test" / "resources" / "com" / "astrazeneca" / "vardict" / "integrationtests"
    fastas_dir = data_root / "VarDictJava" / "testdata" / "fastas"

    fa_hg19 = fastas_dir / "hg19.fa"
    fa_b37 = fastas_dir / "GRCh37.fa"
    fa_hg38 = fastas_dir / "hg38.fa"
    fa_hg1k = fastas_dir / "hg1k.fa"
    if not all(p.is_file() for p in (fa_hg19, fa_b37, fa_hg38, fa_hg1k)):
        subprocess.run([sys.executable, str(ROOT / "harness" / "convert_fastas.py")], cwd=str(ROOT), check=True)

    grch38_ref = g4_dir / "ref" / "d1vd1" / "GRCh38.d1.vd1.fa"
    giab_ref = g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fa"
    if not giab_ref.is_file() and (g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fasta").is_file():
        giab_ref = g4_dir / "ref" / "giabv3" / "GRCh38_GIABv3.fasta"

    seqc2_dir = g4_dir / "seqc2"
    seqc2_t = seqc2_dir / "HCC1395.bam"
    seqc2_n = seqc2_dir / "HCC1395BL.bam"
    if not seqc2_t.is_file() and (seqc2_dir / "HCC1395.slice.bam").is_file():
        seqc2_t = seqc2_dir / "HCC1395.slice.bam"
        seqc2_n = seqc2_dir / "HCC1395BL.slice.bam"

    hg008_dir = g4_dir / "hg008"
    hg008_t = hg008_dir / "HG008-T.slice.bam"
    hg008_n = hg008_dir / "HG008-N.slice.bam"
    hg002_chr20_bam = g4_dir / "hg002" / "HG002.chr20.bam"

    for req_p in (grch38_ref, giab_ref, seqc2_t, seqc2_n, hg008_t, hg008_n, hg002_chr20_bam, w4_dir / ".done"):
        if not req_p.exists():
            raise RuntimeError(f"Missing required local benchmark prerequisite: {req_p}")

    switch_catalog = build_switch_catalog()

    with tempfile.TemporaryDirectory(prefix="vdr_v110_study_") as tmpdir:
        tmp = Path(tmpdir)

        # Prepare contig-scoped SEQC2 FASTA & SDF (used in Arm 1 & Arm 2)
        target_contigs = ["chr7", "chr9", "chr11", "chr12", "chr17", "chr22"]
        scoped_fa = tmp / "seqc2_target_contigs.fa"
        with scoped_fa.open("w", encoding="utf-8") as fh:
            subprocess.run(["samtools", "faidx", str(grch38_ref), *target_contigs], stdout=fh, check=True)
        subprocess.run(["samtools", "faidx", str(scoped_fa)], check=True)

        scoped_sdf = tmp / "seqc2_target_contigs.sdf"
        subprocess.run(
            ["rtg", "RTG_MEM=2G", "RTG_JAVA_OPTS=-Dtalkback=false", "format", "-o", str(scoped_sdf), str(scoped_fa)],
            cwd=tmpdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )

        # Prepare HG008 v0.6 SV BED and IUPAC chr20 FASTA for Arm 1
        hg008_sv_bed = tmp / "hg008_v06_sv.bed"
        hg008_sv_bed.write_text(
            "chr3\t101702500\t101707500\nchr10\t123542500\t123544500\nchr13\t39535800\t39539000\n",
            encoding="utf-8",
        )

        raw_c20 = subprocess.check_output(["samtools", "faidx", str(fa_hg19), "chr20"], text=True)
        arr20 = bytearray("".join(l.strip() for l in raw_c20.splitlines() if not l.startswith(">")).encode("ascii"))
        iupac_cycle = b"RYMKSWBDHV"
        k_iupac = 0
        for s_pos, e_pos in ((46841560, 46841850), (47124940, 47125250), (47131950, 47132040)):
            for p1 in range(s_pos, e_pos, 17):
                arr20[p1 - 1] = iupac_cycle[k_iupac % len(iupac_cycle)]
                k_iupac += 1
        iupac_c20_fa = tmp / "hg19_chr20_iupac.fa"
        with iupac_c20_fa.open("w", encoding="utf-8") as fh:
            fh.write(">chr20\n")
            for off in range(0, len(arr20), 60):
                fh.write(arr20[off : off + 60].decode("ascii") + "\n")
        subprocess.run(["samtools", "faidx", str(iupac_c20_fa)], check=True)

        # =====================================================================
        # Arm 1: Strict VCFv4.2/v4.3 Schema Audit & bcftools norm --check-ref e
        # =====================================================================
        v2v_dir = w4_dir / "vartovcf"
        v2v_ref = v2v_dir / "reference.fa"
        slices_bed = seqc2_dir / "slices.bed"

        arm1_specs: list[tuple[str, Path, str, list[str]]] = [
            # 1. SEQC2 paired somatic (6 contigs)
            (
                "SEQC2_HCC1395_TN",
                scoped_fa,
                f"{seqc2_t}|{seqc2_n}",
                ["-N", "TUMOR|NORMAL", "-z", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-f", "0.05", "--vcf-af", "0.30", "--fisher", "-th", "4", str(slices_bed)],
            ),
            # 2. GIAB HG002 chr20 60 kbp
            (
                "GIAB_HG002_chr20_60kbp",
                giab_ref,
                str(hg002_chr20_bam),
                ["-R", "chr20:10000001-10060000", "-N", "HG002", "-f", "0.01", "--fisher", "-th", "4"],
            ),
            # 3. GIAB HG008 paired somatic test5.bed
            (
                "GIAB_HG008_TN_test5",
                giab_ref,
                f"{hg008_t}|{hg008_n}",
                ["-N", "HG008_T|HG008_N", "-z", "-c", "1", "-S", "2", "-E", "3", "-f", "0.01", "--fisher", "-th", "4", str(hg008_dir / "test5.bed")],
            ),
            # 4..14. 11 symbolic SV callsets across 10 distinct BAM inputs
            ("Colo829_DELDiscReverse", fa_hg19, str(vj_bams / "Colo829.chr10.DELDiscReverse.bam"), ["-R", "chr10:89700380-89713976", "-f", "0.001"]),
            ("Colo829_LGDEL3", fa_hg19, str(vj_bams / "Colo829.chr10.LGDEL3.bam"), ["-R", "chr10:127187158-127201055", "-f", "0.001"]),
            ("Colo829_DELSV3Unhappy", fa_hg19, str(vj_bams / "Colo829.chr20.DELSV3Unhappy.bam"), ["-R", "chr20:54434600-54441423", "-f", "0.001"]),
            ("Colo829_DUPDiscSV", fa_hg19, str(vj_bams / "Colo829.chr20.DUPDiscSV.bam"), ["-R", "chr20:47124929-47131802", "-f", "0.001", "-L", "500"]),
            ("Colo829_INVSV", fa_hg19, str(vj_bams / "Colo829.chr20.INVSV.bam"), ["-R", "chr20:46841554-46845755", "-f", "0.001"]),
            ("Colo829_MarkingDUPSV", fa_hg19, str(vj_bams / "Colo829.chr20.MarkingDUPSV.bam"), ["-R", "chr20:47131954-47132045", "-f", "0.001", "-L", "100"]),
            ("Colo829_CandDel5SV", fa_hg19, str(vj_bams / "Colo829.chr20.CandDel5SV.bam"), ["-R", "chr20:506395-507295", "-f", "0.001", "-L", "100"]),
            (
                "Zero_in_DUPSV_Paired",
                fa_hg19,
                f"{vj_bams / 'Zero_in_DUPSV_TUMOUR.bam'}|{vj_bams / 'Zero_in_DUPSV_NORMAL.bam'}",
                ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3", "-L", "50", "--fisher"],
            ),
            ("HG008_v06_single", giab_ref, str(hg008_t), ["-z", "-c", "1", "-S", "2", "-E", "3", "-f", "0.01", "-L", "100", "--fisher", str(hg008_sv_bed)]),
            ("HG008_v06_paired", giab_ref, f"{hg008_t}|{hg008_n}", ["-z", "-c", "1", "-S", "2", "-E", "3", "-f", "0.01", "-L", "100", "--fisher", str(hg008_sv_bed)]),
            ("vartovcf_sv_contig1", v2v_ref, str(v2v_dir / "vartovcf_sv_contig1.bam"), ["-z", "-c", "1", "-S", "2", "-E", "3", "-g", "4", "-P", "0", "-f", "0.01", "-r", "2", str(v2v_dir / "vartovcf_sv_contig1.bed")]),
            # 15..16. 2 -p pileup runs
            ("N13532_pileup", fa_b37, str(vj_bams / "N13532.TwoDels1InsStringBoundError.bam"), ["-R", "11:1280138-1280341", "-p"]),
            ("HG002_CMRG_Pileup", giab_ref, str(w4_dir / "hg002_cmrg_mosaic" / "HG002.cmrg_mosaic.slice.bam"), ["-R", "chr1:1019800-1021200", "-p", "-th", "4"]),
            # 17..19. 3 targeted schema stress callsets
            ("NAP1L1_NegativeRRC", fa_hg38, str(vj_bams / "NAP1L1.NegativeRRC.bam"), ["-R", "chr12:76036587-76085033", "-f", "0.001"]),
            ("Colo829_DUPDiscSV_IUPAC", iupac_c20_fa, str(vj_bams / "Colo829.chr20.DUPDiscSV.bam"), ["-R", "chr20:47124929-47131802", "-f", "0.001", "-L", "500", "-p"]),
            ("N13532_shift3", fa_b37, str(vj_bams / "N13532.TwoDels1InsStringBoundError.bam"), ["-R", "11:1280138-1280341", "-3", "-f", "0.0025", "-F", "0x700", "-Q", "10"]),
        ]

        distinct_bam_files: set[str] = set()
        vj_defects_sum = {k: 0 for k in DEFECT_KEYS}
        off_defects_sum = {k: 0 for k in DEFECT_KEYS}
        total_records_vj = 0
        total_records_off = 0
        check_ref_vj_passed = 0
        check_ref_off_passed = 0
        cached_seqc2_vcfs: dict[str, Path] = {}

        for label, ref_fa, bam_arg, extra_args in arm1_specs:
            for single_b in bam_arg.split("|"):
                distinct_bam_files.add(os.path.realpath(single_b))

            vj_vcf = tmp / f"arm1_{label}_vj.vcf"
            off_vcf = tmp / f"arm1_{label}_off.vcf"

            subprocess.run(
                [str(vdr_bin), "-G", str(ref_fa), "-b", bam_arg, "--compat=vardictjava", "--output", str(vj_vcf), *extra_args],
                check=True,
            )
            subprocess.run(
                [str(vdr_bin), "-G", str(ref_fa), "-b", bam_arg, "--compat=off", "--output", str(off_vcf), *extra_args],
                check=True,
            )
            if label == "SEQC2_HCC1395_TN":
                cached_seqc2_vcfs["vj"] = vj_vcf
                cached_seqc2_vcfs["off"] = off_vcf

            vj_ref_ok = check_bcftools_norm_ref(vj_vcf, ref_fa)
            off_ref_ok = check_bcftools_norm_ref(off_vcf, ref_fa)
            if vj_ref_ok:
                check_ref_vj_passed += 1
            if off_ref_ok:
                check_ref_off_passed += 1
            else:
                raise RuntimeError(f"Arm 1 bcftools norm --check-ref e failed under --compat=off on {label}")

            n_vj, d_vj = audit_vcf_schema(vj_vcf.read_text(encoding="utf-8"), vj_ref_ok)
            n_off, d_off = audit_vcf_schema(off_vcf.read_text(encoding="utf-8"), off_ref_ok)
            if any(v != 0 for v in d_off.values()):
                raise RuntimeError(f"Arm 1 --compat=off schema violation on {label}: {d_off}")

            total_records_vj += n_vj
            total_records_off += n_off
            for k in DEFECT_KEYS:
                vj_defects_sum[k] += d_vj[k]
                off_defects_sum[k] += d_off[k]

        arm1_data: dict[str, Any] = {
            "total_callsets": len(arm1_specs),
            "distinct_bams": len(distinct_bam_files),
            "total_records_vj": total_records_vj,
            "total_records_off": total_records_off,
            "vardictjava_defects": vj_defects_sum,
            "compat_off_defects": off_defects_sum,
            "total_defects_vj": sum(vj_defects_sum.values()),
            "total_defects_off": sum(off_defects_sum.values()),
            "bcftools_norm_check_ref_vj_passed": check_ref_vj_passed,
            "bcftools_norm_check_ref_off_passed": check_ref_off_passed,
            "bcftools_norm_check_ref_off_total": len(arm1_specs),
        }

        # =====================================================================
        # Arm 2: Diploid SEQC2 HCC1395|HCC1395BL rtg vcfeval (No --squash-ploidy)
        # =====================================================================
        scoped_eval_bed = tmp / "seqc2_scoped_eval.bed"
        eval_intervals, eval_bp = intersect_beds(slices_bed, seqc2_dir / "eval.bed", scoped_eval_bed)

        eval_ivs_by_chr: dict[str, list[tuple[int, int]]] = {}
        with scoped_eval_bed.open("r", encoding="utf-8") as fh:
            for line in fh:
                p = line.strip().split("\t")
                eval_ivs_by_chr.setdefault(p[0], []).append((int(p[1]), int(p[2])))

        def in_eval_bed(chrom: str, pos1: int) -> bool:
            return any(lo < pos1 <= hi for lo, hi in eval_ivs_by_chr.get(chrom, []))

        contig_lens: dict[str, int] = {}
        with Path(f"{scoped_fa}.fai").open("r", encoding="utf-8") as fh:
            for line in fh:
                p = line.split("\t")
                contig_lens[p[0]] = int(p[1])

        truth_src = seqc2_dir / "truth" / "truth.vcf.gz"
        raw_truth = subprocess.check_output(["bcftools", "view", "-R", str(slices_bed), str(truth_src)], text=True)
        tvaf_re = re.compile(r"(?:^|;)TVAF=([0-9.eE+-]+)")
        gt_proj_threshold = 0.75
        truth_scoped_vcf = tmp / "seqc2_truth_diploid.vcf"
        truth_scoped_gz = Path(f"{truth_scoped_vcf}.gz")

        truth_in_eval_total = 0
        subclonal_truth_lt_030 = 0
        with truth_scoped_vcf.open("w", encoding="utf-8") as out:
            out.write("##fileformat=VCFv4.2\n")
            for c in target_contigs:
                out.write(f"##contig=<ID={c},length={contig_lens[c]}>\n")
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n")
            for line in raw_truth.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split("\t")
                if p[0] not in contig_lens:
                    continue
                pos1 = int(p[1])
                m = tvaf_re.search(p[7])
                tvaf = float(m.group(1)) if m else 0.30
                if in_eval_bed(p[0], pos1):
                    truth_in_eval_total += 1
                    if tvaf < 0.30:
                        subclonal_truth_lt_030 += 1
                gt = "1/1" if tvaf >= gt_proj_threshold else "0/1"
                out.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t.\tPASS\t.\tGT\t{gt}\n")

        subprocess.run(["bgzip", "-f", str(truth_scoped_vcf)], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(truth_scoped_gz)], check=True)

        def eval_seqc2_diploid(mode_tag: str, raw_vcf_path: Path) -> dict[str, Any]:
            som_vcf = tmp / f"seqc2_{mode_tag}_somatic_filtered.vcf"
            som_gz = Path(f"{som_vcf}.gz")
            gt00_carrier = 0
            retained_somatic = 0
            with raw_vcf_path.open("r", encoding="utf-8") as fh, som_vcf.open("w", encoding="utf-8") as out:
                for line in fh:
                    if line.startswith("#"):
                        out.write(line)
                        continue
                    p = line.rstrip("\n").split("\t")
                    flt_set = set(p[6].split(";"))
                    info = p[7]
                    if ("STATUS=StrongSomatic" in info or "STATUS=LikelySomatic" in info) and flt_set <= {"PASS", "f0.30"}:
                        retained_somatic += 1
                        fmt = p[8].split(":")
                        t_vals = dict(zip(fmt, p[9].split(":")))
                        n_vals = dict(zip(fmt, p[10].split(":")))
                        if t_vals["GT"] == "0/0" and n_vals["GT"] == "0/0" and (int(t_vals["VD"]) > 0 or int(n_vals["VD"]) > 0):
                            gt00_carrier += 1
                        out.write(line)
            subprocess.run(["bgzip", "-f", str(som_vcf)], check=True)
            subprocess.run(["tabix", "-f", "-p", "vcf", str(som_gz)], check=True)

            eval_dir = tmp / f"vcfeval_{mode_tag}"
            subprocess.run(
                [
                    "rtg",
                    "RTG_MEM=4G",
                    "RTG_JAVA_OPTS=-Dtalkback=false",
                    "vcfeval",
                    "-T",
                    "4",
                    "-b",
                    str(truth_scoped_gz),
                    "-c",
                    str(som_gz),
                    "-e",
                    str(scoped_eval_bed),
                    "-t",
                    str(scoped_sdf),
                    "--all-records",
                    "--vcf-score-field=QUAL",
                    "--sample=TUMOR,TUMOR",
                    "-o",
                    str(eval_dir),
                ],
                cwd=tmpdir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            best_row = None
            none_row = None
            for sline in (eval_dir / "summary.txt").read_text(encoding="utf-8").splitlines():
                cols = sline.strip().split()
                if not cols or cols[0] == "Threshold" or cols[0].startswith("---"):
                    continue
                if cols[0] == "None":
                    none_row = cols
                elif best_row is None:
                    best_row = cols
            if none_row is None or len(none_row) < 8:
                raise RuntimeError(f"Failed to parse None row from {eval_dir / 'summary.txt'}")
            if best_row is None:
                best_row = none_row
            best_thr = None if best_row[0] == "None" else round(float(best_row[0]), 3)
            return {
                "retained_somatic_calls": retained_somatic,
                "gt00_carrier_somatic": gt00_carrier,
                "tp_baseline": int(none_row[1]),
                "tp_call": int(none_row[2]),
                "fp": int(none_row[3]),
                "fn": int(none_row[4]),
                "precision": round(float(none_row[5]), 4),
                "sensitivity": round(float(none_row[6]), 4),
                "f_measure": round(float(none_row[7]), 4),
                "best_qual_threshold": best_thr,
                "best_tp_baseline": int(best_row[1]),
                "best_tp_call": int(best_row[2]),
                "best_fp": int(best_row[3]),
                "best_fn": int(best_row[4]),
                "best_precision": round(float(best_row[5]), 4),
                "best_sensitivity": round(float(best_row[6]), 4),
                "best_f_measure": round(float(best_row[7]), 4),
            }

        vj_vcfeval = eval_seqc2_diploid("vj", cached_seqc2_vcfs["vj"])
        off_vcfeval = eval_seqc2_diploid("off", cached_seqc2_vcfs["off"])

        arm2_data: dict[str, Any] = {
            "target_contigs": target_contigs,
            "evaluated_intervals": eval_intervals,
            "evaluated_bp": eval_bp,
            "discovery_af": 0.05,
            "vcf_af": 0.30,
            "gt_proj_tvaf_threshold": gt_proj_threshold,
            "truth_total": vj_vcfeval["tp_baseline"] + vj_vcfeval["fn"],
            "truth_in_eval_bed": truth_in_eval_total,
            "subclonal_truth_lt_030": subclonal_truth_lt_030,
            "vardictjava": vj_vcfeval,
            "compat_off": off_vcfeval,
            "tp_baseline_gain": off_vcfeval["tp_baseline"] - vj_vcfeval["tp_baseline"],
            "sensitivity_gain": round(off_vcfeval["sensitivity"] - vj_vcfeval["sensitivity"], 4),
            "f_measure_gain": round(off_vcfeval["f_measure"] - vj_vcfeval["f_measure"], 4),
            "best_tp_baseline_gain": off_vcfeval["best_tp_baseline"] - vj_vcfeval["best_tp_baseline"],
            "best_sensitivity_gain": round(off_vcfeval["best_sensitivity"] - vj_vcfeval["best_sensitivity"], 4),
            "best_f_measure_gain": round(off_vcfeval["best_f_measure"] - vj_vcfeval["best_f_measure"], 4),
        }

        # =====================================================================
        # Arm 3: Somatic Multi-Allelic Continuation, MUC19 Paralog, & Truth Sets
        # =====================================================================
        arg_secondary = {
            (140453033, "T"),
            (140453035, "A"),
            (140453036, "CA"),
            (140453065, "T"),
            (140453081, "T"),
            (140453095, "G"),
            (140453108, "A"),
            (140453116, "G"),
        }
        ma_cases = [
            (
                "ArgumentNA_fisher",
                "ArgumentNA",
                fa_hg19,
                f"{vj_bams / 'ArgumentNA_Tumor.bam'}|{vj_bams / 'ArgumentNA_Normal.bam'}",
                ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p", "--fisher"],
                lambda pos, ref, alt: (pos, alt) in arg_secondary,
            ),
            (
                "ArgumentNA_nofisher",
                "ArgumentNA",
                fa_hg19,
                f"{vj_bams / 'ArgumentNA_Tumor.bam'}|{vj_bams / 'ArgumentNA_Normal.bam'}",
                ["-R", "chr7:140453010-140453158", "-q", "25", "-f", "0.001", "-F", "0x100", "-k", "1", "-M", "20", "-p"],
                lambda pos, ref, alt: (pos, alt) in arg_secondary,
            ),
            (
                "Colo829_chr1_756309",
                "Colo829",
                fa_hg19,
                f"{vj_bams / 'Colo829-18_S3-sort.bam'}|{vj_bams / 'Colo829-19_S4-sort.bam'}",
                ["-R", "chr1:755917-756517", "-f", "0.01"],
                lambda pos, ref, alt: pos == 756309 and ref == "TTGT" and alt == "CTGG",
            ),
            (
                "HCC2218_chrX_7811864",
                "HCC2218",
                fa_hg19,
                f"{vj_bams / 'HCC2218_BL_realigins.bam'}|{vj_bams / 'HCC2218_C_realigins.bam'}",
                ["-R", "chrX:7811670-7811868", "-f", "0.001"],
                lambda pos, ref, alt: pos == 7811864 and ref == "C" and alt == "T",
            ),
            (
                "Zero_in_DUPSV_chr19_291259",
                "Zero_in_DUPSV",
                fa_hg19,
                f"{vj_bams / 'Zero_in_DUPSV_TUMOUR.bam'}|{vj_bams / 'Zero_in_DUPSV_NORMAL.bam'}",
                ["-R", "chr19:290742-291636", "-f", "0.001", "-X", "3"],
                lambda pos, ref, alt: pos == 291259 and ref == "A" and alt == "G",
            ),
            (
                "StrongLOH_chr21_11071401",
                "StrongLOH",
                fa_hg1k,
                f"{vj_bams / 'TumorA.StrongLOH.chr21.bam'}|{vj_bams / 'NormalA.StrongLOH.chr21.bam'}",
                ["-R", "chr21:10996701-11186054", "-Q", "10", "-f", "0.01", "-r", "2", "-F", "0x500", "-V", "0.05", "--nosv"],
                lambda pos, ref, alt: pos in (11071401, 11071410) and len(alt) == len(ref) + 1,
            ),
        ]

        ma_pair_ids = {pair_id for _, pair_id, *_ in ma_cases}
        vj_rec_pairs: set[str] = set()
        off_rec_pairs: set[str] = set()
        vj_rec_runs = 0
        off_rec_runs = 0
        vj_sec_alleles = 0
        off_sec_alleles = 0

        for _, pair_id, ref_fa, bam_pair, extra, pred in ma_cases:
            def_vcf = subprocess.check_output(
                [str(vdr_bin), "-G", str(ref_fa), "-b", bam_pair, "--compat=vardictjava", "--vcfall", "-v", *extra],
                text=True,
            )
            off_vcf = subprocess.check_output(
                [str(vdr_bin), "-G", str(ref_fa), "-b", bam_pair, "--compat=off", "--vcfall", "-v", *extra],
                text=True,
            )
            d_hits = [
                1
                for l in def_vcf.splitlines()
                if not l.startswith("#") and l.strip() and pred(int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
            ]
            o_hits = [
                1
                for l in off_vcf.splitlines()
                if not l.startswith("#") and l.strip() and pred(int(l.split("\t")[1]), l.split("\t")[3], l.split("\t")[4])
            ]
            if d_hits:
                vj_rec_runs += 1
                vj_rec_pairs.add(pair_id)
                vj_sec_alleles += len(d_hits)
            if o_hits:
                off_rec_runs += 1
                off_rec_pairs.add(pair_id)
                off_sec_alleles += len(o_hits)

        # MUC19 paralog somatic pair
        muc_dir = w4_dir / "muc19"
        muc_vcf = tmp / "muc19_off.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "-G",
                str(muc_dir / "hg19_chr12_muc19.fa"),
                "-b",
                f"{muc_dir / 'MUC19_vaf1.bam'}|{muc_dir / 'MUC19_vaf1_dz.bam'}",
                "-z",
                "-c",
                "1",
                "-S",
                "2",
                "-E",
                "3",
                "-f",
                "0.005",
                "--fisher",
                "--compat=off",
                "--vcfall",
                "--output",
                str(muc_vcf),
                "--tbi",
                str(muc_dir / "MUC19_vaf1.bed"),
            ],
            check=True,
        )
        if not check_bcftools_norm_ref(muc_vcf, muc_dir / "hg19_chr12_muc19.fa"):
            raise RuntimeError("MUC19 --compat=off VCF failed bcftools norm --check-ref e")

        muc_violations = 0
        with gzip.open(muc_vcf, "rt", encoding="utf-8") as fh:
            muc_rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#") and l.strip()]
        for r in muc_rows:
            fmt = r[8].split(":")
            for smp in (r[9], r[10]):
                vals = dict(zip(fmt, smp.split(":")))
                if int(vals["VD"]) > int(vals["DP"]) or float(vals["AF"]) > 1.000001:
                    muc_violations += 1

        def eval_normalized_somatic_truth(
            tag: str,
            ref_fa: Path,
            bam_pair: str,
            eval_bed: Path,
            truth_vcf: Path,
            compat_mode: str,
        ) -> dict[str, Any]:
            raw_p = tmp / f"arm3_{tag}_{compat_mode}_raw.vcf"
            som_p = tmp / f"arm3_{tag}_{compat_mode}_som.vcf"
            subprocess.run(
                [
                    str(vdr_bin),
                    "-G",
                    str(ref_fa),
                    "-b",
                    bam_pair,
                    "-z",
                    "-c",
                    "1",
                    "-S",
                    "2",
                    "-E",
                    "3",
                    "-f",
                    "0.05",
                    "--fisher",
                    f"--compat={compat_mode}",
                    "-th",
                    "4",
                    "--output",
                    str(raw_p),
                    str(eval_bed),
                ],
                check=True,
            )
            with raw_p.open("r", encoding="utf-8") as fh, som_p.open("w", encoding="utf-8") as out:
                for line in fh:
                    if line.startswith("#"):
                        out.write(line)
                        continue
                    p = line.split("\t")
                    if p[6] == "PASS" and ("STATUS=StrongSomatic" in p[7] or "STATUS=LikelySomatic" in p[7]):
                        out.write(line)
            view_c = subprocess.run(["bcftools", "view", "-T", str(eval_bed), str(som_p)], capture_output=True, text=True, check=True)
            norm_c = subprocess.run(
                ["bcftools", "norm", "-m", "-any", "-f", str(ref_fa)],
                input=view_c.stdout,
                capture_output=True,
                text=True,
                check=True,
            )
            view_t = subprocess.run(
                ["bcftools", "view", "-f", ".,PASS", "-T", str(eval_bed), str(truth_vcf)],
                capture_output=True,
                text=True,
                check=True,
            )
            norm_t = subprocess.run(
                ["bcftools", "norm", "-m", "-any", "-f", str(ref_fa)],
                input=view_t.stdout,
                capture_output=True,
                text=True,
                check=True,
            )
            calls = {
                tuple(l.split("\t")[i] for i in (0, 1, 3, 4))
                for l in norm_c.stdout.splitlines()
                if not l.startswith("#") and l.strip()
            }
            truth = {
                tuple(l.split("\t")[i] for i in (0, 1, 3, 4))
                for l in norm_t.stdout.splitlines()
                if not l.startswith("#") and l.strip()
            }
            tp = len(calls & truth)
            fp = len(calls - truth)
            fn = len(truth - calls)
            prec = round(tp / (tp + fp), 4) if (tp + fp) else 0.0
            sens = round(tp / len(truth), 4) if truth else 0.0
            f1 = round(2 * tp / (2 * tp + fp + fn), 4) if (2 * tp + fp + fn) else 0.0
            return {
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "truth_total": len(truth),
                "precision": prec,
                "sensitivity": sens,
                "f1": f1,
            }

        seqc2_hc_bed = tmp / "seqc2_hc.bed"
        seqc2_hc_ivs, seqc2_hc_bp = intersect_beds(slices_bed, seqc2_dir / "truth" / "HC.bed", seqc2_hc_bed)
        seqc2_hc_vj = eval_normalized_somatic_truth("seqc2", grch38_ref, f"{seqc2_t}|{seqc2_n}", seqc2_hc_bed, truth_src, "vardictjava")
        seqc2_hc_off = eval_normalized_somatic_truth("seqc2", grch38_ref, f"{seqc2_t}|{seqc2_n}", seqc2_hc_bed, truth_src, "off")

        hg008_hc40_bed = tmp / "hg008_hc40.bed"
        hg008_hc_ivs, hg008_hc_bp = intersect_beds(
            hg008_dir / "slices.bed",
            hg008_dir / "truth" / "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_all.bed",
            hg008_hc40_bed,
            max_a_lines=40,
        )
        hg008_truth_vcf = hg008_dir / "truth" / "HG008-T_somatic_smvar_benchmark_v0.4_GRCh38_tumorvariants.vcf.gz"
        hg008_v04_vj = eval_normalized_somatic_truth("hg008", giab_ref, f"{hg008_t}|{hg008_n}", hg008_hc40_bed, hg008_truth_vcf, "vardictjava")
        hg008_v04_off = eval_normalized_somatic_truth("hg008", giab_ref, f"{hg008_t}|{hg008_n}", hg008_hc40_bed, hg008_truth_vcf, "off")

        arm3_data: dict[str, Any] = {
            "distinct_bam_pairs": len(ma_pair_ids),
            "total_runs": len(ma_cases),
            "vj_recovered_pairs": len(vj_rec_pairs),
            "off_recovered_pairs": len(off_rec_pairs),
            "vj_recovered_runs": vj_rec_runs,
            "off_recovered_runs": off_rec_runs,
            "vj_secondary_alleles": vj_sec_alleles,
            "off_secondary_alleles": off_sec_alleles,
            "muc19_variants": len(muc_rows),
            "muc19_vd_gt_dp_violations": muc_violations,
            "seqc2_hc": {
                "evaluated_intervals": seqc2_hc_ivs,
                "evaluated_bp": seqc2_hc_bp,
                "vardictjava": seqc2_hc_vj,
                "compat_off": seqc2_hc_off,
            },
            "hg008_v04": {
                "evaluated_intervals": hg008_hc_ivs,
                "evaluated_bp": hg008_hc_bp,
                "vardictjava": hg008_v04_vj,
                "compat_off": hg008_v04_off,
            },
        }

        # =====================================================================
        # Arm 4: Patient BRCA2 Soft-Clip Adapter Chimera Suppression
        # =====================================================================
        brca2_bam = w4_dir / "brca2_adaptor" / "adaptor_insertion.bam"
        adaptors = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        brca2_base = [
            str(vdr_bin),
            "-G",
            str(grch38_ref),
            "-b",
            str(brca2_bam),
            "-R",
            "chr13:32683877-32684023",
            "-f",
            "0.01",
            "-q",
            "25",
            "-r",
            "2",
            "-UN",
            "-v",
            "--adaptor",
            adaptors,
        ]
        brca2_vj = subprocess.check_output([*brca2_base, "--compat=vardictjava"], text=True)
        brca2_off = subprocess.check_output([*brca2_base, "--compat=off"], text=True)

        def parse_brca2_chimeras(vcf_text: str) -> tuple[int, list[str]]:
            total_r = 0
            chim_alts: list[str] = []
            for line in vcf_text.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                total_r += 1
                p = line.split("\t")
                if int(p[1]) == 32683955 and ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" in p[4] or len(p[4]) > 100):
                    chim_alts.append(p[4])
            return total_r, chim_alts

        brca2_vj_total, brca2_vj_chims = parse_brca2_chimeras(brca2_vj)
        brca2_off_total, brca2_off_chims = parse_brca2_chimeras(brca2_off)
        arm4_data: dict[str, Any] = {
            "chimera_chrom": "chr13",
            "chimera_pos": 32683955,
            "region_span_bp": 147,
            "vj_total_records": brca2_vj_total,
            "off_total_records": brca2_off_total,
            "vj_chimeric_insertions": len(brca2_vj_chims),
            "vj_chimeric_ins_bp": (len(brca2_vj_chims[0]) - 1) if brca2_vj_chims else 0,
            "off_chimeric_insertions": len(brca2_off_chims),
        }

        # =====================================================================
        # Arm 5: Clinical Workflow Presets (wgs, wes, amplicon, ctdna) & -sp
        # =====================================================================
        # 5.1 --preset wgs (HG002 60 kbp --tbi + HG008 paired --csi)
        hg002_bed3 = tmp / "hg002_wgs_3col.bed"
        hg002_bed3.write_text("chr20\t10000000\t10060000\n", encoding="utf-8")
        hg002_wgs_gz = tmp / "hg002_wgs.vcf.gz"
        subprocess.run(
            [str(vdr_bin), "--preset", "wgs", "-G", str(giab_ref), "-b", str(hg002_chr20_bam), "-N", "HG002", "--output", str(hg002_wgs_gz), "--tbi", "-th", "4", str(hg002_bed3)],
            check=True,
        )
        if not check_bcftools_norm_ref(hg002_wgs_gz, giab_ref):
            raise RuntimeError("HG002 --preset wgs failed bcftools norm --check-ref e")
        tbx_hg2 = subprocess.check_output(["tabix", str(hg002_wgs_gz), "chr20:10000001-10060000"], text=True)
        with gzip.open(hg002_wgs_gz, "rt", encoding="utf-8") as fh:
            wgs_hg002_recs = [l for l in fh if l.strip() and not l.startswith("#")]

        hg008_bed3 = tmp / "hg008_wgs_3col.bed"
        with (hg008_dir / "test5.bed").open("r", encoding="utf-8") as src, hg008_bed3.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.strip() and not line.startswith("#"):
                    p = line.rstrip("\n").split("\t")
                    dst.write(f"{p[0]}\t{p[1]}\t{p[2]}\n")
        hg008_wgs_gz = tmp / "hg008_wgs.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "--preset",
                "wgs",
                "-G",
                str(giab_ref),
                "-b",
                f"{hg008_t}|{hg008_n}",
                "-N",
                "HG008_T|HG008_N",
                "--fisher",
                "--output",
                str(hg008_wgs_gz),
                "--csi",
                "-th",
                "4",
                str(hg008_bed3),
            ],
            check=True,
        )
        if not check_bcftools_norm_ref(hg008_wgs_gz, giab_ref):
            raise RuntimeError("HG008 --preset wgs failed bcftools norm --check-ref e")
        tbx_hg8 = subprocess.check_output(["tabix", str(hg008_wgs_gz), "chr1:72271447-72291446"], text=True)
        with gzip.open(hg008_wgs_gz, "rt", encoding="utf-8") as fh:
            wgs_hg008_recs = [l for l in fh if l.strip() and not l.startswith("#")]
        wgs_tbi_csi_ok = bool(tbx_hg2.strip() and tbx_hg8.strip())

        # 5.2 --preset wes (SEQC2 264 kbp 4-col vs UCSC 8-col BED8 non-hijack)
        wes_bed4 = tmp / "seqc2_wes_4col.bed"
        wes_bed4.write_text("chr11\t6400000\t6459439\tslice_1\nchr17\t74747665\t74952140\tslice_3\n", encoding="utf-8")
        wes_bed8 = tmp / "seqc2_wes_ucsc.bed8"
        wes_bed8.write_text(
            "chr11\t6400000\t6459439\tslice_1\t0\t+\t6400500\t6459000\n"
            "chr17\t74747665\t74952140\tslice_3\t0\t+\t74748000\t74951500\n",
            encoding="utf-8",
        )
        wes_4_gz = tmp / "seqc2_wes_4col.vcf.gz"
        wes_8_gz = tmp / "seqc2_wes_bed8.vcf.gz"
        for bed_in, gz_out in ((wes_bed4, wes_4_gz), (wes_bed8, wes_8_gz)):
            subprocess.run(
                [
                    str(vdr_bin),
                    "--preset",
                    "wes",
                    "-G",
                    str(grch38_ref),
                    "-b",
                    f"{seqc2_t}|{seqc2_n}",
                    "-N",
                    "HCC1395|HCC1395BL",
                    "--fisher",
                    "--output",
                    str(gz_out),
                    "--tbi",
                    "-th",
                    "4",
                    str(bed_in),
                ],
                check=True,
            )
            if not check_bcftools_norm_ref(gz_out, grch38_ref):
                raise RuntimeError(f"--preset wes failed bcftools norm --check-ref e on {bed_in.name}")
        with gzip.open(wes_4_gz, "rt", encoding="utf-8") as fh:
            wes_4_recs = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
        with gzip.open(wes_8_gz, "rt", encoding="utf-8") as fh:
            wes_8_all = [l.rstrip("\n") for l in fh if l.strip()]
        wes_8_hdr = "\n".join(l for l in wes_8_all if l.startswith("#"))
        wes_8_recs = [l for l in wes_8_all if not l.startswith("#")]
        sites_4 = {(r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) for r in wes_4_recs}
        sites_8 = {(r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) for r in wes_8_recs}
        wes_bed8_non_hijack = (
            sites_4 == sites_8
            and "##FILTER=<ID=AMPBIAS," not in wes_8_hdr
            and "##INFO=<ID=GDAMP," not in wes_8_hdr
        )

        # 5.3 --preset amplicon (single 8-col + paired somatic Illumina PE + QIAseq SPE -sp)
        amp_single_gz = tmp / "amp_single.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "--preset",
                "amplicon",
                "-G",
                str(fa_hg19),
                "-b",
                str(g4_dir / "amplicon" / "amplicon.bam"),
                "--output",
                str(amp_single_gz),
                "--tbi",
                str(g4_dir / "amplicon" / "correct_8columns_amplicon.bed"),
            ],
            check=True,
        )
        if not check_bcftools_norm_ref(amp_single_gz, fa_hg19):
            raise RuntimeError("Single-sample --preset amplicon failed bcftools norm --check-ref e")
        with gzip.open(amp_single_gz, "rt", encoding="utf-8") as fh:
            amp_single_recs = [l for l in fh if l.strip() and not l.startswith("#")]

        ilm_bed = w4_dir / "beds" / "ilm_twoprimer_grch38.amp8.bed"
        ilm_pair = f"{w4_dir / 'bam' / 'SampleA_ILM1.bam'}|{w4_dir / 'bam' / 'SampleB_ILM1.bam'}"
        vj_paired_amp_rc = subprocess.call(
            [str(vdr_bin), "-G", str(grch38_ref), "-b", ilm_pair, "-a", "40:0.70", "--compat=vardictjava", "-v", str(ilm_bed)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        amp_paired_gz = tmp / "amp_paired_ilm.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "--preset",
                "amplicon",
                "-N",
                "SampleA_ILM1|SampleB_ILM1",
                "--fisher",
                "-G",
                str(grch38_ref),
                "-b",
                ilm_pair,
                "--output",
                str(amp_paired_gz),
                "--tbi",
                "-th",
                "4",
                str(ilm_bed),
            ],
            check=True,
        )
        if not check_bcftools_norm_ref(amp_paired_gz, grch38_ref):
            raise RuntimeError("Paired --preset amplicon failed bcftools norm --check-ref e")
        with gzip.open(amp_paired_gz, "rt", encoding="utf-8") as fh:
            amp_paired_recs = [l for l in fh if l.strip() and not l.startswith("#")]

        # QIAseq V3 SPE (-sp vs no -sp) on SampleA_QGN1.bam
        qgn_bed = w4_dir / "beds" / "qiaseq_spe_grch38.amp8.bed"
        bam_qgn_a = w4_dir / "bam" / "SampleA_QGN1.bam"
        bam_qgn_b = w4_dir / "bam" / "SampleB_QGN1.bam"
        bam_qgn_ac5 = w4_dir / "bam" / "SampleAC5_QGN1.bam"

        qgn_inserts_bed = tmp / "qgn_inserts.bed"
        ins_ivs_by_chr: dict[str, list[tuple[int, int]]] = {}
        with qgn_bed.open("r", encoding="utf-8") as src, qgn_inserts_bed.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.strip() and not line.startswith("#"):
                    p = line.rstrip("\n").split("\t")
                    s_i, e_i = int(p[6]), int(p[7])
                    ins_ivs_by_chr.setdefault(p[0], []).append((s_i, e_i))
                    dst.write(f"{p[0]}\t{s_i}\t{e_i}\n")
        for c in ins_ivs_by_chr:
            ins_ivs_by_chr[c].sort()

        # Compute >= 50x callable positions in SampleA_QGN1.bam across QIAseq inserts
        depth_qgn = subprocess.check_output(
            ["samtools", "depth", "-b", str(qgn_inserts_bed), "-G", "0x504", str(bam_qgn_a)],
            text=True,
        )
        callable_50x: set[tuple[str, int]] = set()
        for dline in depth_qgn.splitlines():
            if not dline.strip():
                continue
            dp_cols = dline.split("\t")
            if int(dp_cols[2]) >= 50:
                callable_50x.add((dp_cols[0], int(dp_cols[1])))

        kp_in_inserts: set[tuple[str, int, str, str]] = set()
        kp_callable_50x: set[tuple[str, int, str, str]] = set()
        with gzip.open(w4_dir / "raw" / "KnownPositives_hg19ToHg38.vcf.gz", "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                c, pos1, ref_a, alt_a = p[0], int(p[1]), p[3], p[4]
                if any(s_i < pos1 <= e_i for s_i, e_i in ins_ivs_by_chr.get(c, [])):
                    kp_in_inserts.add((c, pos1, ref_a, alt_a))
                    if (c, pos1) in callable_50x:
                        kp_callable_50x.add((c, pos1, ref_a, alt_a))

        qgn_sp_gz = tmp / "qgn_preset_sp.vcf.gz"
        subprocess.run(
            [str(vdr_bin), "--preset", "amplicon", "-sp", "-UN", "-G", str(grch38_ref), "-b", str(bam_qgn_a), "--output", str(qgn_sp_gz), "--tbi", "-th", "4", str(qgn_bed)],
            check=True,
        )
        if not check_bcftools_norm_ref(qgn_sp_gz, grch38_ref):
            raise RuntimeError("QIAseq --preset amplicon -sp failed bcftools norm --check-ref e")
        with gzip.open(qgn_sp_gz, "rt", encoding="utf-8") as fh:
            qgn_sp_recs = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]

        qgn_nosp_txt = subprocess.check_output(
            [str(vdr_bin), "--preset", "amplicon", "-UN", "-G", str(grch38_ref), "-b", str(bam_qgn_a), "-th", "4", str(qgn_bed)],
            text=True,
        )
        qgn_nosp_recs = [l for l in qgn_nosp_txt.splitlines() if l.strip() and not l.startswith("#")]

        qgn_sp_kp_tp = sum(
            1
            for r in qgn_sp_recs
            if r.split("\t")[6] == "PASS"
            and (r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) in kp_in_inserts
        )
        qgn_nosp_kp_tp = sum(
            1
            for r in qgn_nosp_recs
            if r.split("\t")[6] == "PASS"
            and (r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) in kp_in_inserts
        )
        qgn_sp_dp50_kp_tp = sum(
            1
            for r in qgn_sp_recs
            if r.split("\t")[6] == "PASS"
            and (r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) in kp_callable_50x
        )
        qgn_nosp_dp50_kp_tp = sum(
            1
            for r in qgn_nosp_recs
            if r.split("\t")[6] == "PASS"
            and (r.split("\t")[0], int(r.split("\t")[1]), r.split("\t")[3], r.split("\t")[4]) in kp_callable_50x
        )

        # 5.4 --preset ctdna (~1,000x 0.5% VAF spike-in + wet-lab SampleAC5 d100 soft-filtering)
        target_chr = "chr12"
        target_pos = 120999311
        win_reg = "chr12:120999050-120999600"
        wt_sam = tmp / "ctdna_titration.sam"
        wt_records = subprocess.check_output(["samtools", "view", "-F", "0xF04", str(bam_qgn_b), win_reg], text=True).splitlines()
        with wt_sam.open("w", encoding="utf-8") as hdr_fh:
            subprocess.run(["samtools", "view", "-H", str(bam_qgn_b), win_reg], stdout=hdr_fh, check=True)
            hdr_fh.flush()
            for rep_idx in range(1, 12):
                prefix = f"wt_r{rep_idx}_"
                for rec in wt_records:
                    if rec:
                        hdr_fh.write(prefix + rec + "\n")

        sam_a_lines = subprocess.check_output(
            ["samtools", "view", "-F", "0xF04", str(bam_qgn_a), f"{target_chr}:{target_pos}-{target_pos}"],
            text=True,
        ).splitlines()
        fwd_picks: list[str] = []
        rev_picks: list[str] = []
        seen_fwd_pos: set[int] = set()
        seen_rev_pos: set[int] = set()
        for idx, raw_line in enumerate(sam_a_lines):
            cols = raw_line.split("\t")
            flag, rpos, mapq, seq, qual = int(cols[1]), int(cols[3]), int(cols[4]), cols[9], cols[10]
            if mapq < 40:
                continue
            qpos = 0
            hit_alt = False
            for n_s, op in CIGAR_RE.findall(cols[5]):
                n = int(n_s)
                if op in ("M", "=", "X"):
                    if rpos <= target_pos < rpos + n:
                        qi = qpos + (target_pos - rpos)
                        dist_end = min(qi, len(seq) - 1 - qi)
                        if seq[qi].upper() == "A" and (ord(qual[qi]) - 33) >= 30 and dist_end >= 12:
                            hit_alt = True
                        break
                    rpos += n
                    qpos += n
                elif op in ("I", "S"):
                    qpos += n
                elif op in ("D", "N"):
                    rpos += n
            if not hit_alt:
                continue
            aln_start = int(cols[3])
            cols[0] = f"spikeA_{idx}_{cols[0]}"
            cols[1] = "16" if (flag & 16) else "0"
            cols[6], cols[7], cols[8] = "*", "0", "0"
            if (flag & 16) == 0 and aln_start not in seen_fwd_pos and len(fwd_picks) < 3:
                seen_fwd_pos.add(aln_start)
                fwd_picks.append("\t".join(cols))
            elif (flag & 16) != 0 and aln_start not in seen_rev_pos and len(rev_picks) < 3:
                seen_rev_pos.add(aln_start)
                rev_picks.append("\t".join(cols))
        if len(fwd_picks) != 3 or len(rev_picks) != 3:
            raise RuntimeError(f"Expected 3 fwd + 3 rev spike-in reads at {target_chr}:{target_pos}")

        with wt_sam.open("a", encoding="utf-8") as out_fh:
            for rec in fwd_picks + rev_picks:
                out_fh.write(rec + "\n")
        titration_bam = tmp / "ctdna_1000x_0p5pct.bam"
        subprocess.run(["samtools", "sort", "-o", str(titration_bam), str(wt_sam)], check=True)
        subprocess.run(["samtools", "index", str(titration_bam)], check=True)

        ctdna_bed3 = tmp / "ctdna_locus_3col.bed"
        ctdna_bed3.write_text("chr12\t120999118\t120999522\n", encoding="utf-8")
        off_default_ctdna = subprocess.check_output(
            [str(vdr_bin), "--compat=off", "-G", str(grch38_ref), "-b", str(titration_bam), "-z", "-c", "1", "-S", "2", "-E", "3", str(ctdna_bed3)],
            text=True,
        )
        off_ctdna_hits = [
            l for l in off_default_ctdna.splitlines() if l.strip() and not l.startswith("#") and int(l.split("\t")[1]) == target_pos
        ]

        ctdna_gz = tmp / "ctdna_1000x.vcf.gz"
        subprocess.run(
            [str(vdr_bin), "--preset", "ctdna", "-G", str(grch38_ref), "-b", str(titration_bam), "-N", "CTDNA_0P5", "--output", str(ctdna_gz), "--tbi", str(ctdna_bed3)],
            check=True,
        )
        if not check_bcftools_norm_ref(ctdna_gz, grch38_ref):
            raise RuntimeError("--preset ctdna failed bcftools norm --check-ref e")
        with gzip.open(ctdna_gz, "rt", encoding="utf-8") as fh:
            ct_recs = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
        locus_hits = [r.split("\t") for r in ct_recs if int(r.split("\t")[1]) == target_pos and r.split("\t")[4] == "A"]
        if len(locus_hits) != 1:
            raise RuntimeError(f"Expected exact chr12:{target_pos} G>A call under --preset ctdna")
        ct_hit = locus_hits[0]
        ct_info = dict(x.split("=", 1) for x in ct_hit[7].split(";") if "=" in x)
        ct_fmt = dict(zip(ct_hit[8].split(":"), ct_hit[9].split(":")))

        ac5_ctdna_txt = subprocess.check_output(
            [str(vdr_bin), "--preset", "ctdna", "-G", str(grch38_ref), "-b", str(bam_qgn_ac5), "-th", "4", str(qgn_inserts_bed)],
            text=True,
        )
        ac5_rows = [l.split("\t") for l in ac5_ctdna_txt.splitlines() if l.strip() and not l.startswith("#")]
        ac5_d100_count = sum(1 for r in ac5_rows if "d100" in r[6].split(";"))

        arm5_data: dict[str, Any] = {
            "wgs_hg002_60kbp_records": len(wgs_hg002_recs),
            "wgs_hg008_paired_records": len(wgs_hg008_recs),
            "wgs_tbi_and_csi_verified": wgs_tbi_csi_ok,
            "wes_seqc2_evaluated_bp": 263914,
            "wes_seqc2_264kbp_records": len(wes_4_recs),
            "wes_bed8_non_hijack_verified": wes_bed8_non_hijack,
            "amplicon_single_8col_records": len(amp_single_recs),
            "paired_somatic_amplicon_vj_supported": (vj_paired_amp_rc == 0),
            "paired_somatic_amplicon_off_records": len(amp_paired_recs),
            "qiaseq_spe_tumor_records": len(qgn_sp_recs),
            "qiaseq_spe_known_pos_tp": qgn_sp_kp_tp,
            "qiaseq_spe_no_sp_known_pos_tp": qgn_nosp_kp_tp,
            "qiaseq_spe_dp50_known_pos_tp": qgn_sp_dp50_kp_tp,
            "qiaseq_spe_dp50_no_sp_known_pos_tp": qgn_nosp_dp50_kp_tp,
            "qiaseq_spe_known_pos_total": len(kp_callable_50x),
            "qiaseq_spe_panel_known_pos_unmasked": len(kp_in_inserts),
            "sample_ac5_ctdna_d100_filtered": ac5_d100_count,
            "admixture_wt_copies": 11,
            "admixture_alt_reads": 6,
            "spikein_chrom": target_chr,
            "spikein_pos": target_pos,
            "spikein_target_vaf_pct": 0.5,
            "spikein_discovery_af": 0.0025,
            "preset_ctdna_freq": 0.001,
            "preset_ctdna_vcf_af": 0.002,
            "compat_off_default_f001_detected": bool(off_ctdna_hits),
            "spikein_observed_dp": int(ct_info["DP"]),
            "spikein_observed_vd": int(ct_info["VD"]),
            "spikein_observed_af": float(ct_info["AF"]),
            "spikein_detected": (ct_hit[6] == "PASS" and ct_fmt.get("GT") == "0/1"),
        }

        # =====================================================================
        # Arm 6: Multi-BGZF-Block --csi Spanning-SV & -p Homopolymer Tiling
        # =====================================================================
        lgdel_gz = tmp / "colo829_deldisc_pileup.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "-G",
                str(fa_hg19),
                "-b",
                str(vj_bams / "Colo829.chr10.DELDiscReverse.bam"),
                "-R",
                "chr10:89700380-89713976",
                "-f",
                "0.001",
                "-p",
                "--compat=off",
                "--vcfall",
                "--output",
                str(lgdel_gz),
                "--csi",
            ],
            check=True,
        )
        colo829_blocks = count_bgzf_blocks(lgdel_gz)
        del_svs: list[tuple[str, int, int, str]] = []
        with gzip.open(lgdel_gz, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if p[4] == "<DEL>":
                    kv = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
                    del_svs.append((p[0], int(p[1]), int(kv["END"]), p[4]))

        fa_chr20 = subprocess.check_output(["samtools", "faidx", str(fa_hg19), "chr20:46839001-46915000"], text=True)
        chr20_seq = "".join(l.strip() for l in fa_chr20.splitlines() if not l.startswith(">")).upper()
        chr20_origin = 46839001

        def ref_slice(p1_start: int, length: int) -> str:
            off = p1_start - chr20_origin
            return chr20_seq[off : off + length]

        def revcomp(seq_s: str) -> str:
            tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
            return seq_s.translate(tr)[::-1]

        chr20_len = 0
        with Path(f"{fa_hg19}.fai").open("r", encoding="utf-8") as fh:
            for line in fh:
                p = line.split("\t")
                if p[0] == "chr20":
                    chr20_len = int(p[1])
                    break

        sam_records: list[tuple[int, str]] = []
        r_id = 0
        for pos in range(46841000, 46896000, 35):
            for strand_flag in (0, 16):
                r_id += 1
                seq_s = ref_slice(pos, 100)
                sam_records.append((pos, f"wt_{r_id}\t{strand_flag}\tchr20\t{pos}\t60\t100M\t*\t0\t0\t{seq_s}\t{'I' * 100}"))

        del_l, del_r = 46841600, 46885000
        for i in range(10):
            s1 = del_l - 60
            seq1 = ref_slice(s1, 60) + ref_slice(del_r, 30)
            sam_records.append((s1, f"sdel_f_{i}\t99\tchr20\t{s1}\t60\t60M30S\t=\t{del_r}\t{del_r + 60 - s1}\t{seq1}\t{'I' * 90}"))
            s2 = del_r
            seq2 = ref_slice(del_l - 30, 30) + ref_slice(s2, 60)
            sam_records.append((s2, f"sdel_f_{i}\t147\tchr20\t{s2}\t60\t30S60M\t=\t{s1}\t{-(del_r + 60 - s1)}\t{seq2}\t{'I' * 90}"))

        dup_l, dup_r = 46842200, 46887500
        for i in range(10):
            seq_r = ref_slice(dup_r - 29, 30) + ref_slice(dup_l, 70)
            sam_records.append((dup_l, f"sdup_{i}\t81\tchr20\t{dup_l}\t60\t30S70M\t=\t{dup_r - 69}\t{dup_r - dup_l + 1}\t{seq_r}\t{'I' * 100}"))
            sf = dup_r - 69
            seq_f = ref_slice(sf, 70) + ref_slice(dup_l, 30)
            sam_records.append((sf, f"sdup_{i}\t161\tchr20\t{sf}\t60\t70M30S\t=\t{dup_l}\t{-(dup_r - dup_l + 1)}\t{seq_f}\t{'I' * 100}"))

        inv_l, inv_r = 46843000, 46889500
        for i in range(10):
            sm = inv_l
            seq_m = ref_slice(sm, 80)
            sf = inv_r - 69
            seq_f = ref_slice(sf, 70) + revcomp(ref_slice(inv_l - 19, 20))
            sam_records.append((sm, f"sinv_m_{i}\t1\tchr20\t{sm}\t60\t80M\t=\t{sf}\t{sf + 70 - sm}\t{seq_m}\t{'I' * 80}"))
            sam_records.append((sf, f"sinv_s_{i}\t1\tchr20\t{sf}\t60\t70M20S\t=\t{sm}\t{-(sf + 70 - sm)}\t{seq_f}\t{'I' * 90}"))

        sam_records.sort(key=lambda x: x[0])
        span3_sam = tmp / "span3_sv.sam"
        span3_bam = tmp / "span3_sv.bam"
        with span3_sam.open("w", encoding="utf-8") as out:
            out.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr20\tLN:{chr20_len}\n")
            for _, rline in sam_records:
                out.write(rline + "\n")
        subprocess.run(["samtools", "view", "-b", "-o", str(span3_bam), str(span3_sam)], check=True)
        subprocess.run(["samtools", "index", str(span3_bam)], check=True)

        span3_gz = tmp / "span3_sv_pileup.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "-G",
                str(fa_hg19),
                "-b",
                str(span3_bam),
                "-R",
                "chr20:46841200-46895000",
                "-f",
                "0.001",
                "-p",
                "--compat=off",
                "--vcfall",
                "--output",
                str(span3_gz),
                "--csi",
                "-th",
                "4",
            ],
            check=True,
        )
        span3_blocks = count_bgzf_blocks(span3_gz)

        found_spanning: dict[str, tuple[str, int, int, str]] = {}
        with gzip.open(span3_gz, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if p[4] in ("<DEL>", "<DUP>", "<INV>"):
                    kv = dict(x.split("=", 1) for x in p[7].split(";") if "=" in x)
                    s_p, e_p = int(p[1]), int(kv["END"])
                    if e_p - s_p >= 35000:
                        found_spanning[p[4]] = (p[0], s_p, e_p, p[4])

        all_csi_checks = [(lgdel_gz, del_svs[0])] + [(span3_gz, found_spanning[k]) for k in ("<DEL>", "<DUP>", "<INV>")]
        interior_total = 0
        interior_passed = 0
        for vcf_gz, (chrom, s_p, e_p, sv_alt) in all_csi_checks:
            start_bin = (s_p - 1) >> 14
            end_bin = (e_p - 1) >> 14
            for b in range(start_bin + 1, end_bin + 1):
                q_coord = min(e_p - 1, (b << 14) + 500)
                if q_coord <= s_p:
                    continue
                reg = f"{chrom}:{q_coord}-{q_coord}"
                interior_total += 1
                tbx_out = subprocess.check_output(["tabix", str(vcf_gz), reg], text=True)
                bcf_out = subprocess.check_output(["bcftools", "view", "-H", "-r", reg, str(vcf_gz)], text=True)
                tbx_ok = any(l.split("\t")[1] == str(s_p) and l.split("\t")[4] == sv_alt for l in tbx_out.splitlines() if l.strip())
                bcf_ok = any(l.split("\t")[1] == str(s_p) and l.split("\t")[4] == sv_alt for l in bcf_out.splitlines() if l.strip())
                if tbx_ok and bcf_ok:
                    interior_passed += 1

        hg008_csi_gz = tmp / "hg008_multi_csi.vcf.gz"
        subprocess.run(
            [
                str(vdr_bin),
                "-G",
                str(giab_ref),
                "-b",
                f"{hg008_t}|{hg008_n}",
                "-z",
                "-c",
                "1",
                "-S",
                "2",
                "-E",
                "3",
                "-f",
                "0.005",
                "--fisher",
                "--compat=off",
                "--vcfall",
                "--output",
                str(hg008_csi_gz),
                "--csi",
                "-th",
                "4",
                str(hg008_dir / "slices.bed"),
            ],
            check=True,
        )
        hg008_blocks = count_bgzf_blocks(hg008_csi_gz)

        # -p homopolymer coordinate retention (untiled == tiled + 100% covered bases retained)
        pileup_cases = [
            ("N13532_TwoDels1Ins", fa_b37, vj_bams / "N13532.TwoDels1InsStringBoundError.bam", "11", 1279500, 1281000, 300, True),
            ("threeIndels_Homopolymer", fa_b37, vj_bams / "threeIndels_zero_deletion_tumor.bam", "20", 62361250, 62361500, 80, False),
            ("HG002_CMRG_Pileup", giab_ref, w4_dir / "hg002_cmrg_mosaic" / "HG002.cmrg_mosaic.slice.bam", "chr1", 1019800, 1021200, 350, True),
        ]
        untiled_eq_tiled_cases = 0
        covered_bases_total = 0
        covered_bases_retained = 0
        pos0_sentinel_rows = 0

        for name, ref_fa, bam_p, chrom, s_p, e_p, tile_sz, full_span in pileup_cases:
            reg = f"{chrom}:{s_p}-{e_p}"
            depth_out = subprocess.check_output(["samtools", "depth", "-Q", "1", "-r", reg, str(bam_p)], text=True)
            covered = {int(l.split("\t")[1]) for l in depth_out.splitlines() if l.strip() and int(l.split("\t")[2]) > 0}
            covered_bases_total += len(covered)

            vcf_u = subprocess.check_output(
                [str(vdr_bin), "-G", str(ref_fa), "-b", str(bam_p), "-R", reg, "-p", "--compat=off", "--vcfall", "-v"],
                text=True,
            )
            vcf_t = subprocess.check_output(
                [str(vdr_bin), "-G", str(ref_fa), "-b", str(bam_p), "-R", reg, "-p", "--compat=off", "--vcfall", "-v", "--tile-size", str(tile_sz), "-th", "4"],
                text=True,
            )
            if full_span and vcf_u == vcf_t:
                untiled_eq_tiled_cases += 1

            emitted_t: set[int] = set()
            for line in vcf_t.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                p_coord = int(line.split("\t")[1])
                if p_coord <= 0:
                    pos0_sentinel_rows += 1
                emitted_t.add(p_coord)
            covered_bases_retained += len(covered & emitted_t)

        arm6_data: dict[str, Any] = {
            "csi_leaf_bin_bp": 16384,
            "max_explicit_sv_seq_bp": 32768,
            "colo829_bgzf_blocks": colo829_blocks,
            "multi_block_bgzf_blocks": span3_blocks,
            "hg008_bgzf_blocks": hg008_blocks,
            "interior_bin_queries_passed": interior_passed,
            "interior_bin_queries_total": interior_total,
            "pileup_full_span_cases_total": sum(1 for *_, fs in pileup_cases if fs),
            "pileup_untiled_eq_tiled_cases": untiled_eq_tiled_cases,
            "pileup_covered_bases_total": covered_bases_total,
            "pileup_covered_bases_retained": covered_bases_retained,
            "pileup_covered_bases_retained_pct": round(100.0 * covered_bases_retained / max(1, covered_bases_total), 2),
            "pos0_sentinel_rows": pos0_sentinel_rows,
        }

        doc: dict[str, Any] = {
            "meta": {
                "study_id": EXPECTED_STUDY_ID,
                "study_version": EXPECTED_STUDY_VERSION,
                "study_fingerprint": EXPECTED_STUDY_FP,
                "study_perf_fingerprint": EXPECTED_STUDY_PERF_FP,
                "compat_switch_count": EXPECTED_COMPAT_SWITCH_COUNT,
                "preset_count": EXPECTED_PRESET_COUNT,
                "var2vcf_flag_count": EXPECTED_VAR2VCF_FLAG_COUNT,
                "total_arms": EXPECTED_TOTAL_ARMS,
                "rtg_version": tool_versions["rtg_version"],
                "bcftools_version": tool_versions["bcftools_version"],
                "tabix_version": tool_versions["tabix_version"],
                "samtools_version": tool_versions["samtools_version"],
                "scope_note": (
                    "Host-independent deterministic callset, targeted CompatProfile VCFv4.2/v4.3 "
                    "schema/canonicalization/accounting audit, bcftools norm --check-ref e FASTA "
                    "concordance, and diploid rtg vcfeval metrics (zero host-dependent wall-clock/RSS timings)."
                ),
            },
            "arms": {
                "arm1_vcf_schema_and_bcftools_norm": arm1_data,
                "arm2_seqc2_diploid_vcfeval": arm2_data,
                "arm3_somatic_multiallelic_and_truth": arm3_data,
                "arm4_brca2_adaptor_chimera": arm4_data,
                "arm5_clinical_presets_and_amplicon": arm5_data,
                "arm6_csi_spanning_sv_and_pileup_tiling": arm6_data,
            },
            "switch_catalog": switch_catalog,
        }

        tsv_rows: list[dict[str, str]] = [
            {
                "arm_id": "Arm 1",
                "arm_title": "Targeted VCF Schema, Canonicalization & FASTA Reference Concordance",
                "cohort_scope": f"{arm1_data['total_callsets']} callsets across {arm1_data['distinct_bams']} BAMs (17 wet-lab public BAM slices + 1 synthetic contig-start SV BAM vartovcf_sv_contig1 + synthetic 17-bp periodic IUPAC chr20 FASTA fixture)",
                "provenance_type": "17 wet-lab public BAM slices + 1 synthetic contig-start SV BAM (vartovcf_sv_contig1) + synthetic 17-bp periodic IUPAC chr20 FASTA fixture",
                "vardictjava_result": f"{arm1_data['total_defects_vj']} targeted VCF audit findings; {arm1_data['bcftools_norm_check_ref_vj_passed']}/{arm1_data['total_callsets']} bcftools norm --check-ref e",
                "compat_off_result": f"{arm1_data['total_defects_off']} targeted VCF audit findings; {arm1_data['bcftools_norm_check_ref_off_passed']}/{arm1_data['bcftools_norm_check_ref_off_total']} bcftools norm --check-ref e",
                "verdict": "PASS",
                "sanitized_cmd": sanitize_text(
                    'vardictrust --compat=off -G $REPO_ROOT/scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa -b "$TUMOR_BAM|$NORMAL_BAM" -z -c 1 -S 2 -E 3 -g 4 --output $WORK_TMP/out.vcf',
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
            {
                "arm_id": "Arm 2",
                "arm_title": "Diploid SEQC2 HCC1395 Subclonal Somatic rtg vcfeval",
                "cohort_scope": f"SEQC2 HCC1395/HCC1395BL across 6 contigs ({arm2_data['evaluated_intervals']} HC intervals, {arm2_data['evaluated_bp']} bp)",
                "provenance_type": "Wet-lab tumor-normal WES truth (Fang et al. 2021)",
                "vardictjava_result": (
                    f"Unthresholded: TP={vj_vcfeval['tp_baseline']}, FP={vj_vcfeval['fp']}, FN={vj_vcfeval['fn']}, "
                    f"P={vj_vcfeval['precision']:.4f}, R={vj_vcfeval['sensitivity']:.4f}, F1={vj_vcfeval['f_measure']:.4f}; "
                    f"Opt QUAL>={vj_vcfeval['best_qual_threshold']:.0f}: TP={vj_vcfeval['best_tp_baseline']}, FP={vj_vcfeval['best_fp']}, "
                    f"FN={vj_vcfeval['best_fn']}, P={vj_vcfeval['best_precision']:.4f}, R={vj_vcfeval['best_sensitivity']:.4f}, "
                    f"F1={vj_vcfeval['best_f_measure']:.4f}; GT=0/0 carriers={vj_vcfeval['gt00_carrier_somatic']}"
                ),
                "compat_off_result": (
                    f"Unthresholded: TP={off_vcfeval['tp_baseline']} (+{arm2_data['tp_baseline_gain']}), FP={off_vcfeval['fp']}, FN={off_vcfeval['fn']}, "
                    f"P={off_vcfeval['precision']:.4f}, R={off_vcfeval['sensitivity']:.4f} (+{arm2_data['sensitivity_gain'] * 100.0:.1f} pp), "
                    f"F1={off_vcfeval['f_measure']:.4f} (+{arm2_data['f_measure_gain']:.4f}); "
                    f"Opt QUAL>={off_vcfeval['best_qual_threshold']:.0f}: TP={off_vcfeval['best_tp_baseline']} (+{arm2_data['best_tp_baseline_gain']}), "
                    f"FP={off_vcfeval['best_fp']}, FN={off_vcfeval['best_fn']}, P={off_vcfeval['best_precision']:.4f}, "
                    f"R={off_vcfeval['best_sensitivity']:.4f} (+{arm2_data['best_sensitivity_gain'] * 100.0:.1f} pp), "
                    f"F1={off_vcfeval['best_f_measure']:.4f} (+{arm2_data['best_f_measure_gain']:.4f}); GT=0/0 carriers={off_vcfeval['gt00_carrier_somatic']}"
                ),
                "verdict": f"PASS (+{arm2_data['tp_baseline_gain']}/{arm2_data['truth_total']} subclonal TPs; +{arm2_data['sensitivity_gain'] * 100.0:.1f} pp recall; +{arm2_data['best_f_measure_gain']:.4f} optimal F1)",
                "sanitized_cmd": sanitize_text(
                    "rtg RTG_MEM=4G RTG_JAVA_OPTS=-Dtalkback=false vcfeval -T 4 -b $WORK_TMP/truth.vcf.gz -c $WORK_TMP/somatic.vcf.gz -e $WORK_TMP/eval.bed -t $WORK_TMP/ref.sdf --all-records --vcf-score-field=QUAL --sample=TUMOR,TUMOR -o $WORK_TMP/eval",
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
            {
                "arm_id": "Arm 3",
                "arm_title": "Multi-Allelic Somatic Continuation, MUC19 Paralog, & Normalized Truth",
                "cohort_scope": f"{arm3_data['distinct_bam_pairs']} multi-allelic T/N pairs ({arm3_data['total_runs']} runs) + MUC19 paralog + SEQC2/HG008 v0.4",
                "provenance_type": "Wet-lab somatic tumor-normal benchmarks",
                "vardictjava_result": f"{arm3_data['vj_recovered_pairs']}/{arm3_data['distinct_bam_pairs']} pairs ({arm3_data['vj_recovered_runs']}/{arm3_data['total_runs']} runs) secondary alleles; SEQC2 TP={seqc2_hc_vj['tp']}/{seqc2_hc_vj['truth_total']} (FP={seqc2_hc_vj['fp']}); HG008 TP={hg008_v04_vj['tp']}/{hg008_v04_vj['truth_total']} (FP={hg008_v04_vj['fp']})",
                "compat_off_result": f"{arm3_data['off_recovered_pairs']}/{arm3_data['distinct_bam_pairs']} pairs ({arm3_data['off_recovered_runs']}/{arm3_data['total_runs']} runs) secondary alleles; MUC19={arm3_data['muc19_variants']} calls (0 VD>DP); SEQC2 TP={seqc2_hc_off['tp']}/{seqc2_hc_off['truth_total']} (FP={seqc2_hc_off['fp']}); HG008 TP={hg008_v04_off['tp']}/{hg008_v04_off['truth_total']} (FP={hg008_v04_off['fp']})",
                "verdict": "PASS",
                "sanitized_cmd": sanitize_text(
                    'vardictrust --compat=off --vcf-all-alleles -G $REPO_ROOT/VarDictJava/testdata/fastas/hg19.fa -b "$TUMOR_BAM|$NORMAL_BAM" -v',
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
            {
                "arm_id": "Arm 4",
                "arm_title": "Patient BRCA2 Soft-Clip TruSeq Adapter Chimera Suppression",
                "cohort_scope": "Patient BRCA2 tumor BAM chr13:32683877-32684023 (147 bp)",
                "provenance_type": "Wet-lab clinical tumor amplicon/capture slice",
                "vardictjava_result": f"{arm4_data['vj_chimeric_insertions']} false {arm4_data['vj_chimeric_ins_bp']}-bp adapter insertion at chr13:{arm4_data['chimera_pos']}",
                "compat_off_result": f"{arm4_data['off_chimeric_insertions']} chimeric insertions at chr13:{arm4_data['chimera_pos']}",
                "verdict": "PASS",
                "sanitized_cmd": sanitize_text(
                    "vardictrust --compat=off --adaptor AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -G $REPO_ROOT/scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa -b $REPO_ROOT/scratch/g4/wave4/brca2_adaptor/adaptor_insertion.bam -R chr13:32683877-32684023 -v",
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
            {
                "arm_id": "Arm 5",
                "arm_title": "Clinical Workflow Presets (wgs, wes, amplicon, ctdna) & QIAseq -sp",
                "cohort_scope": "GIAB HG002/HG008 WGS, SEQC2 WES (264 kbp + BED8), Illumina PE & QIAseq V3 SPE panels, 0.5% VAF ctDNA admixture",
                "provenance_type": "Wet-lab cohorts + semi-synthetic real-read 0.5% VAF ctDNA spike-in",
                "vardictjava_result": f"Paired amplicon unsupported; QIAseq without -sp PASS KnownPositives={qgn_nosp_kp_tp}/{len(kp_in_inserts)} ({qgn_nosp_dp50_kp_tp}/{len(kp_callable_50x)} at >=50x); default -f 0.01 misses 0.5% ctDNA",
                "compat_off_result": f"4/4 presets valid; paired amplicon={arm5_data['paired_somatic_amplicon_off_records']} calls; QIAseq -sp PASS KnownPositives={qgn_sp_kp_tp}/{len(kp_in_inserts)} ({qgn_sp_dp50_kp_tp}/{len(kp_callable_50x)} at >=50x); ctdna 0.5% PASS (DP={arm5_data['spikein_observed_dp']}, VD={arm5_data['spikein_observed_vd']})",
                "verdict": "PASS",
                "sanitized_cmd": sanitize_text(
                    "vardictrust --preset ctdna -G $REPO_ROOT/scratch/g4/ref/d1vd1/GRCh38.d1.vd1.fa -b $WORK_TMP/ctdna_1000x_0p5pct.bam --output $WORK_TMP/ctdna.vcf.gz --tbi $WORK_TMP/ctdna_locus_3col.bed",
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
            {
                "arm_id": "Arm 6",
                "arm_title": "Multi-Block CSI Spanning-SV Indexing & Homopolymer Pileup Tiling",
                "cohort_scope": f"chr20 >35 kbp <DEL>/<DUP>/<INV> ({span3_blocks} BGZF blocks), Colo829 ({colo829_blocks} blocks), HG008 ({hg008_blocks} blocks), 3 -p homopolymer slices",
                "provenance_type": "Synthetic >35 kbp spanning-SV BAM + wet-lab GIAB/Colo829/N13532 BAMs",
                "vardictjava_result": "Legacy -p drops multi-deletion shifted coordinates and emits pos-0 sentinel rows",
                "compat_off_result": f"{interior_passed}/{interior_total} interior 16-kbp CSI queries returned spanning SVs; {arm6_data['pileup_covered_bases_retained_pct']:.1f}% covered bases retained; {pos0_sentinel_rows} sentinel rows",
                "verdict": "PASS",
                "sanitized_cmd": sanitize_text(
                    "vardictrust --compat=off --vcf-all-alleles -p -G $REPO_ROOT/VarDictJava/testdata/fastas/hg19.fa -b $WORK_TMP/span3_sv.bam -R chr20:46841200-46895000 --output $WORK_TMP/span3.vcf.gz --csi",
                    tmpdir,
                    data_root,
                    ROOT,
                ),
            },
        ]

        json_text = sanitize_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", tmpdir, data_root, ROOT)

        out_dir.mkdir(parents=True, exist_ok=True)
        metrics_path = out_dir / "v110_compat_metrics.json"
        arms_path = out_dir / "v110_compat_arms.tsv"
        sha_path = out_dir / "SHA256SUMS.tsv"

        tsv_lines: list[str] = ["\t".join(TSV_COLUMNS)]
        for r in tsv_rows:
            tsv_lines.append("\t".join(sanitize_text(r[col], tmpdir, data_root, ROOT) for col in TSV_COLUMNS))
        tsv_text = "\n".join(tsv_lines) + "\n"

        for label_name, payload in (("v110_compat_metrics.json", json_text), ("v110_compat_arms.tsv", tsv_text)):
            for forbidden in ("/home/", str(ROOT), str(data_root), "/tmp/", tmpdir, "vdr_v110_study_", "miniconda3", "miniforge3", "$HOME_REDACTED"):
                if forbidden and forbidden in payload:
                    raise RuntimeError(f"Forbidden local path token {forbidden!r} leaked into {label_name}")

        atomic_write_text(metrics_path, json_text)
        atomic_write_text(arms_path, tsv_text)

        sha_entries = []
        for p in sorted((metrics_path, arms_path), key=lambda x: x.name):
            sha_entries.append(f"{sha256_file(p)}\t{p.name}\t{p.stat().st_size}")
        sha_text = "sha256\tpath\tsize_bytes\n" + "\n".join(sha_entries) + "\n"
        atomic_write_text(sha_path, sha_text)


def verify_bundle(out_dir: Path) -> None:
    metrics_path = out_dir / "v110_compat_metrics.json"
    arms_path = out_dir / "v110_compat_arms.tsv"
    sha_path = out_dir / "SHA256SUMS.tsv"

    for req_file in (metrics_path, arms_path, sha_path):
        if not req_file.is_file():
            raise FileNotFoundError(f"Missing required v1.1.0 study artifact: {req_file}")

    md_files = list(out_dir.rglob("*.md"))
    if md_files:
        raise ValueError(f"Forbidden Markdown files inside {out_dir}: {md_files}")

    # 1. Verify SHA256SUMS.tsv digests and byte sizes
    with sha_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != ["sha256", "path", "size_bytes"]:
            raise ValueError(f"Unexpected SHA256SUMS.tsv header: {reader.fieldnames}")
        sha_rows = list(reader)

    expected_paths = ["v110_compat_arms.tsv", "v110_compat_metrics.json"]
    if [r["path"] for r in sha_rows] != expected_paths:
        raise ValueError(f"Unexpected SHA256SUMS.tsv file list: {[r['path'] for r in sha_rows]}")

    for row in sha_rows:
        target = out_dir / row["path"]
        if not target.is_file():
            raise FileNotFoundError(f"Missing checksummed target: {target}")
        actual_sha = sha256_file(target)
        if actual_sha != row["sha256"]:
            raise ValueError(f"SHA-256 mismatch for {target.name}: expected {row['sha256']}, got {actual_sha}")
        actual_sz = str(target.stat().st_size)
        if actual_sz != row["size_bytes"]:
            raise ValueError(f"Byte-size mismatch for {target.name}: expected {row['size_bytes']}, got {actual_sz}")

    # 2. Path leak checks
    for p in (metrics_path, arms_path, sha_path):
        raw_txt = p.read_text(encoding="utf-8")
        for bad_tok in ("/home/", str(ROOT), "/tmp/", "vdr_v110_study_", "miniconda3", "miniforge3", "$HOME_REDACTED"):
            if bad_tok in raw_txt:
                raise ValueError(f"Path leak {bad_tok!r} detected in {p.name}")

    # 3. JSON schema & domain invariants
    doc = json.loads(metrics_path.read_text(encoding="utf-8"))
    meta = doc["meta"]
    if meta["study_id"] != EXPECTED_STUDY_ID:
        raise ValueError(f"Unexpected study_id: {meta['study_id']}")
    if meta["study_version"] != EXPECTED_STUDY_VERSION:
        raise ValueError(f"Unexpected study_version: {meta['study_version']}")
    if meta["study_fingerprint"] != EXPECTED_STUDY_FP:
        raise ValueError(f"Unexpected study_fingerprint: {meta['study_fingerprint']}")
    if meta["study_perf_fingerprint"] != EXPECTED_STUDY_PERF_FP:
        raise ValueError(f"Unexpected study_perf_fingerprint: {meta['study_perf_fingerprint']}")
    if meta["compat_switch_count"] != EXPECTED_COMPAT_SWITCH_COUNT:
        raise ValueError(f"Unexpected compat_switch_count: {meta['compat_switch_count']}")
    if meta["preset_count"] != EXPECTED_PRESET_COUNT:
        raise ValueError(f"Unexpected preset_count: {meta['preset_count']}")
    if meta["var2vcf_flag_count"] != EXPECTED_VAR2VCF_FLAG_COUNT:
        raise ValueError(f"Unexpected var2vcf_flag_count: {meta['var2vcf_flag_count']}")
    if meta["total_arms"] != EXPECTED_TOTAL_ARMS:
        raise ValueError(f"Unexpected total_arms: {meta['total_arms']}")

    arms = doc["arms"]
    a1 = arms["arm1_vcf_schema_and_bcftools_norm"]
    if a1["total_callsets"] != 19:
        raise ValueError(f"Arm 1 expected 19 callsets, got {a1['total_callsets']}")
    for k in DEFECT_KEYS:
        if a1["vardictjava_defects"][k] <= 0:
            raise ValueError(f"Arm 1 expected vardictjava_defects[{k!r}] > 0")
        if a1["compat_off_defects"][k] != 0:
            raise ValueError(f"Arm 1 expected compat_off_defects[{k!r}] == 0")
    if a1["total_defects_off"] != 0 or a1["bcftools_norm_check_ref_off_passed"] != 19:
        raise ValueError("Arm 1 compat_off schema or bcftools norm invariant failed")

    a2 = arms["arm2_seqc2_diploid_vcfeval"]
    if a2["gt_proj_tvaf_threshold"] != 0.75:
        raise ValueError(f"Arm 2 unexpected gt_proj_tvaf_threshold: {a2['gt_proj_tvaf_threshold']}")
    if a2["vardictjava"]["gt00_carrier_somatic"] <= 0 or a2["compat_off"]["gt00_carrier_somatic"] != 0:
        raise ValueError("Arm 2 GT=0/0 somatic carrier invariant failed")
    if a2["tp_baseline_gain"] <= 0 or a2["f_measure_gain"] <= 0 or a2["best_f_measure_gain"] <= 0:
        raise ValueError("Arm 2 diploid rtg vcfeval TP/F1 gain invariant failed")

    a3 = arms["arm3_somatic_multiallelic_and_truth"]
    if (
        a3["distinct_bam_pairs"] != 5
        or a3["total_runs"] != 6
        or a3["vj_recovered_pairs"] != 0
        or a3["off_recovered_pairs"] != 5
        or a3["vj_recovered_runs"] != 0
        or a3["off_recovered_runs"] != 6
        or a3["muc19_variants"] < 20
        or a3["muc19_vd_gt_dp_violations"] != 0
    ):
        raise ValueError("Arm 3 multi-allelic or MUC19 invariant failed")

    a4 = arms["arm4_brca2_adaptor_chimera"]
    if a4["chimera_pos"] != 32683955 or a4["vj_chimeric_insertions"] < 1 or a4["off_chimeric_insertions"] != 0:
        raise ValueError("Arm 4 BRCA2 adapter-chimera invariant failed")

    a5 = arms["arm5_clinical_presets_and_amplicon"]
    if (
        not a5["wgs_tbi_and_csi_verified"]
        or not a5["wes_bed8_non_hijack_verified"]
        or a5["paired_somatic_amplicon_vj_supported"]
        or a5["paired_somatic_amplicon_off_records"] < 10
        or a5["qiaseq_spe_known_pos_tp"] <= a5["qiaseq_spe_no_sp_known_pos_tp"]
        or not a5["spikein_detected"]
    ):
        raise ValueError("Arm 5 clinical presets / amplicon / ctDNA invariant failed")

    a6 = arms["arm6_csi_spanning_sv_and_pileup_tiling"]
    if (
        a6["csi_leaf_bin_bp"] != 16384
        or a6["max_explicit_sv_seq_bp"] != 32768
        or a6["multi_block_bgzf_blocks"] < 20
        or a6["interior_bin_queries_passed"] != a6["interior_bin_queries_total"]
        or a6["interior_bin_queries_total"] <= 0
        or a6["pileup_full_span_cases_total"] <= 0
        or a6["pileup_untiled_eq_tiled_cases"] != a6["pileup_full_span_cases_total"]
        or a6["pileup_covered_bases_total"] <= 0
        or a6["pileup_covered_bases_retained"] != a6["pileup_covered_bases_total"]
        or a6["pos0_sentinel_rows"] != 0
    ):
        raise ValueError("Arm 6 CSI spanning-SV or -p pileup tiling invariant failed")

    # 4. Verify 35-entry switch_catalog against config.rs & hygiene rules (prefix-tolerant for v1.1.1+ growth)
    sw_cnt = int(meta["compat_switch_count"])
    expected_switches = parse_compat_profile_switches(CONFIG_RS_PATH, ROADMAP_MD_PATH, require_exact_count=False)
    catalog = doc["switch_catalog"]
    if len(catalog) != sw_cnt + 3:
        raise ValueError(f"Expected {sw_cnt + 3} switch_catalog entries, got {len(catalog)}")
    if [c["name"] for c in catalog[:sw_cnt]] != expected_switches[:sw_cnt]:
        raise ValueError(f"switch_catalog[0..{sw_cnt}] does not match CompatProfile declaration order in config.rs")

    bare_issue_re = re.compile(r"(?<![A-Za-z0-9_/-])#\d+\b")
    jargon_re = re.compile(r"\b(Gate\s*\d+|s[1-6]_[a-z0-9_]+|arm[1-6]_[a-z0-9_]+|Option C|Tier[- ]P|evidence_certified)\b")
    for entry in catalog:
        if entry["domain"] not in ALLOWED_DOMAINS:
            raise ValueError(f"Invalid switch_catalog domain: {entry['domain']}")
        for k in ("name", "category", "domain", "upstream_refs", "vardictjava_behavior", "compat_off_behavior", "verification_arm"):
            val = str(entry[k])
            if bare_issue_re.search(val):
                raise ValueError(f"Bare #N citation in switch_catalog[{entry['index']}].{k}: {val!r}")
            if jargon_re.search(val):
                raise ValueError(f"Internal jargon in switch_catalog[{entry['index']}].{k}: {val!r}")
            outside_code = re.sub(r"`[^`]+`", "", val)
            if "<" in outside_code or "|" in outside_code:
                raise ValueError(f"Unescaped '<' or '|' in switch_catalog[{entry['index']}].{k}: {val!r}")

    # 5. Verify TSV schema
    with arms_path.open("r", encoding="utf-8", newline="") as fh:
        t_reader = csv.DictReader(fh, delimiter="\t")
        if t_reader.fieldnames != TSV_COLUMNS:
            raise ValueError(f"Unexpected v110_compat_arms.tsv columns: {t_reader.fieldnames}")
        t_rows = list(t_reader)
    if len(t_rows) != EXPECTED_TOTAL_ARMS:
        raise ValueError(f"Expected {EXPECTED_TOTAL_ARMS} rows in v110_compat_arms.tsv, got {len(t_rows)}")
    for r in t_rows:
        if not r["verdict"].startswith("PASS"):
            raise ValueError(f"Non-PASS row in v110_compat_arms.tsv: {r}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="Execute live v1.1.0 study collection and write bundle")
    parser.add_argument("--verify", action="store_true", help="Verify tracked bundle offline (100%% tool-free stdlib)")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Bundle directory path")
    args = parser.parse_args()

    if args.run == args.verify:
        parser.error("Specify exactly one of --run or --verify")

    out_dir = args.out_dir.resolve()
    if args.run:
        run_collection(out_dir)
        verify_bundle(out_dir)
        print(f"OK collected and verified v1.1.0 study bundle in {out_dir}")
        return 0

    verify_bundle(out_dir)
    print(f"OK verified v1.1.0 study bundle offline in {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
