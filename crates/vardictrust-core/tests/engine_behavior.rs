//! First-class functional engine behavior integration suite (`34` tests).
//!
//! Exercises CRAM decoding, Fisher exact test tables/VCFs, VCF filter edge cases,
//! structural variant calling, read/pileup filtering flags, hostile inputs, and
//! multi-thread tile determinism entirely in-process via
//! `vardictrust_core::run_with_writer` with zero `samtools` or `vardictrust`
//! subprocesses.

use std::fs;
use std::sync::atomic::Ordering;

use noodles::core::Position;
use noodles::sam::alignment::{
    record::cigar::{op::Kind, Op},
    record_buf::{Cigar, RecordBuf},
};
use tempfile::TempDir;
use vardictrust_core::{
    config::{extract_sample_name, Config},
    types::VardictError,
};

mod common;
use common::{
    data_root, ensure_fastas_converted, extract_variants, indexed_bam, run_core_var, run_core_vcf,
    snv_record, workspace_root, write_indexed_cram, write_ref,
};

fn base_micro_config() -> Config {
    let root = workspace_root();
    let bam = root
        .join("fixtures/micro/micro.bam")
        .to_string_lossy()
        .into_owned();
    Config {
        reference: root
            .join("fixtures/micro/micro.ref.fa")
            .to_string_lossy()
            .into_owned(),
        sample: extract_sample_name(&bam, None),
        bam,
        threads: 1,
        region: Some("chr1:40-260".to_string()),
        ..Config::default()
    }
}

fn normalize_var_rows(text: &str, expected_cols: usize) -> Vec<String> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut cols: Vec<&str> = l.split('\t').collect();
            if !cols.is_empty() {
                cols[0] = "SAMPLE";
            }
            assert_eq!(
                cols.len(),
                expected_cols,
                ".var row must emit {expected_cols} columns"
            );
            cols.join("\t")
        })
        .collect()
}

// 1. mod cram (1 test)
mod cram {
    use super::*;

    #[test]
    fn soft_masked_lowercase_ref_and_omitted_nm_matches_bam() {
        ensure_fastas_converted();
        let droot = data_root();
        let ref_fa = droot.join("VarDictJava/testdata/fastas/hg19.fa");
        let bam_path = droot.join(
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CRISPR.bam",
        );

        let tmp_dir = tempfile::tempdir().expect("tempdir");
        let cram_path = tmp_dir.path().join("crispr.cram");
        write_indexed_cram(&bam_path, &ref_fa, &cram_path, true, Some("chr20"))
            .expect("pure-Rust CRAM 3.0 encode with stripped NM");

        let cfg_bam = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: bam_path.to_string_lossy().into_owned(),
            sample: "COLO829".to_string(),
            region: Some("chr20:50454281-50455444".to_string()),
            freq: 0.001,
            af_threshold: 0.001,
            threads: 1,
            ..Config::default()
        };

        let mut cfg_cram = cfg_bam.clone();
        cfg_cram.bam = cram_path.to_string_lossy().into_owned();

        let vcf_bam = run_core_vcf(&cfg_bam).expect("BAM VCF");
        let vcf_cram = run_core_vcf(&cfg_cram).expect("CRAM VCF");

        let bam_vars = extract_variants(&vcf_bam);
        let cram_vars = extract_variants(&vcf_cram);
        assert_eq!(bam_vars.len(), 8, "expected 8 real COLO829 variants");
        assert_eq!(
            bam_vars, cram_vars,
            "soft-masked hg19 CRAM with omitted NM must produce byte-identical VCF variants to BAM"
        );

        let mut cfg_tiled = cfg_cram.clone();
        cfg_tiled.tile_size = Some(400);
        let vcf_tiled = run_core_vcf(&cfg_tiled).expect("tiled CRAM VCF");
        assert_eq!(
            bam_vars,
            extract_variants(&vcf_tiled),
            "tiled CRAM run (--tile-size 400) must match untiled BAM output"
        );
    }
}

// 2. mod fisher (1 test)
mod fisher {
    use super::*;

    #[test]
    fn java_exact_test_vcf_and_intermediate_table() {
        ensure_fastas_converted();
        let droot = data_root();

        // (A) Simple --fisher: 38-column .var + VCF PMEAN/PSTD/QSTD/SBF/ODDRATIO
        let mut cfg_simple = base_micro_config();
        cfg_simple.fisher = true;
        let simple_var = run_core_var(&cfg_simple).expect("simple --fisher .var");
        let simple_rows = normalize_var_rows(&simple_var, 38);
        assert!(!simple_rows.is_empty(), "simple --fisher must emit rows");
        let simple_vcf = run_core_vcf(&cfg_simple).expect("simple --fisher VCF");
        let simple_vars = extract_variants(&simple_vcf);
        assert!(
            !simple_vars.is_empty()
                && simple_vars.iter().all(|v| {
                    v.contains("PMEAN=")
                        && v.contains("PSTD=")
                        && v.contains("QSTD=")
                        && v.contains("SBF=")
                        && v.contains("ODDRATIO=")
                }),
            "simple --fisher VCF must emit PMEAN, PSTD, QSTD, SBF, and ODDRATIO"
        );

        // (B) Somatic --fisher: 61-column .var oracle comparison + VCF SSF/SOR
        let ref_hg19 = droot.join("VarDictJava/testdata/fastas/hg19.fa");
        let bam_dir =
            droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
        let bam_pair = format!(
            "{}|{}",
            bam_dir.join("ArgumentNA_Tumor.bam").display(),
            bam_dir.join("ArgumentNA_Normal.bam").display()
        );
        let mut cfg_somatic = Config {
            reference: ref_hg19.to_string_lossy().into_owned(),
            sample: extract_sample_name(&bam_pair, None),
            bam: bam_pair,
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            column_for_chromosome: Some(0),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            region: Some("chr7:140453010-140453158:testbed".to_string()),
            base_qual: 25.0,
            flag_filter: 0x100,
            freq: 0.001,
            af_threshold: 0.001,
            realign: 1,
            minmatch: 20,
            pileup: true,
            fisher: true,
            ..Config::default()
        };
        let somatic_var = run_core_var(&cfg_somatic).expect("somatic --fisher .var");
        let expected_path = droot.join(
            "VarDictJava/testdata/integrationtestcases/Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt",
        );
        let expected_raw = fs::read_to_string(&expected_path).expect("read somatic fisher oracle");
        let expected_body: String = expected_raw
            .lines()
            .skip(1)
            .collect::<Vec<&str>>()
            .join("\n");
        let actual_rows = normalize_var_rows(&somatic_var, 61);
        let expected_rows = normalize_var_rows(&expected_body, 61);
        assert_eq!(actual_rows, expected_rows);

        cfg_somatic.pileup = false;
        let somatic_vcf = run_core_vcf(&cfg_somatic).expect("somatic --fisher VCF");
        let somatic_vars = extract_variants(&somatic_vcf);
        assert!(
            !somatic_vars.is_empty()
                && somatic_vars
                    .iter()
                    .all(|v| v.contains(";SSF=") && v.contains(";SOR=")),
            "somatic --fisher VCF must emit SSF and SOR"
        );

        // (C) Amplicon --fisher: 40-column .var table
        let tmp = tempfile::tempdir().expect("tempdir");
        let amp_bed = tmp.path().join("amp.bed");
        fs::write(&amp_bed, "chr1\t81\t140\tmicro_amp\t0\t+\t81\t140\n").expect("write amp bed");
        let mut cfg_amp = base_micro_config();
        cfg_amp.region = None;
        cfg_amp.bed = Some(amp_bed.to_string_lossy().into_owned());
        cfg_amp.zero_based = true;
        cfg_amp.col_chr = 0;
        cfg_amp.col_start = 1;
        cfg_amp.col_end = 2;
        cfg_amp.col_seg_start = 1;
        cfg_amp.col_seg_end = 2;
        cfg_amp.col_gene = Some(3);
        cfg_amp.column_for_chromosome = Some(0);
        cfg_amp.omit_end_tag = true;
        cfg_amp.amplicon = Some("30:0.5".to_string());
        cfg_amp.freq = 0.01;
        cfg_amp.af_threshold = 0.01;
        cfg_amp.min_reads = 1;
        cfg_amp.fisher = true;
        let amp_var = run_core_var(&cfg_amp).expect("amplicon --fisher .var");
        let amp_rows = normalize_var_rows(&amp_var, 40);
        assert!(
            !amp_rows.is_empty(),
            "amplicon --fisher must emit 40-col rows"
        );

        // (D) R-mode (145, 494, 1, 0) exact 5-decimal midpoint regression check (73/320 = 0.228125):
        // glibc libm (0x3fcd333333333333) rounds p_less to 0.22812; musl libm (0x3fcd333333333334)
        // is 1 ULP higher across transcendental dhyper evaluation and rounds p_less to 0.22813.
        let fe_r = vardictrust_core::stats::fisher::FisherExact::new(145, 494, 1, 0);
        assert_eq!(fe_r.get_p_value().to_bits(), 0.22812_f64.to_bits());
        assert_eq!(fe_r.get_p_value_greater().to_bits(), 1.0_f64.to_bits());
        assert_eq!(fe_r.get_odds_ratio(), "0");
        #[cfg(not(target_env = "musl"))]
        assert_eq!(fe_r.get_p_value_less().to_bits(), 0.22812_f64.to_bits());
        #[cfg(target_env = "musl")]
        assert_eq!(fe_r.get_p_value_less().to_bits(), 0.22813_f64.to_bits());
    }
}

// 3. mod vcf_filters (1 test)
mod vcf_filters {
    use super::*;

    #[test]
    fn zero_allele_frequency_perl_falsy_fallback_vs_zero_decimal() {
        ensure_fastas_converted();
        let droot = data_root();
        let ref_fa = droot.join("VarDictJava/testdata/fastas/GRCh37.fa");
        let bam = droot.join(
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_20156208_20234556.bam",
        );

        let mut cfg_f0 = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            sample: "NA24385".to_string(),
            bam: bam.to_string_lossy().into_owned(),
            region: Some("20:20156208-20157500".to_string()),
            freq: 0.0,
            af_threshold: 0.0,
            af_threshold_truthy_zero: false, // -f 0
            threads: 1,
            ..Config::default()
        };
        let mut cfg_f0_dec = cfg_f0.clone();
        cfg_f0_dec.af_threshold_truthy_zero = true; // -f 0.0

        let stdout_f0 = run_core_vcf(&cfg_f0).expect("VCF with -f 0");
        let stdout_f0_dec = run_core_vcf(&cfg_f0_dec).expect("VCF with -f 0.0");

        assert!(
            stdout_f0.contains("##FILTER=<ID=f0.02,Description=\"Allele frequency < 0.02\">"),
            "-f 0 header must emit f0.02 filter description"
        );
        assert!(
            stdout_f0_dec.contains("##FILTER=<ID=f0,Description=\"Allele frequency < 0\">"),
            "-f 0.0 header must emit f0 filter description"
        );

        let vars_f0 = extract_variants(&stdout_f0);
        let vars_f0_dec = extract_variants(&stdout_f0_dec);
        assert_eq!(vars_f0.len(), vars_f0_dec.len());
        assert!(
            vars_f0
                .iter()
                .any(|v| v.contains("f0.02") && v.contains("\t0/0:")),
            "sub-2% AF candidate in NA24385 must receive f0.02 filter and GT=0/0 under -f 0"
        );
        assert!(
            !vars_f0_dec.iter().any(|v| v.contains("f0.02")),
            "-f 0.0 must not apply f0.02 filter to sub-2% AF candidates"
        );
        let _ = &mut cfg_f0;
    }
}

// 4. mod structural_variants (3 tests)
mod structural_variants {
    use super::*;

    #[test]
    fn duplication_5prime_deletion_and_out_of_bounds_deletion() {
        ensure_fastas_converted();
        let droot = data_root();
        let cases = [
            (
                "VarDictJava/testdata/fastas/hg19.fa",
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.MarkingDUPSV.bam",
                "chr20:47131954-47132045:testbed",
                "VarDictJava/testdata/integrationtestcases/Simple;hg19.fa;Colo829.chr20.MarkingDUPSV.bam;chr20;47131954-47132045;-f 0.001;.txt",
                Some("C/+196"),
            ),
            (
                "VarDictJava/testdata/fastas/hg19.fa",
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CandDel5SV.bam",
                "chr20:506503-507295:testbed",
                "VarDictJava/testdata/integrationtestcases/Simple;hg19.fa;Colo829.chr20.CandDel5SV.bam;chr20;506503-507295;-f 0.001;.txt",
                Some("T/-632"),
            ),
            (
                "VarDictJava/testdata/fastas/hg38.fa",
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/DelOutOfBound.chr1_150578117_150578617.bam",
                "chr1:150578117-150578617:testbed",
                "VarDictJava/testdata/integrationtestcases/Simple;hg38.fa;DelOutOfBound.bam;chr1;150578117-150578617;-f0.001.txt",
                Some("C/-752"),
            ),
        ];

        for (fa_rel, bam_rel, region, oracle_rel, required_tok) in cases {
            let bam = droot.join(bam_rel).to_string_lossy().into_owned();
            let cfg = Config {
                reference: droot.join(fa_rel).to_string_lossy().into_owned(),
                sample: extract_sample_name(&bam, None),
                bam,
                threads: 1,
                zero_based: true,
                omit_end_tag: true,
                column_for_chromosome: Some(0),
                col_chr: 0,
                col_start: 1,
                col_end: 2,
                col_seg_start: 1,
                col_seg_end: 2,
                col_gene: Some(3),
                region: Some(region.to_string()),
                freq: 0.001,
                af_threshold: 0.001,
                ..Config::default()
            };
            let var_out = run_core_var(&cfg).expect("SV fixture .var");
            if let Some(tok) = required_tok {
                assert!(
                    var_out.contains(tok),
                    "expected token {tok} in {bam_rel} output"
                );
            }
            let expected_raw = fs::read_to_string(droot.join(oracle_rel)).expect("read oracle");
            let expected_body: String = expected_raw
                .lines()
                .skip(1)
                .collect::<Vec<&str>>()
                .join("\n");
            let mut actual_rows = normalize_var_rows(&var_out, 36);
            let mut expected_rows = normalize_var_rows(&expected_body, 36);
            actual_rows.sort();
            expected_rows.sort();
            assert_eq!(
                actual_rows, expected_rows,
                "SV parity mismatch on {bam_rel}"
            );
        }
    }

    #[test]
    fn deldupvar_boundary_inclusion() {
        ensure_fastas_converted();
        let droot = data_root();
        let bam = droot
            .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/DelOutOfBound.chr1_150578117_150578617.bam")
            .to_string_lossy()
            .into_owned();
        let mut cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg38.fa")
                .to_string_lossy()
                .into_owned(),
            sample: extract_sample_name(&bam, None),
            bam,
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            column_for_chromosome: Some(0),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            region: Some("chr1:150578117-150578617:testbed".to_string()),
            freq: 0.001,
            af_threshold: 0.001,
            deldupvar: false,
            ..Config::default()
        };
        let without_flag = run_core_var(&cfg).expect("without --deldupvar");
        cfg.deldupvar = true;
        let with_flag = run_core_var(&cfg).expect("with --deldupvar");

        let rows_without = normalize_var_rows(&without_flag, 36);
        let rows_with = normalize_var_rows(&with_flag, 36);
        assert_eq!(rows_without.len(), 7);
        assert_eq!(rows_with.len(), 6);
        assert!(
            without_flag.contains("C/-752") && !with_flag.contains("C/-752"),
            "--deldupvar must filter out boundary-spanning deletion C/-752 whose start precedes region start"
        );
    }

    #[test]
    fn nosv_suppresses_structural_variants() {
        ensure_fastas_converted();
        let droot = data_root();
        let bam = droot
            .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.INVSV.bam")
            .to_string_lossy()
            .into_owned();
        let mut cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg19.fa")
                .to_string_lossy()
                .into_owned(),
            sample: "Colo829".to_string(),
            bam,
            threads: 1,
            zero_based: true,
            region: Some("chr20:46841555-46845755:testbed".to_string()),
            freq: 0.001,
            af_threshold: 0.001,
            no_sv: false,
            ..Config::default()
        };
        let vcf_with_sv = run_core_vcf(&cfg).expect("VCF with SV enabled");
        assert!(
            vcf_with_sv.contains("<INV>"),
            "INVSV fixture must emit <INV> when SV calling is enabled"
        );

        cfg.no_sv = true;
        let vcf_nosv = run_core_vcf(&cfg).expect("VCF with --nosv");
        assert!(
            !vcf_nosv.contains("<INV>")
                && !vcf_nosv.contains("<DEL>")
                && !vcf_nosv.contains("<DUP>"),
            "--nosv must suppress <DEL>/<DUP>/<INV> structural variants"
        );
    }
}

// 5. mod filtering_and_pileup (13 tests)
mod filtering_and_pileup {
    use super::*;

    #[test]
    fn min_matches_threshold() {
        let mut cfg_filtered = base_micro_config();
        cfg_filtered.minmatch = 200;
        let vars_filtered = extract_variants(&run_core_vcf(&cfg_filtered).expect("-M 200"));
        assert!(
            vars_filtered.is_empty(),
            "-M 200 should filter out ~80bp reads and produce 0 variants"
        );

        let mut cfg_pass = base_micro_config();
        cfg_pass.minmatch = 50;
        let vars_pass = extract_variants(&run_core_vcf(&cfg_pass).expect("-M 50"));
        assert!(
            !vars_pass.is_empty(),
            "-M 50 should retain ~80bp reads and emit variants"
        );
    }

    #[test]
    fn output_mapq_cap() {
        ensure_fastas_converted();
        let droot = data_root();
        let cfg_unfiltered = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/GRCh37.fa")
                .to_string_lossy()
                .into_owned(),
            bam: droot
                .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/NA24385_chr_20_20156208_20234556.bam")
                .to_string_lossy()
                .into_owned(),
            sample: "NA24385".to_string(),
            region: Some("20:20156208-20157000".to_string()),
            out_mapq: 0.0,
            threads: 1,
            ..Config::default()
        };
        let mut cfg_filtered = cfg_unfiltered.clone();
        cfg_filtered.out_mapq = 65.0;

        let vars_unfiltered = extract_variants(&run_core_vcf(&cfg_unfiltered).expect("-O 0"));
        let vars_filtered = extract_variants(&run_core_vcf(&cfg_filtered).expect("-O 65"));
        assert!(
            vars_filtered.len() < vars_unfiltered.len(),
            "-O 65 must filter low-MAPQ candidates relative to -O 0"
        );
        assert!(
            !vars_filtered.is_empty(),
            "-O 65 must still retain high-AF / high-MAPQ variants"
        );
    }

    #[test]
    fn move_indels_to_3_prime() {
        let mut cfg = base_micro_config();
        cfg.region = Some("chr2:110-230".to_string());
        cfg.three_prime = true;
        let vars = extract_variants(&run_core_vcf(&cfg).expect("-3"));
        assert!(
            !vars.is_empty(),
            "-3 must call variants on chr2 indel region"
        );
    }

    #[test]
    fn reference_extension_window() {
        let root = workspace_root();
        let mut cfg = base_micro_config();
        cfg.region = None;
        cfg.bed = Some(
            root.join("fixtures/micro/micro.bed")
                .to_string_lossy()
                .into_owned(),
        );
        cfg.zero_based = true;
        cfg.col_chr = 0;
        cfg.col_start = 1;
        cfg.col_end = 2;
        cfg.col_seg_start = 1;
        cfg.col_seg_end = 2;
        cfg.col_gene = Some(3);
        cfg.column_for_chromosome = Some(0);
        cfg.omit_end_tag = true;
        cfg.extension = 50;
        let vars = extract_variants(&run_core_vcf(&cfg).expect("-x 50"));
        assert!(
            !vars.is_empty(),
            "-x 50 must expand regions and call variants"
        );
    }

    #[test]
    fn duplicate_read_filtering() {
        let root = workspace_root();
        let cfg_nodedup = Config {
            reference: root
                .join("fixtures/pr402/ref.fa")
                .to_string_lossy()
                .into_owned(),
            bam: root
                .join("fixtures/pr402/pr402.bam")
                .to_string_lossy()
                .into_owned(),
            sample: "test".to_string(),
            region: Some("chr7:1-140".to_string()),
            freq: 0.01,
            af_threshold: 0.01,
            omit_end_tag: false,
            dedup: false,
            threads: 1,
            ..Config::default()
        };
        let vars_nodedup = extract_variants(&run_core_vcf(&cfg_nodedup).expect("without -t"));
        assert_eq!(vars_nodedup.len(), 1);
        assert!(vars_nodedup[0].contains("DP=30;END=30;VD=30"));

        let mut cfg_dedup = cfg_nodedup.clone();
        cfg_dedup.dedup = true;
        let vars_dedup = extract_variants(&run_core_vcf(&cfg_dedup).expect("with -t"));
        assert_eq!(
            vars_dedup.len(),
            0,
            "with -t, duplicate reads at chr7:1 are removed and 0 variants called"
        );

        // Issue G Part 1: when -t removes duplicate reads but a variant still passes
        // (e.g. 2 identical single-end reads at pos 20 + 1 at pos 21 + 1 at pos 22, all
        // carrying SNV at pos 40 -> 1 of 4 reads is a duplicate, leaving DP=3, VD=3),
        // VarDictJava v1.8.3 snapshots totalReads/duplicateReads before the read loop
        // (CigarParser.java:78,162-163), so DUPRATE is always "0" in .var and VCF.
        let tmp = tempfile::tempdir().expect("tempdir");
        let fa = write_ref(&tmp, "chr1", 200, true);
        let bam_path = tmp.path().join("dedup_snv.bam");
        let make_read_with_snv_at_40 = |name: &str, start: usize| -> RecordBuf {
            let snv_offset = 40 - start;
            let mut seq = vec![b'A'; 50];
            seq[snv_offset] = b'C';
            let cigar: Cigar = [Op::new(Kind::Match, 50)].into_iter().collect();
            RecordBuf::builder()
                .set_name(name)
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start).unwrap())
                .set_cigar(cigar)
                .set_sequence(seq.into())
                .set_quality_scores(vec![40u8; 50].into())
                .build()
        };
        let recs = vec![
            make_read_with_snv_at_40("r1_dup_a", 20),
            make_read_with_snv_at_40("r1_dup_b", 20),
            make_read_with_snv_at_40("r2_uniq", 21),
            make_read_with_snv_at_40("r3_uniq", 22),
        ];
        indexed_bam(&bam_path, "chr1", 200, &recs).expect("write dedup_snv.bam");
        let cfg_dup_rate = Config {
            reference: fa.to_string_lossy().into_owned(),
            bam: bam_path.to_string_lossy().into_owned(),
            sample: "S".to_string(),
            region: Some("chr1:1-150".to_string()),
            freq: 0.01,
            af_threshold: 0.01,
            dedup: true,
            threads: 1,
            ..Config::default()
        };
        let var_out = run_core_var(&cfg_dup_rate).expect(".var with -t");
        let var_cols: Vec<&str> = var_out
            .lines()
            .next()
            .expect("one .var row")
            .split('\t')
            .collect();
        assert_eq!(var_cols[7], "3", "DP must reflect 1 removed duplicate read");
        assert_eq!(var_cols[34], "0", "36-col .var Duprate must be 0");

        let vcf_out = extract_variants(&run_core_vcf(&cfg_dup_rate).expect("VCF with -t"));
        assert_eq!(vcf_out.len(), 1);
        assert!(
            vcf_out[0].contains(";DUPRATE=0;"),
            "VCF INFO must contain DUPRATE=0, got: {}",
            vcf_out[0]
        );
    }

    #[test]
    fn region_segment_column_indexing() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let bed_a = tmp.path().join("se_a.bed");
        let bed_b = tmp.path().join("se_b.bed");
        fs::write(
            &bed_a,
            "chr1\t40\t120\tmicro_chr1\t0\t+\t40\t120\t0\t40\t80\tmicro_chr1\n",
        )
        .expect("write bed_a");
        fs::write(
            &bed_b,
            "chr1\t40\t120\tmicro_chr1\t0\t+\t40\t120\t0\t90\t120\tmicro_chr1\n",
        )
        .expect("write bed_b");

        let mut cfg_a = base_micro_config();
        cfg_a.region = None;
        cfg_a.bed = Some(bed_a.to_string_lossy().into_owned());
        cfg_a.zero_based = true;
        cfg_a.col_chr = 0;
        cfg_a.col_start = 1;
        cfg_a.col_end = 2;
        cfg_a.col_seg_start = 9;
        cfg_a.col_seg_end = 10;
        cfg_a.col_gene = Some(3);
        cfg_a.column_for_chromosome = Some(0);
        cfg_a.omit_end_tag = true;
        cfg_a.extension = 0;

        let mut cfg_b = cfg_a.clone();
        cfg_b.bed = Some(bed_b.to_string_lossy().into_owned());

        let vars_a = extract_variants(&run_core_vcf(&cfg_a).expect("segment 40..80"));
        let vars_b = extract_variants(&run_core_vcf(&cfg_b).expect("segment 90..120"));
        assert_eq!(vars_a.len(), 0, "segment 40..80 must exclude chr1:101");
        assert_eq!(vars_b.len(), 1, "segment 90..120 must include chr1:101");
        assert!(vars_b[0].contains("chr1\t101\t"));
    }

    #[test]
    fn unique_first_alignment_dedup() {
        ensure_fastas_converted();
        let droot = data_root();
        let bam = droot
            .join(
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/L861Q.bam",
            )
            .to_string_lossy()
            .into_owned();
        let mut cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg19.fa")
                .to_string_lossy()
                .into_owned(),
            sample: extract_sample_name(&bam, None),
            bam,
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            column_for_chromosome: Some(0),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            region: Some("chr7:55259400-55259600:testbed".to_string()),
            unique_mode_second_in_pair_enabled: true,
            ..Config::default()
        };
        let var_un = run_core_var(&cfg).expect("--UN .var");
        cfg.unique_mode_second_in_pair_enabled = false;
        let var_default = run_core_var(&cfg).expect("default .var");
        assert!(!var_un.is_empty());
        assert_ne!(
            var_un, var_default,
            "--UN must deduplicate overlapping mate pairs relative to default"
        );
    }

    #[test]
    fn pileup_mode_all_positions() {
        let mut cfg = base_micro_config();
        cfg.region = Some("chr1:100-105".to_string());
        cfg.pileup = true;
        let var_out = run_core_var(&cfg).expect("-p .var");
        let rows: Vec<&str> = var_out.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(
            rows.len(),
            6,
            "-p on chr1:100-105 (1-based) must emit all 6 positions"
        );
    }

    #[test]
    fn splice_junction_counting() {
        ensure_fastas_converted();
        let droot = data_root();
        let bam = droot
            .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam")
            .to_string_lossy()
            .into_owned();
        let cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg38.fa")
                .to_string_lossy()
                .into_owned(),
            sample: extract_sample_name(&bam, None),
            bam,
            threads: 1,
            zero_based: true,
            region: Some("chr8:90534611-90620070:testbed".to_string()),
            freq: 0.1,
            af_threshold: 0.1,
            splice: true,
            ..Config::default()
        };
        let splice_out = run_core_var(&cfg).expect("-i splice table");
        let rows: Vec<&str> = splice_out
            .lines()
            .filter(|l| !l.trim().is_empty())
            .collect();
        assert_eq!(rows.len(), 16, "-i must emit 16 splice junction rows");
    }

    #[test]
    fn amplicon_edge_assignment() {
        let tmp = tempfile::tempdir().expect("tempdir");
        let amp_bed = tmp.path().join("amp.bed");
        fs::write(&amp_bed, "chr1\t81\t140\tmicro_amp\t0\t+\t81\t140\n").expect("write amp bed");
        let mut cfg = base_micro_config();
        cfg.region = None;
        cfg.bed = Some(amp_bed.to_string_lossy().into_owned());
        cfg.zero_based = true;
        cfg.col_chr = 0;
        cfg.col_start = 1;
        cfg.col_end = 2;
        cfg.col_seg_start = 1;
        cfg.col_seg_end = 2;
        cfg.col_gene = Some(3);
        cfg.column_for_chromosome = Some(0);
        cfg.omit_end_tag = true;
        cfg.amplicon = Some("30:0.5".to_string());
        cfg.freq = 0.01;
        cfg.af_threshold = 0.01;
        cfg.min_reads = 1;

        let vcf_out = run_core_vcf(&cfg).expect("amplicon VCF");
        assert!(
            vcf_out.contains("GDAMP=1") || vcf_out.contains("TLAMP=1"),
            "amplicon output must contain GDAMP/TLAMP annotations"
        );
    }

    #[test]
    fn paired_amplicon_vcf_rejected_while_intermediate_var_runs_first_bam() {
        let dir = TempDir::new().expect("tempdir");
        let fa = write_ref(&dir, "chr1", 500, true);

        // First BAM carries an SNV at chr1:100 (50bp at 80..129); second BAM carries an SNV at
        // chr1:102 (50bp at 82..131) so both satisfy auto-detected `10:0.95` bounds on `80..131`.
        let first_bam = dir.path().join("tumor_amp.bam");
        let second_bam = dir.path().join("normal_amp.bam");
        indexed_bam(
            &first_bam,
            "chr1",
            500,
            &[snv_record("t0", 80), snv_record("t1", 80)],
        )
        .expect("write first_bam");
        indexed_bam(
            &second_bam,
            "chr1",
            500,
            &[snv_record("n0", 82), snv_record("n1", 82)],
        )
        .expect("write second_bam");

        let amp_bed = dir.path().join("paired_amp.bed");
        fs::write(&amp_bed, "chr1\t79\t131\tAMP1\t0\t+\t80\t131\n").expect("write 8-col amp bed");

        let single_cfg = Config {
            reference: fa.to_string_lossy().into_owned(),
            sample: "TUMOR_AMP".to_string(),
            bam: first_bam.to_string_lossy().into_owned(),
            bed: Some(amp_bed.to_string_lossy().into_owned()),
            zero_based: true,
            zero_based_defined: true,
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            column_for_chromosome: Some(0),
            omit_end_tag: true,
            amplicon: Some("30:0.5".to_string()),
            freq: 0.01,
            af_threshold: 0.01,
            ..Config::default()
        };
        let expected_single_var = run_core_var(&single_cfg).expect("single-BAM amplicon .var");
        let single_rows = normalize_var_rows(&expected_single_var, 38);
        assert_eq!(single_rows.len(), 1);
        assert!(expected_single_var.contains("\tchr1\t100\t100\tA\tC\t"));

        // Exercise both explicit `-a 30:0.5` and 8-column BED amplicon auto-detection (`amplicon: None`).
        for amp_opt in [Some("30:0.5".to_string()), None] {
            let mut paired_cfg = single_cfg.clone();
            paired_cfg.bam = format!("{}|{}", first_bam.display(), second_bam.display());
            paired_cfg.amplicon = amp_opt.clone();

            // 1. VCF mode (`emit_intermediate == false`) rejects paired + amplicon (#378)
            let vcf_err = run_core_vcf(&paired_cfg)
                .expect_err("paired + amplicon VCF must fail with Unsupported (#378)");
            match vcf_err {
                VardictError::Unsupported(msg) => {
                    assert!(
                        msg.contains("not supported in VCF mode (#378)"),
                        "unexpected #378 message for amplicon={amp_opt:?}: {msg}"
                    );
                }
                other => panic!("expected VardictError::Unsupported for #378, got: {other:?}"),
            }

            // 2. Intermediate `.var` mode (`emit_intermediate == true`) succeeds on the first BAM only
            let paired_var =
                run_core_var(&paired_cfg).expect("paired + amplicon intermediate .var must run");
            let paired_rows = normalize_var_rows(&paired_var, 38);
            assert_eq!(paired_rows.len(), 1);
            assert!(
                paired_var.contains("\tchr1\t100\t100\tA\tC\t"),
                "first BAM SNV at chr1:100 must be called in .var mode"
            );
            assert!(
                !paired_var.contains("\tchr1\t102\t102\t"),
                "second BAM SNV at chr1:102 must NOT be called in AmpliconMode (#378)"
            );
            assert_eq!(
                paired_var, expected_single_var,
                "paired amplicon .var must match single first-BAM amplicon .var byte-for-byte for amplicon={amp_opt:?}"
            );
        }

        // 3. Direct library `run_with_writer` validates `Config::sample_regexp` and resolves
        //    samples across single-sample, BED, `-R`, paired VCF, and paired `.var` while
        //    preserving `Config::default()` (`"SAMPLE"`) when `sample_regexp: None`.
        let t_cohort_dir = dir.path().join("cohort").join("T_LIB");
        let n_cohort_dir = dir.path().join("cohort").join("N_LIB");
        fs::create_dir_all(&t_cohort_dir).expect("mkdir T_LIB");
        fs::create_dir_all(&n_cohort_dir).expect("mkdir N_LIB");
        let t_cohort_bam = t_cohort_dir.join("t.bam");
        let n_cohort_bam = n_cohort_dir.join("n.bam");
        fs::copy(&first_bam, &t_cohort_bam).expect("copy t_cohort_bam");
        fs::copy(
            first_bam.with_extension("bam.bai"),
            t_cohort_bam.with_extension("bam.bai"),
        )
        .expect("copy t_cohort_bam.bai");
        fs::copy(&second_bam, &n_cohort_bam).expect("copy n_cohort_bam");
        fs::copy(
            second_bam.with_extension("bam.bai"),
            n_cohort_bam.with_extension("bam.bai"),
        )
        .expect("copy n_cohort_bam.bai");
        let plain_bed = dir.path().join("plain.bed");
        fs::write(&plain_bed, "chr1\t79\t131\tGENE1\n").expect("write plain.bed");

        let default_sample_cfg = Config {
            reference: fa.to_string_lossy().into_owned(),
            bam: t_cohort_bam.to_string_lossy().into_owned(),
            region: Some("chr1:80-131".to_string()),
            threads: 1,
            ..Config::default()
        };
        let mut default_buf = Vec::new();
        vardictrust_core::run_with_writer(&default_sample_cfg, &mut default_buf)
            .expect("run_with_writer default sample");
        assert!(
            String::from_utf8_lossy(&default_buf)
                .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"),
            "Config::default() with sample_regexp: None must keep SAMPLE"
        );

        let mut bad_regex_cfg = default_sample_cfg.clone();
        bad_regex_cfg.sample_regexp = Some("([unclosed".to_string());
        let mut bad_buf = Vec::new();
        let bad_err = vardictrust_core::run_with_writer(&bad_regex_cfg, &mut bad_buf)
            .expect_err("malformed Config::sample_regexp must fail fast");
        assert!(
            matches!(bad_err, VardictError::Parse(_)),
            "expected VardictError::Parse for malformed sample_regexp, got {bad_err:?}"
        );

        for pat in [r"/cohort\/([^/]+)\//", r"//cohort/([^/]+)//"] {
            let mut single_r_cfg = default_sample_cfg.clone();
            single_r_cfg.sample_regexp = Some(pat.to_string());
            let mut out_r = Vec::new();
            vardictrust_core::run_with_writer(&single_r_cfg, &mut out_r)
                .expect("single -R sample_regexp");
            assert!(
                String::from_utf8_lossy(&out_r)
                    .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT_LIB\n"),
                "single -R failed for {pat}"
            );

            let single_bed_cfg = Config {
                reference: fa.to_string_lossy().into_owned(),
                bam: t_cohort_bam.to_string_lossy().into_owned(),
                bed: Some(plain_bed.to_string_lossy().into_owned()),
                zero_based: true,
                zero_based_defined: true,
                col_chr: 0,
                col_start: 1,
                col_end: 2,
                col_seg_start: 1,
                col_seg_end: 2,
                col_gene: Some(3),
                column_for_chromosome: Some(0),
                sample_regexp: Some(pat.to_string()),
                threads: 1,
                ..Config::default()
            };
            let mut out_bed = Vec::new();
            vardictrust_core::run_with_writer(&single_bed_cfg, &mut out_bed)
                .expect("single BED sample_regexp");
            assert!(
                String::from_utf8_lossy(&out_bed)
                    .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT_LIB\n"),
                "single BED failed for {pat}"
            );

            let mut paired_vcf_cfg = single_bed_cfg.clone();
            paired_vcf_cfg.bam = format!("{}|{}", t_cohort_bam.display(), n_cohort_bam.display());
            let mut out_paired_vcf = Vec::new();
            vardictrust_core::run_with_writer(&paired_vcf_cfg, &mut out_paired_vcf)
                .expect("paired VCF sample_regexp");
            assert!(
                String::from_utf8_lossy(&out_paired_vcf).contains(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT_LIB\tN_LIB\n"
                ),
                "paired VCF failed for {pat}"
            );

            let paired_var_cfg = Config {
                reference: fa.to_string_lossy().into_owned(),
                bam: format!("{}|{}", t_cohort_bam.display(), n_cohort_bam.display()),
                region: Some("chr1:80-131".to_string()),
                emit_intermediate: true,
                sample_regexp: Some(pat.to_string()),
                threads: 1,
                ..Config::default()
            };
            let mut out_paired_var = Vec::new();
            vardictrust_core::run_with_writer(&paired_var_cfg, &mut out_paired_var)
                .expect("paired -R .var sample_regexp");
            let paired_var_str = String::from_utf8_lossy(&out_paired_var);
            assert!(
                paired_var_str.lines().any(|l| l.starts_with("T_LIB\t")),
                "paired -R .var failed for {pat}: {paired_var_str}"
            );
        }
    }

    #[test]
    fn qratio_threshold_calculation() {
        let mut cfg_pass = base_micro_config();
        cfg_pass.freq = 0.01;
        cfg_pass.af_threshold = 0.01;
        cfg_pass.qratio = 0.0;
        let vars_pass = extract_variants(&run_core_vcf(&cfg_pass).expect("-o 0.0"));
        assert!(!vars_pass.is_empty(), "-o 0.0 must allow variant calls");

        let mut cfg_filter = cfg_pass.clone();
        cfg_filter.qratio = 100.0;
        let vars_filter = extract_variants(&run_core_vcf(&cfg_filter).expect("-o 100.0"));
        assert!(
            vars_filter.is_empty(),
            "-o 100.0 must filter out variants with qratio < 100.0"
        );
    }

    #[test]
    fn realignment_chimeric_adapter_and_crispr_options() {
        let mut cfg = base_micro_config();
        cfg.sample = "CUSTOM_SAMPLE".to_string();
        cfg.fisher = true;
        cfg.column_header = true;
        cfg.debug = true;
        cfg.unique_mode_alignment_enabled = true;
        cfg.include_n = true;
        cfg.freq = 0.01;
        cfg.af_threshold = 0.01;
        cfg.min_reads = 2;
        cfg.min_reads_bias = 1;
        cfg.base_qual = 20.0;
        cfg.map_qual = Some(10);
        cfg.max_mismatches = 5;
        cfg.read_pos = 3.0;
        cfg.extension_nuc = 2;
        cfg.max_region_len = 1000;
        cfg.indel_size = 40;
        cfg.flag_filter = 0x504;
        cfg.sv_read_len = 250;
        cfg.sv_flanking = 80;
        cfg.sv_min_af = 2;
        cfg.sv_min_len = 500;
        cfg.lofreq = 0.02;
        cfg.no_sv = true;
        cfg.chimeric = true;
        cfg.mfreq = 0.2;
        cfg.nmfreq = 0.08;
        cfg.crispr_cutting_site = Some(100);
        cfg.crispr_flanking_size = Some(5);
        cfg.tile_size = Some(500);
        cfg.adaptor = vec!["AGATCGGAAGAGC".to_string()];
        cfg.realign = 0;
        cfg.qratio = 1.5;

        let var_out = run_core_var(&cfg).expect("complex operational options");
        let lines: Vec<&str> = var_out.lines().collect();
        assert!(!lines.is_empty(), "must emit table lines");
        assert!(
            lines[0].starts_with("Sample\t"),
            "with column_header and emit_intermediate, first line must be column header"
        );
        assert!(var_out.contains("CUSTOM_SAMPLE"));

        // E2E forward + reverse `-T` (`trim_bases_after`) 5'-distance read trimming
        let trim_dir = TempDir::new().expect("trim tempdir");
        let trim_fa = write_ref(&trim_dir, "chr1", 500, true);
        let cigar_30m: Cigar = [(30usize, Kind::Match)]
            .iter()
            .map(|&(l, k)| Op::new(k, l))
            .collect();
        let mut trim_recs = Vec::new();
        for idx in 0..4 {
            // Forward 30M @ 100..=129: read_pos_inc=9 (pos 109, dist=9 <= 15) and read_pos_inc=19 (pos 119, dist=19 > 15)
            let mut fwd_seq = vec![b'A'; 30];
            fwd_seq[9] = b'C';
            fwd_seq[19] = b'G';
            trim_recs.push(
                RecordBuf::builder()
                    .set_name(format!("fwd_var_{idx}"))
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(100usize).unwrap())
                    .set_cigar(cigar_30m.clone())
                    .set_sequence(fwd_seq.into())
                    .set_quality_scores(vec![40u8; 30].into())
                    .build(),
            );
            trim_recs.push(
                RecordBuf::builder()
                    .set_name(format!("fwd_ref_{idx}"))
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(100usize).unwrap())
                    .set_cigar(cigar_30m.clone())
                    .set_sequence(vec![b'A'; 30].into())
                    .set_quality_scores(vec![40u8; 30].into())
                    .build(),
            );
            // Reverse 30M @ 200..=229 (5' of sequenced read is right end):
            // read_pos_inc=9 (pos 209, dist=30-9=21 > 15) and read_pos_inc=20 (pos 220, dist=30-20=10 <= 15)
            let mut rev_seq = vec![b'A'; 30];
            rev_seq[9] = b'T';
            rev_seq[20] = b'C';
            trim_recs.push(
                RecordBuf::builder()
                    .set_name(format!("rev_var_{idx}"))
                    .set_flags(noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED)
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(200usize).unwrap())
                    .set_cigar(cigar_30m.clone())
                    .set_sequence(rev_seq.into())
                    .set_quality_scores(vec![40u8; 30].into())
                    .build(),
            );
            trim_recs.push(
                RecordBuf::builder()
                    .set_name(format!("rev_ref_{idx}"))
                    .set_flags(noodles::sam::alignment::record::Flags::REVERSE_COMPLEMENTED)
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(200usize).unwrap())
                    .set_cigar(cigar_30m.clone())
                    .set_sequence(vec![b'A'; 30].into())
                    .set_quality_scores(vec![40u8; 30].into())
                    .build(),
            );
        }
        trim_recs.sort_by_key(|r| r.alignment_start());
        let trim_bam = trim_dir.path().join("trim_e2e.bam");
        indexed_bam(&trim_bam, "chr1", 500, &trim_recs).expect("write trim_e2e.bam");

        let mut untrimmed_cfg = Config {
            reference: trim_fa.to_string_lossy().into_owned(),
            bam: trim_bam.to_string_lossy().into_owned(),
            sample: "TRIM".to_string(),
            region: Some("chr1:90-240".to_string()),
            max_mismatches: 5,
            threads: 1,
            ..Config::default()
        };
        let untrimmed_vcf = run_core_vcf(&untrimmed_cfg).expect("untrimmed vcf");
        for pos in [109, 119, 209, 220] {
            assert!(
                untrimmed_vcf.contains(&format!("chr1\t{pos}\t")),
                "expected untrimmed VCF to contain pos {pos}, got:\n{untrimmed_vcf}"
            );
        }

        untrimmed_cfg.trim_bases_after = 15;
        let trimmed_vcf = run_core_vcf(&untrimmed_cfg).expect("trimmed -T 15 vcf");
        for kept in [109, 220] {
            assert!(
                trimmed_vcf.contains(&format!("chr1\t{kept}\t")),
                "expected -T 15 VCF to retain pos {kept}, got:\n{trimmed_vcf}"
            );
        }
        for dropped in [119, 209] {
            assert!(
                !trimmed_vcf.contains(&format!("chr1\t{dropped}\t")),
                "expected -T 15 VCF to trim pos {dropped}, got:\n{trimmed_vcf}"
            );
        }
    }
}

// 6. mod hostile_inputs (6 tests)
mod hostile_inputs {
    use super::*;

    #[rustfmt::skip]
    #[allow(clippy::type_complexity)]
    const BAD_CIGARS: &[(&str, usize, &[(usize, Kind)])] = &[
        ("del-at-chr-start",   1,   &[(10, Kind::Match), (200, Kind::Deletion), (10, Kind::Match)]),
        ("leading-D",          50,  &[(2, Kind::Deletion), (10, Kind::Match)]),
        ("trailing-D",         70,  &[(10, Kind::Match), (3, Kind::Deletion)]),
        ("giant-ins",          90,  &[(5, Kind::Match), (100, Kind::Insertion), (5, Kind::Match)]),
        ("mixed-mid",          120, &[(1, Kind::Match), (1, Kind::Insertion), (1, Kind::Deletion), (1, Kind::Match)]),
        ("h-before-s",         140, &[(5, Kind::HardClip), (5, Kind::SoftClip), (10, Kind::Match)]),
        ("del-past-ref-end",   985, &[(10, Kind::Match), (500, Kind::Deletion), (10, Kind::Match)]),
    ];

    fn make_read(pos: usize, ops: &[(usize, Kind)]) -> RecordBuf {
        let cigar: Cigar = ops.iter().map(|&(l, k)| Op::new(k, l)).collect();
        let rl = cigar.read_length();
        assert!(rl > 0, "zero-read-length CIGARs are unspecifiable");
        RecordBuf::builder()
            .set_name(format!("r{pos}"))
            .set_flags(noodles::sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(pos).unwrap())
            .set_cigar(cigar)
            .set_sequence(vec![b'A'; rl].into())
            .set_quality_scores(vec![40; rl].into())
            .build()
    }

    fn bed_config(
        dir: &TempDir,
        bed_lines: &str,
        fasta: &std::path::Path,
        bam: &std::path::Path,
    ) -> Config {
        let bed = dir.path().join("case.bed");
        fs::write(&bed, bed_lines).expect("write bed");
        Config {
            reference: fasta.to_string_lossy().into_owned(),
            bam: bam.to_string_lossy().into_owned(),
            sample: "HOSTILE".to_string(),
            bed: Some(bed.to_string_lossy().into_owned()),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            column_for_chromosome: Some(0),
            omit_end_tag: true,
            threads: 1,
            ..Config::default()
        }
    }

    #[test]
    fn empty_bam_produces_empty_output() {
        let dir = TempDir::new().expect("tempdir");
        let empty_bam = dir.path().join("empty.bam");
        indexed_bam(&empty_bam, "chr1", 1000, &[]).expect("write empty bam");
        let fa = write_ref(&dir, "chr1", 1000, true);
        let cfg = bed_config(&dir, "chr1\t100\t200\tG1\n", &fa, &empty_bam);

        let vcf = run_core_vcf(&cfg).expect("empty BAM must succeed cleanly");
        assert!(vcf.contains("#CHROM"));
        assert!(extract_variants(&vcf).is_empty());
    }

    #[test]
    fn malformed_cigar_returns_clean_result() {
        let dir = TempDir::new().expect("tempdir");
        let bad_bam = dir.path().join("bad.bam");
        let recs: Vec<_> = BAD_CIGARS
            .iter()
            .map(|(_, pos, ops)| make_read(*pos, ops))
            .collect();
        indexed_bam(&bad_bam, "chr1", 1000, &recs).expect("write bad cigar bam");
        let fa = write_ref(&dir, "chr1", 1000, true);
        let cfg = bed_config(&dir, "chr1\t50\t999\tG2\n", &fa, &bad_bam);

        let vcf = run_core_vcf(&cfg).expect("malformed CIGARs must not crash or error");
        assert!(vcf.contains("#CHROM"));
    }

    #[test]
    fn out_of_bounds_bed_region_clamps_cleanly() {
        let dir = TempDir::new().expect("tempdir");
        let empty_bam = dir.path().join("empty.bam");
        indexed_bam(&empty_bam, "chr1", 1000, &[]).expect("write empty bam");
        let fa = write_ref(&dir, "chr1", 1000, true);
        let cfg = bed_config(
            &dir,
            "chr1\t300\t100\tBAD\nchr1\t900\t5000\tOOB\n",
            &fa,
            &empty_bam,
        );
        // Either succeeds after clamping OOB region or returns a clean Parse/Other error, never panicking.
        let res = run_core_vcf(&cfg);
        assert!(
            res.is_ok() || matches!(res, Err(VardictError::Parse(_) | VardictError::Other(_))),
            "unexpected result on OOB BED: {res:?}"
        );
    }

    #[test]
    fn missing_fasta_index_returns_error() {
        let dir = TempDir::new().expect("tempdir");
        let empty_bam = dir.path().join("empty.bam");
        indexed_bam(&empty_bam, "chr1", 1000, &[]).expect("write empty bam");
        let fa_no_idx = write_ref(&dir, "chr1", 1000, false);
        let cfg = bed_config(&dir, "chr1\t100\t200\tG4\n", &fa_no_idx, &empty_bam);

        let err = run_core_vcf(&cfg).expect_err("missing .fai must return Err");
        assert!(
            matches!(err, VardictError::Io(_) | VardictError::Other(_)),
            "expected Io or Other error for missing .fai, got: {err:?}"
        );
    }

    #[test]
    fn absent_chromosome_skips_cleanly() {
        let dir = TempDir::new().expect("tempdir");
        let chrz_bam = dir.path().join("chrz.bam");
        indexed_bam(
            &chrz_bam,
            "chrZ",
            500,
            &[make_read(50, &[(15, Kind::Match)])],
        )
        .expect("write chrZ bam");
        let fa = write_ref(&dir, "chr1", 1000, true);
        let cfg = bed_config(&dir, "chrZ\t50\t200\tG5\n", &fa, &chrz_bam);

        let res = run_core_vcf(&cfg);
        assert!(
            res.is_ok() || matches!(res, Err(VardictError::Io(_))),
            "absent chromosome must return clean result without panicking, got: {res:?}"
        );
    }

    #[test]
    fn validation_error_rate_threshold_enforced() {
        let root = workspace_root();
        let orig_bam = root.join("fixtures/pr402/pr402.bam");
        let ref_fa = root.join("fixtures/pr402/ref.fa");
        let bed_path = root.join("fixtures/pr402/pr402.bed");

        let mut reader = noodles::bam::io::reader::Builder
            .build_from_path(&orig_bam)
            .expect("open fixture bam");
        let header = reader.read_header().expect("read bam header");

        let mut records: Vec<RecordBuf> = Vec::new();
        let mut bad_template: Option<RecordBuf> = None;
        for rec_res in reader.record_bufs(&header) {
            let rec = rec_res.expect("read record");
            if rec.name().map(|n| n == "bad_read").unwrap_or(false) {
                bad_template = Some(rec.clone());
            }
            records.push(rec);
        }
        let template = bad_template.expect("fixture must contain bad_read");
        for i in 1..=10 {
            let mut dup = template.clone();
            *dup.name_mut() = Some(format!("bad_read_{i}").into());
            records.push(dup);
        }
        records.sort_by_key(|r| r.alignment_start().map(|p| p.get()).unwrap_or(usize::MAX));

        let dir = TempDir::new().expect("tempdir");
        let bam_11 = dir.path().join("sabotaged_11.bam");
        indexed_bam(&bam_11, "chr7", 600, &records).expect("write indexed 11-bad-read BAM");

        let cfg1 = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: orig_bam.to_string_lossy().into_owned(),
            sample: "SAMPLE".to_string(),
            bed: Some(bed_path.to_string_lossy().into_owned()),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            column_for_chromosome: Some(0),
            omit_end_tag: true,
            threads: 1,
            ..Config::default()
        };
        let mut cfg11 = cfg1.clone();
        cfg11.bam = bam_11.to_string_lossy().into_owned();
        cfg11.exception_counter = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));

        let out1 = run_core_vcf(&cfg1).expect("1 bad read must succeed");
        let out11 = run_core_vcf(&cfg11).expect("11 bad reads (>10 threshold) must succeed");

        assert_eq!(
            cfg1.exception_counter.load(Ordering::Relaxed),
            1,
            "expected 1 caught CIGAR modification exception on 1-bad-read BAM"
        );
        assert_eq!(
            cfg11.exception_counter.load(Ordering::Relaxed),
            11,
            "expected 11 caught CIGAR modification exceptions (>10 ceiling) on 11-bad-read BAM"
        );
        assert_eq!(
            out1, out11,
            "VCF output must be byte-identical between 1 bad read and 11 bad reads"
        );
    }
}

// 7. mod concurrency (4 tests)
mod concurrency {
    use super::*;

    fn snv_pair(pos: usize) -> Vec<RecordBuf> {
        (0..2)
            .map(|i| snv_record(&format!("s{pos}r{i}"), pos - 20))
            .collect()
    }

    fn insertion_pair(pos: usize) -> Vec<RecordBuf> {
        let cigar: Cigar = [
            (40usize, Kind::Match),
            (10, Kind::Insertion),
            (10, Kind::Match),
        ]
        .iter()
        .map(|&(l, k)| Op::new(k, l))
        .collect();
        (0..2)
            .map(|i| {
                RecordBuf::builder()
                    .set_name(format!("i{pos}r{i}"))
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(pos - 20).unwrap())
                    .set_cigar(cigar.clone())
                    .set_sequence(
                        [vec![b'A'; 40], vec![b'T'; 10], vec![b'A'; 10]]
                            .concat()
                            .into(),
                    )
                    .set_quality_scores(vec![40u8; 60].into())
                    .build()
            })
            .collect()
    }

    fn deletion_pair(pos: usize, del_len: usize) -> Vec<RecordBuf> {
        let cigar: Cigar = [
            (20usize, Kind::Match),
            (del_len, Kind::Deletion),
            (30, Kind::Match),
        ]
        .iter()
        .map(|&(l, k)| Op::new(k, l))
        .collect();
        (0..2)
            .map(|i| {
                RecordBuf::builder()
                    .set_name(format!("d{pos}r{i}"))
                    .set_flags(noodles::sam::alignment::record::Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::try_from(pos - 20).unwrap())
                    .set_cigar(cigar.clone())
                    .set_sequence(vec![b'A'; 50].into())
                    .set_quality_scores(vec![40u8; 50].into())
                    .build()
            })
            .collect()
    }

    #[test]
    fn vcf_output_byte_identical_across_thread_counts() {
        let dir = TempDir::new().expect("tempdir");
        let fa = write_ref(&dir, "chr1", 2000, true);

        let mut recs = Vec::new();
        for pos in (150..=1000).step_by(80) {
            recs.extend(snv_pair(pos));
        }
        recs.extend(insertion_pair(600));
        recs.sort_by_key(|r| r.alignment_start());
        let bam = dir.path().join("det.bam");
        indexed_bam(&bam, "chr1", 2000, &recs).expect("write indexed bam");

        let bed = dir.path().join("det.bed");
        fs::write(&bed, "chr1\t100\t1100\tDET\n").expect("write bed");

        let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
        let mut counts: Vec<usize> = [1usize, 2, 8].iter().map(|t| (*t).min(cores)).collect();
        counts.dedup();

        let runs: Vec<Vec<u8>> = counts
            .iter()
            .map(|&t| {
                let bam_str = bam.to_string_lossy().into_owned();
                let cfg = Config {
                    reference: fa.to_string_lossy().into_owned(),
                    sample: extract_sample_name(&bam_str, None),
                    bam: bam_str,
                    bed: Some(bed.to_string_lossy().into_owned()),
                    col_chr: 0,
                    col_start: 1,
                    col_end: 2,
                    col_seg_start: 1,
                    col_seg_end: 2,
                    col_gene: Some(3),
                    column_for_chromosome: Some(0),
                    omit_end_tag: true,
                    threads: t,
                    tile_size: Some(300),
                    emit_intermediate: false,
                    ..Config::default()
                };
                let mut buf = Vec::new();
                vardictrust_core::run_with_writer(&cfg, &mut buf).expect("run_with_writer");
                buf
            })
            .collect();

        let s = String::from_utf8_lossy(&runs[0]);
        assert!(s.contains("CHROM"), "no VCF header");
        let rows = s.lines().filter(|l| l.starts_with("chr1\t")).count();
        assert!(
            rows >= 12,
            "expected >=12 variant rows spanning multiple tiles, got {rows}"
        );
        for pos in [150, 950] {
            assert!(
                s.contains(&format!("chr1\t{pos}\t")),
                "tile-spanning variant @{pos} missing"
            );
        }
        for pair in runs.windows(2) {
            assert_eq!(pair[0], pair[1], "VCF bytes differ across thread counts");
        }
    }

    #[test]
    fn amplicon_multithreaded_deterministic() {
        let dir = TempDir::new().expect("tempdir");
        let fa = write_ref(&dir, "chr1", 2000, true);

        let mut recs = Vec::new();
        for pos in [150usize, 450, 750, 1050] {
            recs.extend(snv_pair(pos));
        }
        recs.sort_by_key(|r| r.alignment_start());
        let bam = dir.path().join("amp_det.bam");
        indexed_bam(&bam, "chr1", 2000, &recs).expect("write indexed bam");

        // 4 non-overlapping amplicon segments matching the 50bp reads at pos-20..pos+29
        let bed = dir.path().join("amp_det.bed");
        fs::write(
            &bed,
            "chr1\t130\t179\tAMP1\t0\t+\t135\t174\n\
             chr1\t430\t479\tAMP2\t0\t+\t435\t474\n\
             chr1\t730\t779\tAMP3\t0\t+\t735\t774\n\
             chr1\t1030\t1079\tAMP4\t0\t+\t1035\t1074\n",
        )
        .expect("write amplicon bed");

        let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
        let mut counts: Vec<usize> = [1usize, 2, 4].iter().map(|t| (*t).min(cores)).collect();
        counts.dedup();

        for emit_intermediate in [true, false] {
            let runs: Vec<Vec<u8>> = counts
                .iter()
                .map(|&t| {
                    let bam_str = bam.to_string_lossy().into_owned();
                    let cfg = Config {
                        reference: fa.to_string_lossy().into_owned(),
                        sample: extract_sample_name(&bam_str, None),
                        bam: bam_str,
                        bed: Some(bed.to_string_lossy().into_owned()),
                        col_chr: 0,
                        col_start: 1,
                        col_end: 2,
                        col_seg_start: 1,
                        col_seg_end: 2,
                        col_gene: Some(3),
                        column_for_chromosome: Some(0),
                        amplicon: Some("10:0.95".to_string()),
                        omit_end_tag: true,
                        threads: t,
                        emit_intermediate,
                        ..Config::default()
                    };
                    let mut buf = Vec::new();
                    vardictrust_core::run_with_writer(&cfg, &mut buf).expect("run_with_writer");
                    buf
                })
                .collect();

            let s = String::from_utf8_lossy(&runs[0]);
            let data_rows = s.lines().filter(|l| !l.starts_with('#')).count();
            assert!(
                data_rows >= 4,
                "expected >=4 amplicon variant rows (emit_intermediate={emit_intermediate}), got {data_rows}"
            );
            for pair in runs.windows(2) {
                assert_eq!(
                    pair[0], pair[1],
                    "amplicon output differs across thread counts (emit_intermediate={emit_intermediate})"
                );
            }
        }
    }

    #[test]
    fn tiled_vs_untiled_var_and_vcf_conservation_simple_and_somatic() {
        let dir = TempDir::new().expect("tempdir");
        let fa = write_ref(&dir, "chr1", 2000, true);

        let mut tumor_recs = Vec::new();
        for pos in (150..=1000).step_by(80) {
            tumor_recs.extend(snv_pair(pos));
        }
        // Boundary-straddling 10-bp insertions (`pos + 19` = 250, 350, 600) across
        // tile_size = 150 (250, 400, 550, 700, 850, 1000), 250 (350, 600, 850), and 500 (600),
        // alongside deletion_pair(500, 3).
        tumor_recs.extend(insertion_pair(231));
        tumor_recs.extend(insertion_pair(331));
        tumor_recs.extend(insertion_pair(581));
        tumor_recs.extend(deletion_pair(500, 3));
        tumor_recs.sort_by_key(|r| r.alignment_start());
        let tumor_bam = dir.path().join("tumor.bam");
        indexed_bam(&tumor_bam, "chr1", 2000, &tumor_recs).expect("write tumor bam");

        let ref_cigar: Cigar = [(50usize, Kind::Match)]
            .iter()
            .map(|&(l, k)| Op::new(k, l))
            .collect();
        let mut normal_recs = Vec::new();
        for pos in (150..=1000).step_by(80).chain([231, 331, 500, 581]) {
            for i in 0..2 {
                normal_recs.push(
                    RecordBuf::builder()
                        .set_name(format!("n{pos}r{i}"))
                        .set_flags(noodles::sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(0)
                        .set_alignment_start(Position::try_from(pos - 20).unwrap())
                        .set_cigar(ref_cigar.clone())
                        .set_sequence(vec![b'A'; 50].into())
                        .set_quality_scores(vec![40u8; 50].into())
                        .build(),
                );
            }
        }
        normal_recs.extend(snv_pair(310));
        normal_recs.sort_by_key(|r| r.alignment_start());
        let normal_bam = dir.path().join("normal.bam");
        indexed_bam(&normal_bam, "chr1", 2000, &normal_recs).expect("write normal bam");

        let bed = dir.path().join("multi_tile.bed");
        fs::write(&bed, "chr1\t100\t1100\tDET\n").expect("write bed");

        let base_cfg = Config {
            reference: fa.to_string_lossy().into_owned(),
            sample: "TUMOR".to_string(),
            bam: tumor_bam.to_string_lossy().into_owned(),
            bed: Some(bed.to_string_lossy().into_owned()),
            col_chr: 0,
            col_start: 1,
            col_end: 2,
            col_seg_start: 1,
            col_seg_end: 2,
            col_gene: Some(3),
            column_for_chromosome: Some(0),
            omit_end_tag: true,
            threads: 1,
            ..Config::default()
        };

        // Compute the 4 untiled (tile_size = None) baselines once before the tile-size loop.
        let mut simple_untiled = base_cfg.clone();
        simple_untiled.tile_size = None;
        let mut untiled_simple_var = normalize_var_rows(
            &run_core_var(&simple_untiled).expect("simple untiled .var"),
            36,
        );
        assert!(
            untiled_simple_var.len() >= 15,
            "expected >=15 simple .var rows across tiles, got {}",
            untiled_simple_var.len()
        );
        untiled_simple_var.sort();
        let untiled_simple_vcf =
            extract_variants(&run_core_vcf(&simple_untiled).expect("simple untiled .vcf"));

        let somatic_bam = format!("{}|{}", tumor_bam.display(), normal_bam.display());
        let mut somatic_untiled = base_cfg.clone();
        somatic_untiled.bam = somatic_bam;
        somatic_untiled.sample = "TUMOR|NORMAL".to_string();
        somatic_untiled.tile_size = None;
        let mut untiled_somatic_var = normalize_var_rows(
            &run_core_var(&somatic_untiled).expect("somatic untiled .var"),
            55,
        );
        assert!(
            untiled_somatic_var.len() >= 15,
            "expected >=15 somatic .var rows across tiles, got {}",
            untiled_somatic_var.len()
        );
        untiled_somatic_var.sort();
        let untiled_somatic_vcf =
            extract_variants(&run_core_vcf(&somatic_untiled).expect("somatic untiled .vcf"));

        for tile_size in [150u32, 250, 500] {
            let mut simple_tiled = simple_untiled.clone();
            simple_tiled.tile_size = Some(tile_size);

            let mut simple_var_tiled =
                normalize_var_rows(&run_core_var(&simple_tiled).expect("simple tiled .var"), 36);
            simple_var_tiled.sort();
            assert_eq!(
                simple_var_tiled, untiled_simple_var,
                "simple .var row set must be identical between tiled ({tile_size}bp) and untiled runs"
            );

            let simple_vcf_tiled =
                extract_variants(&run_core_vcf(&simple_tiled).expect("simple tiled .vcf"));
            assert_eq!(
                simple_vcf_tiled, untiled_simple_vcf,
                "simple .vcf variant lines must be identical between tiled ({tile_size}bp) and untiled runs"
            );

            let mut somatic_tiled = somatic_untiled.clone();
            somatic_tiled.tile_size = Some(tile_size);

            let mut somatic_var_tiled = normalize_var_rows(
                &run_core_var(&somatic_tiled).expect("somatic tiled .var"),
                55,
            );
            somatic_var_tiled.sort();
            assert_eq!(
                somatic_var_tiled, untiled_somatic_var,
                "somatic .var row set must be identical between tiled ({tile_size}bp) and untiled runs"
            );

            let somatic_vcf_tiled =
                extract_variants(&run_core_vcf(&somatic_tiled).expect("somatic tiled .vcf"));
            assert_eq!(
                somatic_vcf_tiled, untiled_somatic_vcf,
                "somatic .vcf variant lines must be identical between tiled ({tile_size}bp) and untiled runs"
            );
        }

        // Somatic RefSiteMode (`CompactOut` + `KeepWhere` + `hydrate_compact_ref_site` +
        // `combine_analysis` `DropAll`) E2E equivalence against `KeepAll` (`debug = true`):
        //   1. Pure ref-only tumor site right after a `p - 1` insertion (p = 250 after +GG at 249)
        //      + normal variant at p = 250, plus pure ref-only tumor site at q = 360 hydrated via
        //      `hydrate_compact_ref_site`.
        //   2. BAM1 sub-threshold allele (AF = 1/40 = 0.025 < freq = 0.05) + BAM2 germline/LOH call
        //      at r = 470.
        //   3. `combine_analysis` low-AF rescue (`combined AF = 5/120 <= freq = 0.05`) on `-12`
        //      deletion at d = 600 (3 explicit `25M12D25M` reads in tumor + 2 `25M6S` soft-clips
        //      in normal rescued as Germline).
        let mut lcg = 0x1234_5678u32;
        let ref_bytes: Vec<u8> = (0..1200)
            .map(|_| {
                lcg = lcg.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[((lcg >> 24) & 3) as usize]
            })
            .collect();
        let ref_fa = dir.path().join("somatic_refmode.fa");
        fs::write(
            &ref_fa,
            format!(">chr1\n{}\n", std::str::from_utf8(&ref_bytes).unwrap()),
        )
        .expect("write somatic_refmode.fa");
        fs::write(
            ref_fa.with_extension("fa.fai"),
            "chr1\t1200\t6\t1200\t1201\n",
        )
        .expect("write somatic_refmode.fa.fai");

        let slice_1b =
            |pos1: usize, len: usize| -> Vec<u8> { ref_bytes[pos1 - 1..pos1 - 1 + len].to_vec() };
        let alt_of = |pos1: usize| -> u8 {
            match ref_bytes[pos1 - 1] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            }
        };
        let mk_read = |name: String, start1: usize, ops: &[(usize, Kind)], seq: Vec<u8>| {
            let rl = seq.len();
            let cigar: Cigar = ops.iter().map(|&(l, k)| Op::new(k, l)).collect();
            RecordBuf::builder()
                .set_name(name)
                .set_flags(noodles::sam::alignment::record::Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(start1).unwrap())
                .set_cigar(cigar)
                .set_sequence(seq.into())
                .set_quality_scores(vec![40u8; rl].into())
                .build()
        };

        let mut t_recs = Vec::new();
        let mut n_recs = Vec::new();

        // Case 1a: p = 250 right after +GG insertion at 249 in tumor; normal SNV at 250
        for i in 0..4 {
            let mut seq = slice_1b(230, 20);
            seq.extend_from_slice(b"GG");
            seq.extend_from_slice(&slice_1b(250, 20));
            t_recs.push(mk_read(
                format!("t_ins249_{i}"),
                230,
                &[(20, Kind::Match), (2, Kind::Insertion), (20, Kind::Match)],
                seq,
            ));
        }
        for i in 0..6 {
            t_recs.push(mk_read(
                format!("t_ref250_{i}"),
                230,
                &[(40, Kind::Match)],
                slice_1b(230, 40),
            ));
        }
        for i in 0..6 {
            let mut seq = slice_1b(230, 40);
            seq[20] = alt_of(250);
            n_recs.push(mk_read(
                format!("n_snv250_{i}"),
                230,
                &[(40, Kind::Match)],
                seq,
            ));
        }
        for i in 0..4 {
            n_recs.push(mk_read(
                format!("n_ref250_{i}"),
                230,
                &[(40, Kind::Match)],
                slice_1b(230, 40),
            ));
        }

        // Case 1b: q = 360 pure ref-only tumor site (hydrated from CompactRefSite) + normal SNV
        for i in 0..10 {
            t_recs.push(mk_read(
                format!("t_ref360_{i}"),
                340,
                &[(40, Kind::Match)],
                slice_1b(340, 40),
            ));
        }
        for i in 0..6 {
            let mut seq = slice_1b(340, 40);
            seq[20] = alt_of(360);
            n_recs.push(mk_read(
                format!("n_snv360_{i}"),
                340,
                &[(40, Kind::Match)],
                seq,
            ));
        }
        for i in 0..4 {
            n_recs.push(mk_read(
                format!("n_ref360_{i}"),
                340,
                &[(40, Kind::Match)],
                slice_1b(340, 40),
            ));
        }

        // Case 2: r = 470 sub-threshold tumor allele (1/40 = 0.025 < freq 0.05) + normal SNV (10/20)
        {
            let mut seq = slice_1b(450, 40);
            seq[20] = alt_of(470);
            t_recs.push(mk_read(
                "t_sub470".to_string(),
                450,
                &[(40, Kind::Match)],
                seq,
            ));
        }
        for i in 0..39 {
            t_recs.push(mk_read(
                format!("t_ref470_{i}"),
                450,
                &[(40, Kind::Match)],
                slice_1b(450, 40),
            ));
        }
        for i in 0..10 {
            let mut seq = slice_1b(450, 40);
            seq[20] = alt_of(470);
            n_recs.push(mk_read(
                format!("n_snv470_{i}"),
                450,
                &[(40, Kind::Match)],
                seq,
            ));
        }
        for i in 0..10 {
            n_recs.push(mk_read(
                format!("n_ref470_{i}"),
                450,
                &[(40, Kind::Match)],
                slice_1b(450, 40),
            ));
        }

        // Case 3: d = 600 12-bp deletion rescued via combine_analysis with combined AF = 5/120 <= 0.05
        for i in 0..3 {
            let mut seq = slice_1b(575, 25);
            seq.extend_from_slice(&slice_1b(612, 25));
            t_recs.push(mk_read(
                format!("t_del600_{i}"),
                575,
                &[(25, Kind::Match), (12, Kind::Deletion), (25, Kind::Match)],
                seq,
            ));
        }
        for i in 0..17 {
            t_recs.push(mk_read(
                format!("t_ref600_{i}"),
                575,
                &[(50, Kind::Match)],
                slice_1b(575, 50),
            ));
        }
        for i in 0..2 {
            let mut seq = slice_1b(575, 25);
            seq.extend_from_slice(&slice_1b(612, 6));
            n_recs.push(mk_read(
                format!("n_sc600_{i}"),
                575,
                &[(25, Kind::Match), (6, Kind::SoftClip)],
                seq,
            ));
        }
        for i in 0..98 {
            n_recs.push(mk_read(
                format!("n_ref600_{i}"),
                575,
                &[(50, Kind::Match)],
                slice_1b(575, 50),
            ));
        }

        t_recs.sort_by_key(|r| r.alignment_start());
        n_recs.sort_by_key(|r| r.alignment_start());
        let t_bam = dir.path().join("refmode_tumor.bam");
        let n_bam = dir.path().join("refmode_normal.bam");
        indexed_bam(&t_bam, "chr1", 1200, &t_recs).expect("write refmode_tumor.bam");
        indexed_bam(&n_bam, "chr1", 1200, &n_recs).expect("write refmode_normal.bam");

        let refmode_cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!("{}|{}", t_bam.display(), n_bam.display()),
            sample: "TUMOR|NORMAL".to_string(),
            region: Some("chr1:200-700".to_string()),
            freq: 0.05,
            af_threshold: 0.05,
            min_reads: 2,
            threads: 1,
            ..Config::default()
        };
        let mut keepall_cfg = refmode_cfg.clone();
        keepall_cfg.debug = true;

        let refmode_vcf = run_core_vcf(&refmode_cfg).expect("somatic RefSiteMode vcf");
        let keepall_vcf = run_core_vcf(&keepall_cfg).expect("somatic KeepAll debug vcf");
        assert_eq!(
            refmode_vcf, keepall_vcf,
            "somatic RefSiteMode VCF must match KeepAll (debug=true) byte-for-byte"
        );

        let refmode_var = run_core_var(&refmode_cfg).expect("somatic RefSiteMode var");
        let keepall_var_raw = run_core_var(&keepall_cfg).expect("somatic KeepAll debug var");
        let refmode_rows: Vec<String> = refmode_var
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(str::to_string)
            .collect();
        let keepall_rows_55: Vec<String> = keepall_var_raw
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let cols: Vec<&str> = l.split('\t').collect();
                assert_eq!(cols.len(), 56, "debug somatic .var row must emit 56 cols");
                cols[..55].join("\t")
            })
            .collect();
        assert_eq!(
            refmode_rows, keepall_rows_55,
            "somatic RefSiteMode .var 55-col table must match KeepAll (debug=true) identically"
        );
        assert!(
            refmode_rows.iter().any(|r| r.contains("\t250\t250\t")),
            "expected normal variant at p=250 right after p-1 insertion in somatic .var, got:\n{refmode_var}"
        );
        assert!(
            refmode_rows.iter().any(|r| r.contains("\t360\t360\t")),
            "expected hydrated CompactRefSite variant at q=360 in somatic .var, got:\n{refmode_var}"
        );
        assert!(
            refmode_rows
                .iter()
                .any(|r| r.contains("\t470\t470\t") && r.contains("\t40\t1\t")),
            "expected sub-threshold BAM1 allele (DP=40, VD=1) paired at r=470, got:\n{refmode_var}"
        );
        assert!(
            refmode_rows
                .iter()
                .any(|r| r.contains("\t599\t") && r.contains("\tGermline\t")),
            "expected combine_analysis low-AF deletion at 599..611 rescued as Germline, got:\n{refmode_var}"
        );
    }

    #[test]
    fn tiled_vs_untiled_sv_and_splice_conservation() {
        ensure_fastas_converted();
        let droot = data_root();

        // 1. Symbolic `<DEL>` and `<INV>` SV BAM fixtures: `tile_size = Some(200)` vs `None`
        let sv_cases = [
            (
                "VarDictJava/testdata/fastas/hg19.fa",
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr10.DELDiscReverse.bam",
                "chr10:89712340-89713976:testbed",
                "\t<DEL>\t",
            ),
            (
                "VarDictJava/testdata/fastas/hg19.fa",
                "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.INVSV.bam",
                "chr20:46842754-46845755:testbed",
                "\t<INV>\t",
            ),
        ];

        for (fa_rel, bam_rel, region, expected_sv_tag) in sv_cases {
            let bam = droot.join(bam_rel).to_string_lossy().into_owned();
            let untiled_cfg = Config {
                reference: droot.join(fa_rel).to_string_lossy().into_owned(),
                sample: extract_sample_name(&bam, None),
                bam,
                threads: 1,
                zero_based: true,
                zero_based_defined: true,
                region: Some(region.to_string()),
                freq: 0.001,
                af_threshold: 0.001,
                tile_size: None,
                ..Config::default()
            };
            let mut tiled_cfg = untiled_cfg.clone();
            tiled_cfg.tile_size = Some(200);
            tiled_cfg.threads = 4;

            let untiled_vcf = run_core_vcf(&untiled_cfg).expect("untiled SV VCF");
            let tiled_vcf = run_core_vcf(&tiled_cfg).expect("tiled (200bp) SV VCF");
            assert!(
                untiled_vcf.contains(expected_sv_tag),
                "expected symbolic SV {expected_sv_tag:?} in untiled VCF for {bam_rel}"
            );
            assert_eq!(
                tiled_vcf, untiled_vcf,
                "tiled (tile_size=Some(200)) SV VCF must match untiled None byte-for-byte on {bam_rel}"
            );

            let mut untiled_var =
                normalize_var_rows(&run_core_var(&untiled_cfg).expect("untiled SV .var"), 36);
            let mut tiled_var = normalize_var_rows(
                &run_core_var(&tiled_cfg).expect("tiled (200bp) SV .var"),
                36,
            );
            untiled_var.sort();
            tiled_var.sort();
            assert_eq!(
                tiled_var, untiled_var,
                "tiled (tile_size=Some(200)) SV .var row set must match untiled None on {bam_rel}"
            );
        }

        // 2. `-i` splice junction output: `tile_size = Some(200)` vs `None` byte-for-byte
        let splice_bam = droot
            .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Tcov_zero_SV_DEL_candidate3.bam")
            .to_string_lossy()
            .into_owned();
        let splice_untiled_cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg38.fa")
                .to_string_lossy()
                .into_owned(),
            sample: extract_sample_name(&splice_bam, None),
            bam: splice_bam,
            threads: 1,
            zero_based: true,
            zero_based_defined: true,
            region: Some("chr8:90534611-90620070:testbed".to_string()),
            freq: 0.1,
            af_threshold: 0.1,
            splice: true,
            tile_size: None,
            ..Config::default()
        };
        let mut splice_tiled_cfg = splice_untiled_cfg.clone();
        splice_tiled_cfg.tile_size = Some(200);
        splice_tiled_cfg.threads = 4;

        let untiled_splice = run_core_var(&splice_untiled_cfg).expect("untiled -i splice");
        let tiled_splice = run_core_var(&splice_tiled_cfg).expect("tiled (200bp) -i splice");
        let splice_rows: Vec<&str> = untiled_splice
            .lines()
            .filter(|l| !l.trim().is_empty())
            .collect();
        assert_eq!(
            splice_rows.len(),
            16,
            "expected 16 splice junction rows on Tcov_zero_SV_DEL_candidate3.bam"
        );
        assert_eq!(
            tiled_splice, untiled_splice,
            "tiled (tile_size=Some(200)) -i splice output must match untiled None byte-for-byte"
        );
    }
}

// 8. mod parity_edge_cases (5 tests)
mod parity_edge_cases {
    use super::*;
    use noodles::sam::alignment::record::Flags;
    use vardictrust_core::variations::Variant;

    fn lcg_ref_bytes(len: usize) -> Vec<u8> {
        let mut s = 0x9E37_79B9_7F4A_7C15_u64;
        let bases = [b'A', b'C', b'G', b'T'];
        (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                bases[((s >> 33) & 3) as usize]
            })
            .collect()
    }

    fn write_custom_ref(dir: &TempDir, tag: &str, name: &str, seq: &[u8]) -> std::path::PathBuf {
        let fa = dir.path().join(format!("{tag}.fa"));
        let mut content = Vec::with_capacity(name.len() + seq.len() + 4);
        content.push(b'>');
        content.extend_from_slice(name.as_bytes());
        content.push(b'\n');
        content.extend_from_slice(seq);
        content.push(b'\n');
        fs::write(&fa, &content).expect("write custom fasta");
        let len = seq.len();
        fs::write(
            fa.with_extension("fa.fai"),
            format!("{name}\t{len}\t{}\t{len}\t{}\n", name.len() + 2, len + 1),
        )
        .expect("write custom fai");
        fa
    }

    fn alt_base(b: u8) -> u8 {
        match b {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        }
    }

    fn mk_rec(
        name: &str,
        flags: Flags,
        pos_1b: usize,
        ops: &[(usize, Kind)],
        seq: Vec<u8>,
        qual: Vec<u8>,
        mate: Option<(usize, i32)>,
    ) -> RecordBuf {
        let cigar: Cigar = ops.iter().map(|&(l, k)| Op::new(k, l)).collect();
        let (mate_pos_1b, tlen) = mate.map_or((None, 0), |(mp, tl)| (Some(mp), tl));
        let mut b = RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(pos_1b).expect("pos"))
            .set_mapping_quality(noodles::sam::alignment::record::MappingQuality::new(60).unwrap())
            .set_cigar(cigar)
            .set_sequence(seq.into())
            .set_quality_scores(qual.into())
            .set_template_length(tlen);
        if let Some(mp) = mate_pos_1b {
            b = b
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(Position::try_from(mp).expect("mate pos"));
        }
        b.build()
    }

    #[test]
    fn gap1_split_tile_multi_bam_active_limit_and_tail_validation() {
        let dir = tempfile::tempdir().expect("tempdir");
        let ref_len = 605_000usize;
        let ref_seq = lcg_ref_bytes(ref_len);
        let ref_fa = write_custom_ref(&dir, "gap1_ref", "chr1", &ref_seq);

        let slice_1b = |pos_1b: usize, len: usize| ref_seq[pos_1b - 1..pos_1b - 1 + len].to_vec();
        let snv_at = |read_pos_1b: usize, len: usize, var_pos_1b: usize| {
            let mut s = slice_1b(read_pos_1b, len);
            let off = var_pos_1b - read_pos_1b;
            s[off] = alt_base(s[off]);
            s
        };

        // --- (A) Simple mode `b1:b2:b3` over 600 kbp (6 tiles of 100 kbp, 3 prescan bands):
        // b2 has 0 reads in Tile 0 (1..=100_000) and passing reads in Tile 3 (320_000);
        // b3 has SNV reads in Tile 0 (pos = 200, var at 225).
        let b1 = dir.path().join("g1_b1.bam");
        let b2 = dir.path().join("g1_b2.bam");
        let b3 = dir.path().join("g1_b3.bam");
        let b2_empty = dir.path().join("g1_b2_empty.bam");
        let b4_unindexed = dir.path().join("g1_b4_no_bai.bam");

        indexed_bam(
            &b1,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "b1_t0_1",
                    Flags::empty(),
                    200,
                    &[(50, Kind::Match)],
                    slice_1b(200, 50),
                    vec![40; 50],
                    None,
                ),
                mk_rec(
                    "b1_t3_1",
                    Flags::empty(),
                    320_000,
                    &[(50, Kind::Match)],
                    slice_1b(320_000, 50),
                    vec![40; 50],
                    None,
                ),
            ],
        )
        .expect("b1");

        indexed_bam(
            &b2,
            "chr1",
            ref_len,
            &[mk_rec(
                "b2_t3_only",
                Flags::empty(),
                320_100,
                &[(50, Kind::Match)],
                slice_1b(320_100, 50),
                vec![40; 50],
                None,
            )],
        )
        .expect("b2");

        indexed_bam(
            &b3,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "b3_t0_v1",
                    Flags::empty(),
                    200,
                    &[(50, Kind::Match)],
                    snv_at(200, 50, 225),
                    vec![38; 50],
                    None,
                ),
                mk_rec(
                    "b3_t0_v2",
                    Flags::REVERSE_COMPLEMENTED,
                    202,
                    &[(50, Kind::Match)],
                    snv_at(202, 50, 225),
                    vec![41; 50],
                    None,
                ),
                mk_rec(
                    "b3_t0_v3",
                    Flags::empty(),
                    204,
                    &[(50, Kind::Match)],
                    snv_at(204, 50, 225),
                    vec![39; 50],
                    None,
                ),
            ],
        )
        .expect("b3");

        indexed_bam(&b2_empty, "chr1", ref_len, &[]).expect("b2_empty");
        indexed_bam(&b4_unindexed, "chr1", ref_len, &[]).expect("b4_unindexed");
        fs::remove_file(b4_unindexed.with_extension("bam.bai")).expect("remove b4 bai");

        let simple_untiled = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!("{}:{}:{}", b1.display(), b2.display(), b3.display()),
            sample: "G1".to_string(),
            region: Some("chr1:1-600000".to_string()),
            threads: 1,
            tile_size: None,
            ..Config::default()
        };
        let expected_simple_vcf = run_core_vcf(&simple_untiled).expect("untiled simple vcf");
        let expected_simple_var = run_core_var(&simple_untiled).expect("untiled simple var");
        assert!(
            expected_simple_vcf.contains("chr1\t225\t"),
            "untiled b1:b2:b3 must retain b3 SNV at pos 225, got:\n{expected_simple_vcf}"
        );

        for threads in [1usize, 4] {
            let mut tiled_cfg = simple_untiled.clone();
            tiled_cfg.tile_size = Some(100_000);
            tiled_cfg.threads = threads;
            let mut vcf_buf = Vec::new();
            vardictrust_core::run_with_writer(&tiled_cfg, &mut vcf_buf).expect("tiled simple vcf");
            assert_eq!(
                String::from_utf8(vcf_buf).unwrap(),
                expected_simple_vcf,
                "tiled simple b1:b2:b3 (threads={threads}) must match untiled VCF"
            );
            let mut var_cfg = tiled_cfg.clone();
            var_cfg.emit_intermediate = true;
            let mut var_buf = Vec::new();
            vardictrust_core::run_with_writer(&var_cfg, &mut var_buf).expect("tiled simple var");
            assert_eq!(
                String::from_utf8(var_buf).unwrap(),
                expected_simple_var,
                "tiled simple b1:b2:b3 (threads={threads}) must match untiled .var"
            );
        }

        // When b2 is empty across the ENTIRE parent region (`b1:b2_empty:b3`), active_bam_limit == 2
        // stops before b3 on both untiled and tiled runs.
        let mut stopped_untiled = simple_untiled.clone();
        stopped_untiled.bam = format!("{}:{}:{}", b1.display(), b2_empty.display(), b3.display());
        let stopped_vcf = run_core_vcf(&stopped_untiled).expect("stopped untiled");
        assert!(
            extract_variants(&stopped_vcf).is_empty(),
            "empty b2 across parent region must stop before b3"
        );
        for threads in [1usize, 4] {
            let mut stopped_tiled = stopped_untiled.clone();
            stopped_tiled.tile_size = Some(100_000);
            stopped_tiled.threads = threads;
            let mut out = Vec::new();
            vardictrust_core::run_with_writer(&stopped_tiled, &mut out).expect("stopped tiled");
            assert_eq!(String::from_utf8(out).unwrap(), stopped_vcf);
        }

        // --- (B) AR-9: Trailing `b4` (`idx = 3 >= limit = 2`) with missing `.bai` still fails fast (#187)
        let mut missing_bai_cfg = stopped_untiled.clone();
        missing_bai_cfg.bam = format!(
            "{}:{}:{}:{}",
            b1.display(),
            b2_empty.display(),
            b3.display(),
            b4_unindexed.display()
        );
        missing_bai_cfg.tile_size = Some(100_000);
        missing_bai_cfg.threads = 2;
        let mut err_buf = Vec::new();
        assert!(
            vardictrust_core::run_with_writer(&missing_bai_cfg, &mut err_buf).is_err(),
            "trailing unindexed BAM beyond active_bam_limit must still fail fast (#187)"
        );

        // AR-9: Out-of-region SV & large-indel `partial_pipeline` still inspect full `bam_paths`
        // (`idx >= active_bam_limit`): region `chr1:1-2100` (`tile_size = Some(1000)`),
        // `sv_b1` has forward discordant+5'-matched split reads at pos 500 -> mate window ~10001
        // (outside loaded reference!) AND a 5'-clipped 400 bp `realignlgdel` anchor at pos 2020 -> 2470,
        // `sv_b2` has 0 reads in 1..=2100 (`active_bam_limit = 2`) and ref reads at 2480 & 10030,
        // `sv_b3` (`idx = 2 >= 2`) has 0 reads in 1..=2100 and reverse mate/3'-clip reads at 2470 & 10001!
        let sv_b1 = dir.path().join("g1_sv_b1.bam");
        let sv_b2 = dir.path().join("g1_sv_b2.bam");
        let sv_b3 = dir.path().join("g1_sv_b3.bam");
        let mut sv_b1_recs = Vec::new();
        for i in 0..4usize {
            let mut seq = slice_1b(500 + i, 50 - i);
            seq.extend_from_slice(&slice_1b(10_001, 20));
            sv_b1_recs.push(mk_rec(
                &format!("sv_f_{i}"),
                Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED,
                500 + i,
                &[(50 - i, Kind::Match), (20, Kind::SoftClip)],
                seq,
                vec![40; 70 - i],
                Some((10_001, 9551 - i as i32)),
            ));
        }
        for i in 0..3usize {
            let mut seq = slice_1b(2020 + i, 50 - i);
            seq.extend_from_slice(&slice_1b(2470, 20));
            sv_b1_recs.push(mk_rec(
                &format!("lgdel_5p_{i}"),
                Flags::empty(),
                2020 + i,
                &[(50 - i, Kind::Match), (20, Kind::SoftClip)],
                seq,
                vec![40; 70 - i],
                None,
            ));
        }
        sv_b1_recs.sort_by_key(|r| r.alignment_start());
        indexed_bam(&sv_b1, "chr1", ref_len, &sv_b1_recs).expect("sv_b1");
        indexed_bam(
            &sv_b2,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "sv_b2_lgdel_ref",
                    Flags::empty(),
                    2480,
                    &[(50, Kind::Match)],
                    slice_1b(2480, 50),
                    vec![40; 50],
                    None,
                ),
                mk_rec(
                    "sv_b2_mate_ref",
                    Flags::empty(),
                    10_030,
                    &[(50, Kind::Match)],
                    slice_1b(10_030, 50),
                    vec![40; 50],
                    None,
                ),
            ],
        )
        .expect("sv_b2");
        let mut sv_b3_recs = Vec::new();
        for i in 0..3usize {
            let mut seq = slice_1b(2050, 20);
            seq.extend_from_slice(&slice_1b(2470, 50));
            sv_b3_recs.push(mk_rec(
                &format!("lgdel_3p_{i}"),
                Flags::REVERSE_COMPLEMENTED,
                2470,
                &[(20, Kind::SoftClip), (50, Kind::Match)],
                seq,
                vec![40; 70],
                None,
            ));
        }
        for i in 0..4usize {
            let mut seq = slice_1b(530, 20);
            seq.extend_from_slice(&slice_1b(10_001, 50));
            sv_b3_recs.push(mk_rec(
                &format!("sv_r_{i}"),
                Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::REVERSE_COMPLEMENTED,
                10_001,
                &[(20, Kind::SoftClip), (50, Kind::Match)],
                seq,
                vec![40; 70],
                Some((500 + i, -(9551 - i as i32))),
            ));
        }
        sv_b3_recs.sort_by_key(|r| r.alignment_start());
        indexed_bam(&sv_b3, "chr1", ref_len, &sv_b3_recs).expect("sv_b3");

        let sv_untiled_cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!(
                "{}:{}:{}",
                sv_b1.display(),
                sv_b2.display(),
                sv_b3.display()
            ),
            sample: "G1SV".to_string(),
            region: Some("chr1:1-2100".to_string()),
            freq: 0.01,
            af_threshold: 0.01,
            threads: 1,
            tile_size: None,
            ..Config::default()
        };
        let sv_untiled_vcf = run_core_vcf(&sv_untiled_cfg).expect("sv untiled vcf");
        let mut sv_without_b3_cfg = sv_untiled_cfg.clone();
        sv_without_b3_cfg.bam = format!("{}:{}", sv_b1.display(), sv_b2.display());
        let sv_without_b3_vcf = run_core_vcf(&sv_without_b3_cfg).expect("sv without b3 vcf");
        assert_ne!(
            sv_untiled_vcf, sv_without_b3_vcf,
            "sv_b3 out-of-region reads at 2470 & 10001 must alter SV/indel calls vs b1:b2 alone"
        );
        let mut sv_tiled_cfg = sv_untiled_cfg.clone();
        sv_tiled_cfg.tile_size = Some(1_000);
        sv_tiled_cfg.threads = 2;
        let mut sv_tiled_buf = Vec::new();
        vardictrust_core::run_with_writer(&sv_tiled_cfg, &mut sv_tiled_buf).expect("sv tiled vcf");
        assert_eq!(
            String::from_utf8(sv_tiled_buf).unwrap(),
            sv_untiled_vcf,
            "tiled SV & realign_indels partial_pipeline must still inspect full bam_paths beyond active_bam_limit"
        );

        // --- (C) Somatic mode (`bam1 = t1:t2:t3`, `bam2 = n1:n2:n3`) across 600 kbp:
        // Both t2 and n2 have 0 reads in Tile 0 and passing reads in Tile 3;
        // t3 has tumor SNV reads at pos 225 in Tile 0, and n3 has normal SNV reads at pos 225
        // (turning Status into Germline and increasing normal DP!).
        let n1 = dir.path().join("g1_n1.bam");
        let n2 = dir.path().join("g1_n2.bam");
        let n3 = dir.path().join("g1_n3.bam");
        fs::copy(&b1, &n1).unwrap();
        fs::copy(b1.with_extension("bam.bai"), n1.with_extension("bam.bai")).unwrap();
        fs::copy(&b2, &n2).unwrap();
        fs::copy(b2.with_extension("bam.bai"), n2.with_extension("bam.bai")).unwrap();
        fs::copy(&b3, &n3).unwrap();
        fs::copy(b3.with_extension("bam.bai"), n3.with_extension("bam.bai")).unwrap();

        let som_untiled = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!(
                "{}:{}:{}|{}:{}:{}",
                b1.display(),
                b2.display(),
                b3.display(),
                n1.display(),
                n2.display(),
                n3.display()
            ),
            sample: "TUM|NRM".to_string(),
            region: Some("chr1:1-600000".to_string()),
            threads: 1,
            tile_size: None,
            ..Config::default()
        };
        let expected_som_vcf = run_core_vcf(&som_untiled).expect("untiled somatic vcf");
        let expected_som_var = run_core_var(&som_untiled).expect("untiled somatic var");
        assert!(
            expected_som_vcf.contains("chr1\t225\t") && expected_som_vcf.contains("STATUS=Germline;"),
            "untiled somatic t1:t2:t3|n1:n2:n3 must retain both t3 and n3 at pos 225 as Germline, got:\n{expected_som_vcf}"
        );
        for threads in [1usize, 4] {
            let mut som_tiled = som_untiled.clone();
            som_tiled.tile_size = Some(100_000);
            som_tiled.threads = threads;
            let mut vcf_buf = Vec::new();
            vardictrust_core::run_with_writer(&som_tiled, &mut vcf_buf).expect("tiled somatic vcf");
            assert_eq!(
                String::from_utf8(vcf_buf).unwrap(),
                expected_som_vcf,
                "tiled somatic t1:t2:t3|n1:n2:n3 (threads={threads}) must match untiled VCF"
            );
            let mut var_cfg = som_tiled.clone();
            var_cfg.emit_intermediate = true;
            let mut var_buf = Vec::new();
            vardictrust_core::run_with_writer(&var_cfg, &mut var_buf).expect("tiled somatic var");
            assert_eq!(
                String::from_utf8(var_buf).unwrap(),
                expected_som_var,
                "tiled somatic t1:t2:t3|n1:n2:n3 (threads={threads}) must match untiled .var"
            );
        }
    }

    #[test]
    fn gap2_amplicon_bam_feed_uses_region_bounds_without_halo() {
        let dir = tempfile::tempdir().expect("tempdir");
        let ref_len = 2_000usize;
        let ref_seq = lcg_ref_bytes(ref_len);
        let ref_fa = write_custom_ref(&dir, "gap2_ref", "chr1", &ref_seq);
        let slice_1b = |pos_1b: usize, len: usize| ref_seq[pos_1b - 1..pos_1b - 1 + len].to_vec();
        let snv_at = |read_pos_1b: usize, len: usize, var_pos_1b: usize| {
            let mut s = slice_1b(read_pos_1b, len);
            let off = var_pos_1b - read_pos_1b;
            s[off] = alt_base(s[off]);
            s
        };

        // 8-column 0-based amplicon BED -> 1-based [600, 700], insert [600, 700]
        let amp_bed = dir.path().join("gap2_amp.bed");
        fs::write(&amp_bed, "chr1\t599\t700\tAMP1\t0\t+\t599\t700\n").expect("write amp bed");

        // (A) Multi-BAM `-a` (`b1:b2:b3`): `b2` has a read ONLY in `[region.start - 500, region.start - 1]`
        // (`200..=249`) and 0 reads overlapping `[600, 700]`, so `SamView.java:20-21` stops before `b3`.
        let b1 = dir.path().join("g2_b1.bam");
        let b2_halo = dir.path().join("g2_b2_halo.bam");
        let b3_snv = dir.path().join("g2_b3_snv.bam");

        indexed_bam(
            &b1,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "a_ref_1",
                    Flags::empty(),
                    600,
                    &[(101, Kind::Match)],
                    slice_1b(600, 101),
                    vec![40; 101],
                    None,
                ),
                mk_rec(
                    "a_ref_2",
                    Flags::REVERSE_COMPLEMENTED,
                    600,
                    &[(101, Kind::Match)],
                    slice_1b(600, 101),
                    vec![40; 101],
                    None,
                ),
            ],
        )
        .expect("g2_b1");

        indexed_bam(
            &b2_halo,
            "chr1",
            ref_len,
            &[mk_rec(
                "b2_upstream_halo_only",
                Flags::empty(),
                200,
                &[(50, Kind::Match)],
                slice_1b(200, 50),
                vec![40; 50],
                None,
            )],
        )
        .expect("g2_b2_halo");

        indexed_bam(
            &b3_snv,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "b3_amp_v1",
                    Flags::empty(),
                    600,
                    &[(101, Kind::Match)],
                    snv_at(600, 101, 630),
                    vec![39; 101],
                    None,
                ),
                mk_rec(
                    "b3_amp_v2",
                    Flags::REVERSE_COMPLEMENTED,
                    601,
                    &[(100, Kind::Match)],
                    snv_at(601, 100, 630),
                    vec![41; 100],
                    None,
                ),
                mk_rec(
                    "b3_amp_v3",
                    Flags::empty(),
                    602,
                    &[(99, Kind::Match)],
                    snv_at(602, 99, 630),
                    vec![40; 99],
                    None,
                ),
            ],
        )
        .expect("g2_b3_snv");

        let amp_multi_cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: format!(
                "{}:{}:{}",
                b1.display(),
                b2_halo.display(),
                b3_snv.display()
            ),
            sample: "AMP_MULTI".to_string(),
            bed: Some(amp_bed.to_string_lossy().into_owned()),
            amplicon: Some("10:0.95".to_string()),
            zero_based: true,
            zero_based_defined: true,
            freq: 0.01,
            af_threshold: 0.01,
            threads: 1,
            ..Config::default()
        };
        let multi_vcf = run_core_vcf(&amp_multi_cfg).expect("amplicon multi-bam vcf");
        assert!(
            extract_variants(&multi_vcf).is_empty(),
            "amplicon b2 with reads only in upstream 500bp halo must stop before b3 (SamView.java:20-21), got:\n{multi_vcf}"
        );

        // (B) `-a -t` (`dedup`): an SE upstream-halo read at `pos = 595` (`3M`, `595..=597` < `600`)
        // must NOT enter `SamView.queryOverlapping(chr, 600, 700)` and must NOT deduplicate a valid
        // amplicon read starting at `pos = 595` (`106M`, `595..=700`, overlap 101/106 = 0.9528 >= 0.95)!
        let b_dedup = dir.path().join("g2_dedup.bam");
        indexed_bam(
            &b_dedup,
            "chr1",
            ref_len,
            &[
                mk_rec(
                    "halo_short_non_ov",
                    Flags::empty(),
                    595,
                    &[(3, Kind::Match)],
                    slice_1b(595, 3),
                    vec![40; 3],
                    None,
                ),
                mk_rec(
                    "amp_snv_pos595",
                    Flags::empty(),
                    595,
                    &[(106, Kind::Match)],
                    snv_at(595, 106, 630),
                    vec![39; 106],
                    None,
                ),
                mk_rec(
                    "amp_snv_pos596",
                    Flags::REVERSE_COMPLEMENTED,
                    596,
                    &[(105, Kind::Match)],
                    snv_at(596, 105, 630),
                    vec![41; 105],
                    None,
                ),
                mk_rec(
                    "amp_ref_pos600",
                    Flags::empty(),
                    600,
                    &[(101, Kind::Match)],
                    slice_1b(600, 101),
                    vec![40; 101],
                    None,
                ),
            ],
        )
        .expect("g2_dedup");

        let mut amp_dedup_cfg = amp_multi_cfg.clone();
        amp_dedup_cfg.bam = b_dedup.to_string_lossy().into_owned();
        amp_dedup_cfg.dedup = true;
        amp_dedup_cfg.min_reads = 2;
        let dedup_vcf = run_core_vcf(&amp_dedup_cfg).expect("amplicon -t vcf");
        let dedup_vars = extract_variants(&dedup_vcf);
        assert_eq!(
            dedup_vars.len(),
            1,
            "upstream non-overlapping read 595..=597 must not -t dedup amplicon read at pos=595, got:\n{dedup_vcf}"
        );
        assert!(dedup_vars[0].starts_with("chr1\t630\t"));
    }

    #[test]
    fn gap3_crispr_filter_counts_eqx_cigar_ops() {
        let dir = tempfile::tempdir().expect("tempdir");
        let ref_fa = write_ref(&dir, "chr1", 400, true);
        let mut seq_125 = vec![b'A'; 150];
        seq_125[75] = b'C'; // 1-based pos 200
        let mut seq_126 = vec![b'A'; 150];
        seq_126[74] = b'C'; // 1-based pos 200

        let eqx_bam = dir.path().join("gap3_eqx.bam");
        indexed_bam(
            &eqx_bam,
            "chr1",
            400,
            &[
                // Boundary read at pos = 75: 75 + 150 - 200 = 25 (NOT > 25 -> must be skipped!)
                mk_rec(
                    "eqx_skip_75",
                    Flags::empty(),
                    75,
                    &[
                        (125, Kind::SequenceMatch),
                        (1, Kind::SequenceMismatch),
                        (24, Kind::SequenceMatch),
                    ],
                    {
                        let mut s = vec![b'A'; 150];
                        s[125] = b'G';
                        s
                    },
                    vec![40; 150],
                    None,
                ),
                // Target read 1 at pos = 125: 75=1X74= -> rlen3 = 150, 200-125=75>25 && 125+150-200=75>25
                mk_rec(
                    "eqx_keep_125",
                    Flags::empty(),
                    125,
                    &[
                        (75, Kind::SequenceMatch),
                        (1, Kind::SequenceMismatch),
                        (74, Kind::SequenceMatch),
                    ],
                    seq_125,
                    vec![38; 150],
                    None,
                ),
                // Target read 2 at pos = 126: 74=1X75= -> rlen3 = 150, 200-126=74>25 && 126+150-200=76>25
                mk_rec(
                    "eqx_keep_126",
                    Flags::REVERSE_COMPLEMENTED,
                    126,
                    &[
                        (74, Kind::SequenceMatch),
                        (1, Kind::SequenceMismatch),
                        (75, Kind::SequenceMatch),
                    ],
                    seq_126,
                    vec![41; 150],
                    None,
                ),
            ],
        )
        .expect("gap3_eqx.bam");

        for realign in [0, 1] {
            let cfg = Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                bam: eqx_bam.to_string_lossy().into_owned(),
                sample: "CRISPR_EQX".to_string(),
                region: Some("chr1:100-300".to_string()),
                crispr_cutting_site: Some(200),
                crispr_flanking_size: Some(25),
                realign,
                read_pos: -100.0,
                min_reads: 2,
                freq: 0.01,
                af_threshold: 0.01,
                threads: 1,
                ..Config::default()
            };
            let var_out = run_core_var(&cfg).expect("crispr =/X .var");
            let rows = normalize_var_rows(&var_out, 37); // 36 + 1 CRISPR column
            assert_eq!(
                rows.len(),
                1,
                "CRISPR -J 200 -j 25 -k {realign} must retain 75=1X74= variant at pos 200 and skip pos 75, got:\n{var_out}"
            );
            assert!(
                rows[0].contains("\tchr1\t200\t200\tA\tC\t"),
                "expected SNV A->C at pos 200 under -k {realign}, got: {}",
                rows[0]
            );
        }
    }

    #[test]
    fn gap4_tiled_splice_prescan_halo_max_rl_and_boundary_dedup() {
        let dir = tempfile::tempdir().expect("tempdir");
        let ref_len = 265_000usize;
        let ref_seq = lcg_ref_bytes(ref_len);
        let ref_fa = write_custom_ref(&dir, "gap4_ref", "chr1", &ref_seq);
        let slice_1b = |pos_1b: usize, len: usize| ref_seq[pos_1b - 1..pos_1b - 1 + len].to_vec();

        let rc12 = |pos_1b: usize| -> Vec<u8> {
            let mut v = slice_1b(pos_1b, 12);
            vardictrust_core::util::reverse_complement(&mut v);
            v
        };

        // Region `chr1:1001-261000`, `tile_size = Some(1000)` (260 tiles, 2 `250 kbp` prescan bands!):
        // Tile 0: start=1001..=2000 (halo 1001..=2500)
        // Tile 1: start=2001..=3000 (halo 1501..=3500)
        let mut recs = Vec::new();

        // 1. Seed initial `max_read_length = 150` at pos = 1100 in Tile 0.
        recs.push(mk_rec(
            "t0_seed_150m",
            Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED,
            1100,
            &[(150, Kind::Match)],
            slice_1b(1100, 150),
            vec![40; 150],
            Some((20, 150)),
        ));

        // 2. Boundary `-i -t` (`dedup`) test at `pos = 1400 < tile1.halo_start (1501)`:
        //    - `early_non_ov`: SE read at pos = 1400, 50M (`1400..=1449`, ends < 1501!)
        //    - `spliced_dup`: SE read at pos = 1400, 40M600N40M200N40M (M=120 <= 150; introns `1439-2039` in Tile 0
        //      and `2079-2279` owned by Tile 1!)
        //    - `spliced_unique`: SE read at pos = 1405, 35M600N40M200N40M (M=115 <= 150; introns `1439-2039` and `2079-2279`, count = 1!)
        recs.push(mk_rec(
            "early_non_ov_1400",
            Flags::empty(),
            1400,
            &[(50, Kind::Match)],
            slice_1b(1400, 50),
            vec![40; 50],
            None,
        ));
        {
            let mut seq = slice_1b(1400, 40);
            seq.extend_from_slice(&slice_1b(2040, 40));
            seq.extend_from_slice(&slice_1b(2280, 40));
            recs.push(mk_rec(
                "spliced_dup_1400",
                Flags::empty(),
                1400,
                &[
                    (40, Kind::Match),
                    (600, Kind::Skip),
                    (40, Kind::Match),
                    (200, Kind::Skip),
                    (40, Kind::Match),
                ],
                seq,
                vec![40; 120],
                None,
            ));
        }
        {
            let mut seq = slice_1b(1405, 35);
            seq.extend_from_slice(&slice_1b(2040, 40));
            seq.extend_from_slice(&slice_1b(2280, 40));
            recs.push(mk_rec(
                "spliced_unique_1405",
                Flags::empty(),
                1405,
                &[
                    (35, Kind::Match),
                    (600, Kind::Skip),
                    (40, Kind::Match),
                    (200, Kind::Skip),
                    (40, Kind::Match),
                ],
                seq,
                vec![40; 115],
                None,
            ));
        }

        // 3. AR-2 halo `MaxRlStep` alignment test at `pos = 1600` inside `[tile1.halo_start (1501), tile1.start (2001))`:
        //    All three reads are paired with distinct mapped `PNEXT >= 10` so `-t` keeps all three:
        //    - Read 0 (`pos = 1600`, `50M`, `read_end = 1650`): overlaps halo `1501..=3500`, NOT core `2001..=3000`.
        //    - Read 1 (`pos = 1600`, `12S40M400N40M300N40M12S`, `[MIS] = 144 <= 150`, `read_end = 2420`,
        //      introns `1639-2039` and `2079-2379` owned by Tile 1):
        //      5' `12S` is RC of unique 12-mer at `ref_pos = 2000` (`|1600 - 2000| = 400`).
        //      Before Read 2, `max_read_length == 150` (`2 * 150 = 300 <= 400`), so `is_chimeric_seed_hit` is FALSE
        //      and `12S...12S` drops Read 1! (If Read 2's `new_max = 300` step misfires early because `eff_start`
        //      used `2001` instead of `1501`, `2 * 300 = 600 > 400` strips 5' `12S` and emits phantom `2079-2379`!)
        //    - Read 2 (`pos = 1600`, `300M`, `read_end = 1900`): overlaps halo `1501..=3500`, NOT core `2001..=3000`,
        //      raising `max_read_length` to `300` AFTER Read 1!
        recs.push(mk_rec(
            "halo_r0_50m",
            Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED,
            1600,
            &[(50, Kind::Match)],
            slice_1b(1600, 50),
            vec![40; 50],
            Some((25, 50)),
        ));
        {
            let mut seq = rc12(2000);
            seq.extend_from_slice(&slice_1b(1600, 40));
            seq.extend_from_slice(&slice_1b(2040, 40));
            seq.extend_from_slice(&slice_1b(2380, 40));
            seq.extend_from_slice(b"AAAAAAAAAAAA");
            recs.push(mk_rec(
                "halo_r1_chimeric_12s_12s",
                Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED,
                1600,
                &[
                    (12, Kind::SoftClip),
                    (40, Kind::Match),
                    (400, Kind::Skip),
                    (40, Kind::Match),
                    (300, Kind::Skip),
                    (40, Kind::Match),
                    (12, Kind::SoftClip),
                ],
                seq,
                vec![40; 144],
                Some((35, 820)),
            ));
        }
        recs.push(mk_rec(
            "halo_r2_300m",
            Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED,
            1600,
            &[(300, Kind::Match)],
            slice_1b(1600, 300),
            vec![40; 300],
            Some((45, 300)),
        ));

        recs.sort_by_key(|r| r.alignment_start());
        let splice_bam = dir.path().join("gap4_splice.bam");
        indexed_bam(&splice_bam, "chr1", ref_len, &recs).expect("gap4_splice.bam");

        // Decoy normal BAM for `-b "t|n" -i` routing test (`VarDictLauncher.java:54-55` -> `bam1` only).
        let decoy_normal = dir.path().join("gap4_normal_decoy.bam");
        {
            let mut dseq = slice_1b(2500, 50);
            dseq.extend_from_slice(&slice_1b(2700, 50));
            indexed_bam(
                &decoy_normal,
                "chr1",
                ref_len,
                &[mk_rec(
                    "decoy_intron",
                    Flags::empty(),
                    2500,
                    &[(50, Kind::Match), (150, Kind::Skip), (50, Kind::Match)],
                    dseq,
                    vec![40; 100],
                    None,
                )],
            )
            .expect("decoy_normal");
        }

        let untiled_cfg = Config {
            reference: ref_fa.to_string_lossy().into_owned(),
            bam: splice_bam.to_string_lossy().into_owned(),
            sample: "SPLICE_S".to_string(),
            sample_explicit: true,
            region: Some("chr1:1001-261000".to_string()),
            splice: true,
            dedup: true,
            realign: 1,
            threads: 1,
            tile_size: None,
            ..Config::default()
        };
        let expected_splice = run_core_var(&untiled_cfg).expect("untiled -i -t");
        assert_eq!(
            expected_splice.trim(),
            "SPLICE_S\tchr1\t1439-2039\t1\nSPLICE_S\tchr1\t2079-2279\t1",
            "untiled -i -t must count 1439-2039 and 2079-2279 once (deduping pos=1400) and drop 1639-2039 / 2079-2379 (max_rl=150 before halo_r2_300m)"
        );

        for bam_spec in [
            splice_bam.to_string_lossy().into_owned(),
            format!("{}|{}", splice_bam.display(), decoy_normal.display()),
        ] {
            for threads in [1usize, 4] {
                let mut tiled_cfg = untiled_cfg.clone();
                tiled_cfg.bam = bam_spec.clone();
                tiled_cfg.tile_size = Some(1_000);
                tiled_cfg.threads = threads;
                let mut out = Vec::new();
                vardictrust_core::run_with_writer(&tiled_cfg, &mut out).expect("tiled -i -t");
                assert_eq!(
                    String::from_utf8(out).unwrap(),
                    expected_splice,
                    "tiled -i -t (threads={threads}, bam={bam_spec}) must match untiled byte-for-byte"
                );
            }
        }
    }

    #[test]
    fn gap5_vcf_chrom_spaces_unicode_linesep_and_debug_half_even() {
        // --- (A) Single-sample `-N "My Sample"` + hostile Unicode line separators (`\u{0085}`, `\u{2028}`, `\u{2029}`)
        let mut single_cfg = base_micro_config();
        single_cfg.sample = "My Sample".to_string();
        single_cfg.sample_explicit = true;
        let single_vcf = run_core_vcf(&single_cfg).expect("single vcf with space");
        let single_chrom = single_vcf
            .lines()
            .find(|l| l.starts_with("#CHROM\t"))
            .expect("#CHROM line");
        assert_eq!(
            single_chrom,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMy Sample"
        );
        let single_vars = extract_variants(&single_vcf);
        assert!(!single_vars.is_empty());
        assert!(
            single_vars[0].contains("\tSAMPLE=My_Sample;TYPE="),
            "INFO/SAMPLE must translate spaces to underscores (var2vcf_valid.pl:44,236), got: {}",
            single_vars[0]
        );

        let mut unicode_single_cfg = single_cfg.clone();
        unicode_single_cfg.sample = "My Sample\u{0085}A\u{2028}B\u{2029}C".to_string();
        let unicode_single_vcf = run_core_vcf(&unicode_single_cfg).expect("unicode single vcf");
        let unicode_chrom_lines: Vec<&str> = unicode_single_vcf
            .lines()
            .filter(|l| l.starts_with("#CHROM\t"))
            .collect();
        assert_eq!(
            unicode_chrom_lines,
            vec!["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMy Sample_A_B_C"],
            "Unicode NEL/LS/PS must be replaced with '_' while preserving ASCII space in #CHROM"
        );
        assert!(
            extract_variants(&unicode_single_vcf)[0].contains("\tSAMPLE=My_Sample_A_B_C;TYPE="),
            "INFO/SAMPLE must replace both ASCII space and Unicode whitespace with '_'"
        );

        // --- (B) Paired somatic `-N "Tumor A|Normal B"` and `-N "Tumor A"` -> `Tumor A-match`
        let dir = tempfile::tempdir().expect("tempdir");
        let ref_fa = write_ref(&dir, "chr1", 400, true);
        let t_bam = dir.path().join("t_space.bam");
        let n_bam = dir.path().join("n_space.bam");
        let mk_snv_125 = |name: &str, flags: Flags, start_1b: usize, q: u8| {
            let mut seq = vec![b'A'; 50];
            seq[125 - start_1b] = b'C';
            mk_rec(
                name,
                flags,
                start_1b,
                &[(50, Kind::Match)],
                seq,
                vec![q; 50],
                None,
            )
        };
        indexed_bam(
            &t_bam,
            "chr1",
            400,
            &[
                mk_snv_125("t1", Flags::empty(), 100, 38),
                mk_snv_125("t2", Flags::REVERSE_COMPLEMENTED, 102, 40),
                mk_snv_125("t3", Flags::empty(), 104, 41),
                mk_snv_125("t4", Flags::REVERSE_COMPLEMENTED, 106, 39),
            ],
        )
        .expect("t_bam");
        indexed_bam(
            &n_bam,
            "chr1",
            400,
            &[
                mk_rec(
                    "n1",
                    Flags::empty(),
                    100,
                    &[(50, Kind::Match)],
                    vec![b'A'; 50],
                    vec![40; 50],
                    None,
                ),
                mk_rec(
                    "n2",
                    Flags::REVERSE_COMPLEMENTED,
                    102,
                    &[(50, Kind::Match)],
                    vec![b'A'; 50],
                    vec![40; 50],
                    None,
                ),
            ],
        )
        .expect("n_bam");

        for (sample_arg, expected_tail, expected_info_sample) in [
            (
                "Tumor A|Normal B",
                "\tFORMAT\tTumor A\tNormal B",
                ";SAMPLE=Tumor_A;TYPE=",
            ),
            (
                "Tumor A",
                "\tFORMAT\tTumor A\tTumor A-match",
                ";SAMPLE=Tumor_A;TYPE=",
            ),
            (
                "Tumor A\u{2028}X|Normal B\u{2029}Y\u{0085}Z",
                "\tFORMAT\tTumor A_X\tNormal B_Y_Z",
                ";SAMPLE=Tumor_A_X;TYPE=",
            ),
        ] {
            let som_cfg = Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                bam: format!("{}|{}", t_bam.display(), n_bam.display()),
                sample: sample_arg.to_string(),
                sample_explicit: true,
                region: Some("chr1:80-200".to_string()),
                min_reads: 2,
                freq: 0.01,
                af_threshold: 0.01,
                threads: 1,
                ..Config::default()
            };
            let som_vcf = run_core_vcf(&som_cfg).expect("paired space vcf");
            let som_chrom = som_vcf
                .lines()
                .find(|l| l.starts_with("#CHROM\t"))
                .expect("paired #CHROM");
            assert!(
                som_chrom.ends_with(expected_tail),
                "expected paired #CHROM ending with {expected_tail:?} for -N {sample_arg:?}, got: {som_chrom:?}"
            );
            let som_vars = extract_variants(&som_vcf);
            assert!(!som_vars.is_empty());
            assert!(
                som_vars[0].contains(expected_info_sample),
                "expected {expected_info_sample:?} in paired INFO, got: {}",
                som_vars[0]
            );
        }

        // --- (C) Gap 5b: `Variant::debug_variants_content` (`Variant.java:268-282`) `round_half_even`
        // including `1.0 / 20_000.0` (`5e-5` `DecimalFormat("0.0000")` underflow tie `M-D1`) -> `"0.0000"`
        let mut v = Variant::new();
        v.description = "C".to_string();
        v.fwd = 1;
        v.rev = 0;
        v.frequency = 1.0 / 20_000.0;
        v.bias = "0;1".to_string();
        v.mean_position = 12.25; // exact binary half-way tie -> HALF_EVEN rounds to 12.2 (not 12.3!)
        v.pstd = true;
        v.mean_quality = 30.25; // exact binary half-way tie -> HALF_EVEN rounds to 30.2 (not 30.3!)
        v.qstd = true;
        v.hi_freq = 1.0 / 20_000.0;
        v.mean_mapping_quality = 60.0;
        v.hq_lq_ratio = 1.0 / 2_000.0; // 5e-4 underflow tie for "0.000" -> "0.000" (not "0.001"!)
        assert_eq!(
            v.debug_variants_content(),
            "C:1:F-1:R-0:0.0000:0;1:12.2:1:30.2:1:0.0000:60.0:0.000"
        );
    }
}
