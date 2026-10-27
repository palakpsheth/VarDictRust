//! First-class functional engine behavior integration suite (`25` tests).
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

// 5. mod filtering_and_pileup (12 tests)
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
            out_mapq: 0,
            threads: 1,
            ..Config::default()
        };
        let mut cfg_filtered = cfg_unfiltered.clone();
        cfg_filtered.out_mapq = 65;

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

// 7. mod concurrency (1 test)
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
}
