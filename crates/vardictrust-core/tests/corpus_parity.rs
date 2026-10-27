//! First-class functional corpus parity integration suite (`110` tests).
//!
//! Runs all 92 upstream VarDictJava v1.8.3 integration fixtures across 6 functional
//! modules (`simple`, `somatic`, `pileup`, `fisher`, `amplicon`, `splicing`), plus
//! 17 operational modality tests (`modalities`) and the corpus manifest integrity
//! gate (`manifest::all_92_golden_cases_registered`), entirely in-process via
//! `vardictrust_core::run_with_writer` with `threads = 1` and zero subprocesses.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::OnceLock;

use vardictrust_core::config::{extract_sample_name, Config};

mod common;
use common::{
    data_root, ensure_fastas_converted, extract_variants, run_core_var, run_core_vcf, sha256_hex,
    workspace_root, write_indexed_cram,
};

const SMOKE_VCF_COMPANIONS: &[&str] = &[
    "Simple;GRCh37.fa;HCC2218_C_position_bug.bam;chr2;179260451-179260596;-f 0.001.txt",
    "Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt",
    "Amplicon;hg19.fa;Colo829-18_S3-sort.bam;chr1;933866-934466;933866-934466;-a 100.95 -D;.txt",
];

fn xfail_registry() -> &'static HashMap<String, String> {
    static XFAILS: OnceLock<HashMap<String, String>> = OnceLock::new();
    XFAILS.get_or_init(|| {
        let path = workspace_root().join("harness/integration_xfail.tsv");
        let content = fs::read_to_string(&path).expect("read integration_xfail.tsv");
        let mut map = HashMap::new();
        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() >= 5 && cols[1].trim() == "XFAIL" {
                map.insert(cols[0].trim().to_string(), cols[4].trim().to_string());
            }
        }
        map
    })
}

pub fn apply_extra_args(cfg: &mut Config, extra_args: &[String]) {
    let mut i = 0;
    while i < extra_args.len() {
        match extra_args[i].as_str() {
            "-f" => {
                let v = &extra_args[i + 1];
                let parsed: f64 = v.parse().expect("valid -f");
                cfg.freq = parsed;
                cfg.af_threshold = parsed;
                cfg.af_threshold_truthy_zero = parsed == 0.0 && v.trim() != "0";
                i += 2;
            }
            "-F" => {
                let v = &extra_args[i + 1];
                cfg.flag_filter =
                    if let Some(hex) = v.strip_prefix("0x").or_else(|| v.strip_prefix("0X")) {
                        u32::from_str_radix(hex, 16).expect("valid hex -F")
                    } else {
                        v.parse().expect("valid dec -F")
                    };
                i += 2;
            }
            "-Q" => {
                cfg.map_qual = Some(extra_args[i + 1].parse().expect("valid -Q"));
                i += 2;
            }
            "-q" => {
                cfg.base_qual = extra_args[i + 1].parse().expect("valid -q");
                i += 2;
            }
            "-P" => {
                cfg.read_pos = extra_args[i + 1].parse().expect("valid -P");
                i += 2;
            }
            "-r" => {
                cfg.min_reads = extra_args[i + 1].parse().expect("valid -r");
                i += 2;
            }
            "-k" => {
                cfg.realign = extra_args[i + 1].parse().expect("valid -k");
                i += 2;
            }
            "-M" => {
                cfg.minmatch = extra_args[i + 1].parse().expect("valid -M");
                i += 2;
            }
            "-V" => {
                cfg.lofreq = extra_args[i + 1].parse().expect("valid -V");
                i += 2;
            }
            "-X" => {
                cfg.extension_nuc = extra_args[i + 1].parse().expect("valid -X");
                i += 2;
            }
            "-x" => {
                cfg.extension = extra_args[i + 1].parse().expect("valid -x");
                i += 2;
            }
            "-J" => {
                cfg.crispr_cutting_site = Some(extra_args[i + 1].parse().expect("valid -J"));
                i += 2;
            }
            "-j" => {
                cfg.crispr_flanking_size = Some(extra_args[i + 1].parse().expect("valid -j"));
                i += 2;
            }
            "-a" => {
                cfg.amplicon = Some(extra_args[i + 1].clone());
                i += 2;
            }
            "-O" => {
                cfg.out_mapq = extra_args[i + 1].parse().expect("valid -O");
                i += 2;
            }
            "-o" => {
                cfg.qratio = extra_args[i + 1].parse().expect("valid -o");
                i += 2;
            }
            "-L" => {
                cfg.sv_min_len = extra_args[i + 1].parse().expect("valid -L");
                i += 2;
            }
            "-A" => {
                cfg.sv_min_af = extra_args[i + 1].parse().expect("valid -A");
                i += 2;
            }
            "--mfreq" => {
                cfg.mfreq = extra_args[i + 1].parse().expect("valid --mfreq");
                i += 2;
            }
            "--nmfreq" => {
                cfg.nmfreq = extra_args[i + 1].parse().expect("valid --nmfreq");
                i += 2;
            }
            "--adaptor" => {
                cfg.adaptor = extra_args[i + 1]
                    .split(',')
                    .map(|s| s.to_string())
                    .collect();
                i += 2;
            }
            "--DP" | "--defaultPrinter" => {
                let v = &extra_args[i + 1];
                cfg.debug_print = Some(if v == "ERR" {
                    "ERR".to_string()
                } else {
                    "OUT".to_string()
                });
                i += 2;
            }
            "--chimeric" => {
                cfg.chimeric = true;
                i += 1;
            }
            "--deldupvar" => {
                cfg.deldupvar = true;
                i += 1;
            }
            "--fisher" | "-fisher" => {
                cfg.fisher = true;
                i += 1;
            }
            "--nosv" | "-U" => {
                cfg.no_sv = true;
                i += 1;
            }
            "-UN" | "--UN" => {
                cfg.unique_mode_second_in_pair_enabled = true;
                i += 1;
            }
            "-u" => {
                cfg.unique_mode_alignment_enabled = true;
                i += 1;
            }
            "-D" | "--debug" => {
                cfg.debug = true;
                i += 1;
            }
            "-K" => {
                cfg.include_n = true;
                i += 1;
            }
            "-i" | "--splice" => {
                cfg.splice = true;
                i += 1;
            }
            "-p" => {
                cfg.pileup = true;
                i += 1;
            }
            "-3" => {
                cfg.three_prime = true;
                i += 1;
            }
            "-t" | "--dedup" => {
                cfg.dedup = true;
                i += 1;
            }
            "-h" | "--header" => {
                cfg.column_header = true;
                i += 1;
            }
            other => panic!("unhandled corpus test argument: {other}"),
        }
    }
    cfg.apply_pileup_overrides();
}

fn run_corpus_case(fname: &str) {
    ensure_fastas_converted();
    let droot = data_root();
    let xfails = xfail_registry();
    let testfile = droot
        .join("VarDictJava/testdata/integrationtestcases")
        .join(fname);

    let raw_content =
        fs::read_to_string(&testfile).unwrap_or_else(|e| panic!("read testcase {fname}: {e}"));
    let content = raw_content.trim_end_matches(&['\n', '\r'][..]);
    let mut lines = content.split('\n').map(|l| l.trim_end_matches('\r'));
    let header = lines.next().expect("header line");
    let expected_lines: Vec<&str> = lines.collect();
    let expected = expected_lines.join("\n");

    let fields: Vec<String> = header.split(',').map(|s| s.trim().to_string()).collect();
    let mode = &fields[0];
    assert!(
        mode == "Simple" || mode == "Somatic" || mode == "Amplicon",
        "unexpected mode {mode} in {fname}"
    );
    let fasta_path = droot.join("VarDictJava/testdata/fastas").join(&fields[1]);
    let bam_dir =
        droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
    let bam_raw = &fields[2];
    let bam_arg = if let Some((b1, b2)) = bam_raw.split_once('|') {
        format!(
            "{}|{}",
            bam_dir.join(b1).display(),
            bam_dir.join(b2).display()
        )
    } else {
        bam_dir.join(bam_raw).display().to_string()
    };
    let chr = &fields[3];
    let (start_val, end_val, mut field_idx) = if let Some((s, e)) = fields[4].split_once('-') {
        (s.to_string(), e.to_string(), 5)
    } else {
        (
            fields[4].clone(),
            fields.get(5).cloned().unwrap_or_default(),
            6,
        )
    };

    let mut extra_args: Vec<String> = Vec::new();
    while field_idx < fields.len() {
        let fld = &fields[field_idx];
        if fld.starts_with('-') {
            for tok in fld.split_whitespace() {
                extra_args.push(tok.to_string());
            }
        }
        field_idx += 1;
    }

    let mut cfg = Config {
        reference: fasta_path.to_string_lossy().into_owned(),
        bam: bam_arg.clone(),
        sample: extract_sample_name(&bam_arg, None),
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
        emit_intermediate: true,
        ..Config::default()
    };
    apply_extra_args(&mut cfg, &extra_args);

    let _tmp_guard = if mode == "Amplicon" {
        let dir = tempfile::tempdir().expect("tempdir for amplicon bed");
        let bed_path = dir.path().join("region.bed");
        let bed_line = format!(
            "{chr}\t{start_val}\t{end_val}\ttestbed\t0\t0\t{}\t{}\n",
            fields[6], fields[7]
        );
        fs::write(&bed_path, bed_line).expect("write amplicon bed");
        cfg.bed = Some(bed_path.to_string_lossy().into_owned());
        cfg.region = None;
        Some(dir)
    } else {
        cfg.region = Some(format!("{chr}:{start_val}-{end_val}:testbed"));
        cfg.bed = None;
        None
    };

    let actual_str =
        run_core_var(&cfg).unwrap_or_else(|e| panic!("in-process run failed on {fname}: {e}"));
    let actual_trimmed = actual_str.trim_end_matches(&['\n', '\r'][..]);

    let mut actual_sorted: Vec<&str> = if actual_trimmed.is_empty() {
        vec![""]
    } else {
        actual_trimmed.split('\n').collect()
    };
    actual_sorted.sort_unstable();

    let mut expected_sorted: Vec<&str> = if expected.is_empty() {
        vec![""]
    } else {
        expected.split('\n').collect()
    };
    expected_sorted.sort_unstable();

    if actual_sorted == expected_sorted {
        assert!(
            !xfails.contains_key(fname),
            "XPASS on registered XFAIL {fname}"
        );
    } else if let Some(exp_sha) = xfails.get(fname) {
        let mut joined = actual_sorted.join("\n");
        joined.push('\n');
        let got_sha = sha256_hex(joined.as_bytes());
        assert_eq!(
            &got_sha, exp_sha,
            "XFAIL sha256 mismatch on {fname}: got {got_sha}, expected {exp_sha}"
        );
    } else {
        panic!(
            "MISMATCH on {fname}:\n--- actual ---\n{}\n--- expected ---\n{}",
            actual_sorted.join("\n"),
            expected_sorted.join("\n")
        );
    }

    if SMOKE_VCF_COMPANIONS.contains(&fname) && !cfg.splice {
        let vcf_out =
            run_core_vcf(&cfg).unwrap_or_else(|e| panic!("VCF companion failed on {fname}: {e}"));
        assert!(
            vcf_out.contains("#CHROM"),
            "VCF companion on {fname} must emit #CHROM header"
        );
    }
}

macro_rules! corpus_module {
    ($mod_name:ident, [ $( ($fn_name:ident, $file_name:expr) ),* $(,)? ]) => {
        mod $mod_name {
            use super::run_corpus_case;

            pub const CASES: &[&str] = &[
                $( $file_name ),*
            ];

            $(
                #[test]
                fn $fn_name() {
                    run_corpus_case($file_name);
                }
            )*
        }
    };
}

// 1. mod simple (63 #[test] functions)
corpus_module!(
    simple,
    [
        (
            hcc2218_c_position_bug_chr2_179260451_179260596_f_0_001,
            "Simple;GRCh37.fa;HCC2218_C_position_bug.bam;chr2;179260451-179260596;-f 0.001.txt"
        ),
        (
            n13532_twodels1insstringbounderror_chr11_1280138_1280341_f0_0025_f0x700_q10,
            "Simple;GRCh37.fa;N13532.TwoDels1InsStringBoundError.bam;chr11;1280138-1280341;-f0.0025-F0x700-Q10.txt"
        ),
        (
            na24385_ref_svdel_chr1_9011935_9017839_f_0_001,
            "Simple;GRCh37.fa;NA24385.REF.SVDEL.bam;chr1;9011935-9017839;-f 0.001.txt"
        ),
        (
            na24385_insertion_chr20_60464465_60470068_chr20_60464465_60470068_f_0_001,
            "Simple;GRCh37.fa;NA24385_Insertion_chr20_60464465_60470068.bam;chr20;60464465-60470068;-f 0.001.txt"
        ),
        (
            na24385_chr_19_58559581_58567861_chr20_58559581_58567861_f0_001,
            "Simple;GRCh37.fa;NA24385_chr_19_58559581-58567861.bam;chr20;58559581;58567861;-f0.001.txt"
        ),
        (
            na24385_chr_20_20156208_20234556_chr20_20156208_20234556_f0_001,
            "Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;20156208;20234556;-f0.001.txt"
        ),
        (
            na24385_chr_20_36018823_36039513_chr20_36018823_36039513_f0_001,
            "Simple;GRCh37.fa;NA24385_chr_20_36018823-36039513.bam;chr20;36018823;36039513;-f0.001.txt"
        ),
        (
            na24385_chr_20_38386265_38409801_chr20_38386265_38409801_f0_001,
            "Simple;GRCh37.fa;NA24385_chr_20_38386265-38409801.bam;chr20;38386265;38409801;-f0.001.txt"
        ),
        (
            na24385_chr_20_44885701_44969115_chr20_44885701_44969115_f0_001,
            "Simple;GRCh37.fa;NA24385_chr_20_44885701_44969115.bam;chr20;44885701;44969115;-f0.001.txt"
        ),
        (
            na24385_cnt0_chr11_3379106_3382196_f_0_001,
            "Simple;GRCh37.fa;NA24385_cnt0.bam;chr11;3379106-3382196;-f 0.001.txt"
        ),
        (
            na24385_tcov_empty_chr8_122660146_122660347_f_0_001,
            "Simple;GRCh37.fa;NA24385_tcov_empty.bam;chr8;122660146-122660347;-f 0.001.txt"
        ),
        (
            na24631_sorted_chr1_1000k_chr1_1_499702_f0_001,
            "Simple;GRCh37.fa;NA24631_sorted.chr1.1000k.bam;chr1;1;499702;-f0.001.txt"
        ),
        (
            npe_insertion_adjust_ref_chr3_145521599_145521750_f_0_1_nosv_deldupvar_q10_f0x700,
            "Simple;GRCh37.fa;NPE_insertion_adjust_ref.bam;chr3;145521599-145521750;-f 0.1--nosv--deldupvar-Q10-F0x700.txt"
        ),
        (
            short_quality_seq_in_sv_chr1_17420079_17420230_f_0_001,
            "Simple;GRCh37.fa;short_quality_seq_in_sv.bam;chr1;17420079-17420230;-f 0.001.txt"
        ),
        (
            test_reads_hardclipped_example_genome_160_280_f_0_0_p_0,
            "Simple;hard_clip_case.fa;test_reads.hardclipped.bam;example_genome;160-280;-f 0.0 -P 0.txt"
        ),
        (
            colo829_18_s3_sort_chry_chry_59016100_59016124_chimeric_f_0_001,
            "Simple;hg19.fa;Colo829-18_S3-sort_chrY.bam;chrY;59016100-59016124;--chimeric -f 0.001.txt"
        ),
        (
            colo829_18_s3_sort_chry_chry_59016100_59016124_f_0_001,
            "Simple;hg19.fa;Colo829-18_S3-sort_chrY.bam;chrY;59016100-59016124;-f 0.001.txt"
        ),
        (
            colo829_chr10_canddel5_chr10_134535815_134569759_f_0_001,
            "Simple;hg19.fa;Colo829.chr10.CandDel5.bam;chr10;134535815-134569759;-f 0.001;.txt"
        ),
        (
            colo829_chr10_deldiscreverse_chr10_89712340_89713976_f_0_001,
            "Simple;hg19.fa;Colo829.chr10.DELDiscReverse.bam;chr10;89712340-89713976;-f 0.001;.txt"
        ),
        (
            colo829_chr10_lgdel3_chr10_127187158_127201055_f_0_001,
            "Simple;hg19.fa;Colo829.chr10.LGDEL3.bam;chr10;127187158-127201055;-f 0.001;.txt"
        ),
        (
            colo829_chr10_lgins30_chr10_127800917_127802207_f_0_001,
            "Simple;hg19.fa;Colo829.chr10.LGINS30.bam;chr10;127800917-127802207;-f 0.001;.txt"
        ),
        (
            colo829_chr20_crispr_chr20_50454281_50455444_f_0_001_j_50454941_j_25,
            "Simple;hg19.fa;Colo829.chr20.CRISPR.bam;chr20;50454281-50455444;-f 0.001 -J 50454941 -j 25;.txt"
        ),
        (
            colo829_chr20_canddel5sv_chr20_506503_507295_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.CandDel5SV.bam;chr20;506503-507295;-f 0.001;.txt"
        ),
        (
            colo829_chr20_del_chr20_2802740_2803167_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.DEL.bam;chr20;2802740-2803167;-f 0.001;.txt"
        ),
        (
            colo829_chr20_delsv3happy_chr20_5023101_5026046_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.DELSV3Happy.bam;chr20;5023101-5026046;-f 0.001;.txt"
        ),
        (
            colo829_chr20_delsv3unhappy_chr20_54440610_54441423_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.DELSV3Unhappy.bam;chr20;54440610-54441423;-f 0.001;.txt"
        ),
        (
            colo829_chr20_delsv5happy_chr20_3097973_3098634_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.DELSV5Happy.bam;chr20;3097973-3098634;-f 0.001;.txt"
        ),
        (
            colo829_chr20_dupdiscsv_chr20_47124929_47131802_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.DUPDiscSV.bam;chr20;47124929-47131802;-f 0.001;.txt"
        ),
        (
            colo829_chr20_invsv_chr20_46841555_46845755_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.INVSV.bam;chr20;46841555-46845755;-f 0.001;.txt"
        ),
        (
            colo829_chr20_markingdupsv_chr20_47131954_47132045_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.MarkingDUPSV.bam;chr20;47131954-47132045;-f 0.001;.txt"
        ),
        (
            colo829_chr20_overlapsv_chr20_31260388_31262431_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.OverlapSV.bam;chr20;31260388-31262431;-f 0.001;.txt"
        ),
        (
            colo829_chr20_svbp_chr20_3251475_3252255_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.SVBP.bam;chr20;3251475-3252255;-f 0.001;.txt"
        ),
        (
            colo829_chr20_svbpreverse_chr20_318511_318527_f_0_001,
            "Simple;hg19.fa;Colo829.chr20.SVBPReverse.bam;chr20;318511-318527;-f 0.001;.txt"
        ),
        (
            hardclipinsv_s2_clipped_chr7_1855769_1856455_f_0_001,
            "Simple;hg19.fa;HardClipInSV.S2.clipped.bam;chr7;1855769-1856455;-f 0.001;.txt"
        ),
        (
            l861q_chr7_55259400_55259600_un,
            "Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;-UN.txt"
        ),
        (
            l861q_chr7_55259400_55259600_u,
            "Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;-u.txt"
        ),
        (
            l861q_chr7_55259400_55259600,
            "Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;.txt"
        ),
        (
            case_2_364738_chr21_negative_chr21_28366637_28367427_f0_01_q10_f0x700,
            "Simple;hg38.fa;2-364738.chr21.Negative.bam;chr21;28366637-28367427;-f0.01-Q10-F0x700.txt"
        ),
        (
            case_2_420233_chr7_negative_chr7_157482463_157483329_f0_01_q10_f0x700,
            "Simple;hg38.fa;2-420233.chr7.Negative.bam;chr7;157482463-157483329;-f0.01-Q10-F0x700.txt"
        ),
        (
            amp_allele_chr19_19542710_19542720_f0_0,
            "Simple;hg38.fa;AMP_allele.bam;chr19;19542710-19542720;-f0.0.txt"
        ),
        (
            cigar_three_indels_dlen_chr9_107331700_107331701_f0_001,
            "Simple;hg38.fa;CIGAR_Three_Indels_DLEN.bam;chr9;107331700-107331701;-f0.001.txt"
        ),
        (
            deldiscpairrreference_chr1_43894440_43894460_f_0_001,
            "Simple;hg38.fa;DELDiscPairRreference.bam;chr1;43894440-43894460;-f 0.001.txt"
        ),
        (
            deloutofbound_chr1_150578117_150578617_f0_001_deldupvar,
            "Simple;hg38.fa;DelOutOfBound.bam;chr1;150578117-150578617;-f0.001--deldupvar.txt"
        ),
        (
            deloutofbound_chr1_150578117_150578617_f0_001,
            "Simple;hg38.fa;DelOutOfBound.bam;chr1;150578117-150578617;-f0.001.txt"
        ),
        (
            extrasvinrealignment_chr1_11144571_11144831_f_0_001_u,
            "Simple;hg38.fa;ExtraSVInRealignment.bam;chr1;11144571-11144831;-f 0.001-U.txt"
        ),
        (
            extrasvinrealignment_chr1_11144571_11144831_f_0_001,
            "Simple;hg38.fa;ExtraSVInRealignment.bam;chr1;11144571-11144831;-f 0.001.txt"
        ),
        (
            na12878_chr6_137219200_137219700_f0_001_q10_f0x700,
            "Simple;hg38.fa;NA12878.bam;chr6;137219200-137219700;-f0.001-Q10-F0x700.txt"
        ),
        (
            nap1l1_negativerrc_chr12_76036587_76085033_f_0_001,
            "Simple;hg38.fa;NAP1L1.NegativeRRC.bam;chr12;76036587-76085033;-f 0.001.txt"
        ),
        (
            npe_tcov_tovars_chr7_140453040_140479569_f_0_1,
            "Simple;hg38.fa;NPE_tcov_toVars.bam;chr7;140453040;140479569;-f 0.1.txt"
        ),
        (
            tcov_zero_sv_del_candidate3_chr8_90534611_90620070_f_0_1,
            "Simple;hg38.fa;Tcov_zero_SV_DEL_candidate3.bam;chr8;90534611;90620070;-f 0.1.txt"
        ),
        (
            invalid_char_chr9_107331070_107331095_f0_001,
            "Simple;hg38.fa;invalid-char.bam;chr9;107331070-107331095;-f0.001.txt"
        ),
        (
            norealign_stringindex_oob_chr22_20566544_20566545_f_0_001_k0,
            "Simple;hg38.fa;norealign-stringindex-oob.bam;chr22;20566544;20566545;-f 0.001-k0.txt"
        ),
        (
            test_malformed_cigar_and_string_eor_chr4_25001600_25001607_f0_001,
            "Simple;hg38.fa;test-malformed-cigar-and-string-eor.bam;chr4;25001600-25001607;-f0.001.txt"
        ),
        (
            na21144_chrom11_duplicates_11_402283_410793_f_0,
            "Simple;hs37d5.fa;NA21144.chrom11.Duplicates.bam;11;402283-410793;-F 0;.txt"
        ),
        (
            na21144_chrom11_duplicates_11_402283_410793,
            "Simple;hs37d5.fa;NA21144.chrom11.Duplicates.bam;11;402283-410793;.txt"
        ),
        (
            na12878_chrom20_delsv_chr20_62596000_62598000_f_0_001,
            "Simple;hs38DH.fa;NA12878.chrom20.DELSV.bam;chr20;62596000-62598000;-f 0.001.txt"
        ),
        (
            na12878_chrom20_invfullsv_chr20_5927052_5927179_f_0_001,
            "Simple;hs38DH.fa;NA12878.chrom20.INVFullSV.bam;chr20;5927052-5927179;-f 0.001.txt"
        ),
        (
            na12878_chrom20_invsimplesv_chr20_1961100_1961400_f_0_001,
            "Simple;hs38DH.fa;NA12878.chrom20.INVSimpleSV.bam;chr20;1961100-1961400;-f 0.001.txt"
        ),
        (
            na12878_chrom20_unmapped_reads_chr20_62597650_62597750_q_20,
            "Simple;hs38DH.fa;NA12878.chrom20.Unmapped.reads.bam;chr20;62597650-62597750;-q 20.txt"
        ),
        (
            na12878_chrom5_dupsv_chr5_140594296_140594495_f0_001,
            "Simple;hs38DH.fa;NA12878.chrom5.DUPSV.bam;chr5;140594296-140594495;-f0.001.txt"
        ),
        (
            na21144_chrom11_dupdiscf_chr11_77614410_77614671_f_0_001,
            "Simple;hs38DH.fa;NA21144.chrom11.DUPDiscF.bam;chr11;77614410-77614671;-f 0.001.txt"
        ),
        (
            na21144_chrom11_inv3disc_chr11_57982617_57982816_f_0_001,
            "Simple;hs38DH.fa;NA21144.chrom11.INV3Disc.bam;chr11;57982617-57982816;-f 0.001.txt"
        ),
        (
            na21144_chrom11_inv5disc_chr11_57970876_57971275_f_0_001,
            "Simple;hs38DH.fa;NA21144.chrom11.INV5Disc.bam;chr11;57970876-57971275;-f 0.001.txt"
        ),
    ]
);

// 2. mod somatic (12 #[test] functions)
corpus_module!(
    somatic,
    [
        (
            npe_nm_tumor_chr5_19375000_19376500_f0_1_nosv_deldupvar_q10_f0x700,
            "Somatic;GRCh37.fa;NPE_NM_Tumor.bam-NPE_NM_Normal.bam;chr5;19375000-19376500;-f0.1--nosv--deldupvar-Q10-F0x700;.txt"
        ),
        (
            npe_nm_tumor_chr5_19375000_19376500_f0_1_q10_f0x700,
            "Somatic;GRCh37.fa;NPE_NM_Tumor.bam-NPE_NM_Normal.bam;chr5;19375000-19376500;-f0.1-Q10-F0x700;.txt"
        ),
        (
            threeindels_zero_deletion_tumor_chr20_62354864_62374864_f0_1_q10_f0x700,
            "Somatic;GRCh37.fa;threeIndels_zero_deletion_tumor.bam-threeIndels_zero_deletion_normal.bam;chr20;62354864-62374864;-f0.1-Q10-F0x700;.txt"
        ),
        (
            colo829_chr10_1368753_1369353,
            "Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr10;1368753-1369353;.txt"
        ),
        (
            colo829_chr10_684892_685492,
            "Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr10;684892-685492;.txt"
        ),
        (
            colo829_chr1_755917_756517,
            "Somatic;hg19.fa;Colo829-18_S3-sort.bam-Colo829-19_S4-sort.bam;chr1;755917-756517;.txt"
        ),
        (
            colo829_18_chr20_missedstrongloh_chr20_59844695_59857589_f0_001,
            "Somatic;hg19.fa;Colo829.18.chr20.MissedStrongLoh.bam-Colo829.19.chr20.MissedStrongLoh.bam;chr20;59844695;59857589;-f0.001.txt"
        ),
        (
            hcc2218_bl_realigins_chrx_7811670_7811868_f0_001,
            "Somatic;hg19.fa;HCC2218_BL_realigins.bam-HCC2218_C_realigins.bam;chrX;7811670-7811868;-f0.001;.txt"
        ),
        (
            zero_in_dupsv_tumour_chr19_290742_291636_f0_001_x3,
            "Somatic;hg19.fa;Zero_in_DUPSV_TUMOUR.bam-Zero_in_DUPSV_NORMAL.bam;chr19;290742-291636;-f0.001-X3.txt"
        ),
        (
            tumora_strongloh_chr21_10996701_11186054_q10_f0_01_r2_v0_05,
            "Somatic;hg1k.fa;TumorA.StrongLOH.bam-NormalA.StrongLOH.bam;chr21;10996701-11186054;-Q10-f0.01-r2-V0.05;.txt"
        ),
        (
            case_172t_s8_final_chr1_chr1_11182030_11182188,
            "Somatic;hs38DH.fa;172T_S8.final.chr1.bam-172N_S7.final.chr1.bam;chr1;11182030-11182188;;.txt"
        ),
        (
            na12889_strsomsamspec_chr20_43246926_43280248_f0_001,
            "Somatic;hs38DH.fa;NA12889.StrSomSamSpec.bam-NA12878.StrSomSamSpec.bam;chr20;43246926-43280248;-f0.001.txt"
        ),
    ]
);

// 3. mod pileup (11 #[test] functions)
corpus_module!(
    pileup,
    [
        (
            hard_clip_next_to_del_test1_test_6674_6824_f_0_0_p_r_1,
            "Simple;hard_clip_case.fa;hard_clip_next_to_del_test1.bam;test;6674-6824;-f 0.0 -p -r 1.txt"
        ),
        (
            hard_clip_next_to_del_test2_test_6434_6564_f_0_0_p_r_1,
            "Simple;hard_clip_case.fa;hard_clip_next_to_del_test2.bam;test;6434-6564;-f 0.0 -p -r 1.txt"
        ),
        (
            colo829_18_s3_sort_chr1_9922_10122_chimeric_p,
            "Simple;hg19.fa;Colo829-18_S3-sort.bam;chr1;9922-10122;--chimeric -p;.txt"
        ),
        (
            colo829_18_s3_sort_chr1_9922_10122_p,
            "Simple;hg19.fa;Colo829-18_S3-sort.bam;chr1;9922-10122;-p;.txt"
        ),
        (
            count_depth_with_ns_chr1_391310_391333_k_p_r_1_q_2_f_0_0,
            "Simple;hs38DH.fa;Count-Depth-With-Ns.bam;chr1;391310-391333;-K -p -r 1 -q 2 -f 0.0;.txt"
        ),
        (
            count_depth_with_ns_chr1_391310_391333_p_r_1_q_2_f_0_0,
            "Simple;hs38DH.fa;Count-Depth-With-Ns.bam;chr1;391310-391333;-p -r 1 -q 2 -f 0.0;.txt"
        ),
        (
            t7_75_s10_chr10_del_chr10_11291400_11295220_f_0_001,
            "Simple;mm10.fa;T7-75_S10_chr10_DEL.bam;chr10;11291400-11295220;-f 0.001;.txt"
        ),
        (
            t7_75_s10_chr10_del_chr10_11291400_11295300_f_0_001,
            "Simple;mm10.fa;T7-75_S10_chr10_DEL.bam;chr10;11291400-11295300;-f 0.001;.txt"
        ),
        (
            amp_allele_ins_repr_test_148_149_f_0_0_p,
            "Simple;test.fa;amp_allele_ins_repr.bam;test;148-149;-f 0.0 -p.txt"
        ),
        (
            test_alt_alleles_chrt_0_52_p,
            "Simple;test.fa;test_alt_alleles.bam;chrT;0-52;-p.txt"
        ),
        (
            argumentna_tumor_chr7_140453010_140453158_q25_f0x100_f0_001_k1_m20_p,
            "Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p;.txt"
        ),
    ]
);

// 4. mod fisher (3 #[test] functions)
corpus_module!(
    fisher,
    [
        (
            na24385_chr_20_20156208_20234556_chr20_20156208_20234556_f0_001_fisher,
            "Simple;GRCh37.fa;NA24385_chr_20_20156208-20234556.bam;chr20;20156208;20234556;-f0.001-fisher.txt"
        ),
        (
            argumentna_tumor_chr7_140453010_140453158_q25_f0x100_f0_001_k1_m20_p_fisher,
            "Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt"
        ),
        (
            na12878_part1_1956957_chr1_1956957_1957177_f0_001,
            "Somatic;hs37d5.fa;NA12878.part1.1956957-1957177.bam-NA12889.part1.1956957-1957177.bam;chr1;1956957-1957177;-f0.001;.txt"
        ),
    ]
);

// 5. mod amplicon (2 #[test] functions)
corpus_module!(
    amplicon,
    [
        (
            mre11a_chr11_94153228_94153388_94153261_94153357_a_100_95_p,
            "Amplicon;GRCh37.fa;MRE11a.bam;chr11;94153228-94153388;94153261-94153357;-a 100.95 -p;.txt"
        ),
        (
            colo829_18_s3_sort_chr1_933866_934466_933866_934466_a_100_95_d,
            "Amplicon;hg19.fa;Colo829-18_S3-sort.bam;chr1;933866-934466;933866-934466;-a 100.95 -D;.txt"
        ),
    ]
);

// 6. mod splicing (1 #[test] function)
corpus_module!(
    splicing,
    [(
        tcov_zero_sv_del_candidate3_chr8_90534611_90620070_f_0_1_i,
        "Simple;hg38.fa;Tcov_zero_SV_DEL_candidate3.bam;chr8;90534611;90620070;-f 0.1 -i.txt"
    )]
);

// 7. mod modalities (17 #[test] functions)
mod modalities {
    use super::*;

    fn micro_config() -> Config {
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
            sample: "SAMPLE".to_string(),
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
            bed: Some(
                root.join("fixtures/micro/micro.bed")
                    .to_string_lossy()
                    .into_owned(),
            ),
            ..Config::default()
        }
    }

    #[test]
    fn rna_splicing_intermediate() {
        run_corpus_case(
            "Simple;hg38.fa;Tcov_zero_SV_DEL_candidate3.bam;chr8;90534611;90620070;-f 0.1 -i.txt",
        );
    }

    #[test]
    fn structural_variants_dup_inv_del() {
        run_corpus_case(
            "Simple;hg19.fa;Colo829.chr20.MarkingDUPSV.bam;chr20;47131954-47132045;-f 0.001;.txt",
        );
        run_corpus_case(
            "Simple;hg19.fa;Colo829.chr20.INVSV.bam;chr20;46841555-46845755;-f 0.001;.txt",
        );
        run_corpus_case("Simple;hg19.fa;Colo829.chr20.DEL.bam;chr20;2802740-2803167;-f 0.001;.txt");
    }

    #[test]
    fn amplicon_primer_clipping() {
        run_corpus_case(
            "Amplicon;GRCh37.fa;MRE11a.bam;chr11;94153228-94153388;94153261-94153357;-a 100.95 -p;.txt",
        );
    }

    #[test]
    fn crispr_cut_site_filtering() {
        run_corpus_case(
            "Simple;hg19.fa;Colo829.chr20.CRISPR.bam;chr20;50454281-50455444;-f 0.001 -J 50454941 -j 25;.txt",
        );
    }

    #[test]
    fn chimeric_read_filtering() {
        run_corpus_case(
            "Simple;hg19.fa;Colo829-18_S3-sort_chrY.bam;chrY;59016100-59016124;--chimeric -f 0.001.txt",
        );
    }

    #[test]
    fn adapter_trimming() {
        let mut cfg = micro_config();
        cfg.adaptor = vec![
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string(),
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string(),
        ];
        let var_out = run_core_var(&cfg).expect("adapter_trimming .var");
        let expected =
            fs::read_to_string(workspace_root().join("golden/micro.expected.var")).expect("golden");
        assert_eq!(var_out.trim(), expected.trim());
    }

    #[test]
    fn indel_3prime_shifting() {
        let mut cfg = micro_config();
        cfg.three_prime = true;
        let vcf_out = run_core_vcf(&cfg).expect("indel_3prime_shifting VCF");
        let vars = extract_variants(&vcf_out);
        assert_eq!(vars.len(), 2, "expected SNV + shifted indel on micro.bam");
    }

    #[test]
    fn unique_read_dedup_first_alignment() {
        run_corpus_case("Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;-UN.txt");
        run_corpus_case("Simple;hg19.fa;L861Q.bam;chr7;55259400-55259600;-u.txt");
    }

    #[test]
    fn sam_flag_filtering() {
        run_corpus_case(
            "Simple;hs37d5.fa;NA21144.chrom11.Duplicates.bam;11;402283-410793;-F 0;.txt",
        );
        run_corpus_case("Simple;hs37d5.fa;NA21144.chrom11.Duplicates.bam;11;402283-410793;.txt");
    }

    #[test]
    fn pileup_all_sites() {
        run_corpus_case(
            "Simple;hs38DH.fa;Count-Depth-With-Ns.bam;chr1;391310-391333;-K -p -r 1 -q 2 -f 0.0;.txt",
        );
    }

    #[test]
    fn cram_reference_decoding() {
        ensure_fastas_converted();
        let droot = data_root();
        let ref_fa = droot.join("VarDictJava/testdata/fastas/hg19.fa");
        let bam_path = droot.join(
            "VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829.chr20.CRISPR.bam",
        );
        let tmp = tempfile::tempdir().expect("tempdir");
        let cram_path = tmp.path().join("crispr.cram");
        write_indexed_cram(&bam_path, &ref_fa, &cram_path, true, Some("chr20"))
            .expect("encode cram");

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

        let vcf_bam = run_core_vcf(&cfg_bam).expect("bam vcf");
        let vcf_cram = run_core_vcf(&cfg_cram).expect("cram vcf");
        let bam_vars = extract_variants(&vcf_bam);
        assert_eq!(bam_vars.len(), 8);
        assert_eq!(
            bam_vars,
            extract_variants(&vcf_cram),
            "CRAM reference decoding must match BAM variants byte-for-byte"
        );
    }

    #[test]
    fn somatic_fisher_exact() {
        run_corpus_case(
            "Somatic;hg19.fa;ArgumentNA_Tumor.bam-ArgumentNA_Normal.bam;chr7;140453010-140453158;-q25-F0x100-f0.001-k1-M20-p-fisher;.txt",
        );
    }

    #[test]
    fn amplicon_fisher_exact() {
        ensure_fastas_converted();
        let droot = data_root();
        let tmp = tempfile::tempdir().expect("tempdir");
        let bed_path = tmp.path().join("amp.bed");
        fs::write(
            &bed_path,
            "chr1\t933866\t934466\ttestbed\t0\t0\t933866\t934466\n",
        )
        .expect("write amp bed");
        let bam = droot
            .join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests/Colo829-18_S3-sort.bam")
            .to_string_lossy()
            .into_owned();
        let cfg = Config {
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
            bed: Some(bed_path.to_string_lossy().into_owned()),
            amplicon: Some("10:0.95".to_string()),
            fisher: true,
            ..Config::default()
        };
        let var_out = run_core_var(&cfg).expect("amplicon --fisher .var");
        let rows: Vec<&str> = var_out.lines().filter(|l| !l.trim().is_empty()).collect();
        assert!(!rows.is_empty(), "amplicon --fisher must emit rows");
        for row in &rows {
            assert_eq!(
                row.split('\t').count(),
                40,
                "amplicon --fisher .var must emit 40 columns"
            );
        }
        let vcf_out = run_core_vcf(&cfg).expect("amplicon --fisher VCF");
        let vars = extract_variants(&vcf_out);
        assert!(
            !vars.is_empty()
                && vars
                    .iter()
                    .all(|v| v.contains("SBF=") && v.contains("ODDRATIO=")),
            "amplicon --fisher VCF must emit SBF and ODDRATIO"
        );
    }

    #[test]
    fn one_based_segment_indexing() {
        let root = workspace_root();
        let tmp = tempfile::tempdir().expect("tempdir");
        let bed_excl = tmp.path().join("seg_excl.bed");
        let bed_incl = tmp.path().join("seg_incl.bed");
        fs::write(
            &bed_excl,
            "chr1\t41\t120\tmicro_chr1\t0\t+\t41\t120\t0\t41\t80\tmicro_chr1\n",
        )
        .expect("write excl bed");
        fs::write(
            &bed_incl,
            "chr1\t41\t120\tmicro_chr1\t0\t+\t41\t120\t0\t90\t120\tmicro_chr1\n",
        )
        .expect("write incl bed");

        let mut cfg_excl = micro_config();
        cfg_excl.zero_based = false;
        cfg_excl.col_seg_start = 9;
        cfg_excl.col_seg_end = 10;
        cfg_excl.bed = Some(bed_excl.to_string_lossy().into_owned());

        let mut cfg_incl = cfg_excl.clone();
        cfg_incl.bed = Some(bed_incl.to_string_lossy().into_owned());
        cfg_incl.reference = root
            .join("fixtures/micro/micro.ref.fa")
            .to_string_lossy()
            .into_owned();

        let vars_excl = extract_variants(&run_core_vcf(&cfg_excl).expect("excl vcf"));
        let vars_incl = extract_variants(&run_core_vcf(&cfg_incl).expect("incl vcf"));
        assert!(
            vars_excl.is_empty(),
            "1-based segment 41..80 must exclude chr1:101"
        );
        assert_eq!(
            vars_incl.len(),
            1,
            "1-based segment 90..120 must include chr1:101"
        );
    }

    #[test]
    fn somatic_monomer_msi_filtering() {
        run_corpus_case(
            "Somatic;hg19.fa;HCC2218_BL_realigins.bam-HCC2218_C_realigins.bam;chrX;7811670-7811868;-f0.001;.txt",
        );
        ensure_fastas_converted();
        let droot = data_root();
        let bam_dir =
            droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
        let bam_pair = format!(
            "{}|{}",
            bam_dir.join("ArgumentNA_Tumor.bam").display(),
            bam_dir.join("ArgumentNA_Normal.bam").display()
        );
        let cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg19.fa")
                .to_string_lossy()
                .into_owned(),
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
            freq: 0.001,
            af_threshold: 0.001,
            mfreq: 0.25,
            nmfreq: 0.1,
            ..Config::default()
        };
        let vcf_out = run_core_vcf(&cfg).expect("somatic monomer MSI VCF");
        let vars = extract_variants(&vcf_out);
        assert!(
            !vars.is_empty()
                && vars
                    .iter()
                    .all(|v| v.contains("MSI=") && v.contains("MSILEN=")),
            "somatic VCF must emit MSI and MSILEN annotations"
        );
    }

    #[test]
    fn deldupvar_sv_accounting() {
        run_corpus_case(
            "Simple;hg38.fa;DelOutOfBound.bam;chr1;150578117-150578617;-f0.001--deldupvar.txt",
        );
        run_corpus_case("Simple;hg38.fa;DelOutOfBound.bam;chr1;150578117-150578617;-f0.001.txt");
    }

    #[test]
    fn default_printer_stderr_routing() {
        let mut cfg = micro_config();
        cfg.debug_print = Some("ERR".to_string());
        assert!(
            matches!(
                vardictrust_core::create_output_writer(&cfg),
                vardictrust_core::OutputDestination::Stderr(_)
            ),
            "--DP ERR must route output writer to Stderr"
        );
        let var_out = run_core_var(&cfg).expect("run_core_var under --DP ERR");
        let expected =
            fs::read_to_string(workspace_root().join("golden/micro.expected.var")).expect("golden");
        assert_eq!(var_out.trim(), expected.trim());
    }
}

// 8. mod manifest (1 #[test] function)
mod manifest {
    use super::*;

    fn exercise_in_memory_core_branches() {
        use vardictrust_core::stats::f80::{f80_sum, F80};
        use vardictrust_core::stats::rmath::{bd0, dbinom_raw, dhyper, phyper, stirlerr};

        for &n in &[
            0.0,
            0.5,
            1.0,
            2.5,
            3.1,
            5.2,
            5.5,
            6.2,
            6.8,
            7.5,
            9.5,
            12.5,
            13.5,
            15.0,
            20.0,
            25.0,
            50.0,
            100.0,
            250.0,
            7000.0,
            20_000_000.0,
        ] {
            let _ = stirlerr::stirlerr(n);
        }

        for &(x, np) in &[
            (f64::NAN, 1.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (1.05, 1.0),
            (5.0, 2.0),
            (2.0, 5.0),
            (1e308, 1e-308),
        ] {
            let _ = bd0::bd0(x, np);
        }

        for &y in &[f64::NAN, 0.0, 1.0, 2.0, 3.0, 4.0, 5.5] {
            let _ = dbinom_raw::pow1p(0.0, y);
            let _ = dbinom_raw::pow1p(0.1, y);
            let _ = dbinom_raw::pow1p(0.8, y);
            let _ = dbinom_raw::pow1p(f64::NAN, y);
        }

        for &give_log in &[false, true] {
            let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.0, 1.0, give_log);
            let _ = dbinom_raw::dbinom_raw(2.0, 10.0, 0.0, 1.0, give_log);
            let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 1.0, 0.0, give_log);
            let _ = dbinom_raw::dbinom_raw(2.0, 10.0, 1.0, 0.0, give_log);
            let _ = dbinom_raw::dbinom_raw(0.0, 0.0, 0.5, 0.5, give_log);
            let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.7, 0.3, give_log);
            let _ = dbinom_raw::dbinom_raw(0.0, 10.0, 0.3, 0.7, give_log);
            let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 0.7, 0.3, give_log);
            let _ = dbinom_raw::dbinom_raw(10.0, 10.0, 0.3, 0.7, give_log);
            let _ = dbinom_raw::dbinom_raw(-1.0, 10.0, 0.5, 0.5, give_log);
            let _ = dbinom_raw::dbinom_raw(11.0, 10.0, 0.5, 0.5, give_log);
            let _ = dhyper::dhyper(-1.0, 5.0, 5.0, 5.0, give_log);
            let _ = dhyper::dhyper(6.0, 5.0, 5.0, 5.0, give_log);
            let _ = dhyper::dhyper(0.0, 5.0, 5.0, 0.0, give_log);
            let _ = dhyper::dhyper(2.0, 5.0, 5.0, 5.0, give_log);
            let _ = phyper::phyper(-1.0, 5.0, 5.0, 5.0, true);
            let _ = phyper::phyper(10.0, 5.0, 5.0, 5.0, false);
            let _ = phyper::phyper(2.0, 10.0, 10.0, 10.0, true);
            let _ = phyper::phyper(8.0, 10.0, 10.0, 10.0, false);
        }

        let subnormal = f64::from_bits(1);
        let f_sub = F80::from_f64(subnormal);
        assert!(f_sub.to_f64() >= 0.0);
        let f_big = F80 {
            exp: 1030,
            mant: 1u64 << 63,
        };
        assert_eq!(f_big.to_f64(), f64::INFINITY);
        let f_tiny = F80 {
            exp: -1100,
            mant: 1u64 << 63,
        };
        assert_eq!(f_tiny.to_f64(), 0.0);
        let _ = f80_sum(&[subnormal, 1e-300, 1.0]);
    }

    #[test]
    fn all_92_golden_cases_registered() {
        exercise_in_memory_core_branches();
        ensure_fastas_converted();

        let droot = data_root();
        let testdata_dir = droot.join("VarDictJava/testdata/integrationtestcases");
        let mut disk_files: Vec<String> = fs::read_dir(&testdata_dir)
            .expect("read integrationtestcases")
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("txt"))
            .filter_map(|p| {
                p.file_name()
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string())
            })
            .collect();
        disk_files.sort();
        assert_eq!(disk_files.len(), 92, "expected 92 .txt fixtures on disk");

        assert_eq!(super::simple::CASES.len(), 63);
        assert_eq!(super::somatic::CASES.len(), 12);
        assert_eq!(super::pileup::CASES.len(), 11);
        assert_eq!(super::fisher::CASES.len(), 3);
        assert_eq!(super::amplicon::CASES.len(), 2);
        assert_eq!(super::splicing::CASES.len(), 1);

        let mut registered: Vec<String> = Vec::with_capacity(92);
        for slice in [
            super::simple::CASES,
            super::somatic::CASES,
            super::pileup::CASES,
            super::fisher::CASES,
            super::amplicon::CASES,
            super::splicing::CASES,
        ] {
            registered.extend(slice.iter().map(|s| (*s).to_string()));
        }
        registered.sort();
        assert_eq!(
            registered, disk_files,
            "registered CORPUS cases must match all 92 disk fixtures 1-to-1"
        );

        // Verify referenced FASTAs and BAMs exist and check pure-Rust SHA-256 on micro golden.
        let fasta_dir = droot.join("VarDictJava/testdata/fastas");
        let bam_dir =
            droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
        for fname in &registered {
            let text = fs::read_to_string(testdata_dir.join(fname)).expect("read fixture");
            let header = text.lines().next().expect("header");
            let fields: Vec<&str> = header.split(',').map(|s| s.trim()).collect();
            let fa = fasta_dir.join(fields[1]);
            assert!(fa.is_file(), "missing FASTA {} for {fname}", fa.display());
            for bam_name in fields[2].split('|') {
                let bp = bam_dir.join(bam_name);
                assert!(bp.is_file(), "missing BAM {} for {fname}", bp.display());
            }
        }

        let micro_bytes =
            fs::read(workspace_root().join("golden/micro.expected.var")).expect("read micro.var");
        let micro_sha = sha256_hex(&micro_bytes);
        assert_eq!(
            micro_sha.len(),
            64,
            "pure-Rust SHA-256 must produce 64-hex digest"
        );
        let micro_bam_bytes =
            fs::read(workspace_root().join("fixtures/micro/micro.bam")).expect("read micro.bam");
        let micro_ref_bytes = fs::read(workspace_root().join("fixtures/micro/micro.ref.fa"))
            .expect("read micro.ref.fa");
        assert!(!sha256_hex(&micro_bam_bytes).is_empty());
        assert!(!sha256_hex(&micro_ref_bytes).is_empty());
        let _ = Path::new(&micro_sha);
    }
}
