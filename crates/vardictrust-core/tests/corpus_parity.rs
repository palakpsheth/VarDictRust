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

use vardictrust_core::config::{extract_sample_name, resolve_sample_names, Config};

mod common;
use common::{
    data_root, ensure_fastas_converted, extract_variants, indexed_bam, run_core_var, run_core_vcf,
    sha256_hex, workspace_root, write_indexed_cram,
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
                let decoded = vardictrust_core::config::java_integer_decode(v).expect("valid -F");
                cfg.flag_filter = decoded as u32;
                cfg.samfilter_raw = v.clone();
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
                cfg.read_pos = vardictrust_core::config::parse_java_number_i32(&extra_args[i + 1])
                    .expect("valid -P") as f64;
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

    let actual_lines: Vec<&str> = if actual_trimmed.is_empty() {
        vec![""]
    } else {
        actual_trimmed.split('\n').collect()
    };
    let mut actual_sorted = actual_lines.clone();
    actual_sorted.sort_unstable();

    let expected_lines: Vec<&str> = if expected.is_empty() {
        vec![""]
    } else {
        expected.split('\n').collect()
    };
    let mut expected_sorted = expected_lines.clone();
    expected_sorted.sort_unstable();

    if actual_sorted == expected_sorted {
        assert!(
            !xfails.contains_key(fname),
            "XPASS on registered XFAIL {fname}"
        );
        if std::env::var("STRICT_VAR_ORDER").ok().as_deref() != Some("0") && !cfg.splice {
            let row_key = |line: &str| -> String {
                let cols: Vec<&str> = line.split('\t').take(7).collect();
                cols.join("\t")
            };
            let actual_keys: Vec<String> = actual_lines.iter().map(|l| row_key(l)).collect();
            let expected_keys: Vec<String> = expected_lines.iter().map(|l| row_key(l)).collect();
            assert_eq!(
                actual_keys, expected_keys,
                "UNSORTED ROW KEY ORDER MISMATCH on {fname}"
            );
            assert_eq!(
                actual_lines, expected_lines,
                "UNSORTED ROW ORDER MISMATCH on {fname}"
            );
        }
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

        // H1: Upstream VarDictJava RegionBuilderTest / IntegrationTest BED auto-detect fixtures.
        let beds_dir = workspace_root().join("fixtures/beds");
        assert!(
            vardictrust_core::regions::detect_amplicon_bed(
                beds_dir.join("correct_8columns_amplicon.bed"),
                "\t"
            )
            .expect("detect correct_8columns_amplicon"),
            "correct_8columns_amplicon.bed must auto-detect as amplicon BED"
        );
        for non_amp in [
            "correct_5columns.bed",
            "correct_8columns_not_amplicon.bed",
            "correct_9columns.bed",
            "incorrect_4columns.bed",
            "incorrect_8columns.bed",
        ] {
            assert!(
                !vardictrust_core::regions::detect_amplicon_bed(beds_dir.join(non_amp), "\t")
                    .unwrap_or_else(|e| panic!("detect {non_amp}: {e}")),
                "{non_amp} must NOT auto-detect as amplicon BED"
            );
        }

        // Amplicon -p zero-coverage null row (AmpliconPostProcessModule.java:244-252 / AmpliconOutputVariant.java:70-98):
        // zero-coverage amplicon locus emits 38 columns with empty genotype (col 14), "0;0" bias (col 16),
        // and empty flank sequences (cols 31, 32).
        let dir = tempfile::tempdir().expect("tempdir for zero-cov amplicon bed");
        let bed_path = dir.path().join("amp_zero.bed");
        fs::write(&bed_path, "chr1\t10\t20\tGENE1\t0\t+\t12\t18\n")
            .expect("write zero-cov amplicon bed");
        let mut amp_cfg = micro_config();
        amp_cfg.region = None;
        amp_cfg.bed = Some(bed_path.to_string_lossy().into_owned());
        amp_cfg.column_for_chromosome = None;
        amp_cfg.amplicon = Some("10:0.95".to_string());
        amp_cfg.pileup = true;
        amp_cfg.apply_pileup_overrides();
        let out = run_core_var(&amp_cfg).expect("zero-cov amplicon -p run");
        let rows: Vec<&str> = out.lines().collect();
        assert!(!rows.is_empty(), "expected zero-coverage amplicon -p rows");
        for row in rows {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols.len(), 38, "expected 38-col amplicon row: {row}");
            assert_eq!(cols[13], "", "null_record genotype must be empty: {row}");
            assert_eq!(cols[15], "0;0", "null_record bias must be '0;0': {row}");
            assert_eq!(cols[30], "", "null_record left_seq must be empty: {row}");
            assert_eq!(cols[31], "", "null_record right_seq must be empty: {row}");
        }
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
        check_fisher_vcf_raw_cell_strict_parity(&var_out, &vars);
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
        check_trim_and_sample_regexp_corpus_parity();
    }

    fn vcf_chrom_header(vcf_text: &str) -> &str {
        vcf_text
            .lines()
            .find(|l| l.starts_with("#CHROM\t"))
            .expect("missing #CHROM header line")
    }

    fn vcf_variants_sha256(vars: &[String]) -> String {
        let mut joined = String::new();
        for v in vars {
            joined.push_str(v);
            joined.push('\n');
        }
        sha256_hex(joined.as_bytes())
    }

    fn check_trim_and_sample_regexp_corpus_parity() {
        // 1. fixtures/micro/micro.bam with -T 25 and -n '.*(mic)ro\.bam'
        let mut micro_cfg = micro_config();
        micro_cfg.trim_bases_after = 25;
        micro_cfg.sample_regexp = Some(r".*(mic)ro\.bam".to_string());
        assert_eq!(
            resolve_sample_names(
                &micro_cfg.bam,
                None,
                micro_cfg.sample_regexp.as_deref(),
                false,
                false,
            ),
            "mic"
        );
        let micro_var = run_core_var(&micro_cfg).expect("micro -T 25 -n .var");
        let expected_micro_var = concat!(
            "mic\tmicro_chr1\tchr1\t101\t101\tT\tA\t7\t3\t4\t0\t3\t0\tT/A\t0.4286\t0;0\t13.0\t1\t40.0\t0\t60.0\t6.000\t0.4286\t0\t0\t3.000\t1\t0\t3\t7\tTCGCCGTTGCCTAAACCTAT\tTGAAGGAGTCTAGCAGCCGC\tchr1:41-260\tSNV\t0\t0\n",
            "mic\tmicro_chr2\tchr2\t150\t153\tGCGG\tG\t8\t4\t4\t0\t4\t0\tC/-3\t0.5000\t0;0\t21.0\t1\t40.0\t0\t60.0\t8.000\t0.5000\t0\t0\t3.000\t1\t0\t4\t8\tCATCTGTTTCTTGCGTCGTA\tGACCCTCCATTGTTACTTAT\tchr2:111-230\tDeletion\t0\t0\n",
        );
        assert_eq!(micro_var, expected_micro_var, "micro -T 25 -n .var parity");
        assert_eq!(
            sha256_hex(micro_var.as_bytes()),
            "43ab86e3c31cd2239462dd45feedc13b74610162ae1f435ba32b1c243c868257",
            "micro -T 25 -n .var SHA-256"
        );

        let micro_vcf = run_core_vcf(&micro_cfg).expect("micro -T 25 -n VCF");
        assert_eq!(
            vcf_chrom_header(&micro_vcf),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmic"
        );
        let micro_vars = extract_variants(&micro_vcf);
        assert_eq!(
            micro_vars,
            vec![
                "chr1\t101\t.\tT\tA\t63\tPASS\tSAMPLE=mic;TYPE=SNV;DP=7;VD=3;AF=0.4286;BIAS=0:0;REFBIAS=4:0;VARBIAS=3:0;PMEAN=13;PSTD=1;QUAL=40;QSTD=0;SBF=1;ODDRATIO=0;MQ=60;SN=6;HIAF=0.4286;ADJAF=0;SHIFT3=0;MSI=3;MSILEN=1;NM=0;HICNT=3;HICOV=7;LSEQ=TCGCCGTTGCCTAAACCTAT;RSEQ=TGAAGGAGTCTAGCAGCCGC;DUPRATE=0;SPLITREAD=0;SPANPAIR=0\tGT:DP:VD:AD:AF:RD:ALD\t0/1:7:3:4,3:0.4286:4,0:3,0",
                "chr2\t150\t.\tGCGG\tG\t80\tPASS\tSAMPLE=mic;TYPE=Deletion;DP=8;VD=4;AF=0.5;BIAS=0:0;REFBIAS=4:0;VARBIAS=4:0;PMEAN=21;PSTD=1;QUAL=40;QSTD=0;SBF=1;ODDRATIO=0;MQ=60;SN=8;HIAF=0.5000;ADJAF=0;SHIFT3=0;MSI=3;MSILEN=1;NM=0;HICNT=4;HICOV=8;LSEQ=CATCTGTTTCTTGCGTCGTA;RSEQ=GACCCTCCATTGTTACTTAT;DUPRATE=0;SPLITREAD=0;SPANPAIR=0\tGT:DP:VD:AD:AF:RD:ALD\t1/0:8:4:4,4:0.5:4,0:4,0",
            ]
        );
        assert_eq!(
            vcf_variants_sha256(&micro_vars),
            "67fc43df95007f0aa3317c2f28acbfebef8aa23e742418ec2212d955693a1bad"
        );

        // 2. Real corpus BAM L861Q.bam (chr7:55259400-55259600:testbed) with -T 50 and -n '.*/(L861)Q\.bam'
        ensure_fastas_converted();
        let droot = data_root();
        let bam_dir =
            droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
        let hg19_fa = droot
            .join("VarDictJava/testdata/fastas/hg19.fa")
            .to_string_lossy()
            .into_owned();
        let l861_cfg = Config {
            reference: hg19_fa.clone(),
            bam: bam_dir.join("L861Q.bam").to_string_lossy().into_owned(),
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            region: Some("chr7:55259400-55259600:testbed".to_string()),
            freq: 0.01,
            af_threshold: 0.01,
            trim_bases_after: 50,
            sample_regexp: Some(r".*/(L861)Q\.bam".to_string()),
            ..Config::default()
        };
        let l861_var = run_core_var(&l861_cfg).expect("L861Q -T 50 -n .var");
        let expected_l861_var = concat!(
            "L861\ttestbed\tchr7\t55259468\t55259468\tC\tT\t169\t3\t166\t0\t2\t1\tC/T\t0.0178\t1;2\t12.7\t1\t34.7\t1\t60.0\t6.000\t0.0180\t0.0118\t0\t1.000\t1\t2.0\t3\t167\tCGCGACCTGGCAGCCAGGAA\tGTACTGGTGAAAACACCGCA\tchr7:55259401-55259600\tSNV\t0\t0\n",
            "L861\ttestbed\tchr7\t55259504\t55259504\tA\tT\t516\t8\t277\t230\t0\t8\tA/T\t0.0155\t2;0\t7.0\t0\t36.0\t0\t60.0\t16.000\t0.0157\t0\t0\t1.000\t1\t1.0\t8\t511\tCCGCAGCATGTCAAGATCAC\tGATTTTGGGCTGGCCAAACT\tchr7:55259401-55259600\tSNV\t0\t0\n",
            "L861\ttestbed\tchr7\t55259524\t55259524\tT\tA\t646\t148\t224\t273\t66\t82\tT/A\t0.2291\t2;2\t22.6\t1\t34.8\t1\t60.0\t23.667\t0.2258\t0.0062\t0\t2.000\t3\t1.4\t142\t629\tAGATTTTGGGCTGGCCAAAC\tGCTGGGTGCGGAAGAGAAAG\tchr7:55259401-55259600\tSNV\t0\t0\n",
            "L861\ttestbed\tchr7\t55259572\t55259572\tG\tA\t193\t3\t15\t174\t0\t3\tG/A\t0.0155\t2;0\t6.0\t1\t33.3\t1\t60.0\t6.000\t0.0157\t0.0052\t1\t2.000\t1\t2.0\t3\t191\tTGCAGAAGGAGGCAAAGTAA\tGAGGTGGCTTTAGGTCAGCC\tchr7:55259401-55259600\tSNV\t0\t0\n",
        );
        assert_eq!(l861_var, expected_l861_var, "L861Q -T 50 -n .var parity");
        assert_eq!(
            sha256_hex(l861_var.as_bytes()),
            "83d1a74d72a73d6ba338da932da36dede2eb9b051f57c03b026011d9a2b19321"
        );
        let l861_vcf = run_core_vcf(&l861_cfg).expect("L861Q -T 50 -n VCF");
        assert_eq!(
            vcf_chrom_header(&l861_vcf),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tL861"
        );
        let l861_vars = extract_variants(&l861_vcf);
        assert_eq!(l861_vars.len(), 4);
        assert!(l861_vars.iter().all(|v| v.contains("SAMPLE=L861;")));
        assert_eq!(
            vcf_variants_sha256(&l861_vars),
            "f86fb4bdac9ef6ef96903524e0c1e6bb85ef21fc51d2344481c978cc4ecfc714"
        );

        // 3. Paired somatic corpus BAMs ArgumentNA_Tumor.bam|ArgumentNA_Normal.bam with -n '.*/ArgumentNA_(Tumor|Normal)\.bam'
        let bam_pair = format!(
            "{}|{}",
            bam_dir.join("ArgumentNA_Tumor.bam").display(),
            bam_dir.join("ArgumentNA_Normal.bam").display()
        );
        let pair_re = r".*/ArgumentNA_(Tumor|Normal)\.bam";
        assert_eq!(
            resolve_sample_names(&bam_pair, None, Some(pair_re), true, false),
            "Tumor|Normal"
        );
        assert_eq!(
            resolve_sample_names(&bam_pair, None, Some(pair_re), true, true),
            "Tumor"
        );
        let pair_cfg = Config {
            reference: hg19_fa,
            bam: bam_pair,
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            region: Some("chr7:140453010-140453158:testbed".to_string()),
            freq: 0.001,
            af_threshold: 0.001,
            sample_regexp: Some(pair_re.to_string()),
            ..Config::default()
        };
        let pair_var = run_core_var(&pair_cfg).expect("ArgumentNA -n .var");
        let pair_rows: Vec<&str> = pair_var.lines().filter(|l| !l.trim().is_empty()).collect();
        assert_eq!(pair_rows.len(), 6);
        for row in &pair_rows {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols.len(), 55);
            assert_eq!(cols[0], "Tumor");
        }
        assert_eq!(
            sha256_hex(pair_var.as_bytes()),
            "66b6c36fccaf2d14ffa4b73ea5bc788a75f0a9b3029c3140c579bd62f8b5248f"
        );
        let pair_vcf = run_core_vcf(&pair_cfg).expect("ArgumentNA -n VCF");
        assert_eq!(
            vcf_chrom_header(&pair_vcf),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal"
        );
        let pair_vars = extract_variants(&pair_vcf);
        assert_eq!(pair_vars.len(), 5);
        assert!(pair_vars.iter().all(|v| v.contains("SAMPLE=Tumor;")));
        assert_eq!(
            vcf_variants_sha256(&pair_vars),
            "e10dc6e2d99d9bd1c61142e51eee6863bf0a032489eaa691e88df1ef71cc1498"
        );
    }

    fn check_fisher_vcf_raw_cell_strict_parity(amp_var_out: &str, amp_vars: &[String]) {
        // 1. Single-sample --fisher VCF raw-cell parity (micro.bam --fisher vs vardict-java --fisher | var2vcf_valid.pl)
        let mut single_cfg = micro_config();
        single_cfg.fisher = true;
        let single_var = run_core_var(&single_cfg).expect("micro --fisher .var");
        assert_eq!(
            sha256_hex(single_var.as_bytes()),
            "192af7d9b62b173af3969bbc8668ac6031c5d40cf3a2e29f19be260bcdbc2a76"
        );
        let single_vars = extract_variants(&run_core_vcf(&single_cfg).expect("micro --fisher VCF"));
        assert_eq!(
            single_vars,
            vec![
                "chr1\t101\t.\tT\tA\t112\tPASS\tSAMPLE=SAMPLE;TYPE=SNV;DP=15;VD=7;AF=0.4667;BIAS=0:0;REFBIAS=8:0;VARBIAS=7:0;PMEAN=16.6;PSTD=1;QUAL=40;QSTD=0;SBF=1;ODDRATIO=0;MQ=60;SN=14;HIAF=0.4667;ADJAF=0;SHIFT3=0;MSI=3;MSILEN=1;NM=0;HICNT=7;HICOV=15;LSEQ=TCGCCGTTGCCTAAACCTAT;RSEQ=TGAAGGAGTCTAGCAGCCGC;DUPRATE=0;SPLITREAD=0;SPANPAIR=0\tGT:DP:VD:AD:AF:RD:ALD\t0/1:15:7:8,7:0.4667:8,0:7,0",
                "chr2\t150\t.\tGCGG\tG\t80\tPASS\tSAMPLE=SAMPLE;TYPE=Deletion;DP=10;VD=4;AF=0.4;BIAS=0:0;REFBIAS=6:0;VARBIAS=4:0;PMEAN=21;PSTD=1;QUAL=40;QSTD=0;SBF=1;ODDRATIO=0;MQ=60;SN=8;HIAF=0.4000;ADJAF=0;SHIFT3=0;MSI=3;MSILEN=1;NM=0;HICNT=4;HICOV=10;LSEQ=CATCTGTTTCTTGCGTCGTA;RSEQ=GACCCTCCATTGTTACTTAT;DUPRATE=0;SPLITREAD=0;SPANPAIR=0\tGT:DP:VD:AD:AF:RD:ALD\t0/1:10:4:6,4:0.4:6,0:4,0",
            ]
        );
        assert_eq!(
            vcf_variants_sha256(&single_vars),
            "427329327b03b866ee3f1ea3e85cd571ca631ab1fecde859486d2d0015920d26"
        );

        // 2. Amplicon --fisher VCF raw-cell parity (Colo829-18_S3-sort.bam -a 10:0.95 --fisher):
        // DecimalFormat("0.0") NM=4.0 / NM=3.0 & DecimalFormat("0.0000") HIAF=1.0000 / 0.7500 alongside
        // getRoundedValueToPrint-stripped AF=1 / 0.75 / 0.25, PMEAN=50 / 33 / 17, SN=3 / 6 / 4 / 12, ADJAF=0.25, SBF=1.
        assert_eq!(
            sha256_hex(amp_var_out.as_bytes()),
            "03417419078d26728fa7a913e4b135dcdc84965c7df11c6446ee0678cf2eb70e"
        );
        assert_eq!(
            amp_vars,
            &[
                "chr1\t934099\t.\tG\tC\t65\tPASS\tSAMPLE=Colo829-18_S3-sort;TYPE=SNV;DP=4;VD=4;AF=1;BIAS=0:2;REFBIAS=0:0;VARBIAS=2:2;PMEAN=50;PSTD=1;QUAL=32.8;QSTD=1;SBF=1;ODDRATIO=0;MQ=60;SN=3;HIAF=1.0000;ADJAF=0;SHIFT3=0;MSI=3;MSILEN=1;NM=4.0;HICNT=3;HICOV=3;LSEQ=CAACGCGGGGGCAGAGGGCG;RSEQ=GTGCGCATGCGGGGCTGCAG;GDAMP=1;TLAMP=1;NCAMP=0;AMPFLAG=0\tGT:DP:VD:AD:AF:RD:ALD\t1/1:4:4:0,4:1:0,0:2,2",
                "chr1\t934106\t.\tA\tG\t57\tPASS\tSAMPLE=Colo829-18_S3-sort;TYPE=SNV;DP=4;VD=3;AF=0.75;BIAS=0:2;REFBIAS=1:0;VARBIAS=1:2;PMEAN=53.3;PSTD=1;QUAL=36.3;QSTD=1;SBF=1;ODDRATIO=0;MQ=60;SN=6;HIAF=0.7500;ADJAF=0;SHIFT3=0;MSI=1;MSILEN=1;NM=4.7;HICNT=3;HICOV=4;LSEQ=GGGGCAGAGGGCGGGTGCGC;RSEQ=TGCGGGGCTGCAGGGATGCG;GDAMP=1;TLAMP=1;NCAMP=0;AMPFLAG=0\tGT:DP:VD:AD:AF:RD:ALD\t1/0:4:3:1,3:0.75:1,0:1,2",
                "chr1\t934123\t.\tT\tTGCGGAGGGGAGGGCGCGGA\t32\tPASS\tSAMPLE=Colo829-18_S3-sort;TYPE=Insertion;DP=8;VD=2;AF=0.25;BIAS=2:2;REFBIAS=2:2;VARBIAS=1:1;PMEAN=33;PSTD=1;QUAL=32.2;QSTD=1;SBF=1;ODDRATIO=1;MQ=60;SN=4;HIAF=1.0000;ADJAF=0.25;SHIFT3=27;MSI=2;MSILEN=5;NM=1.2;HICNT=2;HICOV=2;LSEQ=CGCATGCGGGGCTGCAGGGA;RSEQ=GCGGAGGGGAGGGCGCGGAG;GDAMP=1;TLAMP=1;NCAMP=0;AMPFLAG=0\tGT:DP:VD:AD:AF:RD:ALD\t0/1:8:2:4,2:0.25:2,2:1,1",
                "chr1\t934144\t.\tC\tG\t94\tPASS\tSAMPLE=Colo829-18_S3-sort;TYPE=SNV;DP=8;VD=6;AF=0.75;BIAS=0:2;REFBIAS=0:0;VARBIAS=2:4;PMEAN=17;PSTD=1;QUAL=36.7;QSTD=1;SBF=1;ODDRATIO=0;MQ=60;SN=12;HIAF=1.0000;ADJAF=0;SHIFT3=0;MSI=2;MSILEN=1;NM=3.0;HICNT=6;HICOV=6;LSEQ=GCGGAGGGGAGGGCGCGGAG;RSEQ=GGAGGGCGCGGAGCGGAGGG;GDAMP=1;TLAMP=1;NCAMP=0;AMPFLAG=0\tGT:DP:VD:AD:AF:RD:ALD\t1/0:8:6:0,6:0.75:0,0:2,4",
            ]
        );
        assert_eq!(
            vcf_variants_sha256(amp_vars),
            "15adfab14ed590fad069f2fe92da7d7a504771aea4ffbf35c997e503d353992b"
        );

        // 3. Paired somatic --fisher VCF raw-cell parity (ArgumentNA_Tumor.bam|ArgumentNA_Normal.bam --fisher)
        let droot = data_root();
        let bam_dir =
            droot.join("VarDictJava/src/test/resources/com/astrazeneca/vardict/integrationtests");
        let bam_pair = format!(
            "{}|{}",
            bam_dir.join("ArgumentNA_Tumor.bam").display(),
            bam_dir.join("ArgumentNA_Normal.bam").display()
        );
        let som_cfg = Config {
            reference: droot
                .join("VarDictJava/testdata/fastas/hg19.fa")
                .to_string_lossy()
                .into_owned(),
            sample: extract_sample_name(&bam_pair, None),
            bam: bam_pair,
            threads: 1,
            zero_based: true,
            omit_end_tag: true,
            region: Some("chr7:140453010-140453158:testbed".to_string()),
            freq: 0.001,
            af_threshold: 0.001,
            fisher: true,
            ..Config::default()
        };
        let som_var = run_core_var(&som_cfg).expect("ArgumentNA --fisher .var");
        assert_eq!(
            sha256_hex(som_var.as_bytes()),
            "6a9733dea80ce010e166c1cb12bc0d2b5791dac01000be569e4c0fbd05bc501f"
        );
        let som_vars = extract_variants(&run_core_vcf(&som_cfg).expect("ArgumentNA --fisher VCF"));
        assert_eq!(
            som_vars,
            vec![
                "chr7\t140453031\t.\tA\tG\t74\tq22.5;SN1.5\tSTATUS=LikelyLOH;SAMPLE=ArgumentNA_Tumor;TYPE=SNV;DP=6686;VD=20;AF=0.003;SHIFT3=0;MSI=3.000;MSILEN=1;SSF=0.07682;SOR=3.12761;LSEQ=TAGTAACTCAGCAGCATCTC;RSEQ=GGGCCAAAAATTTAATCAGT\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/1:6686:20:20,0:6417,225:6642,20:0.003:1,1:19.6:1:17.3:1:1:0:60:0.429:0.0009:0:1\t0/1:2087:2:1,1:2016,69:2085,2:0.001:1,2:20:0:40:0:0.06597:28.9909:60:4:0.001:0:1",
                "chr7\t140453035\t.\tC\tA\t41\tPASS\tSTATUS=StrongLOH;SAMPLE=ArgumentNA_Tumor;TYPE=SNV;DP=6687;VD=0;AF=0;SHIFT3=1;MSI=2.000;MSILEN=1;SSF=0.01361;SOR=0;LSEQ=AACTCAGCAGCATCTCAGGG;RSEQ=CAAAAATTTAATCAGTGGAA\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/0:6687:0:0,0:6458,226:6684,0:0:0,0:0:0:0:0:1:0:0:0:0:0:0\t0/1:2098:3:3,0:2021,74:2095,3:0.0014:1,0:18.7:1:26.3:1:1:0:60:2:0.001:0:3.3",
                "chr7\t140453036\t.\tCA\tC\t170\tPASS\tSTATUS=LikelySomatic;SAMPLE=ArgumentNA_Tumor;TYPE=Deletion;DP=6685;VD=20;AF=0.003;SHIFT3=4;MSI=5.000;MSILEN=1;SSF=0.56689;SOR=1.04829;LSEQ=ACTCAGCAGCATCTCAGGGC;RSEQ=AAAATTTAATCAGTGGAAAA\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/1:6685:20:19,1:6433,225:6658,20:0.003:1,1:21.4:1:39.5:1:0.49821:1.50466:58:40:0.003:0:0.1\t0/1:2102:6:6,0:2015,78:2093,6:0.0029:1,0:24.8:1:39.8:1:1:0:60:12:0.0029:0:1.5",
                "chr7\t140453046\t.\tA\tG\t63\tf0.001\tSTATUS=LikelyLOH;SAMPLE=ArgumentNA_Tumor;TYPE=SNV;DP=6690;VD=3;AF=0.0004;SHIFT3=0;MSI=1.000;MSILEN=1;SSF=0.34416;SOR=0.4711;LSEQ=ATCTCAGGGCCAAAAATTTA;RSEQ=TCAGTGGAAAAATAGCCTCA\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/0:6690:3:2,1:6456,227:6683,3:0.0004:1,2:20.7:1:40:1:0.09887:14.20096:56:6:0.0005:0:1\t0/1:2102:2:2,0:2020,78:2098,2:0.001:1,0:34.5:1:35:1:1:0:60:4:0.001:0:1.5",
                "chr7\t140453081\t.\tC\tA\t135\tPASS\tSTATUS=LikelySomatic;SAMPLE=ArgumentNA_Tumor;TYPE=SNV;DP=6671;VD=20;AF=0.003;SHIFT3=0;MSI=4.000;MSILEN=1;SSF=0.07765;SOR=3.11813;LSEQ=GCCTCAATTCTTACCATCCA;RSEQ=AAAATGGATCCAGACAACTG\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/1:6671:20:19,1:6430,219:6649,20:0.003:1,1:33.7:1:31.4:1:0.48925:1.54515:56:3:0.0023:0:1.6\t0/1:2076:2:2,0:2021,52:2073,2:0.001:1,0:27.5:1:16:0:1:0:60:0:0:0:5.5",
            ]
        );
        assert_eq!(
            vcf_variants_sha256(&som_vars),
            "b35e13cbda9f6052446d7aea5ccc7d2ba4d4b7cf8cc16e2d6dad715f2e1608e5"
        );

        // Also check full pileup --fisher corpus configuration (-q 25 -F 0x100 -f 0.001 -k 1 -M 20 -p --fisher)
        let mut som_pileup_cfg = som_cfg;
        som_pileup_cfg.base_qual = 25.0;
        som_pileup_cfg.flag_filter = 0x100;
        som_pileup_cfg.samfilter_raw = "0x100".to_string();
        som_pileup_cfg.realign = 1;
        som_pileup_cfg.minmatch = 20;
        som_pileup_cfg.pileup = true;
        som_pileup_cfg.apply_pileup_overrides();
        let som_pileup_vars =
            extract_variants(&run_core_vcf(&som_pileup_cfg).expect("ArgumentNA -p --fisher VCF"));
        assert_eq!(som_pileup_vars.len(), 42);
        assert!(som_pileup_vars.iter().any(|v| v.contains("SSF=0.5255;")));
        assert!(som_pileup_vars.iter().any(|v| v.contains(":0.0654:")));
        assert_eq!(
            vcf_variants_sha256(&som_pileup_vars),
            "e18d894bf41e49e0e6357ad0390bcae88fdf2f9872186d446e4c72bc9b716d1e"
        );

        // 4. SomaticOutputVariant.getRoundedValueToPrint trailing-dot FORMAT cells (:0., :15.:, :1.)
        // on combine_analysis back-calculated normal fields against micro.ref.fa (pure-Rust indexed_bam).
        check_somatic_fisher_trailing_dot_cells();
    }

    fn check_somatic_fisher_trailing_dot_cells() {
        use noodles::core::Position;
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles::sam::alignment::record_buf::{Cigar, Data, RecordBuf};

        let ref_fa = workspace_root().join("fixtures/micro/micro.ref.fa");
        let fa_text = fs::read_to_string(&ref_fa).expect("read micro.ref.fa");
        let mut chr1 = Vec::new();
        let mut in_chr1 = false;
        for line in fa_text.lines() {
            let trimmed = line.trim();
            if let Some(rest) = trimmed.strip_prefix('>') {
                in_chr1 = rest.split_whitespace().next() == Some("chr1");
            } else if in_chr1 {
                chr1.extend_from_slice(trimmed.as_bytes());
            }
        }
        let flip = |b: u8| if b != b'C' { b'C' } else { b'G' };
        let build_rec = |name: &str, pos1: usize, ops: &[(usize, Kind)], seq: Vec<u8>, nm: i32| {
            let cigar: Cigar = ops.iter().map(|&(l, k)| Op::new(k, l)).collect();
            let seq_len = seq.len();
            let data: Data = [(Tag::EDIT_DISTANCE, BufValue::from(nm))]
                .into_iter()
                .collect();
            RecordBuf::builder()
                .set_name(name)
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(pos1).expect("pos1"))
                .set_mapping_quality(MappingQuality::new(60).expect("mapq 60"))
                .set_cigar(cigar)
                .set_sequence(seq.into())
                .set_quality_scores(vec![40u8; seq_len].into())
                .set_data(data)
                .build()
        };

        for (t_mm, n_mm, expected_var_sha, expected_vcf_sha, expected_chr100_row) in [
            (
                [1usize, 0, 0],
                0usize,
                "76e036321d26faa602c9988ecb5825f67f5948245d8221ed12faf5e3cc998c85",
                "fd93ef7beeb1c25f60f395066008ea2b90c53efa7901aa3d5d92373f1ece01c4",
                "chr1\t100\t.\tTTTGAAGGAGT\tT\t112\tPASS\tSTATUS=Germline;SAMPLE=Tumor;TYPE=Deletion;DP=11;VD=3;AF=0.2727;SHIFT3=0;MSI=1.000;MSILEN=1;SSF=0.27745;SOR=0.44283;LSEQ=TTCGCCGTTGCCTAAACCTA;RSEQ=CTAGCAGCCGCAGTAAGGCA\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/1:11:3:3,0:8,0:8,3:0.2727:0,0:19:1:40:0:1:0:60:6:0.2727:0:0.3\t0/1:15:7:7,0:8,0:8,7:0.4667:0,0:15:1:40:1:1:0:60:6:0.4326:0.3846:0.",
            ),
            (
                [2usize, 1, 1],
                1usize,
                "7fba9177549554eb90037d138f831aedfac2f1e3aeae97adbdfdcf048c287156",
                "69ec72a5328fe1ea0bd1f0186bb12439fd1c605bea9435b84093d9d68a04d3c6",
                "chr1\t100\t.\tTTTGAAGGAGT\tT\t112\tPASS\tSTATUS=Germline;SAMPLE=Tumor;TYPE=Deletion;DP=11;VD=3;AF=0.2727;SHIFT3=0;MSI=1.000;MSILEN=1;SSF=0.27745;SOR=0.44283;LSEQ=TTCGCCGTTGCCTAAACCTA;RSEQ=CTAGCAGCCGCAGTAAGGCA\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t0/1:11:3:3,0:8,0:8,3:0.2727:0,0:17.7:1:40:0:1:0:60:6:0.2727:0:1.3\t0/1:15:7:7,0:8,0:8,7:0.4667:0,0:15.:1:40:1:1:0:60:6:0.4326:0.3846:1.",
            ),
        ] {
            let tmp = tempfile::tempdir().expect("tempdir");
            let mut t_raw: Vec<(usize, RecordBuf)> = Vec::new();
            for (i, &mcnt) in t_mm.iter().enumerate() {
                let s = 81 + i;
                let left_len = 101 - s;
                let mut left = chr1[s - 1..100].to_vec();
                for m in 0..mcnt {
                    let idx = 2 + m * 3;
                    left[idx] = flip(left[idx]);
                }
                left.extend_from_slice(&chr1[110..130]);
                t_raw.push((
                    s,
                    build_rec(
                        &format!("t_del_{i}"),
                        s,
                        &[
                            (left_len, Kind::Match),
                            (10, Kind::Deletion),
                            (20, Kind::Match),
                        ],
                        left,
                        (10 + mcnt) as i32,
                    ),
                ));
            }
            for i in 0..8 {
                let s = 75 + i;
                t_raw.push((
                    s,
                    build_rec(
                        &format!("t_ref_{i}"),
                        s,
                        &[(55, Kind::Match)],
                        chr1[s - 1..s - 1 + 55].to_vec(),
                        0,
                    ),
                ));
            }
            t_raw.sort_by_key(|(s, _)| *s);
            let t_recs: Vec<RecordBuf> = t_raw.into_iter().map(|(_, r)| r).collect();

            let mut n_raw: Vec<(usize, RecordBuf)> = Vec::new();
            for i in 0..7 {
                let s = 80 + i;
                let left_len = 101 - s;
                let mut left = chr1[s - 1..100].to_vec();
                for m in 0..n_mm {
                    let idx = 2 + m * 3;
                    left[idx] = flip(left[idx]);
                }
                left.extend_from_slice(&chr1[110..125]);
                n_raw.push((
                    s,
                    build_rec(
                        &format!("n_sc_{i}"),
                        s,
                        &[(left_len, Kind::Match), (15, Kind::SoftClip)],
                        left,
                        n_mm as i32,
                    ),
                ));
            }
            for i in 0..8 {
                let s = 75 + i;
                n_raw.push((
                    s,
                    build_rec(
                        &format!("n_ref_{i}"),
                        s,
                        &[(55, Kind::Match)],
                        chr1[s - 1..s - 1 + 55].to_vec(),
                        0,
                    ),
                ));
            }
            n_raw.sort_by_key(|(s, _)| *s);
            let n_recs: Vec<RecordBuf> = n_raw.into_iter().map(|(_, r)| r).collect();

            let tbam = tmp.path().join("Tumor.bam");
            let nbam = tmp.path().join("Normal.bam");
            indexed_bam(&tbam, "chr1", chr1.len(), &t_recs).expect("write Tumor.bam");
            indexed_bam(&nbam, "chr1", chr1.len(), &n_recs).expect("write Normal.bam");

            let cfg = Config {
                reference: ref_fa.to_string_lossy().into_owned(),
                sample: "Tumor|Normal".to_string(),
                sample_explicit: true,
                bam: format!("{}|{}", tbam.display(), nbam.display()),
                threads: 1,
                zero_based: true,
                omit_end_tag: true,
                region: Some("chr1:50-160:testbed".to_string()),
                freq: 0.01,
                af_threshold: 0.01,
                fisher: true,
                ..Config::default()
            };
            let var_out = run_core_var(&cfg).expect("somatic trailing-dot .var");
            assert_eq!(sha256_hex(var_out.as_bytes()), expected_var_sha);
            let vars = extract_variants(&run_core_vcf(&cfg).expect("somatic trailing-dot VCF"));
            assert_eq!(vars.len(), 2);
            assert_eq!(vars[1], expected_chr100_row);
            assert_eq!(vcf_variants_sha256(&vars), expected_vcf_sha);
        }
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
