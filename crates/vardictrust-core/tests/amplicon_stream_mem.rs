#![cfg(target_os = "linux")]
//! -p streaming memory pin (audit-amplicon-expand follow-up).
//!
//! The -p amplicon path emits a (null/ref) row at EVERY insert base, so the
//! old per-base `BTreeMap` materialized one entry per base before any row was
//! written (measured: 902 MB VmHWM for a 2 Mbp insert in a debug build). The
//! streaming path must emit all rows with flat working memory.
//!
//! Process-wide VmHWM cannot pin this in-suite: parallel test threads pollute
//! it (measured 315 MB of unrelated-test noise in one full-suite run vs
//! <64 MB isolated for this exact stream). Following tests/cram_stack.rs,
//! the stream therefore runs in an isolated child process (re-executed test
//! binary, `VDR_STREAM_CHILD=1`), which self-reports rows, wall time and its
//! own VmHWM; the parent asserts all three.

use std::collections::HashSet;
use std::process::Command;

use vardictrust_core::config::Config;
use vardictrust_core::postprocess::amplicon::{
    amplicon_position_stream, amplicon_post_process_streaming, AmpliconVariant,
};
use vardictrust_core::regions::{AmpliconInfo, Region};
use vardictrust_core::types::{Pos1, VardictError};

const INSERT_LEN: u32 = 500_000;
// Old per-base BTreeMap uses ~225 MB VmHWM for a 500 kbp insert (debug build);
// the streaming path uses ~20 MB and must stay well below 64 MB.
const PEAK_KB_CAP: u64 = 64_000;
const WALL_SECS_CAP: u64 = 60;

fn peak_kb() -> Option<u64> {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines().find_map(|l| {
                l.strip_prefix("VmHWM:")?
                    .split_whitespace()
                    .next()?
                    .parse()
                    .ok()
            })
        })
}

/// Child side: stream the -p null rows for the whole synthetic insert and
/// self-report `ROWS=<n> ELAPSED_MS=<ms> PEAK_KB=<hwm>`.
#[test]
fn pileup_stream_child() {
    if std::env::var("VDR_STREAM_CHILD").is_err() {
        return;
    }
    let config = Config {
        pileup: true,
        emit_intermediate: true,
        ..Default::default()
    };
    let region = Region {
        chr: "chr1".to_string(),
        start: Pos1(1),
        end: Pos1(INSERT_LEN),
        gene: "TEST".to_string(),
        amplicon: Some(AmpliconInfo {
            insert_start: Pos1(1),
            insert_end: Pos1(INSERT_LEN),
        }),
    };
    let segment = vec![region.clone()];
    let vars_per_amplicon = vec![vardictrust_core::javamap::JavaHashMap::new()];

    let t0 = std::time::Instant::now();
    let mut rows = 0u64;
    {
        let mut count = |_av: AmpliconVariant| -> Result<(), VardictError> {
            rows += 1;
            Ok(())
        };
        amplicon_post_process_streaming(
            &config,
            &region,
            &vars_per_amplicon,
            amplicon_position_stream(&segment),
            &HashSet::new(),
            &mut count,
        )
        .expect("streaming is infallible here");
    }
    let elapsed_ms = t0.elapsed().as_millis();
    let peak = peak_kb().expect("/proc/self/status VmHWM");
    println!("ROWS={rows} ELAPSED_MS={elapsed_ms} PEAK_KB={peak}");
}

fn run_child() -> std::process::Output {
    let exe = std::env::current_exe().expect("current test exe");
    Command::new(exe)
        .arg("pileup_stream_child")
        .arg("--exact")
        .arg("--nocapture")
        .env("VDR_STREAM_CHILD", "1")
        .output()
        .expect("spawn child")
}

fn parse_field(stdout: &str, key: &str) -> u64 {
    stdout
        .split_whitespace()
        .find_map(|tok| {
            tok.strip_prefix(key)
                .and_then(|v| v.strip_prefix('='))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or_else(|| panic!("child stdout missing {key}=...: {stdout}"))
}

#[test]
fn pileup_stream_huge_interval_bounded_child() {
    let child = run_child();
    assert!(
        child.status.success(),
        "stream child failed ({}): {}",
        child.status,
        String::from_utf8_lossy(&child.stderr)
    );
    let stdout = String::from_utf8_lossy(&child.stdout).into_owned();
    assert_eq!(
        parse_field(&stdout, "ROWS"),
        u64::from(INSERT_LEN),
        "every insert base must yield its -p row: {stdout}"
    );
    assert!(
        parse_field(&stdout, "ELAPSED_MS") < WALL_SECS_CAP * 1000,
        "2 Mbp -p stream took too long: {stdout}"
    );
    assert!(
        parse_field(&stdout, "PEAK_KB") < PEAK_KB_CAP,
        "child peak RSS over the 2 Mbp -p stream exceeds the cap (old map: ~902 MB): {stdout}"
    );
}
