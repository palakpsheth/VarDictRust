//! CRAM worker-stack regression pin (T3D-03).
//!
//! Root cause (verified by a symbolized debug backtrace of the release engine
//! on the HG008 CRAM slices): decoding a high-depth CRAM slice nests ~97
//! frames deep inside `noodles-cram 0.100.0` — container/slice iterators,
//! `Block::decode`, `name_tokenizer` (10k-name blocks), striped `rans_nx16`
//! chunk decodes — which overflows the 2 MiB default worker-thread
//! stack (`thread '<unknown>' has overflowed its stack`). No frame in
//! `vardictrust-core` is oversized, so there is nothing to Box on our side;
//! the `RUST_MIN_STACK=67108864` override in `harness/cram_twin.sh` is
//! required until upstream shrinks its decoder frames.
//!
//! This test pins the workaround with measured thresholds: the tile-0 CRAM
//! query MUST die on a small stack and MUST succeed with an identical record
//! count on a 64 MiB stack. Each case runs in a child process so the fatal
//! overflow cannot abort the test harness. If the small-stack child ever
//! starts succeeding, upstream likely fixed the frames and the workaround
//! should be re-evaluated (the test fails loudly in that case by design).
//!
//! The small stack is profile-dependent (measured on the reference host):
//! release needs somewhere in (2 MiB, 8 MiB] — the 2 MiB rayon-worker default
//! overflows deterministically — while debug (unoptimized frames) needs only
//! (256 KiB, 1 MiB]. Both are far below the 64 MiB
//! workaround, so the pin is robust in either profile.
//!
//! Skipped when the scratch HG008 slices are absent (CI has no such data;
//! `harness/cram_twin.sh` covers the full-path gate where scratch exists).
//! Synthetic CRAMs do NOT reproduce this: htslib emits gzip-only blocks for
//! small/synthetic inputs and reserves striped-rANS name blocks for real
//! high-depth data (verified by block-codec profiling during diagnosis).

use std::path::{Path, PathBuf};
use std::process::Command;

use vardictrust_core::config::Config;
use vardictrust_core::types::Pos0;

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../")
}

// Tile 0 of scratch/g4/hg008/engine.bed; measured: 1805 records from
// HG008-T.slice.cram at >=8 MiB stacks, SIGABRT at 2 MiB (release).
const TILE_CHR: &str = "chr1";
const TILE_START: u32 = 1743056;
const TILE_END: u32 = 1744086;
const TILE_T_COUNT: usize = 1805;

fn scratch_present(root: &Path) -> bool {
    root.join("scratch/g4/hg008/HG008-T.slice.cram").exists()
        && root
            .join("scratch/g4/hg008/HG008-T.slice.cram.crai")
            .exists()
        && root.join("scratch/g4/ref/giabv3/GRCh38_GIABv3.fa").exists()
        && root
            .join("scratch/g4/ref/giabv3/GRCh38_GIABv3.fa.fai")
            .exists()
}

fn query_tile0(cram: &Path, reference: &str) -> usize {
    let config = Config {
        reference: reference.to_string(),
        ..Default::default()
    };
    vardictrust_core::bam::query_region(cram, TILE_CHR, Pos0(TILE_START), Pos0(TILE_END), &config)
        .expect("tile-0 CRAM query")
        .len()
}

/// Child side: decode tile 0 on a thread with `VDR_CRam_STACK` bytes of stack.
/// On overflow the process dies by signal (never returns); on success it
/// prints `COUNT=<n>` and exits 0.
#[test]
fn cram_stack_child() {
    if std::env::var("VDR_CRam_CHILD").is_err() {
        return;
    }
    let stack: usize = std::env::var("VDR_CRam_STACK")
        .expect("VDR_CRam_STACK")
        .parse()
        .expect("parse VDR_CRam_STACK");
    let root = workspace_root();
    let cram = root.join("scratch/g4/hg008/HG008-T.slice.cram");
    let reference = root
        .join("scratch/g4/ref/giabv3/GRCh38_GIABv3.fa")
        .to_str()
        .unwrap()
        .to_string();
    let h = std::thread::Builder::new()
        .name("cram-stack-child".to_string())
        .stack_size(stack)
        .spawn(move || query_tile0(&cram, &reference))
        .expect("spawn decode thread");
    let n = h.join().expect("decode thread survived");
    println!("COUNT={n}");
}

fn run_child(stack: usize) -> std::process::Output {
    let exe = std::env::current_exe().expect("current test exe");
    Command::new(exe)
        .arg("cram_stack_child")
        .arg("--exact")
        .arg("--nocapture")
        .env("VDR_CRam_CHILD", "1")
        .env("VDR_CRam_STACK", stack.to_string())
        // The child sets its own explicit stack size; the inherited
        // RUST_MIN_STACK must not interfere either way.
        .env_remove("RUST_MIN_STACK")
        .output()
        .expect("spawn child")
}

#[test]
fn cram_worker_stack_thresholds() {
    let root = workspace_root();
    if !scratch_present(&root) {
        println!("SKIP: scratch HG008 CRAM slices absent (CI has no such data)");
        return;
    }

    // Small stack: the decode MUST die. Release uses the 2 MiB rayon-worker
    // default (measured: overflows); debug uses 256 KiB (measured: overflows;
    // 2 MiB fits unoptimized frames, so it cannot pin there).
    let small: usize = if cfg!(debug_assertions) {
        256 * 1024
    } else {
        2 * 1024 * 1024
    };
    let died = run_child(small);
    assert!(
        !died.status.success(),
        "{small}B child unexpectedly SUCCEEDED ({}). If upstream shrank the \
         noodles-cram decoder frames, re-evaluate the RUST_MIN_STACK=64M \
         workaround in harness/cram_twin.sh — it may now be removable.",
        String::from_utf8_lossy(&died.stdout)
    );
    // Death must be by signal (stack overflow), not a wrong-reason exit
    // (missing .crai, config error, ...): on Unix a signal death has no
    // exit code and the runtime prints "overflowed its stack" to stderr.
    #[cfg(unix)]
    assert!(
        (died.status.code().is_none()
            || String::from_utf8_lossy(&died.stderr).contains("overflowed its stack"))
            && !String::from_utf8_lossy(&died.stderr).contains("tile-0 CRAM query"),
        "{small}B child failed for the wrong reason ({}): stdout={} stderr={}. Expected death-by-signal (stack overflow).",
        died.status,
        String::from_utf8_lossy(&died.stdout),
        String::from_utf8_lossy(&died.stderr)
    );
    println!("PINNED: {small}B worker stack fails ({})", died.status);

    // Workaround stack: the decode MUST succeed with the measured count.
    let big = run_child(64 * 1024 * 1024);
    assert!(
        big.status.success(),
        "64 MiB child failed ({}): {}",
        big.status,
        String::from_utf8_lossy(&big.stderr)
    );
    let stdout = String::from_utf8_lossy(&big.stdout);
    assert!(
        stdout.contains(&format!("COUNT={TILE_T_COUNT}")),
        "64 MiB child decoded an unexpected record count: {stdout}"
    );
    println!("PINNED: 64 MiB worker stack decodes {TILE_T_COUNT} records");
}
