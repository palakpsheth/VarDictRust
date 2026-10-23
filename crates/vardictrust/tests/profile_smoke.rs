#![cfg(any(
    feature = "profile-cpu",
    feature = "profile-dhat",
    feature = "profile-heap"
))]

use std::process::Command;

fn root() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[cfg(feature = "profile-cpu")]
#[test]
fn cpu_profile_files_exist_and_dir_probe_fails_fast() {
    let dir = tempfile::tempdir().unwrap();
    let ok = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_CPU", dir.path())
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-z",
            "-th",
            "2",
            "fixtures/micro/micro.bed",
        ])
        .output()
        .unwrap();
    assert!(ok.status.success());
    // ms-long run => 0 samples is legal at 499 Hz: assert existence/validity, never non-emptiness.
    assert!(dir.path().join("cpu-collapsed.txt").exists());
    assert!(dir.path().join("cpu-profile.pb").exists());
    assert!(dir.path().join("cpu-meta.json").exists());
    assert!(dir.path().join("cpu-flamegraph.svg").exists());

    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_CPU", "/nonexistent/readonly/dir")
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
}

#[cfg(feature = "profile-dhat")]
#[test]
fn dhat_profile_files_exist_and_dir_probe_fails_fast() {
    let dir = tempfile::tempdir().unwrap();
    let ok = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_DHAT", dir.path())
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-z",
            "-th",
            "2",
            "fixtures/micro/micro.bed",
        ])
        .output()
        .unwrap();
    assert!(ok.status.success());
    assert!(dir.path().join("dhat-heap.json").exists());

    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_DHAT", "/nonexistent/readonly/dir")
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
}

#[cfg(feature = "profile-heap")]
#[test]
fn heap_profile_files_exist_and_dir_probe_fails_fast() {
    let dir = tempfile::tempdir().unwrap();
    let ok = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_HEAP", dir.path())
        .args([
            "-G",
            "fixtures/micro/micro.ref.fa",
            "-b",
            "fixtures/micro/micro.bam",
            "-c",
            "1",
            "-S",
            "2",
            "-E",
            "3",
            "-g",
            "4",
            "-z",
            "-th",
            "2",
            "fixtures/micro/micro.bed",
        ])
        .output()
        .unwrap();
    assert!(ok.status.success());
    let gz_exists = dir.path().join("heap.pb.gz").exists();
    let pb_exists = dir.path().join("heap.pb").exists();
    assert!(gz_exists || pb_exists, "heap protobuf file must exist");
    assert!(dir.path().join("heap-flamegraph.svg").exists());

    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_HEAP", "/nonexistent/readonly/dir")
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
}
