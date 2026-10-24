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

    // Minor 9: Use a regular file as parent so create_dir_all fails with ENOTDIR even when run as root.
    let bad_parent = tempfile::NamedTempFile::new().unwrap();
    let bad_dir = bad_parent.path().join("sub");
    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_CPU", &bad_dir)
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
    let stderr = String::from_utf8_lossy(&bad.stderr);
    assert!(
        stderr.contains("VARDICT_PROFILE_CPU") || stderr.contains(bad_dir.to_str().unwrap()),
        "stderr must name the variable or path, got: {stderr}"
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

    let bad_parent = tempfile::NamedTempFile::new().unwrap();
    let bad_dir = bad_parent.path().join("sub");
    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_DHAT", &bad_dir)
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
    let stderr = String::from_utf8_lossy(&bad.stderr);
    assert!(
        stderr.contains("VARDICT_PROFILE_DHAT") || stderr.contains(bad_dir.to_str().unwrap()),
        "stderr must name the variable or path, got: {stderr}"
    );
}

// Minor 5: When profile-dhat is enabled, dhat owns the global allocator, disabling jemalloc heap profiling.
#[cfg(all(feature = "profile-heap", not(feature = "profile-dhat")))]
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

    let bad_parent = tempfile::NamedTempFile::new().unwrap();
    let bad_dir = bad_parent.path().join("sub");
    let bad = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_HEAP", &bad_dir)
        .args(["--version"])
        .output()
        .unwrap();
    assert!(
        !bad.status.success(),
        "unwritable profile dir must fail before the run"
    );
    let stderr = String::from_utf8_lossy(&bad.stderr);
    assert!(
        stderr.contains("VARDICT_PROFILE_HEAP") || stderr.contains(bad_dir.to_str().unwrap()),
        "stderr must name the variable or path, got: {stderr}"
    );
}

#[test]
fn compiled_out_warning_when_features_disabled_and_no_files_created() {
    let _dir = tempfile::tempdir().unwrap();

    #[cfg(not(feature = "profile-cpu"))]
    {
        let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
            .current_dir(root())
            .env("VARDICT_PROFILE_CPU", _dir.path())
            .args(["--version"])
            .output()
            .unwrap();
        assert!(out.status.success());
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("VARDICT_PROFILE_CPU is set, but the feature is not compiled in"),
            "expected compiled-out CPU warning, got: {stderr}"
        );
        assert_eq!(
            std::fs::read_dir(_dir.path()).unwrap().count(),
            0,
            "no files must be created when feature is compiled out"
        );
    }

    #[cfg(not(feature = "profile-dhat"))]
    {
        let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
            .current_dir(root())
            .env("VARDICT_PROFILE_DHAT", _dir.path())
            .args(["--version"])
            .output()
            .unwrap();
        assert!(out.status.success());
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("VARDICT_PROFILE_DHAT is set, but the feature is not compiled in"),
            "expected compiled-out DHAT warning, got: {stderr}"
        );
        assert_eq!(
            std::fs::read_dir(_dir.path()).unwrap().count(),
            0,
            "no files must be created when feature is compiled out"
        );
    }

    #[cfg(not(feature = "profile-heap"))]
    {
        let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
            .current_dir(root())
            .env("VARDICT_PROFILE_HEAP", _dir.path())
            .args(["--version"])
            .output()
            .unwrap();
        assert!(out.status.success());
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("VARDICT_PROFILE_HEAP is set, but the feature is not compiled in"),
            "expected compiled-out HEAP warning, got: {stderr}"
        );
        assert_eq!(
            std::fs::read_dir(_dir.path()).unwrap().count(),
            0,
            "no files must be created when feature is compiled out"
        );
    }
}

#[test]
fn version_and_help_with_no_env_vars_creates_no_files() {
    let dir = tempfile::tempdir().unwrap();
    let version = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(dir.path())
        .env_remove("VARDICT_PROFILE_CPU")
        .env_remove("VARDICT_PROFILE_HEAP")
        .env_remove("VARDICT_PROFILE_DHAT")
        .args(["--version"])
        .output()
        .unwrap();
    assert!(version.status.success());
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        0,
        "running --version with no env vars must create no files"
    );

    let help = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(dir.path())
        .env_remove("VARDICT_PROFILE_CPU")
        .env_remove("VARDICT_PROFILE_HEAP")
        .env_remove("VARDICT_PROFILE_DHAT")
        .args(["--help"])
        .output()
        .unwrap();
    assert!(help.status.success());
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        0,
        "running --help with no env vars must create no files"
    );
}

#[test]
fn empty_profile_env_var_fails_fast_with_actionable_error() {
    for var in [
        "VARDICT_PROFILE_CPU",
        "VARDICT_PROFILE_HEAP",
        "VARDICT_PROFILE_DHAT",
    ] {
        let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
            .current_dir(root())
            .env(var, "")
            .args(["--version"])
            .output()
            .unwrap();
        assert!(
            !out.status.success(),
            "empty {var} must fail fast before running"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(&format!("{var} must be a non-empty directory path")),
            "expected actionable error for empty {var}, got: {stderr}"
        );
    }
}

#[cfg(all(feature = "profile-heap", feature = "profile-dhat"))]
#[test]
fn heap_dhat_allocator_mismatch_warning() {
    let dir = tempfile::tempdir().unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .current_dir(root())
        .env("VARDICT_PROFILE_HEAP", dir.path())
        .args(["--version"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("VARDICT_PROFILE_HEAP requested, but global allocator does not match; disabling heap profiling"),
        "expected allocator mismatch warning, got: {stderr}"
    );
}
