use std::path::Path;
use std::process::Command;

fn workspace_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
}

#[test]
fn test_dp_err_routes_to_stderr() {
    let bin = env!("CARGO_BIN_EXE_vardictrust");
    let out = Command::new(bin)
        .current_dir(workspace_root())
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
            "fixtures/micro/micro.bed",
            "--DP",
            "ERR",
        ])
        .output()
        .expect("run vardictrust with --DP ERR");

    if !out.status.success() {
        panic!(
            "command failed: {}\nstdout: {}\nstderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    assert!(out.status.success());
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty under --DP ERR, but got {} bytes",
        out.stdout.len()
    );
    assert!(
        !out.stderr.is_empty(),
        "stderr must contain variant output under --DP ERR, but was empty"
    );
}

#[test]
fn test_default_printer_long_alias_routes_to_stderr() {
    let bin = env!("CARGO_BIN_EXE_vardictrust");
    let out = Command::new(bin)
        .current_dir(workspace_root())
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
            "fixtures/micro/micro.bed",
            "--defaultPrinter",
            "ERR",
        ])
        .output()
        .expect("run vardictrust with --defaultPrinter ERR");

    if !out.status.success() {
        panic!(
            "command failed: {}\nstdout: {}\nstderr: {}",
            out.status,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    assert!(out.status.success());
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty under --defaultPrinter ERR, but got {} bytes",
        out.stdout.len()
    );
    assert!(
        !out.stderr.is_empty(),
        "stderr must contain variant output under --defaultPrinter ERR, but was empty"
    );
}
