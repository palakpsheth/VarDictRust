#[test]
fn version_credits_original_and_discloses_ai() {
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_vardictrust"))
        .arg("--version")
        .output()
        .expect("run vardictrust");
    assert!(out.status.success());
    let s = String::from_utf8_lossy(&out.stdout);
    for needle in [
        concat!("vardictrust ", env!("CARGO_PKG_VERSION")),
        "VarDictJava v1.8.3",
        "0886261",
        "009e017",
        "Lai",
        "gkw227",
        "AstraZeneca-NGS",
        "AI disclosure",
        "MIT",
    ] {
        assert!(s.contains(needle), "--version lacks {needle:?}:\n{s}");
    }
}
