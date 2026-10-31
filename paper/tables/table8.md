### Table 8. Supply-Chain Security (H4), Static musl Binary Equivalence, Functional Coverage & Divergence Summary

| Source Cell ID | Dimension / Criterion | VarDictJava v1.8.3 (0886261) | vardictcpp (ae93fe39dc) | VarDictRust v1.0.6 (5865a2a0b24a8c17) |
| :--- | :--- | :--- | :--- | :--- |
| s5_h4_maintenance_security_matrix | Upstream Repository Maintenance Status | Archived Read-Only (EOL) | Preprint Prototype (10 Divergences) | Actively Maintained (MIT, v1.0.6) |
| s5_h4_maintenance_security_matrix | Runtime Dependency & Binary Footprint | JVM + R + Perl (102 conda packages) | C++17 + htslib dynamic libs | 4.46 MiB GNU / 4.56 MiB Static musl (statically linked) |
| s5_h4_maintenance_security_matrix | Supply-Chain Attestation & Vulnerability Audit | Unpinned ecosystem dependencies | No SLSA provenance | SLSA Level 3 + 0 RustSec CVEs |
| s2_tierp_mode_parity_sweep | 13-Mode Functional Coverage & CRAM Support | 12/13 (CRAM fails without htsjdk classpath) | 12/13 (Amplicon CPP-3 & Splice -i unsupported) | 13/13 Modes (100% MATCH) |
| s5_h4_maintenance_security_matrix | Empirical Divergence Classes vs Oracle | Reference Oracle (Baseline) | CPP-1..10 (+309,569 boundary dupes on HG002 Arm A) | 0 Divergences on Valid Inputs (#![forbid(unsafe_code)]) |

*Static musl binary equivalence verified via harness/bin_equiv.sh and certified in evidence_certified/musl-ci.json.*
