### Table 8. Supply-Chain Security, Static musl Binary Equivalence, Functional Coverage & Divergence Summary

| Engineering & Security Dimension | VarDictJava v1.8.3 (0886261) | vardictcpp (ae93fe39dc) | VarDictRust v1.0.6 (5865a2a0b24a8c17) |
| :--- | :--- | :--- | :--- |
| Upstream Repository Maintenance Status | Archived Read-Only (EOL) | Preprint Prototype (10 Divergences) | Actively Maintained (MIT, v1.0.6) |
| Runtime Dependency & Binary Footprint | JVM + R + Perl (102 conda packages) | C++17 + htslib dynamic libs | 4.46 MiB GNU / 4.56 MiB Static musl (statically linked) |
| Supply-Chain Attestation & Vulnerability Audit | Unpinned ecosystem dependencies | No SLSA provenance | SLSA Level 3 + 0 RustSec CVEs |
| 13-Mode Functional Coverage & CRAM Support | 12/13 (CRAM fails without htsjdk classpath) | 12/13 (Amplicon VCF & Splice -i unsupported) | 13/13 Modes (100% Exact Match) |
| Empirical Divergence Classes vs. Oracle | Reference Oracle (Baseline) | 10 Cataloged Divergence Classes (+309,569 chunk-boundary duplicates on unchunked HG002 WGS) | 0 Divergences on Valid Inputs (#![forbid(unsafe_code)]) |

*Static musl binary equivalence verified across all regression and multimodal benchmark suites.*
