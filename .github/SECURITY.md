# Security Policy: VarDictRust

The VarDictRust project takes the security, memory safety, and integrity of scientific genomics workflows seriously. As an idiomatic Rust rewrite engineered to process large, untrusted Next-Generation Sequencing (NGS) files (BAM, CRAM, FASTA, BED), VarDictRust enforces `#![forbid(unsafe_code)]` in its core variant calling engine (`crates/vardict-core`) to eliminate entire classes of memory corruption vulnerabilities by construction.

---

## Supported Versions

Only the latest active major and minor releases receive security updates and vulnerability patches.

| Version | Supported | Security Maintenance |
| :--- | :--- | :--- |
| **`1.0.x`** | :white_check_mark: | Active support (patches for critical/high vulnerabilities and supply-chain advisories) |
| `< 1.0.0` | :x: | Unsupported (development / pre-release milestones) |

---

## Reporting a Vulnerability

If you discover a security vulnerability or potential exploit in VarDictRust, **please do not disclose it in a public GitHub issue, pull request, or discussion forum.** Public disclosure before a patch is available puts clinical and research pipelines at risk.

### Preferred Method: GitHub Private Vulnerability Reporting
Please report vulnerabilities privately using GitHub's native Security Advisory interface:
1. Navigate to the **Security** tab of the VarDictRust repository on GitHub.
2. Under "Reporting", click **Report a vulnerability** (or use the direct link: `https://github.com/<owner>/VarDictRust/security/advisories/new`).
3. Fill out the report with:
   - A clear description of the vulnerability and its potential impact (e.g. denial of service, memory exhaustion, arbitrary execution, data tampering).
   - Minimal reproduction steps (including minimal FASTA/BED/BAM inputs or synthetic fixtures).
   - Any proposed remediation or patches.

Reports submitted through Private Vulnerability Reporting are confidential, encrypted, and accessible only to project maintainers and the reporter until mutual coordination and disclosure.

---

## Response & Remediation SLA

We commit to the following Service Level Agreement (SLA) for confirmed reports:

- **Initial Acknowledgement:** Within **48 hours** of report receipt.
- **Triage & Reproduction:** Within **5 business days** confirming vulnerability severity (CVSS score).
- **Patch Development & Testing:** Coordinated within a private security advisory branch with full regression test verification.
- **Release & Coordinated Disclosure:**
  - **Critical / High Severity:** Patch release within **14–30 calendar days**.
  - **Medium / Low Severity:** Patch release within **60 calendar days**.
  - A formal GitHub Security Advisory (GHSA) and CVE identifier will be published concurrently with the patch release.

---

## Scope & Non-Vulnerabilities

### In Scope
- Memory safety violations, buffer overflows, or panics in the variant calling engine when parsing untrusted BAM/CRAM/FASTA/BED inputs.
- Unbounded memory consumption (OOM) or infinite loops triggered by malformed CIGAR strings, truncated alignments, or aberrant sequence records.
- Cryptographic or supply-chain vulnerabilities in dependencies (monitored via `cargo-audit` and `cargo-deny`).
- CI/CD workflow security, GitHub Actions permission escalation, or release provenance tampering.

### Out of Scope (Scientific Parity Heuristics)
- Differences in variant calling sensitivity, specificity, or allele frequency thresholds that reflect intentional adherence to the **VarDictJava v1.8.3 reference implementation** (the project's prime directive).
- Known upstream VarDictJava algorithmic quirks or heuristics (which are documented in [`DIVERGENCES.md`](../DIVERGENCES.md) and governed post-1.0 by the `--fix <id>` opt-in flags per [`ROADMAP.md`](../ROADMAP.md)).
