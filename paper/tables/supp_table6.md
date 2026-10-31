### Supplementary Table S6 (Table S6). 13-Mode Functional Parity Matrix and CPP-1..10 vardictcpp Divergence Taxonomy

| Source Cell ID | Mode / Divergence ID | VarDictJava Rows / Behavior | VarDictRust Rows / Behavior | vardictcpp Rows / Behavior | Verification Verdict & Root Cause |
| :--- | :--- | :--- | :--- | :--- | :--- |
| s2_tierp_mode_parity_sweep | simple_f0.01 | 228 | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | simple_f0_mrd | 248 | 248 (FP=0, FN=0) | 248 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | simple_chunk50000 | 228 | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | simple_fisher | 650 | 650 (FP=0, FN=0) | 650 (FP=1, FN=1) | MATCH |
| s2_tierp_mode_parity_sweep | targeted_panel_chr17 | 83 | 83 (FP=0, FN=0) | 83 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | sv_deletion_regions | 3 | 3 (FP=0, FN=0) | 3 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | low_coverage | 146 | 146 (FP=0, FN=0) | 146 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | amplicon_c17amp | 145 | 145 (FP=0, FN=0) | 145 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | amplicon_fisher | 145 | 145 (FP=0, FN=0) | 145 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | somatic_f0.01_big | 664 | 664 (FP=0, FN=0) | 664 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | somatic_f0_mrd_c17 | 1396 | 1,396 (FP=0, FN=0) | 1,396 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | somatic_fisher_c17 | 527 | 527 (FP=0, FN=0) | 527 (FP=0, FN=0) | MATCH |
| s2_tierp_mode_parity_sweep | P5_cram_twin | 0_NoClassDefFoundError | 228 (FP=0, FN=0) | 228 (FP=0, FN=0) | MATCH_BAM_TWIN |
| s5_h4_maintenance_security_matrix | CPP-1: VCF Header & INFO/END Schema | Standard VCF 4.2 + END tag | Exact v1.8.3 header & END parity | Malformed/omitted END & header lines | Schema mismatch in downstream parsers |
| s5_h4_maintenance_security_matrix | CPP-2: CRISPR (-J/-j) & Structural Modes | Supported | Supported (100% MATCH) | Unimplemented / drops flags | Verified on extended-coverage-tier2 |
| s5_h4_maintenance_security_matrix | CPP-3: Amplicon (-a) Output Format | 38/40-col .var + VCF | Exact 38/40-col .var + VCF MATCH | Emits raw TSV instead of VCF | Breaks compare_vcf.py & bcftools |
| s5_h4_maintenance_security_matrix | CPP-4: Concatenated CLI Flag Parsing | Supported | Supported (60-flag parity) | Exits with code 120 on bundled flags | CLI wrapper incompatibility |
| s5_h4_maintenance_security_matrix | CPP-5: Chunk Boundary Deduplication (--chunk) | Single-region (or OOM) | 100 kbp halo tiling + cross-tile dedup | Emits +309,569 duplicate boundary calls on HG002 Arm A | Uncoordinated 50 kbp chunk splitting |
| s5_h4_maintenance_security_matrix | CPP-6: Missing ##contig Header Lines | Emits ##contig from .fai | Emits exact ##contig lines from .fai | Omits ##contig header lines | Requires header injection for rtg vcfeval |
| s5_h4_maintenance_security_matrix | CPP-7: GRCh37 / b37 Chromosome Naming | Supported | Supported (NA24385-b37 EQUAL) | Hardcoded chr-prefix assumptions | Tested on na24385-b37-truth |
| s5_h4_maintenance_security_matrix | CPP-8: BAM/CRAM Reference & Decode | BAM (CRAM needs classpath) | Native BAM + CRAM (-G) zero-diff | Fails or diverges on CRAM twins | Verified on hg008-cram-twin & Ultima CRAM |
| s5_h4_maintenance_security_matrix | CPP-9: MNP Splitting & Somatic SV Drops | Reference pileup & realigner | Exact parity with VarDictJava | Differs on 394 HG002 Arm B calls & drops SVs | Realignment & tie-break divergence |
| s5_h4_maintenance_security_matrix | CPP-10: RNA-Seq Splice Mode (-i) | Supported (12,668 rows on HG002 mRNA) | Supported (12,668 rows 100% MATCH) | Unimplemented (SKIP) | Verified on 51.8 GiB HG002 STAR BAM |

*Every divergence class was empirically observed and cataloged during the 60-cell GCP evaluation.*
