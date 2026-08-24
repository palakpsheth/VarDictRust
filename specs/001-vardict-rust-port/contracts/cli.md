# CLI Interface Contract: VarDictRust

This document specifies the command-line interface contract for the `vardictrust` binary.

## 1. Parameters

The CLI supports all standard short options of VarDictJava v1.8.3, plus long options and subcommands.

| Short | Long | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `-G` | `--ref` | `Path` | *Required* | Path to the reference genome FASTA file. |
| `-b` | `--bam` | `String` | *Required* | Path to BAM/CRAM file. For somatic mode, specify as `"tumor.bam\|normal.bam"`. |
| `-f` | `--af` | `Float` | `0.01` | Minimum allele frequency threshold. |
| `-N` | `--sample` | `String` | `"SAMPLE"` | Sample name. For somatic mode, specify as `"TUMOR\|NORMAL"`. |
| `-R` | `--region` | `String` | `None` | Process specific region (format: `chr:start-end`). |
| `-c` | `--chrom-col`| `Int` | `1` | BED column number for chromosome (1-based). |
| `-S` | `--start-col`| `Int` | `2` | BED column number for region start (1-based). |
| `-E` | `--end-col` | `Int` | `3` | BED column number for region end (1-based). |
| `-g` | `--name-col` | `Int` | `4` | BED column number for region gene name (1-based). |
| `-t` | `--threads`  | `Int` | Max Cores| Number of threads to use for parallel tiling. |
| `-k` | `--realignment`| `Int` | `1` | Enable local realignment (1) or disable (0). |
| *None*| `--fisher` | `Flag` | `false` | Enable calculation of Fisher p-values and odds ratios. |
| *None*| `--emit-intermediate`| `Flag`| `false`| Emit the legacy intermediate `.var` table format instead of VCF. |

## 2. Expected Outputs

### VCF Output (Default)
Emitted to standard output (`stdout`). Must strictly conform to VCF v4.2 format, including contig lengths fetched from reference `.fai` index, correct format fields, genotypes, and somatic statuses.

### Legacy Intermediate Table (`--emit-intermediate`)
Emitted to standard output (`stdout`) when requested. Column counts must match the oracle modes:
- Simple Mode: **36** columns (without `--fisher`) or **38** columns (with `--fisher`).
- Amplicon Mode: **38** columns (without `--fisher`) or **40** columns (with `--fisher`).
- Somatic Mode: **55** columns (without `--fisher`) or **61** columns (with `--fisher`).

## 3. Exit Status Codes

- `0`: Variant calling completed successfully.
- `1`: Execution error (e.g. file missing, invalid formatting).
- `2`: Unsupported argument or config option.
