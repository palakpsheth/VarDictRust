#!/usr/bin/env bash
# =============================================================================
# packaging/nf-core/prepare_nfcore_pr.sh
# Automates updating the Nextflow nf-core/modules vardictrust module after a
# Bioconda release is live on anaconda.org for both linux-64 and linux-aarch64,
# including the v1.1.0+ pure-MIT meta.yml licence transition and automatic
# Seqera Wave container URI synchronization into README.md.
# =============================================================================
set -euo pipefail

if [ -n "${HOME:-}" ] && [ -d "$HOME/.local/bin" ]; then
  case ":${PATH:-}:" in
    *":$HOME/.local/bin:"*) ;;
    *) export PATH="${PATH:+$PATH:}$HOME/.local/bin" ;;
  esac
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/harness/lib/portable.sh"
IN_TREE_MODULE_DIR="${IN_TREE_MODULE_DIR:-$SCRIPT_DIR/modules/vardictrust}"
README_MD="${README_MD:-$ROOT_DIR/README.md}"

update_nfcore_module_metadata() {
  local env_yml_path="$1"
  local meta_yml_path="$2"
  local ver="$3"

  python3 - "$env_yml_path" "$meta_yml_path" "$ver" <<'PY'
import pathlib, re, sys

env_path = pathlib.Path(sys.argv[1])
meta_path = pathlib.Path(sys.argv[2])
ver = sys.argv[3]
ver_tuple = tuple(int(x) for x in ver.split("."))

env_txt, count = re.subn(
    r"bioconda::vardictrust=[0-9]+\.[0-9]+\.[0-9]+",
    f"bioconda::vardictrust={ver}",
    env_path.read_text(encoding="utf-8"),
)
if count != 1:
    sys.exit(f"ERROR: expected 1 match for bioconda::vardictrust in {env_path}, got {count}")
env_path.write_text(env_txt, encoding="utf-8")

if ver_tuple >= (1, 1, 0):
    meta_txt = meta_path.read_text(encoding="utf-8")
    pat = re.compile(
        r'(^\s*-\s*"vardictrust":\s*\n(?:(?!^\s*-\s*"[A-Za-z0-9_-]+":)[^\n]*\n)*?'
        r"^(\s*)licence:\s*\n)"
        r"(\s*-\s*\"[^\n]+\"\s*\n)+",
        re.MULTILINE,
    )

    def _repl(m: re.Match[str]) -> str:
        item_indent = m.group(2) + "  "
        return f'{m.group(1)}{item_indent}- "MIT"\n'

    meta_txt, lic_count = pat.subn(_repl, meta_txt, count=1)
    if lic_count != 1:
        sys.exit(
            f"ERROR: expected 1 vardictrust licence block in {meta_path}, got {lic_count}"
        )
    if "GPL" in meta_txt:
        sys.exit(f"ERROR: unexpected residual GPL token in {meta_path} for pure-MIT v{ver}")
    if not re.search(
        r'^\s*-\s*"htslib":[\s\S]*?^\s*licence:\s*\n\s*-\s*"MIT"',
        meta_txt,
        re.MULTILINE,
    ):
        sys.exit(f"ERROR: htslib MIT licence entry missing or corrupted in {meta_path}")
    meta_path.write_text(meta_txt, encoding="utf-8")
PY
}

sync_readme_wave_uris() {
  local main_nf_path="$1"
  local meta_yml_path="$2"
  local readme_path="$3"

  python3 - "$main_nf_path" "$meta_yml_path" "$readme_path" <<'PY'
import pathlib, re, sys

main_nf = pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")
meta_yml = pathlib.Path(sys.argv[2]).read_text(encoding="utf-8")
readme_file = pathlib.Path(sys.argv[3])
readme_txt = readme_file.read_text(encoding="utf-8")

nf_sing_m = re.search(
    r"https://community-cr-prod\.seqera\.io/docker/registry/v2/blobs/sha256/[0-9a-f]{2}/[0-9a-f]{64}/data",
    main_nf,
)
nf_dock_m = re.search(
    r"(?<!oras://)community\.wave\.seqera\.io/library/htslib_vardictrust:([0-9a-f]{16})",
    main_nf,
)
if not nf_sing_m or not nf_dock_m:
    sys.exit(
        f"ERROR: could not extract Singularity HTTPS and Docker Wave URIs from {sys.argv[1]}"
    )
nf_sing_https = nf_sing_m.group(0)
nf_docker_amd64_tag = nf_dock_m.group(1)

containers_idx = meta_yml.find("\ncontainers:\n")
if containers_idx == -1 and not meta_yml.startswith("containers:\n"):
    sys.exit(f"ERROR: missing containers: section in {sys.argv[2]}")
containers_block = meta_yml[containers_idx:] if containers_idx != -1 else meta_yml

def _extract_arch_field(section: str, arch: str, field_pat: str) -> str:
    sec_m = re.search(
        rf"^\s{{2}}{re.escape(section)}:\s*\n((?:^\s{{4,}}[^\n]*\n)+)",
        containers_block,
        re.MULTILINE,
    )
    if not sec_m:
        sys.exit(f"ERROR: missing containers.{section} block in {sys.argv[2]}")
    arch_m = re.search(
        rf"^\s{{4}}{re.escape(arch)}:\s*\n((?:^\s{{6,}}[^\n]*\n)+)",
        sec_m.group(1),
        re.MULTILINE,
    )
    if not arch_m:
        sys.exit(f"ERROR: missing containers.{section}.{arch} block in {sys.argv[2]}")
    val_m = re.search(field_pat, arch_m.group(1), re.MULTILINE)
    if not val_m:
        sys.exit(
            f"ERROR: could not parse {field_pat!r} under containers.{section}.{arch} in {sys.argv[2]}"
        )
    return val_m.group(1)

meta_docker_amd64 = _extract_arch_field(
    "docker",
    "linux/amd64",
    r"^\s*name:\s*community\.wave\.seqera\.io/library/htslib_vardictrust:([0-9a-f]{16})\s*$",
)
meta_docker_arm64 = _extract_arch_field(
    "docker",
    "linux/arm64",
    r"^\s*name:\s*community\.wave\.seqera\.io/library/htslib_vardictrust:([0-9a-f]{16})\s*$",
)
meta_sing_amd64 = _extract_arch_field(
    "singularity",
    "linux/amd64",
    r"^\s*name:\s*oras://community\.wave\.seqera\.io/library/htslib_vardictrust:([0-9a-f]{16})\s*$",
)
meta_sing_amd64_https = _extract_arch_field(
    "singularity",
    "linux/amd64",
    r"^\s*https:\s*(https://community-cr-prod\.seqera\.io/docker/registry/v2/blobs/sha256/[0-9a-f]{2}/[0-9a-f]{64}/data)\s*$",
)
meta_sing_arm64 = _extract_arch_field(
    "singularity",
    "linux/arm64",
    r"^\s*name:\s*oras://community\.wave\.seqera\.io/library/htslib_vardictrust:([0-9a-f]{16})\s*$",
)

if nf_docker_amd64_tag != meta_docker_amd64:
    sys.exit(
        f"ERROR: main.nf Docker tag ({nf_docker_amd64_tag}) != meta.yml containers.docker.linux/amd64 ({meta_docker_amd64})"
    )
if nf_sing_https != meta_sing_amd64_https:
    sys.exit(
        f"ERROR: main.nf Singularity HTTPS URL ({nf_sing_https}) != meta.yml containers.singularity.linux/amd64.https ({meta_sing_amd64_https})"
    )

readme_txt, c_dock_hdr = re.subn(
    r"(# Docker / Podman \(linux/amd64:\s*)[0-9a-f]{16}(\s*\|\s*linux/arm64:\s*)[0-9a-f]{16}(\))",
    rf"\g<1>{meta_docker_amd64}\g<2>{meta_docker_arm64}\g<3>",
    readme_txt,
)
readme_txt, c_dock_uri = re.subn(
    r"(?<!oras://)community\.wave\.seqera\.io/library/htslib_vardictrust:[0-9a-f]{16}",
    f"community.wave.seqera.io/library/htslib_vardictrust:{nf_docker_amd64_tag}",
    readme_txt,
)
readme_txt, c_sing_hdr = re.subn(
    r"(# Apptainer / Singularity ORAS \(linux/amd64:\s*)[0-9a-f]{16}(\s*\|\s*linux/arm64:\s*)[0-9a-f]{16}(\))",
    rf"\g<1>{meta_sing_amd64}\g<2>{meta_sing_arm64}\g<3>",
    readme_txt,
)
readme_txt, c_sing_oras = re.subn(
    r"oras://community\.wave\.seqera\.io/library/htslib_vardictrust:[0-9a-f]{16}",
    f"oras://community.wave.seqera.io/library/htslib_vardictrust:{meta_sing_amd64}",
    readme_txt,
)
readme_txt, _ = re.subn(
    r"https://community-cr-prod\.seqera\.io/docker/registry/v2/blobs/sha256/[0-9a-f]{2}/[0-9a-f]{64}/data",
    nf_sing_https,
    readme_txt,
)

if c_dock_hdr != 1 or c_sing_hdr != 1 or c_dock_uri < 1 or c_sing_oras < 1:
    sys.exit(
        f"ERROR: failed to synchronize Wave container URIs in {readme_file} "
        f"(dock_hdr={c_dock_hdr}, dock_uri={c_dock_uri}, sing_hdr={c_sing_hdr}, sing_oras={c_sing_oras})"
    )

readme_file.write_text(readme_txt, encoding="utf-8")
PY
}

VERSION="${1:-1.0.10}"
if [ "${1:-}" = "--hermetic-sync" ]; then
  if [ $# -ne 6 ]; then
    echo "Usage: prepare_nfcore_pr.sh --hermetic-sync <VERSION> <ENV_YML> <META_YML> <MAIN_NF> <README_MD>" >&2
    exit 1
  fi
  HERMETIC_VER="$2"
  if [[ ! "$HERMETIC_VER" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "ERROR: Invalid VERSION '${HERMETIC_VER}' (expected strict SemVer X.Y.Z)." >&2
    exit 1
  fi
  update_nfcore_module_metadata "$3" "$4" "$HERMETIC_VER"
  sync_readme_wave_uris "$5" "$4" "$6"
  exit 0
fi

: "${NFCORE_MODULES_DIR:?Set NFCORE_MODULES_DIR to a trusted local nf-core/modules checkout path}"
NFCORE_MODULE_PATH="$NFCORE_MODULES_DIR/modules/nf-core/vardictrust"

if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "ERROR: Invalid VERSION '${VERSION}' (expected strict SemVer X.Y.Z)." >&2
  exit 1
fi
ANACONDA_API="https://api.anaconda.org/package/bioconda/vardictrust"

for tool in curl python3 nf-core nextflow nf-test docker; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: Required tool '$tool' not found on PATH ($PATH)" >&2
    exit 1
  fi
done

echo "Checking Bioconda availability for vardictrust=${VERSION} (linux-64 + linux-aarch64)..."
ANACONDA_JSON="$(curl --proto '=https' --tlsv1.2 -fsSL --connect-timeout 10 --max-time 30 --retry 3 "$ANACONDA_API")"

if ! python3 -c '
import json, sys
data = json.loads(sys.stdin.read())
ver = sys.argv[1]
versions = data.get("versions", [])
if ver not in versions:
    print(f"Version {ver} not in anaconda.org versions: {versions}", file=sys.stderr)
    sys.exit(1)
subdirs = {
    f.get("attrs", {}).get("subdir")
    for f in data.get("files", [])
    if f.get("version") == ver
}
required = {"linux-64", "linux-aarch64"}
missing = required - subdirs
if missing:
    print(f"Version {ver} is missing required architectures on anaconda.org: {sorted(missing)} (found: {sorted(subdirs)})", file=sys.stderr)
    sys.exit(1)
print(f"Verified bioconda::vardictrust={ver} available for {sorted(subdirs)}.")
' "$VERSION" <<<"$ANACONDA_JSON"; then
  echo "================================================================================" >&2
  echo "ERROR: bioconda::vardictrust=${VERSION} is not yet live for both linux-64 and" >&2
  echo "linux-aarch64 on anaconda.org ($ANACONDA_API)." >&2
  echo "Wait for the bioconda/bioconda-recipes PR to merge and finish uploading before" >&2
  echo "generating Wave containers with nf-core modules containers create." >&2
  echo "================================================================================" >&2
  exit 1
fi

if [[ ! -d "$NFCORE_MODULES_DIR/.git" ]]; then
  echo "ERROR: nf-core/modules clone not found at $NFCORE_MODULES_DIR" >&2
  echo "Set NFCORE_MODULES_DIR=/path/to/nf-core-modules and re-run." >&2
  exit 1
fi

if [[ -n "$(git -C "$NFCORE_MODULES_DIR" status --porcelain)" ]]; then
  echo "ERROR: $NFCORE_MODULES_DIR has uncommitted changes. Commit or stash first." >&2
  exit 1
fi

# Sync in-tree module into nf-core/modules checkout and update environment.yml + (v1.1.0+) pure-MIT licence
rm -rf "$NFCORE_MODULE_PATH"
mkdir -p "$NFCORE_MODULE_PATH"
cp -a "$IN_TREE_MODULE_DIR/." "$NFCORE_MODULE_PATH/"

update_nfcore_module_metadata "$NFCORE_MODULE_PATH/environment.yml" "$NFCORE_MODULE_PATH/meta.yml" "$VERSION"

# Remove previous Wave conda-lock files so containers create writes fresh lockfiles only
rm -f "$NFCORE_MODULE_PATH/.conda-lock/"*.txt

(
  cd "$NFCORE_MODULES_DIR"
  export NXF_VER="${NXF_VER:-25.10.4}"
  echo "Provisioning multi-arch Seqera Wave containers for vardictrust=${VERSION}..."
  nf-core modules containers create vardictrust

  echo "Running nf-test --update-snapshot and stability verification (--profile docker)..."
  nf-test test modules/nf-core/vardictrust/tests/main.nf.test --update-snapshot --profile docker
  nf-test test modules/nf-core/vardictrust/tests/main.nf.test --profile docker

  echo "Running nf-core modules lint and nextflow v2 syntax lint..."
  nf-core modules lint vardictrust
  NXF_SYNTAX_PARSER=v2 nextflow lint modules/nf-core/vardictrust/
)

# Sync updated Wave container URLs, .conda-lock files, and snapshots back in-tree
rm -rf "$IN_TREE_MODULE_DIR"
mkdir -p "$IN_TREE_MODULE_DIR"
cp -a "$NFCORE_MODULE_PATH/." "$IN_TREE_MODULE_DIR/"

# Synchronize multi-arch Seqera Wave container URIs from main.nf + meta.yml into README.md
sync_readme_wave_uris "$IN_TREE_MODULE_DIR/main.nf" "$IN_TREE_MODULE_DIR/meta.yml" "$README_MD"

echo "Successfully updated $IN_TREE_MODULE_DIR, $NFCORE_MODULE_PATH, and $README_MD to v${VERSION}."
