#!/usr/bin/env bash
# =============================================================================
# packaging/nf-core/prepare_nfcore_pr.sh
# Automates updating the Nextflow nf-core/modules vardictrust module after a
# Bioconda release is live on anaconda.org for both linux-64 and linux-aarch64.
# =============================================================================
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN_TREE_MODULE_DIR="$SCRIPT_DIR/modules/vardictrust"
: "${NFCORE_MODULES_DIR:?Set NFCORE_MODULES_DIR to a trusted local nf-core/modules checkout path}"
NFCORE_MODULE_PATH="$NFCORE_MODULES_DIR/modules/nf-core/vardictrust"

VERSION="${1:-1.0.5}"
ANACONDA_API="https://api.anaconda.org/package/bioconda/vardictrust"

for tool in curl python3 nf-core nextflow nf-test docker; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: Required tool '$tool' not found on PATH ($PATH)" >&2
    exit 1
  fi
done

echo "Checking Bioconda availability for vardictrust=${VERSION} (linux-64 + linux-aarch64)..."
ANACONDA_JSON="$(curl -fsSL --connect-timeout 10 --max-time 30 --retry 3 "$ANACONDA_API")"

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

# Sync in-tree module into nf-core/modules checkout and update environment.yml version
rm -rf "$NFCORE_MODULE_PATH"
mkdir -p "$NFCORE_MODULE_PATH"
cp -a "$IN_TREE_MODULE_DIR/." "$NFCORE_MODULE_PATH/"

sed -i -E "s/bioconda::vardictrust=[0-9]+\.[0-9]+\.[0-9]+/bioconda::vardictrust=${VERSION}/" "$NFCORE_MODULE_PATH/environment.yml"

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

echo "Successfully updated packaging/nf-core/modules/vardictrust and $NFCORE_MODULE_PATH to v${VERSION}."
