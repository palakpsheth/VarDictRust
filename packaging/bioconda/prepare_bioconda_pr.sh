#!/usr/bin/env bash
# =============================================================================
# packaging/bioconda/prepare_bioconda_pr.sh
# Automates preparing a pull request to github.com/bioconda/bioconda-recipes.
# Computes the release tarball SHA256 checksum and updates meta.yaml.
# Fails loudly per Rule 0/0b if the tag does not exist on GitHub.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/harness/lib/portable.sh"
RECIPE_DIR="$SCRIPT_DIR/recipes/vardictrust"
META_YAML="$RECIPE_DIR/meta.yaml"

VERSION="${1:-1.0.9}"
if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "ERROR: Invalid VERSION '${VERSION}' (expected strict SemVer X.Y.Z)." >&2
  exit 1
fi
TARBALL_URL="https://github.com/palakpsheth/VarDictRust/archive/refs/tags/v${VERSION}.tar.gz"

echo "Preparing Bioconda recipe for VarDictRust v${VERSION}..."
echo "Fetching release tarball from GitHub: ${TARBALL_URL}"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

if ! curl --proto '=https' --tlsv1.2 -LsSf --connect-timeout 10 --max-time 60 --retry 3 "$TARBALL_URL" -o "$TMP_DIR/release.tar.gz"; then
  echo "================================================================================" >&2
  echo "ERROR: Could not fetch release tarball from GitHub!" >&2
  echo "URL: $TARBALL_URL" >&2
  echo "Per Rule 0b: We NEVER use a local git archive fallback because local gzip" >&2
  echo "compression produces a different SHA256 than GitHub's release tarball." >&2
  echo "Action: Ensure tag v${VERSION} is pushed to GitHub before preparing Bioconda PR:" >&2
  echo "  git push origin v${VERSION}" >&2
  echo "================================================================================" >&2
  exit 1
fi

SHA256="$(portable_sha256 "$TMP_DIR/release.tar.gz")"
echo "Computed live GitHub release SHA256: ${SHA256}"

# Update sha256, version, and reset build number to 0 in meta.yaml
python3 - "$META_YAML" "$VERSION" "$SHA256" <<'PY'
import pathlib, re, sys
path = pathlib.Path(sys.argv[1])
ver, sha = sys.argv[2], sys.argv[3]
txt = path.read_text(encoding="utf-8")
for pat, repl, desc in [
    (r"sha256: [0-9a-fA-F]+", f"sha256: {sha}", "sha256"),
    (r'\{% set version = "[^"]*" %\}', f'{{% set version = "{ver}" %}}', "version"),
    (r"(?m)^([ \t]*number:[ \t]*)[0-9]+", r"\g<1>0", "build.number"),
]:
    txt, count = re.subn(pat, repl, txt)
    if count != 1:
        sys.exit(f"ERROR: expected 1 match for {desc} in {path}, got {count}")
path.write_text(txt, encoding="utf-8")
PY

echo "Updated $META_YAML with version ${VERSION} and sha256 ${SHA256}."
echo ""
echo "Next steps to publish on Bioconda:"
echo "1. Fork https://github.com/bioconda/bioconda-recipes"
echo "2. Clone your fork: git clone git@github.com:<your-user>/bioconda-recipes.git"
echo "3. Branch: cd bioconda-recipes && git checkout -b add-vardictrust"
echo "4. Copy recipe: cp -r $RECIPE_DIR recipes/vardictrust"
echo "5. Commit: git add recipes/vardictrust && git commit -m 'add vardictrust ${VERSION}'"
echo "6. Push and PR: git push origin add-vardictrust and open PR against master"
echo "7. Comment '@BiocondaBot please build' and await human review from @bioconda/core"
