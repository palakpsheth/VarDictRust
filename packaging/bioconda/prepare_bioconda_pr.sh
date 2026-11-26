#!/usr/bin/env bash
# =============================================================================
# packaging/bioconda/prepare_bioconda_pr.sh
# Automates preparing a pull request to github.com/bioconda/bioconda-recipes.
# Computes the release tarball SHA256 checksum and updates meta.yaml
# (including the v1.1.0+ pure-MIT license & 32-switch/preset metadata update).
# Fails loudly per Rule 0/0b if the tag does not exist on GitHub.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$ROOT_DIR/harness/lib/portable.sh"
RECIPE_DIR="$SCRIPT_DIR/recipes/vardictrust"
META_YAML="$RECIPE_DIR/meta.yaml"

if [[ "${1:-}" == -* ]]; then
  set -- "" "$@"
fi
VERSION="${1:-1.1.0}"
if [ $# -gt 0 ]; then
  shift
fi

OFFLINE_SHA256=""
while [ $# -gt 0 ]; do
  case "$1" in
    --version)
      if [ $# -lt 2 ]; then
        echo "ERROR: --version requires a SemVer X.Y.Z argument." >&2
        exit 1
      fi
      VERSION="$2"
      shift 2
      ;;
    --offline-sha256)
      if [ $# -lt 2 ]; then
        echo "ERROR: --offline-sha256 requires a 64-character hex SHA256 digest." >&2
        exit 1
      fi
      OFFLINE_SHA256="$2"
      shift 2
      ;;
    --meta-yaml)
      if [ $# -lt 2 ]; then
        echo "ERROR: --meta-yaml requires a target meta.yaml file path." >&2
        exit 1
      fi
      META_YAML="$2"
      shift 2
      ;;
    *)
      echo "ERROR: Unknown argument '$1' (usage: prepare_bioconda_pr.sh [VERSION] [--offline-sha256 <HEX64>] [--meta-yaml <PATH>])." >&2
      exit 1
      ;;
  esac
done

if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
  echo "ERROR: Invalid VERSION '${VERSION}' (expected strict SemVer X.Y.Z)." >&2
  exit 1
fi

if [ ! -f "$META_YAML" ]; then
  echo "ERROR: Target recipe file not found: $META_YAML" >&2
  exit 1
fi

update_bioconda_meta_yaml() {
  local meta_path="$1"
  local ver="$2"
  local sha="$3"

  python3 - "$meta_path" "$ver" "$sha" <<'PY'
import pathlib, re, sys

path = pathlib.Path(sys.argv[1])
ver = sys.argv[2]
sha = sys.argv[3].lower()
if not re.fullmatch(r"[0-9a-f]{64}", sha):
    sys.exit(f"ERROR: invalid 64-hex SHA256 digest: {sys.argv[3]!r}")

ver_tuple = tuple(int(x) for x in ver.split("."))
txt = path.read_text(encoding="utf-8")

replacements = [
    (r"sha256: [0-9a-fA-F]+", f"sha256: {sha}", "sha256"),
    (r'\{% set version = "[^"]*" %\}', f'{{% set version = "{ver}" %}}', "version"),
    (r"(?m)^([ \t]*number:[ \t]*)[0-9]+", r"\g<1>0", "build.number"),
]

if ver_tuple >= (1, 1, 0):
    v11_description = (
        "\\g<1>"
        "    VarDictRust is a fast, memory-safe Rust somatic and germline variant caller\n"
        "    that consolidates the VarDictJava v1.8.3 + R + Perl pipeline into a single\n"
        "    MIT-licensed binary with bounded-memory WGS tiling, default byte-identical\n"
        "    VarDictJava v1.8.3 parity on valid inputs (--compat=vardictjava), native\n"
        "    multithreaded BGZF and Tabix/CSI VCF output, opt-in --compat=off correctness,\n"
        "    SV, and SAM =/X (--eqx) alignment fixes (32 switches), clinical workflow\n"
        "    presets (--preset wgs|wes|amplicon|ctdna), --amplicon-single-primer support,\n"
        "    and 18/18 var2vcf post-filtering including --vcf-all-alleles multi-allelic\n"
        "    VCF output.\n"
    )
    replacements.extend(
        [
            (r"(?m)^([ \t]*license:[ \t]*).+$", r"\g<1>MIT", "about.license"),
            (r"(?m)^([ \t]*license_family:[ \t]*).+$", r"\g<1>MIT", "about.license_family"),
            (
                r"(?ms)^([ \t]*description:[ \t]*\|[ \t]*\n)(?:[ \t]{4,}[^\n]*\n)+",
                v11_description,
                "about.description",
            ),
        ]
    )

for pat, repl, desc in replacements:
    txt, count = re.subn(pat, repl, txt)
    if count != 1:
        sys.exit(f"ERROR: expected 1 match for {desc} in {path}, got {count}")

if ver_tuple >= (1, 1, 0) and "GPL" in txt:
    sys.exit(f"ERROR: unexpected residual GPL token in {path} for pure-MIT v{ver}")

path.write_text(txt, encoding="utf-8")
PY
}

if [ -n "$OFFLINE_SHA256" ]; then
  if [[ ! "$OFFLINE_SHA256" =~ ^[0-9a-fA-F]{64}$ ]]; then
    echo "ERROR: Invalid --offline-sha256 '${OFFLINE_SHA256}' (expected 64 hex characters)." >&2
    exit 1
  fi
  SHA256="$(printf '%s' "$OFFLINE_SHA256" | tr '[:upper:]' '[:lower:]')"
  echo "Preparing Bioconda recipe for VarDictRust v${VERSION} (hermetic offline sha256)..."
else
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
fi

# Update sha256, version, build number, and (for v1.1.0+) pure-MIT metadata in meta.yaml
update_bioconda_meta_yaml "$META_YAML" "$VERSION" "$SHA256"

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
