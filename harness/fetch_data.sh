#!/usr/bin/env bash
# Usage: fetch_data.sh [--allow-legacy-md5] <group>...  (refs seqc2 hg008 hg002). Resumable; sha256-verified; appends benchmarks/MANIFEST.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/portable.sh"

ALLOW_LEGACY_MD5=0
groups=()
for arg in "$@"; do
  case "$arg" in
    --allow-legacy-md5)
      ALLOW_LEGACY_MD5=1
      ;;
    -*)
      echo "Unknown option: $arg" >&2
      exit 1
      ;;
    *)
      groups+=("$arg")
      ;;
  esac
done

if [ "${#groups[@]}" -eq 0 ]; then
  echo "Usage: $0 [--allow-legacy-md5] <group>...  (refs seqc2 hg008 hg002)" >&2
  exit 1
fi

DATA=${DATA_DIR:-scratch/g4}
MAN=${MANIFEST:-benchmarks/MANIFEST}
ASSETS=${ASSETS_TSV:-harness/assets.tsv}

mkdir -p "$DATA" "$(dirname "$MAN")"
[ -f "$MAN" ] || printf '# name\turl\tsha256\tchecksum_source\tbytes\tfetched\n' > "$MAN"

want=" ${groups[*]} "
while IFS=$'\t' read -r group name url md5 sha256 dest extra; do
  case "$group" in ''|\#*) continue ;; esac
  [[ "$want" == *" $group "* ]] || continue
  out="$DATA/$dest"
  mkdir -p "$(dirname "$out")"

  # Step A (Fail-closed check FIRST):
  if [ -z "$sha256" ] || [ "$sha256" = "-" ]; then
    if [ "$ALLOW_LEGACY_MD5" -ne 1 ]; then
      echo "fetch: missing required trusted SHA-256 for $name (use --allow-legacy-md5 to accept legacy MD5)" >&2
      exit 1
    fi
  fi

  # Step B (Digest-bound markers $out.sha256 vs $out.legacy_md5):
  if [ -n "$sha256" ] && [ "$sha256" != "-" ]; then
    # Fast marker shortcut ONLY when .sha256 marker matches expected sha256 and target (or ref index) exists
    if [ -f "$out.sha256" ] && [ "$(cat "$out.sha256" 2>/dev/null)" = "$sha256" ] && { [ -f "$out" ] || { [ "$group" = "refs" ] && [ -s "${out}.fai" ]; }; }; then
      echo "have: $name"
      continue
    fi

    # Pre-migration / missing .sha256 marker promotion:
    if [ -f "$out" ]; then
      got=$(portable_sha256 "$out")
      if [ "$got" = "$sha256" ]; then
        printf '%s\n' "$sha256" > "$out.sha256"
        touch "$out.done"
        if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "upstream" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
        fi
        echo "verified: $name"
        continue
      elif [ "$group" = "refs" ] && [ -f "$out.done" ] && [ -s "${out}.fai" ]; then
        echo "have (transformed ref): $name"
        continue
      fi
      # For any other file where got != sha256, fall through to re-fetch/verify with stage_verified_sha256
    fi

    echo "fetch: $name"
    if ! stage_verified_sha256 "$url" "$sha256" "$out"; then
      echo "fetch: failed to fetch/verify $name" >&2
      exit 1
    fi
    printf '%s\n' "$sha256" > "$out.sha256"
    touch "$out.done"
    if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$sha256" "upstream" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
    fi
  else
    # Legacy MD5 row (sha256 = "-", already passed ALLOW_LEGACY_MD5 == 1)
    if [ -z "$md5" ] || [ "$md5" = "-" ]; then
      echo "fetch: missing required checksum for $name" >&2
      exit 1
    fi

    # Fast marker shortcut ONLY when .legacy_md5 marker matches expected md5 and target exists
    if [ -f "$out.legacy_md5" ] && [ "$(cat "$out.legacy_md5" 2>/dev/null)" = "$md5" ] && [ -f "$out" ]; then
      echo "have: $name"
      continue
    fi

    if [ -f "$out" ]; then
      got=$(portable_md5 "$out")
      if [ "$got" = "$md5" ]; then
        printf '%s\n' "$md5" > "$out.legacy_md5"
        touch "$out.done"
        if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "legacy-md5" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
        fi
        echo "verified: $name"
        continue
      fi
    fi

    echo "fetch: $name"
    part="${out}.part"
    rm -f "$part"
    if ! curl -sSLf --retry 5 --retry-delay 5 -o "$part" "$url"; then
      echo "fetch: download failed for $url" >&2
      rm -f "$part"
      exit 1
    fi
    got=$(portable_md5 "$part")
    if [ "$got" != "$md5" ]; then
      echo "fetch: MD5 MISMATCH $name got=$got want=$md5" >&2
      rm -f "$part"
      exit 1
    fi
    mv -f "$part" "$out"
    printf '%s\n' "$md5" > "$out.legacy_md5"
    touch "$out.done"
    if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "legacy-md5" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
    fi
  fi
done < "$ASSETS"
