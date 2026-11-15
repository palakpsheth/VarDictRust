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

  if [ -f "$out.done" ]; then
    echo "have: $name"
    continue
  fi

  # In-place promotion for existing files missing .done
  if [ -f "$out" ]; then
    if [ -n "$sha256" ] && [ "$sha256" != "-" ]; then
      got=$(portable_sha256 "$out")
      if [ "$got" = "$sha256" ]; then
        echo "verified: $name"
        touch "$out.done"
        if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "upstream" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
        fi
        continue
      fi
    elif [ "$ALLOW_LEGACY_MD5" -eq 1 ] && [ -n "$md5" ] && [ "$md5" != "-" ]; then
      got=$(portable_md5 "$out")
      if [ "$got" = "$md5" ]; then
        echo "verified: $name"
        touch "$out.done"
        if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
          printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "legacy-md5" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
        fi
        continue
      fi
    fi
  fi

  if [ -z "$sha256" ] || [ "$sha256" = "-" ]; then
    if [ "$ALLOW_LEGACY_MD5" -ne 1 ]; then
      echo "fetch: missing required trusted SHA-256 for $name (use --allow-legacy-md5 to accept legacy MD5)" >&2
      exit 1
    fi
    if [ -z "$md5" ] || [ "$md5" = "-" ]; then
      echo "fetch: missing required checksum for $name" >&2
      exit 1
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
    touch "$out.done"
    if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "legacy-md5" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
    fi
  else
    echo "fetch: $name"
    if ! stage_verified_sha256 "$url" "$sha256" "$out"; then
      echo "fetch: failed to fetch/verify $name" >&2
      exit 1
    fi
    touch "$out.done"
    if ! awk -F'\t' -v n="$name" '$1 == n {found=1; exit} END{exit !found}' "$MAN"; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$sha256" "upstream" "$(portable_filesize "$out")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
    fi
  fi
done < "$ASSETS"
