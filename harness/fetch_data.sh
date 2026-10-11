#!/usr/bin/env bash
# Usage: fetch_data.sh <group>...  (refs seqc2 hg008 hg002). Resumable; md5-verified; appends benchmarks/MANIFEST.
set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 <group>...  (refs seqc2 hg008 hg002)" >&2
  exit 1
fi


DATA=${DATA_DIR:-scratch/g4}
MAN=${MANIFEST:-benchmarks/MANIFEST}
ASSETS=${ASSETS_TSV:-harness/assets.tsv}

mkdir -p "$DATA" "$(dirname "$MAN")"
[ -f "$MAN" ] || printf '# name\turl\tmd5\tmd5_source\tbytes\tfetched\n' > "$MAN"

want=" $* "
while IFS=$'\t' read -r group name url md5 dest; do
  case "$group" in ''|\#*) continue ;; esac
  [[ "$want" == *" $group "* ]] || continue
  out="$DATA/$dest"
  mkdir -p "$(dirname "$out")"
  if [ -f "$out.done" ]; then
    echo "have: $name"
    continue
  fi
  echo "fetch: $name"
  curl -fL --retry 5 --retry-delay 30 -C - -o "$out" "$url"
  got=$(md5sum "$out" | cut -d' ' -f1)
  if [ -n "$md5" ] && [ "$md5" != "-" ]; then
    if [ "$got" != "$md5" ]; then
      echo "fetch: MD5 MISMATCH $name got=$got want=$md5" >&2
      rm -f "$out"
      exit 1
    fi
    src=upstream
  else
    src=computed
  fi
  grep -qF "$name" "$MAN" || printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$got" "$src" "$(stat -c %s "$out")" "$(date -Is)" >> "$MAN"
  touch "$out.done"
done < "$ASSETS"
