#!/usr/bin/env bash
# Usage: fetch_data.sh [--allow-legacy-md5] <group>...  (refs seqc2 hg008 hg002).
# Atomic download; SHA-256 verified (fail-closed unless --allow-legacy-md5); appends benchmarks/MANIFEST.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/portable.sh"

allow_legacy_md5=0
groups=()
for arg in "$@"; do
  case "$arg" in
    --allow-legacy-md5)
      allow_legacy_md5=1
      ;;
    -h|--help)
      echo "Usage: $0 [--allow-legacy-md5] <group>...  (refs seqc2 hg008 hg002)"
      exit 0
      ;;
    -*)
      echo "fetch_data.sh: unknown option '$arg'" >&2
      echo "Usage: $0 [--allow-legacy-md5] <group>...  (refs seqc2 hg008 hg002)" >&2
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

record_manifest() {
  local m_name="$1" m_url="$2" m_digest="$3" m_src="$4" m_file="$5"
  if ! awk -F'\t' -v n="$m_name" -v u="$m_url" -v d="$m_digest" -v s="$m_src" \
    '!/^#/ && $1 == n && $2 == u && $3 == d && $4 == s { found=1; exit } END { exit !found }' "$MAN"; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$m_name" "$m_url" "$m_digest" "$m_src" "$(portable_filesize "$m_file")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$MAN"
  fi
}

# Verify unpacked outputs against pinned SHA-256 digests ONLY when a compressed `refs`
# archive was removed after `stage_g4.sh refs` extracted its contents.
refs_transformed_verified() {
  local r_group="$1" r_dest="$2" r_out="$3" r_sha256="$4"
  [ "$r_group" = "refs" ] || return 1
  [ ! -f "$r_out" ] || return 1
  case "$r_dest" in
    ref/d1vd1/GRCh38.d1.vd1_BWA.tar.gz)
      [ "$r_sha256" = "a354f980bc141e6ed8cb1b51aba48b989a41b8fb280b09e3f20fd66df97825c6" ] || return 1
      local bwa_dir
      bwa_dir="$(dirname "$r_out")"
      [ -f "$bwa_dir/GRCh38.d1.vd1.fa.amb" ] && [ "$(portable_sha256 "$bwa_dir/GRCh38.d1.vd1.fa.amb")" = "3e953413da3e09c29e7f99b9cb944da7a42f0402f110fe88dcade46e23ef665b" ] \
        && [ -f "$bwa_dir/GRCh38.d1.vd1.fa.ann" ] && [ "$(portable_sha256 "$bwa_dir/GRCh38.d1.vd1.fa.ann")" = "7e63f266c4e6aa8fb3883141aec0a79403372c9855d0b3453d7055540798443d" ] \
        && [ -f "$bwa_dir/GRCh38.d1.vd1.fa.bwt" ] && [ "$(portable_sha256 "$bwa_dir/GRCh38.d1.vd1.fa.bwt")" = "b7c850d47b4f1bd02100cce58f5f3da4cbd0742761fd5749f2dbbad3101e7991" ] \
        && [ -f "$bwa_dir/GRCh38.d1.vd1.fa.pac" ] && [ "$(portable_sha256 "$bwa_dir/GRCh38.d1.vd1.fa.pac")" = "681f32f51da26b273ff60a07c86885c51b5bc432cb42ef02242a810b1ae8670a" ] \
        && [ -f "$bwa_dir/GRCh38.d1.vd1.fa.sa" ]  && [ "$(portable_sha256 "$bwa_dir/GRCh38.d1.vd1.fa.sa")"  = "381a1ef763a73fc8a3aa20f9777653492c648d6c20dd18c4704ec111bd88f217" ]
      ;;
    ref/giabv3/GRCh38_GIABv3.fasta.gz)
      [ "$r_sha256" = "b37cb9a9c0f8716fdb13e2c5c60f08cb8be40de39c3088176f1f62d9e8eeeed8" ] || return 1
      [ -f "${r_out%.gz}" ] && [ "$(portable_sha256 "${r_out%.gz}")" = "0713026dfdd61fb1cd68be4b9689a9e5c55c8708a3843248b0b574506fd0f0e9" ]
      ;;
    *)
      return 1
      ;;
  esac
}

want=" ${groups[*]} "
while IFS=$'\t' read -r -u 3 group name url md5 sha256 dest extra; do
  case "$group" in ''|\#*) continue ;; esac
  [[ "$want" == *" $group "* ]] || continue
  if [ -z "${dest:-}" ] || [ -n "${extra:-}" ]; then
    echo "fetch: malformed assets.tsv row for ${name:-<unknown>} (expected 6 columns)" >&2
    exit 1
  fi

  # Enforce fail-closed checksum policy BEFORE inspecting any completion markers!
  use_legacy_md5=0
  if [ -z "$sha256" ] || [ "$sha256" = "-" ]; then
    if [ "$allow_legacy_md5" -ne 1 ]; then
      echo "fetch: missing required SHA-256 checksum for $name (pass --allow-legacy-md5 to allow legacy MD5)" >&2
      exit 1
    fi
    if [[ ! "$md5" =~ ^[0-9a-f]{32}$ ]]; then
      echo "fetch: missing or invalid legacy MD5 checksum for $name" >&2
      exit 1
    fi
    use_legacy_md5=1
  elif [[ ! "$sha256" =~ ^[0-9a-f]{64}$ ]]; then
    echo "fetch: invalid SHA-256 checksum '$sha256' for $name" >&2
    exit 1
  fi

  out="$DATA/$dest"
  mkdir -p "$(dirname "$out")"

  if [ "$use_legacy_md5" -eq 0 ]; then
    # Always rehash raw cached files before accepting them; never trust sidecars alone.
    if [ -f "$out" ] && [ "$(portable_sha256 "$out")" = "$sha256" ]; then
      printf '%s\n' "$sha256" > "$out.sha256"
      rm -f "$out.legacy_md5"
      [ -f "$out.done" ] || touch "$out.done"
      record_manifest "$name" "$url" "$sha256" "upstream" "$out"
      echo "have: $name"
      continue
    fi

    if [ ! -f "$out" ] && { [ -f "$out.sha256" ] || [ -f "$out.done" ]; } && refs_transformed_verified "$group" "$dest" "$out" "$sha256"; then
      printf '%s\n' "$sha256" > "$out.sha256"
      rm -f "$out.legacy_md5"
      [ -f "$out.done" ] || touch "$out.done"
      echo "have: $name"
      continue
    fi

    rm -f "$out" "$out.done" "$out.sha256" "$out.legacy_md5"
    echo "fetch: $name"
    stage_verified_sha256 "$url" "$sha256" "$out"
    record_manifest "$name" "$url" "$sha256" "upstream" "$out"
    printf '%s\n' "$sha256" > "$out.sha256"
    touch "$out.done"
  else
    # Always rehash legacy-MD5 cached files before accepting them; never trust .legacy_md5 alone.
    if [ -f "$out" ] && [ "$(portable_md5 "$out")" = "$md5" ]; then
      printf '%s\n' "$md5" > "$out.legacy_md5"
      [ -f "$out.done" ] || touch "$out.done"
      record_manifest "$name" "$url" "$md5" "legacy-md5" "$out"
      echo "have: $name"
      continue
    fi

    rm -f "$out" "$out.done" "$out.sha256" "$out.legacy_md5" "$out.part"
    echo "fetch: $name"
    if ! curl -sSLf --retry 5 --retry-delay 5 -o "$out.part" "$url"; then
      rm -f "$out.part"
      echo "fetch: download failed for $name ($url)" >&2
      exit 1
    fi
    got="$(portable_md5 "$out.part")"
    if [ "$got" != "$md5" ]; then
      rm -f "$out.part"
      echo "fetch: MD5 MISMATCH $name got=$got want=$md5" >&2
      exit 1
    fi
    mv -f "$out.part" "$out"
    record_manifest "$name" "$url" "$got" "legacy-md5" "$out"
    printf '%s\n' "$md5" > "$out.legacy_md5"
    touch "$out.done"
  fi
done 3< "$ASSETS"
