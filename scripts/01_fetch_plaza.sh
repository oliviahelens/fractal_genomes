#!/usr/bin/env bash
# Download PLAZA gene families + species metadata for all instances used.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=data/raw/plaza; mkdir -p "$OUT/species"
F=https://ftp.psb.ugent.be/pub/plaza
get(){ [ -s "$2" ] || curl -sfL --retry 4 -m 900 -o "$2" "$1"; }
get $F/plaza_public_dicots_05/GeneFamilies/genefamily_data.HOMFAM.csv.gz   $OUT/dicots_HOMFAM.csv.gz
get $F/plaza_public_monocots_05/GeneFamilies/genefamily_data.HOMFAM.csv.gz $OUT/plaza_public_monocots_05_HOMFAM.csv.gz
get $F/plaza_basal_01/GeneFamilies/genefamily_data.HOMFAM.csv.gz           $OUT/plaza_basal_01_HOMFAM.csv.gz
get $F/plaza_gymno_01/GeneFamilies/genefamily_data.hom.csv.gz              $OUT/plaza_gymno_01_HOMFAM.csv.gz
for inst in plaza_public_dicots_05 plaza_public_monocots_05; do
  for sp in $(zcat $OUT/$( [ $inst = plaza_public_dicots_05 ] && echo dicots || echo $inst )_HOMFAM.csv.gz 2>/dev/null | grep -v '^#' | cut -f2 | sort -u); do
    get $F/$inst/SpeciesInformation/annot_sources.$sp.csv.gz $OUT/species/$sp.csv.gz || true
  done
done
