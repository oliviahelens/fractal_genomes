#!/usr/bin/env bash
# Download one protein sequence per gene (PLAZA "selected transcript") for every species in the matrix.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT=data/raw/proteomes; mkdir -p $OUT
F=https://ftp.psb.ugent.be/pub/plaza
get(){ [ -s "$2" ] || curl -sfL --retry 4 -m 1800 -o "$2" "$1" || echo "FAILED $1"; }
tail -n +2 data/processed/species_instance.tsv | while IFS=$'\t' read sp inst; do
  case $inst in
    dicots)   get $F/plaza_public_dicots_05/Fasta/proteome.selected_transcript.$sp.fasta.gz   $OUT/$sp.fa.gz ;;
    monocots) get $F/plaza_public_monocots_05/Fasta/proteome.selected_transcript.$sp.fasta.gz $OUT/$sp.fa.gz ;;
    basal)    get $F/plaza_basal_01/Fasta/proteome.selected_transcript.$sp.fasta.gz           $OUT/$sp.fa.gz ;;
    gymno)    get $F/plaza_gymno_01/Fasta/proteome.$sp.csv.gz                                  $OUT/$sp.gymno.csv.gz ;;
  esac
done
