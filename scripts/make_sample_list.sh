#!/bin/bash
# Build samples_to_align.txt: one line per sample
#   animal_id   R1_path   R2_path
# Uses i7 and i5 from FASTQ filenames to look up animal in exp_design.txt.
# Run from project root (parent of data/ and exp_design.txt).

set -e
DATA_DIR="${DATA_DIR:-data}"
EXP_DESIGN="${EXP_DESIGN:-exp_design.txt}"
OUT="${OUT:-samples_to_align.txt}"

mkdir -p "$(dirname "$OUT")"

echo "animal	R1	R2" > "$OUT"

for R1 in "${DATA_DIR}"/*-R1.fastq.gz; do
  [ -f "$R1" ] || continue
  R2="${R1/-R1.fastq.gz/-R2.fastq.gz}"
  [ -f "$R2" ] || { echo "Missing R2 for $R1" >&2; continue; }
  # Filename like ...-P141-AATCGTTA-ACGTTATT-R1.fastq.gz -> i7=AATCGTTA, i5=ACGTTATT (last two fields)
  base=$(basename "$R1" -R1.fastq.gz)
  i7_i5=$(echo "$base" | awk -F'-' '{print $(NF-1)"\t"$NF}')
  i7=$(echo "$i7_i5" | cut -f1)
  i5=$(echo "$i7_i5" | cut -f2)
  animal=$(awk -v i7="$i7" -v i5="$i5" 'NR>1 && $12==i7 && $13==i5 {print $1; exit}' "$EXP_DESIGN")
  if [ -z "$animal" ]; then
    echo "No animal for i7=$i7 i5=$i5 ($R1)" >&2
    continue
  fi
  echo "${animal}	${R1}	${R2}" >> "$OUT"
done

echo "Wrote $(wc -l < "$OUT" | tr -d ' ') lines (including header) to $OUT"
