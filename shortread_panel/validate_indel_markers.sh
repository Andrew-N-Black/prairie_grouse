#!/bin/bash
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -J z_indel_validate
#SBATCH -o logs/14h_validate_%x_%A_%a.log
# =============================================================================
# validate_indel_markers.sh -- check the pilot Z-indel markers (step 22b)
#   (a) in higher-coverage short-read data (raw/Shotgun): each bird's reads are
#       aligned to the marker alleles only, and reads spanning each indel
#       junction are counted, giving a genotype per bird and marker. Females
#       (one Z) should show one allele; a second allele in a female points to a
#       W copy, a duplication or contamination.
#   (b) against other references: every allele is aligned to the rock ptarmigan
#       (female; has Z and W) and chicken genomes and to the grouse female W
#       sequence, reporting the best-hit chromosome and any hit to the W.
#
# Usage (from the folder holding 14h_*.sh and 14h_*.py; logs/ must exist):
#   bash validate_indel_markers.sh submit     # prep -> per-bird array -> summary
# or step by step:
#   sbatch validate_indel_markers.sh prep
#   sbatch --array=1-N validate_indel_markers.sh map      (N = birds in samples.tsv)
#   bash   validate_indel_markers.sh summarize
# Outputs ($OUT):
#   pilot.genotypes.tsv   one row per bird x marker: read counts per allele, call, status
#   pilot.summary.tsv     per marker x species: n expected / unexpected / two alleles
#   pilot.refcheck.tsv    per allele x reference: best-hit chromosome, W hits
# =============================================================================
set -euo pipefail
MODE=${1:-submit}

# ---- edit -------------------------------------------------------------------
BASE=$CLUSTER_SCRATCH/GROUSE/prairie_grouse
STEP14G=$BASE/14g_z_indel_markers                 # output folder of step 22b
OUT=$BASE/14h_validate_indel_markers
FASTQ_DIR=$CLUSTER_SCRATCH/GROUSE/grouse_asm/raw/Shotgun
# sample <tab> species <tab> sex (M/F/NA), header line; birds not listed get NA
META=$BASE/config/shotgun_metadata.tsv
REF=$CLUSTER_SCRATCH/GROUSE/grouse_asm/ref
PTARMIGAN=$REF/GCF_023343835.1_bLagMut1_primary_genomic.fna
PTARMIGAN_MAP=$REF/GCF_023343835.1_bLagMut1_primary_assembly_report.txt
CHICKEN=$REF/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna
CHICKEN_MAP=$REF/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_assembly_report.txt
PILOT="NW_026294758.1_17651301 NW_026294758.1_30763712 NW_026294758.1_29441058
NW_026294758.1_41082646 NW_026294813.1_15223534 NW_026294758.1_15609875
NW_026294758.1_16879095 NW_026294813.1_3493004 NW_026294758.1_12563775
NW_026294758.1_32464091"
CONDA_ENV=zindel
THREADS=${SLURM_NTASKS:-8}
# -----------------------------------------------------------------------------

if [[ -n $CONDA_ENV ]]; then
  module load conda 2>/dev/null || module load anaconda 2>/dev/null || true
  eval "$(conda shell.bash hook)"; set +u; conda activate "$CONDA_ENV"; set -u
fi
for t in minimap2 samtools python3; do command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
PY=$HERE/validate_indel_markers.py; [[ -f $PY ]] || PY=$HERE/shortread_panel/validate_indel_markers.py
[[ -f $PY ]] || { echo "ERROR: $PY not found; run from the folder holding both scripts" >&2; exit 1; }
mkdir -p "$OUT"/{counts,bam,paf} "$HERE/logs"

case $MODE in
# ---- submit: chain prep -> array -> summarize ---------------------------------
submit)
  # sample sheet first, so the array size is known
  bash "$0" samples
  N=$(($(wc -l < "$OUT/samples.tsv") - 1))
  j1=$(sbatch --parsable -q standby "$0" prep)
  j2=$(sbatch --parsable -q standby --array=1-$N -d afterok:$j1 "$0" map)
  j3=$(sbatch --parsable -q standby -n 1 --mem=4G -t 0:30:00 -d afterok:$j2 "$0" summarize)
  echo "submitted: prep $j1, map array $j2 (1-$N), summarize $j3"
  ;;

# ---- samples: build the sample sheet from FASTQ names + metadata ------------
samples)
  python3 - "$FASTQ_DIR" "$META" "$OUT/samples.tsv" <<'EOF'
import sys, os, re, glob, csv
fq, meta, out = sys.argv[1:]
M = {}
if os.path.exists(meta):
    for r in csv.DictReader(open(meta), delimiter="\t"):
        M[r["sample"]] = r
rows = []
for r1 in sorted(glob.glob(os.path.join(fq, "*_R1_001.fastq.gz"))):
    b = os.path.basename(r1)
    s = re.sub(r"_(S\d+|SRR\d+)_R1_001\.fastq\.gz$", "", b)
    r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    if not os.path.exists(r2):
        continue
    m = M.get(s, {})
    rows.append([s, m.get("species", "NA") or "NA", m.get("sex", "NA") or "NA", r1, r2])
with open(out, "w") as f:
    f.write("sample\tspecies\tsex\tR1\tR2\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
na = [r[0] for r in rows if r[1] == "NA"]
print("samples: %d birds; %d without species in %s%s" % (len(rows), len(na), meta,
      (": " + " ".join(na)) if na else ""))
EOF
  ;;

# ---- prep: alleles + junctions, reference check -----------------------------
prep)
  cd "$OUT"
  printf "%s\n" $PILOT > pilot_markers.txt
  python3 "$PY" alleles --markers pilot_markers.txt --haps "$STEP14G/haplotypes.tsv" \
    --pcr-dir "$STEP14G/pcr" --out-prefix pilot
  samtools faidx pilot.alleles.fa
  echo "alleles: $(grep -c '>' pilot.alleles.fa) for $(wc -l < pilot_markers.txt) markers"
  # whole allele sequences (product + 400-bp flanks) against other references;
  # two presets: sensitive short-seed (divergent copies) and asm20 (long alleles)
  for spec in "ptarmigan:$PTARMIGAN" "chicken:$CHICKEN" "grouseW:$STEP14G/pcr/W_all.fa"; do
    lab=${spec%%:*}; fa=${spec#*:}
    [[ -s $fa ]] || { echo "WARN: $fa missing, skipping $lab" >&2; continue; }
    { minimap2 -c -x sr -k 11 -w 3 -N 50 -t "$THREADS" "$fa" pilot.alleles.fa 2>/dev/null
      minimap2 -c -x asm20 -N 50 -t "$THREADS" "$fa" pilot.alleles.fa 2>/dev/null; } > paf/$lab.paf
  done
  python3 "$PY" refcheck --alleles-fa pilot.alleles.fa --out pilot.refcheck.tsv --refs \
    "ptarmigan=paf/ptarmigan.paf=$PTARMIGAN_MAP" "chicken=paf/chicken.paf=$CHICKEN_MAP" "grouseW=paf/grouseW.paf="
  ;;

# ---- map: one bird per array task ----------------------------------------------
map)
  cd "$OUT"
  i=${SLURM_ARRAY_TASK_ID:?run as an array}
  read -r s sp sex r1 r2 < <(sed -n "$((i + 1))p" samples.tsv)
  # reads aligned to the marker alleles only; unmapped reads are dropped
  minimap2 -ax sr -t "$THREADS" pilot.alleles.fa "$r1" "$r2" 2> logs_map_$s.txt \
    | samtools view -b -F 4 -o bam/$s.unsorted.bam -
  samtools sort -o bam/$s.bam bam/$s.unsorted.bam && rm bam/$s.unsorted.bam && samtools index bam/$s.bam
  samtools view bam/$s.bam \
    | python3 "$PY" genotype --junctions pilot.junctions.tsv --sample "$s" --out counts/$s.counts.tsv
  rm -f logs_map_$s.txt
  echo "$s ($sp, $sex): $(awk 'NR>1{n+=$3} END{print n+0}' counts/$s.counts.tsv) junction-spanning reads"
  ;;

# ---- summarize -------------------------------------------------------------------
summarize)
  cd "$OUT"
  python3 "$PY" summarize --samples samples.tsv --alleles pilot.alleles.tsv --counts-dir counts --out-prefix pilot
  echo; echo "== per marker x species =="; (column -t pilot.summary.tsv 2>/dev/null || cat pilot.summary.tsv)
  echo; echo "== reference check (allele, reference, best chromosome, W hits, verdict) =="
  cut -f1,2,3,10,11 pilot.refcheck.tsv | (column -t 2>/dev/null || cat)
  ;;
*) echo "unknown mode $MODE (submit | samples | prep | map | summarize)" >&2; exit 1 ;;
esac
