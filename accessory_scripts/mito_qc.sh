#!/bin/bash
# =============================================================================
# accessory_scripts/mito_qc.sh — QA/QC of the MitoHiFi mitogenomes
# (run after accessory_scripts/mitohifi_array.sh). Login node is fine: the data
# are tiny (one ~16.7 kb sequence and a few thousand reads per bird).
#
#   bash accessory_scripts/mito_qc.sh qc            # table for all birds
#   bash accessory_scripts/mito_qc.sh trim SAMPLE   # fix one over-long assembly
#   bash accessory_scripts/mito_qc.sh collect       # one multi-FASTA of passing birds
#
# qc writes mitogenomes/mito_qc.tsv, one row per bird:
#   length_bp, ref_len_bp, len_vs_ref    should be ~1.00 (control-region indels
#                                        allow a few hundred bp)
#   circular                             MitoHiFi's circularisation flag for the
#                                        chosen contig
#   CDS / rRNA / tRNA                    expect 13 / 2 / 22
#   first_feature                        expect trnF (tRNA-Phe): standard rotation
#   pct_id_ref, ref_cov                  identity to the reference and the fraction
#                                        of it covered (expect cov >= 0.98)
#   depth_min / depth_med                HiFi read depth along the mitogenome; a
#                                        window near 0 means a misassembly or a
#                                        NUMT-derived segment
#   PASS                                 all of the above within limits
# It also writes mitogenomes/mito_pairwise_identity.tsv: all-vs-all identity.
# Use it to check the clustering. Within species expect >99%; a bird
# clustering with the "wrong" species is either a sample swap or
# mitochondrial introgression, and needs following up either way.
# =============================================================================
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
MITO_DIR="${PROJECT_DIR}/mitogenomes"
REF_DIR="${MITO_DIR}/reference"
SIF="${PROJECT_DIR}/containers/mitohifi.sif"
APPTAINER_BIN="$(command -v apptainer || command -v singularity)"
cx() { "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$SIF" "$@"; }

REF_LEN=$(grep -v '>' "${REF_DIR}/ref.fasta" | tr -d '\n' | wc -c)
seqlen() { grep -v '>' "$1" | tr -d '\n\r' | wc -c; }

case "${1:-qc}" in
# ------------------------------------------------------------------------------
qc)
  OUT="${MITO_DIR}/mito_qc.tsv"
  printf "sample\tlength_bp\tref_len_bp\tlen_vs_ref\tcircular\tCDS\trRNA\ttRNA\tfirst_feature\tpct_id_ref\tref_cov\tdepth_min\tdepth_med\tPASS\n" > "$OUT"
  for D in "$MITO_DIR"/*_F*/; do
    S=$(basename "$D"); F="${D}final_mitogenome.fasta"; GB="${D}final_mitogenome.gb"
    [[ -s "$F" ]] || { printf "%s\tMISSING\n" "$S" >> "$OUT"; continue; }
    L=$(seqlen "$F")
    RATIO=$(awk -v a="$L" -v b="$REF_LEN" 'BEGIN{printf "%.3f", a/b}')
    # circular flag of the chosen contig (MitoHiFi contigs_stats.tsv)
    CIRC=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if(tolower($i) ~ /circular/) c=i; next}
                       c && NR==2 {print $c; exit}' "${D}contigs_stats.tsv" 2>/dev/null || true)
    CIRC=${CIRC:-NA}
    # feature counts and first feature from the GenBank file
    CDS=$(grep -cE '^     CDS ' "$GB" 2>/dev/null || echo 0)
    RRNA=$(grep -cE '^     rRNA ' "$GB" 2>/dev/null || echo 0)
    TRNA=$(grep -cE '^     tRNA ' "$GB" 2>/dev/null || echo 0)
    FIRST=$(awk '/^     (gene|tRNA|rRNA|CDS) /{f=1} f && /\/(gene|product)=/{gsub(/.*="|"$/,""); print; exit}' "$GB" 2>/dev/null || echo NA)
    # identity / coverage against the reference
    read -r PID COV < <(cx minimap2 -x asm10 -c --secondary=no "${REF_DIR}/ref.fasta" "$F" 2>/dev/null \
        | awk -v rl="$REF_LEN" '{m+=$10; a+=$11; c+=$9-$8} END{if(a) printf "%.2f %.3f\n", 100*m/a, c/rl; else print "NA NA"}')
    # read depth along the assembly (the pre-filtered mito reads)
    READS="${D}mito_candidate_reads.fa"; DMIN=NA; DMED=NA
    if [[ -s "$READS" ]]; then
      read -r DMIN DMED < <(cx bash -c "minimap2 -t 4 -ax map-hifi --secondary=no '$F' '$READS' 2>/dev/null \
          | samtools sort -o - 2>/dev/null | samtools depth -a - " \
          | awk '{d[NR]=$3} END{n=asort(d); printf "%d %d\n", d[1], d[int((n+1)/2)]}')
    fi
    PASS=$(awk -v r="$RATIO" -v c="$CDS" -v rr="$RRNA" -v t="$TRNA" -v cov="$COV" -v dm="$DMIN" -v f="$FIRST" \
      'BEGIN{ok = (r>0.97 && r<1.05 && c==13 && rr==2 && t==22 && cov>=0.98 && (dm=="NA" || dm>=5));
             print (ok ? "PASS" : "CHECK")}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$S" "$L" "$REF_LEN" "$RATIO" "$CIRC" "$CDS" "$RRNA" "$TRNA" "$FIRST" "$PID" "$COV" "$DMIN" "$DMED" "$PASS" >> "$OUT"
  done
  column -t -s$'\t' "$OUT" 2>/dev/null || cat "$OUT"

  # all-vs-all identity matrix (long format)
  ALL="${MITO_DIR}/_all_for_pairwise.fa"; > "$ALL"
  for F in "$MITO_DIR"/*_F*/final_mitogenome.fasta; do
    S=$(basename "$(dirname "$F")"); echo ">$S" >> "$ALL"; grep -v '>' "$F" >> "$ALL"
  done
  cx minimap2 -x asm10 -c -X "$ALL" "$ALL" 2>/dev/null \
    | awk 'BEGIN{OFS="\t"; print "a","b","pct_identity","aligned_bp"} {print $1,$6,sprintf("%.2f",100*$10/$11),$11}' \
    > "${MITO_DIR}/mito_pairwise_identity.tsv"
  echo; echo "QC table: ${OUT}"; echo "Pairwise identity: ${MITO_DIR}/mito_pairwise_identity.tsv"
  ;;
# ------------------------------------------------------------------------------
trim)
  # Over-long assembly = more than one copy of the circle in a row. Find where the
  # reference start (tRNA-Phe) and end land in the assembly, cut out exactly one
  # copy, then let MitoHiFi re-annotate it in contig mode (-c).
  S="${2:?usage: $0 trim SAMPLE   e.g. STGR_F5480}"
  D="${MITO_DIR}/${S}"; F="${D}/final_mitogenome.fasta"
  [[ -s "$F" ]] || { echo "ERROR: ${F} not found"; exit 1; }
  W="${D}/trimmed"; mkdir -p "$W"; cd "$W"
  cp -n "$F" "${D}/final_mitogenome.untrimmed.fasta"
  # best single hit of the whole reference on the assembly
  HIT=$(cx minimap2 -x asm10 -c --secondary=no "$F" "${REF_DIR}/ref.fasta" 2>/dev/null \
        | awk -v rl="$REF_LEN" '($4-$3)/rl >= 0.98 {print $8, $9, $5, ($4-$3)/rl; exit}')
  [[ -n "$HIT" ]] || { echo "No single hit covers >=98% of the reference — inspect ${D}/contigs_stats.tsv and potential_contigs/ instead"; exit 1; }
  read -r TS TE STRAND QCOV <<< "$HIT"
  echo ">>> one full copy at ${TS}-${TE} (${STRAND}; covers $(awk -v q="$QCOV" 'BEGIN{printf "%.1f", 100*q}')% of reference)"
  NAME=$(head -1 "$F" | sed 's/^>//; s/ .*//')
  cx samtools faidx "$F"
  cx samtools faidx $( [[ "$STRAND" == "-" ]] && echo "-i" ) "$F" "${NAME}:$((TS+1))-${TE}" \
     | sed "1s/.*/>${S}_mito_trimmed/" > one_copy.fa
  echo ">>> trimmed length: $(seqlen one_copy.fa) bp (reference ${REF_LEN})"
  cx mitohifi.py -c one_copy.fa -f "${REF_DIR}/ref.fasta" -g "${REF_DIR}/ref.gb" -t 4 -o 2 -a animal
  # promote the trimmed result (untrimmed copy kept alongside)
  for X in fasta gb gff; do
    [[ -s "final_mitogenome.${X}" ]] && cp "final_mitogenome.${X}" "${D}/final_mitogenome.${X}"
  done
  echo ">>> ${S}: final_mitogenome.* replaced; rerun '$0 qc'"
  ;;
# ------------------------------------------------------------------------------
collect)
  # One multi-FASTA of the PASS birds, for NCBI / phylogenetics. Not added to
  # the nuclear assemblies or the pangenome.
  OUT="${MITO_DIR}/mitogenomes_pass.fasta"; > "$OUT"
  awk -F'\t' 'NR>1 && $NF=="PASS"{print $1}' "${MITO_DIR}/mito_qc.tsv" | while read -r S; do
    echo ">${S}_mitogenome" >> "$OUT"; grep -v '>' "${MITO_DIR}/${S}/final_mitogenome.fasta" >> "$OUT"
  done
  echo "$(grep -c '>' "$OUT") mitogenomes -> ${OUT}"
  ;;
*) echo "usage: $0 {qc|trim SAMPLE|collect}"; exit 1 ;;
esac
