#!/bin/bash
# =============================================================================
# 12d_mito_tree.sh — maximum-likelihood tree of the mitogenomes (report Figure 17C)
#
# Run after 12c_mito_qc.sh qc + collect. Small job: login node or a short standby job.
#   bash 12d_mito_tree.sh            # fetch outgroup, align, infer the tree
#   Rscript 12e_plot_mito_tree_ggtree.R $CLUSTER_SCRATCH/GROUSE/grouse_asm/mitogenomes/tree/mito_tree.contree
#
# Inputs : mitogenomes/mitogenomes_pass.fasta   (12c_mito_qc.sh collect; the 23 birds)
#          mitogenomes/reference/ref.fasta      (the GRPC reference MW574394.1, from
#                                                12_mitohifi_array.sh setup)
#          rock ptarmigan mitogenome PV239887.1 (outgroup; fetched from NCBI once)
# Method (as in report Appendix S1): MAFFT v7 FFT-NS-i (--maxiterate 1000);
#          IQ-TREE 3 with ModelFinder (-m MFP) and 1,000 ultrafast bootstraps,
#          rooted on the ptarmigan (-o).
# Output : mitogenomes/tree/mito_aln.fasta, mito_tree.{treefile,contree,iqtree,log}
#
# NOTE: this script was written from the methods text after the analysis; the tree
# in report v80 was inferred by hand with these settings. Re-running it should
# reproduce the report tree (identical topology; supports within UFBoot noise).
# =============================================================================
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
MITO_DIR="${PROJECT_DIR}/mitogenomes"
OUT_DIR="${MITO_DIR}/tree"
OUTGROUP_ACC="${OUTGROUP_ACC:-PV239887.1}"         # rock ptarmigan, Lagopus muta
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
SEED="${SEED:-12345}"

module load biocontainers mafft >/dev/null 2>&1 || true
module load biocontainers iqtree >/dev/null 2>&1 || true
MAFFT_BIN="${MAFFT_BIN:-$(command -v mafft || true)}"
IQTREE_BIN="${IQTREE_BIN:-$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)}"
[[ -n "$MAFFT_BIN" ]] || { echo "ERROR: mafft not found (set MAFFT_BIN)"; exit 1; }
[[ -n "$IQTREE_BIN" ]] || { echo "ERROR: iqtree not found (set IQTREE_BIN; IQ-TREE 3 was used for the report)"; exit 1; }

PASS="${MITO_DIR}/mitogenomes_pass.fasta"
REF="${MITO_DIR}/reference/ref.fasta"
[[ -s "$PASS" ]] || { echo "ERROR: $PASS missing; run: bash 12c_mito_qc.sh collect"; exit 1; }
[[ -s "$REF" ]]  || { echo "ERROR: $REF missing; run: bash 12_mitohifi_array.sh setup"; exit 1; }
mkdir -p "$OUT_DIR"

OG="${OUT_DIR}/${OUTGROUP_ACC}.fasta"
if [[ ! -s "$OG" ]]; then
    echo ">>> fetching outgroup ${OUTGROUP_ACC} from NCBI"
    curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${OUTGROUP_ACC}&rettype=fasta&retmode=text" \
        | awk -v id="$OUTGROUP_ACC" 'NR==1 && /^>/ {print ">" id; next} {print}' > "${OG}.part"
    grep -q '^>' "${OG}.part" || { echo "ERROR: NCBI fetch failed (run on a login node)"; rm -f "${OG}.part"; exit 1; }
    mv "${OG}.part" "$OG"
fi

# one record each; headers trimmed to the accession / sample name
IN="${OUT_DIR}/mito_input.fasta"
{
    cat "$PASS"
    awk 'NR==1 && /^>/ {split(substr($0,2), a, " "); print ">" a[1]; next} {print}' "$REF"
    cat "$OG"
} > "$IN"
echo ">>> $(grep -c '^>' "$IN") sequences (23 birds + GRPC reference + outgroup expected)"

echo ">>> MAFFT FFT-NS-i"
"$MAFFT_BIN" --thread "$THREADS" --retree 2 --maxiterate 1000 "$IN" > "${OUT_DIR}/mito_aln.fasta"

echo ">>> IQ-TREE ($("$IQTREE_BIN" --version 2>&1 | grep -m1 -io 'version [0-9.]*'))"
"$IQTREE_BIN" -s "${OUT_DIR}/mito_aln.fasta" -m MFP -B 1000 -o "$OUTGROUP_ACC" \
    -T "$THREADS" --seed "$SEED" --prefix "${OUT_DIR}/mito_tree" -redo --quiet

echo ">>> best model: $(grep -m1 'Best-fit model' "${OUT_DIR}/mito_tree.iqtree" 2>/dev/null || grep -m1 'Best-fit model' "${OUT_DIR}/mito_tree.log")"
echo ">>> tree: ${OUT_DIR}/mito_tree.contree"
echo "    plot: Rscript 12e_plot_mito_tree_ggtree.R ${OUT_DIR}/mito_tree.contree"
