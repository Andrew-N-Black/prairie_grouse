#!/bin/bash
# =============================================================================
# accessory/z_hic_centromere.sh -- is the 41-44-Mb low-diversity segment of the Z the centromere?
# Hi-C contact map and centromere profiles (trans contacts, Rabl symmetry, insulation) along chr_Z of
# Z-bearing haplotypes, from the Hi-C alignments step 05 made to the final assemblies
# (qc_ptarmigan/<SAMPLE>.<HAP>.hic2final.cram). The rock ptarmigan segment (NC_064472.1:41.1-44.3 Mb) is
# placed on each assembly's chr_Z with minimap2 and marked.
#
#   sbatch --array=0-2 accessory/z_hic_centromere.sh            # from the repo root; one task per haplotype
#   HAPS="LEPC_F5541_hap2 ..." sbatch --array=0-N ...           # other haplotypes (males: both carry a Z)
# Output: ${OUT}/<SP>_<SAMPLE>_<HAP>.z_hic.png / .tsv
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -J z_hic
#SBATCH -o logs/%x_%A_%a.out
set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ASM="${ASM:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm}"
FINAL_DIR="${FINAL_DIR:-${ASM}/final_ptarmigan}"
QC_DIR="${QC_DIR:-${ASM}/qc_ptarmigan}"
PTARMIGAN="${PTARMIGAN:-${ASM}/ref/GCF_023343835.1_bLagMut1_primary_genomic.fna}"
SEGMENT="${SEGMENT:-NC_064472.1:41100000-44300000}"
OUT="${OUT:-${ASM}/z_hic_centromere}"
ENV="${ENV:-${CLUSTER_SCRATCH}/GROUSE/z_markers/env}"          # samtools, minimap2, python + numpy/matplotlib
HAPS="${HAPS:-LEPC_F5540_hap1 GRPC_F5597_hap1 STGR_F5457_hap1}"  # males (a Z in both haplotypes)
BIN="${BIN:-250000}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PATH="${ENV}/bin:${PATH}"
mkdir -p "$OUT"

read -r -a H <<< "$HAPS"
P="${H[${SLURM_ARRAY_TASK_ID:-0}]}"                 # e.g. LEPC_F5540_hap1
SAMPLE=$(cut -d_ -f2 <<< "$P"); HAP=$(cut -d_ -f3 <<< "$P")
FA="${FINAL_DIR}/${P}.pseudo_chr.fasta"
CRAM="${QC_DIR}/${SAMPLE}.${HAP}.hic2final.cram"
for f in "$FA" "$CRAM" "$PTARMIGAN"; do [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
[[ -s "${FA}.fai" ]] || samtools faidx "$FA"
[[ -s "${CRAM}.crai" ]] || samtools index -@ "$THREADS" "$CRAM"
ZLEN=$(awk '$1=="chr_Z"{print $2}' "${FA}.fai")
[[ -n "$ZLEN" ]] || { echo "ERROR: no chr_Z in $FA (W-bearing or Z-less haplotype?)"; exit 1; }

# place the ptarmigan low-diversity segment on this chr_Z (span of the hits on chr_Z, mapping quality >= 20)
[[ -s "${OUT}/segment.fa" ]] || samtools faidx "$PTARMIGAN" "$SEGMENT" > "${OUT}/segment.fa"
samtools faidx "$FA" chr_Z > "${OUT}/${P}.chrZ.fa"
MARK=$(minimap2 -x asm20 -t "$THREADS" "${OUT}/${P}.chrZ.fa" "${OUT}/segment.fa" 2>/dev/null \
       | awk '$12>=20 {if (s=="" || $8<s) s=$8; if ($9>e) e=$9} END{if (s!="") print s"-"e}')
rm -f "${OUT}/${P}.chrZ.fa"
echo ">>> ${P}: chr_Z ${ZLEN} bp; ptarmigan ${SEGMENT} -> chr_Z:${MARK:-not placed}"

samtools view -@ "$THREADS" -F 0x904 -f 64 -q 10 -T "$FA" "$CRAM" chr_Z \
    | python3 "${REPO_DIR}/accessory/z_hic_centromere.py" --length "$ZLEN" --bin "$BIN" \
        ${MARK:+--mark "$MARK"} --title "${P} chr_Z" --out "${OUT}/${P}.z_hic"
echo ">>> ${OUT}/${P}.z_hic.png"
