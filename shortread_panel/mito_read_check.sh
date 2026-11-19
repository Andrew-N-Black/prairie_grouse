#!/bin/bash
# =============================================================================
# 12f_mito_read_check.sh — which mitochondrial lineage do a bird's own reads carry?
#
# For birds whose assembled mitogenome is suspect (e.g. GRPC/STGR placed in the LEPC clade),
# count the LEPC vs GRPC/STGR alleles directly in the bird's reads at LEPC-diagnostic sites:
#   ~1.0 LEPC  -> the bird carries an LEPC mitogenome (real; compare with nuclear ancestry)
#   ~0.0 LEPC  -> its own reads are GRPC/STGR: the assembly came from someone else's reads
#   in between -> two mitogenomes in the reads (contamination / index hopping / mixed sample)
#
# 1. Diagnostic sites: the 23 HiFi/Hi-C mitogenomes (12c collect) are aligned to the GRPC
#    reference (MW574394.1); a site is diagnostic when all 7 LEPC share one base and every
#    GRPC + STGR (and the reference) another.
# 2. Reads: the bird's raw Illumina FASTQs (first MAX_PAIRS pairs; FASTQ_DIR/<ID>_S*_R1/R2_001.fastq.gz)
#    are mapped to MW574394.1 (minimap2 -x sr); without FASTQs, the unmapped reads of its panel CRAM
#    (no mitochondrial contig there) are used. Alleles are counted at the diagnostic sites
#    (mapping and base quality >= 20).
#
# USAGE (Gautschi; uses the step-14 conda env):
#   bash 12f_mito_read_check.sh F5594 F5436 F5440 F5596 F5547 F5462      # login node: a few birds
#   sbatch -A dewoody -p cpu -q standby -t 04:00:00 -c 16 --mem=16G -o logs/%x_%j.out \
#          -J mito_reads 12f_mito_read_check.sh F5594 F5436 F5440 F5596 F5547 F5462
# Output: ${OUT}/mito_read_check.tsv (one row per bird) and per-site counts per bird.
# =============================================================================
set -euo pipefail
unset LD_PRELOAD

ASM="${ASM:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm}"
MITO_DIR="${MITO_DIR:-${ASM}/mitogenomes}"
MREF="${MREF:-${MITO_DIR}/reference/ref.fasta}"                 # MW574394.1 (GRPC)
MITOS="${MITOS:-${MITO_DIR}/mitogenomes_pass.fasta}"            # 23 checked mitogenomes
ZM="${ZM:-${CLUSTER_SCRATCH}/GROUSE/z_markers}"
ENV="${ENV:-${ZM}/env}"
LEPC_REF="${LEPC_REF:-$(cat "${ZM}/ref/lepc_ref.path")}"       # reference the CRAMs were written with
CRAM_DIR="${CRAM_DIR:-/scratch/gautschi/blackan/GROUSE/nexus/crams_4.66x}"
FASTQ_DIR="${FASTQ_DIR:-${ASM}/raw/Shotgun}"                  # <ID>_S*_R1_001.fastq.gz / _R2_ (preferred)
MAX_PAIRS="${MAX_PAIRS:-10000000}"                            # read pairs used per bird (mito is high-copy)
OUT="${OUT:-${MITO_DIR}/read_check}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PATH="${ENV}/bin:${PATH}"

[[ $# -gt 0 ]] || { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
for f in "$MREF" "$MITOS" "$LEPC_REF"; do [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
mkdir -p "$OUT"

# ---- 1. diagnostic sites ----------------------------------------------------
if [[ ! -s "${OUT}/diagnostic_sites.tsv" ]]; then
    minimap2 -c -x asm5 --cs "$MREF" "$MITOS" 2>/dev/null > "${OUT}/mitos_vs_ref.paf"
    python3 - "${OUT}/mitos_vs_ref.paf" "$MREF" "${OUT}/diagnostic_sites.tsv" <<'PY'
import re, sys
paf, ref_fa, out = sys.argv[1:4]
ref = "".join(l.strip() for l in open(ref_fa) if not l.startswith(">")).upper()
calls, cov = {}, {}
for line in open(paf):
    t = line.rstrip("\n").split("\t")
    q = t[0]
    cs = [x[5:] for x in t[12:] if x.startswith("cs:Z:")][0]
    pos = int(t[7])
    calls.setdefault(q, {}); cov.setdefault(q, [])
    s0 = pos
    for op, val in re.findall(r"([:*+\-~])([0-9a-z]+)", cs):
        if op == ":":
            pos += int(val)
        elif op == "*":
            calls[q][pos] = val[1].upper(); pos += 1
        elif op == "-":
            for k in range(len(val)):
                calls[q][pos + k] = "-"
            pos += len(val)
    cov[q].append((s0, pos))
def base(q, p):
    if p in calls[q]:
        return calls[q][p]
    return ref[p] if any(s <= p < e for s, e in cov[q]) else None
sp = {q: ("LEPC" if q.startswith("LEPC") else "OTHER") for q in calls}
n = {"LEPC": sum(v == "LEPC" for v in sp.values()), "OTHER": sum(v == "OTHER" for v in sp.values())}
rows = []
for p in range(len(ref)):
    b = {q: base(q, p) for q in calls}
    if any(v is None or v == "-" for v in b.values()):
        continue
    L = {b[q] for q in calls if sp[q] == "LEPC"}
    O = {b[q] for q in calls if sp[q] == "OTHER"} | {ref[p]}
    if len(L) == 1 and len(O) == 1 and L != O:
        rows.append((p + 1, ref[p], L.pop(), O.pop()))
with open(out, "w") as fh:
    fh.write("pos\tref\tLEPC\tother\n")
    for r in rows:
        fh.write("%d\t%s\t%s\t%s\n" % r)
sys.stderr.write("[sites] %d mitogenomes (%d LEPC, %d GRPC/STGR); %d LEPC-diagnostic sites\n"
                 % (len(calls), n["LEPC"], n["OTHER"], len(rows)))
PY
fi
awk 'NR>1{print "MW_REF\t"$1-1"\t"$1}' "${OUT}/diagnostic_sites.tsv" > "${OUT}/diag.bed"
MCONTIG=$(awk '/^>/{sub(/^>/,""); print $1; exit}' "$MREF")
sed -i "s/^MW_REF/${MCONTIG}/" "${OUT}/diag.bed"
[[ -s "${MREF}.fai" ]] || samtools faidx "$MREF"

# ---- 2. per bird --------------------------------------------------------------
echo -e "ID\treads_used\tmito_reads\tdiag_sites_covered\tmedian_depth\treads_LEPC\treads_other\tfrac_LEPC\tsites_mixed" > "${OUT}/mito_read_check.tsv"
for ID in "$@"; do
    B="${OUT}/${ID}.mito.bam"
    R1=$(ls "${FASTQ_DIR}/${ID}"_S*_R1_001.fastq.gz 2>/dev/null | head -1 || true)
    if [[ ! -s "${B}.bai" ]]; then
        if [[ -n "$R1" ]]; then                      # raw reads: first MAX_PAIRS pairs, paired mapping
            R2="${R1/_R1_001/_R2_001}"
            minimap2 -ax sr -t "$THREADS" "$MREF" \
                <(zcat "$R1" | head -n $((MAX_PAIRS * 4))) <(zcat "$R2" | head -n $((MAX_PAIRS * 4))) 2>/dev/null \
                | samtools view -b -F 4 -q 20 - | samtools sort -o "$B" - && samtools index "$B"
            echo "${MAX_PAIRS}_pairs_FASTQ" > "${OUT}/${ID}.n_unmapped"
        else                                         # fallback: unmapped reads kept in the panel CRAM
            C="${CRAM_DIR}/${ID}.cram"
            [[ -s "$C" ]] || { echo "skip $ID: no FASTQ in $FASTQ_DIR and no $C"; continue; }
            samtools fastq -f 4 -T "$LEPC_REF" "$C" 2>/dev/null > "${OUT}/${ID}.unmapped.fq"
            minimap2 -ax sr -t "$THREADS" "$MREF" "${OUT}/${ID}.unmapped.fq" 2>/dev/null \
                | samtools view -b -F 4 -q 20 - | samtools sort -o "$B" - && samtools index "$B"
            echo "$(awk 'END{print NR/4}' "${OUT}/${ID}.unmapped.fq")_unmapped_CRAM_reads" > "${OUT}/${ID}.n_unmapped"
            rm -f "${OUT}/${ID}.unmapped.fq"
        fi
    fi
    samtools mpileup -l "${OUT}/diag.bed" -q 20 -Q 20 -B -f "$MREF" "$B" 2>/dev/null > "${OUT}/${ID}.diag.pileup"
    python3 - "${OUT}/diagnostic_sites.tsv" "${OUT}/${ID}.diag.pileup" "$ID" "$(cat "${OUT}/${ID}.n_unmapped")" \
        "$(samtools view -c "$B")" >> "${OUT}/mito_read_check.tsv" <<'PY'
import re, statistics, sys
sites, pile, bird, n_un, n_mito = sys.argv[1:6]
d = {}
for l in list(open(sites))[1:]:
    p, r, L, O = l.split()
    d[int(p)] = (r, L, O)
def counts(bases, ref):
    out = {}
    bases = re.sub(r"\^.", "", bases).replace("$", "")
    i = 0
    while i < len(bases):
        ch = bases[i]
        if ch in "+-":
            m = re.match(r"[+-](\d+)", bases[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        b = ref if ch in ".," else ch.upper()
        out[b] = out.get(b, 0) + 1; i += 1
    return out
nl = no = mixed = 0; deps = []
for l in open(pile):
    t = l.split("\t")
    p = int(t[1])
    if p not in d:
        continue
    r, L, O = d[p]
    c = counts(t[4], t[2].upper())
    a, b = c.get(L, 0), c.get(O, 0)
    if a + b == 0:
        continue
    nl += a; no += b; deps.append(a + b)
    if min(a, b) / (a + b) > 0.1:
        mixed += 1
tot = nl + no
print("\t".join(map(str, [bird, n_un, n_mito, len(deps), int(statistics.median(deps)) if deps else 0,
                          nl, no, "%.3f" % (nl / tot) if tot else "NA", mixed])))
PY
    tail -1 "${OUT}/mito_read_check.tsv"
done
echo "table: ${OUT}/mito_read_check.tsv"
