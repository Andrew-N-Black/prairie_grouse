#!/usr/bin/env python3
"""
validate_indel_markers.py -- check pilot Z-indel markers (step 22b) in
higher-coverage short-read data, and check which chromosome each marker
matches in other references (Z vs W).

Sub-commands (called by validate_indel_markers.sh):

  alleles    for each marker, cut the in-silico PCR product (+ flanks) from
             the Z haplotypes, keep one sequence per band-size class, and find
             the indel breakpoints ("junctions") in each allele
  genotype   count reads that span each allele's junctions in one sample
             (SAM on stdin, from minimap2 -ax sr against alleles.fa)
  summarize  combine per-sample counts into genotype calls and a per-marker,
             per-species summary
  refcheck   summarize alignments of each allele to other references
             (rock ptarmigan, chicken): best-hit chromosome, any hit to the W

A read supports an allele only if it aligns across a junction with
>= --anchor bp on both sides and few mismatches. Mobile-element insertions
attract reads from every copy of the element in the genome; requiring the
junction between the insertion and its unique flank removes those.

Standard library only; uses samtools for sequence extraction.
"""
import argparse, collections, csv, difflib, glob, os, re, subprocess, sys

SPECIES = ["LEPC", "GRPC", "STGR"]


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def faidx(fa, region):
    out = subprocess.run(["samtools", "faidx", fa, region], capture_output=True, text=True, check=True).stdout
    return "".join(out.splitlines()[1:]).upper()


def contig_len(fa, ctg):
    fai = fa + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", fa], check=True)
    for line in open(fai):
        t = line.split("\t")
        if t[0] == ctg:
            return int(t[1])
    raise KeyError(ctg)


# ----------------------------------------------------------------- alleles
def cmd_alleles(a):
    markers = [l.split()[0] for l in open(a.markers) if l.strip() and not l.startswith("#")]
    haps = list(csv.DictReader(open(a.haps), delimiter="\t"))
    prods = collections.defaultdict(list)       # marker -> [(species, hap, ctg, s, e, strand, fasta)]
    for h in haps:
        bed = os.path.join(a.pcr_dir, "Z.%s.bed" % h["hap_id"])
        if not os.path.exists(bed):
            continue
        hits = collections.defaultdict(list)
        for line in open(bed):
            t = line.rstrip("\n").split("\t")
            if t[3] in markers and int(t[2]) - int(t[1]) <= a.max_len:
                hits[t[3]].append(t)
        for m, ts in hits.items():
            if len(ts) == 1:                    # single product only
                t = ts[0]
                prods[m].append((h["species"], h["hap_id"], t[0], int(t[1]), int(t[2]), t[5], h["z_fasta"]))
    fa_out = open(a.out_prefix + ".alleles.fa", "w")
    jt = open(a.out_prefix + ".junctions.tsv", "w")
    jt.write("allele\tmarker\tsize\tspecies\tjunction\n")
    info = open(a.out_prefix + ".alleles.tsv", "w")
    info.write("marker\tallele\tsize\tspecies\tn_haps\trepresentative\n")
    for m in markers:
        P = prods.get(m, [])
        if not P:
            sys.stderr.write("WARN: no products for %s\n" % m)
            continue
        # cluster product lengths into size classes (<= --size-tol bp apart)
        P.sort(key=lambda p: p[4] - p[3])
        classes = []
        for p in P:
            L = p[4] - p[3]
            if classes and L - classes[-1][-1][4] + classes[-1][-1][3] <= a.size_tol:
                classes[-1].append(p)
            else:
                classes.append([p])
        alleles = []
        for k, cl in enumerate(classes):
            lens = sorted(p[4] - p[3] for p in cl)
            med = lens[len(lens) // 2]
            rep = min(cl, key=lambda p: abs((p[4] - p[3]) - med))
            sp, hap, ctg, s, e, strand, fa = rep
            clen = contig_len(fa, ctg)
            fs, fe = max(0, s - a.flank), min(clen, e + a.flank)
            seq = faidx(fa, "%s:%d-%d" % (ctg, fs + 1, fe))
            lf, rf = s - fs, fe - e
            if strand == "-":
                seq, lf, rf = rc(seq), rf, lf
            spp = sorted(set(p[0] for p in cl), key=SPECIES.index)
            name = "%s__a%d_%dbp_%s" % (m, k + 1, med, "-".join(spp))
            counts = collections.Counter(p[0] for p in cl)
            alleles.append(dict(name=name, seq=seq, lf=lf, prod=seq[lf:len(seq) - rf], size=med, spp=spp))
            info.write("%s\t%s\t%d\t%s\t%s\t%s\n" % (m, name, med, ",".join("%s:%d" % (x, counts[x]) for x in spp),
                                                     len(cl), hap))
        # junctions: indel breakpoints of every allele against every other
        junc = collections.defaultdict(set)
        for i in range(len(alleles)):
            for j in range(i + 1, len(alleles)):
                A, B = alleles[i], alleles[j]           # A shorter or equal
                if len(A["prod"]) > len(B["prod"]):
                    A, B = B, A
                sm = difflib.SequenceMatcher(None, A["prod"], B["prod"], autojunk=False)
                blocks = sm.get_matching_blocks()
                for b1, b2 in zip(blocks, blocks[1:]):
                    ga = b2.a - (b1.a + b1.size)
                    gb = b2.b - (b1.b + b1.size)
                    if abs(gb - ga) >= a.min_indel:
                        junc[A["name"]].add(A["lf"] + b1.a + b1.size)
                        junc[B["name"]].add(B["lf"] + b1.b + b1.size)
                        junc[B["name"]].add(B["lf"] + b2.b)
                        if ga > 0:
                            junc[A["name"]].add(A["lf"] + b2.a)
        for A in alleles:
            fa_out.write(">%s\n%s\n" % (A["name"], A["seq"]))
            pts = sorted(junc[A["name"]])
            kept = [p for k, p in enumerate(pts) if k == 0 or p - pts[k - 1] > 5]
            for p in kept:
                jt.write("%s\t%s\t%d\t%s\t%d\n" % (A["name"], m, A["size"], ",".join(A["spp"]), p))
            if not kept:
                sys.stderr.write("WARN: no junction found for %s\n" % A["name"])
    fa_out.close(); jt.close(); info.close()


# ---------------------------------------------------------------- genotype
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def ref_blocks(pos, cigar):
    """aligned reference intervals (0-based) and largest deletion inside them"""
    p, blocks = pos, []
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=X":
            blocks.append((p, p + n)); p += n
        elif op in "DN":
            blocks.append((p, p + n, "D")); p += n
    return blocks


def spans(blocks, j, anchor):
    """read aligns continuously (no deletion > 5 bp) over [j-anchor, j+anchor)"""
    lo, hi = j - anchor, j + anchor
    cov = [b for b in blocks if b[1] > lo and b[0] < hi]
    if not cov or cov[0][0] > lo or cov[-1][1] < hi:
        return False
    for b in cov:
        if len(b) == 3 and b[1] - b[0] > 5:
            return False
    return all(cov[k][0] <= cov[k - 1][1] for k in range(1, len(cov)))


def cmd_genotype(a):
    J = collections.defaultdict(list)
    for r in csv.DictReader(open(a.junctions), delimiter="\t"):
        J[r["allele"]].append(int(r["junction"]))
    support = collections.defaultdict(set)
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        t = line.split("\t")
        flag, ref, pos, mapq, cigar = int(t[1]), t[2], int(t[3]) - 1, int(t[4]), t[5]
        if flag & 0x904 or ref not in J or mapq < a.min_mapq or cigar == "*":
            continue
        nm = next((int(x[5:]) for x in t[11:] if x.startswith("NM:i:")), 99)
        if nm > a.max_nm:
            continue
        bl = ref_blocks(pos, cigar)
        if any(spans(bl, j, a.anchor) for j in J[ref]):
            support[ref].add(t[0])
    with open(a.out, "w") as f:
        f.write("sample\tallele\treads\n")
        for al in sorted(J):
            f.write("%s\t%s\t%d\n" % (a.sample, al, len(support.get(al, ()))))


# --------------------------------------------------------------- summarize
def cmd_summarize(a):
    meta = {r["sample"]: r for r in csv.DictReader(open(a.samples), delimiter="\t")}
    alle = {r["allele"]: r for r in csv.DictReader(open(a.alleles), delimiter="\t")}
    counts = collections.defaultdict(dict)       # (sample, marker) -> allele -> reads
    for f in sorted(glob.glob(os.path.join(a.counts_dir, "*.counts.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            counts[(r["sample"], alle[r["allele"]]["marker"])][r["allele"]] = int(r["reads"])
    rows, summ = [], collections.defaultdict(collections.Counter)
    for (s, m), c in sorted(counts.items()):
        tot = sum(c.values())
        sp = meta.get(s, {}).get("species", "NA") or "NA"
        sex = meta.get(s, {}).get("sex", "NA") or "NA"
        called = [al for al, n in c.items() if n >= a.min_reads and n >= a.min_frac * tot]
        if tot < a.min_reads:
            call, status = "no_data", "no_data"
        else:
            call = "/".join(sorted(alle[x]["size"] + "bp" for x in called))
            exp = [al for al in c if sp in alle[al]["species"]]
            if len(called) == 1:
                status = ("one_allele" if sp not in SPECIES else       # hybrids, unknown species
                          "expected" if called[0] in exp else "unexpected")
            else:
                status = "two_alleles"
                if sex == "F":
                    status = "two_alleles_FEMALE"    # one Z: second band = W copy, duplication or contamination
        rows.append(dict(sample=s, species=sp, sex=sex, marker=m, total_reads=tot, call=call, status=status,
                         **{"reads_" + alle[al]["size"] + "bp": n for al, n in c.items()}))
        summ[(m, sp)][status] += 1
    cols = ["sample", "species", "sex", "marker", "total_reads", "call", "status"]
    extra = sorted({k for r in rows for k in r if k.startswith("reads_")}, key=lambda k: int(k[6:-2]))
    with open(a.out_prefix + ".genotypes.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols + extra, delimiter="\t", restval=0)
        w.writeheader(); w.writerows(rows)
    with open(a.out_prefix + ".summary.tsv", "w") as f:
        f.write("marker\tspecies\tn\texpected\tunexpected\ttwo_alleles\ttwo_alleles_FEMALE\tone_allele\tno_data\n")
        for (m, sp), c in sorted(summ.items()):
            f.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (m, sp, sum(c.values()), c["expected"], c["unexpected"],
                                                              c["two_alleles"], c["two_alleles_FEMALE"], c["one_allele"],
                                                              c["no_data"]))
    sys.stderr.write("wrote %s.genotypes.tsv and %s.summary.tsv\n" % (a.out_prefix, a.out_prefix))


# ---------------------------------------------------------------- refcheck
def read_chrmap(path):
    """accession -> chromosome name, from a 2+ column table; takes the column
    that looks like an accession (NC_/NW_/CM...) and the next short name."""
    m = {}
    if not path or not os.path.exists(path):
        return m
    if path.endswith("assembly_report.txt"):        # NCBI: col1 name, col7 RefSeq accession
        for line in open(path):
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) > 6 and t[6] not in ("na", ""):
                m[t[6]] = t[0]
        return m
    for line in open(path):
        t = [x.strip() for x in re.split(r"\t|,", line.rstrip("\n")) if x.strip()]
        acc = [x for x in t if re.match(r"^(NC|NW|NT|CM|JA)\w*\.\d+$", x)]
        name = [x for x in t if re.match(r"^(chr_?)?([0-9]+|[ZW]|MT|LG\w+)$", x, re.I)]
        if acc and name:
            m[acc[0]] = name[0]
    return m


def cmd_refcheck(a):
    sizes = {}
    for line in open(a.alleles_fa):
        if line.startswith(">"):
            cur = line[1:].strip(); sizes[cur] = 0
        else:
            sizes[cur] += len(line.strip())
    out = []
    for spec in a.refs:                      # label=paf=chrmap
        label, paf, cmap = (spec.split("=") + ["", ""])[:3]
        names = read_chrmap(cmap)
        hits = collections.defaultdict(list)
        if os.path.exists(paf):
            for line in open(paf):
                t = line.split("\t")
                ql, qs, qe, nm, bl = int(t[1]), int(t[2]), int(t[3]), int(t[9]), int(t[10])
                if bl == 0:
                    continue
                chrom = names.get(t[5], t[5])
                hits[t[0]].append((nm * (qe - qs) / bl, nm / bl, (qe - qs) / ql, chrom, t[5], int(t[7])))
        for q in sizes:
            H = sorted(hits.get(q, []), reverse=True)
            best = H[0] if H else None
            isW = lambda c: re.search(r"(^|_)W$", c, re.I)
            isZ = lambda c: re.search(r"(^|_)Z$", c, re.I)
            w = [h for h in H if isW(h[3]) and h[1] >= a.min_ident and h[2] >= a.min_cov]
            out.append(dict(
                allele=q, reference=label,
                best_chrom=best[3] if best else "no_hit", best_acc=best[4] if best else "",
                best_pos=best[5] if best else "", best_ident="%.3f" % best[1] if best else "",
                best_qcov="%.2f" % best[2] if best else "",
                n_hits=len(H), Z_hits=sum(1 for h in H if isZ(h[3])), W_hits=len(w),
                verdict=("W_HIT" if w else "Z_best" if best and isZ(best[3]) else
                         "no_hit" if not best else "best_not_Z")))
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader(); w.writerows(out)
    v = collections.Counter((r["reference"], r["verdict"]) for r in out)
    sys.stderr.write("refcheck: " + "; ".join("%s %s %d" % (k[0], k[1], n) for k, n in sorted(v.items())) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("alleles")
    p.add_argument("--markers", required=True); p.add_argument("--haps", required=True)
    p.add_argument("--pcr-dir", required=True); p.add_argument("--out-prefix", required=True)
    p.add_argument("--flank", type=int, default=400); p.add_argument("--size-tol", type=int, default=3)
    p.add_argument("--min-indel", type=int, default=8); p.add_argument("--max-len", type=int, default=3000)
    p = sub.add_parser("genotype")
    p.add_argument("--junctions", required=True); p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True); p.add_argument("--anchor", type=int, default=20)
    p.add_argument("--max-nm", type=int, default=4); p.add_argument("--min-mapq", type=int, default=1)
    p = sub.add_parser("summarize")
    p.add_argument("--samples", required=True); p.add_argument("--alleles", required=True)
    p.add_argument("--counts-dir", required=True); p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-reads", type=int, default=3); p.add_argument("--min-frac", type=float, default=0.15)
    p = sub.add_parser("refcheck")
    p.add_argument("--alleles-fa", required=True); p.add_argument("--out", required=True)
    p.add_argument("--refs", nargs="+", required=True, help="label=alignments.paf=chr_map.tsv")
    p.add_argument("--min-ident", type=float, default=0.8); p.add_argument("--min-cov", type=float, default=0.5)
    a = ap.parse_args()
    {"alleles": cmd_alleles, "genotype": cmd_genotype, "summarize": cmd_summarize, "refcheck": cmd_refcheck}[a.cmd](a)


if __name__ == "__main__":
    main()
