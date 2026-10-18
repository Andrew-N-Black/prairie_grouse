#!/usr/bin/env python3
"""
Helper for 14_z_markers.sh (species-diagnostic Z-chromosome markers).

  discover   Merge per-haplotype variant calls against the LEPC reference Z (paftools call) into
             one table of sites with each haplotype's allele, and keep the sites that separate
             species: every called haplotype of species A carries one allele, every called
             haplotype of species B another. Writes candidates.tsv, amplicons.fa (window around
             each candidate, for the uniqueness/W checks), primer3 input, and sites.txt (SNPs, for
             the ANGSD panel check).
  hits       Parse minimap2 PAF of the amplicons against the whole LEPC genome and against the
             female W sequences -> number of genomic copies and W matches per candidate.
  primers    Parse primer3_core output.
  panel      Merge ANGSD .mafs.gz of every panel group into per-candidate allele frequencies and
             write the final ranked table.

Haplotype statuses at a site: A = carries the variant allele, R = reference allele (the position
lies in a region covered by exactly one alignment of that haplotype and no variant is called there),
O = a different variant overlaps (multi-allelic or messy; treated as uncalled), . = not covered.
"""
import argparse
import bisect
import csv
import glob
import gzip
import os
import sys
from collections import defaultdict


def species_of(hap):
    return hap.split("_")[0]


# ------------------------------------------------------------------------------------------ discover
def read_call(path):
    """paftools call output -> (regions{ctg:[(s,e)]}, variants{key:(s,e)}) ; depth-1 variants only."""
    regions, var = defaultdict(list), {}
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if t[0] == "R":
                regions[t[1]].append((int(t[2]), int(t[3])))
            elif t[0] == "V" and int(t[4]) == 1:          # covered by exactly one contig
                s, e = int(t[2]), int(t[3])
                var[(t[1], s, e, t[6].upper(), t[7].upper())] = (s, e)
    for c in regions:
        regions[c].sort()
    return regions, var


def covered(regs, s, e):
    """True if [s, e) (or the insertion point s for s == e) lies inside one depth-1 region."""
    if not regs:
        return False
    lo, hi = (s - 1, s + 1) if s == e else (s, e)
    i = bisect.bisect_right(regs, (lo, float("inf"))) - 1
    return i >= 0 and regs[i][0] <= lo and regs[i][1] >= hi


def cmd_discover(a):
    haps, calls = [], {}
    for f in sorted(glob.glob(os.path.join(a.calls, "*.call.txt"))):
        h = os.path.basename(f)[:-len(".call.txt")]
        haps.append(h)
        calls[h] = read_call(f)
    if not haps:
        sys.exit("no *.call.txt in " + a.calls)
    species = sorted({species_of(h) for h in haps})
    by_sp = {sp: [h for h in haps if species_of(h) == sp] for sp in species}
    sys.stderr.write("[discover] %d haplotypes: %s\n" % (len(haps), ", ".join(
        "%s %d" % (sp, len(v)) for sp, v in by_sp.items())))

    sites = set()
    for h in haps:
        sites.update(calls[h][1].keys())
    # per-haplotype sorted variant starts, to spot a different variant at the same position
    starts = {h: defaultdict(list) for h in haps}
    for h in haps:
        for (c, s, e, r, al) in calls[h][1]:
            starts[h][c].append((s, e))
        for c in starts[h]:
            starts[h][c].sort()
    # every variable position (any haplotype), for flank cleanliness and primer3 exclusions
    varpos = defaultdict(list)
    for (c, s, e, r, al) in sites:
        varpos[c].append((s, max(e, s + 1)))
    for c in varpos:
        varpos[c].sort()

    def status(h, site):
        c, s, e, r, al = site
        if site in calls[h][1]:
            return "A"
        lst = starts[h].get(c, [])
        i = bisect.bisect_left(lst, (s - 50, 0))
        while i < len(lst) and lst[i][0] < max(e, s + 1):
            if lst[i][1] > s or lst[i][0] == s:        # another variant overlaps this site
                return "O"
            i += 1
        return "R" if covered(calls[h][0].get(c), s, e) else "."

    pairs = [(x, y) for i, x in enumerate(species) for y in species[i + 1:]]
    out = []
    for site in sorted(sites):
        st = {h: status(h, site) for h in haps}
        alle = {}
        ok = True
        for sp in species:
            called = [st[h] for h in by_sp[sp] if st[h] in "AR"]
            need = max(a.min_copies, int(a.min_called * len(by_sp[sp]) + 0.9999))
            if len(called) < need:
                ok = False
                break
            alle[sp] = set(called)
        if not ok:
            continue
        diag = [x + "-" + y for x, y in pairs
                if len(alle[x]) == 1 and len(alle[y]) == 1 and alle[x] != alle[y]]
        if not diag:
            continue
        c, s, e, r, al = site
        kind = "SNP" if len(r) == 1 and len(al) == 1 and "-" not in (r, al) else "INDEL"
        size = max(len(r.replace("-", "")), len(al.replace("-", "")))
        # distance to the nearest other variable site
        vp = varpos[c]
        j = bisect.bisect_left(vp, (s, 0))
        left = s - max([q[1] for q in vp[max(0, j - 5):j] if q[1] <= s] or [-10 ** 9])
        right = min([q[0] for q in vp[j:j + 6] if q[0] >= max(e, s + 1)] or [10 ** 9]) - max(e, s + 1)
        out.append({
            "id": "%s:%d:%s>%s" % (c, s + 1, r, al), "chrom": c, "pos": s + 1, "end": max(e, s + 1),
            "ref": r, "alt": al, "type": kind, "size": size, "diagnostic_pairs": ",".join(diag),
            **{"allele_" + sp: "".join(sorted(alle[sp])) for sp in species},
            **{"n_" + sp: sum(st[h] in "AR" for h in by_sp[sp]) for sp in species},
            "clean_left": min(left, 99999), "clean_right": min(right, 99999),   # bp to the nearest other variant
        })
    sys.stderr.write("[discover] %d variable sites, %d separate at least one species pair\n"
                     % (len(sites), len(out)))
    os.makedirs(a.out, exist_ok=True)
    fields = list(out[0].keys()) if out else ["id"]
    with open(os.path.join(a.out, "candidates.raw.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    # amplicon windows for the copy-number / W checks (bed, 0-based), and variable positions
    with open(os.path.join(a.out, "windows.bed"), "w") as fh:
        for o in out:
            s0 = max(0, o["pos"] - 1 - a.flank)
            fh.write("%s\t%d\t%d\t%s\n" % (o["chrom"], s0, o["end"] + a.flank, o["id"]))
    with open(os.path.join(a.out, "varpos.bed"), "w") as fh:
        for c in sorted(varpos):
            for s, e in varpos[c]:
                fh.write("%s\t%d\t%d\n" % (c, s, e))
    for p in [x + "-" + y for x, y in pairs]:
        n = sum(p in o["diagnostic_pairs"].split(",") for o in out)
        sys.stderr.write("    %-10s %7d sites (%d SNP, %d indel)\n" % (p, n,
                         sum(p in o["diagnostic_pairs"].split(",") and o["type"] == "SNP" for o in out),
                         sum(p in o["diagnostic_pairs"].split(",") and o["type"] == "INDEL" for o in out)))


# ------------------------------------------------------------------------------------------ hits
def read_fa(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def good_hits(paf, min_id, min_cov):
    n = defaultdict(int)
    if not os.path.exists(paf):
        return n
    with open(paf) as fh:
        for line in fh:
            t = line.split("\t")
            qlen, qs, qe, match, blen = int(t[1]), int(t[2]), int(t[3]), int(t[9]), int(t[10])
            if blen and match / blen >= min_id and (qe - qs) / qlen >= min_cov:
                n[t[0]] += 1
    return n


def cmd_hits(a):
    seqs = read_fa(a.amplicons)
    gen = good_hits(a.genome_paf, a.min_id, a.min_cov)
    w = good_hits(a.w_paf, a.min_id, a.min_cov)
    vp = defaultdict(list)
    with open(a.varpos) as fh:
        for line in fh:
            c, s, e = line.split()
            vp[c].append((int(s), int(e)))
    rows = list(csv.DictReader(open(a.cand), delimiter="\t"))
    win = {}
    with open(a.windows) as fh:
        for line in fh:
            c, s, e, i = line.split()
            win[i] = (c, int(s), int(e))
    p3 = open(a.primer3_in, "w")
    for r in rows:
        i = r["id"]
        sq = seqs.get(i, "")
        r["masked_frac"] = "%.2f" % (sum(ch.islower() for ch in sq) / len(sq)) if sq else "NA"
        r["N_frac"] = "%.2f" % (sq.upper().count("N") / len(sq)) if sq else "NA"
        r["genome_copies"] = gen.get(i, 0)
        r["W_hits"] = w.get(i, 0)
        c, s0, e0 = win[i]
        others = [(s - s0, e - s) for s, e in vp[c][bisect.bisect_left(vp[c], (s0, 0)):
                                                     bisect.bisect_left(vp[c], (e0, 0))]
                  if not (s == int(r["pos"]) - 1)]
        r["other_var_in_window"] = len(others)
        if sq and r["N_frac"] != "NA" and float(r["N_frac"]) == 0:
            tgt = int(r["pos"]) - 1 - s0
            tlen = max(1, int(r["end"]) - int(r["pos"]) + 1)
            p3.write("SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\nSEQUENCE_TARGET=%d,%d\n" % (i, sq.upper(), max(0, tgt - 2), tlen + 4))
            if others:
                p3.write("SEQUENCE_EXCLUDED_REGION=%s\n" % " ".join("%d,%d" % (o, max(1, l)) for o, l in others if o >= 0))
            p3.write("=\n")
    p3.close()
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)


# ------------------------------------------------------------------------------------------ primers
def cmd_primers(a):
    res, cur = {}, {}
    with open(a.p3out) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "=":
                if "SEQUENCE_ID" in cur:
                    res[cur["SEQUENCE_ID"]] = cur
                cur = {}
            elif "=" in line:
                k, v = line.split("=", 1)
                cur[k] = v
    rows = list(csv.DictReader(open(a.cand), delimiter="\t"))
    for r in rows:
        p = res.get(r["id"], {})
        r["primer_F"] = p.get("PRIMER_LEFT_0_SEQUENCE", "")
        r["primer_R"] = p.get("PRIMER_RIGHT_0_SEQUENCE", "")
        r["product_bp"] = p.get("PRIMER_PAIR_0_PRODUCT_SIZE", "")
        r["Tm_F"] = p.get("PRIMER_LEFT_0_TM", "")[:5]
        r["Tm_R"] = p.get("PRIMER_RIGHT_0_TM", "")[:5]
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    sys.stderr.write("[primers] primer pairs for %d of %d candidates\n" % (sum(bool(r["primer_F"]) for r in rows), len(rows)))


# ------------------------------------------------------------------------------------------ panel
def cmd_panel(a):
    rows = list(csv.DictReader(open(a.cand), delimiter="\t"))
    freq, nind, groups = defaultdict(dict), defaultdict(dict), []
    for f in sorted(glob.glob(os.path.join(a.mafs_dir, "*.mafs.gz"))):
        g = os.path.basename(f)[:-len(".mafs.gz")]
        if g.endswith(".tmp"):                     # unfinished ANGSD run
            continue
        groups.append(g)
        with gzip.open(f, "rt") as fh:
            hdr = fh.readline().split()
            ic, ip = hdr.index("chromo"), hdr.index("position")
            iq = hdr.index("knownEM") if "knownEM" in hdr else hdr.index("unknownEM")
            im = hdr.index("minor")
            iN = hdr.index("nInd") if "nInd" in hdr else None
            for line in fh:
                t = line.split()
                key = (t[ic], int(t[ip]))
                freq[key][g] = (t[im].upper(), float(t[iq]))
                nind[key][g] = int(t[iN]) if iN is not None else -1
    for r in rows:
        key = (r["chrom"], int(r["pos"]))
        for g in groups:
            m = freq.get(key, {}).get(g)
            # sites.txt gives minor = alt, so knownEM is the alt-allele frequency
            r["altfreq_" + g] = "%.3f" % m[1] if m and m[0] == r["alt"] else ("NA" if not m else "minor!=alt")
            r["nInd_" + g] = nind.get(key, {}).get(g, "")
        best = []
        for p in r["diagnostic_pairs"].split(","):
            x, y = p.split("-")
            fx, fy = r.get("altfreq_%s_allo" % x, "NA"), r.get("altfreq_%s_allo" % y, "NA")
            try:
                d = abs(float(fx) - float(fy))
                r["dAF_" + p] = "%.3f" % d
                if d >= a.min_daf:
                    best.append(p)
            except ValueError:
                r["dAF_" + p] = "NA"
        r["panel_confirmed_pairs"] = ",".join(best)
    cols = list(rows[0].keys())
    for p in ("LEPC-GRPC", "GRPC-LEPC", "GRPC-STGR", "LEPC-STGR"):
        if "dAF_" + p not in cols and any("dAF_" + p in r for r in rows):
            cols.append("dAF_" + p)
    rows.sort(key=lambda r: (-(len(r["panel_confirmed_pairs"].split(",")) if r["panel_confirmed_pairs"] else 0),
                             r["type"] != "INDEL"))
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", restval="NA")
        wr.writeheader()
        wr.writerows(rows)
    for p in sorted({p for r in rows for p in r["diagnostic_pairs"].split(",")}):
        n = sum(p in r["panel_confirmed_pairs"].split(",") for r in rows)
        sys.stderr.write("[panel] %-10s %d markers with allele-frequency difference >= %.2f in allopatric birds\n"
                         % (p, n, a.min_daf))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd")
    s = sp.add_parser("discover")
    s.add_argument("--calls", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--min-called", type=float, default=0.8)
    s.add_argument("--min-copies", type=int, default=4)
    s.add_argument("--flank", type=int, default=150)
    s = sp.add_parser("hits")
    for x in ("cand", "amplicons", "windows", "varpos", "genome-paf", "w-paf", "primer3-in", "out"):
        s.add_argument("--" + x, required=True)
    s.add_argument("--min-id", type=float, default=0.85)
    s.add_argument("--min-cov", type=float, default=0.8)
    s = sp.add_parser("primers")
    s.add_argument("--cand", required=True)
    s.add_argument("--p3out", required=True)
    s.add_argument("--out", required=True)
    s = sp.add_parser("panel")
    s.add_argument("--cand", required=True)
    s.add_argument("--mafs-dir", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--min-daf", type=float, default=0.95)
    a = ap.parse_args()
    if not a.cmd:
        ap.error("choose a subcommand")
    {"discover": cmd_discover, "hits": cmd_hits, "primers": cmd_primers, "panel": cmd_panel}[a.cmd](a)


if __name__ == "__main__":
    main()
