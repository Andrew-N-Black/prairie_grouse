#!/usr/bin/env python3
"""
22b_z_indel_markers.py -- find species-diagnostic Z-linked indels that can be
scored by band length on an agarose gel, design primers for them, and check
every primer pair by in-silico PCR.

Called by 22b_z_indel_markers.sh, one sub-command per step:

  candidates  parse paftools-call output for every Z-bearing haplotype, cluster
              indels across haplotypes, and keep windows where species differ
              consistently in net length (Z-linked length polymorphisms)
  primers     write a Primer3 input for each candidate (target = indel window,
              other variable sites excluded), run primer3_core, parse results
  report      combine in-silico PCR products (seqkit amplicon) on every Z
              haplotype, the W of each female and the LEPC reference genome
              into a ranked marker table

Standard library only (no pandas), so it runs with the system python3.
"""
import argparse, bisect, collections, csv, gzip, os, subprocess, sys

SPECIES = ["LEPC", "GRPC", "STGR"]
PAIRS = [("LEPC", "GRPC"), ("GRPC", "STGR"), ("LEPC", "STGR")]


# ----------------------------------------------------------------------------- io
def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_haps(path):
    """haplotypes.tsv: hap_id  species  z_fasta  calls(paftools call output)"""
    haps = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            haps.append(row)
    return haps


def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with opn(path) as f:
        for line in f:
            if line.startswith(">"):
                if name and (keep is None or name in keep):
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name and (keep is None or name in keep):
        seqs[name] = "".join(buf)
    return seqs


# ----------------------------------------------------------------------- calls
def parse_calls(path):
    """paftools.js call output.
    R ctg start end                    regions covered by exactly one alignment
    V ctg start end depth mapq ref alt qname qstart qend strand
    Coordinates are 0-based, half-open on the reference."""
    cov = collections.defaultdict(list)
    var = []
    with opn(path) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if t[0] == "R":
                cov[t[1]].append((int(t[2]), int(t[3])))
            elif t[0] == "V" and t[4] == "1":
                ref = "" if t[6] == "-" else t[6].upper()
                alt = "" if t[7] == "-" else t[7].upper()
                var.append((t[1], int(t[2]), int(t[3]), ref, alt))
    for c in cov:
        cov[c].sort()
    covs = {c: ([x[0] for x in v], v) for c, v in cov.items()}
    byctg = collections.defaultdict(list)
    for v in sorted(var, key=lambda v: (v[0], v[1])):
        byctg[v[0]].append(v)
    vidx = {c: ([v[1] for v in vs], vs) for c, vs in byctg.items()}
    return covs, vidx, var


def covered(covs, ctg, s, e):
    """True if [s,e) lies inside a single depth-1 covered interval."""
    if ctg not in covs:
        return False
    starts, iv = covs[ctg]
    i = bisect.bisect_right(starts, s) - 1
    return i >= 0 and iv[i][1] >= e


def vars_in(vidx, ctg, s, e):
    if ctg not in vidx:
        return []
    starts, vs = vidx[ctg]
    return vs[bisect.bisect_left(starts, s):bisect.bisect_left(starts, e)]


def is_tandem(seq, max_unit=6, min_copies=2):
    """Indel sequence made of a short repeat unit (STR/homopolymer): unstable
    and prone to length homoplasy, so excluded by default."""
    if not seq:
        return False
    for u in range(1, max_unit + 1):
        if len(seq) % u == 0 and len(seq) // u >= min_copies and seq == seq[:u] * (len(seq) // u):
            return True
    return False


def cmd_candidates(a):
    haps = read_haps(a.haps)
    ref = read_fasta(a.ref_z)
    calls, allvar = {}, []
    for h in haps:
        cov, vidx, var = parse_calls(h["calls"])
        calls[h["hap_id"]] = (h["species"], cov, vidx)
        allvar += [(v[0], v[1], v[2], len(v[4]) - len(v[3]), v[3] or v[4]) for v in var]
    # every variable site (any haplotype) -> used later to exclude primer positions
    with open(a.out_prefix + ".all_variable_sites.tsv", "w") as f:
        for v in sorted(set((x[0], x[1], x[2]) for x in allvar)):
            f.write("%s\t%d\t%d\n" % v)
    # indels >= min_indel in at least one haplotype
    ind = sorted(set((c, s, e) for c, s, e, d, sq in allvar if abs(d) >= a.min_indel))
    # cluster across haplotypes
    clusters = []
    for c, s, e in ind:
        if clusters and clusters[-1][0] == c and s - clusters[-1][2] <= a.cluster_gap:
            clusters[-1][2] = max(clusters[-1][2], e)
        else:
            clusters.append([c, s, e])
    nsp = collections.Counter(h["species"] for h in haps)
    out = []
    for c, s, e in clusters:
        ws, we = max(0, s - a.pad), e + a.pad
        lens, seqs = collections.defaultdict(list), []
        for hid, (sp, cov, vidx) in calls.items():
            if not covered(cov, c, ws, we):
                continue
            net = 0
            for vc, vs, ve, vr, vt in vars_in(vidx, c, ws, we):
                if len(vr) != len(vt):
                    net += len(vt) - len(vr)
                    if abs(len(vt) - len(vr)) >= a.min_indel:
                        seqs.append(vr or vt)
            lens[sp].append(net)
        # enough haplotypes called in every species compared
        ok = {sp: len(lens[sp]) >= max(a.min_hap, a.min_frac * nsp[sp]) for sp in SPECIES}
        res = {}
        for p, q in PAIRS:
            if not (ok[p] and ok[q]):
                continue
            lp, lq = lens[p], lens[q]
            if max(lp) - min(lp) > a.max_within or max(lq) - min(lq) > a.max_within:
                continue
            gap = min(abs(x - y) for x in lp for y in lq)
            if gap >= a.min_diff:
                res[p + "_" + q] = gap
        if not res:
            continue
        tandem = any(is_tandem(sq) for sq in seqs)
        if tandem and not a.keep_tandem:
            continue
        refwin = ref[c][ws:we] if c in ref else ""
        soft = sum(ch.islower() for ch in refwin) / max(1, len(refwin))
        out.append(dict(
            cand_id="%s_%d" % (c, s), contig=c, win_start=ws, win_end=we,
            pairs=",".join(res), min_gap=min(res.values()),
            three_way=int(len(res) == 3),
            net_LEPC=summ(lens["LEPC"]), net_GRPC=summ(lens["GRPC"]), net_STGR=summ(lens["STGR"]),
            n_LEPC=len(lens["LEPC"]), n_GRPC=len(lens["GRPC"]), n_STGR=len(lens["STGR"]),
            softmask_frac="%.2f" % soft, tandem=int(tandem)))
    write_tsv(a.out_prefix + ".candidates.tsv", out)
    sys.stderr.write("candidates: %d windows (%s)\n" % (
        len(out), ", ".join("%s %d" % (k, sum(k in r["pairs"] for r in out)) for k in ["LEPC_GRPC", "GRPC_STGR", "LEPC_STGR"])))


def summ(v):
    if not v:
        return "NA"
    return str(v[0]) if min(v) == max(v) else "%d..%d" % (min(v), max(v))


def write_tsv(path, rows, cols=None):
    with open(path, "w", newline="") as f:
        if not rows:
            f.write("" if cols is None else "\t".join(cols) + "\n")
            return
        w = csv.DictWriter(f, fieldnames=cols or list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --------------------------------------------------------------------- primers
def cmd_primers(a):
    ref = read_fasta(a.ref_z)
    cands = list(csv.DictReader(open(a.candidates), delimiter="\t"))
    cands = [r for r in cands if float(r["softmask_frac"]) <= a.max_soft]
    var = collections.defaultdict(list)
    for line in open(a.varsites):
        c, s, e = line.split()
        var[c].append((int(s), int(e)))
    vstart = {c: [x[0] for x in v] for c, v in var.items()}
    # per species pair: rank three-way markers first, then the largest length
    # difference; keep the best candidate within any --min-spacing bp, up to
    # --max-per-pair (so the scarcer LEPC-GRPC markers are not crowded out)
    chosen = {}
    for p, q in PAIRS:
        key = p + "_" + q
        pool = sorted((r for r in cands if key in r["pairs"].split(",")),
                      key=lambda r: (-int(r["three_way"]), -int(r["min_gap"])))
        got = []
        for r in pool:
            c, x = r["contig"], int(r["win_start"])
            if all(c != g["contig"] or abs(x - int(g["win_start"])) >= a.min_spacing for g in got):
                got.append(r)
            if len(got) >= a.max_per_pair:
                break
        for r in got:
            chosen[r["cand_id"]] = r
    chosen = list(chosen.values())
    boulder = []
    kept = []
    for r in chosen:
        c, ws, we = r["contig"], int(r["win_start"]), int(r["win_end"])
        ts, te = max(0, ws - a.flank), min(len(ref[c]), we + a.flank)
        tmpl = ref[c][ts:te].upper()
        excl = []
        lo, hi = bisect.bisect_left(vstart.get(c, []), ts), bisect.bisect_left(vstart.get(c, []), te)
        for s, e in var[c][lo:hi]:
            if not (ws <= s < we):
                excl.append("%d,%d" % (s - ts, max(1, e - s)))
        boulder.append("\n".join([
            "SEQUENCE_ID=" + r["cand_id"],
            "SEQUENCE_TEMPLATE=" + tmpl,
            "SEQUENCE_TARGET=%d,%d" % (ws - ts, we - ws),
            ] + (["SEQUENCE_EXCLUDED_REGION=" + " ".join(excl)] if excl else []) + [
            "PRIMER_TASK=generic", "PRIMER_PICK_LEFT_PRIMER=1", "PRIMER_PICK_RIGHT_PRIMER=1",
            "PRIMER_OPT_SIZE=21", "PRIMER_MIN_SIZE=18", "PRIMER_MAX_SIZE=25",
            "PRIMER_OPT_TM=60.0", "PRIMER_MIN_TM=57.0", "PRIMER_MAX_TM=63.0",
            "PRIMER_MIN_GC=35", "PRIMER_MAX_GC=65", "PRIMER_MAX_POLY_X=4",
            "PRIMER_PRODUCT_SIZE_RANGE=%s" % a.product_range,
            "PRIMER_NUM_RETURN=1", "PRIMER_EXPLAIN_FLAG=1", "="]))
        kept.append(r)
        r["_tmpl"] = tmpl
    p3in = a.out_prefix + ".primer3_input.txt"
    open(p3in, "w").write("\n".join(boulder) + "\n")
    res = subprocess.run(["primer3_core", p3in], capture_output=True, text=True, check=True).stdout
    rec, rows = {}, []
    info = {r["cand_id"]: r for r in kept}
    tmpls = {}
    for line in res.splitlines():
        if line == "=":
            if rec.get("PRIMER_PAIR_NUM_RETURNED", "0") != "0":
                r = dict(info[rec["SEQUENCE_ID"]])
                r.update(left=rec["PRIMER_LEFT_0_SEQUENCE"], right=rec["PRIMER_RIGHT_0_SEQUENCE"],
                         ref_product=rec["PRIMER_PAIR_0_PRODUCT_SIZE"],
                         tm_left="%.1f" % float(rec["PRIMER_LEFT_0_TM"]),
                         tm_right="%.1f" % float(rec["PRIMER_RIGHT_0_TM"]))
                ls = int(rec["PRIMER_LEFT_0"].split(",")[0])
                re_ = int(rec["PRIMER_RIGHT_0"].split(",")[0])
                tmpls[r["cand_id"]] = r["_tmpl"][ls:re_ + 1]
                rows.append(r)
            rec = {}
        elif "=" in line:
            k, v = line.split("=", 1)
            rec[k] = v
    write_tsv(a.out_prefix + ".primers.tsv", rows)
    with open(a.out_prefix + ".amplicons.fa", "w") as f:        # LEPC-reference products
        for k, v in tmpls.items():
            f.write(">%s\n%s\n" % (k, v))
    with open(a.out_prefix + ".primers_seqkit.tsv", "w") as f:   # name  F  R
        for r in rows:
            f.write("%s\t%s\t%s\n" % (r["cand_id"], r["left"], r["right"]))
    sys.stderr.write("primers: designed for %d of %d candidates\n" % (len(rows), len(kept)))


# ---------------------------------------------------------------------- report
def read_amplicons(path, max_len=3000):
    """seqkit amplicon --bed output: seq start end primer_name score strand"""
    hits = collections.defaultdict(list)
    if os.path.exists(path):
        for line in open(path):
            t = line.split("\t")
            n = int(t[2]) - int(t[1])
            if n <= max_len:
                hits[t[3]].append(n)
    return hits


def cmd_report(a):
    rows = list(csv.DictReader(open(a.primers), delimiter="\t"))
    haps = read_haps(a.haps)
    zprod = {h["hap_id"]: read_amplicons(os.path.join(a.pcr_dir, "Z." + h["hap_id"] + ".bed")) for h in haps}
    wfiles = [f for f in os.listdir(a.pcr_dir) if f.startswith("W.")]
    wprod = [read_amplicons(os.path.join(a.pcr_dir, f)) for f in wfiles]
    gprod = read_amplicons(os.path.join(a.pcr_dir, "refgenome.bed"))
    # W gametolog homology: product aligned to any female W (minimap2 PAF)
    whom = collections.Counter()
    wpaf = os.path.join(a.pcr_dir, "W_homology.paf")
    if os.path.exists(wpaf):
        for line in open(wpaf):
            t = line.split("\t")
            qlen, qs, qe, nm, bl = int(t[1]), int(t[2]), int(t[3]), int(t[9]), int(t[10])
            if bl and nm / bl >= a.w_min_ident and (qe - qs) / qlen >= a.w_min_cov:
                whom[t[0]] += 1
    out = []
    for r in rows:
        m = r["cand_id"]
        sizes = collections.defaultdict(list)
        multi = 0
        for h in haps:
            p = zprod[h["hap_id"]].get(m, [])
            if len(p) == 1:
                sizes[h["species"]].append(p[0])
            elif len(p) > 1:
                multi += 1
        w_hits = sum(len(w.get(m, [])) for w in wprod) + whom[m]
        g_hits = len(gprod.get(m, []))
        # species band sizes and separation from in-silico PCR on the haplotypes
        band = {sp: summ(sorted(sizes[sp])) for sp in SPECIES}
        sep = {}
        for p, q in PAIRS:
            if sizes[p] and sizes[q]:
                gap = min(abs(x - y) for x in sizes[p] for y in sizes[q])
                hi = max(sizes[p] + sizes[q])
                sep[p + "_" + q] = (gap, gap / hi)
        resolvable = [k for k, (g, f) in sep.items() if g >= a.gel_min_bp and f >= a.gel_min_frac]
        amplified = sum(len(v) for v in sizes.values())
        pass_all = int(amplified >= a.min_amplified_frac * len(haps) and multi == 0
                       and w_hits == 0 and g_hits == 1 and bool(resolvable))
        out.append(dict(
            marker=m, contig=r["contig"], win_start=r["win_start"], win_end=r["win_end"],
            pairs_gel=",".join(resolvable) or "none", pass_all=pass_all,
            band_LEPC=band["LEPC"], band_GRPC=band["GRPC"], band_STGR=band["STGR"],
            gap_LEPC_GRPC=sep.get("LEPC_GRPC", ("NA",))[0], gap_GRPC_STGR=sep.get("GRPC_STGR", ("NA",))[0],
            gap_LEPC_STGR=sep.get("LEPC_STGR", ("NA",))[0],
            haps_amplified="%d/%d" % (amplified, len(haps)), haps_multi_product=multi,
            W_products=w_hits, refgenome_products=g_hits,
            left=r["left"], right=r["right"], tm_left=r["tm_left"], tm_right=r["tm_right"],
            ref_product=r["ref_product"], softmask_frac=r["softmask_frac"]))
    out.sort(key=lambda d: (-d["pass_all"], -len(d["pairs_gel"].split(",")) if d["pairs_gel"] != "none" else 0,
                            -int(d["gap_LEPC_GRPC"]) if d["gap_LEPC_GRPC"] != "NA" else 0))
    write_tsv(a.out_prefix + ".markers.tsv", out)
    n = collections.Counter()
    for d in out:
        if d["pass_all"]:
            for k in d["pairs_gel"].split(","):
                n[k] += 1
    sys.stderr.write("report: %d markers pass all checks (%s)\n" % (
        sum(d["pass_all"] for d in out), ", ".join("%s %d" % kv for kv in sorted(n.items()))))


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("candidates")
    c.add_argument("--haps", required=True); c.add_argument("--ref-z", required=True)
    c.add_argument("--out-prefix", required=True)
    c.add_argument("--min-indel", type=int, default=10, help="smallest indel considered (bp)")
    c.add_argument("--min-diff", type=int, default=15, help="min length difference between species (bp)")
    c.add_argument("--max-within", type=int, default=2, help="max length range within a species (bp)")
    c.add_argument("--cluster-gap", type=int, default=20)
    c.add_argument("--pad", type=int, default=25, help="bp added either side of an indel cluster")
    c.add_argument("--min-frac", type=float, default=0.8, help="min fraction of a species' haplotypes called")
    c.add_argument("--min-hap", type=int, default=4)
    c.add_argument("--keep-tandem", action="store_true", help="keep STR/homopolymer indels")
    p = sub.add_parser("primers")
    p.add_argument("--candidates", required=True); p.add_argument("--ref-z", required=True)
    p.add_argument("--varsites", required=True); p.add_argument("--out-prefix", required=True)
    p.add_argument("--flank", type=int, default=250)
    p.add_argument("--product-range", default="150-450")
    p.add_argument("--max-soft", type=float, default=0.2)
    p.add_argument("--min-spacing", type=int, default=10000, help="min distance between candidates (bp)")
    p.add_argument("--max-per-pair", type=int, default=150,
                   help="top-ranked candidates per species pair sent to Primer3 and in-silico PCR")
    r = sub.add_parser("report")
    r.add_argument("--primers", required=True); r.add_argument("--haps", required=True)
    r.add_argument("--pcr-dir", required=True); r.add_argument("--out-prefix", required=True)
    r.add_argument("--gel-min-bp", type=int, default=15)
    r.add_argument("--gel-min-frac", type=float, default=0.05, help="min difference as fraction of larger band")
    r.add_argument("--min-amplified-frac", type=float, default=0.9)
    r.add_argument("--w-min-ident", type=float, default=0.8, help="W homology: min identity of product alignment")
    r.add_argument("--w-min-cov", type=float, default=0.5, help="W homology: min fraction of product aligned")
    a = ap.parse_args()
    {"candidates": cmd_candidates, "primers": cmd_primers, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
