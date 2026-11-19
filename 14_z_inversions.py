#!/usr/bin/env python3
"""Step 06 analysis: large inversions on the Z chromosome from the pangenome graph.

Inputs (written by 06_z_inversions.sh):
  untangle.tsv     odgi untangle of every Z path against one reference Z path
  z_paths.fa       sequences of every Z path (odgi paths -f), for gaps and path lengths
  ref.txt          the reference path name
  ref_to_ptarmigan.paf (optional)  reference Z aligned to the rock ptarmigan Z, for coordinates
Python 3 standard library only.

Method
  1. For each query path, every untangled segment is forward or reverse relative to the
     reference (odgi's `inv` = fraction of the segment traversed inverted; > 0.5 = reverse).
     The path's majority orientation (by bp) is "normal"; segments in the other orientation
     are inverted. Inverted segments closer than --merge bp on the reference are merged, and
     blocks spanning >= --min-inv bp of the reference are kept as inversion calls.
  2. A call's breakpoints are flagged when a gap (run of >= 10 N) or the end of the path lies
     within --flank bp: scaffolding gaps, RagTag joins and assembly-piece ends are where a
     mis-oriented scaffold would create a false inversion.
  3. Calls from different haplotypes with >= 50% reciprocal overlap on the reference are one
     event. A haplotype is "called" for an event if its untangled segments cover >= --cover of
     the event's span; called haplotypes that are not carriers are non-carriers. The reference
     haplotype is a non-carrier by definition.
  4. Events are classified like the structural variants of step 05: diagnostic for species X
     (X fixed for one state, both others fixed for the other), fixed between a species pair, or
     polymorphic / private. A species must have >= --min-called-frac of its Z haplotypes called
     (haplotypes whose Z path length is >= half the median; female haplotypes holding only a
     fragment of the Z are not counted).
"""
import argparse, bisect, math, os, re
from collections import defaultdict

SPECIES = ["GRPC", "LEPC", "STGR"]
PAIRS = [("GRPC", "LEPC"), ("GRPC", "STGR"), ("LEPC", "STGR")]


def hap_of(path):                      # PanSN species.sample#hap#contig -> species.sample#hap
    return "#".join(path.split("#")[:2])


def sp_of(path):
    return path.split(".")[0]


def read_paths(fa):
    """Path lengths and N-run intervals (>= 10 N) per path."""
    lens, gaps, name, pos, run = {}, defaultdict(list), None, 0, None
    def close():
        nonlocal run
        if run is not None and pos - run >= 10:
            gaps[name].append((run, pos))
        run = None
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    close(); lens[name] = pos
                name, pos, run = line[1:].split()[0], 0, None
                continue
            s = line.strip()
            for m in re.finditer(r"[^Nn]+|[Nn]+", s):
                seg = m.group(0)
                if seg[0] in "Nn":
                    if run is None:
                        run = pos
                else:
                    close()
                pos += len(seg)
    if name is not None:
        close(); lens[name] = pos
    return lens, gaps


def read_untangle(path):
    rows, cols = [], None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").lstrip("#").split("\t")
            if cols is None and "query.name" in f:
                cols = {c: i for i, c in enumerate(f)}
                continue
            if cols is None or len(f) < len(cols):
                continue
            rows.append(dict(q=f[cols["query.name"]], qs=int(f[cols["query.start"]]), qe=int(f[cols["query.end"]]),
                             r=f[cols["ref.name"]], rs=int(f[cols["ref.start"]]), re=int(f[cols["ref.end"]]),
                             inv=float(f[cols["inv"]])))
    if cols is None:
        raise SystemExit(f"ERROR: no header with query.name in {path}; check odgi untangle output format")
    return rows


def union_len(iv):
    tot, end = 0, -1
    for a, b in sorted(iv):
        if b > end:
            tot += b - max(a, end); end = b
    return tot


def covered(iv, a, b):
    return union_len([(max(x, a), min(y, b)) for x, y in iv if y > a and x < b])


def load_liftover(paf):
    """Blocks (qs, qe, ts, te, strand, alnlen) from the reference Z -> ptarmigan Z PAF."""
    if not paf or not os.path.exists(paf):
        return None
    bl = []
    for line in open(paf):
        f = line.split("\t")
        if len(f) >= 12 and int(f[11]) >= 1:
            bl.append((int(f[2]), int(f[3]), int(f[7]), int(f[8]), f[4], int(f[10])))
    return bl or None


def lift(bl, p):
    best = None
    for qs, qe, ts, te, st, al in bl:
        if qs <= p <= qe and (best is None or al > best[5]):
            best = (qs, qe, ts, te, st, al)
    if best is None:
        return None
    qs, qe, ts, te, st, _ = best
    return ts + (p - qs) if st == "+" else te - (p - qs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--min-inv", type=int, default=50000)
    ap.add_argument("--merge", type=int, default=10000)
    ap.add_argument("--flank", type=int, default=10000)
    ap.add_argument("--cover", type=float, default=0.8)
    ap.add_argument("--min-called-frac", type=float, default=0.7)
    ap.add_argument("--segment", default="41100000-44300000", help="ptarmigan Z interval to test for overlap")
    a = ap.parse_args()
    D = a.dir
    ref = open(os.path.join(D, "ref.txt")).read().split()[0]
    lens, gaps = read_paths(os.path.join(D, "z_paths.fa"))
    rows = read_untangle(os.path.join(D, "untangle.tsv"))
    bl = load_liftover(os.path.join(D, "ref_to_ptarmigan.paf"))
    seg_lo, seg_hi = (int(x) for x in a.segment.split("-"))

    # haplotypes carrying the Z (total path length >= half the median), per species
    hap_len = defaultdict(int)
    for p, L in lens.items():
        hap_len[hap_of(p)] += L
    med = sorted(hap_len.values())[len(hap_len) // 2]
    zhaps = {h for h, L in hap_len.items() if L >= 0.5 * med}
    n_sp = {s: sum(1 for h in zhaps if sp_of(h) == s) for s in SPECIES}
    need = {s: max(1, math.ceil(a.min_called_frac * n_sp[s])) for s in SPECIES}

    # 1. inversion calls per query path
    by_path = defaultdict(list)
    for r in rows:
        if r["r"] == ref:
            by_path[r["q"]].append(r)
    cover_iv = defaultdict(list)                       # haplotype -> reference intervals covered
    calls = []
    for q, segs in by_path.items():
        for s in segs:
            cover_iv[hap_of(q)].append((min(s["rs"], s["re"]), max(s["rs"], s["re"])))
        fwd = sum(s["qe"] - s["qs"] for s in segs if s["inv"] <= 0.5)
        rev = sum(s["qe"] - s["qs"] for s in segs if s["inv"] > 0.5)
        normal_rev = rev > fwd
        opp = sorted([s for s in segs if (s["inv"] > 0.5) != normal_rev], key=lambda s: min(s["rs"], s["re"]))
        blocks = []
        for s in opp:
            rs, re_ = min(s["rs"], s["re"]), max(s["rs"], s["re"])
            if blocks and rs - blocks[-1]["re"] <= a.merge:
                b = blocks[-1]
                b["re"] = max(b["re"], re_); b["qs"] = min(b["qs"], s["qs"]); b["qe"] = max(b["qe"], s["qe"])
            else:
                blocks.append(dict(rs=rs, re=re_, qs=s["qs"], qe=s["qe"]))
        for b in blocks:
            if b["re"] - b["rs"] < a.min_inv:
                continue
            g = gaps.get(q, []); L = lens.get(q, 0)
            def near(p):
                return p < a.flank or L - p < a.flank or any(x - a.flank <= p <= y + a.flank for x, y in g)
            calls.append(dict(hap=hap_of(q), sp=sp_of(q), path=q, rs=b["rs"], re=b["re"], qs=b["qs"], qe=b["qe"],
                              near_gap=near(b["qs"]) or near(b["qe"])))

    # 2. cluster calls into events (>= 50% reciprocal overlap on the reference)
    events = []
    for c in sorted(calls, key=lambda c: c["rs"]):
        for e in events:
            ov = min(e["re"], c["re"]) - max(e["rs"], c["rs"])
            if ov > 0 and ov >= 0.5 * (c["re"] - c["rs"]) and ov >= 0.5 * (e["re"] - e["rs"]):
                e["calls"].append(c); e["rs"] = min(e["rs"], c["rs"]); e["re"] = max(e["re"], c["re"]); break
        else:
            events.append(dict(rs=c["rs"], re=c["re"], calls=[c]))

    # 3. carriers, called haplotypes and classification
    cover_iv[hap_of(ref)].append((0, lens.get(ref, 10**12)))
    out_ev = []
    for i, e in enumerate(sorted(events, key=lambda e: e["rs"]), start=1):
        span = e["re"] - e["rs"]
        carriers = {c["hap"] for c in e["calls"]} & zhaps
        called = {h for h in zhaps if h in carriers or covered(cover_iv[h], e["rs"], e["re"]) >= a.cover * span}
        n = {s: sum(1 for h in called if sp_of(h) == s) for s in SPECIES}
        k = {s: sum(1 for h in carriers if sp_of(h) == s) for s in SPECIES}
        ok = all(n[s] >= need[s] for s in SPECIES)
        fr = {s: (k[s] / n[s] if n[s] else float("nan")) for s in SPECIES}
        st = {s: (1 if fr[s] == 1 else 0 if fr[s] == 0 else None) for s in SPECIES}
        cls, fixed = "not_classified", []
        if ok:
            diag = [x for x in SPECIES if st[x] is not None and all(st[o] == 1 - st[x] for o in SPECIES if o != x)]
            fixed = [f"{p}-{q}" for p, q in PAIRS if st[p] is not None and st[q] is not None and st[p] != st[q]]
            car_sp = [s for s in SPECIES if k[s] > 0]
            if diag:
                cls = f"diagnostic_{diag[0]}"
            elif fixed:
                cls = "fixed_" + ",".join(fixed)
            elif len(car_sp) == 1:
                cls = f"private_{car_sp[0]}"
            else:
                cls = "polymorphic"
        clean = sum(1 for c in e["calls"] if not c["near_gap"] and c["hap"] in carriers)
        evid = "multiple_clean" if clean >= 2 else "single_clean" if clean == 1 else "breakpoints_at_gaps_only"
        pt = (lift(bl, e["rs"]), lift(bl, e["re"])) if bl else (None, None)
        in_seg = ""
        if pt[0] is not None and pt[1] is not None:
            lo, hi = sorted(pt)
            in_seg = "yes" if lo < seg_hi and hi > seg_lo else "no"
        out_ev.append(dict(id=f"ZINV{i:03d}", rs=e["rs"], re=e["re"], span=span, pt=pt, n=n, k=k, fr=fr, cls=cls,
                           fixed=fixed, clean=clean, evid=evid, in_seg=in_seg,
                           carriers=sorted(carriers), calls=e["calls"]))

    # ---- outputs ----
    with open(os.path.join(D, "inversion_calls.tsv"), "w") as f:
        f.write("haplotype\tspecies\tpath\tref_start\tref_end\tref_span\tpath_start\tpath_end\tbreakpoint_near_gap_or_end\n")
        for c in sorted(calls, key=lambda c: (c["rs"], c["hap"])):
            f.write(f"{c['hap']}\t{c['sp']}\t{c['path']}\t{c['rs']}\t{c['re']}\t{c['re'] - c['rs']}\t{c['qs']}\t{c['qe']}\t"
                    f"{'yes' if c['near_gap'] else 'no'}\n")
    with open(os.path.join(D, "inversion_events.tsv"), "w") as f:
        f.write("event\tref_start\tref_end\tspan_bp\tptarmigan_start\tptarmigan_end\tin_41-44Mb_segment\t"
                + "\t".join(f"carriers_{s}\tcalled_{s}\tfreq_{s}" for s in SPECIES)
                + "\tclass\tclean_carrier_calls\tevidence\tcarrier_haplotypes\n")
        for e in out_ev:
            f.write(f"{e['id']}\t{e['rs']}\t{e['re']}\t{e['span']}\t{e['pt'][0] if e['pt'][0] is not None else 'NA'}\t"
                    f"{e['pt'][1] if e['pt'][1] is not None else 'NA'}\t{e['in_seg'] or 'NA'}\t"
                    + "\t".join(f"{e['k'][s]}\t{e['n'][s]}\t{e['fr'][s]:.2f}" for s in SPECIES)
                    + f"\t{e['cls']}\t{e['clean']}\t{e['evid']}\t{','.join(e['carriers'])}\n")
    # validation candidates: classified events that separate species, plus any clean event >= 500 kb
    with open(os.path.join(D, "validation_candidates.tsv"), "w") as f:
        f.write("event\tclass\tevidence\tcarrier_path\tcarrier_start\tcarrier_end\tbreakpoint_near_gap_or_end\t"
                "noncarrier_example\tref_start\tref_end\n")
        for e in out_ev:
            if not (e["cls"].startswith(("diagnostic", "fixed")) or (e["evid"] != "breakpoints_at_gaps_only" and e["span"] >= 500000)):
                continue
            non = sorted(h for h in zhaps if h not in e["carriers"])
            for c in sorted(e["calls"], key=lambda c: c["near_gap"]):
                f.write(f"{e['id']}\t{e['cls']}\t{e['evid']}\t{c['path']}\t{c['qs']}\t{c['qe']}\t"
                        f"{'yes' if c['near_gap'] else 'no'}\t{non[0] if non else 'NA'}\t{e['rs']}\t{e['re']}\n")

    cnt = defaultdict(int)
    for e in out_ev:
        cnt[e["cls"].split("_")[0] if not e["cls"].startswith("fixed") else "fixed"] += 1
    sep = [e for e in out_ev if e["cls"].startswith(("diagnostic", "fixed"))]
    lines = [f"Reference Z path: {ref} ({lens.get(ref, 0):,} bp)",
             "Z haplotypes per species: " + ", ".join(f"{s} {n_sp[s]}" for s in SPECIES)
             + f" (classification needs >= {a.min_called_frac:.0%} called per species)",
             f"Inversion calls >= {a.min_inv:,} bp: {len(calls)} in {len({c['hap'] for c in calls})} haplotypes; "
             f"{sum(1 for c in calls if not c['near_gap'])} with both breakpoints away from gaps and path ends",
             f"Events: {len(out_ev)} ({', '.join(f'{k} {v}' for k, v in sorted(cnt.items()))})",
             f"Events separating species (diagnostic or fixed between a pair): {len(sep)}"]
    for e in sep:
        pt = (f"; ptarmigan {min(e['pt']) / 1e6:.2f}-{max(e['pt']) / 1e6:.2f} Mb" if e["pt"][0] is not None and e["pt"][1] is not None else "")
        lines.append(f"  {e['id']} ref {e['rs'] / 1e6:.2f}-{e['re'] / 1e6:.2f} Mb ({e['span'] / 1e3:,.0f} kb){pt}: {e['cls']}; "
                     + ", ".join(f"{s} {e['k'][s]}/{e['n'][s]}" for s in SPECIES)
                     + f"; {e['evid']}" + (f"; overlaps 41-44 Mb segment: {e['in_seg']}" if e["in_seg"] else ""))
    if bl is None:
        lines.append("(no ref_to_ptarmigan.paf: ptarmigan coordinates and segment overlap not computed)")
    open(os.path.join(D, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
