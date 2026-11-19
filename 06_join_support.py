#!/usr/bin/env python3
# =============================================================================
# 11_join_support.py — Hi-C support for every RagTag join
#
# RagTag builds each pseudo-chromosome by concatenating yahs scaffolds in the
# order the rock ptarmigan reference implies, inserting a 100 bp gap (gap type
# "align_genus") at each join. Those joins are assertions the REFERENCE makes;
# the grouse Hi-C data has not been consulted. This script asks, for every one
# of them, whether the Hi-C agrees.
#
# METHOD
#   One streaming pass over the Hi-C alignment builds a binned contact matrix
#   over the chr_* sequences (default 100 kb bins). Everything else is computed
#   from that matrix, so each haplotype is read exactly once.
#
#   For a join at coordinate J, with flanking windows of W bins either side:
#       obs = contacts between the left window and the right window
#       exp = what those same bin pairs would carry if the sequence were
#             continuous, read off a distance-decay curve P(d) built ONLY from
#             bin pairs lying inside a single scaffold (never across a join)
#       bg  = inter-chromosomal contact density over the same number of bin
#             pairs, i.e. what two unlinked sequences would show
#
#       support = (obs - bg) / (exp - bg)
#
#   ~1.0 means the Hi-C supports the join as strongly as genuinely contiguous
#   sequence at the same separation. ~0.0 means the two sides look no more
#   linked than sequences on different chromosomes.
#
# WHY THE CONTROLS MATTER
#   A threshold on its own is worth little, so the script measures its own
#   null and alternative from the same data:
#     - POSITIVE controls: pseudo-joins placed at interior positions of large
#       intact scaffolds. These are contiguous by construction and should score
#       near 1.0. If they do not, the statistic is miscalibrated and the verdicts
#       are not trustworthy — the script says so rather than reporting anyway.
#     - NEGATIVE control: the inter-chromosomal background, which is subtracted.
#   The control distribution is written out and summarised, so the separation
#   between real joins and controls is visible rather than assumed.
#
# INPUT (per haplotype)
#   --agp        ragtag/<sample>.<hap>/ragtag.scaffold.agp
#   --rename-map ragtag/<sample>.<hap>/<sample>.<hap>.rename_map.tsv
#   --fai        final/<SPECIES>_<SAMPLE>_<HAP>.pseudo_chr.fasta.fai
#   SAM records on stdin (samtools view of the Hi-C alignment to the FINAL
#   assembly, which is what qc/<sample>.<hap>.hic2final.cram holds).
#
# OUTPUT
#   <out-prefix>.joins.tsv     one row per RagTag join, with the verdict
#   <out-prefix>.controls.tsv  the positive-control distribution
#   a summary to stdout
# =============================================================================

import argparse
import os
import sys

import numpy as np


# ----------------------------------------------------------------- parsing --
def read_fai(path, chrom_only=True):
    seqs = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            if chrom_only and not f[0].startswith("chr_"):
                continue
            seqs.append((f[0], int(f[1])))
    return seqs


def read_rename_map(path):
    m = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    m[f[0]] = f[1]
    return m


def read_agp(path, rename):
    """
    -> {final_seq_name: [(comp_name, start, end), ...]} for W components only,
    in object coordinates, with the object renamed to its chr_* name.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[4] != "W":
                continue
            obj = rename.get(f[0], f[0])
            out.setdefault(obj, []).append((f[5], int(f[1]), int(f[2])))
    for k in out:
        out[k].sort(key=lambda t: t[1])
    return out


# ------------------------------------------------------------ bin geometry --
class Bins:
    def __init__(self, seqs, bin_size):
        self.bin_size = bin_size
        self.offset = {}
        self.nbins_of = {}
        self.seq_of_bin = []
        n = 0
        for name, length in seqs:
            nb = max(1, -(-length // bin_size))      # ceil
            self.offset[name] = n
            self.nbins_of[name] = nb
            self.seq_of_bin.extend([name] * nb)
            n += nb
        self.n = n

    def bin_of(self, seq, pos):
        off = self.offset.get(seq)
        if off is None:
            return None
        b = off + (pos - 1) // self.bin_size
        return b if b < off + self.nbins_of[seq] else None


def build_matrix(bins, stream, min_mapq):
    """One pass over SAM records -> symmetric binned contact matrix."""
    M = np.zeros((bins.n, bins.n), dtype=np.int32)
    kept = skipped = 0
    for line in stream:
        if line[0] == "@":
            continue
        f = line.split("\t", 9)
        if len(f) < 9:
            continue
        flag = int(f[1])
        # primary, both mates mapped, count each pair once via READ1
        if flag & 0x904 or flag & 0x8 or not (flag & 0x40):
            continue
        if int(f[4]) < min_mapq:
            continue
        rname, pos, rnext, pnext = f[2], int(f[3]), f[6], int(f[7])
        if rnext == "=":
            rnext = rname
        b1 = bins.bin_of(rname, pos)
        b2 = bins.bin_of(rnext, pnext)
        if b1 is None or b2 is None:
            skipped += 1
            continue
        M[b1, b2] += 1
        if b1 != b2:
            M[b2, b1] += 1
        kept += 1
    return M, kept, skipped


# -------------------------------------------------------------- statistics --
def decay_curve(M, bins, components, max_d):
    """
    Mean contacts per bin pair as a function of separation, using ONLY bin
    pairs that fall inside one scaffold. Never crosses a join, so it is an
    independent expectation for what a contiguous join should look like.
    """
    tot = np.zeros(max_d + 1, dtype=np.float64)
    cnt = np.zeros(max_d + 1, dtype=np.float64)
    for seq, comps in components.items():
        off = bins.offset.get(seq)
        if off is None:
            continue
        for _, start, end in comps:
            b0 = off + (start - 1) // bins.bin_size
            b1 = off + (end - 1) // bins.bin_size
            if b1 - b0 < 2:
                continue
            sub = M[b0:b1 + 1, b0:b1 + 1]
            k = sub.shape[0]
            for d in range(1, min(max_d, k - 1) + 1):
                diag = np.diagonal(sub, offset=d)
                tot[d] += diag.sum()
                cnt[d] += diag.size
    with np.errstate(invalid="ignore", divide="ignore"):
        curve = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    return curve, cnt


def trans_density(M, bins):
    """Mean contacts per bin pair between DIFFERENT chr_* sequences."""
    total = 0.0
    pairs = 0.0
    names = list(bins.offset)
    for i, a in enumerate(names):
        oa, na = bins.offset[a], bins.nbins_of[a]
        for b in names[i + 1:]:
            ob, nb = bins.offset[b], bins.nbins_of[b]
            total += M[oa:oa + na, ob:ob + nb].sum()
            pairs += na * nb
    return (total / pairs) if pairs else 0.0


def score_boundary(M, curve, bg, b_left_end, w):
    """
    Contacts across a boundary sitting at bin b_left_end, over w bins either
    side. Returns (obs, exp, bg_total, ratio).

    The boundary bin ITSELF is excluded from both windows. RagTag gaps are only
    100 bp, so the bin holding the end of the left component almost always also
    holds the start of the right one; counting it on the left leaks that
    component's own short-range contacts across the junction and inflates the
    score. Leaving it out costs one bin of signal and removes the bias — in
    testing it was the difference between a false join scoring 0.37 and 0.00.
    """
    l0, l1 = b_left_end - w, b_left_end          # [b-w, b-1]
    r0, r1 = b_left_end + 1, b_left_end + 1 + w  # [b+1, b+w]
    if l0 < 0 or r1 > M.shape[0]:
        return None
    obs = float(M[l0:l1, r0:r1].sum())
    exp = 0.0
    for i in range(l0, l1):
        d = np.arange(r0, r1) - i
        vals = curve[np.clip(d, 0, len(curve) - 1)]
        exp += np.nansum(vals)
    npairs = w * w
    bg_total = bg * npairs
    denom = exp - bg_total
    ratio = (obs - bg_total) / denom if denom > 1e-9 else float("nan")
    return obs, exp, bg_total, ratio


# -------------------------------------------------------------------- main --
def main():
    p = argparse.ArgumentParser(description="Hi-C support for RagTag joins.")
    p.add_argument("--agp", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--rename-map", default="")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--assembly", default="", help="label for the output rows")
    p.add_argument("--bin-size", type=int, default=100000)
    p.add_argument("--window-bins", type=int, default=10,
                   help="flank size in bins either side of a join (default 10)")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-component", type=int, default=1000000,
                   help="only score joins where BOTH flanking components are at "
                        "least this long; shorter ones carry too little signal")
    p.add_argument("--supported-at", type=float, default=0.40,
                   help="call a join supported at this FRACTION of the positive-"
                        "control median support (default 0.40)")
    p.add_argument("--unsupported-at", type=float, default=0.15,
                   help="call a join unsupported at or below this fraction of "
                        "the control median (default 0.15)")
    p.add_argument("--max-controls", type=int, default=400)
    args = p.parse_args()

    label = args.assembly or os.path.basename(args.out_prefix)
    seqs = read_fai(args.fai)
    if not seqs:
        sys.exit(f"ERROR: no chr_* sequences in {args.fai}")
    rename = read_rename_map(args.rename_map)
    comps = read_agp(args.agp, rename)

    bins = Bins(seqs, args.bin_size)
    print(f">>> {label}: {len(seqs)} chromosomes, {bins.n} bins of "
          f"{args.bin_size // 1000} kb", flush=True)

    M, kept, skipped = build_matrix(bins, sys.stdin, args.min_mapq)
    print(f"    {kept:,} read pairs binned ({skipped:,} outside chr_* sequences)",
          flush=True)
    if kept == 0:
        sys.exit("ERROR: no usable read pairs — check the alignment and --min-mapq")

    w = args.window_bins
    curve, curve_n = decay_curve(M, bins, comps, max_d=4 * w + 2)
    bg = trans_density(M, bins)
    print(f"    inter-chromosomal background: {bg:.2f} contacts per bin pair",
          flush=True)
    # Robust sparsity guard: judge the decay curve over the whole window range
    # rather than at a single separation, which is noisy and can be empty.
    near = curve[1:w + 1]
    if np.all(np.isnan(near)) or np.nanmedian(near) <= bg:
        sys.exit("ERROR: the distance-decay curve over the first %d bins does "
                 "not rise above the inter-chromosomal background (median "
                 "%.2f vs %.2f). The contact data are too sparse at this bin "
                 "size; try a larger --bin-size." %
                 (w, float(np.nanmedian(near)) if not np.all(np.isnan(near)) else float('nan'), bg))

    # ---------------- positive controls: interior points of large scaffolds --
    controls = []
    rng = np.random.default_rng(0)
    for seq, cl in comps.items():
        off = bins.offset.get(seq)
        if off is None:
            continue
        for name, start, end in cl:
            if end - start + 1 < args.min_component * 2:
                continue
            b0 = off + (start - 1) // args.bin_size
            b1 = off + (end - 1) // args.bin_size
            lo, hi = b0 + 2 * w, b1 - 2 * w
            if hi <= lo:
                continue
            for b in rng.choice(np.arange(lo, hi), size=min(8, hi - lo),
                                replace=False):
                r = score_boundary(M, curve, bg, int(b), w)
                if r and np.isfinite(r[3]):
                    controls.append((seq, name, int(b), r[3]))
    if len(controls) > args.max_controls:
        idx = rng.choice(len(controls), args.max_controls, replace=False)
        controls = [controls[i] for i in idx]

    with open(f"{args.out_prefix}.controls.tsv", "w") as fh:
        fh.write("assembly\tsequence\tscaffold\tbin\tsupport\n")
        for s, n, b, r in controls:
            fh.write(f"{label}\t{s}\t{n}\t{b}\t{r:.3f}\n")

    if not controls:
        sys.exit("ERROR: no positive controls could be placed — scaffolds are "
                 "too short relative to --window-bins. Verdicts withheld.")
    cvals = np.array([c[3] for c in controls])
    c_med = float(np.median(cvals))
    c_lo = float(np.percentile(cvals, 5))
    print(f"    positive controls (contiguous by construction): n={len(cvals)}, "
          f"median support {c_med:.2f}, 5th pct {c_lo:.2f}", flush=True)

    # Thresholds are fractions of the control median rather than absolutes, so
    # the call adapts to each library's contact density instead of relying on a
    # number tuned elsewhere.
    sup_cut = args.supported_at * c_med
    uns_cut = args.unsupported_at * c_med
    print(f"    thresholds: supported >= {sup_cut:.2f}, unsupported <= {uns_cut:.2f} "
          f"({args.supported_at:g} and {args.unsupported_at:g} x control median)",
          flush=True)

    # The controls are contiguous by construction, so any that fall below the
    # thresholds are misclassifications. Reporting that rate turns the verdicts
    # from an assertion into a measurement the reader can weigh.
    fp_unsup = float(np.mean(cvals <= uns_cut))
    fp_notsup = float(np.mean(cvals < sup_cut))
    print(f"    control misclassification at these thresholds: "
          f"{100 * fp_unsup:.1f}% of known-contiguous positions would be called "
          f"'unsupported', {100 * fp_notsup:.1f}% would fail 'supported'",
          flush=True)
    if fp_unsup > 0.05:
        print("    WARNING: more than 5% of contiguous controls score as "
              "unsupported, so individual 'unsupported' calls here carry a "
              "real false-positive rate. Treat them as candidates to inspect, "
              "not conclusions.", flush=True)

    calibrated = 0.60 <= c_med <= 1.40
    if not calibrated:
        print("    WARNING: controls do not centre near 1.0, so the support "
              "statistic is miscalibrated here. Rows are written with verdict "
              "'uncalibrated' and should not be acted on.", flush=True)

    # ------------------------------------------------- score the real joins --
    rows = []
    for seq, cl in comps.items():
        off = bins.offset.get(seq)
        if off is None or len(cl) < 2:
            continue
        for i in range(len(cl) - 1):
            lname, lstart, lend = cl[i]
            rname_, rstart, rend = cl[i + 1]
            lsize, rsize = lend - lstart + 1, rend - rstart + 1
            small = lsize < args.min_component or rsize < args.min_component
            b = off + (lend - 1) // args.bin_size
            res = score_boundary(M, curve, bg, b, w)
            if res is None:
                continue
            obs, exp, bgt, ratio = res
            if not calibrated:
                verdict = "uncalibrated"
            elif small:
                verdict = "not_scored_small"
            elif not np.isfinite(ratio):
                verdict = "indeterminate"
            elif ratio >= sup_cut:
                verdict = "supported"
            elif ratio <= uns_cut:
                verdict = "unsupported"
            else:
                verdict = "ambiguous"
            rows.append((label, seq, lend, lname, rname_, lsize, rsize,
                         obs, exp, bgt, ratio, verdict))

    with open(f"{args.out_prefix}.joins.tsv", "w") as fh:
        fh.write("assembly\tsequence\tjoin_pos\tleft_scaffold\tright_scaffold\t"
                 "left_bp\tright_bp\tobs\texp\tbackground\tsupport\tverdict\n")
        for r in rows:
            fh.write("%s\t%s\t%d\t%s\t%s\t%d\t%d\t%.0f\t%.1f\t%.1f\t%.3f\t%s\n" % r)

    scored = [r for r in rows if r[11] in ("supported", "unsupported", "ambiguous")]
    print(f"    {len(rows)} RagTag joins; {len(scored)} scored "
          f"(both flanks >= {args.min_component // 1000} kb)")
    for v in ("supported", "ambiguous", "unsupported", "not_scored_small",
              "indeterminate", "uncalibrated"):
        n = sum(1 for r in rows if r[11] == v)
        if n:
            print(f"      {v:<18} {n}")
    if scored:
        worst = sorted(scored, key=lambda r: r[10])[:8]
        print("    weakest joins:")
        for r in worst:
            print(f"      {r[1]:<9} at {r[2]:>11,}  support {r[10]:+.2f}  "
                  f"{r[3]} | {r[4]}  ({r[5]/1e6:.1f} / {r[6]/1e6:.1f} Mb)  {r[11]}")
    print(f">>> wrote {args.out_prefix}.joins.tsv and .controls.tsv")


if __name__ == "__main__":
    main()
