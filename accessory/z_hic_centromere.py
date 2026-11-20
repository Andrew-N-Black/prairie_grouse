#!/usr/bin/env python3
"""Where is the Z centromere? Hi-C contacts along chr_Z of the final (ptarmigan-ordered) assemblies.
Run by z_hic_centromere.sh; reads `samtools view` of the Hi-C CRAM restricted to chr_Z on stdin.

Three centromere signals, each per bin along the Z:
  1. trans fraction   share of a bin's contacts that go to other chromosomes. Centromeres (and telomeres)
                      of different chromosomes cluster in many nuclei, so centromeric bins gain trans contacts.
  2. Rabl symmetry    in a Rabl configuration the two arms fold back at the centromere, so positions
                      equidistant from it contact each other: for each candidate position c, the mean
                      observed/expected contact between c-k and c+k (k = KMIN..KMAX) peaks at the centromere
                      (the anti-diagonal "X" in the heatmap).
  3. insulation dip   contacts crossing each position at 1-5 Mb separation (observed/expected); arm boundaries
                      and the centromere often insulate.
Output: <out>.png (heatmap of log2 observed/expected + the three profiles, with the marked segment) and
<out>.tsv (profiles per bin)."""
import argparse, sys, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr_Z")
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--bin", type=int, default=250_000)
    ap.add_argument("--minq", type=int, default=10)
    ap.add_argument("--kmin", type=float, default=2e6)
    ap.add_argument("--kmax", type=float, default=15e6)
    ap.add_argument("--mark", default="", help="start-end (bp) of the low-diversity segment on this assembly")
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    nb = a.length // a.bin + 1
    M = np.zeros((nb, nb)); trans = np.zeros(nb); cis = np.zeros(nb)
    n = 0
    for line in sys.stdin:
        t = line.split("\t", 9)
        flag, rname, pos, mapq, rnext, pnext = int(t[1]), t[2], int(t[3]), int(t[4]), t[6], int(t[7])
        if rname != a.chrom or mapq < a.minq or not flag & 1 or flag & 0x904:
            continue
        if not flag & 64:                       # count each pair once (read 1)
            continue
        i = pos // a.bin
        if rnext == "=" or rnext == a.chrom:
            j = pnext // a.bin
            if 0 <= j < nb and abs(pnext - pos) >= 1000:     # drop self-ligation / short-range pairs
                M[i, j] += 1; M[j, i] += 1; cis[i] += 1; cis[j] += 1
        elif rnext != "*":
            trans[i] += 1
        n += 1
    sys.stderr.write(f"[z_hic] {n} read-1 records; {int(M.sum() / 2)} cis pairs, {int(trans.sum())} trans\n")
    # observed / expected by distance
    E = np.zeros_like(M)
    for d in range(nb):
        v = np.diagonal(M, d).mean() if nb - d > 0 else 0
        idx = np.arange(nb - d)
        E[idx, idx + d] = v; E[idx + d, idx] = v
    with np.errstate(divide="ignore", invalid="ignore"):
        OE = np.where(E > 0, M / E, np.nan)
    tf = trans / np.maximum(trans + cis, 1)
    kmin, kmax = int(a.kmin // a.bin), int(a.kmax // a.bin)
    rabl = np.full(nb, np.nan)
    for c in range(nb):
        ks = [k for k in range(kmin, kmax + 1) if c - k >= 0 and c + k < nb]
        if len(ks) >= max(3, (kmax - kmin) // 3):
            rabl[c] = np.nanmean([np.log2(OE[c - k, c + k]) if OE[c - k, c + k] > 0 else np.nan for k in ks])
    d1, d2 = int(1e6 // a.bin), int(5e6 // a.bin)
    ins = np.full(nb, np.nan)
    for x in range(nb):
        vals = [OE[i, j] for i in range(max(0, x - d2), x) for j in range(x, min(nb, x + d2))
                if d1 <= j - i <= d2 and np.isfinite(OE[i, j])]
        if len(vals) > 10:
            ins[x] = np.log2(np.mean(vals))
    with open(a.out + ".tsv", "w") as o:
        o.write("bin_start\tcis\ttrans\ttrans_fraction\trabl_log2OE\tinsulation_log2OE\n")
        for b in range(nb):
            o.write(f"{b * a.bin}\t{int(cis[b])}\t{int(trans[b])}\t{tf[b]:.4f}\t{rabl[b]:.4f}\t{ins[b]:.4f}\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(nb) * a.bin / 1e6
    fig = plt.figure(figsize=(6.5, 9.5))
    gs = fig.add_gridspec(4, 1, height_ratios=[3.2, 1, 1, 1], hspace=0.35)
    ax = fig.add_subplot(gs[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log2(OE)
    im = ax.imshow(L, cmap="RdBu_r", vmin=-2, vmax=2, extent=(0, x[-1] + a.bin / 1e6, x[-1] + a.bin / 1e6, 0),
                   interpolation="nearest")
    fig.colorbar(im, ax=ax, fraction=0.04, label="log2 observed/expected")
    ax.set_title(a.title or a.chrom, loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel(f"{a.chrom} (Mb)"); ax.set_ylabel("Mb")
    mk = [float(v) / 1e6 for v in a.mark.split("-")] if a.mark else None
    if mk:
        for v in mk:
            ax.axvline(v, color="black", lw=0.6, ls=(0, (3, 2))); ax.axhline(v, color="black", lw=0.6, ls=(0, (3, 2)))
    for k, (y, lab) in enumerate([(tf, "Trans fraction"), (rabl, "Rabl symmetry\n(log2 O/E)"),
                                  (ins, "Insulation\n(log2 O/E, 1-5 Mb)")]):
        axp = fig.add_subplot(gs[k + 1], sharex=None)
        axp.plot(x, y, color="black", lw=0.9)
        if mk:
            axp.axvspan(mk[0], mk[1], color="#BBBBBB", alpha=0.5, lw=0)
        axp.set_ylabel(lab, fontsize=8); axp.set_xlim(0, x[-1])
        axp.spines[["top", "right"]].set_visible(False)
    axp.set_xlabel(f"{a.chrom} (Mb); grey = low-diversity segment")
    fig.savefig(a.out + ".png", dpi=200, bbox_inches="tight")
    best = int(np.nanargmax(rabl)) if np.isfinite(rabl).any() else -1
    msg = f"[z_hic] Rabl peak at {best * a.bin / 1e6:.2f} Mb; trans-fraction peak at {np.argmax(tf) * a.bin / 1e6:.2f} Mb"
    if mk:
        msg += f"; marked segment {mk[0]:.2f}-{mk[1]:.2f} Mb"
    print(msg)


if __name__ == "__main__":
    main()
