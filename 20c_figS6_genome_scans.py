#!/usr/bin/env python3
"""Fig. S6 (step 20): FST, dxy and pi along the genome for each species pair.
For every pair three tracks: FST and dxy (chromosomes alternating between the colours of the two
species compared) and pi of each of the two species (species colours; alternate chromosomes on a
grey band). 500-kb windows summed from pixy counts; windows excluded by 20_pixy_summary.py are
left out. Autosomes and Z FST/dxy use all 23 birds; Z pi uses the 17 males.
usage: python3 20c_figS6_genome_scans.py <summary dir> [<out prefix>]
needs <summary dir>/{windows.tsv, groups.tsv}"""
import sys, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SUM = sys.argv[1]
OUTP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(SUM, "FigS6_genome_scans")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = {"GRPC": "#DAA520", "LEPC": "#A52A2A", "STGR": "#000000"}          # house species colours
PAIRS = ["GRPC-LEPC", "LEPC-STGR", "GRPC-STGR"]
PLAB = {"GRPC-LEPC": "LEPC–GRPC", "LEPC-STGR": "LEPC–STGR", "GRPC-STGR": "GRPC–STGR"}
C = ["dxy_d", "dxy_c", "p1_d", "p1_c", "p2_d", "p2_c"]
BIN = 500_000

W = pd.read_csv(os.path.join(SUM, "windows.tsv"), sep="\t")
W = W[W.kept]
G = pd.read_csv(os.path.join(SUM, "groups.tsv"), sep="\t")


def ckey(c):
    m = re.fullmatch(r"chr_(\d+)", c)
    return int(m.group(1)) if m else 99


def binned(df):
    df = df.assign(bin=(df.window_pos_1 - 1) // BIN)
    s = df.groupby(["pair", "chromosome", "bin"])[C + ["sites"]].sum().reset_index()
    s = s[s.sites >= 100_000]
    p1, p2 = s.p1_d / s.p1_c, s.p2_d / s.p2_c
    return s.assign(dxy=s.dxy_d / s.dxy_c, pi1=p1, pi2=p2, fst=1 - ((p1 + p2) / 2) / (s.dxy_d / s.dxy_c))


a = binned(W[W.chromosome.str.match(r"chr_(\d+|Z)$")])                     # all birds
zm = binned(W[W.chromosome == "chr_Z_males"]).assign(chromosome="chr_Z")    # Z, males only (for pi)
chroms = sorted(a.chromosome.unique(), key=ckey)
L = a.groupby("chromosome").bin.max() + 1
off, x = {}, 0
for c in chroms:
    off[c] = x
    x += L[c] * BIN + 6e6
for d in (a, zm):
    d["x"] = (d.bin * BIN + BIN / 2 + d.chromosome.map(off)) / 1e6

DXY_MAX = np.nanpercentile(a.dxy * 1e3, 99.8) * 1.05
PI_MAX = np.nanpercentile(pd.concat([a.pi1, a.pi2]) * 1e3, 99.8) * 1.05
big = {"chr_" + str(i) for i in range(1, 9)} | {"chr_10", "chr_15", "chr_20", "chr_Z"}

fig = plt.figure(figsize=(7.2, 11))
gsS = fig.add_gridspec(11, 1, height_ratios=[1, 1, 1, 0.35, 1, 1, 1, 0.35, 1, 1, 1], hspace=0.25)
rows = [0, 1, 2, 4, 5, 6, 8, 9, 10]                       # blank rows 3 and 7 separate the pairs
axes = []
for i in rows:
    axes.append(fig.add_subplot(gsS[i], sharex=axes[0] if axes else None))
for ax in axes[:-1]:
    ax.tick_params(labelbottom=False)
for p, pair in enumerate(PAIRS):
    s1, s2 = pair.split("-")
    spA, spB = ("LEPC", s1 if s1 != "LEPC" else s2) if "LEPC" in (s1, s2) else (s1, s2)
    d = a[a.pair == pair]
    dz = zm[zm.pair == pair]
    gA = G[(G.pair == pair) & (G.group == "autosomes")].iloc[0]
    for t, stat in enumerate(("fst", "dxy", "pi")):
        ax = axes[3 * p + t]
        for k, c in enumerate(chroms):
            if c == "chr_Z":
                ax.axvspan(off[c] / 1e6 - 2, (off[c] + L[c] * BIN) / 1e6 + 2, color="#F3E3E3", lw=0, zorder=0)
            elif stat == "pi" and k % 2:
                ax.axvspan(off[c] / 1e6, (off[c] + L[c] * BIN) / 1e6, color="#EFEFEF", lw=0, zorder=0)
            dc = d[d.chromosome == c]
            if stat == "fst":
                ax.scatter(dc.x, dc.fst, s=1.4, lw=0, color=COL[spA] if k % 2 == 0 else COL[spB], rasterized=True)
            elif stat == "dxy":
                ax.scatter(dc.x, (dc.dxy * 1e3).clip(upper=DXY_MAX), s=1.4, lw=0,
                           color=COL[spA] if k % 2 == 0 else COL[spB], rasterized=True)
            else:
                src = dz if c == "chr_Z" else dc
                for sp, col in ((s1, "pi1"), (s2, "pi2")):
                    ax.scatter(src.x, (src[col] * 1e3).clip(upper=PI_MAX), s=1.2, lw=0, color=COL[sp],
                               alpha=0.8, rasterized=True)
        if stat == "fst":
            ax.set_ylim(-0.02, 1.02)
            ax.set_ylabel(f"$F_{{ST}}$", fontsize=7)
            ax.axhline(gA.fst_hudson, color="#333333", lw=0.6, ls=(0, (3, 2)), xmax=off["chr_Z"] / x * 0.995)
            ax.set_title(PLAB[pair], loc="left", fontsize=8, fontweight="bold", pad=2)
        elif stat == "dxy":
            ax.set_ylim(0, DXY_MAX)
            ax.set_ylabel("$d_{XY}$\n(×10$^{-3}$)", fontsize=7)
            ax.axhline(gA.dxy * 1e3, color="#333333", lw=0.6, ls=(0, (3, 2)), xmax=off["chr_Z"] / x * 0.995)
        else:
            ax.set_ylim(0, PI_MAX)
            ax.set_ylabel("π\n(×10$^{-3}$)", fontsize=7)
            ax.legend([Line2D([], [], ls="", marker="o", ms=3, color=COL[s]) for s in (spA, spB)], [spA, spB],
                      frameon=False, loc="upper left", ncol=2, fontsize=6, handletextpad=0.1,
                      borderaxespad=0.1, columnspacing=0.8)
        ax.set_xlim(-3, x / 1e6)
        ax.tick_params(labelsize=6)
axes[-1].set_xticks([(off[c] + L[c] * BIN / 2) / 1e6 for c in chroms if c in big])
axes[-1].set_xticklabels([c.replace("chr_", "") for c in chroms if c in big])
axes[-1].set_xlabel("Chromosome (500-kb windows)")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUTP}.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {OUTP}.png/.pdf")
