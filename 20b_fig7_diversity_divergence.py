#!/usr/bin/env python3
"""Figure 7 (step 07): FST along the genome, pi and divergence by chromosome class, and the Z in
rock ptarmigan coordinates. Run after 07_pixy_summary.py.
usage: python3 07_fig7_diversity_divergence.py <summary dir> [<out prefix>]
needs <summary dir>/{windows.tsv, groups.tsv, ref_to_ptarmigan.paf}"""
import sys, os
SUM = sys.argv[1]; OUTP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(SUM, "Fig7_diversity_divergence")
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False})
COL = {"GRPC": "#DAA520", "LEPC": "#A52A2A", "STGR": "#000000"}
PCOL = {"GRPC-LEPC": "#1F6F8B", "LEPC-STGR": "#7F7F7F"}
PLAB = {"GRPC-LEPC": "LEPC–GRPC", "LEPC-STGR": "LEPC–STGR"}
C = ["dxy_d", "dxy_c", "p1_d", "p1_c", "p2_d", "p2_c"]
W = pd.read_csv(os.path.join(SUM, "windows.tsv"), sep="\t"); W = W[W.kept]
G = pd.read_csv(os.path.join(SUM, "groups.tsv"), sep="\t")
# lift Z window midpoints to ptarmigan coordinates through the reference-Z -> ptarmigan-Z alignment
paf = pd.read_csv(os.path.join(SUM, "ref_to_ptarmigan.paf"), sep="\t", header=None, usecols=range(12))
paf = paf[(paf[11] >= 5) & ((paf[3] - paf[2]) >= 10000)].sort_values(10, ascending=False)
def lift(p):
    for r in paf.itertuples(index=False):
        if r[2] <= p < r[3]:
            f = (p - r[2]) / (r[3] - r[2])
            return r[7] + f * (r[8] - r[7]) if r[4] == "+" else r[8] - f * (r[8] - r[7])
    return np.nan
Z = W[W.chromosome.isin(["chr_Z", "chr_Z_males"])].copy()
mids = {p: lift(p + 25000) for p in Z.window_pos_1.unique()}
Z["pt"] = Z.window_pos_1.map(mids)

def stats(d):
    dxy = d.dxy_d / d.dxy_c; p1 = d.p1_d / d.p1_c; p2 = d.p2_d / d.p2_c
    return d.assign(dxy=dxy, pi1=p1, pi2=p2, fst=1 - ((p1 + p2) / 2) / dxy)

def ckey(c):
    m = re.fullmatch(r"chr_(\d+)", c); return int(m.group(1)) if m else 99

BIN = 500_000
a = W[(W.chromosome.str.match(r"chr_(\d+|Z)$")) & (W.sites >= 5000)].copy()
a["bin"] = (a.window_pos_1 - 1) // BIN
a = stats(a.groupby(["pair", "chromosome", "bin"])[C + ["sites"]].sum().reset_index())
a = a[a.sites >= 100_000]
chroms = sorted(a.chromosome.unique(), key=ckey)
L = a.groupby("chromosome").bin.max() + 1
off, x = {}, 0
for c in chroms:
    off[c] = x; x += L[c] * BIN + 6e6
a["x"] = (a.bin * BIN + BIN / 2 + a.chromosome.map(off)) / 1e6

fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(4, 6, height_ratios=[1, 1, 1.15, 1.15], hspace=0.55, wspace=1.1)

# A: Fst along the genome
for r, pair in enumerate(["GRPC-LEPC", "LEPC-STGR"]):
    ax = fig.add_subplot(gs[r, :])
    d = a[a.pair == pair]
    for k, c in enumerate(chroms):
        dc = d[d.chromosome == c]
        if c == "chr_Z":
            ax.axvspan(off[c] / 1e6 - 2, (off[c] + L[c] * BIN) / 1e6 + 2, color="#F3E3E3", lw=0, zorder=0)
        ax.scatter(dc.x, dc.fst, s=1.6, lw=0, color=PCOL[pair] if k % 2 == 0 or c == "chr_Z" else "#B9C9D0" if pair == "GRPC-LEPC" else "#C4C4C4", rasterized=True, zorder=2)
    ax.set_ylim(-0.02, 1.02); ax.set_xlim(-3, x / 1e6)
    ax.set_ylabel(f"$F_{{ST}}$\n{PLAB[pair]}")
    big = {"chr_" + str(i) for i in range(1, 9)} | {"chr_10", "chr_15", "chr_20", "chr_Z"}
    ax.set_xticks([(off[c] + L[c] * BIN / 2) / 1e6 for c in chroms if c in big])
    ax.set_xticklabels([c.replace("chr_", "") for c in chroms if c in big] if r == 1 else [])
    if r == 0:
        ax.text(-0.065, 1.08, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for grp, lab in (("autosomes", "autosomes"), ("Z_males_only" if pair == "GRPC-LEPC" else "Z", "Z")):
        pass
    gA = G[(G.pair == pair) & (G.group == "autosomes")].fst_hudson.iloc[0]
    ax.axhline(gA, color="#333333", lw=0.6, ls=(0, (3, 2)), xmax=off["chr_Z"] / x * 0.995)
ax.set_xlabel("Chromosome (500-kb windows)")

# B: pi by chromosome class
order = [("macro_chr1-8", "Chr 1–8"), ("micro_chr9+", "Chr 9–37"), ("Z_males_only", "Z (males)")]
def pi_of(sp, grp):
    for pair in G.pair.unique():
        s1, s2 = pair.split("-")
        r = G[(G.pair == pair) & (G.group == grp)]
        if s1 == sp: return r.pi_1.iloc[0]
        if s2 == sp: return r.pi_2.iloc[0]
ax = fig.add_subplot(gs[2:, 0:2])
for j, sp in enumerate(["LEPC", "GRPC", "STGR"]):
    ys = [pi_of(sp, g) * 1e3 for g, _ in order]
    ax.plot(np.arange(3) + (j - 1) * 0.15, ys, "o", ms=4.5, color=COL[sp], mec="white", mew=0.6, label=sp, zorder=3)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in order])
ax.set_ylabel("Nucleotide diversity π (×10$^{-3}$)"); ax.set_ylim(0, 5.5); ax.set_xlim(-0.5, 2.5)
ax.legend(frameon=False, loc="upper left", handletextpad=0.1, borderaxespad=0.2)
ax.text(-0.33, 1.03, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")

# C: dxy and da by chromosome class
ax = fig.add_subplot(gs[2:, 2:4])
order2 = [("macro_chr1-8", "Chr 1–8"), ("micro_chr9+", "Chr 9–37"), ("Z", "Z")]
for j, pair in enumerate(["GRPC-LEPC", "LEPC-STGR"]):
    r = [G[(G.pair == pair) & (G.group == g)].iloc[0] for g, _ in order2]
    xs = np.arange(3) + (j - 0.5) * 0.28
    ax.bar(xs, [q.dxy * 1e3 for q in r], width=0.26, color=PCOL[pair], alpha=0.35, lw=0)
    ax.bar(xs, [q.da * 1e3 for q in r], width=0.26, color=PCOL[pair], lw=0)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in order2])
ax.set_ylabel("Divergence (×10$^{-3}$)")
h = [Line2D([], [], marker="s", ls="", ms=6, color=PCOL[p]) for p in PCOL] + \
    [Line2D([], [], marker="s", ls="", ms=6, color="#888888", alpha=0.35), Line2D([], [], marker="s", ls="", ms=6, color="#888888")]
ax.legend(h, [PLAB[p] for p in PCOL] + ["$d_{XY}$ (full bar)", "$d_a$ (net)"], frameon=False, loc="upper left", borderaxespad=0.2, handletextpad=0.1)
ax.set_ylim(0, 12.5)
ax.text(-0.33, 1.03, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")

# D: Z close-up in ptarmigan coordinates
z = Z.dropna(subset=["pt"]).copy()
z["bin"] = (z.pt // BIN)
def zagg(df):
    s = df.groupby(["pair", "bin"])[C].sum().reset_index()
    s = stats(s); s["mb"] = (s.bin * BIN + BIN / 2) / 1e6; return s
zm, za = zagg(z[z.chromosome == "chr_Z_males"]), zagg(z[z.chromosome == "chr_Z"])
ax1 = fig.add_subplot(gs[2, 4:]); ax2 = fig.add_subplot(gs[3, 4:], sharex=ax1)
for ax in (ax1, ax2):
    ax.axvspan(41.11, 44.30, color="#E8E8E8", lw=0, zorder=0)
for sp in ["LEPC", "GRPC", "STGR"]:
    for pair in zm.pair.unique():
        s1, s2 = pair.split("-")
        if sp in (s1, s2):
            d = zm[zm.pair == pair].sort_values("mb")
            ax1.plot(d.mb, (d.pi1 if sp == s1 else d.pi2) * 1e3, color=COL[sp], lw=0.9, label=sp); break
ax1.set_ylabel("π (×10$^{-3}$)\nZ, males"); ax1.set_ylim(0, 4.5)
ax1.tick_params(labelbottom=False)
for pair in ["GRPC-LEPC", "LEPC-STGR"]:
    d = za[za.pair == pair].sort_values("mb")
    ax2.plot(d.mb, d.dxy * 1e3, color=PCOL[pair], lw=0.9, label=PLAB[pair])
ax2.set_ylabel("$d_{XY}$ (×10$^{-3}$)\nZ"); ax2.set_ylim(0, None)
ax2.set_xlabel("Z position (Mb, rock ptarmigan)")
ax2.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0, 1.22), borderaxespad=0, ncol=2, handlelength=1.2)
ax2.set_ylim(0, 14)
ax1.text(-0.2, 1.06, "D", transform=ax1.transAxes, fontsize=10, fontweight="bold")
ax1.text(42.7, ax1.get_ylim()[1] * 0.97, "41–44 Mb", ha="center", va="top", fontsize=6, color="#555555")
for ext in ("png", "pdf"):
    fig.savefig(f"{OUTP}.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {OUTP}.png/.pdf")
