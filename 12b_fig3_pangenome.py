#!/usr/bin/env python3
"""Step 12b: Figure 3A (pangenome growth and core by species) and 3B (repeat and gene content
per haplotype), drawn from the step-12 summary tables.
usage: python3 12b_fig3_pangenome.py <dir with pangenome_histgrowth_{all,GRPC,LEPC,STGR}.tsv,
                                      genome_length_by_haplotype.tsv, per_haplotype_summary.tsv> [out dir]"""
import sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
COL = {"GRPC": "#DAA520", "LEPC": "#A52A2A", "STGR": "#000000", "All": "#8C8C8C"}
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 1.1, "xtick.major.width": 1.1,
                     "ytick.major.width": 1.1})

# ---------------- A: growth ----------------
fig, ax = plt.subplots(figsize=(6.5, 4.5))
for key, lab in [("all", "All"), ("GRPC", "GRPC"), ("LEPC", "LEPC"), ("STGR", "STGR")]:
    d = pd.read_csv(f"{D}/pangenome_histgrowth_{key}.tsv", sep="\t", skiprows=5, header=None, names=["n", "pan", "core"])
    d = d[d.n > 0]
    lw = 2.2 if lab == "All" else 1.6
    ax.plot(d.n, d.pan / 1e6, color=COL[lab], lw=lw, solid_capstyle="round")
    ax.plot(d.n, d.core / 1e6, color=COL[lab], lw=lw, ls=(0, (4, 2)))
ax.set_xlabel("Number of haplotypes")
ax.set_ylabel("Sequence (Mb)")
ax.set_xlim(0, 47); ax.set_ylim(820, 1170)
ax.set_xticks([1, 10, 20, 30, 40, 46])
h = [Line2D([], [], color=COL[k], lw=2.2 if k == "All" else 1.6) for k in ["All", "GRPC", "LEPC", "STGR"]]
l1 = ax.legend(h, ["All (n = 46)", "GRPC (n = 12)", "LEPC (n = 14)", "STGR (n = 20)"], title="Haplotypes",
               frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0), alignment="left")
ax.add_artist(l1)
h2 = [Line2D([], [], color="k", lw=1.6), Line2D([], [], color="k", lw=1.6, ls=(0, (4, 2)))]
ax.legend(h2, ["Pangenome", "Core"], title="Sequence", frameon=False, loc="upper left",
          bbox_to_anchor=(1.0, 0.45), alignment="left")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig3A_pangenome_growth.png", dpi=600)
fig.savefig(f"{OUT}/Fig3A_pangenome_growth.pdf")

# ---------------- B: repeats and genes ----------------
gl = pd.read_csv(f"{D}/genome_length_by_haplotype.tsv", sep="\t")
gl["chrom"] = gl.contig.str.replace(r"_p\d+$", "", regex=True)
z = gl[gl.chrom == "chr_Z"].groupby(["sample", "haplotype"]).path_length_bp.sum()
ph = pd.read_csv(f"{D}/per_haplotype_summary.tsv", sep="\t")
fem = {"F5503", "F5595", "F5596", "F5598", "F5599", "F5600"}
ph["zmb"] = [z.get((s, hh), 0) / 1e6 for s, hh in zip(ph["sample"], ph.haplotype)]
ph["cls"] = np.where(~ph["sample"].isin(fem), "male", np.where(ph.zmb > 60, "fZW", "fnoZ"))
MK = {"male": dict(marker="o", fill=True), "fZW": dict(marker="^", fill=True), "fnoZ": dict(marker="s", fill=False)}
order = ["GRPC", "LEPC", "STGR"]
rng = np.random.default_rng(7)

fig, axs = plt.subplots(1, 2, figsize=(6.5, 4.5))
for ax, col, lab in [(axs[0], "te_pct_masked", "Masked sequence (%)"), (axs[1], "genes_mapped", "Genes transferred")]:
    for i, sp in enumerate(order):
        s = ph[ph.species == sp]
        med = s.loc[s.cls == "male", col].median()        # median of male haplotypes
        ax.plot([i - 0.3, i + 0.3], [med, med], color="#555555", lw=1.6, zorder=2)
        for c, m in MK.items():
            t = s[s.cls == c]
            x = i + rng.uniform(-0.18, 0.18, len(t))
            ax.scatter(x, t[col], s=34, marker=m["marker"], linewidths=1.0, zorder=3,
                       facecolors=COL[sp] if m["fill"] else "white", edgecolors="black" if sp != "STGR" else "#666666")
    ax.set_xticks(range(3)); ax.set_xticklabels(order); ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel(lab)
axs[1].yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v:,.0f}"))
h = [Line2D([], [], ls="", marker="o", mfc="#999999", mec="black", ms=7),
     Line2D([], [], ls="", marker="^", mfc="#999999", mec="black", ms=7),
     Line2D([], [], ls="", marker="s", mfc="white", mec="black", ms=7)]
fig.legend(h, ["Male", "Female, Z + W", "Female, no Z"], frameon=False, loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(f"{OUT}/Fig3B_repeat_gene_content.png", dpi=600)
fig.savefig(f"{OUT}/Fig3B_repeat_gene_content.pdf")
print(ph.groupby(["species", "cls"]).size())
