#!/usr/bin/env python3
"""Fig. S11: reduced Z introgression between LEPC and GRPC (from hybrid_index.sh output).
  A  S vs H, autosomes;  B  S vs H, Z of males;  C  Z S vs autosomal S (females: one Z);
  D  autosomal and Z S of allopatric and range-edge LEPC vs distance to the GRPC range.
Also writes table S10 (one row per bird).
usage: python3 figS11_hybrid_index.py <hybrid_index.tsv> <panel_groups_lepc_grpc.tsv> <out prefix> [mito_panel.tsv]
                                      [qc_excluded_birds.tsv]
Birds listed in qc_excluded_birds.tsv (first column ID) are left out of the figure and table.
mito_panel.tsv (mito_panel.sh) adds the maternal lineage to table S10: LEPC (>= 90% LEPC alleles at the
diagnostic sites), GRPC/STGR (<= 10%) or unassigned."""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

hi, groups, outp = sys.argv[1:4]
mito = sys.argv[4] if len(sys.argv) > 4 else None
excl = set()
if len(sys.argv) > 5:
    excl = {l.split("\t")[0] for l in open(sys.argv[5]) if l.strip() and not l.startswith(("#", "ID\t"))}
COL = {"LEPC": "#A52A2A", "GRPC": "#DAA520"}
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.8})
d = pd.read_csv(hi, sep="\t")
d = d[d.pair == "LEPC-GRPC"]
g = pd.read_csv(groups, sep="\t")[["ID", "county", "km_GRPC_range_map"]]
w = d.pivot_table(index=["ID", "species", "grp", "context", "state", "sex"], columns="class",
                  values=["S", "S_lo", "S_hi", "H", "H_lo", "H_hi", "n_sites"]).reset_index()
w.columns = ["_".join(c).strip("_") for c in w.columns]
w = w.merge(g, on="ID", how="left")
w = w[w.species.isin(["LEPC", "GRPC"]) & ~w.ID.isin(excl)]
print(f"{len(excl)} birds excluded by QC")


def style(r):
    """marker, filled for a bird's group"""
    if r.grp.endswith("_allo"):
        return "o", True
    if r.grp == "LEPC_edge":
        return "s", False
    if r.grp == "GRPC_control":
        return "D", True
    return "^", False


fig, axs = plt.subplots(2, 2, figsize=(7.0, 6.4))
for ax, cls, tag, title in ((axs[0, 0], "autosomes", "A", "Autosomes"), (axs[0, 1], "Z", "B", "Z (males)")):
    ax.plot([0, 0.5, 1, 0], [0, 1, 0, 0], color="#999999", lw=0.7, zorder=0)
    t = w.dropna(subset=[f"S_{cls}", f"H_{cls}"])
    for _, r in t.sort_values("grp").iterrows():
        mk, fill = style(r)
        ax.scatter(r[f"S_{cls}"], r[f"H_{cls}"], s=13, marker=mk, linewidths=0.6, zorder=2,
                   facecolors=COL[r.species] if fill else "white", edgecolors=COL[r.species])
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Hybrid index, S (GRPC ancestry)"); ax.set_ylabel("Interspecific heterozygosity, H")
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold")
    ax.text(-0.2, 1.04, tag, transform=ax.transAxes, fontsize=11, fontweight="bold")
h = [Line2D([], [], ls="", marker="o", mfc=COL["LEPC"], mec=COL["LEPC"], ms=4.5),
     Line2D([], [], ls="", marker="s", mfc="white", mec=COL["LEPC"], ms=4.5),
     Line2D([], [], ls="", marker="^", mfc="white", mec=COL["LEPC"], ms=4.5),
     Line2D([], [], ls="", marker="o", mfc=COL["GRPC"], mec=COL["GRPC"], ms=4.5),
     Line2D([], [], ls="", marker="^", mfc="white", mec=COL["GRPC"], ms=4.5),
     Line2D([], [], ls="", marker="D", mfc=COL["GRPC"], mec=COL["GRPC"], ms=4)]
axs[0, 0].legend(h, ["LEPC reference", "LEPC range edge", "LEPC contact zone", "GRPC reference",
                     "GRPC contact zone", "GRPC South Dakota"], frameon=False, fontsize=6.5, loc="upper right",
                 handletextpad=0.2)
# C: Z vs autosomes
ax = axs[1, 0]
ax.plot([0, 1], [0, 1], color="#999999", lw=0.7, ls=(0, (3, 2)), zorder=0)
t = w.dropna(subset=["S_Z", "S_autosomes"])
for sx, mk in (("M", "o"), ("F", "s")):
    u = t[t.sex == sx]
    ax.scatter(u.S_autosomes, u.S_Z, s=13, marker=mk, c=[COL[s] for s in u.species], alpha=0.85, lw=0, zorder=2)
ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
ax.set_xlabel("Autosomal S"); ax.set_ylabel("Z S")
ax.legend([Line2D([], [], ls="", marker="o", mfc="#666666", mec="none", ms=4.5),
           Line2D([], [], ls="", marker="s", mfc="#666666", mec="none", ms=4.5)],
          ["Male (two Z)", "Female (one Z)"], frameon=False, fontsize=6.5, loc="upper left")
ax.text(-0.2, 1.04, "C", transform=ax.transAxes, fontsize=11, fontweight="bold")
# D: distance gradient in LEPC outside the contact zone
ax = axs[1, 1]
L = w[w.grp.isin(["LEPC_allo", "LEPC_edge"])].dropna(subset=["km_GRPC_range_map"])
rng = np.random.default_rng(3)
jit = lambda n: np.exp(rng.normal(0, 0.04, n))
ax.scatter(L.km_GRPC_range_map * jit(len(L)), L.S_autosomes, s=12, marker="o", color=COL["LEPC"], lw=0,
           alpha=0.8, label="Autosomes")
Lz = L.dropna(subset=["S_Z"])
ax.scatter(Lz.km_GRPC_range_map * jit(len(Lz)), Lz.S_Z, s=12, marker="o", facecolors="white",
           edgecolors="black", lw=0.6, label="Z")
ax.set_xscale("log"); ax.set_xlim(12, 800)
ax.set_xticks([20, 50, 100, 200, 500]); ax.set_xticklabels(["20", "50", "100", "200", "500"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("Distance to GRPC range (km)"); ax.set_ylabel("S (GRPC ancestry)")
ax.legend(frameon=False, fontsize=6.5, loc="upper right")
ax.text(-0.2, 1.04, "D", transform=ax.transAxes, fontsize=11, fontweight="bold")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{outp}.{ext}", dpi=600 if ext == "png" else None)
# table S10
cols = ["ID", "species", "grp", "context", "state", "county", "sex", "km_GRPC_range_map"]
for cls in ("autosomes", "Z"):
    cols += [f"{v}_{cls}" for v in ("n_sites", "S", "S_lo", "S_hi", "H", "H_lo", "H_hi")]
tab = w[cols].copy()
if mito:
    m = pd.read_csv(mito, sep="\t")[["ID", "mito_reads", "diag_sites_covered", "frac_LEPC", "sites_mixed"]]
    m["mt_lineage"] = np.select([m.frac_LEPC >= 0.9, m.frac_LEPC <= 0.1], ["LEPC", "GRPC/STGR"], "unassigned")
    tab = tab.merge(m, on="ID", how="left")
for c in tab.columns:
    if tab[c].dtype == float and not c.startswith(("n_sites", "km")):
        tab[c] = tab[c].round(4)
tab.sort_values(["species", "grp", "state", "county", "ID"]).to_csv(f"{outp}_tableS10.tsv", sep="\t", index=False)
print(f"wrote {outp}.png/pdf and {outp}_tableS10.tsv ({len(tab)} birds)")
