#!/usr/bin/env python3
"""Summarize pixy pi/dxy/fst windows (step 07) for STGR, GRPC and LEPC.

All values are recomputed from pixy's raw counts (count_diffs / count_comparisons),
which is pixy's recommended way to combine windows:
  dxy      = sum(diffs between species) / sum(comparisons between species)
  pi       = sum(diffs within species)  / sum(comparisons within species)
  Fst      = 1 - mean(pi1, pi2) / dxy            (Hudson, as in pixy --fst_type hudson)
  da       = dxy - mean(pi1, pi2)                 (net divergence)

Outputs (in --out):
  windows.tsv           every window x species pair (dxy, Fst, da, pi1, pi2, sites)
  per_chromosome.tsv    the same per chromosome (chr_Z_males = Z with males only)
  groups.tsv            autosomes, macro (chr_1-chr_8), micro (chr_9+), Z, Z males only
  z_vs_autosomes.tsv    Z / autosome ratios with window-bootstrap 95% CIs
  pi_by_species.tsv     pi per species: autosomes, Z, Z males only; Z/A ratio (males-only Z)
                        with bootstrap CI, against the neutral expectations 0.75 (equal
                        breeding sex ratio) and 1.125 (maximum, few breeding males)
  genome_scan.{png,pdf} dxy and Fst along the genome for each species pair
  summary.txt           the main numbers in words
"""
import argparse, glob, os, re
import numpy as np
import pandas as pd

SP = ["STGR", "GRPC", "LEPC"]
PAIRS = [("GRPC", "LEPC"), ("LEPC", "STGR"), ("GRPC", "STGR")]
MACRO = {f"chr_{i}" for i in range(1, 9)}


def chrom_key(c):
    c = c.replace("_males", "")
    m = re.fullmatch(r"chr_(\d+)", c)
    return (0, int(m.group(1))) if m else (1, {"chr_Z": 0, "chr_W": 1}.get(c, 2))


def load(pixy_dir, stat, males):
    files = [f for f in glob.glob(os.path.join(pixy_dir, f"*_{stat}.txt")) if ("_males_" in f) == males]
    if not files:
        return pd.DataFrame()
    df = pd.concat((pd.read_csv(f, sep="\t", na_values=["NA"]) for f in files), ignore_index=True)
    if males:
        df["chromosome"] = df["chromosome"] + "_males"
    return df


def windows(pixy_dir, males):
    pi, dxy = load(pixy_dir, "pi", males), load(pixy_dir, "dxy", males)
    if dxy.empty:
        return pd.DataFrame()
    key = ["chromosome", "window_pos_1", "window_pos_2"]
    pi = pi[key + ["pop", "no_sites", "count_diffs", "count_comparisons"]]
    out = []
    for a, b in PAIRS:
        d = dxy[((dxy.pop1 == a) & (dxy.pop2 == b)) | ((dxy.pop1 == b) & (dxy.pop2 == a))]
        d = d[key + ["no_sites", "count_diffs", "count_comparisons"]].rename(
            columns={"count_diffs": "dxy_d", "count_comparisons": "dxy_c", "no_sites": "sites"})
        for sp, tag in ((a, "p1"), (b, "p2")):
            p = pi[pi["pop"] == sp].drop(columns=["pop", "no_sites"]).rename(
                columns={"count_diffs": f"{tag}_d", "count_comparisons": f"{tag}_c"})
            d = d.merge(p, on=key, how="left")
        d.insert(0, "pair", f"{a}-{b}")
        out.append(d)
    w = pd.concat(out, ignore_index=True).fillna({c: 0 for c in ["dxy_d", "dxy_c", "p1_d", "p1_c", "p2_d", "p2_c", "sites"]})
    return w


def stats(d_d, d_c, p1_d, p1_c, p2_d, p2_c):
    with np.errstate(divide="ignore", invalid="ignore"):
        dxy = d_d / d_c
        pi1, pi2 = p1_d / p1_c, p2_d / p2_c
        piw = (pi1 + pi2) / 2
        return dxy, 1 - piw / dxy, dxy - piw, pi1, pi2


COLS = ["dxy_d", "dxy_c", "p1_d", "p1_c", "p2_d", "p2_c"]


def agg(df):
    s = df[COLS].sum()
    dxy, fst, da, pi1, pi2 = stats(*[s[c] for c in COLS])
    return pd.Series({"windows": len(df), "sites": int(df["sites"].sum()), "dxy": dxy, "fst_hudson": fst,
                      "da": da, "pi_1": pi1, "pi_2": pi2})


def group_of(c):
    if c == "chr_Z_males":
        return ["Z_males_only"]
    if c == "chr_Z":
        return ["Z"]
    if c == "chr_W":
        return []
    return ["autosomes", "macro_chr1-8" if c in MACRO else "micro_chr9+"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pixy", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-sites", type=int, default=5000)
    ap.add_argument("--boot", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(1)

    w = pd.concat([windows(a.pixy, False), windows(a.pixy, True)], ignore_index=True)
    if w.empty:
        raise SystemExit("ERROR: no pixy output found")
    w["chr_order"] = w["chromosome"].map(chrom_key)
    w = w.sort_values(["pair", "chr_order", "window_pos_1"]).drop(columns="chr_order")
    dxy, fst, da, pi1, pi2 = stats(*[w[c].to_numpy(float) for c in COLS])
    w = w.assign(dxy=dxy, fst_hudson=fst, da=da, pi_1=pi1, pi_2=pi2)
    w.to_csv(os.path.join(a.out, "windows.tsv"), sep="\t", index=False, float_format="%.6g")

    # per chromosome
    pc = w.groupby(["pair", "chromosome"], sort=False).apply(agg).reset_index()
    pc.to_csv(os.path.join(a.out, "per_chromosome.tsv"), sep="\t", index=False, float_format="%.6g")

    # groups
    rows = []
    for (pair, chrom), d in w.groupby(["pair", "chromosome"], sort=False):
        for g in group_of(chrom):
            rows.append(d.assign(group=g))
    gw = pd.concat(rows, ignore_index=True)
    grp = gw.groupby(["pair", "group"]).apply(agg).reset_index()
    grp.to_csv(os.path.join(a.out, "groups.tsv"), sep="\t", index=False, float_format="%.6g")

    # Z vs autosomes: ratio of aggregated values, bootstrap over windows within each group
    zr = []
    for pair in [f"{x}-{y}" for x, y in PAIRS]:
        for zgroup in ("Z", "Z_males_only"):
            Z = gw[(gw.pair == pair) & (gw.group == zgroup) & (gw.sites >= a.min_sites)][COLS].to_numpy(float)
            A = gw[(gw.pair == pair) & (gw.group == "autosomes") & (gw.sites >= a.min_sites)][COLS].to_numpy(float)
            if len(Z) < 10 or len(A) < 10:
                continue
            obs = [np.array(stats(*x.sum(0))) for x in (Z, A)]
            bz = np.stack([np.array(stats(*Z[rng.integers(0, len(Z), len(Z))].sum(0))) for _ in range(a.boot)])
            ba = np.stack([np.array(stats(*A[rng.integers(0, len(A), len(A))].sum(0))) for _ in range(a.boot)])
            sp1, sp2 = pair.split("-")
            for i, name in enumerate(["dxy", "fst_hudson", "da", f"pi_{sp1}", f"pi_{sp2}"]):
                with np.errstate(divide="ignore", invalid="ignore"):
                    r = bz[:, i] / ba[:, i]
                lo, hi = np.nanpercentile(r, [2.5, 97.5])
                zr.append({"pair": pair, "z_group": zgroup, "stat": name, "Z": obs[0][i], "autosomes": obs[1][i],
                           "ratio_Z_A": obs[0][i] / obs[1][i], "ci95_lo": lo, "ci95_hi": hi,
                           "Z_windows": len(Z), "A_windows": len(A)})
    zr = pd.DataFrame(zr)
    zr.to_csv(os.path.join(a.out, "z_vs_autosomes.tsv"), sep="\t", index=False, float_format="%.6g")

    # pi per species; Z/A uses the males-only Z (females are hemizygous)
    pis = []
    for males in (False, True):
        p = load(a.pixy, "pi", males)
        if not p.empty:
            pis.append(p)
    pi = pd.concat(pis, ignore_index=True)
    pi = pi[(pi.no_sites >= a.min_sites) & (pi.chromosome != "chr_W")]
    pi["group"] = pi.chromosome.map(lambda c: "Z_males_only" if c == "chr_Z_males" else "Z" if c == "chr_Z" else "autosomes")
    prow = []
    for sp in SP:
        d = pi[pi["pop"] == sp]
        obs = {g: d[d.group == g][["count_diffs", "count_comparisons"]].to_numpy(float) for g in ("autosomes", "Z", "Z_males_only")}
        val = {g: (x[:, 0].sum() / x[:, 1].sum() if len(x) else np.nan) for g, x in obs.items()}
        r = {"species": sp, "pi_autosomes": val["autosomes"], "pi_Z_all": val["Z"], "pi_Z_males": val["Z_males_only"],
             "A_windows": len(obs["autosomes"]), "Z_males_windows": len(obs["Z_males_only"])}
        Z, A = obs["Z_males_only"], obs["autosomes"]
        if len(Z) >= 10 and len(A) >= 10:
            bz = np.array([Z[i].sum(0) for i in (rng.integers(0, len(Z), len(Z)) for _ in range(a.boot))])
            ba = np.array([A[i].sum(0) for i in (rng.integers(0, len(A), len(A)) for _ in range(a.boot))])
            rb = (bz[:, 0] / bz[:, 1]) / (ba[:, 0] / ba[:, 1])
            r["Z_A_ratio"] = val["Z_males_only"] / val["autosomes"]
            r["ci95_lo"], r["ci95_hi"] = np.nanpercentile(rb, [2.5, 97.5])
            r["boot_frac_below_0.75"] = float(np.mean(rb < 0.75))
            r["vs_neutral"] = ("below 0.75" if r["ci95_hi"] < 0.75 else "above 1.125" if r["ci95_lo"] > 1.125
                               else "within 0.75-1.125" if r["ci95_lo"] >= 0.75 and r["ci95_hi"] <= 1.125
                               else "overlaps 0.75" if r["ci95_lo"] < 0.75 else "overlaps 1.125")
        prow.append(r)
    ps = pd.DataFrame(prow)
    ps.to_csv(os.path.join(a.out, "pi_by_species.tsv"), sep="\t", index=False, float_format="%.6g")

    # words
    with open(os.path.join(a.out, "summary.txt"), "w") as fo:
        fo.write(f"pixy summary: {w[~w.chromosome.str.endswith('_males')].groupby('pair').sites.sum().iloc[0]:,} "
                 f"callable sites per pair (genome, all birds); windows with >= {a.min_sites} sites used for ratios\n\n")
        for pair in [f"{x}-{y}" for x, y in PAIRS]:
            g = grp[grp.pair == pair].set_index("group")
            sp1, sp2 = pair.split("-")
            fo.write(f"{pair}\n")
            for gname in ("autosomes", "macro_chr1-8", "micro_chr9+", "Z", "Z_males_only"):
                if gname in g.index:
                    r = g.loc[gname]
                    fo.write(f"  {gname:<14} dxy {r.dxy:.5f}  Fst {r.fst_hudson:.3f}  da {r.da:.5f}  "
                             f"pi {sp1} {r.pi_1:.5f}  pi {sp2} {r.pi_2:.5f}  ({int(r.windows)} windows)\n")
            for _, r in zr[(zr.pair == pair) & (zr.stat.isin(["dxy", "fst_hudson", "da"]))].iterrows():
                fo.write(f"  Z/A {r.stat:<10} ({r.z_group}): {r.ratio_Z_A:.2f} (95% CI {r.ci95_lo:.2f}-{r.ci95_hi:.2f})\n")
            fo.write("\n")
        fo.write("Nucleotide diversity per species (Z/A uses the males-only Z; neutral expectation 0.75 with an\n"
                 "equal breeding sex ratio, up to 1.125 when few males breed, as on leks)\n")
        for _, r in ps.iterrows():
            line = (f"  {r.species}: pi autosomes {r.pi_autosomes:.5f}  Z all birds {r.pi_Z_all:.5f}  "
                    f"Z males {r.pi_Z_males:.5f}")
            if "Z_A_ratio" in r and pd.notna(r.get("Z_A_ratio")):
                line += (f"  Z/A {r.Z_A_ratio:.2f} (95% CI {r.ci95_lo:.2f}-{r.ci95_hi:.2f}; {r.vs_neutral}; "
                         f"{int(r.Z_males_windows)} Z windows)")
            fo.write(line + "\n")
        fo.write("  (GRPC males-only Z = one bird, F5597: two Z copies, so its Z pi is noisy)\n")
    print(open(os.path.join(a.out, "summary.txt")).read())

    # genome scan
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g = w[(~w.chromosome.str.endswith("_males")) & (w.chromosome != "chr_W") & (w.sites >= a.min_sites)].copy()
    chroms = sorted(g.chromosome.unique(), key=chrom_key)
    lens = g.groupby("chromosome").window_pos_2.max()
    off, x0 = {}, 0
    for c in chroms:
        off[c] = x0
        x0 += lens[c] + 5e6
    g["x"] = g.window_pos_1.astype(float) + g.chromosome.map(off)
    fig, axes = plt.subplots(6, 1, figsize=(12, 11), sharex=True)
    for i, (sa, sb) in enumerate(PAIRS):
        d = g[g.pair == f"{sa}-{sb}"]
        for j, (col, lab) in enumerate((("dxy", "dxy"), ("fst_hudson", "Fst"))):
            ax = axes[2 * i + j]
            for k, c in enumerate(chroms):
                dc = d[d.chromosome == c]
                col_c = "#c0392b" if c == "chr_Z" else ("#4d4d4d" if k % 2 == 0 else "#a6a6a6")
                ax.scatter(dc.x / 1e6, dc[col], s=1.5, color=col_c, rasterized=True, linewidths=0)
            ax.set_ylabel(f"{lab}\n{sa}–{sb}", fontsize=8)
            ax.tick_params(labelsize=7)
            if col == "fst_hudson":
                ax.set_ylim(min(0, np.nanmin(d[col])), 1)
    ticks = [off[c] / 1e6 + lens[c] / 2e6 for c in chroms]
    labels = [c.replace("chr_", "") for c in chroms]
    axes[-1].set_xticks([t for t, l in zip(ticks, labels) if l in {"1", "2", "3", "4", "5", "6", "7", "8", "10", "15", "20", "Z"}])
    axes[-1].set_xticklabels([l for l in labels if l in {"1", "2", "3", "4", "5", "6", "7", "8", "10", "15", "20", "Z"}], fontsize=7)
    axes[-1].set_xlabel("Chromosome (reference: pangenome reference paths)", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"genome_scan.{ext}"), dpi=300)
    print(f"wrote {a.out}/genome_scan.png")


if __name__ == "__main__":
    main()
