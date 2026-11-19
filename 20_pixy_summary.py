#!/usr/bin/env python3
"""Summarize pixy pi/dxy/fst windows (step 07) for STGR, GRPC and LEPC.

All values are recomputed from pixy's raw counts (count_diffs / count_comparisons),
which is pixy's recommended way to combine windows:
  dxy      = sum(diffs between species) / sum(comparisons between species)
  pi       = sum(diffs within species)  / sum(comparisons within species)
  Fst      = 1 - mean(pi1, pi2) / dxy            (Hudson, as in pixy --fst_type hudson)
  da       = dxy - mean(pi1, pi2)                 (net divergence)

Window filters (applied to every statistic and every species pair alike):
  - fewer than --min-sites callable sites
  - divergence outliers: pi of any species or dxy of any pair > --max-div (default 0.02,
    about five times the genome-wide mean). Such windows (chromosome ends, the first
    ~3 Mb of the reference Z, the smallest microchromosomes) have high pi even within
    males and are taken to be collapsed repeats or paralogous read mapping. They are
    listed in excluded_windows.tsv; set --max-div 1 to keep them.

Outputs (in --out):
  windows.tsv           every window x species pair, with a 'kept' column
  excluded_windows.tsv  windows removed by the filters, and why
  per_chromosome.tsv    per chromosome (chr_Z_males = Z with males only)
  groups.tsv            autosomes, macro (chr_1-chr_8), micro (chr_9+), Z, Z males only
  z_vs_autosomes.tsv    Z / autosome ratios with block-bootstrap 95% CIs (1-Mb blocks)
  pi_by_species.tsv     pi per species: autosomes, macro, Z, Z males only; Z/A ratio
                        (males-only Z) with bootstrap CI, against the neutral range for a
                        ZW system: 0.75 with an equal breeding sex ratio, falling to
                        9/16 = 0.5625 when few males breed (lek mating), because the Z
                        spends two-thirds of its time in males:
                        Ne_Z/Ne_A = 9(Nm+Nf) / (4(2Nm+4Nf))
  genome_scan.{png,pdf} dxy and Fst along the genome for each species pair
  summary.txt           the main numbers in words

Input: --pixy <dir of pixy *_pi/_dxy.txt>, or --windows <a windows.tsv written by this
script> to recompute everything with other filters without the raw pixy files.
"""
import argparse, glob, os, re
import numpy as np
import pandas as pd

SP = ["STGR", "GRPC", "LEPC"]
PAIRS = [("GRPC", "LEPC"), ("LEPC", "STGR"), ("GRPC", "STGR")]
MACRO = {f"chr_{i}" for i in range(1, 9)}
COLS = ["dxy_d", "dxy_c", "p1_d", "p1_c", "p2_d", "p2_c"]
KEY = ["chromosome", "window_pos_1", "window_pos_2"]


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


def windows_from_pixy(pixy_dir, males):
    pi, dxy = load(pixy_dir, "pi", males), load(pixy_dir, "dxy", males)
    if dxy.empty:
        return pd.DataFrame()
    pi = pi[KEY + ["pop", "count_diffs", "count_comparisons"]]
    out = []
    for a, b in PAIRS:
        d = dxy[((dxy.pop1 == a) & (dxy.pop2 == b)) | ((dxy.pop1 == b) & (dxy.pop2 == a))]
        d = d[KEY + ["no_sites", "count_diffs", "count_comparisons"]].rename(
            columns={"count_diffs": "dxy_d", "count_comparisons": "dxy_c", "no_sites": "sites"})
        for sp, tag in ((a, "p1"), (b, "p2")):
            p = pi[pi["pop"] == sp].drop(columns=["pop"]).rename(
                columns={"count_diffs": f"{tag}_d", "count_comparisons": f"{tag}_c"})
            d = d.merge(p, on=KEY, how="left")
        d.insert(0, "pair", f"{a}-{b}")
        out.append(d)
    return pd.concat(out, ignore_index=True).fillna({c: 0 for c in COLS + ["sites"]})


def stats(d_d, d_c, p1_d, p1_c, p2_d, p2_c):
    with np.errstate(divide="ignore", invalid="ignore"):
        dxy = d_d / d_c
        pi1, pi2 = p1_d / p1_c, p2_d / p2_c
        piw = (pi1 + pi2) / 2
        return dxy, 1 - piw / dxy, dxy - piw, pi1, pi2


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


def blocks(df, cols, size):
    """Sum window counts into consecutive blocks of `size` bp per chromosome (bootstrap units,
    so that linked neighbouring windows are resampled together)."""
    b = df.assign(_b=(df.window_pos_1 - 1) // size)
    return b.groupby(["chromosome", "_b"])[cols].sum().to_numpy(float)


def boot_ratio(Z, A, f, n, rng):
    """Bootstrap windows within Z and within A; f maps summed counts to a statistic."""
    bz = np.array([f(Z[rng.integers(0, len(Z), len(Z))].sum(0)) for _ in range(n)])
    ba = np.array([f(A[rng.integers(0, len(A), len(A))].sum(0)) for _ in range(n)])
    with np.errstate(divide="ignore", invalid="ignore"):
        return bz / ba


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pixy")
    src.add_argument("--windows")
    ap.add_argument("--samples", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-sites", type=int, default=5000)
    ap.add_argument("--max-div", type=float, default=0.02)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--block", type=int, default=1000000, help="bootstrap block size (bp); windows are summed into blocks")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(1)

    # ---------------- windows and filters ----------------
    if a.pixy:
        w = pd.concat([windows_from_pixy(a.pixy, False), windows_from_pixy(a.pixy, True)], ignore_index=True)
        if w.empty:
            raise SystemExit("ERROR: no pixy output found")
    else:
        w = pd.read_csv(a.windows, sep="\t")
    w = w[["pair"] + KEY + ["sites"] + COLS].copy()
    w["chr_order"] = w["chromosome"].map(chrom_key)
    w = w.sort_values(["pair", "chr_order", "window_pos_1"]).drop(columns="chr_order")
    dxy, fst, da, pi1, pi2 = stats(*[w[c].to_numpy(float) for c in COLS])
    w = w.assign(dxy=dxy, fst_hudson=fst, da=da, pi_1=pi1, pi_2=pi2)
    w["base"] = w.chromosome.str.replace("_males", "", regex=False)
    low = w.groupby(["base", "window_pos_1"]).sites.transform("min") < a.min_sites
    hi_row = (w[["dxy", "pi_1", "pi_2"]].max(axis=1) > a.max_div)
    hi = hi_row.groupby([w.base, w.window_pos_1]).transform("any")   # any pair, all-birds or males-only
    w["kept"] = ~(low | hi)
    w["excluded_because"] = np.where(low, f"sites<{a.min_sites}", np.where(hi, f"pi_or_dxy>{a.max_div}", ""))
    w.drop(columns="base").to_csv(os.path.join(a.out, "windows.tsv"), sep="\t", index=False, float_format="%.6g")
    ex = w[~w.kept].drop_duplicates(KEY)[KEY + ["excluded_because"]]
    ex.to_csv(os.path.join(a.out, "excluded_windows.tsv"), sep="\t", index=False)
    allw = w.drop_duplicates(KEY)
    exs = (allw.assign(g=allw.chromosome.map(lambda c: "Z" if c.startswith("chr_Z") else "W" if c == "chr_W" else
                                               "macro" if c in MACRO else "micro"))
               .groupby(["g", "excluded_because"]).size().unstack(fill_value=0))
    k = w[w.kept].copy()

    # ---------------- per chromosome and groups ----------------
    pc = k.groupby(["pair", "chromosome"], sort=False).apply(agg).reset_index()
    pc.to_csv(os.path.join(a.out, "per_chromosome.tsv"), sep="\t", index=False, float_format="%.6g")
    gw = pd.concat([d.assign(group=g) for (_, c), d in k.groupby(["pair", "chromosome"], sort=False)
                    for g in group_of(c)], ignore_index=True)
    grp = gw.groupby(["pair", "group"]).apply(agg).reset_index()
    grp.to_csv(os.path.join(a.out, "groups.tsv"), sep="\t", index=False, float_format="%.6g")

    # ---------------- Z vs autosomes, per pair ----------------
    zr = []
    names = ["dxy", "fst_hudson", "da"]
    for pair in [f"{x}-{y}" for x, y in PAIRS]:
        A = blocks(gw[(gw.pair == pair) & (gw.group == "autosomes")], COLS, a.block)
        for zgroup in ("Z", "Z_males_only"):
            Z = blocks(gw[(gw.pair == pair) & (gw.group == zgroup)], COLS, a.block)
            nwin = int(((gw.pair == pair) & (gw.group == zgroup)).sum())
            if len(Z) < 10 or len(A) < 10:
                continue
            oz, oa = np.array(stats(*Z.sum(0))), np.array(stats(*A.sum(0)))
            rb = boot_ratio(Z, A, lambda s: np.array(stats(*s))[:3], a.boot, rng)
            for i, nm in enumerate(names):
                lo, hi_ = np.nanpercentile(rb[:, i], [2.5, 97.5])
                zr.append({"pair": pair, "z_group": zgroup, "stat": nm, "Z": oz[i], "autosomes": oa[i],
                           "ratio_Z_A": oz[i] / oa[i], "ci95_lo": lo, "ci95_hi": hi_,
                           "Z_windows": nwin, "Z_blocks": len(Z), "A_blocks": len(A)})
    zr = pd.DataFrame(zr)
    zr.to_csv(os.path.join(a.out, "z_vs_autosomes.tsv"), sep="\t", index=False, float_format="%.6g")

    # ---------------- pi per species ----------------
    def sp_counts(sp, chroms):
        for pair in [f"{x}-{y}" for x, y in PAIRS]:
            s1, s2 = pair.split("-")
            if sp in (s1, s2):
                tag = "p1" if sp == s1 else "p2"
                d = k[(k.pair == pair) & (k.chromosome.isin(chroms))]
                return blocks(d, [f"{tag}_d", f"{tag}_c"], a.block)
    auto = [c for c in k.chromosome.unique() if re.fullmatch(r"chr_\d+", c)]
    prow = []
    for sp in SP:
        A, M = sp_counts(sp, auto), sp_counts(sp, list(MACRO))
        Za, Zm = sp_counts(sp, ["chr_Z"]), sp_counts(sp, ["chr_Z_males"])
        v = lambda x: x[:, 0].sum() / x[:, 1].sum() if len(x) else np.nan
        r = {"species": sp, "pi_autosomes": v(A), "pi_macro_chr1-8": v(M), "pi_Z_all": v(Za), "pi_Z_males": v(Zm),
             "A_blocks": len(A), "Z_males_blocks": len(Zm)}
        if len(Zm) >= 10:
            for ref, R in (("A", A), ("macro", M)):
                rb = boot_ratio(Zm, R, lambda s: s[0] / s[1], a.boot, rng)
                lo, hi_ = np.nanpercentile(rb, [2.5, 97.5])
                r[f"Z_{ref}_ratio"] = v(Zm) / v(R)
                r[f"Z_{ref}_ci95_lo"], r[f"Z_{ref}_ci95_hi"] = lo, hi_
                if ref == "A":
                    r["boot_frac_below_0.5625"] = float(np.mean(rb < 0.5625))
                    r["vs_neutral"] = ("below 0.5625 (lower than any neutral sex ratio)" if hi_ < 0.5625
                                       else "above 0.75" if lo > 0.75
                                       else "within 0.5625-0.75" if lo >= 0.5625 and hi_ <= 0.75
                                       else "overlaps 0.5625" if lo < 0.5625 else "overlaps 0.75")
        prow.append(r)
    ps = pd.DataFrame(prow)
    ps.to_csv(os.path.join(a.out, "pi_by_species.tsv"), sep="\t", index=False, float_format="%.6g")

    # ---------------- words ----------------
    with open(os.path.join(a.out, "summary.txt"), "w") as fo:
        fo.write(f"pixy summary. Windows kept: {int(allw.kept.sum())} of {len(allw)}; excluded by filter:\n")
        fo.write("  " + exs.to_string().replace("\n", "\n  ") + "\n\n")
        for pair in [f"{x}-{y}" for x, y in PAIRS]:
            g = grp[grp.pair == pair].set_index("group")
            sp1, sp2 = pair.split("-")
            fo.write(f"{pair}\n")
            for gname in ("autosomes", "macro_chr1-8", "micro_chr9+", "Z", "Z_males_only"):
                if gname in g.index:
                    r = g.loc[gname]
                    fo.write(f"  {gname:<14} dxy {r.dxy:.5f}  Fst {r.fst_hudson:.3f}  da {r.da:.5f}  "
                             f"pi {sp1} {r.pi_1:.5f}  pi {sp2} {r.pi_2:.5f}  ({int(r.windows)} windows)\n")
            for _, r in zr[zr.pair == pair].iterrows():
                fo.write(f"  Z/A {r.stat:<10} ({r.z_group}): {r.ratio_Z_A:.2f} (95% CI {r.ci95_lo:.2f}-{r.ci95_hi:.2f})\n")
            fo.write("\n")
        fo.write("Nucleotide diversity per species, Z = males only. Neutral Z/A expectation for a ZW system:\n"
                 "0.75 with an equal breeding sex ratio, falling to 0.5625 when few males breed (leks)\n")
        for _, r in ps.iterrows():
            line = (f"  {r.species}: pi autosomes {r.pi_autosomes:.5f}  chr1-8 {r['pi_macro_chr1-8']:.5f}  "
                    f"Z all birds {r.pi_Z_all:.5f}  Z males {r.pi_Z_males:.5f}")
            if pd.notna(r.get("Z_A_ratio")):
                line += (f"\n      Z/A {r.Z_A_ratio:.2f} (95% CI {r.Z_A_ci95_lo:.2f}-{r.Z_A_ci95_hi:.2f}; {r.vs_neutral}); "
                         f"Z/chr1-8 {r.Z_macro_ratio:.2f} ({r.Z_macro_ci95_lo:.2f}-{r.Z_macro_ci95_hi:.2f})")
            fo.write(line + "\n")
        fo.write("  (GRPC males-only Z = one bird, F5597: two Z copies, so its Z pi is noisy)\n")
    print(open(os.path.join(a.out, "summary.txt")).read())

    # ---------------- genome scan ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g = k[(~k.chromosome.str.endswith("_males")) & (k.chromosome != "chr_W")].copy()
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
            for kk, c in enumerate(chroms):
                dc = d[d.chromosome == c]
                col_c = "#c0392b" if c == "chr_Z" else ("#4d4d4d" if kk % 2 == 0 else "#a6a6a6")
                ax.scatter(dc.x / 1e6, dc[col], s=1.5, color=col_c, rasterized=True, linewidths=0)
            ax.set_ylabel(f"{lab}\n{sa}–{sb}", fontsize=8)
            ax.tick_params(labelsize=7)
            if col == "fst_hudson":
                ax.set_ylim(min(0, np.nanmin(d[col])), 1)
    show = {"1", "2", "3", "4", "5", "6", "7", "8", "10", "15", "20", "Z"}
    ticks = [(off[c] + lens[c] / 2) / 1e6 for c in chroms]
    labels = [c.replace("chr_", "") for c in chroms]
    axes[-1].set_xticks([t for t, l in zip(ticks, labels) if l in show])
    axes[-1].set_xticklabels([l for l in labels if l in show], fontsize=7)
    axes[-1].set_xlabel("Chromosome (reference: pangenome reference paths)", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"genome_scan.{ext}"), dpi=300)
    print(f"wrote {a.out}/genome_scan.png")


if __name__ == "__main__":
    main()
