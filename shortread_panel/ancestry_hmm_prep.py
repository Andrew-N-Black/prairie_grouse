#!/usr/bin/env python3
"""Helpers for ancestry_hmm.sh (local ancestry with Ancestry_HMM; Corbett-Detig & Nielsen 2017).

  sites      ancestry-informative sites from hybrid_index.sh output: allele frequencies of the LEPC and GRPC
             reference groups re-estimated from the genotype likelihoods (red/*.npz) with the estimate-step groups
             (panel_groups_lepc_grpc.tsv) minus QC-excluded birds; sites with dAF >= --daf, one per --thin bp;
             major/minor bases from the ANGSD .mafs.gz files. Writes sites.tsv and sites.bed.
  count      major/minor read counts of one bird at the sites, from `samtools mpileup` output on stdin.
  build      Ancestry_HMM input + sample file for one admixed group and chromosome class (A or Z): reference
             panels as expected allele counts (frequency x 2n), Morgan distances from a uniform rate, bird read
             counts; ploidy 2, except females on the Z (1).
  summarize  posterior files -> expected GRPC ancestry per site and bird, 1-Mb windows averaged over birds,
             per-bird genome-wide estimates; table + figure.
"""
import argparse, glob, gzip, os, re, sys
import numpy as np


def em_freq(gl, n_iter=30):
    f = np.full(gl.shape[0], 0.2)
    for _ in range(n_iter):
        pri = np.stack([(1 - f) ** 2, 2 * f * (1 - f), f ** 2], axis=1)[:, None, :]
        post = gl * pri
        post /= post.sum(axis=2, keepdims=True) + 1e-300
        f = np.clip((post[:, :, 1] + 2 * post[:, :, 2]).mean(axis=1) / 2, 1e-4, 1 - 1e-4)
    return f


def read_tsv(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def cmd_sites(a):
    excl = {r["ID"] for r in read_tsv(a.exclude)} if a.exclude else set()
    grp = {r["ID"]: r["grp"] for r in read_tsv(a.groups) if r["ID"] not in excl}
    zset = set(a.z.split(","))
    ch, po, f0, f1, n0, n1 = [], [], [], [], [], []
    for fn in sorted(glob.glob(os.path.join(a.red, "*.npz"))):
        z = np.load(fn)
        if len(z["pos"]) == 0:
            continue
        ids = list(z["ids"])
        m0 = np.array([grp.get(b) == "LEPC_allo" for b in ids]); m1 = np.array([grp.get(b) == "GRPC_allo" for b in ids])
        gl = z["gl"].astype(np.float32)
        ch.append(z["chrom"]); po.append(z["pos"])
        f0.append(em_freq(gl[:, m0, :])); f1.append(em_freq(gl[:, m1, :]))
        n0.append(int(m0.sum())); n1.append(int(m1.sum()))
    ch = np.concatenate(ch); po = np.concatenate(po); f0 = np.concatenate(f0); f1 = np.concatenate(f1)
    n0, n1 = n0[0], n1[0]
    daf = np.abs(f0 - f1)
    keep = np.where(daf >= a.daf)[0]
    # thin: best site per bin
    key = np.array([f"{c}:{p // a.thin}" for c, p in zip(ch[keep], po[keep])])
    order = np.lexsort((-daf[keep], key))
    _, first = np.unique(key[order], return_index=True)
    sel = np.sort(keep[order[first]])
    # alleles from the ANGSD mafs files
    want = {(ch[i], int(po[i])) for i in sel}
    alle = {}
    for fn in glob.glob(os.path.join(a.mafs, "*.mafs.gz")):
        with gzip.open(fn, "rt") as fh:
            h = fh.readline().rstrip("\n").split("\t")
            ic, ip, ima, imi = h.index("chromo"), h.index("position"), h.index("major"), h.index("minor")
            for line in fh:
                t = line.split("\t")
                k = (t[ic], int(t[ip]))
                if k in want:
                    alle[k] = (t[ima], t[imi])
    sel = [i for i in sel if (ch[i], int(po[i])) in alle]
    sel.sort(key=lambda i: (ch[i], int(po[i])))
    with open(a.out + ".tsv", "w") as o, open(a.out + ".bed", "w") as b:
        o.write("chrom\tpos\tmajor\tminor\tclass\tf_LEPC\tf_GRPC\tn_LEPC\tn_GRPC\n")
        for i in sel:
            c, p = ch[i], int(po[i]); ma, mi = alle[(c, p)]
            o.write(f"{c}\t{p}\t{ma}\t{mi}\t{'Z' if c in zset else 'A'}\t{f0[i]:.4f}\t{f1[i]:.4f}\t{n0}\t{n1}\n")
            b.write(f"{c}\t{p - 1}\t{p}\n")
    nz = sum(ch[i] in zset for i in sel)
    print(f"[sites] {len(sel)} sites (dAF >= {a.daf}, 1 per {a.thin} bp): {len(sel) - nz} autosomal, {nz} Z; "
          f"reference birds {n0} LEPC, {n1} GRPC")


def pile_counts(bases, ref):
    out = {}
    bases = re.sub(r"\^.", "", bases).replace("$", "")
    i = 0
    while i < len(bases):
        c = bases[i]
        if c in "+-":
            m = re.match(r"[+-](\d+)", bases[i:]); i += len(m.group(0)) + int(m.group(1)); continue
        b = ref if c in ".," else c.upper()
        out[b] = out.get(b, 0) + 1; i += 1
    return out


def cmd_count(a):
    al = {}
    for r in read_tsv(a.sites):
        al[(r["chrom"], int(r["pos"]))] = (r["major"], r["minor"])
    n = 0
    with open(a.out, "w") as o:
        for line in sys.stdin:
            t = line.rstrip("\n").split("\t")
            k = (t[0], int(t[1]))
            if k not in al or len(t) < 5:
                continue
            c = pile_counts(t[4], t[2].upper()); ma, mi = al[k]
            if c.get(ma, 0) + c.get(mi, 0):
                o.write(f"{t[0]}\t{t[1]}\t{c.get(ma, 0)}\t{c.get(mi, 0)}\n"); n += 1
    print(f"[count] {n} sites with reads -> {a.out}")


def cmd_build(a):
    excl = {r["ID"] for r in read_tsv(a.exclude)} if a.exclude else set()
    sex = {r["ID"]: r["sex"] for r in read_tsv(a.sex)}
    birds = [r["ID"] for r in read_tsv(a.groups) if r["grp"] in a.grp.split(",") and r["ID"] not in excl]
    if a.cls == "Z":
        birds = [b for b in birds if sex.get(b) in ("M", "F")]
    birds = [b for b in birds if os.path.exists(os.path.join(a.counts, b + ".counts"))]
    sites = [r for r in read_tsv(a.sites) if r["class"] == a.cls]
    idx = {(r["chrom"], int(r["pos"])): i for i, r in enumerate(sites)}
    cnt = np.zeros((len(sites), len(birds), 2), np.int32)
    for j, b in enumerate(birds):
        for line in open(os.path.join(a.counts, b + ".counts")):
            c, p, x, y = line.split()
            i = idx.get((c, int(p)))
            if i is not None:
                cnt[i, j] = (int(x), int(y))
    # Morgans per bp: per-scaffold rate from --rate-table (scaffold<TAB>cM/Mb, e.g. the chicken rate of the homologous
    # chromosome) where given, else the uniform --rate-cm-mb. Z: x 2/3, because only the two Z copies in males recombine
    # (sex-averaged rate per Z copy), if the table gives the male (= only) Z map rate.
    rt = {}
    if a.rate_table:
        for line in open(a.rate_table):
            if line.strip() and not line.startswith("#"):
                c, v = line.split()[:2]
                rt[c] = float(v)
    zf = 2 / 3 if a.cls == "Z" else 1
    rate_of = lambda c: rt.get(c, a.rate_cm_mb) / 100 / 1e6 * zf
    keep = cnt.sum(axis=(1, 2)) > 0
    with open(a.out + ".input", "w") as o:
        prev = (None, 0)
        for i, r in enumerate(sites):
            if not keep[i]:
                continue
            c, p = r["chrom"], int(r["pos"])
            d = 0.5 if c != prev[0] else (p - prev[1]) * rate_of(c)
            prev = (c, p)
            n0, n1 = 2 * int(r["n_LEPC"]), 2 * int(r["n_GRPC"])
            a0 = int(round((1 - float(r["f_LEPC"])) * n0)); a1 = int(round((1 - float(r["f_GRPC"])) * n1))
            cols = [c, str(p), str(a0), str(n0 - a0), str(a1), str(n1 - a1), f"{d:.3e}"]
            cols += [f"{x}\t{y}" for x, y in cnt[i]]
            o.write("\t".join(cols) + "\n")
    with open(a.out + ".samples", "w") as o:
        for b in birds:
            o.write(f"{b}\t{1 if (a.cls == 'Z' and sex.get(b) == 'F') else 2}\n")
    print(f"[build] {a.grp} {a.cls}: {len(birds)} birds, {int(keep.sum())} sites -> {a.out}.input / .samples")


def cmd_summarize(a):
    import pandas as pd, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for run in a.runs:
        name = os.path.basename(run.rstrip("/"))
        for fn in sorted(glob.glob(os.path.join(run, "*.posterior"))):
            bird = os.path.basename(fn)[:-len(".posterior")]
            d = pd.read_csv(fn, sep="\t")
            states = d.columns[2:]
            g = np.array([int(s.split(",")[1]) for s in states])           # GRPC copies in each state
            ploidy = np.array([sum(int(v) for v in s.split(",")) for s in states])[0]
            exp = (d[states].values * g).sum(1) / ploidy
            rows.append(pd.DataFrame({"run": name, "bird": bird, "chrom": d.iloc[:, 0], "pos": d.iloc[:, 1],
                                      "grpc": exp}))
    if not rows:
        sys.exit("no posterior files found")
    x = pd.concat(rows)
    x["win"] = x.pos // a.window
    bird = x.groupby(["run", "bird"]).grpc.mean().rename("mean_grpc").reset_index()
    win = x.groupby(["run", "chrom", "win"]).agg(n_sites=("pos", "nunique"), mean_grpc=("grpc", "mean"),
                                                  max_bird=("grpc", "max")).reset_index()
    win["start"] = win.win * a.window
    bird.to_csv(a.out + "_per_bird.tsv", sep="\t", index=False)
    win.drop(columns="win").to_csv(a.out + "_windows.tsv", sep="\t", index=False)
    runs = list(win.run.unique())
    fig, axs = plt.subplots(len(runs), 1, figsize=(10, 2.2 * len(runs)), squeeze=False)
    for ax, r in zip(axs[:, 0], runs):
        w = win[win.run == r].copy()
        order = w.groupby("chrom").start.max().sort_values(ascending=False).index
        off, xs = {}, 0
        for c in order:
            off[c] = xs; xs += w[w.chrom == c].start.max() + a.window * 3
        w["x"] = (w.start + w.chrom.map(off)) / 1e6
        isz = w.chrom.isin(a.z.split(","))
        ax.scatter(w.x[~isz], w.mean_grpc[~isz], s=3, color="#A52A2A", lw=0)
        ax.scatter(w.x[isz], w.mean_grpc[isz], s=3, color="#000000", lw=0)
        ax.set_ylabel("Mean GRPC\nancestry"); ax.set_title(r, loc="left", fontsize=9, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    axs[-1, 0].set_xlabel("Scaffolds, longest first (Mb; Z in black)")
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=300)
    s = win.assign(Z=win.chrom.isin(a.z.split(","))).groupby(["run", "Z"]).mean_grpc.describe()
    print(s.round(4).to_string())


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("sites")
    for k in ("red", "mafs", "groups", "z", "out"):
        s.add_argument("--" + k, required=True)
    s.add_argument("--exclude"); s.add_argument("--daf", type=float, default=0.8); s.add_argument("--thin", type=int, default=5000)
    c = sp.add_parser("count"); c.add_argument("--sites", required=True); c.add_argument("--out", required=True)
    b = sp.add_parser("build")
    for k in ("sites", "counts", "groups", "sex", "grp", "cls", "out"):
        b.add_argument("--" + k, required=True)
    b.add_argument("--exclude"); b.add_argument("--rate-cm-mb", type=float, default=3.11)
    b.add_argument("--rate-table", help="scaffold<TAB>cM/Mb (sex-averaged; Z: male map)")
    m = sp.add_parser("summarize")
    m.add_argument("--runs", nargs="+", required=True); m.add_argument("--z", required=True)
    m.add_argument("--out", required=True); m.add_argument("--window", type=int, default=1_000_000)
    a = ap.parse_args()
    {"sites": cmd_sites, "count": cmd_count, "build": cmd_build, "summarize": cmd_summarize}[a.cmd](a)


if __name__ == "__main__":
    main()
