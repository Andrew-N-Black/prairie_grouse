#!/usr/bin/env python3
"""Hybrid index (S) and interspecific heterozygosity (H) from genotype likelihoods
(shortread_panel/; driven by hybrid_index.sh).

  reduce    one ANGSD chunk (beagle GLs, all panel birds) -> candidate ancestry-informative sites:
            per-site allele frequencies in the allopatric reference groups (EM on the GLs) and the
            GLs of every bird at sites where a reference pair differs by >= --pre-daf. Saved as .npz.
  estimate  all chunks -> final markers (dAF >= --daf between the two reference groups, thinned to
            one per --thin bp), then per bird and per chromosome class (autosomes, Z):
              diploid (autosomes; Z of males): maximum-likelihood S and H (HIest model of
                Fitzpatrick 2012, summed over genotype uncertainty). Ancestry classes at a locus:
                both copies species 2 (p = S - H/2), one of each (H), both species 1 (1 - S - H/2);
                each copy's allele is drawn from that species' reference frequency.
              haploid (Z of females, one Z): S only; per site, L = S*P(reads|sp2) + (1-S)*P(reads|sp1)
                using the GLs of the two homozygous genotypes as the allele likelihoods.
            Reference birds are scored leave-one-out: marker choice and frequencies are recomputed
            without them. 95% intervals from a bootstrap of 1-Mb blocks.
  plot      triangle plots (S vs H) for autosomes and the male Z, and Z vs autosomal S.

S is the fraction of ancestry from species 2 of --pair (e.g. LEPC-GRPC: S = GRPC ancestry)."""
import argparse, glob, gzip, os, sys
import multiprocessing as mp
import numpy as np


def em_freq(gl, n_iter=30):
    """gl: sites x inds x 3 (genotype likelihoods of 0/1/2 minor alleles). Returns ML minor-allele
    frequency per site and the posterior expected dosage per bird."""
    f = np.full(gl.shape[0], 0.2)
    for _ in range(n_iter):
        pri = np.stack([(1 - f) ** 2, 2 * f * (1 - f), f ** 2], axis=1)[:, None, :]
        post = gl * pri
        post /= post.sum(axis=2, keepdims=True) + 1e-300
        dos = post[:, :, 1] + 2 * post[:, :, 2]
        f = np.clip(dos.mean(axis=1) / 2, 1e-4, 1 - 1e-4)
    return f, dos


# --------------------------------------------------------------------------- reduce
def cmd_reduce(a):
    ids = [l.strip() for l in open(a.ids) if l.strip()]
    grp = {}
    for i, l in enumerate(open(a.groups)):
        t = l.rstrip("\n").split("\t")
        if i == 0:
            ci, cg = t.index("ID"), t.index("grp")
            continue
        grp[t[ci]] = t[cg]
    refs = a.refs.split(",")
    idx = {g: np.array([i for i, b in enumerate(ids) if grp.get(b) == g]) for g in refs}
    chrom, pos, rows = [], [], []
    op = gzip.open if a.beagle.endswith(".gz") else open
    with op(a.beagle, "rt") as fh:
        hdr = fh.readline().split()
        n = (len(hdr) - 3) // 3
        assert n == len(ids), f"beagle has {n} birds, id list {len(ids)}"
        for line in fh:
            t = line.split()
            c, p = t[0].rsplit("_", 1)
            chrom.append(c); pos.append(int(p))
            rows.append(np.asarray(t[3:], dtype=np.float32))
    if not rows:
        np.savez_compressed(a.out, chrom=np.array([], dtype="U40"), pos=np.array([], dtype=np.int64),
                            gl=np.zeros((0, len(ids), 3), np.float16), freqs=np.zeros((0, len(refs))),
                            refs=np.array(refs), ids=np.array(ids))
        print("[reduce] 0 sites"); return
    gl = np.stack(rows).reshape(len(rows), len(ids), 3)
    gl /= gl.sum(axis=2, keepdims=True)
    fr = np.stack([em_freq(gl[:, idx[g], :])[0] for g in refs], axis=1)
    pairs = [(i, j) for i in range(len(refs)) for j in range(i + 1, len(refs))]
    daf = np.max(np.stack([np.abs(fr[:, i] - fr[:, j]) for i, j in pairs], axis=1), axis=1)
    k = daf >= a.pre_daf
    np.savez_compressed(a.out, chrom=np.array(chrom)[k], pos=np.array(pos)[k], gl=gl[k].astype(np.float16),
                        freqs=fr[k], refs=np.array(refs), ids=np.array(ids))
    print(f"[reduce] {len(rows)} SNPs, {k.sum()} with a reference-pair dAF >= {a.pre_daf}")


# --------------------------------------------------------------------------- likelihoods
def site_terms(gl, p1, p2):
    """Per site, likelihood of the reads under each ancestry class (11, 12, 22) for one bird.
    gl: sites x 3; p1, p2: minor-allele frequency in species 1 and 2."""
    def geno(pa, pb):    # P(genotype 0/1/2) when one copy is drawn with pa and the other with pb
        return np.stack([(1 - pa) * (1 - pb), pa * (1 - pb) + pb * (1 - pa), pa * pb], axis=1)
    return np.stack([(gl * geno(p1, p1)).sum(1), (gl * geno(p1, p2)).sum(1), (gl * geno(p2, p2)).sum(1)], axis=1)


def em_mix(c, W, n_iter=500, tol=1e-6):
    """ML mixture weights. c: sites x K (likelihood of the reads under each ancestry class);
    W: reps x sites (site weights; 1 for the point estimate, bootstrap counts otherwise).
    Returns reps x K weights. For the diploid model K = 3 (classes 11, 12, 22), so
    H = w12 and S = w22 + w12/2; the HIest constraints hold automatically."""
    K = c.shape[1]
    pi = np.full((W.shape[0], K), 1.0 / K)
    tot = W.sum(1, keepdims=True)
    for _ in range(n_iter):
        m = c[None, :, :] * pi[:, None, :]                      # reps x sites x K
        m /= m.sum(2, keepdims=True) + 1e-300
        new = (W[:, :, None] * m).sum(1) / tot
        if np.max(np.abs(new - pi)) < tol:
            pi = new; break
        pi = new
    return pi


def fit_diploid(c, W):
    pi = em_mix(c, W)
    return pi[:, 2] + pi[:, 1] / 2, pi[:, 1]                     # S, H


def fit_haploid(gl, p1, p2, W):
    """Female Z: allele likelihoods from the GLs of the two homozygous genotypes; S = weight of species 2."""
    la, lb = gl[:, 0], gl[:, 2]                                   # P(reads | major), P(reads | minor)
    c = np.stack([(1 - p1) * la + p1 * lb, (1 - p2) * la + p2 * lb], axis=1)
    return em_mix(c, W)[:, 1]


# --------------------------------------------------------------------------- estimate
def thin(chrom, pos, score, bp):
    """keep the highest-scoring site in each bp-long bin per chromosome"""
    key = np.array([f"{c}:{p // bp}" for c, p in zip(chrom, pos)])
    order = np.lexsort((-score, key))
    _, first = np.unique(key[order], return_index=True)
    return np.sort(order[first])


_CTX = {}


def _one_bird(b_i):
    """S and H (with block-bootstrap intervals) for one bird; reads the shared arrays from _CTX."""
    X = _CTX; a = X["a"]; bird = X["ids"][b_i]
    p1, p2 = X["f1"].copy(), X["f2"].copy()
    use = np.ones(len(p1), bool)
    if X["m1"][b_i] or X["m2"][b_i]:      # leave-one-out: frequencies and marker choice without the bird
        if X["m1"][b_i]:
            p1 = np.clip((X["D1"] - X["d1"][:, X["r1"][b_i]]) / (2 * (X["n1"] - 1)), 1e-4, 1 - 1e-4)
        else:
            p2 = np.clip((X["D2"] - X["d2"][:, X["r2"][b_i]]) / (2 * (X["n2"] - 1)), 1e-4, 1 - 1e-4)
        use = np.abs(p1 - p2) >= a.daf
    g = X["G"][:, b_i, :]
    info = use & (np.ptp(g, axis=1) > 1e-6)          # flat GLs (no reads) carry no information
    sx = X["sexd"].get(bird, {}).get("sex", "?")
    if X["cls"] == "Z" and sx not in ("M", "F"):
        return None
    hap = X["cls"] == "Z" and sx == "F"
    gi, q1, q2, bi = g[info], p1[info], p2[info], X["binv"][info]
    if len(gi) < a.min_sites:
        return None
    rng = np.random.default_rng(1000 + b_i); nb = X["nb"]
    W = np.vstack([np.ones(len(gi))] +
                  [np.bincount(rng.integers(0, nb, nb), minlength=nb)[bi].astype(float) for _ in range(a.boot)])
    if hap:
        Sv = fit_haploid(gi, q1, q2, W); Hv = np.full(len(Sv), np.nan)
    else:
        Sv, Hv = fit_diploid(site_terms(gi, q1, q2), W)
    S, H = Sv[0], Hv[0]
    slo, shi = np.percentile(Sv[1:], [2.5, 97.5])
    hlo, hhi = (np.nan, np.nan) if hap else np.percentile(Hv[1:], [2.5, 97.5])
    md = X["meta"].get(bird, {}); sm = X["sexd"].get(bird, {})
    return [bird, md.get("common", ""), md.get("grp", ""), md.get("GROUP", ""), md.get("STATE_NAME", ""),
            sm.get("county", ""), sx, X["pair"], X["cls"], "haploid" if hap else "diploid", str(len(gi))] + \
           [f"{v:.4f}" if np.isfinite(v) else "NA" for v in (S, slo, shi, H, hlo, hhi)]


def cmd_estimate(a):
    files = sorted(glob.glob(os.path.join(a.indir, "*.npz")))
    assert files, "no .npz chunks in " + a.indir
    ch, po, gl, fr = [], [], [], []
    for f in files:
        z = np.load(f)
        if len(z["pos"]) == 0:
            continue
        ch.append(z["chrom"]); po.append(z["pos"]); gl.append(z["gl"].astype(np.float32)); fr.append(z["freqs"])
        refs, ids = list(z["refs"]), list(z["ids"])
    ch = np.concatenate(ch); po = np.concatenate(po); gl = np.concatenate(gl); fr = np.concatenate(fr)
    zset = set(a.z.split(","))
    isz = np.array([c in zset for c in ch])
    meta = {}
    for i, l in enumerate(open(a.groups)):
        t = l.rstrip("\n").split("\t")
        if i == 0:
            hdr = t; continue
        meta[t[hdr.index("ID")]] = dict(zip(hdr, t))
    sexd = {}
    for i, l in enumerate(open(a.sex)):
        t = l.rstrip("\n").split("\t")
        if i == 0:
            hs = t; continue
        sexd[t[0]] = dict(zip(hs, t))
    out = open(a.out, "w")
    cols = ["ID", "species", "grp", "context", "state", "county", "sex", "pair", "class", "model", "n_sites",
            "S", "S_lo", "S_hi", "H", "H_lo", "H_hi"]
    out.write("\t".join(cols) + "\n")
    summ = []
    for pair in a.pairs.split(","):
        g1, g2 = [f"{s}_allo" for s in pair.split("-")]
        i1, i2 = refs.index(g1), refs.index(g2)
        m1 = np.array([meta.get(b, {}).get("grp") == g1 for b in ids])
        m2 = np.array([meta.get(b, {}).get("grp") == g2 for b in ids])
        # full-panel reference frequencies and per-bird dosages (for leave-one-out)
        for cls, sel in (("autosomes", ~isz), ("Z", isz)):
            S_idx = np.where(sel)[0]
            G = gl[S_idx]
            f1, d1 = em_freq(G[:, m1, :]); f2, d2 = em_freq(G[:, m2, :])
            daf = np.abs(f1 - f2)
            keep0 = daf >= a.daf
            kk = thin(ch[S_idx][keep0], po[S_idx][keep0], daf[keep0], a.thin)
            sites = np.where(keep0)[0][kk]
            blocks = np.array([f"{c}:{p // a.block}" for c, p in zip(ch[S_idx][sites], po[S_idx][sites])])
            ub, binv = np.unique(blocks, return_inverse=True)
            summ.append(f"{pair}\t{cls}\t{len(sites)} markers (dAF >= {a.daf}, 1 per {a.thin} bp) in {len(ub)} blocks; "
                        f"{int(m1.sum())} {g1}, {int(m2.sum())} {g2} reference birds")
            n1, n2 = m1.sum(), m2.sum()
            D1, D2 = d1[sites].sum(1), d2[sites].sum(1)
            r1 = np.cumsum(m1) - 1; r2 = np.cumsum(m2) - 1            # bird -> column in d1 / d2
            _CTX.update(dict(ids=ids, m1=m1, m2=m2, f1=f1[sites], f2=f2[sites], D1=D1, D2=D2,
                             d1=d1[sites], d2=d2[sites], r1=r1, r2=r2, n1=n1, n2=n2, G=G[sites], binv=binv,
                             nb=len(ub), cls=cls, pair=pair, meta=meta, sexd=sexd, a=a))
            with mp.get_context("fork").Pool(a.threads) as pool:
                for row in pool.imap(_one_bird, range(len(ids)), chunksize=4):
                    if row:
                        out.write("\t".join(row) + "\n")
            sys.stderr.write(summ[-1] + "\n")
    out.close()
    open(a.out.replace(".tsv", "") + ".markers.txt", "w").write("\n".join(summ) + "\n")


# --------------------------------------------------------------------------- plot
def cmd_plot(a):
    import pandas as pd, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"GRPC": "#DAA520", "LEPC": "#A52A2A", "STGR": "#000000", "STGR / GRPC": "#7F7F7F"}
    d = pd.read_csv(a.tsv, sep="\t")
    plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False})
    for pair, dp in d.groupby("pair"):
        s1, s2 = pair.split("-")
        # only birds labelled as one of the two species (or a hybrid of them); a third species has no
        # meaning on this axis
        dp = dp[dp.species.astype(str).apply(lambda x: any(s in x for s in (s1, s2)))]
        fig, axs = plt.subplots(1, 3, figsize=(10, 3.4))
        for ax, (cls, title) in zip(axs[:2], [("autosomes", "Autosomes"), ("Z", "Z (males)")]):
            t = dp[(dp["class"] == cls) & (dp.model == "diploid")]
            ax.plot([0, 0.5, 1, 0], [0, 1, 0, 0], color="#999999", lw=0.8)
            for sp, ts in t.groupby("species"):
                for ctx, mk in (("Allopatric", "o"), ("Edge", "s"), ("Sympatric", "^")):
                    u = ts[ts.context == ctx]
                    if len(u) == 0:
                        continue
                    ax.scatter(u.S, u.H, s=14, marker=mk, facecolors=COL.get(sp, "#7F7F7F") if ctx == "Allopatric" else "none",
                               edgecolors=COL.get(sp, "#7F7F7F"), linewidths=0.7,
                               label=f"{sp} {'range edge' if ctx == 'Edge' else ctx.lower()}")
            ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
            ax.set_xlabel(f"Hybrid index (S, {s2} ancestry)"); ax.set_ylabel("Interspecific heterozygosity (H)")
            ax.set_title(title, loc="left", fontsize=9, fontweight="bold")
        axs[0].legend(frameon=False, fontsize=6.5, loc="upper right")
        w = dp.pivot_table(index=["ID", "species", "sex"], columns="class", values="S").reset_index().dropna()
        ax = axs[2]
        ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, ls=(0, (3, 2)))
        for sp, ts in w.groupby("species"):
            for sx, mk in (("M", "o"), ("F", "s")):
                u = ts[ts.sex == sx]
                ax.scatter(u.autosomes, u.Z, s=14, marker=mk, color=COL.get(sp, "#7F7F7F"), alpha=0.85, lw=0)
        ax.set_xlabel("Autosomal S"); ax.set_ylabel("Z S (females: one Z)")
        ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
        ax.set_title("Z vs autosomes", loc="left", fontsize=9, fontweight="bold")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{a.prefix}_{pair}.{ext}", dpi=300)
        print(f"wrote {a.prefix}_{pair}.png/pdf")


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reduce")
    r.add_argument("--beagle", required=True); r.add_argument("--ids", required=True)
    r.add_argument("--groups", required=True); r.add_argument("--out", required=True)
    r.add_argument("--refs", default="LEPC_allo,GRPC_allo,STGR_allo")
    r.add_argument("--pre-daf", type=float, default=0.5)
    e = sp.add_parser("estimate")
    e.add_argument("--indir", required=True); e.add_argument("--groups", required=True)
    e.add_argument("--sex", required=True); e.add_argument("--out", required=True)
    e.add_argument("--z", required=True, help="comma-separated Z scaffolds")
    e.add_argument("--pairs", default="LEPC-GRPC,GRPC-STGR")
    e.add_argument("--daf", type=float, default=0.8)
    e.add_argument("--thin", type=int, default=5000)
    e.add_argument("--block", type=int, default=1_000_000)
    e.add_argument("--boot", type=int, default=100)
    e.add_argument("--min-sites", type=int, default=50)
    e.add_argument("--threads", type=int, default=8)
    p = sp.add_parser("plot")
    p.add_argument("--tsv", required=True); p.add_argument("--prefix", required=True)
    a = ap.parse_args()
    {"reduce": cmd_reduce, "estimate": cmd_estimate, "plot": cmd_plot}[a.cmd](a)


if __name__ == "__main__":
    main()
