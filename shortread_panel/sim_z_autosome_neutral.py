#!/usr/bin/env python3
"""Neutral expectation for foreign (GRPC) ancestry on the Z vs the autosomes in a LEPC population
that receives GRPC migrants (shortread_panel/; runs anywhere, no cluster data needed).

Individual-based Wright-Fisher model with ZW sex determination and local ancestry tracked along
chromosomes:
  * N birds (half female) per generation; each generation a fraction m of the parents are replaced by
    immigrants: m_G pure GRPC and m_L pure LEPC (m_L keeps the deme LEPC-majority, as in the contact
    zone and range edge). A fraction f_fem of GRPC immigrants are female (0 = males only, 1 = females only).
  * Autosomes: N_CHR chromosomes of 1 Morgan, LOCI loci each, crossovers in both sexes.
    Z: one 1-Morgan chromosome; crossovers only in male meiosis (ZW females pass on their single Z
    intact to sons, the W to daughters).
  * Males: two Z copies; females: one.
Summaries per sampled bird (same definitions as hybrid_index.py): S (GRPC fraction), H (fraction of
loci with one ancestry of each) on the autosomes, S on the Z.
Reported per scenario, matching the observed statistics:
  Z/A      ratio of mean GRPC ancestry per Z copy to mean autosomal GRPC ancestry (all sampled birds)
  adm_Z/A  the same for admixed non-F1 birds (autosomal S 0.1-0.4, H < 0.8)
  adm_noZ  fraction of those admixed birds with no foreign Z (Z S < 0.02)
Scenarios also include one with F1 hybrid females sterile (Haldane's rule) and a pulse (no ongoing
migration, generations since admixture).

usage: python3 sim_z_autosome_neutral.py [--reps 20] [--out sim_z_autosome.tsv]"""
import argparse, itertools
import numpy as np

N_CHR, LOCI = 20, 50


def gamete(h0, h1, r, rng):
    """h0, h1: n x L ancestry (0 LEPC, 1 GRPC) of the two copies; recombine along L loci with per-interval
    crossover probability r. Returns n x L."""
    n, L = h0.shape
    start = rng.random(n) < 0.5
    sw = rng.random((n, L - 1)) < r
    phase = np.concatenate([start[:, None], sw], axis=1).cumsum(1) % 2 == 1
    return np.where(phase, h1, h0)


def simulate(N, gens, mG, mL, ffem, sterile_f1, pulse, rng, sample=150):
    L, rA, rZ = LOCI, 1.0 / LOCI, 1.0 / LOCI
    # state: autosomes A[n, 2, N_CHR, L]; Z[n, 2, L] (females: copy 1 unused = -1); sex 0 male 1 female
    sex = (np.arange(N) % 2).astype(np.int8)
    A = np.zeros((N, 2, N_CHR, L), np.int8)
    Z = np.zeros((N, 2, L), np.int8); Z[sex == 1, 1] = -1
    if pulse:                          # a single pulse: a fraction mG of the founders are GRPC
        k = rng.random(N) < mG * (2 if ffem in (0.0, 1.0) else 1)   # same overall GRPC share when one sex only
        if ffem == 1.0:
            k &= sex == 1
        elif ffem == 0.0:
            k &= sex == 0
        A[k] = 1; Z[k & (sex == 0)] = 1; Z[k & (sex == 1), 0] = 1
    for g in range(gens):
        if not pulse and g < gens:     # immigrants replace residents before mating
            u = rng.random(N)
            gi = u < mG; li = (u >= mG) & (u < mG + mL)
            # GRPC immigrant sex follows ffem
            sex[gi] = (rng.random(gi.sum()) < ffem).astype(np.int8)
            sex[li] = (rng.random(li.sum()) < 0.5).astype(np.int8)
            for msk, a in ((gi, 1), (li, 0)):
                A[msk] = a; Z[msk] = a
            Z[sex == 1, 1] = -1
        males, females = np.where(sex == 0)[0], np.where(sex == 1)[0]
        if sterile_f1:                 # Haldane's rule: F1-like females (autosomal H > 0.8) do not breed
            Hf = (A[females, 0] != A[females, 1]).mean(axis=(1, 2))
            females = females[Hf < 0.8]
        if len(males) == 0 or len(females) == 0:
            break
        fa = rng.choice(males, N); mo = rng.choice(females, N)
        nsex = (rng.random(N) < 0.5).astype(np.int8)
        nA = np.empty_like(A)
        for c in range(N_CHR):
            nA[:, 0, c] = gamete(A[fa, 0, c], A[fa, 1, c], rA, rng)
            nA[:, 1, c] = gamete(A[mo, 0, c], A[mo, 1, c], rA, rng)
        nZ = np.full_like(Z, -1)
        nZ[:, 0] = gamete(Z[fa, 0], Z[fa, 1], rZ, rng)        # every chick gets a recombined paternal Z
        sons = nsex == 0
        nZ[sons, 1] = Z[mo[sons], 0]                           # sons also get the mother's single Z
        A, Z, sex = nA, nZ, nsex
    idx = rng.choice(N, min(sample, N), replace=False)
    a = A[idx]; z = Z[idx]; sx = sex[idx]
    S_A = a.mean(axis=(1, 2, 3)); H_A = (a[:, 0] != a[:, 1]).mean(axis=(1, 2))
    S_Z = np.where(sx == 0, z.astype(float).mean(axis=(1, 2)), z[:, 0].mean(1))
    nz = np.where(sx == 0, 2, 1)
    return S_A, H_A, S_Z, nz


def summarise(S_A, H_A, S_Z, nz):
    za = (S_Z * nz).sum() / nz.sum() / max(S_A.mean(), 1e-9)
    adm = (S_A >= 0.1) & (S_A <= 0.4) & (H_A < 0.8)
    if adm.sum():
        az = (S_Z[adm] * nz[adm]).sum() / nz[adm].sum() / S_A[adm].mean()
        nf = (S_Z[adm] < 0.02).mean()
    else:
        az = nf = np.nan
    return S_A.mean(), za, adm.sum(), az, nf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--N", type=int, default=600)
    ap.add_argument("--out", default="sim_z_autosome.tsv")
    a = ap.parse_args()
    rng = np.random.default_rng(42)
    rows = []
    scen = []
    for mG, ffem, ster in itertools.product((0.005, 0.02), (0.0, 0.5, 1.0), (False, True)):
        scen.append(("migration", mG, 0.05, ffem, ster, 60))
    for gens, ffem in itertools.product((2, 5, 10), (0.5, 1.0)):
        scen.append(("pulse", 0.25, 0.0, ffem, False, gens))
    for kind, mG, mL, ffem, ster, gens in scen:
        res = [summarise(*simulate(a.N, gens, mG, mL, ffem, ster, kind == "pulse", rng)) for _ in range(a.reps)]
        res = np.array(res, float)
        m = np.nanmean(res, 0); lo = np.nanpercentile(res, 2.5, 0); hi = np.nanpercentile(res, 97.5, 0)
        rows.append([kind, mG, mL, ffem, int(ster), gens, f"{m[0]:.3f}", f"{m[1]:.2f} ({lo[1]:.2f}-{hi[1]:.2f})",
                     f"{m[2]:.1f}", f"{m[3]:.2f} ({lo[3]:.2f}-{hi[3]:.2f})", f"{m[4]:.2f} ({lo[4]:.2f}-{hi[4]:.2f})"])
        print("\t".join(map(str, rows[-1])), flush=True)
    hdr = ["scenario", "m_GRPC", "m_LEPC", "GRPC_migrants_female", "F1_females_sterile", "generations",
           "mean_auto_S", "Z/A", "n_admixed_per_150", "admixed_Z/A", "admixed_noZ"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
