#!/usr/bin/env python3
"""Step 13 summary: classify pangenome SVs by species and compare the Z with the autosomes.

Reads, per chromosome, <sv-dir>/<CHROM>/<CHROM>.sv.vcf.gz (vg deconstruct + vcfbub),
ref.txt (reference path and its length) and haps_per_species.tsv. Python 3 stdlib only.

The VCFs are vcfwave-decomposed. Each ALT allele whose length differs from REF by >= MIN_SV
is an INS or DEL; a record vcfwave flags as an inversion (INFO INV) with alleles >= MIN_SV is
INV; any other allele where REF and ALT are both >= MIN_SV is COMPLEX. Shorter alleles are
ignored, as are exact duplicate records and alleles carried by no genotyped haplotype (vg
omits the reference bird, so variants private to it have no carrier).

For each SV allele, each species' frequency is computed over its called haplotypes. A
species must have at least MIN_CALLED_FRAC x (haplotypes carrying the chromosome) called.
  diagnostic for X : X fixed for one state (allele present in all or none of X's
                     haplotypes) and both other species fixed for the opposite state
  fixed difference : for a species pair, one fixed present and the other fixed absent
  private to X     : present only in X, not fixed
Each classification is made twice: with all called haplotypes, and after randomly
subsampling EQUAL_N called haplotypes per species ("_eqN"), so that chromosomes or
species with fewer haplotypes (GRPC Z: 7) are not favoured by small samples.
chr_W is counted but not classified (no LEPC W).
"""
import argparse, gzip, math, os, random, re
from collections import Counter, defaultdict

SPECIES = ["GRPC", "LEPC", "STGR"]
PAIRS = [("GRPC", "LEPC"), ("GRPC", "STGR"), ("LEPC", "STGR")]
BINS = [(50, 100), (100, 300), (300, 1000), (1000, 5000), (5000, 10000), (10000, 10**9)]


def bin_label(n):
    for lo, hi in BINS:
        if lo <= n < hi:
            return f"{lo}-{hi - 1}" if hi < 10**9 else f">={lo}"
    return "NA"


def fisher_two_sided(a, b, c, d):
    """Fisher's exact test for [[a, b], [c, d]]; returns (odds ratio, two-sided P)."""
    def lp(x):  # log hypergeometric probability of x in cell a
        r1, c1, n = a + b, a + c, a + b + c + d
        return (math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                - math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1)
                - math.lgamma(c1 - x + 1) - math.lgamma(n - r1 - c1 + x + 1))
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    p0 = lp(a)
    p = sum(math.exp(lp(x)) for x in range(lo, hi + 1) if lp(x) <= p0 + 1e-7)
    orr = (a * d) / (b * c) if b * c else float("inf")
    return orr, min(1.0, p)


def state(freq):
    return 1 if freq == 1 else 0 if freq == 0 else None


def classify(calls):
    """calls: species -> list of 0/1 (carries allele) over called haplotypes. Returns labels."""
    fr = {s: sum(v) / len(v) for s, v in calls.items()}
    st = {s: state(f) for s, f in fr.items()}
    diag = ""
    for x in SPECIES:
        others = [s for s in SPECIES if s != x]
        if st[x] is not None and all(st[o] == 1 - st[x] for o in others):
            diag = x
    fixed = [f"{p}_{q}" for p, q in PAIRS if st[p] is not None and st[q] is not None and st[p] != st[q]]
    carriers = [s for s in SPECIES if fr[s] > 0]
    private = carriers[0] if len(carriers) == 1 and st[carriers[0]] != 1 else ""
    return fr, diag, fixed, private


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sv-dir", required=True)
    ap.add_argument("--chroms", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-sv", type=int, default=50)
    ap.add_argument("--min-called-frac", type=float, default=0.7)
    ap.add_argument("--equal-n", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    eq = f"eq{a.equal_n}"
    os.makedirs(a.out, exist_ok=True)

    chroms = [l.strip() for l in open(a.chroms) if l.strip()]
    per_chrom, sizes = {}, Counter()
    tab = gzip.open(os.path.join(a.out, "sv_table.tsv.gz"), "wt")
    tab.write("\t".join(["chrom", "pos", "ref_len", "alt_len", "type", "size"]
                        + [f"n_{s}" for s in SPECIES] + [f"af_{s}" for s in SPECIES]
                        + ["diagnostic", "fixed_pairs", "private", f"diagnostic_{eq}", f"fixed_pairs_{eq}"]) + "\n")

    for chrom in chroms:
        d = os.path.join(a.sv_dir, chrom)
        vcf = os.path.join(d, f"{chrom}.sv.vcf.gz")
        if not os.path.exists(vcf):
            print(f"WARNING: {vcf} missing, skipped")
            continue
        ref_path, ref_len = open(os.path.join(d, "ref.txt")).read().split()[:2]
        nhap = {r.split("\t")[0]: int(r.split("\t")[1]) for r in open(os.path.join(d, "haps_per_species.tsv")).read().splitlines()[1:]}
        need = {s: max(1, math.ceil(a.min_called_frac * nhap.get(s, 0))) for s in SPECIES}
        classifiable = all(nhap.get(s, 0) > 0 for s in SPECIES)
        c = Counter()
        cols = None
        seen = set()                                  # vcfwave can emit the same variant twice
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("##"):
                    continue
                f = line.rstrip("\n").split("\t")
                if line.startswith("#"):
                    cols = [(i, n.split(".")[0]) for i, n in enumerate(f[9:], start=9)]
                    continue
                ref, alts = f[3], f[4].split(",")
                is_inv = bool(re.search(r'(^|;)INV(=|;|$)', f[7]))
                gi = f[8].split(":").index("GT")
                gts = defaultdict(list)            # species -> allele index per called haplotype
                for i, sp in cols:
                    for al in re.split(r"[|/]", f[i].split(":")[gi]):
                        if al != ".":
                            gts[sp].append(int(al))
                for k, alt in enumerate(alts, start=1):
                    if alt.startswith("<") or alt == "*":
                        continue
                    rl, al_ = len(ref), len(alt)
                    diff = al_ - rl
                    if is_inv and max(rl, al_) >= a.min_sv:
                        typ, size = "INV", max(rl, al_)
                    elif diff >= a.min_sv:
                        typ, size = "INS", diff
                    elif diff <= -a.min_sv:
                        typ, size = "DEL", -diff
                    elif min(rl, al_) >= a.min_sv:
                        typ, size = "COMPLEX", max(rl, al_)
                    else:
                        continue
                    key = (f[1], ref, alt)
                    if key in seen:
                        c["duplicate"] += 1
                        continue
                    seen.add(key)
                    if not any(x == k for s in gts for x in gts[s]):
                        c["no_carrier"] += 1         # carried only by the reference bird (not genotyped) or no one
                        continue
                    c["SV"] += 1; c[typ] += 1
                    sizes[(typ, bin_label(size), "all")] += 1
                    ok = classifiable and all(len(gts[s]) >= need[s] for s in SPECIES)
                    fr, diag, fixed, priv = ({s: float("nan") for s in SPECIES}, "", [], "")
                    diag_e, fixed_e = "", []
                    if ok:
                        c["classified"] += 1
                        calls = {s: [int(x == k) for x in gts[s]] for s in SPECIES}
                        fr, diag, fixed, priv = classify(calls)
                        if all(len(calls[s]) >= a.equal_n for s in SPECIES):
                            c["classified_eq"] += 1
                            sub = {s: rng.sample(calls[s], a.equal_n) for s in SPECIES}
                            _, diag_e, fixed_e, _ = classify(sub)
                        for lab in ([f"diag_{diag}"] if diag else []) + [f"fixed_{p}" for p in fixed] \
                                + ([f"private_{priv}"] if priv else []) \
                                + ([f"diag_{diag_e}_{eq}"] if diag_e else []) + [f"fixed_{p}_{eq}" for p in fixed_e]:
                            c[lab] += 1
                        if diag:
                            c["diag_any"] += 1; sizes[(typ, bin_label(size), "diagnostic")] += 1
                        if diag_e:
                            c[f"diag_any_{eq}"] += 1; sizes[(typ, bin_label(size), f"diagnostic_{eq}")] += 1
                    tab.write("\t".join(map(str, [chrom, f[1], rl, al_, typ, size]
                              + [len(gts[s]) for s in SPECIES] + [f"{fr[s]:.3f}" for s in SPECIES]
                              + [diag, ",".join(fixed), priv, diag_e, ",".join(fixed_e)])) + "\n")
        c["ref_len"] = int(ref_len)
        per_chrom[chrom] = (ref_path, c)
        print(f"{chrom}: {c['SV']:,} SVs ({c['classified']:,} classifiable; skipped {c['duplicate']:,} duplicates, "
              f"{c['no_carrier']:,} with no genotyped carrier), "
              f"diagnostic {c['diag_any']:,} ({c[f'diag_any_{eq}']:,} at n={a.equal_n}), ref {ref_path}")
    tab.close()

    # ---- per chromosome ----
    keys = (["SV", "INS", "DEL", "INV", "COMPLEX", "classified", "classified_eq", "diag_any", f"diag_any_{eq}"]
            + [f"diag_{s}" for s in SPECIES] + [f"diag_{s}_{eq}" for s in SPECIES]
            + [f"fixed_{p}_{q}" for p, q in PAIRS] + [f"fixed_{p}_{q}_{eq}" for p, q in PAIRS]
            + [f"private_{s}" for s in SPECIES])
    with open(os.path.join(a.out, "sv_per_chromosome.tsv"), "w") as f:
        f.write("chrom\treference_path\tref_Mb\tSV_per_Mb\t" + "\t".join(keys) + "\n")
        for chrom, (ref_path, c) in per_chrom.items():
            mb = c["ref_len"] / 1e6
            f.write(f"{chrom}\t{ref_path}\t{mb:.2f}\t{c['SV'] / mb:.1f}\t" + "\t".join(str(c[k]) for k in keys) + "\n")

    # ---- size distribution ----
    with open(os.path.join(a.out, "sv_size_distribution.tsv"), "w") as f:
        f.write("type\tsize_bin\tall\tdiagnostic\tdiagnostic_" + eq + "\n")
        for typ in ["INS", "DEL", "INV", "COMPLEX"]:
            for lo, hi in BINS:
                b = bin_label(lo)
                f.write(f"{typ}\t{b}\t{sizes[(typ, b, 'all')]}\t{sizes[(typ, b, 'diagnostic')]}\t"
                        f"{sizes[(typ, b, 'diagnostic_' + eq)]}\n")

    # ---- Z versus autosomes ----
    auto = [ch for ch in per_chrom if ch not in ("chr_Z", "chr_W")]
    lines, text = [], []
    if "chr_Z" in per_chrom and auto:
        z = per_chrom["chr_Z"][1]
        A = Counter(); [A.update(per_chrom[ch][1]) for ch in auto]
        zmb, amb = z["ref_len"] / 1e6, A["ref_len"] / 1e6
        cats = ([("diag_any", "classified"), (f"diag_any_{eq}", "classified_eq")]
                + [(f"diag_{s}", "classified") for s in SPECIES] + [(f"diag_{s}_{eq}", "classified_eq") for s in SPECIES]
                + [(f"fixed_{p}_{q}", "classified") for p, q in PAIRS]
                + [(f"fixed_{p}_{q}_{eq}", "classified_eq") for p, q in PAIRS])
        hdr = ("category\tZ_n\tZ_classified\tZ_frac\tauto_n\tauto_classified\tauto_frac\t"
               "odds_ratio\tfisher_P\tZ_per_Mb\tauto_per_Mb\tdensity_ratio")
        for cat, den in cats:
            zn, zt, an, at = z[cat], z[den], A[cat], A[den]
            orr, p = fisher_two_sided(zn, zt - zn, an, at - an)
            zr, ar = zn / zmb, an / amb
            lines.append(f"{cat}\t{zn}\t{zt}\t{zn / zt if zt else float('nan'):.4f}\t{an}\t{at}\t"
                         f"{an / at if at else float('nan'):.4f}\t{orr:.2f}\t{p:.3g}\t{zr:.2f}\t{ar:.2f}\t"
                         f"{zr / ar if ar else float('nan'):.2f}")
        with open(os.path.join(a.out, "sv_Z_vs_autosomes.tsv"), "w") as f:
            f.write(hdr + "\n" + "\n".join(lines) + "\n")
        text += [f"Z: {z['SV']:,} SVs on {zmb:.1f} Mb ({z['SV'] / zmb:.1f}/Mb); autosomes: {A['SV']:,} on {amb:.1f} Mb "
                 f"({A['SV'] / amb:.1f}/Mb)"]
        text += ["Z vs autosomes (fraction of classifiable SVs; Fisher's exact test):"] + ["  " + l for l in lines]

    tot = Counter(); [tot.update(c) for _, c in per_chrom.values()]
    text = [f"SV alleles >= {a.min_sv} bp: {tot['SV']:,} (INS {tot['INS']:,}, DEL {tot['DEL']:,}, INV {tot['INV']:,}, COMPLEX {tot['COMPLEX']:,})",
            f"classifiable (all three species called): {tot['classified']:,}; at n = {a.equal_n}: {tot['classified_eq']:,}",
            "diagnostic: " + ", ".join(f"{s} {tot['diag_' + s]:,}" for s in SPECIES)
            + f"; at n = {a.equal_n}: " + ", ".join(f"{s} {tot[f'diag_{s}_{eq}']:,}" for s in SPECIES),
            "fixed differences: " + ", ".join(f"{p}-{q} {tot[f'fixed_{p}_{q}']:,}" for p, q in PAIRS)] + text
    open(os.path.join(a.out, "sv_summary.txt"), "w").write("\n".join(text) + "\n")
    print("\n".join(text))


if __name__ == "__main__":
    main()
