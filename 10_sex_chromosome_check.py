#!/usr/bin/env python3
# =============================================================================
# 10_sex_chromosome_check.py
#
# Tests, quantitatively, whether the six short/low-BUSCO assemblies are short
# because their haplotype lacks most of the Z chromosome.
#
# Three checks, each independent of the others:
#
#   1. CHROMOSOME CONTENT   (needs only final/*.pseudo_chr.fasta.fai)
#      chr_Z, chr_W and chr_MT length per assembly, against the panel median.
#      A haplotype missing the Z shows chr_Z far below the panel's typical
#      ~76 Mb while its partner is normal.
#
#   2. BUSCOs ON chr_Z      (needs only each assembly's own full_table.tsv)
#      How many Complete BUSCOs each assembly places on its own chr_Z. An
#      assembly with a truncated Z should carry correspondingly few. This is
#      the self-contained version of the test and needs no reference run.
#
#   3. MISSING-BUSCO LOCALISATION  (needs the ptarmigan reference BUSCO run,
#      i.e. 08_ptarmigan_reference_busco.sh)
#      Uses the reference run to assign every ortholog to a reference
#      chromosome, then asks where each assembly's MISSING orthologs sit. If
#      the extra missing orthologs in the six outliers are overwhelmingly
#      Z-linked, the Z explanation is established rather than inferred. Check 3
#      is skipped with a note if the reference run is absent.
#
# Pure standard library — no pandas, no plotting. Runs in seconds on a login
# node:
#     python3 10_sex_chromosome_check.py \
#         --final-dir $CLUSTER_SCRATCH/GROUSE/grouse_asm/final_ptarmigan \
#         --busco-dir $CLUSTER_SCRATCH/GROUSE/grouse_asm/qc_ptarmigan/busco \
#         --out-dir   $CLUSTER_SCRATCH/GROUSE/grouse_asm/qc_ptarmigan/sex_chromosome_check
# =============================================================================

import argparse
import glob
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

REF_LABEL_PATTERN = re.compile(r"([A-Z]+)_([A-Za-z0-9.]+)_(hap\d+|ref)")
SEX_SEQS = ("chr_Z", "chr_W", "chr_MT")


# ------------------------------------------------------------------ helpers --
def parse_name(dirname):
    m = REF_LABEL_PATTERN.fullmatch(dirname)
    return m.groups() if m else None


def read_fai(path):
    """sequence -> length."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


def read_fulltable(path):
    """
    busco_id -> (status, sequence). Finds the header line rather than assuming
    it is on line 3, because BUSCO's comment block is not guaranteed to be two
    lines. A duplicated BUSCO appears on several rows; the first is kept, which
    is enough for presence/absence and chromosome assignment.
    """
    out = {}
    with open(path) as fh:
        header_idx = None
        for i, line in enumerate(fh):
            if line.startswith("# Busco id"):
                header_idx = i
                break
        if header_idx is None:
            return None
        for line in fh:                    # file position is after the header
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            bid, status = f[0], f[1]
            seq = f[2].split(":")[0] if len(f) > 2 and f[2] else None
            if bid not in out:
                out[bid] = (status, seq)
    return out


def find_assemblies(busco_dir, final_dir):
    found = []
    for path in sorted(glob.glob(os.path.join(busco_dir, "*"))):
        if not os.path.isdir(path):
            continue
        parsed = parse_name(os.path.basename(path))
        if not parsed:
            continue
        species, sample, hap = parsed
        name = os.path.basename(path)
        tables = sorted(glob.glob(os.path.join(path, "run_*", "full_table.tsv")))
        found.append({
            "name": name, "species": species, "sample": sample, "hap": hap,
            "full_table": tables[0] if tables else None,
            "fai": os.path.join(final_dir, f"{name}.pseudo_chr.fasta.fai"),
        })
    return found


def mb(x):
    return "NA" if x is None else f"{x / 1e6:.2f}"


# --------------------------------------------------------------------- main --
def main():
    p = argparse.ArgumentParser(
        description="Test whether short assemblies are short because they lack the Z.")
    p.add_argument("--final-dir", required=True)
    p.add_argument("--busco-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ref-name", default="LAGMUT_bLagMut1_ref",
                   help="directory name of the reference-genome BUSCO run (check 3)")
    p.add_argument("--z-deficit-frac", type=float, default=0.5,
                   help="flag an assembly whose chr_Z is below this fraction of the "
                        "panel median chr_Z (default 0.5)")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    assemblies = find_assemblies(args.busco_dir, args.final_dir)
    grouse = [a for a in assemblies if a["hap"] != "ref"]
    if not grouse:
        sys.exit(f"ERROR: no <SPECIES>_<SAMPLE>_hapN directories under {args.busco_dir}")
    print(f">>> {len(grouse)} haplotype assemblies "
          f"({len(assemblies) - len(grouse)} reference)")

    # =========================================================== CHECK 1 ====
    print("\n>>> [1/3] Chromosome content from the FASTA index")
    for a in assemblies:
        a["fai_data"] = read_fai(a["fai"]) if os.path.exists(a["fai"]) else None
        if a["fai_data"] is None:
            print(f"    ! {a['name']}: no .fai at {a['fai']}")

    with_fai = [a for a in assemblies if a["fai_data"]]
    for a in with_fai:
        d = a["fai_data"]
        a["total"] = sum(d.values())
        for s in SEX_SEQS:
            a[s] = d.get(s)
        a["autosome"] = sum(v for k, v in d.items()
                            if k.startswith("chr_") and k not in SEX_SEQS)

    z_vals = [a["chr_Z"] for a in with_fai if a["hap"] != "ref" and a.get("chr_Z")]
    z_median = statistics.median(z_vals) if z_vals else None
    if z_median:
        print(f"    panel median chr_Z: {z_median / 1e6:.2f} Mb (n={len(z_vals)})")

    path1 = os.path.join(args.out_dir, "chromosome_content.tsv")
    with open(path1, "w") as fh:
        fh.write("assembly\tspecies\tsample\thap\ttotal_Mb\tautosome_Mb\t"
                 "chr_Z_Mb\tchr_W_Mb\tchr_MT_kb\tchr_Z_vs_median\tz_deficient\n")
        for a in with_fai:
            z = a.get("chr_Z")
            ratio = (z / z_median) if (z and z_median) else None
            a["z_ratio"] = ratio
            a["z_deficient"] = bool(ratio is not None and ratio < args.z_deficit_frac)
            fh.write("\t".join([
                a["name"], a["species"], a["sample"], a["hap"],
                mb(a.get("total")), mb(a.get("autosome")),
                mb(z), mb(a.get("chr_W")),
                "NA" if a.get("chr_MT") is None else f"{a['chr_MT'] / 1e3:.1f}",
                "NA" if ratio is None else f"{ratio:.3f}",
                "yes" if a["z_deficient"] else "no",
            ]) + "\n")
    print(f"    wrote {path1}")

    z_def = [a for a in with_fai if a["hap"] != "ref" and a.get("z_deficient")]
    if z_def:
        print(f"    Z-deficient assemblies (chr_Z < {args.z_deficit_frac:g}x panel median):")
        for a in sorted(z_def, key=lambda v: v["name"]):
            print(f"      {a['name']:<24} chr_Z {mb(a['chr_Z']):>8} Mb "
                  f"({a['z_ratio']:.2f}x)   chr_W {mb(a.get('chr_W')):>8} Mb   "
                  f"total {mb(a['total']):>8} Mb")
    else:
        print("    no assembly has a markedly short chr_Z")

    # =========================================================== CHECK 2 ====
    print("\n>>> [2/3] Complete BUSCOs placed on each assembly's own chr_Z")
    for a in assemblies:
        a["ft"] = read_fulltable(a["full_table"]) if a["full_table"] else None
        if a["ft"] is None and a["full_table"]:
            print(f"    ! {a['name']}: unreadable full_table")

    with_ft = [a for a in assemblies if a["ft"]]
    for a in with_ft:
        per_chr = Counter()
        missing = 0
        for status, seq in a["ft"].values():
            if status == "Missing":
                missing += 1
            elif seq:
                per_chr[seq] += 1
        a["per_chr"] = per_chr
        a["n_missing"] = missing
        a["n_total"] = len(a["ft"])
        a["z_buscos"] = per_chr.get("chr_Z", 0)
        a["w_buscos"] = per_chr.get("chr_W", 0)

    grouse_ft = [a for a in with_ft if a["hap"] != "ref"]
    zb = [a["z_buscos"] for a in grouse_ft]
    zb_median = statistics.median(zb) if zb else 0
    miss = [a["n_missing"] for a in grouse_ft]
    miss_median = statistics.median(miss) if miss else 0
    print(f"    panel median: {zb_median:.0f} Complete BUSCOs on chr_Z, "
          f"{miss_median:.0f} Missing overall")

    path2 = os.path.join(args.out_dir, "busco_by_chromosome.tsv")
    with open(path2, "w") as fh:
        fh.write("assembly\tspecies\tsample\thap\tn_markers\tcomplete_on_chr_Z\t"
                 "complete_on_chr_W\tmissing_total\tmissing_excess_vs_median\t"
                 "chr_Z_Mb\tz_deficient\n")
        for a in sorted(with_ft, key=lambda v: v["name"]):
            excess = a["n_missing"] - miss_median if a["hap"] != "ref" else 0
            a["missing_excess"] = excess
            fh.write("\t".join([
                a["name"], a["species"], a["sample"], a["hap"],
                str(a["n_total"]), str(a["z_buscos"]), str(a["w_buscos"]),
                str(a["n_missing"]), f"{excess:.0f}",
                mb(a.get("chr_Z")), "yes" if a.get("z_deficient") else "no",
            ]) + "\n")
    print(f"    wrote {path2}")

    low_z = [a for a in grouse_ft if a["z_buscos"] < 0.5 * zb_median]
    if low_z:
        print("    assemblies with few BUSCOs on chr_Z:")
        for a in sorted(low_z, key=lambda v: v["name"]):
            print(f"      {a['name']:<24} {a['z_buscos']:>4} on chr_Z "
                  f"(median {zb_median:.0f})   missing {a['n_missing']:>4} "
                  f"(+{a['missing_excess']:.0f} vs median)")
    else:
        print("    every assembly carries a typical number of chr_Z BUSCOs")

    # =========================================================== CHECK 3 ====
    print("\n>>> [3/3] Where the MISSING orthologs sit on the reference genome")
    ref = next((a for a in assemblies if a["name"] == args.ref_name), None)
    if ref is None or not ref.get("ft"):
        print(f"    skipped: no usable BUSCO run for '{args.ref_name}' under {args.busco_dir}")
        print("             Run 08_ptarmigan_reference_busco.sh, then rerun this script.")
    else:
        ref_chr = {bid: seq for bid, (status, seq) in ref["ft"].items()
                   if status != "Missing" and seq}
        n_ref_z = sum(1 for v in ref_chr.values() if v == "chr_Z")
        print(f"    reference run localises {len(ref_chr)} orthologs; "
              f"{n_ref_z} of them sit on chr_Z")

        path3 = os.path.join(args.out_dir, "missing_busco_localisation.tsv")
        with open(path3, "w") as fh:
            fh.write("assembly\tspecies\tsample\thap\tmissing_total\tmissing_localised\t"
                     "missing_on_ref_chr_Z\tpct_of_missing_on_chr_Z\t"
                     "pct_of_ref_chr_Z_orthologs_missing\tz_deficient\n")
            for a in sorted(grouse_ft, key=lambda v: v["name"]):
                missing_ids = [b for b, (s, _) in a["ft"].items() if s == "Missing"]
                loc = [ref_chr[b] for b in missing_ids if b in ref_chr]
                on_z = sum(1 for c in loc if c == "chr_Z")
                a["missing_on_z"] = on_z
                a["missing_loc"] = len(loc)
                pct_missing = (100 * on_z / len(loc)) if loc else 0.0
                pct_z_lost = (100 * on_z / n_ref_z) if n_ref_z else 0.0
                fh.write("\t".join([
                    a["name"], a["species"], a["sample"], a["hap"],
                    str(a["n_missing"]), str(len(loc)), str(on_z),
                    f"{pct_missing:.1f}", f"{pct_z_lost:.1f}",
                    "yes" if a.get("z_deficient") else "no",
                ]) + "\n")
        print(f"    wrote {path3}")

        flagged = [a for a in grouse_ft if a.get("z_deficient") or a["missing_excess"] > 100]
        normal = [a for a in grouse_ft if a not in flagged]
        if flagged and normal:
            f_on_z = statistics.median([a["missing_on_z"] for a in flagged])
            n_on_z = statistics.median([a["missing_on_z"] for a in normal])
            print(f"    median missing orthologs on reference chr_Z:")
            print(f"      flagged assemblies (n={len(flagged)}): {f_on_z:.0f}")
            print(f"      the rest          (n={len(normal)}): {n_on_z:.0f}")
            excess = statistics.median([a["missing_excess"] for a in flagged])
            attributable = (100 * (f_on_z - n_on_z) / excess) if excess else 0
            print(f"    excess missing in flagged assemblies: {excess:.0f} orthologs, "
                  f"of which {f_on_z - n_on_z:.0f} are Z-linked "
                  f"({attributable:.0f}%)")
            print()
            if attributable >= 70:
                print("    VERDICT: the excess missing orthologs are predominantly Z-linked.")
                print("             Loss of Z-linked sequence explains the shortfall.")
            elif attributable >= 30:
                print("    VERDICT: Z-linked loss accounts for part of the shortfall but")
                print("             not most of it. Inspect the per-chromosome breakdown.")
            else:
                print("    VERDICT: the excess missing orthologs are NOT mainly Z-linked.")
                print("             The Z explanation does not hold; look elsewhere.")
        else:
            print("    no clear two-group split to compare — inspect the TSV directly")

    print(f"\n>>> Done. Output under {args.out_dir}")


if __name__ == "__main__":
    main()
