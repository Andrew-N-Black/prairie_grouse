#!/usr/bin/env python3
"""
validate_snp_markers.py -- check the species-diagnostic Z SNPs (step 22)
in the higher-coverage shotgun birds, the same way 14h checks the indels.

Each bird's reads are aligned only to short windows around the SNPs (LEPC
reference); reads with more than a few mismatches are dropped, so reads from
diverged copies (W gametologs, paralogs) do not count. Each bird is then
genotyped at every SNP and compared with its species' expected allele.
Females carry one Z, so a female with two alleles at a SNP flags a W copy, a
duplication or contamination.

Sub-commands:
  regions    read the step-22 SNP table, write merged windows (BED), the site
             list in window coordinates and the expected alleles
  count      per-site base counts for one bird (samtools mpileup on stdin)
  summarize  genotype calls, per-SNP accuracy by species, per-bird Z ancestry
  refcheck   which chromosome each window matches in another reference (PAF)
Standard library only.
"""
import argparse, collections, csv, glob, os, re, sys

SPECIES = ["LEPC", "GRPC", "STGR"]


def pair_species(s):
    """'LEPC_GRPC', 'LEPC-GRPC', 'GRPC vs STGR' -> ['LEPC', 'GRPC'] in order"""
    found = [(m.start(), m.group(0)) for m in re.finditer("|".join(SPECIES), s.upper())]
    return [x for _, x in sorted(found)]


# ------------------------------------------------------------------ regions
def read_table(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip()]
    sep = "\t" if "\t" in lines[0] else None              # tab, else any whitespace
    hdr = lines[0].split(sep)
    return [dict(zip(hdr, l.split(sep))) for l in lines[1:]]


def cmd_regions(a):
    rows = read_table(a.table)
    if not rows:
        sys.exit("empty table")
    wide = all("allele_" + x in rows[0] for x in SPECIES)
    need = [a.col_contig, a.col_pos, a.col_pair] + ([] if wide else [a.col_a, a.col_b])
    for c in need:
        if c not in rows[0]:
            sys.exit("column '%s' not in %s (columns: %s)" % (c, a.table, ", ".join(rows[0])))
    sites = []
    for r in rows:
        if a.col_keep and str(r.get(a.col_keep, "")).strip() not in a.keep_values.split(","):
            continue
        if "type" in r and r["type"] != "SNP":
            continue
        pos = int(float(r[a.col_pos])) - a.pos_offset
        if wide:
            # step-22 format: one row per SNP; allele_<SP> = R (ref) / A (alt);
            # pairs listed in COL_PAIR, e.g. "GRPC-STGR,LEPC-STGR"
            code = {"R": r.get("ref", "").upper(), "A": r.get("alt", "").upper()}
            al = {x: code.get(r["allele_" + x].strip().upper(), "") for x in SPECIES}
            for pr in re.split(r"[,;]", r[a.col_pair]):
                sp = pair_species(pr)
                if len(sp) != 2:
                    continue
                aa, ab = al[sp[0]], al[sp[1]]
                if len(aa) == 1 and len(ab) == 1 and aa in "ACGT" and ab in "ACGT" and aa != ab:
                    sites.append((r[a.col_contig], pos, sp[0], sp[1], aa, ab))
            continue
        sp = pair_species(r[a.col_pair])
        if len(sp) != 2:
            continue
        al_a, al_b = r[a.col_a].strip().upper(), r[a.col_b].strip().upper()
        if not (len(al_a) == 1 and len(al_b) == 1 and al_a in "ACGT" and al_b in "ACGT" and al_a != al_b):
            continue
        sites.append((r[a.col_contig], pos, sp[0], sp[1], al_a, al_b))
    # the same position may be listed for several pairs: keep one row per pair
    sites = sorted(set(sites))
    # merged windows
    bycontig = collections.defaultdict(list)
    for s in sites:
        bycontig[s[0]].append(s[1])
    regions = []
    for c, ps in bycontig.items():
        for p in sorted(set(ps)):
            s, e = max(0, p - a.flank), p + a.flank + 1
            if regions and regions[-1][0] == c and s <= regions[-1][2]:
                regions[-1][2] = max(regions[-1][2], e)
            else:
                regions.append([c, s, e])
    with open(a.out_prefix + ".regions.bed", "w") as f:
        for c, s, e in regions:
            f.write("%s\t%d\t%d\n" % (c, s, e))
    # site positions in window coordinates (samtools faidx names windows "ctg:s+1-e")
    rstart = collections.defaultdict(list)
    for c, s, e in regions:
        rstart[c].append((s, e))
    with open(a.out_prefix + ".sites.tsv", "w") as f, open(a.out_prefix + ".sites_window.bed", "w") as b:
        f.write("site\twindow\twpos\tcontig\tpos\tsp_a\tsp_b\tallele_a\tallele_b\n")
        seen = set()
        for c, p, sa, sb, aa, ab in sites:
            s, e = next((s, e) for s, e in rstart[c] if s <= p < e)
            win = "%s:%d-%d" % (c, s + 1, e)
            site = "%s:%d" % (c, p + 1)
            f.write("%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n" % (site, win, p - s + 1, c, p + 1, sa, sb, aa, ab))
            if (win, p - s) not in seen:
                b.write("%s\t%d\t%d\n" % (win, p - s, p - s + 1)); seen.add((win, p - s))
    pc = collections.Counter("%s_%s" % (s[2], s[3]) for s in sites)
    sys.stderr.write("regions: %d SNP rows (%s) in %d windows (%.1f Mb)\n" % (
        len(sites), ", ".join("%s %d" % kv for kv in sorted(pc.items())), len(regions),
        sum(e - s for _, s, e in regions) / 1e6))


# -------------------------------------------------------------------- count
def cmd_count(a):
    out = open(a.out, "w")
    out.write("window\twpos\tA\tC\tG\tT\n")
    for line in sys.stdin:
        t = line.rstrip("\n").split("\t")
        if len(t) < 5:
            continue
        ref, bases = t[2].upper(), t[4]
        n = collections.Counter()
        i = 0
        while i < len(bases):
            ch = bases[i]
            if ch == "^":
                i += 2; continue
            if ch in "+-":
                j = i + 1
                while j < len(bases) and bases[j].isdigit():
                    j += 1
                i = j + int(bases[i + 1:j]); continue
            if ch in ".,":
                n[ref] += 1
            elif ch.upper() in "ACGT":
                n[ch.upper()] += 1
            i += 1
        out.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (t[0], t[1], n["A"], n["C"], n["G"], n["T"]))
    out.close()


# ---------------------------------------------------------------- summarize
def cmd_summarize(a):
    meta = {r["sample"]: r for r in csv.DictReader(open(a.samples), delimiter="\t")}
    sites = list(csv.DictReader(open(a.sites), delimiter="\t"))
    bysite = collections.defaultdict(list)
    for s in sites:
        bysite[(s["window"], s["wpos"])].append(s)
    per_snp = collections.defaultdict(collections.Counter)   # (site, pair, species) -> status
    bird = collections.defaultdict(lambda: collections.Counter())
    calls_out = open(a.out_prefix + ".bird_ancestry.tsv", "w")
    for f in sorted(glob.glob(os.path.join(a.counts_dir, "*.counts.tsv"))):
        smp = os.path.basename(f)[:-len(".counts.tsv")]
        sp = meta.get(smp, {}).get("species", "NA") or "NA"
        sex = meta.get(smp, {}).get("sex", "NA") or "NA"
        for r in csv.DictReader(open(f), delimiter="\t"):
            n = {b: int(r[b]) for b in "ACGT"}
            tot = sum(n.values())
            for s in bysite.get((r["window"], r["wpos"]), []):
                pair = s["sp_a"] + "_" + s["sp_b"]
                if tot < a.min_reads:
                    st = "no_data"
                else:
                    called = [b for b in "ACGT" if n[b] >= a.min_alt and n[b] >= a.min_frac * tot]
                    other = [b for b in called if b not in (s["allele_a"], s["allele_b"])]
                    if other:
                        st = "third_allele"
                    elif len(called) == 2:
                        st = "two_alleles_FEMALE" if sex == "F" else "two_alleles"
                    elif called == [s["allele_a"]]:
                        st = "allele_" + s["sp_a"]
                    else:
                        st = "allele_" + s["sp_b"]
                    # per-bird ancestry counts by pair (alleles, not genotypes)
                    na_, nb_ = n[s["allele_a"]], n[s["allele_b"]]
                    if na_ + nb_:
                        bird[(smp, pair)][s["sp_a"]] += na_ / (na_ + nb_)
                        bird[(smp, pair)][s["sp_b"]] += nb_ / (na_ + nb_)
                        bird[(smp, pair)]["n"] += 1
                if sp in (s["sp_a"], s["sp_b"]):
                    exp = "allele_" + sp
                    st2 = "expected" if st == exp else ("no_data" if st == "no_data" else
                          st if st.startswith(("two_", "third")) else "unexpected")
                else:
                    st2 = "female_two_alleles" if st == "two_alleles_FEMALE" else "other_species"
                per_snp[(s["site"], pair, sp)][st2] += 1
    # per-bird Z ancestry
    calls_out.write("sample\tspecies\tsex\tpair\tn_snps\tfrac_first\tfrac_second\n")
    for (smp, pair), c in sorted(bird.items()):
        a1, a2 = pair.split("_")
        calls_out.write("%s\t%s\t%s\t%s\t%d\t%.3f\t%.3f\n" % (smp, meta.get(smp, {}).get("species", "NA"),
                        meta.get(smp, {}).get("sex", "NA"), pair, c["n"], c[a1] / c["n"], c[a2] / c["n"]))
    calls_out.close()
    # per-SNP table and pass/fail
    snps = collections.defaultdict(dict)
    for (site, pair, sp), c in per_snp.items():
        snps[(site, pair)][sp] = c
    keep = collections.Counter(); tot = collections.Counter()
    with open(a.out_prefix + ".snp_validation.tsv", "w") as f:
        f.write("site\tpair\tn_first\texp_first\tn_second\texp_second\tunexpected\ttwo_alleles\t"
                "female_two_alleles\tthird_allele\tpass\n")
        for (site, pair), d in sorted(snps.items()):
            a1, a2 = pair.split("_")
            def called(sp):
                c = d.get(sp, collections.Counter())
                return sum(v for k, v in c.items() if k != "no_data"), c["expected"]
            n1, e1 = called(a1); n2, e2 = called(a2)
            allc = collections.Counter()
            for c in d.values():
                allc.update(c)
            fem2 = allc["two_alleles_FEMALE"] + allc["female_two_alleles"]
            ok = (n1 >= a.min_birds and n2 >= a.min_birds and e1 >= a.min_acc * n1 and e2 >= a.min_acc * n2
                  and fem2 == 0 and allc["third_allele"] == 0)
            tot[pair] += 1; keep[pair] += ok
            f.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                site, pair, n1, e1, n2, e2, allc["unexpected"], allc["two_alleles"], fem2,
                allc["third_allele"], int(ok)))
    sys.stderr.write("SNPs passing (>= %.0f%% expected in both species, >= %d birds each, no female "
                     "with two alleles): %s\n" % (100 * a.min_acc, a.min_birds,
                     "; ".join("%s %d of %d" % (p, keep[p], tot[p]) for p in sorted(tot))))


# ----------------------------------------------------------------- refcheck
def cmd_refcheck(a):
    names = {}
    if a.chrmap and os.path.exists(a.chrmap):
        for line in open(a.chrmap):
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) > 6:
                names[t[6]] = t[0]
    best = {}
    for line in open(a.paf):
        t = line.split("\t")
        if int(t[10]) == 0:
            continue
        sc = int(t[9])
        if t[0] not in best or sc > best[t[0]][0]:
            best[t[0]] = (sc, names.get(t[5], t[5]), int(t[9]) / int(t[10]))
    nwin = sum(1 for _ in open(a.regions))
    c = collections.Counter(v[1] for v in best.values())
    with open(a.out, "w") as f:
        f.write("window\tbest_chrom\tidentity\n")
        for k, v in sorted(best.items()):
            f.write("%s\t%s\t%.3f\n" % (k, v[1], v[2]))
    sys.stderr.write("refcheck %s: %d of %d windows aligned; best hits: %s\n" % (
        a.label, len(best), nwin, ", ".join("%s %d" % kv for kv in c.most_common(6))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("regions")
    p.add_argument("--table", required=True); p.add_argument("--out-prefix", required=True)
    p.add_argument("--col-contig", default="contig"); p.add_argument("--col-pos", default="pos")
    p.add_argument("--col-pair", default="pair"); p.add_argument("--col-a", default="allele_a")
    p.add_argument("--col-b", default="allele_b")
    p.add_argument("--col-keep", default="", help="optional filter column (e.g. diagnostic)")
    p.add_argument("--keep-values", default="1,TRUE,True,yes")
    p.add_argument("--pos-offset", type=int, default=1, help="1 if positions are 1-based, 0 if 0-based")
    p.add_argument("--flank", type=int, default=300)
    p = sub.add_parser("count"); p.add_argument("--out", required=True)
    p = sub.add_parser("summarize")
    p.add_argument("--samples", required=True); p.add_argument("--sites", required=True)
    p.add_argument("--counts-dir", required=True); p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-reads", type=int, default=4); p.add_argument("--min-alt", type=int, default=2)
    p.add_argument("--min-frac", type=float, default=0.15); p.add_argument("--min-birds", type=int, default=10)
    p.add_argument("--min-acc", type=float, default=0.95)
    p = sub.add_parser("refcheck")
    p.add_argument("--paf", required=True); p.add_argument("--regions", required=True)
    p.add_argument("--chrmap", default=""); p.add_argument("--label", required=True); p.add_argument("--out", required=True)
    a = ap.parse_args()
    {"regions": cmd_regions, "count": cmd_count, "summarize": cmd_summarize, "refcheck": cmd_refcheck}[a.cmd](a)


if __name__ == "__main__":
    main()
