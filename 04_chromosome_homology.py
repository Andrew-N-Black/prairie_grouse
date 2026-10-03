#!/usr/bin/env python3
"""
Chromosome correspondence between two reference genomes, from a whole-genome
alignment. Companion to 04_chromosome_homology.sh.

WHY THIS EXISTS
    Reference assemblies use two incompatible chromosome-naming conventions:

      homology naming  chrN means "the chromosome the field calls N"
      size-rank naming chrN means "the Nth longest scaffold", no homology claim

    Sanger/VGP curated assemblies (SUPER_1, SUPER_2, ...) use the second. Treating
    those numbers as homology produces chromosome names that silently mean
    different things in different trees, and nothing fails while it happens.

    This script computes the correspondence instead of assuming it, and reports
    the cases that matter for karyotype work:

      one_to_one   query chromosome corresponds to a single target chromosome
      fusion       ONE query chromosome covers TWO OR MORE target chromosomes
                   -> the target split them, or the query lineage fused them
      split        SEVERAL query chromosomes divide ONE target chromosome
                   -> the target fused them, or the query lineage split them
      ambiguous    no target accounts for enough of the alignment to call

    "fusion" and "split" are descriptions of the correspondence, not of history:
    which lineage changed is a question for an outgroup, not for two genomes.

METHOD
    Input is PAF from `minimap2 -x asm20 <target.fa> <query.fa>`. For each query
    chromosome, matching bases (PAF column 10) are summed per target chromosome.
    Summing MATCHES rather than counting alignment records is deliberate: record
    counts are dominated by short repeat hits, which are exactly the alignments
    that cross between non-homologous chromosomes.

    A target is "substantial" for a query chromosome if it holds at least
    --min-frac of that chromosome's total matched bases. Reciprocal coverage is
    also reported: what fraction of the TARGET chromosome's length the shared
    alignment spans, which is what distinguishes a whole-arm correspondence from
    a shared repeat family.

OUTPUTS
    <prefix>.correspondence.tsv   every query x target pair above the threshold
    <prefix>.summary.tsv          one row per query chromosome, with the verdict
    <prefix>.rename_map.tsv       query accession -> homology-based chr_ label,
                                  for step 05 to consume
"""

import argparse
import re
import sys
from collections import defaultdict


# --------------------------------------------------------------------- input --
def read_chr_map(path):
    """accession -> chromosome name, as written by step 03."""
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                out[parts[0]] = parts[1]
    return out


def read_fai(path):
    """sequence name -> length."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = int(f[1])
    return out


def chr_sort_key(name):
    """Numbered chromosomes ascending, then Z, W, MT, then anything else."""
    m = re.fullmatch(r"(\d+)", name)
    if m:
        return (0, int(m.group(1)), "")
    return (1, {"Z": 0, "W": 1, "MT": 2}.get(name, 3), name)


# ------------------------------------------------------------------ analysis --
def parse_paf(paf_path, qmap, tmap, min_block, min_mapq):
    """
    Sum matched bases per (query chromosome, target chromosome).

    Only sequences present in the respective chromosome maps are considered, so
    unplaced scaffolds on either side are ignored rather than contributing noise.
    """
    matched = defaultdict(lambda: defaultdict(int))
    q_total = defaultdict(int)
    # Matched-base-weighted mean position of each query chromosome ON the target.
    # Used to order split suffixes by position rather than by size — see
    # build_rename_labels. Weighted by matched bases so a single stray repeat
    # alignment at the far end of the chromosome cannot move it.
    tpos_num = defaultdict(lambda: defaultdict(float))
    tpos_den = defaultdict(lambda: defaultdict(float))
    skipped_unmapped = 0
    n_records = 0

    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            n_records += 1
            qname, tname = f[0], f[5]
            qchr, tchr = qmap.get(qname), tmap.get(tname)
            if qchr is None or tchr is None:
                skipped_unmapped += 1
                continue
            try:
                n_match = int(f[9])
                block = int(f[10])
                mapq = int(f[11])
            except ValueError:
                continue
            if block < min_block or mapq < min_mapq:
                continue
            matched[qchr][tchr] += n_match
            q_total[qchr] += n_match
            try:
                tmid = (int(f[7]) + int(f[8])) / 2.0
            except ValueError:
                tmid = None
            if tmid is not None:
                tpos_num[qchr][tchr] += tmid * n_match
                tpos_den[qchr][tchr] += n_match

    tpos = {q: {t: (tpos_num[q][t] / tpos_den[q][t])
                for t in tpos_num[q] if tpos_den[q][t] > 0}
            for q in tpos_num}
    return matched, q_total, n_records, skipped_unmapped, tpos


def classify(matched, q_total, t_len, min_frac, min_target_cov,
             high_conc, min_matched_bp):
    """
    Per query chromosome, the target chromosomes it corresponds to, largest first.

    A target must hold >= min_frac of the QUERY chromosome's matched bases, and
    then clear EITHER of two routes:

      A. reciprocal   the alignment spans >= min_target_cov of the TARGET
                      chromosome's length
      B. concentration the query sends >= high_conc of ALL its matched bases to
                      this one target, on at least min_matched_bp of sequence

    Route A alone is not enough, in BOTH directions, and both failures were seen
    on real data rather than imagined:

      Without the reciprocal test, dispersed repeat alignments spread a query
      chromosome's matched bases evenly across every target, so with five targets
      each takes ~20%, clears min_frac, and yields a confident five-way "fusion"
      out of pure noise. A uniform random PAF produces exactly that.

      With ONLY the reciprocal test, a query chromosome that is genuinely just a
      PIECE of a larger target is rejected, because covering a third of a 149 Mb
      chromosome cannot produce high reciprocal coverage by construction. On the
      ptarmigan-vs-chicken run this silently hid both real splits: chicken 2 =
      ptarmigan 3 + 7 (98.0% of its length reconstructed) and chicken 4 =
      ptarmigan 4 + 13 (98.9%) — the second being GGA4p, the exact join the
      grouse Hi-C had independently flagged as unsupported.

    Route B rescues those without readmitting the noise, because the two cases
    differ sharply in concentration: the scattered-noise simulation never sent
    more than 22% of a chromosome's matched bases to any single target, while a
    genuine partial correspondence sends essentially all of them to one.
    min_matched_bp is a floor so that a tiny chromosome with a handful of
    alignments cannot reach 100% by accident.
    """
    calls = {}
    for qchr, targets in matched.items():
        total = q_total[qchr]
        if total <= 0:
            continue
        ranked = sorted(targets.items(), key=lambda kv: -kv[1])
        substantial = []
        for t, b in ranked:
            frac_q = b / total
            tl = t_len.get(t, 0)
            frac_t = (b / tl) if tl else 0.0
            if frac_q < min_frac:
                continue
            by_reciprocal = frac_t >= min_target_cov
            by_concentration = frac_q >= high_conc and b >= min_matched_bp
            if by_reciprocal or by_concentration:
                substantial.append((t, b, frac_q))
        if not substantial:
            # Nothing cleared either route: keep the best hit so the row is
            # still informative, but it will be called ambiguous below.
            t, b = ranked[0]
            substantial = [(t, b, b / total)]
            calls[qchr] = {"total": total, "targets": substantial,
                           "all_ranked": ranked, "passed": False}
            continue
        calls[qchr] = {"total": total, "targets": substantial,
                       "all_ranked": ranked, "passed": True}
    return calls


def assign_verdicts(calls, min_frac):
    """
    one_to_one / fusion / split / ambiguous.

    'split' can only be decided globally — it depends on whether some OTHER query
    chromosome also claims the same target — so it is applied after the per-query
    pass.
    """
    # Which query chromosomes claim each target substantially?
    claimants = defaultdict(list)
    for qchr, info in calls.items():
        if not info.get("passed"):
            continue
        for t, _b, frac in info["targets"]:
            if frac >= min_frac:
                claimants[t].append(qchr)

    verdicts = {}
    for qchr, info in calls.items():
        tg = [(t, b, f) for t, b, f in info["targets"] if f >= min_frac] \
            if info.get("passed") else []
        if not tg:
            verdicts[qchr] = "ambiguous"
        elif len(tg) >= 2:
            verdicts[qchr] = "fusion"
        else:
            only = tg[0][0]
            verdicts[qchr] = "split" if len(claimants.get(only, [])) >= 2 else "one_to_one"
    return verdicts, claimants


def build_rename_labels(calls, verdicts, claimants, q_len, min_frac, tpos=None):
    """
    A chr_ label per query chromosome that says what it corresponds to.

      one_to_one -> chr_<target>                 e.g. chr_1
      fusion     -> chr_<t1>_<t2>                e.g. chr_6_8
      split      -> chr_<target>a / b / c        e.g. chr_2a, chr_2b,
                    lettered by POSITION along the target chromosome, so that
                    chr_2a is the segment nearest the target's start. Ordering
                    these by size instead would read backwards to anyone who
                    expects chr_4a to mean 4p: on the real ptarmigan run the
                    GGA4p microchromosome is the SMALLER of chicken chr4's two
                    pieces but the FIRST along it.
      ambiguous  -> chr_<query>_unresolved

    Sex chromosomes and MT short-circuit to themselves: they are identified by
    name in both assemblies and should never be renamed from alignment evidence.
    """
    labels = {}
    for qchr in calls:
        if qchr in ("Z", "W", "MT"):
            labels[qchr] = f"chr_{qchr}"

    # Split groups: order the claimants by their own length, descending, so the
    # a/b/c suffixes are stable and meaningful rather than dictionary order.
    split_suffix = {}
    for target, qs in claimants.items():
        qs_real = [q for q in qs if verdicts.get(q) == "split"]
        if len(qs_real) < 2:
            continue
        # Primary key: mean position on the TARGET chromosome, ascending.
        # Falls back to descending length when positions are unavailable, which
        # keeps the function usable without a PAF-derived position table.
        if tpos:
            ordered = sorted(qs_real,
                             key=lambda q: (tpos.get(q, {}).get(target, float("inf")),
                                            -q_len.get(q, 0)))
        else:
            ordered = sorted(qs_real, key=lambda q: -q_len.get(q, 0))
        for i, q in enumerate(ordered):
            split_suffix[q] = (target, "abcdefghijklmnop"[i] if i < 16 else str(i))

    for qchr, info in calls.items():
        if qchr in labels:
            continue
        v = verdicts[qchr]
        tg = [t for t, _b, f in info["targets"] if f >= min_frac] \
            if info.get("passed") else []
        if v == "one_to_one":
            labels[qchr] = f"chr_{tg[0]}"
        elif v == "fusion":
            parts = sorted(tg, key=chr_sort_key)
            labels[qchr] = "chr_" + "_".join(parts)
        elif v == "split":
            target, suf = split_suffix.get(qchr, (tg[0], ""))
            labels[qchr] = f"chr_{target}{suf}"
        else:
            labels[qchr] = f"chr_u{qchr}"
    return labels


# ---------------------------------------------------------------------- main --
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="minimap2 PAF, query=new reference, target=chicken")
    p.add_argument("--query-chr-map", required=True, help="accession -> name, new reference")
    p.add_argument("--target-chr-map", required=True, help="accession -> name, chicken")
    p.add_argument("--query-fai", required=True)
    p.add_argument("--target-fai", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--query-label", default="query")
    p.add_argument("--target-label", default="target")
    p.add_argument("--min-frac", type=float, default=0.15,
                   help="target must hold this fraction of a query chromosome's "
                        "matched bases to count as substantial (default 0.15)")
    p.add_argument("--high-concentration", type=float, default=0.80,
                   help="alternative route: a target also counts if the query "
                        "sends at least this fraction of ALL its matched bases "
                        "there (default 0.80). This is what lets a chromosome "
                        "that is only PART of a larger target be recognised, "
                        "which is how splits are detected.")
    p.add_argument("--min-matched-bp", type=int, default=50000,
                   help="floor on matched bases for the concentration route, so "
                        "a tiny chromosome cannot reach 100%% by accident "
                        "(default 50000)")
    p.add_argument("--min-target-cov", type=float, default=0.15,
                   help="AND the alignment must span this fraction of the target "
                        "chromosome's length (default 0.15). This is what rejects "
                        "dispersed repeat hits, which otherwise fabricate "
                        "multi-way fusions out of noise.")
    p.add_argument("--min-block", type=int, default=5000,
                   help="ignore alignment blocks shorter than this (default 5000)")
    p.add_argument("--min-mapq", type=int, default=10)
    args = p.parse_args()

    qmap_acc = read_chr_map(args.query_chr_map)
    tmap_acc = read_chr_map(args.target_chr_map)
    q_fai_acc = read_fai(args.query_fai)
    t_fai_acc = read_fai(args.target_fai)

    if not qmap_acc or not tmap_acc:
        sys.exit("ERROR: a chromosome map is empty — check the inputs.")

    # Lengths keyed by chromosome name rather than accession.
    q_len = {qmap_acc[a]: L for a, L in q_fai_acc.items() if a in qmap_acc}
    t_len = {tmap_acc[a]: L for a, L in t_fai_acc.items() if a in tmap_acc}

    matched, q_total, n_rec, n_skip, tpos = parse_paf(args.paf, qmap_acc, tmap_acc,
                                                      args.min_block, args.min_mapq)
    if not matched:
        sys.exit("ERROR: no alignments survived filtering. Check that the PAF has "
                 "query=new reference and target=chicken, and that both chromosome "
                 "maps use the same accessions as the FASTAs.")

    calls = classify(matched, q_total, t_len, args.min_frac, args.min_target_cov,
                     args.high_concentration, args.min_matched_bp)
    verdicts, claimants = assign_verdicts(calls, args.min_frac)
    labels = build_rename_labels(calls, verdicts, claimants, q_len, args.min_frac,
                                 tpos=tpos)

    # ---------------------------------------------------------------- report --
    print(f">>> PAF records read            : {n_rec}")
    print(f"    skipped (unplaced/unmapped) : {n_skip}")
    print(f"    query chromosomes with data : {len(calls)}")
    print(f"    filters: block >= {args.min_block} bp, MAPQ >= {args.min_mapq}")
    print(f"    substantial: >= {args.min_frac:.0%} of the query chromosome's matched bases,")
    print(f"      AND EITHER >= {args.min_target_cov:.0%} of the target chromosome's length")
    print(f"      OR     >= {args.high_concentration:.0%} of the query's matched bases go to "
          f"that one target (>= {args.min_matched_bp/1000:.0f} kb)")

    corr_path = f"{args.out_prefix}.correspondence.tsv"
    with open(corr_path, "w") as fh:
        fh.write(f"{args.query_label}_chr\t{args.target_label}_chr\tmatched_bp\t"
                 f"frac_of_query_matched\tfrac_of_target_len\tsubstantial\n")
        for qchr in sorted(calls, key=chr_sort_key):
            info = calls[qchr]
            for t, b in info["all_ranked"]:
                frac_q = b / info["total"]
                frac_t = b / t_len[t] if t_len.get(t) else 0.0
                if frac_q < 0.01:
                    continue
                fh.write(f"{qchr}\t{t}\t{b}\t{frac_q:.4f}\t{frac_t:.4f}\t"
                         f"{'yes' if frac_q >= args.min_frac else 'no'}\n")

    summ_path = f"{args.out_prefix}.summary.tsv"
    counts = defaultdict(int)
    with open(summ_path, "w") as fh:
        fh.write(f"{args.query_label}_chr\tlength_bp\tverdict\t{args.target_label}_chr(s)\t"
                 f"fracs\tproposed_label\n")
        for qchr in sorted(calls, key=chr_sort_key):
            info, v = calls[qchr], verdicts[qchr]
            counts[v] += 1
            tg = [(t, f) for t, _b, f in info["targets"] if f >= args.min_frac] \
                if info.get("passed") else []
            if not tg:
                tg = [(info["targets"][0][0], info["targets"][0][2])]
            fh.write(f"{qchr}\t{q_len.get(qchr, 0)}\t{v}\t"
                     f"{'+'.join(t for t, _ in tg)}\t"
                     f"{','.join(f'{f:.2f}' for _, f in tg)}\t{labels[qchr]}\n")

    # Rename map keyed by ACCESSION, which is what step 05 joins on.
    ren_path = f"{args.out_prefix}.rename_map.tsv"
    with open(ren_path, "w") as fh:
        fh.write("# accession\tchromosome_name  (homology-based; consumed by step 05)\n")
        for acc, qchr in sorted(qmap_acc.items(), key=lambda kv: chr_sort_key(kv[1])):
            lab = labels.get(qchr)
            if lab is None:
                continue
            fh.write(f"{acc}\t{lab[4:] if lab.startswith('chr_') else lab}\n")

    # ------------------------------------------------------------ to screen --
    print("")
    print("=" * 72)
    print(f"  {args.query_label} chromosome correspondence to {args.target_label}")
    print("=" * 72)
    # 'matched' is shown because the percentages alone mislead: they are a SHARE
    # of whatever aligned, so a chromosome with 118 bp of evidence can post 61%
    # and sit beside one with 476 kb posting 80%. On the real ptarmigan run that
    # cost a wrong reading of which rows were weak. The absolute number is the
    # one that says whether a percentage means anything.
    print(f"  {'query':<10} {'length':>10} {'matched':>11}  {'verdict':<11} "
          f"{'corresponds to':<18} proposed")
    print(f"  {'-'*10} {'-'*10} {'-'*11}  {'-'*11} {'-'*18} {'-'*12}")
    for qchr in sorted(calls, key=chr_sort_key):
        info, v = calls[qchr], verdicts[qchr]
        tg = [(t, f) for t, _b, f in info["targets"] if f >= args.min_frac] \
            if info.get("passed") else []
        if not tg:
            tg = [(info["targets"][0][0], info["targets"][0][2])]
        desc = " + ".join(f"{t} ({f:.0%})" for t, f in tg)
        if v == "ambiguous":
            desc = "best: " + desc + " — below reciprocal-coverage bar"
        flag = "  <<<" if v in ("fusion", "split") else ""
        matched_bp = info["total"]
        print(f"  {qchr:<10} {q_len.get(qchr,0)/1e6:>8.2f}Mb {matched_bp/1e3:>9.0f}kb  "
              f"{v:<11} {desc:<18} {labels[qchr]}{flag}")

    print("")
    print("  verdict counts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    fusions = [q for q in calls if verdicts[q] == "fusion"]
    splits = [q for q in calls if verdicts[q] == "split"]
    if fusions:
        print("")
        print(f"  FUSION — one {args.query_label} chromosome spans several in {args.target_label}:")
        for q in sorted(fusions, key=chr_sort_key):
            tg = [t for t, _b, f in calls[q]["targets"] if f >= args.min_frac]
            summed = sum(t_len.get(t, 0) for t in tg)
            ql = q_len.get(q, 0)
            print(f"    {args.query_label} {q} = {args.target_label} "
                  f"{' + '.join(sorted(tg, key=chr_sort_key))}")
            # The length check is what separates a real fusion from a shared
            # repeat family: the parts should add up to the whole.
            if ql:
                print(f"        {summed/1e6:.2f} Mb of target summed vs "
                      f"{ql/1e6:.2f} Mb query -> {100.0*summed/ql:.1f}%")
    if splits:
        print("")
        print(f"  SPLIT — several {args.query_label} chromosomes divide one in {args.target_label}:")
        seen = set()
        for q in sorted(splits, key=chr_sort_key):
            t = [tt for tt, _b, f in calls[q]["targets"] if f >= args.min_frac][0]
            if t in seen:
                continue
            seen.add(t)
            members = sorted((x for x in splits
                              if [tt for tt, _b, f in calls[x]["targets"]
                                  if f >= args.min_frac][0] == t), key=chr_sort_key)
            summed = sum(q_len.get(x, 0) for x in members)
            tl = t_len.get(t, 0)
            print(f"    {args.target_label} {t} = {args.query_label} "
                  f"{' + '.join(members)}")
            if tl:
                print(f"        {summed/1e6:.2f} Mb of query summed vs "
                      f"{tl/1e6:.2f} Mb target -> {100.0*summed/tl:.1f}% "
                      f"of its length reconstructed")

    print("")
    print("  Written:")
    print(f"    {summ_path}")
    print(f"    {corr_path}")
    print(f"    {ren_path}")
    print("")
    print("  These are correspondences, not polarised rearrangements: they say")
    print("  the two genomes differ, not which lineage changed. Deciding that")
    print("  needs an outgroup.")


if __name__ == "__main__":
    main()
