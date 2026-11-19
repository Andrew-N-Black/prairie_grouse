#!/usr/bin/env python3
"""Label every gap in the NCBI submission files with its linkage evidence and write
NCBI feature tables (.tbl, assembly_gap features) for table2asn.

Why: NCBI converts runs of N to gaps automatically only when every gap has the same
linkage evidence. These assemblies have two kinds of 100-N gap:
  - YaHS gaps inside Hi-C scaffolds                      -> proximity-ligation
  - RagTag joins ordering scaffolds against the ptarmigan -> align-xgenus
        (+ proximity-ligation when step 06's Hi-C test called the join "supported")
Both are exactly 100 N, so they are told apart by position.

How RagTag joins are found: for every gap row in each haplotype's RagTag AGP, the
40 bp on either side of the gap are taken from RagTag's own output FASTA
(ragtag.scaffold.fasta, AGP coordinates). Every 100-N run in a submission sequence
whose flanks match one of those signatures (either orientation) is that RagTag join;
this works even on chromosomes that were edited after RagTag (misjoin fixes), where
AGP coordinates no longer apply. Every other gap is a YaHS gap.

Inputs: <sub>/files/*.fsa.gz, <sub>/rename_maps/*.tsv, <sub>/haplotype_table.tsv
        <ragtag-dir>/<sample>.<hap>/{ragtag.scaffold.agp, ragtag.scaffold.fasta,
                                     <sample>.<hap>.rename_map.tsv}
        <joins-dir>/<SP>_<sample>_<hap>.joins.tsv   (step 06 verdicts)
Output: <out>/<sample>_hap<n>/<sample>_hap<n>.{fsa,tbl} and <out>/gap_summary.tsv
"""
import argparse, glob, gzip, os, re, sys
from collections import defaultdict

K = 40
GAP_MIN = 10
VERDICTS = {"supported", "unsupported", "ambiguous", "not_scored_small", "not_scored", "uncalibrated"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, head, buf = None, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, head, "".join(buf)
                head = line.rstrip("\n"); name = head[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, head, "".join(buf)


def read_agp_gaps(path):
    gaps = defaultdict(list)            # object -> [(beg, end)] 1-based
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[4] in ("N", "U"):
                gaps[f[0]].append((int(f[1]), int(f[2])))
    return gaps


def read_rmap(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and not line.startswith("#"):
                out[f[0]] = f[1]          # AGP object -> name given by step 05
    return out


def read_verdicts(path, name_to_obj, gaps):
    """(object, gap beg) -> verdict, using the column layout found in the table."""
    if not os.path.exists(path):
        return {}, f"no joins table {path}"
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    if not rows:
        return {}, "empty joins table"
    n = max(len(r) for r in rows)
    col = lambda pred: max(range(n), key=lambda c: sum(1 for r in rows if len(r) > c and pred(r[c])))
    vc = col(lambda v: v in VERDICTS)
    cc = col(lambda v: v.startswith("chr_"))
    best, pc = -1, None
    for c in range(n):
        if c in (vc, cc):
            continue
        hit = 0
        for r in rows:
            if len(r) > c and re.fullmatch(r"-?\d+", r[c] or "") and r[cc] in name_to_obj:
                p = int(r[c])
                hit += any(b - 1 <= p <= e + 1 for b, e in gaps.get(name_to_obj[r[cc]], []))
        if hit > best:
            best, pc = hit, c
    out = {}
    for r in rows:
        if len(r) <= max(vc, cc, pc) or r[cc] not in name_to_obj or not re.fullmatch(r"-?\d+", r[pc] or ""):
            continue
        obj, p = name_to_obj[r[cc]], int(r[pc])
        for b, e in gaps.get(obj, []):
            if b - 1 <= p <= e + 1:
                out[(obj, b)] = r[vc]
    return out, f"verdict col '{header[vc]}', chrom col '{header[cc]}', position col '{header[pc]}' ({best} rows on gaps)"


def signatures(rdir, sample, hap, joins_path):
    agp = os.path.join(rdir, "ragtag.scaffold.agp")
    sfa = os.path.join(rdir, "ragtag.scaffold.fasta")
    rmp = os.path.join(rdir, f"{sample}.{hap}.rename_map.tsv")
    for p in (agp, sfa):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing {p}")
    gaps = read_agp_gaps(agp)
    rmap = read_rmap(rmp) if os.path.exists(rmp) else {}
    name_to_obj = {v: k for k, v in rmap.items()}
    verdict, note = read_verdicts(joins_path, name_to_obj, gaps)
    sig = {}
    n_gap = 0
    for obj, _h, seq in fasta(sfa):
        for b, e in gaps.get(obj, []):
            n_gap += 1
            L, R = seq[max(0, b - 1 - K):b - 1].upper(), seq[e:e + K].upper()
            if len(L) < K or len(R) < K:
                continue
            v = verdict.get((obj, b), "not_in_joins_table")
            sig[(L, R)] = (obj, b, v)
            sig[(rc(R), rc(L))] = (obj, b, v)
    return sig, n_gap, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, help="ncbi_submission directory")
    ap.add_argument("--ragtag-dir", required=True)
    ap.add_argument("--joins-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    table = [l.rstrip("\n").split("\t") for l in open(os.path.join(a.sub, "haplotype_table.tsv"))]
    hdr = {c: i for i, c in enumerate(table[0])}
    summ = open(os.path.join(a.out, "gap_summary.tsv"), "w")
    summ.write("file\tsource\tgaps\tyahs_proximity_ligation\tragtag_hic_supported\tragtag_other\t"
               "ragtag_joins_in_agp\tragtag_joins_found\tverdict_lookup\n")
    cache = {}
    for row in table[1:]:
        sample, sp, fn = row[hdr["sample_name"]], row[hdr["species"]], row[hdr["filename"]]
        stem = fn.replace(".fsa.gz", "")
        if a.only and sample not in a.only.split(","):
            continue
        # source haplotype of every sequence, from the rename map
        srcmap = {}
        for l in open(os.path.join(a.sub, "rename_maps", stem + ".tsv")):
            f = l.rstrip("\n").split("\t")
            if f[0] not in ("new_id", "-"):
                srcmap[f[0]] = f[1]
        wd = os.path.join(a.out, stem); os.makedirs(wd, exist_ok=True)
        counts = defaultdict(int); found = defaultdict(set); notes = {}
        with open(os.path.join(wd, stem + ".fsa"), "w") as fo, open(os.path.join(wd, stem + ".tbl"), "w") as ft:
            for sid, head, seq in fasta(os.path.join(a.sub, "files", fn)):
                fo.write(head + "\n")
                for i in range(0, len(seq), 80):
                    fo.write(seq[i:i + 80] + "\n")
                src = srcmap.get(sid, "")
                runs = [(m.start(), m.end()) for m in re.finditer(r"[Nn]{%d,}" % GAP_MIN, seq)]
                if not runs:
                    continue
                if src not in ("hap1", "hap2"):
                    sys.exit(f"ERROR: {sid} has gaps but no source haplotype in the rename map")
                key = (sp, sample, src)
                if key not in cache:
                    cache[key] = signatures(os.path.join(a.ragtag_dir, f"{sample}.{src}"), sample, src,
                                            os.path.join(a.joins_dir, f"{sp}_{sample}_{src}.joins.tsv"))
                sig, n_agp, note = cache[key]
                notes[src] = (n_agp, note)
                ft.write(f">Feature {sid}\n")
                for s, e in runs:
                    L, R = seq[max(0, s - K):s].upper(), seq[e:e + K].upper()
                    hit = sig.get((L, R))
                    if hit:
                        found[src].add((hit[0], hit[1]))
                        ev = ["proximity-ligation", "align-xgenus"] if hit[2] == "supported" else ["align-xgenus"]
                        counts["ragtag_hic_supported" if hit[2] == "supported" else "ragtag_other"] += 1
                    else:
                        ev = ["proximity-ligation"]
                        counts["yahs"] += 1
                    counts["gaps"] += 1
                    ft.write(f"{s + 1}\t{e}\tassembly_gap\n\t\t\testimated_length\tunknown\n"
                             f"\t\t\tgap_type\twithin scaffold\n")
                    for x in ev:
                        ft.write(f"\t\t\tlinkage_evidence\t{x}\n")
        srcs = sorted(notes)
        summ.write(f"{stem}\t{','.join(srcs)}\t{counts['gaps']}\t{counts['yahs']}\t{counts['ragtag_hic_supported']}\t"
                   f"{counts['ragtag_other']}\t{sum(notes[s][0] for s in srcs)}\t{sum(len(found[s]) for s in srcs)}\t"
                   f"{' | '.join(notes[s][1] for s in srcs)}\n")
        print(f"{stem}: {counts['gaps']} gaps = {counts['yahs']} YaHS, {counts['ragtag_hic_supported']} RagTag Hi-C-supported, "
              f"{counts['ragtag_other']} RagTag other; RagTag joins found {sum(len(found[s]) for s in srcs)} "
              f"of {sum(notes[s][0] for s in srcs)} in the AGP")
    summ.close()


if __name__ == "__main__":
    main()
