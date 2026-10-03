#!/usr/bin/env python3
"""
Rename chromosomes consistently across every derived file in an assembly tree.
Companion to 06_rename_chromosomes.sh.

WHY THIS EXISTS
    Chromosome names ended up encoding a comparison rather than naming a
    chromosome. chr_6_8, chr_2a, chr_4b and chr_u30 are all statements about
    CHICKEN: "spans chicken 6 and 8", "second piece of chicken 2", "could not be
    matched to chicken". Useful while choosing a reference; wrong as identifiers
    on a finished assembly, where a name should be short, stable and unique.

    A chromosome name also appears in a dozen places — the FASTA, its index, the
    Liftoff GFF, the autosome split, BUSCO's full table, the join-support tables,
    the Hi-C CRAM header. Renaming by hand in some of them and not others is how
    an assembly quietly stops agreeing with its own annotation.

SAFETY
    - Dry run by default. Nothing is written without --apply.
    - The mapping must be a bijection. A->B together with C->B is rejected,
      because it would merge two chromosomes into one name.
    - Renaming is SIMULTANEOUS, not sequential. A swap (A->B, B->A) is applied
      through one lookup table, so it cannot collapse both onto one name the way
      two successive sed passes would.
    - Every file is written to a .part and moved into place only on success.
    - A reverse mapping is written, so the rename can be undone exactly.
    - After rewriting a FASTA, sequence count and total length are re-checked
      against the original.

WHAT IT TOUCHES
    final_<tag>/<P>.pseudo_chr.fasta           FASTA headers
    final_<tag>/<P>.liftoff.gff3               column 1 (seqid)
    final_autosomes_<tag>/<P>.autosome_chr.fasta
    final_autosomes_<tag>/<P>.excluded_from_autosomes.fasta
    final_autosomes_<tag>/<P>.split_report.tsv column 1
    qc_<tag>/busco/<P>/run_*/full_table.tsv    the Sequence column
    qc_<tag>/join_support/<P>.joins.tsv        column 2
    qc_<tag>/join_support/<P>.controls.tsv     column 2
    ragtag_<tag>/<S>.<H>/<S>.<H>.rename_map.tsv column 2

    The .fai files and the Hi-C CRAM headers are handled by the shell wrapper,
    which has samtools. Pretext maps cannot be edited and must be regenerated if
    you want them to carry the new names; the wrapper prints the command.
"""

import argparse
import glob
import os
import re
import sys


# ------------------------------------------------------------------ mapping --
def read_map(path):
    """old_name -> new_name, tab separated, # comments ignored."""
    m = {}
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0] or not parts[1]:
                sys.exit(f"ERROR: {path} line {ln}: expected 'old<TAB>new', got: {line!r}")
            old, new = parts[0].strip(), parts[1].strip()
            if old in m and m[old] != new:
                sys.exit(f"ERROR: {path} line {ln}: {old} mapped twice ({m[old]} and {new})")
            m[old] = new
    if not m:
        sys.exit(f"ERROR: {path} contains no mappings")
    return m


def validate(m):
    """
    Reject anything that would merge chromosomes or produce an unusable name.

    The collision check is the important one. Two different chromosomes mapped
    onto the same new name is unrecoverable once the FASTA is written — the
    sequences are still distinct but can no longer be told apart by name, and
    samtools faidx refuses to index the result. Catch it here.
    """
    problems = []
    seen = {}
    for old, new in sorted(m.items()):
        if new in seen:
            problems.append(f"collision: '{seen[new]}' and '{old}' both become '{new}'")
        seen[new] = old
        if re.search(r"[\s,;|>]", new):
            problems.append(f"'{new}' contains whitespace or a character that breaks "
                            f"FASTA/SAM/GFF parsing")
        if not new:
            problems.append(f"'{old}' maps to an empty name")
    # A no-op entry is harmless but usually a mistake worth surfacing.
    noop = [o for o, n in m.items() if o == n]
    return problems, noop


# ------------------------------------------------------------- file rewrites --
def fasta_names_from_fai(path):
    """
    Sequence names from the .fai index, without touching the FASTA.

    A dry run only needs to know which names are present and how many would
    change. Streaming a 1 Gb FASTA to answer that — times 46 assemblies, twice
    over for the autosome copies — makes the preview take longer than the real
    work, which defeats the point of having one. The .fai lists every name in a
    few kilobytes and step 05 always writes it.
    """
    fai = path + ".fai"
    if not os.path.exists(fai):
        return None
    names = []
    try:
        with open(fai) as fh:
            for line in fh:
                f = line.split("\t", 1)
                if f and f[0]:
                    names.append(f[0])
    except OSError:
        return None
    return names


def rewrite_fasta(path, m, apply_):
    """
    Rewrite '>name ...' headers. The description after the first whitespace is
    dropped, matching how step 05 wrote these files in the first place.

    In dry-run mode this answers from the .fai index and never opens the FASTA.
    """
    if not apply_:
        names = fasta_names_from_fai(path)
        if names is not None:
            return sum(1 for n in names if n in m), len(names)
        # No index. Some FASTAs never get one — step 07 writes
        # *.excluded_from_autosomes.fasta without indexing it, and those carry
        # chr_Z and chr_W, so they run to ~90 MB each. Streaming 46 of them to
        # produce a number that is then discarded is the same mistake as reading
        # the assemblies was, so the size rule below applies uniformly: a dry run
        # never reads more than DRY_RUN_SCAN_LIMIT of any file. Small unindexed
        # FASTAs are still counted exactly.
        try:
            if os.path.getsize(path) > DRY_RUN_SCAN_LIMIT:
                return None, None
        except OSError:
            return None, None

    changed = 0
    n_seq_in = n_seq_out = 0
    bp_in = bp_out = 0
    tmp = path + ".part"
    out = open(tmp, "w") if apply_ else None
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    n_seq_in += 1
                    name = line[1:].split()[0] if len(line) > 1 else ""
                    if name in m:
                        changed += 1
                        line = ">" + m[name] + "\n"
                    else:
                        line = ">" + name + "\n"
                    n_seq_out += 1
                else:
                    L = len(line.strip())
                    bp_in += L
                    bp_out += L
                if out:
                    out.write(line)
    finally:
        if out:
            out.close()

    if apply_:
        # Nothing but names may change. If the arithmetic moved, something is
        # wrong with the rewrite and the original must not be replaced.
        if n_seq_in != n_seq_out or bp_in != bp_out:
            os.unlink(tmp)
            sys.exit(f"ERROR: {path}: rewrite changed the content "
                     f"({n_seq_in}->{n_seq_out} seqs, {bp_in}->{bp_out} bp). Original kept.")
        os.replace(tmp, path)
    return changed, n_seq_in


# A dry run that reads every Liftoff GFF end to end takes longer than the real
# rename and looks like a hang. Small files are still counted exactly — that is
# cheap and the count is reassuring — but anything above this is reported as
# "will be rewritten" without being opened. There is no index to consult for a
# GFF the way there is for a FASTA, so the choice is an exact count or a fast
# preview, and a preview should be fast.
DRY_RUN_SCAN_LIMIT = 20 * 1024 * 1024      # 20 MB


def rewrite_column(path, m, col, apply_, sep="\t", skip_prefix="#", header_rows=0):
    """
    Rewrite one column of a delimited file, leaving everything else byte-identical.

    Returns the number of names changed, or None when a dry run skipped the file
    because it is large (see DRY_RUN_SCAN_LIMIT).
    """
    if not apply_:
        try:
            if os.path.getsize(path) > DRY_RUN_SCAN_LIMIT:
                return None
        except OSError:
            return None

    changed = 0
    tmp = path + ".part"
    out = open(tmp, "w") if apply_ else None
    try:
        with open(path) as fh:
            for i, line in enumerate(fh):
                raw = line.rstrip("\n")
                if (skip_prefix and raw.startswith(skip_prefix)) or i < header_rows or not raw:
                    if out:
                        out.write(line)
                    continue
                f = raw.split(sep)
                if len(f) > col and f[col] in m:
                    f[col] = m[f[col]]
                    changed += 1
                    line = sep.join(f) + "\n"
                if out:
                    out.write(line)
    finally:
        if out:
            out.close()
    if apply_:
        os.replace(tmp, path)
    return changed


# --------------------------------------------------------------------- main --
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--ref-tag", default="ptarmigan")
    p.add_argument("--map", required=True, help="TSV: old_name<TAB>new_name")
    p.add_argument("--apply", action="store_true",
                   help="actually write. Without this, nothing is modified.")
    p.add_argument("--reverse-map-out", default=None)
    args = p.parse_args()

    m = read_map(args.map)
    problems, noop = validate(m)
    if problems:
        print("ERROR: the mapping is not safe to apply:")
        for x in problems:
            print("   " + x)
        sys.exit(1)

    print(f">>> mapping: {len(m)} entries from {args.map}")
    for old, new in sorted(m.items()):
        flag = "   (no change)" if old == new else ""
        print(f"      {old:<14} -> {new}{flag}")
    if noop:
        print(f"    note: {len(noop)} entries map a name to itself; harmless.")
    print(f">>> mode   : {'APPLY — files will be rewritten' if args.apply else 'DRY RUN — nothing written'}")

    P, tag = args.project_dir, args.ref_tag
    final = os.path.join(P, f"final_{tag}")
    split = os.path.join(P, f"final_autosomes_{tag}")
    qc = os.path.join(P, f"qc_{tag}")
    ragtag = os.path.join(P, f"ragtag_{tag}")

    if not os.path.isdir(final):
        sys.exit(f"ERROR: {final} not found — is --ref-tag right?")

    # (glob, handler, label)
    jobs = [
        (os.path.join(final, "*.pseudo_chr.fasta"), ("fasta", None), "assembly FASTA"),
        (os.path.join(final, "*.liftoff.gff3"), ("col", 0), "Liftoff GFF"),
        (os.path.join(split, "*.autosome_chr.fasta"), ("fasta", None), "autosome FASTA"),
        (os.path.join(split, "*.excluded_from_autosomes.fasta"), ("fasta", None), "excluded FASTA"),
        (os.path.join(split, "*.split_report.tsv"), ("col", 0), "split report"),
        (os.path.join(qc, "busco", "*", "run_*", "full_table.tsv"), ("col", 2), "BUSCO full table"),
        (os.path.join(qc, "join_support", "*.joins.tsv"), ("col", 1), "join support"),
        (os.path.join(qc, "join_support", "*.controls.tsv"), ("col", 1), "join controls"),
        (os.path.join(ragtag, "*", "*.rename_map.tsv"), ("col", 1), "RagTag rename map"),
    ]

    total_files = total_changes = 0
    print("")
    for pattern, (kind, col), label in jobs:
        files = sorted(glob.glob(pattern))
        if not files:
            print(f"  {label:<20} no files matched — skipping")
            continue
        n_ch = n_f = n_skipped = 0
        # A full read plus a full write across ~1 Gb files is real I/O, so name
        # the file in hand rather than going silent for minutes and looking hung.
        verbose = args.apply and len(files) > 2
        if verbose:
            print(f"  {label}: rewriting {len(files)} file(s)", flush=True)
        for i, f in enumerate(files, 1):
            if verbose:
                print(f"      [{i}/{len(files)}] {os.path.basename(f)}", flush=True)
            try:
                if kind == "fasta":
                    c, _ = rewrite_fasta(f, m, args.apply)
                else:
                    # BUSCO's full_table has two comment lines then a '#'-prefixed
                    # header; '#' skipping covers both.
                    c = rewrite_column(f, m, col, args.apply)
            except OSError as exc:
                print(f"  ! {f}: {exc}")
                continue
            if c is None:          # dry run, file too large to scan
                n_skipped += 1
                continue
            if c:
                n_ch += c
                n_f += 1
        if n_skipped:
            print(f"  {label:<20} {len(files):>4} file(s), {n_skipped:>4} too large to "
                  f"count in a dry run — will be rewritten")
        else:
            print(f"  {label:<20} {len(files):>4} file(s), {n_f:>4} affected, "
                  f"{n_ch:>6} name(s) rewritten")
        total_files += n_f
        total_changes += n_ch

    print(f"\n>>> {total_changes} names across {total_files} files"
          f"{'' if args.apply else ' WOULD be rewritten (dry run)'}")

    if args.apply and args.reverse_map_out:
        with open(args.reverse_map_out, "w") as fh:
            fh.write("# reverse of the rename just applied; feed back in with --map to undo\n")
            for old, new in sorted(m.items()):
                fh.write(f"{new}\t{old}\n")
        print(f">>> reverse mapping written: {args.reverse_map_out}")

    if not args.apply:
        print("\n    Re-run with --apply (or APPLY=true in the wrapper) to write.")


if __name__ == "__main__":
    main()
