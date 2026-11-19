#!/usr/bin/env python3
"""Prepare the 46 haplotype assemblies (23 birds) for NCBI's diploid genome submission portal.

Follows NCBI's instructions for this project (email, Nov 2026):
  - assemblies arranged in pairs, one pair per bird
  - the Z, the W and the mitogenome go in the haplotype submitted as "Haplotype 1"
  - no sequence identifier may contain an accession from the RagTag reference (NC_/NW_/...)
  - chromosome and organelle assignments in the FASTA definition lines

Two modes
  manifest  write birds.tsv (one row per bird: species, sample, both haplotype FASTAs,
            mitogenome FASTA found under --mito-dir, sex from W content) and print the
            sequence-name patterns found, for checking before `prep`
  prep      read birds.tsv and write, per bird, <sample>_hap1.fsa.gz and <sample>_hap2.fsa.gz
            (submission labels), plus rename maps, a per-haplotype table and a report

Rules applied in `prep`
  - chr_<N> -> [location=chromosome] [chromosome=<N>]; several pieces of one chromosome all
    carry the same chromosome value
  - females: the source haplotype with more chr_Z sequence becomes Haplotype 1. Z and W pieces
    in the other haplotype are moved into Haplotype 1 unless >= DUP_COV of their length aligns
    (minimap2 -x asm5) to Haplotype 1's own Z/W, in which case they are treated as duplicated
    sequence and left in Haplotype 2 as unplaced. Without minimap2 nothing is moved.
  - males: hap1 is Haplotype 1; any chr_W piece (mis-assigned sequence) becomes unplaced
  - chr_MT / mitochondrial sequences in the assemblies are dropped; the curated mitogenome is
    added to Haplotype 1 as [location=mitochondrion] [completeness=complete] [topology=circular]
  - everything else is unplaced; terminal Ns are trimmed and sequences < 200 bp dropped
  - new IDs: <sample>_h<1|2>_<chrom>[_<k>] and <sample>_h<1|2>_u<0001>; the run stops if any
    output identifier or definition line contains an accession-like string
"""
import argparse, glob, gzip, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

ORGANISM = {"LEPC": "Tympanuchus pallidicinctus", "GRPC": "Tympanuchus cupido",
            "STGR": "Tympanuchus phasianellus"}
CHROM_RE = re.compile(r"^chr_([0-9]+|Z|W|MT)(?:_.*)?$")
ACC_RE = re.compile(r"\b[A-Z]{1,2}_?\d{5,}(\.\d+)?\b")     # NC_064463.1, NW_026294813.1, CM012345
MIN_LEN = 200
W_FEMALE_MB = 5.0                                          # >= this much chr_W => female


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fa(fh, ident, defline, seq, width=80):
    fh.write(f">{ident} {defline}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def classify(name):
    m = CHROM_RE.match(name)
    return m.group(1) if m else None


def aligned_fraction(queries, targets, threads):
    """Fraction of each query covered by minimap2 asm5 alignments to targets (dict name->seq)."""
    if not queries or not targets or not shutil.which("minimap2"):
        return None
    with tempfile.TemporaryDirectory() as td:
        q, t = os.path.join(td, "q.fa"), os.path.join(td, "t.fa")
        with open(q, "w") as f:
            for n, s in queries.items():
                f.write(f">{n}\n{s}\n")
        with open(t, "w") as f:
            for n, s in targets.items():
                f.write(f">{n}\n{s}\n")
        out = subprocess.run(["minimap2", "-x", "asm5", "-t", str(threads), t, q],
                             capture_output=True, text=True, check=True).stdout
    iv = defaultdict(list)
    for line in out.splitlines():
        f = line.split("\t")
        iv[f[0]].append((int(f[2]), int(f[3])))
    frac = {}
    for n, s in queries.items():
        cov, end = 0, -1
        for a, b in sorted(iv.get(n, [])):
            if b > end:
                cov += b - max(a, end)
                end = b
        frac[n] = cov / len(s) if s else 0.0
    return frac


# ----------------------------------------------------------------------------- manifest
def manifest(a):
    rows, pat = [], defaultdict(int)
    for f in sorted(glob.glob(os.path.join(a.assembly_dir, "*_hap1" + a.suffix))):
        base = os.path.basename(f)[:-len("_hap1" + a.suffix)]
        sp, sample = base.split("_", 1)
        h2 = os.path.join(a.assembly_dir, f"{base}_hap2{a.suffix}")
        mito = ""
        if a.mito_dir:
            hits = sorted(p for p in glob.glob(os.path.join(a.mito_dir, "**", f"*{sample}*"), recursive=True)
                          if re.search(r"\.(fa|fasta|fna)(\.gz)?$", p))
            mito = hits[0] if hits else ""
            if len(hits) > 1:
                print(f"NOTE {sample}: {len(hits)} candidate mitogenome files, using {hits[0]} (check birds.tsv)")
        wlen = 0
        for h in (f, h2):
            if os.path.exists(h):
                for n, s in read_fasta(h).items():
                    pat[re.sub(r"\d+", "#", n)] += 1
                    if classify(n) == "W":
                        wlen += len(s)
        sex = "female" if wlen >= W_FEMALE_MB * 1e6 else "male"
        rows.append([sp, sample, sex, f, h2 if os.path.exists(h2) else "MISSING", mito or "MISSING"])
        print(f"{sp} {sample}: {sex} (chr_W {wlen / 1e6:.1f} Mb), mito {'found' if mito else 'MISSING'}")
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "birds.tsv"), "w") as fh:
        fh.write("species\tsample\tsex\thap1_fasta\thap2_fasta\tmito_fasta\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    # sequence-name patterns across all assemblies, so unexpected names are seen before prep
    if rows:
        print("\nSequence-name patterns across all assemblies (# = digits); anything not chr_<N>/Z/W/MT "
              "becomes an unplaced sequence with a new name:")
        for p, k in sorted(pat.items(), key=lambda x: -x[1]):
            flag = "   <-- contains an accession; will be renamed" if ACC_RE.search(p.replace("#", "1" * 6)) else ""
            print(f"  {k:6d}  {p}{flag}")
    print(f"\nWrote {a.out}/birds.tsv ({len(rows)} birds). Check sex and mito_fasta, then run prep.")


# ----------------------------------------------------------------------------- prep
def prep(a):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(a.out, "birds.tsv"))][1:]
    files = os.path.join(a.out, "files"); maps = os.path.join(a.out, "rename_maps")
    os.makedirs(files, exist_ok=True); os.makedirs(maps, exist_ok=True)
    info = open(os.path.join(a.out, "haplotype_table.tsv"), "w")
    info.write("bird_umbrella_group\tsample_name\tspecies\torganism\tsex\tsubmission_haplotype\t"
               "source_haplotype\tassembly_name\tfilename\tlength_bp\tn_sequences\tchromosomes\t"
               "carries_Z\tcarries_W\tcarries_mito\n")
    report = open(os.path.join(a.out, "prep_report.txt"), "w")
    if not shutil.which("minimap2"):
        print("WARNING: minimap2 not on PATH: female Z/W pieces in the second haplotype will NOT be moved",
              file=sys.stderr)

    for sp, sample, sex, h1f, h2f, mitof in rows:
        if a.only and sample not in a.only.split(","):
            continue
        org = ORGANISM[sp]
        log = [f"== {sp} {sample} ({sex})"]
        src = {"hap1": read_fasta(h1f), "hap2": read_fasta(h2f)}
        zlen = {h: sum(len(s) for n, s in d.items() if classify(n) == "Z") for h, d in src.items()}
        first = "hap1" if sex == "male" else max(zlen, key=zlen.get)
        second = "hap2" if first == "hap1" else "hap1"
        log.append(f"   Haplotype 1 = source {first} (chr_Z {zlen[first] / 1e6:.1f} Mb); "
                   f"Haplotype 2 = source {second} (chr_Z {zlen[second] / 1e6:.1f} Mb)")

        # assignments: (label, name, seq, chrom or None, note)
        out = {1: [], 2: []}
        for lab, h in ((1, first), (2, second)):
            for n, s in src[h].items():
                c = classify(n)
                if c == "MT":
                    log.append(f"   dropped {h}:{n} ({len(s):,} bp) - mitochondrial; curated mitogenome used")
                    continue
                if c == "W" and sex == "male":
                    out[lab].append((n, s, None, "chr_W in a male: unplaced")); continue
                out[lab].append((n, s, c, ""))

        if sex == "female":
            cand = {n: s for n, s, c, _ in out[2] if c in ("Z", "W")}
            h1sex = {n: s for n, s, c, _ in out[1] if c in ("Z", "W")}
            frac = aligned_fraction(cand, h1sex, a.threads)
            keep2 = []
            for n, s, c, note in out[2]:
                if n not in cand:
                    keep2.append((n, s, c, note)); continue
                if frac is None:
                    keep2.append((n, s, None, f"chr_{c} piece in Haplotype 2 left unplaced (not checked: no minimap2)"))
                elif frac[n] >= a.dup_cov:
                    keep2.append((n, s, None, f"chr_{c} piece: {frac[n]:.0%} aligns to Haplotype 1 sex chromosomes; "
                                              f"treated as duplicated, left unplaced in Haplotype 2"))
                else:
                    out[1].append((n, s, c, f"moved from {second} ({frac[n]:.0%} aligns to Haplotype 1)"))
            out[2] = keep2

        for lab in (1, 2):
            tag = f"{sample}_h{lab}"
            fn = f"{sample}_hap{lab}.fsa.gz"
            n_u, nseq, total, used, chroms = 0, 0, 0, defaultdict(int), set()
            mp = open(os.path.join(maps, f"{sample}_hap{lab}.tsv"), "w")
            mp.write("new_id\tsource_haplotype\tsource_name\tlength_bp\tassignment\tnote\n")
            with gzip.open(os.path.join(files, fn), "wt", compresslevel=6) as fh:
                order = sorted(out[lab], key=lambda r: (r[2] is None, *(
                    (0, int(r[2])) if r[2] and r[2].isdigit() else (1, 0)), r[2] or "", -len(r[1])))
                for n, s, c, note in order:
                    s = s.strip("Nn")
                    if len(s) < MIN_LEN:
                        mp.write(f"-\t{first if lab == 1 else second}\t{n}\t{len(s)}\tdropped (< {MIN_LEN} bp)\t{note}\n")
                        continue
                    if c:
                        used[c] += 1
                        ident = f"{tag}_{c}" + (f"_{used[c]}" if used[c] > 1 else "")
                        defl = f"[organism={org}] [isolate={sample}] [location=chromosome] [chromosome={c}]"
                        chroms.add(c)
                    else:
                        n_u += 1
                        ident = f"{tag}_u{n_u:04d}"
                        defl = f"[organism={org}] [isolate={sample}]"
                    if ACC_RE.search(ident + " " + defl):
                        sys.exit(f"ERROR: accession-like string in output: {ident} {defl}")
                    write_fa(fh, ident, defl, s)
                    nseq += 1; total += len(s)
                    srch = (first if lab == 1 else second) if "moved from" not in note else second
                    mp.write(f"{ident}\t{srch}\t{n}\t{len(s)}\t{('chromosome ' + c) if c else 'unplaced'}\t{note}\n")
                if lab == 1:
                    mito = read_fasta(mitof) if mitof not in ("", "MISSING") else {}
                    if len(mito) != 1:
                        sys.exit(f"ERROR {sample}: mitogenome file must hold exactly one sequence ({mitof})")
                    (mn, ms), = mito.items()
                    if not 15000 <= len(ms) <= 18000:
                        sys.exit(f"ERROR {sample}: mitogenome length {len(ms)} outside 15-18 kb")
                    write_fa(fh, f"{tag}_MT", f"[organism={org}] [isolate={sample}] [location=mitochondrion] "
                                              f"[completeness=complete] [topology=circular]", ms.upper())
                    nseq += 1; total += len(ms)
                    mp.write(f"{tag}_MT\tmitogenome\t{mn}\t{len(ms)}\tmitochondrion\tfrom {mitof}\n")
            mp.close()
            cz, cw = "Z" in chroms, "W" in chroms
            info.write(f"{sample}\t{sample}\t{sp}\t{org}\t{sex}\tHaplotype {lab}\t{first if lab == 1 else second}\t"
                       f"{sp}_{sample}.hap{lab}\t{fn}\t{total}\t{nseq}\t{len(chroms)}\t{'yes' if cz else 'no'}\t"
                       f"{'yes' if cw else 'no'}\t{'yes' if lab == 1 else 'no'}\n")
            log.append(f"   Haplotype {lab}: {fn}  {total / 1e6:,.1f} Mb in {nseq} sequences; "
                       f"{len(chroms)} chromosomes{' incl. Z' if cz else ''}{' and W' if cw else ''}"
                       f"{'; + mitogenome' if lab == 1 else ''}; {n_u} unplaced")
        for lab in (1, 2):
            for n, s, c, note in out[lab]:
                if note:
                    log.append(f"   {n} ({len(s):,} bp): {note}")
        print("\n".join(log)); report.write("\n".join(log) + "\n")
    info.close(); report.close()
    print(f"\nWrote {files}/, {maps}/, {a.out}/haplotype_table.tsv and {a.out}/prep_report.txt")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["manifest", "prep"])
    ap.add_argument("--assembly-dir", default="")
    ap.add_argument("--suffix", default=".pseudo_chr.fasta")
    ap.add_argument("--mito-dir", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dup-cov", type=float, default=0.8)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--only", default="", help="comma-separated samples (testing)")
    a = ap.parse_args()
    manifest(a) if a.mode == "manifest" else prep(a)
