#!/usr/bin/env python3
"""TE composition by category for every haplotype, separately for the autosomes, the Z and the W.
Reads the RepeatMasker .out files written by 12_pangenome_build.sh (STAGE=te; one merged library from the
three species' RepeatModeler2 libraries, so all haplotypes are annotated alike) and the final assemblies' .fai.
Masked bp are counted per category after merging overlapping hits within the category.
usage: te_categories.py <manifest> <final_dir> <te_dir> <out.tsv>
Output columns: species sample haplotype chrom_class class_bp category masked_bp"""
import csv, os, re, sys

manifest, final_dir, te_dir, out_tsv = sys.argv[1:5]
AUTO = re.compile(r'^chr_\d+(_p\d+)?$')
CLS = lambda c: 'A' if AUTO.match(c) else ('Z' if re.match(r'^chr_Z(_p\d+)?$', c) else
                                           ('W' if re.match(r'^chr_W(_p\d+)?$', c) else None))


def category(rc):
    """RepeatMasker class/family -> reporting category"""
    c, _, fam = rc.partition('/')
    if c == 'LINE':
        return 'LINE/CR1' if fam.startswith('CR1') else 'LINE/other'
    if c == 'LTR':
        for k in ('ERV1', 'ERVK', 'ERVL'):
            if fam.startswith(k):
                return 'LTR/' + k
        return 'LTR/other'
    if c in ('SINE', 'DNA', 'RC', 'Satellite'):
        return c
    if c in ('Simple_repeat', 'Low_complexity'):
        return 'Simple/low complexity'
    if c in ('rRNA', 'tRNA', 'snRNA', 'scRNA', 'srpRNA', 'ARTEFACT'):
        return None
    return 'Unknown'


def merged_len(xs):
    xs.sort(); tot = 0; cs, ce = xs[0]
    for a, b in xs[1:]:
        if a > ce + 1:
            tot += ce - cs + 1; cs, ce = a, b
        else:
            ce = max(ce, b)
    return tot + ce - cs + 1


with open(manifest) as fh:
    rows = list(csv.DictReader([l for l in fh if not l.startswith('#')], delimiter='\t'))
out = open(out_tsv, 'w')
out.write('species\tsample\thaplotype\tchrom_class\tclass_bp\tcategory\tmasked_bp\n')
for r in rows:
    sp, s = r['species'].strip(), r['sample_id'].strip()
    if not s:
        continue
    for h in ('hap1', 'hap2'):
        stem = f'{sp}_{s}_{h}.pseudo_chr'
        fai = os.path.join(final_dir, stem + '.fasta.fai')
        rmout = os.path.join(te_dir, stem, stem + '.fasta.out')
        if not (os.path.exists(fai) and os.path.exists(rmout)):
            sys.stderr.write(f'missing {fai} or {rmout}\n'); continue
        cls = {}; size = {}
        for line in open(fai):
            f = line.split('\t'); k = CLS(f[0])
            if k:
                cls[f[0]] = k; size[k] = size.get(k, 0) + int(f[1])
        iv = {}
        for line in open(rmout):
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or f[4] not in cls:
                continue
            cat = category(f[10])
            if cat:
                iv.setdefault((cls[f[4]], cat, f[4]), []).append((int(f[5]), int(f[6])))
        tot = {}
        for (k, cat, _), xs in iv.items():
            tot[(k, cat)] = tot.get((k, cat), 0) + merged_len(xs)
        for (k, cat), bp in sorted(tot.items()):
            out.write(f'{sp}\t{s}\t{h}\t{k}\t{size[k]}\t{cat}\t{bp}\n')
        sys.stderr.write(f'{stem}: ' + ', '.join(f'{k} {size[k] / 1e6:.1f} Mb' for k in sorted(size)) + '\n')
out.close()
