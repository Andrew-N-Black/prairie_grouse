# Prairie Grouse Genome Assembly and QC

SLURM/bash pipeline (Purdue RCAC, Gautschi) for **haplotype-resolved reference
genome assembly and quality control** of three prairie grouse species:

| Species | Code | Scientific name |
|---|---|---|
| Lesser Prairie-Chicken | LEPC | *Tympanuchus pallidicinctus* |
| Greater Prairie-Chicken | GRPC | *T. cupido* |
| Sharp-tailed Grouse | STGR | *T. phasianellus* |

**n = 23 individuals × 2 haplotypes = 46 assemblies** (7 LEPC, 6 GRPC,
10 STGR), all assembled from **PacBio HiFi + Hi-C only**. Step 02 still supports ONT
ultra-long reads (`--ul`, for any manifest row whose `ont_ul` is a path rather than `NA`),
but none are used: three birds were sequenced with ONT-UL; F5468's yield was too low,
and F5545 (LEPC) and F5597 (GRPC), first assembled with `--ul`, had Merqury QV 51.7-55.9
against 62.2-64.7 for HiFi-only haplotypes, with no gain in contiguity, so they were
re-assembled without it ([`accessory_scripts/rerun_without_ONT.sh`](accessory_scripts/rerun_without_ONT.sh);
the ONT-based outputs are kept in `grouse_asm/archive_with_ONT/`).

Assemblies are ordered against the **rock ptarmigan** reference
(*Lagopus muta*, bLagMut1, `GCF_023343835.1`), and annotations are transferred
from it.

This repository is the code for **Objective 3** of the USFWS report *Grouse
genomics, October 2026* ("Haplotype-resolved reference genomes"). SLURM
account: `dewoody` (Gautschi; `fnrdewoody` is the Negishi account). Sequence data: NCBI BioProject PRJNA1513026
(embargoed).

| Report item | Produced by |
|---|---|
| Assembly, scaffolding, QA/QC methods | steps 02, 05, 05b, 11 |
| Figure S5 (contiguity: scaffolds, length, N50) | QUAST from step 05 -> `accessory_scripts/plot_assembly_contiguity.R` |
| Figure S6 (BUSCO completeness, 46 assemblies + ptarmigan) | steps 07, 08, 09 |
| Figures 13-14 (BUSCO karyotype and synteny; LEPC F5540, GRPC F5595, STGR F5457 hap1 + ptarmigan) | step 09 |
| Table S3 (sex from Z/W content; with CHD1 PCR, Figure S7) | step 10 |
| Table S4 (assembly chromosome names vs chicken GRCg7b karyotype) | steps 01 + 04 |
| Unsupported RagTag joins (chr_2) | step 11 |
| Mitogenomes and their QC (Appendix S1) | steps 12, 12b, 12c |
| Figure 17C (mitogenome tree) | steps 12d (12e draws it on its own) |
| Table 3, Figure 16 (genome-wide window phylogenies) | step 13 (+ 13b statistics, 13c/13d window recovery); runs on Negishi |
| Figure 17 (autosomal, Z and mitogenome summary trees) | steps 13e (Negishi), 13f (local R) |

The superseded **chicken-guided** version of this pipeline is archived,
unchanged, in [`chicken_reference_chain/`](chicken_reference_chain) for
reference only; nothing in it is part of the current workflow or the report.
Chicken GRCg7b is still downloaded (step 01), but only so step 04 can relate
ptarmigan chromosome numbers to the chicken karyotype.

---

## Why ptarmigan and not chicken

Chicken (GRCg7b) is the conventional reference for galliform work, and the
first version of this pipeline used it. It was replaced because
reference-guided scaffolding does not merely *suggest* a chromosome structure,
it *imposes* one — and wherever the reference and the sample karyotypes differ,
the assembly inherits the reference's version, silently.

Three karyotype differences between chicken and these birds were found by
aligning the two references to each other (step 04) and then testing each
proposed join against the grouse Hi-C data (step 11):

- **Chicken chromosomes 6 and 8 are a single chromosome** in ptarmigan and in
  all three *Tympanuchus* species (ptarmigan chr_5; 65.82 Mb in grouse (median) against
  65.80 + ptarmigan 65.88, i.e. 99.9% of the length accounted for).
- **The short arm of chicken chromosome 4 (GGA4p) is a free-standing
  chromosome** here (ptarmigan chr_13, 19.05 Mb). Chicken carries a fusion that
  these birds do not. The chicken-ordered assemblies were forced to reproduce
  that fusion in 28 of 46 haplotypes, and the Hi-C contradicted it in every one.
- **Chicken chromosome 2 corresponds to two chromosomes** (ptarmigan chr_3 and
  chr_7; 98.0% accounted for), confirmed independently by the AGP, the contact
  maps, and an insulation-score measurement.

Ptarmigan is also simply the closer relative — *Lagopus* and *Tympanuchus* are
both Tetraoninae, *Gallus* is Phasianinae — which shows up directly in
annotation transfer: **a mean of 477 unmapped Liftoff features per haplotype against
ptarmigan (about 21,900 in total) versus 95,451 against chicken**, a 4.3× reduction.

### What chromosome numbers mean here

bLagMut1 is a Sanger/VGP curated assembly and numbers its chromosomes
`SUPER_1, SUPER_2, …` **in descending length order**. That is a size rank, not
a homology statement: **ptarmigan chr_6 is not chicken chromosome 6.**

The assemblies carry the ptarmigan numbering, because a name should identify a
chromosome rather than assert a claim about a third genome. The correspondence
to the chicken karyotype — which is what the galliform literature, gene names
and cytogenetic results all use — is computed by **step 04** and is the thing
to cite when relating these assemblies to published work. Never translate a
chromosome number by assuming the numbers match.

---

## Pipeline

Run in order. Steps 01 and 03 are one-off setup; 02, 05, 07, 11 and 12 are array
jobs over the 23 samples; the rest are single jobs. Steps 12-13 use the finished
assemblies: 12 builds the mitogenomes and their tree, 13 builds nuclear window
phylogenies for Tree House Explorer (THEx) and runs on **Negishi**.

| Step | Script | What it does | Shape |
|---|---|---|---|
| 01 | [`01_download_chicken_reference.sh`](01_download_chicken_reference.sh) | Chicken GRCg7b FASTA + accession→chromosome map. Needed **only** by step 04, to tie ptarmigan chromosome numbers to the chicken karyotype. No annotation is downloaded — none is transferred from chicken. | single, 4 h |
| 02 | [`02_assembly_scaffold_array.sh`](02_assembly_scaffold_array.sh) | **The reference-free half.** HiFiAdapterFilt → hifiasm (Hi-C-phased, `--ul` where ONT exists) → bwa mem Hi-C alignment → yahs scaffolding → Juicebox `.hic`/`.assembly` for manual curation. QUAST after hifiasm and after yahs. Also builds every conda env the later steps use. | array ×23, 48 cpu, 160 G, 10 d |
| 03 | [`03_ptarmigan_reference.sh`](03_ptarmigan_reference.sh) | Rock ptarmigan FASTA + GFF3 + accession→chromosome map, resolving the NCBI directory name at runtime. Prints the karyotype it found and flags size-rank naming. | single, 2 h |
| 04 | [`04_chromosome_homology.sh`](04_chromosome_homology.sh) + [`.py`](04_chromosome_homology.py) | Aligns ptarmigan to chicken (`minimap2 -x asm20`) and works out which chicken chromosome each ptarmigan chromosome corresponds to. **This is where the karyotype findings above come from.** Writes the correspondence table, the summary, and a rename map. | single, 24 cpu, 8 h |
| 05 | [`05_ragtag_liftoff_array.sh`](05_ragtag_liftoff_array.sh) | **The reference-guided half.** RagTag ordering → chromosome naming → QUAST → Liftoff → orientation dot plot → Hi-C realignment + Pretext contact map → tidk → HiFi depth check. Reads step 02's yahs scaffolds; never writes to them. | array ×23, 24 cpu, 96 G, 3 d |
| 05b | [`05b_relabel_male_W.sh`](05b_relabel_male_W.sh) | Males (ZZ) have no W: anything RagTag placed on ptarmigan chr_W in a male haplotype is renamed `scaffold_N` (>= 50 kb) or moved to `unplaced_short` (< 50 kb), with the GFF3 updated and a log written. The six ZW females are never touched. | login node, minutes |
| 06 | [`06_rename_chromosomes.sh`](06_rename_chromosomes.sh) + [`.py`](06_rename_chromosomes.py) | Renames chromosomes consistently across every derived file (FASTA, GFF3, AGP, CRAM, BUSCO tables). **Only needed if step 05 ran with `CHR_NAMING=homology`** — its default produces the final names directly. | single, 4 h |
| 07 | [`07_busco_array.sh`](07_busco_array.sh) | Splits chr_Z / chr_W / chr_MT into their own FASTAs, then BUSCO v5.4.7 genome mode (metaeuk) against `aves_odb10` (8,338 orthologs) on the full assembly. | array ×23, 24 cpu, 2 d |
| 08 | [`08_ptarmigan_reference_busco.sh`](08_ptarmigan_reference_busco.sh) | Runs the ptarmigan reference through the identical BUSCO treatment, so it can be drawn alongside the grouse and so every ortholog gets a reference chromosome. Required by step 10's third check. | single, 24 cpu, 1 d |
| 09 | [`09_busco_plots.sh`](09_busco_plots.sh) + [`.py`](09_busco_plots.py) | BUSCO-Plot-Py figures: per-species completeness barplots, karyoplots, and pairwise horizontal synteny plots. | single, 8 cpu, 4 h |
| 10 | [`10_sex_chromosome_check.sh`](10_sex_chromosome_check.sh) + [`.py`](10_sex_chromosome_check.py) | Tests whether the short, low-BUSCO haplotypes are short because they lack the Z. Three independent checks: chromosome content, BUSCOs on chr_Z, and where the *missing* orthologs sit on the reference. | single, 30 min |
| 11 | [`11_join_support.sh`](11_join_support.sh) + [`.py`](11_join_support.py) | Scores **every join RagTag made** against the Hi-C: observed cross-join contacts versus what contiguous sequence would carry at the same separation, calibrated per library on pseudo-joins inside intact scaffolds. This is what catches a reference imposing a join the data does not support. | array ×23, 8 cpu, 8 h |
| 12 | [`12_mitohifi_array.sh`](12_mitohifi_array.sh) | One circularised, rotated, annotated mitogenome per bird from the HiFi reads (MitoHiFi v3, NUMT-bearing long reads filtered out), replacing the 20-50 kb hifiasm chr_MT concatemers. `setup` fetches the reference mitogenome. | array ×23, 32 cpu, 1 d |
| 12b | [`12b_mito_from_hic.sh`](12b_mito_from_hic.sh) | Mitogenome as a Hi-C read consensus for birds whose HiFi reads are dominated by NUMTs (STGR F5485, F5502). | array (2 rows), 16 cpu, 4 h |
| 12c | [`12c_mito_qc.sh`](12c_mito_qc.sh) | QC table (length, gene content, rotation, identity/coverage vs reference, read depth), `trim`/`rescue` fixes, pairwise identities, and `collect` of the PASS mitogenomes. | login node |
| 12d | [`12d_mito_tree.sh`](12d_mito_tree.sh) | MAFFT FFT-NS-i alignment of the 23 mitogenomes + GRPC reference + rock ptarmigan outgroup, IQ-TREE 3 (ModelFinder, 1,000 UFBoot), rooted on ptarmigan. | login node / short job |
| 12e | [`12e_plot_mito_tree_ggtree.R`](12e_plot_mito_tree_ggtree.R) | Stand-alone mitogenome tree from 12d (ggtree); in the report it is panel C of Figure 17 (13f). | local R |
| 13 | [`13_thex_phylogeny.sh`](13_thex_phylogeny.sh) + [`.py`](13_thex_phylogeny.py) | **Negishi.** Copies the assemblies from Gautschi, projects every haplotype onto the ptarmigan reference (minimap2), builds one alignment per chromosome, then THExBuilder 100-kb windows → trimAl → pairwise filter → IQ-TREE 2 per window, and writes the THEx Tree Viewer file (trees rooted on ptarmigan; `TopologyID` = species-level relationship), gene track and chromosome-length BED. | array → node → array → small job |
| 13b | [`13b_thex_stats.py`](13b_thex_stats.py) | Statistics for Table 3 / Figure 16 with chromosomes as the unit of replication: Z vs autosome rank, Spearman vs chromosome length and gene density, window-level logistic GEE clustered by chromosome. | local, THEx env |
| 13c | [`13c_window_diagnostics.py`](13c_window_diagnostics.py) | Why a region has no trees: per-window reference N/repeat/GC, haplotype coverage, genes and filter fate. | login node |
| 13d | [`13d_thex_rescue.sh`](13d_thex_rescue.sh) + [`.py`](13d_thex_rescue.py) | Re-runs windows lost because one to four haplotypes lack that stretch, without those haplotypes, then reruns step 13's IQ-TREE and summary stages. | node → array → small job |
| 13e | [`13e_summary_trees.sh`](13e_summary_trees.sh) | One summary tree for the autosomes and one for the Z from the step-13 window trees: wASTRAL (weighted by UFBoot and branch length) at haplotype and species level, plus gene concordance factors (IQ-TREE `--gcf`, % of window trees containing each branch). | `setup` once, then one standby job |
| 13f | [`13f_plot_summary_trees.R`](13f_plot_summary_trees.R) | Figure 17: autosomal and Z summary trees (gCF on nodes) beside the 12d mitogenome tree (ggtree + patchwork). | local R |
| 14 | [`14_z_markers.sh`](14_z_markers.sh) + [`.py`](14_z_markers.py) | **Species-diagnostic Z markers.** Finds Z sites that separate species in the 46 assemblies (aligned to the LEPC reference Z; minimap2 + paftools call), checks each ~300-bp window (single copy, no W gametolog match, repeats), designs primers (Primer3), then confirms allele frequencies in the 506-bird panel by species x context (ANGSD, groups in [`panel_groups.tsv`](panel_groups.tsv)). | array -> node -> array -> small job |

### Submitting

```bash
N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)   # 23

bash   01_download_chicken_reference.sh
bash   03_ptarmigan_reference.sh
sbatch 04_chromosome_homology.sh

sbatch --array=0-$((N-1))%6 02_assembly_scaffold_array.sh
# ... then, once every sample has finished:
sbatch --array=0-$((N-1))%6 05_ragtag_liftoff_array.sh

bash   05b_relabel_male_W.sh                    # after 05 (and 06, if used)
sbatch --array=0-$((N-1))%8 07_busco_array.sh
sbatch 08_ptarmigan_reference_busco.sh
sbatch 09_busco_plots.sh                       # after 07 and 08
sbatch 10_sex_chromosome_check.sh              # after 07 and 08
sbatch --array=0-$((N-1))%8 11_join_support.sh  # needs 05 with DO_HIC=true

# mitogenomes (Gautschi)
bash   12_mitohifi_array.sh setup               # once, login node
sbatch --array=0-$((N-1))%8 12_mitohifi_array.sh
sbatch -q standby -t 04:00:00 --array=7,8 12b_mito_from_hic.sh   # STGR F5485, F5502
bash   12c_mito_qc.sh qc && bash 12c_mito_qc.sh collect
bash   12d_mito_tree.sh                           # mitogenome tree; drawn as Figure 17C by 13f

# window phylogenies (Negishi, from a clone of this repo)
bash 13_thex_phylogeny.sh transfer               # assemblies + ptarmigan from Gautschi
bash 13_thex_phylogeny.sh setup                  # THEx conda env on scratch
bash 13_thex_phylogeny.sh submit                 # project -> prep -> iqtree -> summarize
bash 13d_thex_rescue.sh submit                   # optional: recover windows missing 1-4 haplotypes
bash 13e_summary_trees.sh setup && bash 13e_summary_trees.sh submit   # autosome + Z summary trees
Rscript 13f_plot_summary_trees.R autosomes.gcf.cf.tree chrZ.gcf.cf.tree mito_tree.contree   # local, Figure 17
python 13b_thex_stats.py grouse_TreeViewer_input.xlsx chromosome_lengths.bed grouse_thex_stats

# species-diagnostic Z markers (Gautschi)
bash 14_z_markers.sh setup && bash 14_z_markers.sh submit
```

Steps 01, 03 and 04 are independent of 02 and can run while the assemblies are
building. 02 is the long pole — budget ten days per sample.

### Stage toggles

Step 05's stages are individually switchable without editing the file, which is
the convenient way to do a fast structural pass first and add the expensive
stages later:

```bash
sbatch --export=ALL,DO_HIC=false,DO_DEPTH=false --array=0-22%6 05_ragtag_liftoff_array.sh
```

`DO_LIFTOFF`, `DO_QUAST`, `DO_HIC`, `DO_DOTPLOT`, `DO_TIDK`, `DO_DEPTH` and
`CHR_NAMING` all work this way. `DO_HIC` and `DO_DEPTH` are each a full pass
over a sequencing library and together account for most of the walltime;
`DO_HIC=true` is required before step 11.

---

## Output layout

Everything lives under `$CLUSTER_SCRATCH/GROUSE/grouse_asm`:

```
ref/                          both references, their indexes and chromosome maps
conda_envs/                   every tool env, built once by 02
tools/BUSCO-Plot-Py/          plotting library checkout                (09)
hifi_filtered/<sample>/       adapter-filtered HiFi reads               (02)
hifiasm/                      haplotype contigs                        (02)
yahs/                         Hi-C scaffolds + Juicebox .hic/.assembly (02)
qc/                           QUAST at contig and scaffold stage       (02)
qc/chromosome_homology/       ptarmigan↔chicken correspondence         (04)
ragtag_ptarmigan/             RagTag AGP + per-assembly rename map     (05)
liftoff_ptarmigan/            unmapped-feature lists                   (05)
final_ptarmigan/              *.pseudo_chr.fasta, *.liftoff.gff3       (05)
final_autosomes_ptarmigan/    autosome-only and sex/organelle FASTAs   (07)
qc_ptarmigan/quast/           QUAST on the final assemblies            (05)
qc_ptarmigan/depth_check/     HiFi coverage per chromosome             (05)
qc_ptarmigan/busco/           BUSCO runs, 46 assemblies + reference    (07, 08)
qc_ptarmigan/plots/           barplots, karyoplots, synteny            (09)
qc_ptarmigan/sex_chromosome_check/                                     (10)
qc_ptarmigan/join_support/    per-join Hi-C support scores             (11)
mitogenomes/<SPECIES>_<SAMPLE>/ final_mitogenome.{fasta,gb,gff}        (12, 12b)
mitogenomes/mito_qc.tsv, mitogenomes_pass.fasta                        (12c)
mitogenomes/tree/             alignment and IQ-TREE output             (12d)
```

Step 13 runs on Negishi under `$CLUSTER_SCRATCH/GROUSE/grouse_thex`:

```
inputs/                       copies of final_ptarmigan/ and the ptarmigan reference
work/projected/, work/msa/    haplotypes projected onto ptarmigan; one alignment per chromosome
work/thexb/                   THExBuilder windows, trimAl and pairwise-filter output
work/iqtree/<chrom>/          one IQ-TREE run per window (.contree used)
work/thex_inputs/             Tree Viewer .xlsx, chromosome BED, gene GFF, summary tables
```

The deliverable assemblies are
`final_ptarmigan/<SPECIES>_<SAMPLE>_hap{1,2}.pseudo_chr.fasta` with matching
`.liftoff.gff3`.

---

## Manifest

[`assembly_manifest.tsv`](assembly_manifest.tsv) — one row per sample, tab
separated, in this column order:

```
sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
```

`hifi_bams` may be a comma-separated list; `ont_ul` is a path or `NA`. Every
step that loops over samples reads this file and indexes it by
`SLURM_ARRAY_TASK_ID`, so the row order defines the array indices — don't
reorder it between steps.

> **Note:** species names are whitespace-trimmed when read (scripts 02, 05, 07, 11),
> because the LEPC rows once carried a trailing space (`"LEPC "`) that produced
> filenames such as `LEPC _F5545_hap1`. The manifest is now clean as well.

---

## Accessory scripts

[`accessory_scripts/`](accessory_scripts) holds the data-preparation utilities
that run before step 02: ONT ultra-long extraction and QC, HiFi re-sequencing
merges, and parallel decompression. They are run by hand, per dataset, as
needed. It also holds `plot_assembly_contiguity.R`, which draws report
Figure S5 from the final QUAST values (run locally).

---

## Window phylogenies (step 13): design notes

Step 13 moved here from the `grouse_pangenome` repo (it was that repo's step 05) because
it uses only the finished assemblies. The full rationale is in the header of
`13_thex_phylogeny.sh`; in short:

- **Reference = rock ptarmigan**, renamed to the grouse `chr_N` names. It is also the
  outgroup, so it adds no species bias, and no grouse assembly's misjoins enter the
  coordinates. Chromosomes that fewer than half the haplotypes cover (microchromosomes
  30-33, 35, 36, 38) are dropped; W and MT are excluded.
- **THExBuilder 1.0.4 pairwise-filter bug**: it computes a mismatch *count*, not a
  proportion, so its default cutoff masks any mismatch. The script detects this and
  passes 0.15 × sub-window length (15 per 100 bp), and exempts the outgroup.
- **TopologyID** is the species-level relationship, not THEx's topobinner (with 47 tips
  nearly every window is a unique topology).
- Install THEx from GitHub with `pip --no-deps` into a pinned env (`setup`); the
  `ajharris_2374` conda package no longer solves.

## Notes on interpreting the QC

**BUSCO completeness measures the assembly, not the reference.** Reference-guided
ordering cannot change gene content — step 05 only reorders and renames what
hifiasm and yahs produced. So re-running step 07 after changing reference should
reproduce the previous score to within rounding. A score that *moves* means
something went wrong in step 05 (sequence dropped by the length filter, or a
mangled rename) and should be chased before anything downstream is believed.

**Six haplotypes are ~100 Mb short and ~4 points lower in BUSCO.** This is not
an assembly failure: they are the Z-bearing or W-bearing haplotype of a ZW
female, and phasing splits Z and W into separate haplotypes. Step 10 tests this
three independent ways, and step 05's depth check corroborates it from the reads.

**One join is still unresolved.** The junction inside `chr_2` (the chicken
chromosome 3 homolog), splitting it into ~7.4 Mb and ~102.6 Mb pieces, is
*unsupported* by Hi-C in 17 of the 46 assemblies — across all three species and
against both references — while ptarmigan itself carries the chromosome intact.
Against a measured false-positive rate of 0.88%, seeing it 17 times is not
chance (binomial P(≥17 of 46 | p = 0.0088) ≈ 1.6 × 10⁻²³). It is flagged rather than resolved: it may be a real
polymorphism in these species, or a shared assembly artifact, and distinguishing
those needs evidence this pipeline does not produce.
