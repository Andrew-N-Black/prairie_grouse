# Prairie grouse: haplotype-resolved genomes, pangenome and the Z chromosome

SLURM/bash pipeline (Purdue RCAC; Gautschi, step 15 on Negishi) for the manuscript
*Species boundaries in hybridizing prairie grouse are concentrated on the Z chromosome*
(in preparation for *Science*): haplotype-resolved assemblies of 23 birds, a multi-species
pangenome, structural variants, window phylogenies, diversity and divergence, mitogenomes
and species-diagnostic Z markers.

| Species | Code | Scientific name |
|---|---|---|
| Lesser Prairie-Chicken | LEPC | *Tympanuchus pallidicinctus* |
| Greater Prairie-Chicken | GRPC | *T. cupido* |
| Sharp-tailed Grouse | STGR | *T. phasianellus* |

**n = 23 birds × 2 haplotypes = 46 assemblies** (7 LEPC, 6 GRPC, 10 STGR; six females),
PacBio HiFi + Hi-C, ordered against the **rock ptarmigan** (*Lagopus muta*, bLagMut1,
`GCF_023343835.1`). Sequence data: NCBI BioProject PRJNA1513026. SLURM account `dewoody`
(Gautschi), `fnrdewoody` (Negishi).

This repository merges the former `prairie_grouse` (assembly, QC, phylogenies,
mitogenomes, markers) and `grouse_pangenome` (pangenome, structural variants, diversity)
repositories, with the full history of both. Steps are numbered in the order of the
manuscript's Methods; the [old → new table](#renumbering-old--new) maps the previous
numbers.

---

## Layout

| Path | Contents |
|---|---|
| `01_` … `22_` | the numbered pipeline (below) |
| [`shortread_panel/`](shortread_panel) | **unnumbered, not in the manuscript yet**: short-read panel mapping, graph genotyping and introgression tests, marker validation in shotgun birds, panel mitochondrial read check (see [Short-read panel](#short-read-panel-unnumbered)) |
| [`accessory/`](accessory) | data preparation (ONT, HiFi merging, decompression), container pulls, outgroup download, the HiFi-only re-assembly of F5545/F5597 |
| [`ncbi_submission/`](ncbi_submission) | GenBank submission of the 46 assemblies (FASTA prep, per-gap linkage evidence, `.sqn`) |
| [`docs/`](docs) | working notes (`misjoin_fix_notes.md`, dated Oct 2026; its step numbers are the pre-merge ones) |
| [`archive/chicken_reference_chain/`](archive/chicken_reference_chain) | the superseded chicken-guided pipeline, unchanged, for reference only |
| `assembly_manifest.tsv` | one row per bird (see [Manifest](#manifest)) |

---

## Pipeline

Run from the repository root. Every step reads and writes under
`$CLUSTER_SCRATCH/GROUSE/grouse_asm` (the THEx step 15 under `grouse_thex` on Negishi).
**Steps 16–19 are placeholders**: the analyses are in the manuscript but their code has
not yet been added here; the numbers are reserved so that nothing is renumbered later.

### A. References and assembly

| Step | Script | What it does | Shape |
|---|---|---|---|
| 01 | [`01_download_chicken_reference.sh`](01_download_chicken_reference.sh) | Chicken GRCg7b FASTA + accession→chromosome map, used **only** by step 03 to relate ptarmigan chromosomes to the chicken karyotype. | single, 4 h |
| 02 | [`02_ptarmigan_reference.sh`](02_ptarmigan_reference.sh) | Rock ptarmigan FASTA + GFF3 + accession→chromosome map (NCBI directory resolved at runtime); prints the karyotype and flags size-rank naming. | single, 2 h |
| 03 | [`03_chromosome_homology.sh`](03_chromosome_homology.sh) + [`.py`](03_chromosome_homology.py) | Ptarmigan → chicken alignment (`minimap2 -x asm20`); which chicken chromosome each ptarmigan chromosome corresponds to (table S3). | single, 24 cpu, 8 h |
| 04 | [`04_assembly_scaffold_array.sh`](04_assembly_scaffold_array.sh) | **Reference-free half.** HiFiAdapterFilt → hifiasm (Hi-C phasing; v0.19.9-r616) → bwa mem Hi-C → YaHS → Juicebox files; QUAST after hifiasm and YaHS. Builds the conda envs later steps use. | array ×23, 48 cpu, 160 G, 10 d |
| 05 | [`05_ragtag_liftoff_array.sh`](05_ragtag_liftoff_array.sh) | **Reference-guided half.** RagTag ordering on ptarmigan → chromosome naming → QUAST → Liftoff → dot plot → Hi-C realignment + Pretext map → tidk → HiFi depth. Never writes to step 04's output. | array ×23, 24 cpu, 96 G, 3 d |
| 05b | [`05b_relabel_male_W.sh`](05b_relabel_male_W.sh) | Males (ZZ) have no W: sequence RagTag placed on ptarmigan chr_W in a male haplotype becomes `scaffold_N` (≥ 50 kb) or `unplaced_short`. Females untouched. | login node |
| 05c | [`05c_rename_chromosomes.sh`](05c_rename_chromosomes.sh) + [`.py`](05c_rename_chromosomes.py) | *Optional.* Renames chromosomes across every derived file; only needed if step 05 ran with `CHR_NAMING=homology`. | single, 4 h |
| 06 | [`06_join_support.sh`](06_join_support.sh) + [`.py`](06_join_support.py) | Scores **every RagTag join** against the bird's own Hi-C (cross-join contacts vs. contiguous sequence at the same separation, calibrated per library): supported / unsupported / ambiguous / not scored. | array ×23, 8 cpu, 8 h |
| 07 | [`07_break_unsupported_joins.sh`](07_break_unsupported_joins.sh) + [`.py`](07_break_unsupported_joins.py) | Cuts every join step 06 called `unsupported` (989 cuts) and writes the pangenome-input FASTAs (`pangenome_inputs/`), after proving AGP coordinates still match. | single, < 1 h |
| 08 | [`08_chrom_qc.sh`](08_chrom_qc.sh) + [`.py`](08_chrom_qc.py), [`08b_fix_misjoins.py`](08b_fix_misjoins.py) | Chromosome-scale misjoin screen: every haplotype aligned to an anchor assembly (GRPC F5595; `minimap2 -x asm10`), blocks ≥ 500 kb by anchor chromosome; 08b cuts misjoins at the nearest scaffold gap (24 cuts) and relabels the pieces. Notes: [`docs/misjoin_fix_notes.md`](docs/misjoin_fix_notes.md). | array + login node |

### B. Assembly quality and sex chromosomes

| Step | Script | What it does | Shape |
|---|---|---|---|
| 09 | [`09_busco_array.sh`](09_busco_array.sh) | Splits chr_Z / chr_W / chr_MT into their own FASTAs; BUSCO v5.4.7 (`aves_odb10`, 8,338 orthologs) on each assembly. | array ×23, 24 cpu, 2 d |
| 09b | [`09b_ptarmigan_reference_busco.sh`](09b_ptarmigan_reference_busco.sh) | The ptarmigan reference through the identical BUSCO run. | single, 1 d |
| 09c | [`09c_busco_plots.sh`](09c_busco_plots.sh) + [`.py`](09c_busco_plots.py) | BUSCO-Plot-Py completeness barplots, karyoplots and synteny (**figs. S4–S5**). | single, 4 h |
| 10 | [`10_sex_chromosome_check.sh`](10_sex_chromosome_check.sh) + [`.py`](10_sex_chromosome_check.py) | Z/W content per haplotype from chromosome content, Z BUSCOs and where missing orthologs sit (**table S2**; genomic sex, confirmed by CHD1 PCR, fig. S3). | single, 30 min |
| 11 | [`11_assembly_qc_array.sh`](11_assembly_qc_array.sh) | Merqury (k = 21: QV, completeness, spectra-cn) and HiFi depth (hap1 + hap2 together; windows < 0.5× or > 1.5× the median flagged) on the step-07 FASTAs; optional FCS-adaptor. | array ×23 |
| 11b–e | [`11b_collate_qc.sh`](11b_collate_qc.sh), [`11c_qc_report_table.sh`](11c_qc_report_table.sh), [`11d_plot_qc_report.R`](11d_plot_qc_report.R), [`11e_plot_assembly_contiguity.R`](11e_plot_assembly_contiguity.R) | One QC table with REVIEW flags; report table (per haplotype and hap1+hap2); **fig. S2** (completeness, QV, depth); **fig. S1** (contiguity, from QUAST). | login node / local R |

### C. Pangenome and structural variants

| Step | Script | What it does | Shape |
|---|---|---|---|
| 12 | [`12_pangenome_build.sh`](12_pangenome_build.sh) | PanSN renaming; PGGB graph per chromosome (seqwish v0.7.9 from `pggb.sif`), merged with `odgi squeeze`; odgi stats and path lengths; panacus growth/core by species; RepeatModeler per species → RepeatMasker on all 46; per-haplotype and per-species genome, gene and TE tables (**table S4**). | array per chromosome + summary |
| 12b | [`12b_fig3_pangenome.py`](12b_fig3_pangenome.py) | **Fig. S6** (autosomal pangenome growth and core; autosomal repeat and gene content per haplotype) from step 12's `comparison/` tables (`STAGE=autosomes`); writes `FigS6_autosomal_content.png/pdf`. | local |
| 13 | [`13_structural_variants.sh`](13_structural_variants.sh) + [`13_sv_summary.py`](13_sv_summary.py) | `vg deconstruct` (v1.59) per chromosome against one reference path (longest single-sequence LEPC haplotype; STGR for chr_W), `vcfbub -l 0 -a 100000`, `vcfwave -I 1000`; SVs ≥ 50 bp typed INS/DEL/INV/COMPLEX; species-diagnostic SVs and fixed differences (all haplotypes and subsampled to six); Z vs autosomes (Fisher's exact test, density) (**table S8**). | array per chromosome, 96 cpu; summary |
| 14 | [`14_z_inversions.sh`](14_z_inversions.sh) + [`.py`](14_z_inversions.py) | Large Z inversions: every Z path (`odgi paths -f`) aligned to the step-13 reference Z (`minimap2 -x asm10`); inverted blocks ≥ 50 kb clustered into events, classified, flagged near gaps/path ends; optional liftover to ptarmigan Z coordinates (needed by step 20b). | single, 32 cpu, < 1 h |

### D. Phylogenies, population structure, coding sequence and diversity

| Step | Script | What it does | Shape |
|---|---|---|---|
| 15 | [`15_thex_phylogeny.sh`](15_thex_phylogeny.sh) + [`.py`](15_thex_phylogeny.py) | **Negishi.** Every haplotype projected onto ptarmigan (`minimap2 -x asm20`), one alignment per chromosome, THExBuilder 100-kb windows → trimAl → pairwise filter → IQ-TREE 2 per window; Tree Viewer file with species-level `TopologyID` (**Fig. 2**, **Table 1**). | array → node → array → small job |
| 15b | [`15b_thex_stats.py`](15b_thex_stats.py) | Table 1 / Fig. 2 statistics: Z vs autosome rank, Spearman vs chromosome length and gene density (Fig. 2C), window-level GEE. | local |
| 15c | [`15c_window_diagnostics.py`](15c_window_diagnostics.py) | Why a region has no trees: per-window N/repeat/GC, haplotype coverage, filter fate. | login node |
| 15d | [`15d_thex_rescue.sh`](15d_thex_rescue.sh) + [`.py`](15d_thex_rescue.py) | Re-runs windows lost because 1–4 haplotypes lack that stretch, without those haplotypes. | node → array → small job |
| 15e | [`15e_summary_trees.sh`](15e_summary_trees.sh) | Weighted ASTRAL summary trees (autosomes, Z) + gene concordance factors. | standby job |
| 15f | [`15f_plot_summary_trees.R`](15f_plot_summary_trees.R) | **Fig. 3**: autosomal and Z summary trees beside the step-21d mitogenome tree. | local R |
| 16 | *placeholder — code to be added* | **SNP calling from HiFi reads**: pbmm2 (HiFi preset) to the ptarmigan reference, DeepVariant v1.10.0 (PACBIO model) per bird, GLnexus joint genotyping, filters. | |
| 17 | *placeholder — code to be added* | **Population structure**: PLINK PCA of LD-pruned SNPs (autosomes, 23 birds; Z, 17 males) and ADMIXTURE K = 2–6 with cross-validation (**Fig. 4A–D**, **fig. S7**). | |
| 18 | *placeholder — code to be added* | **Sliding-window PCA and diagnostic regions**: WinPCA v1.2 (100-kb windows, 10-kb step), diagnostic score D(w), top-1% windows merged into regions (**Fig. 4E–H**, **fig. S8**, **Table 2**, **table S5**). | |
| 19 | *placeholder — code to be added* | **Coding-sequence divergence and polymorphism**: CDS of all Z genes and 2,000 autosomal genes from each haplotype's Liftoff annotation; Nei–Gojobori dN/dS, πN/πS, fixed differences, McDonald–Kreitman tests; 41–44-Mb segment vs. blocks of consecutive Z genes (**Table 3**, **tables S6–S7**, **fig. S9**). | |
| 20 | [`20_pixy_divergence.sh`](20_pixy_divergence.sh) + [`20_pixy_summary.py`](20_pixy_summary.py) | π, d<sub>XY</sub>, d<sub>a</sub> and Hudson's F<sub>ST</sub>: each bird's HiFi reads mapped (`minimap2 -x map-hifi`) to the step-13 reference paths (STGR W as a decoy); all-sites calling with bcftools in 5-Mb chunks; pixy 50-kb windows (Z also males only). Summary excludes windows with < 5,000 sites or π/d<sub>XY</sub> > 0.02, Z/autosome ratios with 1-Mb block-bootstrap CIs, π vs the neutral ZW range 0.5625–0.75 (**table S9**). | samples → ref → map (array per bird) → call (array per chunk) → summary |
| 20b | [`20b_fig7_diversity_divergence.py`](20b_fig7_diversity_divergence.py) | **Fig. 5** (F<sub>ST</sub> genome scan; π and divergence by chromosome class; the Z in ptarmigan coordinates with the 41–44-Mb segment). Run by step 20's summary stage when step 14's `ref_to_ptarmigan.paf` exists. | local or summary job |
| 20c | [`20c_figS6_genome_scans.py`](20c_figS6_genome_scans.py) | **Fig. S10**: F<sub>ST</sub>, d<sub>XY</sub> and π along the genome for each species pair (species colours; Z π from males). Run by step 20's summary stage. | local or summary job |

### E. Mitogenomes and Z markers

| Step | Script | What it does | Shape |
|---|---|---|---|
| 21 | [`21_mitohifi_array.sh`](21_mitohifi_array.sh) | One circularized, annotated mitogenome per bird from HiFi reads (MitoHiFi v3; NUMT-bearing long reads removed). `setup` fetches the reference mitogenome. | array ×23, 32 cpu, 1 d |
| 21b | [`21b_mito_from_hic.sh`](21b_mito_from_hic.sh) | Mitogenome as a Hi-C read consensus for birds whose HiFi reads are dominated by NUMTs (STGR F5485, F5502). | array (2 rows) |
| 21c | [`21c_mito_qc.sh`](21c_mito_qc.sh) | QC table, `trim`/`rescue` fixes, pairwise identities, `collect` of the passing mitogenomes. | login node |
| 21d | [`21d_mito_tree.sh`](21d_mito_tree.sh) | MAFFT alignment of the 23 mitogenomes + GRPC reference + ptarmigan outgroup; IQ-TREE 3 (ModelFinder, 1,000 UFBoot) (**Fig. 3C**). | short job |
| 21e | [`21e_plot_mito_tree_ggtree.R`](21e_plot_mito_tree_ggtree.R) | Stand-alone mitogenome tree (ggtree). | local R |
| 22 | [`22_z_markers.sh`](22_z_markers.sh) + [`.py`](22_z_markers.py) | **Species-diagnostic Z SNPs** in the 40 Z-bearing haplotypes (aligned to the LEPC reference Z; minimap2 + paftools), single-copy and W-gametolog checks, Primer3 primers (80–200-bp products). Its `submit-panel`, `submit-zbam` and `submit-assign` stages use the short-read panel (ANGSD; [`shortread_panel/panel_groups.tsv`](shortread_panel/panel_groups.tsv)) and are not part of the manuscript. | array → node → array → small job |
| 22b | [`22b_z_indel_markers.sh`](22b_z_indel_markers.sh) + [`.py`](22b_z_indel_markers.py) | **Gel-scorable Z indel markers**: species-diagnostic length differences in the Z-bearing haplotypes, Primer3, in-silico PCR on every Z, the female W and the whole LEPC reference. | single, 16 cpu, 4 h |

### Submitting

```bash
N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)   # 23

# A. references and assembly
bash   01_download_chicken_reference.sh
bash   02_ptarmigan_reference.sh
sbatch 03_chromosome_homology.sh
sbatch --array=0-$((N-1))%6 04_assembly_scaffold_array.sh          # ~10 d per bird
sbatch --array=0-$((N-1))%6 05_ragtag_liftoff_array.sh             # after 04
bash   05b_relabel_male_W.sh
sbatch --array=0-$((N-1))%8 06_join_support.sh                     # needs 05 with DO_HIC=true
DRY_RUN=true bash 07_break_unsupported_joins.sh && sbatch 07_break_unsupported_joins.sh
#   08: see docs/misjoin_fix_notes.md (map array, summarize, 08b_fix_misjoins.py)

# B. quality and sex chromosomes
sbatch --array=0-$((N-1))%8 09_busco_array.sh
sbatch 09b_ptarmigan_reference_busco.sh
sbatch 09c_busco_plots.sh && sbatch 10_sex_chromosome_check.sh
sbatch --array=0-$((N-1))%6 11_assembly_qc_array.sh
bash 11b_collate_qc.sh && bash 11c_qc_report_table.sh                # then 11d / 11e locally

# C. pangenome and structural variants
bash accessory/pull_containers.sh                                   # once: pggb, vg, RepeatMasker images
bash 12_pangenome_build.sh list                                     # prints the per-chromosome array command
STAGE=summary sbatch --export=ALL -t 4-00:00:00 12_pangenome_build.sh
bash 13_structural_variants.sh submit
export PATH=/home/blackan/.conda/envs/zindel/bin:$PATH              # minimap2 for step 14
PTARMIGAN_FA=<ptarmigan.fna> sbatch --export=ALL 14_z_inversions.sh

# D. phylogenies (Negishi) and diversity
bash 15_thex_phylogeny.sh transfer && bash 15_thex_phylogeny.sh setup && bash 15_thex_phylogeny.sh submit
bash 15d_thex_rescue.sh submit                                      # optional
bash 15e_summary_trees.sh setup && bash 15e_summary_trees.sh submit
bash 20_pixy_divergence.sh samples                                  # check samples.tsv
sbatch --export=ALL,STAGE=ref 20_pixy_divergence.sh
bash 20_pixy_divergence.sh map && bash 20_pixy_divergence.sh call  # summary + Fig. 5, fig. S10 run afterwards

# E. mitogenomes and markers
bash 21_mitohifi_array.sh setup && sbatch --array=0-$((N-1))%8 21_mitohifi_array.sh
sbatch -q standby -t 04:00:00 --array=7,8 21b_mito_from_hic.sh
bash 21c_mito_qc.sh qc && bash 21c_mito_qc.sh collect && bash 21d_mito_tree.sh
bash 22_z_markers.sh setup && bash 22_z_markers.sh submit
mkdir -p logs && sbatch 22b_z_indel_markers.sh
```

Step 05's stages can be switched individually (`DO_LIFTOFF`, `DO_QUAST`, `DO_HIC`,
`DO_DOTPLOT`, `DO_TIDK`, `DO_DEPTH`, `CHR_NAMING`), e.g. a fast structural pass first:
`sbatch --export=ALL,DO_HIC=false,DO_DEPTH=false --array=0-22%6 05_ragtag_liftoff_array.sh`.
`DO_HIC=true` is required before step 06.

---

## Manuscript items and the steps that produce them

Figure numbers follow the current draft; script and output file names (e.g. `12b_fig3_…`, `20b_fig7_…`, `20c_figS6_…`) keep their original numbers.

| Manuscript | Steps |
|---|---|
| Assembly, scaffolding, QC (Methods) | 04, 05, 05b, 06, 08, 09, 11 |
| Fig. 1 (sampling map); fig. S1, table S1 (assembly contiguity) | — ; 05 (QUAST), 11e |
| Figs. S4–S5 (karyotype, synteny); table S3 | 09c; 03 |
| Table S2, fig. S3 (sex chromosomes) | 10 (+ CHD1 PCR) |
| Fig. S2 (Merqury, depth) | 11–11d |
| Fig. S6, table S4 (pangenome) | 12, 12b |
| Table S8 (structural variants); Z inversions | 13; 14 |
| Fig. 2, Table 1 (window trees) | 15, 15b–15d |
| Fig. 4, Table 2, figs. S7–S8, table S5 (SNPs, PCA, ADMIXTURE, WinPCA) | 16–18 *(code to be added)* |
| Table 3, tables S6–S7, fig. S9 (coding sequence) | 19 *(code to be added)* |
| Fig. 3 (summary and mitogenome trees) | 15e, 15f, 21d |
| Fig. 5, fig. S10, table S9 (π, d<sub>XY</sub>, F<sub>ST</sub>) | 20, 20b, 20c (+ 14 for ptarmigan coordinates) |
| Fig. S11, table S10 (hybrid index, Z vs autosomal introgression; neutral simulations) | `shortread_panel/`: hybrid_index, figS11_hybrid_index, sim_z_autosome_neutral |
| Mitogenomes | 21–21e |
| Candidate Z-linked markers | 22, 22b |
| GenBank submission | `ncbi_submission/` |

---

## Short-read panel (unnumbered)

[`shortread_panel/`](shortread_panel) holds the analyses that use population short-read
data rather than the 23 assemblies. They are kept but **not numbered**, because they may
or may not be part of this manuscript: introgression across the LEPC–GRPC hybrid zone
cannot be tested without high-depth LEPC from the zone, and those birds cannot be sampled
and sequenced until spring 2027.

| Script | What it does |
|---|---|
| [`shortread_introgression.sh`](shortread_panel/shortread_introgression.sh) + [`shortread_manifest.tsv`](shortread_panel/shortread_manifest.tsv) | Graph GBZ + giraffe indexes, `vg giraffe` → `vg pack` → `vg call` per sample, joint VCF, Dsuite `Dtrios` (genome-wide D) and `Dinvestigate` (windowed D/f<sub>d</sub>). Needs at least one `Outgroup` row (Gunnison sage-grouse, `accessory/download_outgroup.sh`); Dsuite's untreed mode infers the trio from BBAA. Run from the repo root. |
| [`validate_indel_markers.sh`](shortread_panel/validate_indel_markers.sh) + [`.py`](shortread_panel/validate_indel_markers.py) | Genotypes the step-22b indel markers in higher-coverage shotgun birds and checks the alleles against ptarmigan, chicken and the grouse W. |
| [`validate_snp_markers.sh`](shortread_panel/validate_snp_markers.sh) + [`.py`](shortread_panel/validate_snp_markers.py) | The same for the step-22 SNP markers; per-bird Z ancestry. |
| [`mito_read_check.sh`](shortread_panel/mito_read_check.sh) | Which mitochondrial lineage a panel bird's reads carry, at LEPC-diagnostic sites. |
| [`hybrid_index.sh`](shortread_panel/hybrid_index.sh) + [`.py`](shortread_panel/hybrid_index.py) | Hybrid index (S) and interspecific heterozygosity (H) for every panel bird on the autosomes and the Z, from ANGSD genotype likelihoods at markers that differ between allopatric reference groups (dAF ≥ 0.8; reference birds scored leave-one-out; 1-Mb block bootstrap). Females' Z is scored as one copy. Triangle plots and Z vs autosomal S for LEPC–GRPC. Leaves out the 2026 birds listed in [`panel_exclude_2026.txt`](shortread_panel/panel_exclude_2026.txt) (held for the hybrid-zone paper); only the 24 allopatric 2026 GRPC are used, as parental reference birds. Run from `shortread_panel/`: `bash hybrid_index.sh setup`, then `submit`. |
| [`panel_groups.tsv`](shortread_panel/panel_groups.tsv) | Panel birds by species × context, used by step 22's panel stages. |
| [`sim_z_autosome_neutral.py`](shortread_panel/sim_z_autosome_neutral.py) | Neutral expectation for foreign ancestry on the Z vs the autosomes: individual-based Wright–Fisher model with ZW inheritance (Z recombines only in males), ongoing GRPC/LEPC immigration with male-, unbiased or female-only GRPC migrants, F1-female sterility, and single pulses. Reports Z/A ancestry ratio and the fraction of admixed birds with no foreign Z. Runs locally (`python3 sim_z_autosome_neutral.py --reps 20`). |
| [`figS11_hybrid_index.py`](shortread_panel/figS11_hybrid_index.py) | **Fig. S11** (S vs H on the autosomes and the male Z; Z vs autosomal S; S against distance to the GRPC range) and **table S10** from `hybrid_index.tsv` and `panel_groups_lepc_grpc.tsv` (plus maternal lineage from `mito_panel.tsv` when given; birds in `qc_excluded_birds.tsv` left out). |
| [`qc_excluded_birds.tsv`](shortread_panel/qc_excluded_birds.tsv) | The 51 birds excluded from the hybrid-index and mitochondrial results (mixed mitochondrial lineages at ≥ 10 of 47 sites and/or heterozygosity > 3 IQR above their population's median: signs of DNA from two birds). |
| [`mito_panel.sh`](shortread_panel/mito_panel.sh) | Runs `mito_read_check.sh` over all hybrid-index birds in parallel (`prep`, `submit`, `status`) → `mito_panel.tsv`: LEPC vs GRPC/STGR maternal lineage per bird. |
| [`panel_groups_lepc_grpc.tsv`](shortread_panel/panel_groups_lepc_grpc.tsv) | Reference sets for the hybrid-index estimate: LEPC ≥ ~90 km from the GRPC range as the LEPC reference, Edwards, Kiowa, Comanche (KS) and Cheyenne (CO) LEPC (18–70 km) as `LEPC_edge`, Stanley Co. SD GRPC (contact with STGR only) scored as `GRPC_control`, a no-LEPC-contact baseline; `km_GRPC_range_map` from the placeholder range map. |
| [`panel_sex.tsv`](shortread_panel/panel_sex.tsv) | Sex (from Z/autosome read depth), county, state, coordinates and year of each panel bird. |

---

## Why ptarmigan and not chicken

Reference-guided scaffolding imposes the reference's chromosome structure wherever the
karyotypes differ. Aligning chicken and ptarmigan (step 03) and testing each join against
the grouse Hi-C (step 06) found three differences:

- **Chicken chromosomes 6 and 8 are one chromosome** in ptarmigan and all three
  *Tympanuchus* (ptarmigan chr_5; 65.82 Mb median in grouse).
- **The short arm of chicken chromosome 4 (GGA4p) is a free-standing chromosome** here
  (ptarmigan chr_13, 19.05 Mb). Chicken-ordered assemblies reproduced the fusion in 28 of
  46 haplotypes and the Hi-C contradicted it in every one.
- **Chicken chromosome 2 corresponds to two chromosomes** (ptarmigan chr_3 and chr_7).

Ptarmigan is also the closer relative (both Tetraoninae): a mean of 477 unmapped Liftoff
features per haplotype against ptarmigan versus 95,451 in total against chicken.

**Chromosome numbers are size ranks, not homology.** bLagMut1 numbers chromosomes in
descending length, so ptarmigan chr_6 is not chicken chromosome 6. The assemblies carry
the ptarmigan numbers; relate them to the galliform literature through step 03's
correspondence table (table S3), never by assuming the numbers match.

---

## Output layout

Under `$CLUSTER_SCRATCH/GROUSE/grouse_asm`:

```
ref/                          both references, indexes, chromosome maps            (01-02)
conda_envs/                   tool envs                                            (04)
hifi_filtered/, hifiasm/, yahs/   filtered reads, contigs, Hi-C scaffolds           (04)
qc/chromosome_homology/       ptarmigan<->chicken correspondence                   (03)
ragtag_ptarmigan/, liftoff_ptarmigan/                                              (05)
final_ptarmigan/              <SP>_<S>_hap{1,2}.pseudo_chr.fasta, .liftoff.gff3     (05, 05b)
qc_ptarmigan/join_support/    per-join Hi-C verdicts                               (06)
pangenome_inputs/             graph-input FASTAs (unsupported joins and misjoins cut) (07, 08)
qc_ptarmigan/busco/, plots/, sex_chromosome_check/                                  (09-10)
pangenome_qc/                 Merqury, depth, summary tables                       (11)
pangenome/, te_annotation/, comparison/   graphs, repeats, genome/gene/TE tables    (12)
structural_variants/          per-chromosome VCFs, summary/                        (13)
z_inversions/                 events, calls, ref_to_ptarmigan.paf                  (14)
pixy_divergence/              BAMs, all-sites VCF chunks, pixy, summary/           (20)
mitogenomes/                  per bird, mito_qc.tsv, mitogenomes_pass.fasta, tree/ (21)
ncbi_submission/              GenBank files, .sqn                                  (ncbi_submission/)
```

Step 15 runs on Negishi under `$CLUSTER_SCRATCH/GROUSE/grouse_thex`
(`inputs/`, `work/projected/`, `work/msa/`, `work/thexb/`, `work/iqtree/<chrom>/`,
`work/thex_inputs/`). Output folder names are unchanged by the renumbering.

---

## Manifest

[`assembly_manifest.tsv`](assembly_manifest.tsv), one row per bird, tab separated:

```
sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
```

`hifi_bams` may be comma-separated; `ont_ul` is a path or `NA` (all `NA`: no assembly
uses ONT; see `accessory/rerun_without_ONT.sh` for F5545 and F5597). Array steps index
this file by `SLURM_ARRAY_TASK_ID`, so do not reorder it. Species names are
whitespace-trimmed when read.

---

## Design notes

**BUSCO measures the assembly, not the reference.** Step 05 only reorders and renames
what hifiasm and YaHS produced, so a BUSCO score that changes after re-scaffolding means
something went wrong in step 05.

**Six haplotypes are ~100 Mb short.** In all six females the Z and W are on the same
haplotype; the other haplotype lacks the Z (step 10).

**Joins are cut before the graph (step 07).** A graph builder places sequence by
alignment and gains nothing from reference order, while an incorrect join becomes a false
structural variant. The TE and gene comparisons in step 12 still read the uncut
`final_ptarmigan/` assemblies.

**One join remains unresolved**: a junction inside chr_2 is unsupported by Hi-C in 17 of
46 assemblies (binomial P ≈ 1.6 × 10⁻²³ against the 0.88% false-positive rate); it may be
a real polymorphism or a shared assembly artifact.

**Window phylogenies (step 15).** Reference = rock ptarmigan renamed to grouse `chr_N`
names (also the outgroup). Chromosomes covered by fewer than half the haplotypes are
dropped. THExBuilder 1.0.4 computes a mismatch count, not a proportion, for its pairwise
filter; the script passes 0.15 × sub-window length. `TopologyID` is the species-level
relationship. THEx is installed with `pip --no-deps` into a pinned env.

**Pangenome graphs (step 12).** Built per chromosome; graphs from a seqwish development
build were broken and all were rebuilt with seqwish v0.7.9 from `pggb.sif` (`--resume`).
`odgi untangle` aborts on these graphs (sdsl assertion), so step 14 uses whole-chromosome
alignment instead.

**Diversity (step 20).** Windows with π or d<sub>XY</sub> > 0.02 (collapsed repeats or
paralogous mapping; mostly chromosomes 9–37, chromosome ends and the first 3 Mb of the
reference Z) are excluded and listed in `summary/excluded_windows.tsv`. For a ZW system the
neutral Z/A diversity ratio is 0.75 with an equal breeding sex ratio and falls to 0.5625
when few males breed.

**Tooling.** Gautschi has no mosdepth module and its merqury module lacks `meryl`; step 11
fetches release binaries once (`INSTALL_ONLY=true bash 11_assembly_qc_array.sh`). pggb,
odgi, vg and RepeatMasker run from Apptainer images (`accessory/pull_containers.sh`).
panacus is downloaded as a static binary. minimap2 for steps 14 and 22 comes from a conda
env (`export PATH=<env>/bin:$PATH`; never overwrite `PATH`).

---

## Renumbering (old → new)

| Old (`prairie_grouse`) | New | | Old (`grouse_pangenome`) | New |
|---|---|---|---|---|
| 01 chicken reference | 01 | | 01 break unsupported joins | 07 |
| 03 ptarmigan reference | 02 | | `scripts/chrom_qc/` | 08, 08b |
| 04 chromosome homology | 03 | | 02, 02b–02d assembly QC | 11, 11b–11d |
| 02 assembly + scaffolding | 04 | | 03 pangenome build | 12 |
| 05, 05b | 05, 05b | | 05 structural variants | 13 |
| 06 rename chromosomes | 05c | | 06 Z inversions | 14 |
| 11 join support | 06 | | 07 pixy (+ Fig. 7 script) | 20, 20b |
| 07, 08, 09 BUSCO | 09, 09b, 09c | | 04 short-read introgression | `shortread_panel/` |
| 10 sex chromosomes | 10 | | `accessory/ncbi_*` | `ncbi_submission/` |
| `accessory_scripts/plot_assembly_contiguity.R` | 11e | | `accessory/` (containers, outgroup) | `accessory/` |
| 13, 13b–13f THEx | 15, 15b–15f | | | |
| 12, 12b–12e mitogenomes | 21, 21b–21e | | | |
| 12f mito read check | `shortread_panel/` | | | |
| 14 Z markers (+ 14g indel markers) | 22 (+ 22b) | | | |
| 14h, 14i marker validation | `shortread_panel/` | | | |
| `accessory_scripts/` | `accessory/` | | | |
| `chicken_reference_chain/` | `archive/` | | | |

Output directory names (`structural_variants/`, `z_inversions/`, `pixy_divergence/`,
`14g_z_indel_markers/`, etc.) were **not** renamed, so existing results remain valid.
