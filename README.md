# Prairie Grouse Genome Assembly and QC

SLURM/bash pipeline (Purdue RCAC, Gautschi) for **haplotype-resolved reference
genome assembly and quality control** of three prairie grouse species:

| Species | Code | Scientific name |
|---|---|---|
| Lesser Prairie-Chicken | LEPC | *Tympanuchus pallidicinctus* |
| Greater Prairie-Chicken | GRPC | *T. cupido* |
| Sharp-tailed Grouse | STGR | *T. phasianellus* |

**n = 23 individuals × 2 haplotypes = 46 assemblies.** PacBio HiFi + Hi-C
(Omni-C) throughout, with ONT ultra-long reads for three samples.

Assemblies are ordered against the **rock ptarmigan** reference
(*Lagopus muta*, bLagMut1, `GCF_023343835.1`), and annotations are transferred
from it.

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
  all three *Tympanuchus* species (ptarmigan chr_5; 65.85 Mb in grouse against
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
annotation transfer: **22,010 unmapped Liftoff features against ptarmigan
versus 95,451 against chicken**, a 4.3× reduction.

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

Run in order. Steps 01 and 03 are one-off setup; 02, 05, 07 and 11 are array
jobs over the 23 samples; the rest are single jobs.

| Step | Script | What it does | Shape |
|---|---|---|---|
| 01 | [`01_download_chicken_reference.sh`](01_download_chicken_reference.sh) | Chicken GRCg7b FASTA + accession→chromosome map. Needed **only** by step 04, to tie ptarmigan chromosome numbers to the chicken karyotype. No annotation is downloaded — none is transferred from chicken. | single, 4 h |
| 02 | [`02_assembly_scaffold_array.sh`](02_assembly_scaffold_array.sh) | **The reference-free half.** HiFiAdapterFilt → hifiasm (Hi-C-phased, `--ul` where ONT exists) → bwa mem Hi-C alignment → yahs scaffolding → Juicebox `.hic`/`.assembly` for manual curation. QUAST after hifiasm and after yahs. Also builds every conda env the later steps use. | array ×23, 48 cpu, 160 G, 10 d |
| 03 | [`03_ptarmigan_reference.sh`](03_ptarmigan_reference.sh) | Rock ptarmigan FASTA + GFF3 + accession→chromosome map, resolving the NCBI directory name at runtime. Prints the karyotype it found and flags size-rank naming. | single, 2 h |
| 04 | [`04_chromosome_homology.sh`](04_chromosome_homology.sh) + [`.py`](04_chromosome_homology.py) | Aligns ptarmigan to chicken (`minimap2 -x asm20`) and works out which chicken chromosome each ptarmigan chromosome corresponds to. **This is where the karyotype findings above come from.** Writes the correspondence table, the summary, and a rename map. | single, 24 cpu, 8 h |
| 05 | [`05_ragtag_liftoff_array.sh`](05_ragtag_liftoff_array.sh) | **The reference-guided half.** RagTag ordering → chromosome naming → QUAST → Liftoff → orientation dot plot → Hi-C realignment + Pretext contact map → tidk → HiFi depth check. Reads step 02's yahs scaffolds; never writes to them. | array ×23, 24 cpu, 96 G, 3 d |
| 06 | [`06_rename_chromosomes.sh`](06_rename_chromosomes.sh) + [`.py`](06_rename_chromosomes.py) | Renames chromosomes consistently across every derived file (FASTA, GFF3, AGP, CRAM, BUSCO tables). **Only needed if step 05 ran with `CHR_NAMING=homology`** — its default produces the final names directly. | single, 4 h |
| 07 | [`07_busco_array.sh`](07_busco_array.sh) | Splits chr_Z / chr_W / chr_MT into their own FASTAs, then BUSCO v5.4.7 genome mode (metaeuk) against `aves_odb10` (8,338 orthologs) on the full assembly. | array ×23, 24 cpu, 2 d |
| 08 | [`08_ptarmigan_reference_busco.sh`](08_ptarmigan_reference_busco.sh) | Runs the ptarmigan reference through the identical BUSCO treatment, so it can be drawn alongside the grouse and so every ortholog gets a reference chromosome. Required by step 10's third check. | single, 24 cpu, 1 d |
| 09 | [`09_busco_plots.sh`](09_busco_plots.sh) + [`.py`](09_busco_plots.py) | BUSCO-Plot-Py figures: per-species completeness barplots, karyoplots, and pairwise horizontal synteny plots. | single, 8 cpu, 4 h |
| 10 | [`10_sex_chromosome_check.sh`](10_sex_chromosome_check.sh) + [`.py`](10_sex_chromosome_check.py) | Tests whether the short, low-BUSCO haplotypes are short because they lack the Z. Three independent checks: chromosome content, BUSCOs on chr_Z, and where the *missing* orthologs sit on the reference. | single, 30 min |
| 11 | [`11_join_support.sh`](11_join_support.sh) + [`.py`](11_join_support.py) | Scores **every join RagTag made** against the Hi-C: observed cross-join contacts versus what contiguous sequence would carry at the same separation, calibrated per library on pseudo-joins inside intact scaffolds. This is what catches a reference imposing a join the data does not support. | array ×23, 8 cpu, 8 h |

### Submitting

```bash
N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)   # 23

bash   01_download_chicken_reference.sh
bash   03_ptarmigan_reference.sh
sbatch 04_chromosome_homology.sh

sbatch --array=0-$((N-1))%6 02_assembly_scaffold_array.sh
# ... then, once every sample has finished:
sbatch --array=0-$((N-1))%6 05_ragtag_liftoff_array.sh

sbatch --array=0-$((N-1))%8 07_busco_array.sh
sbatch 08_ptarmigan_reference_busco.sh
sbatch 09_busco_plots.sh                       # after 07 and 08
sbatch 10_sex_chromosome_check.sh              # after 07 and 08
sbatch --array=0-$((N-1))%8 11_join_support.sh  # needs 05 with DO_HIC=true
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

`hifi_bams` may be a comma-separated list; `ont_ul` may be empty. Every step
that loops over samples reads this file and indexes it by `SLURM_ARRAY_TASK_ID`,
so the row order defines the array indices — don't reorder it between steps.

---

## Accessory scripts

[`accessory_scripts/`](accessory_scripts) holds the data-preparation utilities
that run before step 02: ONT ultra-long extraction and QC, HiFi re-sequencing
merges, and parallel decompression. They are run by hand, per dataset, as
needed.

---

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
*unsupported* by Hi-C in 15 of the 46 assemblies — across all three species and
against both references — while ptarmigan itself carries the chromosome intact.
Against a measured false-positive rate of 0.88%, seeing it 15 times is not
chance (P ≈ 1.5 × 10⁻³¹). It is flagged rather than resolved: it may be a real
polymorphism in these species, or a shared assembly artifact, and distinguishing
those needs evidence this pipeline does not produce.
