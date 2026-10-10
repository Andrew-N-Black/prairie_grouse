# Prairie Grouse Genomics Pipeline

SLURM/bash pipeline (Purdue RCAC)

**De novo assembly + pangenome + introgression** across three prairie grouse species —
   Lesser Prairie-Chicken (LEPC, *Tympanuchus pallidicinctus*), Greater Prairie-Chicken (GRPC,
   *T. cupido*), and Sharp-tailed Grouse (STGR, *T. phasianellus*).


Scripts are numbered in run order within each project. Step 01 is shared setup (both reference
genomes); 02–04 are the assembly/pangenome/introgression track; 

## Scripts

| Script | Purpose |
|---|---|
| [`01_download_reference_genomes.sh`](01_download_reference_genomes.sh) | Downloads both reference genomes this pipeline needs, once: chicken (*Gallus gallus*, GRCg7b — used by `02`'s RagTag/Liftoff steps) and LEPC (*T. pallidicinctus*, pur_lepc_1.0 — used by `05` onward). Run before either track. |
| [`02_genome_assembly_array.sh`](02_genome_assembly_array.sh) | Array job, one task per sample (23 total). Per sample: HiFiAdapterFilt (adapter trim) → hifiasm (Hi-C-phased, +ONT ultra-long for 3 samples) → yahs (Hi-C scaffolding) → Juicebox `.hic` export → RagTag (reference-guided pseudo-chromosome ordering) → Liftoff (chicken gene annotation transfer) → QC (minimap2 orientation dotplot + Pretext Hi-C contact map). Reads [`assembly_manifest.tsv`](assembly_manifest.tsv). |
| [`02_test_single_sample.sh`](02_test_single_sample.sh) | Same per-sample pipeline as `02_genome_assembly_array.sh`, but for exactly one sample defined via environment variables — run this first to validate the pipeline before committing to the full array. |
| [`03_pangenome_analysis.sh`](03_pangenome_analysis.sh) | Builds one PGGB pangenome graph across all haplotype assemblies from all three species. Compares genome length, transposable element content (RepeatModeler/RepeatMasker), and gene content (Liftoff) across species. |
| [`04_shortread_introgression.sh`](04_shortread_introgression.sh) | Maps population-level short reads to the pangenome graph (`vg giraffe`/`pack`/`call`), merges into a joint VCF, and runs `Dsuite` D-statistics to quantify introgression/hybridization between species. Requires an outgroup sample. Reads [`shortread_manifest.tsv`](shortread_manifest.tsv). |


## Manifests (fill in before running)

- [`assembly_manifest.tsv`](assembly_manifest.tsv) — one row per assembly sample: sample ID, species, raw PacBio HiFi BAM path(s), Hi-C R1/R2, optional ONT ultra-long path.
- [`shortread_manifest.tsv`](shortread_manifest.tsv) — one row per population-resequencing sample: sample ID, species (or `Outgroup`), short-read FASTQ R1/R2.


