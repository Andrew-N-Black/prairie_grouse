# Prairie Grouse Genomics Pipeline

SLURM/bash pipeline (Purdue RCAC, Negishi) for two related projects:

1. **De novo assembly + pangenome + introgression** across three prairie grouse species —
   Lesser Prairie-Chicken (LEPC, *Tympanuchus pallidicinctus*), Greater Prairie-Chicken (GRPC,
   *T. cupido*), and Sharp-tailed Grouse (STGR, *T. phasianellus*).
2. **LEPC population genomics** — a replication/extension of
   [Black et al., PNAS Nexus 2024](https://academic.oup.com/pnasnexus/article/3/8/pgae298/7720645)
   (heterozygosity, PCA, runs of homozygosity), combining the original 468-sample PNAS Nexus
   cohort with new short-read population sequencing, downsampled to matched coverage.

Scripts are numbered in run order within each project. Step 01 is shared setup (both reference
genomes); 02–04 are the assembly/pangenome/introgression track; 05–08 are the LEPC
population-genomics track.

## Scripts

| Script | Purpose |
|---|---|
| [`01_download_reference_genomes.sh`](01_download_reference_genomes.sh) | Downloads both reference genomes this pipeline needs, once: chicken (*Gallus gallus*, GRCg7b — used by `02`'s RagTag/Liftoff steps) and LEPC (*T. pallidicinctus*, pur_lepc_1.0 — used by `05` onward). Run before either track. |
| [`02_genome_assembly_array.sh`](02_genome_assembly_array.sh) | Array job, one task per sample (23 total). Per sample: HiFiAdapterFilt (adapter trim) → hifiasm (Hi-C-phased, +ONT ultra-long for 3 samples) → yahs (Hi-C scaffolding) → Juicebox `.hic` export → RagTag (reference-guided pseudo-chromosome ordering) → Liftoff (chicken gene annotation transfer) → QC (minimap2 orientation dotplot + Pretext Hi-C contact map). Reads [`assembly_manifest.tsv`](assembly_manifest.tsv). |
| [`02_test_single_sample.sh`](02_test_single_sample.sh) | Same per-sample pipeline as `02_genome_assembly_array.sh`, but for exactly one sample defined via environment variables — run this first to validate the pipeline before committing to the full array. |
| [`03_pangenome_analysis.sh`](03_pangenome_analysis.sh) | Builds one PGGB pangenome graph across all haplotype assemblies from all three species. Compares genome length, transposable element content (RepeatModeler/RepeatMasker), and gene content (Liftoff) across species. |
| [`04_shortread_introgression.sh`](04_shortread_introgression.sh) | Maps population-level short reads to the pangenome graph (`vg giraffe`/`pack`/`call`), merges into a joint VCF, and runs `Dsuite` D-statistics to quantify introgression/hybridization between species. Requires an outgroup sample. Reads [`shortread_manifest.tsv`](shortread_manifest.tsv). |
| [`05_combined_alignment_array.sh`](05_combined_alignment_array.sh) | Array job aligning the combined LEPC population cohort (468 PNAS Nexus samples + new short-read samples) to the LEPC reference genome, computing per-sample depth. Final per-sample output is CRAM. |
| [`06_downsample_and_finalize.sh`](06_downsample_and_finalize.sh) | Downsamples the new cohort's CRAMs to match the old cohort's empirical mean depth (coverage-sensitive analyses need matched depth), writes `final_cramlist.txt`. |
| [`07_heterozygosity_array.sh`](07_heterozygosity_array.sh) | Array job: per-sample individual heterozygosity via ANGSD `-doSaf` + `realSFS` (reads CRAM directly), replicating the original PNAS Nexus methodology. |
| [`08_pca_roh.sh`](08_pca_roh.sh) | PCA (`pcangsd`) and runs of homozygosity (`bcftools roh` + the original repo's `rohparser.py`), replicating the original PNAS Nexus methodology genome-wide. |

## Manifests (fill in before running)

- [`assembly_manifest.tsv`](assembly_manifest.tsv) — one row per assembly sample: sample ID, species, raw PacBio HiFi BAM path(s), Hi-C R1/R2, optional ONT ultra-long path.
- [`shortread_manifest.tsv`](shortread_manifest.tsv) — one row per population-resequencing sample: sample ID, species (or `Outgroup`), short-read FASTQ R1/R2.

## CRAM vs. BAM

CRAM is used wherever the consuming tool supports it, to cut storage — reference-based
compression matters most in the LEPC population-genomics track, where ~500 samples all align
to the same reference. Kept as BAM where a tool hard-requires it:

- **`05`–`08` (LEPC track):** the final per-sample alignment (`05`'s quality-filtered output,
  consumed by everything downstream) is CRAM. The internal `bwa mem` → sort → GATK4
  `MarkDuplicates` chain in `05` stays BAM deliberately — GATK4 `MarkDuplicates` has a
  documented, still-open CRAM reliability issue
  ([broadinstitute/gatk#5218](https://github.com/broadinstitute/gatk/issues/5218)) — and only
  converts to CRAM at the quality-filter step immediately after, via plain `samtools view -C`,
  which has no such caveat. ANGSD (`07`, `08`) reads CRAM natively through its htslib backend.
- **`02` (assembly track):** the Hi-C-to-final-assembly BAM used for the Pretext QC step is
  CRAM, since it's only ever read via `samtools view` (piped to PretextMap as SAM text). The
  Hi-C-to-contigs BAM that `yahs` scaffolds from stays BAM — `yahs` only accepts
  BAM/BED/PA5/BIN input, no CRAM support.
- **`03`/`04`:** no BAM/CRAM I/O at all (`odgi`/`pggb` graph files; `vg`'s own GAM/GBZ/pack
  formats) — not applicable.

