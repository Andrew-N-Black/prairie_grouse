# chicken_reference_chain/ (archived, reference only)

The first version of the assembly pipeline, which ordered and annotated the
46 haplotype assemblies against **chicken (GRCg7b)**. It was replaced on
2026-10-03 (commit `660a2df`) by the rock ptarmigan chain in the repository
root, and **none of these scripts produced the assemblies, figures or tables
in the USFWS report.** They are restored here from git history
(`660a2df^`) unchanged, so the chicken-guided runs can be reproduced or
compared; do not submit them as part of the current pipeline.

| File | What it did |
|---|---|
| `01_download_chicken_ref.sh` | GRCg7b FASTA + GFF3 for RagTag/Liftoff |
| `01_download_reference_genomes.sh` | chicken + LEPC reference download (later version of the above) |
| `02_genome_assembly_array.sh` | hifiasm -> yahs -> RagTag (chicken) -> Liftoff (chicken) -> QC, array x23 |
| `02_test_single_sample.sh` | single-sample test of the array script |
| `03_pangenome_analysis.sh` | early 3-species pangenome draft on the chicken-ordered assemblies |
| `04_shortread_introgression.sh` | early short-read mapping to that pangenome |
| `shortread_manifest.tsv` | sample sheet for step 04 |
| `README_chicken_pipeline_original.md` | the README that accompanied this version |

Why chicken was dropped (details in the root README): chicken RagTag
ordering imposed chicken's karyotype -- a GGA4p fusion in 28 of 46
haplotypes that Hi-C contradicted in every case, plus the GGA6/GGA8 and
GGA2 differences -- and Liftoff left 95,451 features unmapped versus 22,010
against ptarmigan. The pangenome (Objective 4) is developed in its own
repository on the ptarmigan-ordered assemblies.

Paths, SLURM account and conda environments are as they were at the time
(`grouse_asm/final`, later `chicken_guided/final`).
