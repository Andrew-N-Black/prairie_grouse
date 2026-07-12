**Objective 1: Temporal comparison of LEPC in New Mexico.**

•	Overview: 

o	Short read sequencing of “old” (2018-2019) vs “new”(2026) LEPC from the same general geographic area in New Mexico. 

o	Thematic map concept for genetic load, heterozygosity, fROH (i.e., genomic susceptibility).

•	Methods:

o	Sample distribution: Map of Targeted Samples 

•	New (n=10)

•	Old (n=10) 

o	Illumina 2x150 shotgun sequencing. We can achieve ~20x DOC by multiplexing 48 samples per lane

o	Align samples to linear reference genome (published LEPC) using dragon ($) or nf-core sarek.

o	Will need to upload to NCBI asap

•	Questions: 

o	Can we actually implement thematic map concept!?

o	Which programs (and reference) to use:

• fROH-bcftools

• Individual heterozygosity-vcftools or angsd

•	Genetic Load / Rxy - VEP and/or GERPP

Objective 2: Assemble tri-species pangenome, compare between species.

•	Overview:

o	Assembly of multiple phased genome assemblies for three grouse species: Lesser Prairie Chicken, Greater Prairie Chicken, and Sharp Tailed Grouse. From these assemblies, we will do a course comparative analysis. 

•	Methods:

•	Sample distribution: Map of Targeted Samples

o	N=23 samples

•	Sequencing

o	n=23 PacBio HiFi

o	n=23 Hi-C

o	n=3 Oxford Nanopore Ultra long reads

•	Assembly

o	Use hifiasm (with --ul) and possibly Verkko for validation of microchromosomes and Sex chromosome.

o	Ragtag for chromosome scaffolding

o	Liftover of chicken annotation to all or just three nanopore assemblies

o	Repeat Masker / EDTA for TE profiling

o	Quast for assembly statistics

o	Minmap cactus or PGGB for comparative analyses

•	Questions:

•	Pangenomic approach?

Objective 3: GRPC vs STGR hybridization

•	Map of Targeted Samples

o	N=29 Sympatric GRPC

o	N=6 Allopatric GRPC (same samples used for long read sequencing)

o	N=20 Sympatric STGR 

o	N=12 Allopatric STGR

o	Additional Samples

•	N=2 putative STGR/GRPC hybrids

•	N=6 GRPC (from Northern DPS)
