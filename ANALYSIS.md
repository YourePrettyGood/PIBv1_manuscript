# Analyses and analysis pipelines for the manuscript titled "Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania"

## Analysis summary

- [Variant filtering](#variant-filtering-and-datasets)
  - "Minimal filters"
  - "MAFLD"
  - "Mappability50"
  - "ROH"
- [Geographic distribution of genetic variation](#geographic-distribution-of-genetic-variation)
  - Novel vs. known variation using dbSNP, 1kGP, and gnomAD
  - Region- or population-specific genetic variation
- [Population structure analysis](#population-structure-analysis)
  - Principal components analysis (PCA)
  - ADMIXTURE
  - Outgroup f<sub>3</sub>
  - Neighbour-Joining tree from pairwise F<sub>ST</sub>
  - Local ancestry inference (RFMix)
  - f<sub>4</sub> statistics
- [Demographic inference](#demographic-inference)
  - Runs of Homozygosity (ROH) analysis
  - SMC++
- [Archaic introgression](#archaic-introgression)
  - PCA projection
  - D statistics
  - f<sub>4</sub> ratio statistics
  - Sprime
    - Phased projections
    - Archaic coverage
    - Match rate mixture models
    - Core haplotypes
- [Selection and adaptive introgression](#selection-and-adaptive-introgression)
- [Functional annotation](#functional-annotation)
- [Massively Parallel Reporter Assay](#massively-parallel-reporter-assay)
  - MPRA oligo design
  - MPRA analysis

## Variant filtering and datasets

### "Minimal filters"

While many different sets of filters have been used for human genomic datasets
over the years, we wanted to empirically evaluate the efficacy of these filters
on our compiled dataset, which comprises sequencing data spanning almost a decade.
As primary measures of efficacy, we calculated the rate of Mendelian segregation
errors in 8 trios across the compiled dataset and compared this to the number of
variant calls retained. These serve as approximate measures of the false positive
and total positive rates, so an efficient set of filters minimize the Mendelian
error rate while maximizing the number of variant calls retained.

We examined a wide variety of filters including:
1. GATK's Variant Quality Score Recalibration (VQSR)
2. A set of sample-independent genome-wide masks ("BEDs")
3. A set of sample-dependent per-site criteria ("Sampledep")
4. A set of sample-dependent per-genotype masks ("GTmasks")

The first three categories were annotated in the VCF as values in the `FILTER`
column, while the fourth category required setting genotypes to missing.
We applied "GTmasks" and annotated the `FILTER` column for categories 2 and 3
using the Nextflow pipeline [`filter_VCF.nf`](/Analysis_Pipelines/filter_VCF.nf).
The config used for this pipeline can be found [here](/Configs/Used/filter_VCF/PIBv1_filter_VCF.config).

A much more detailed explanation of the filter annotation process can be found
[here](/Instructions/filter_VCF_README.md).

Evaluation of the criteria was performed using the Nextflow pipeline
[`VCF_stats.nf`](/Analysis_Pipelines/VCF_stats.nf). The config used for the
"minimal filters" dataset can be found [here](/Configs/Used/VCF_stats/PIBv1_VCF_stats_VQSRpass_Missingness0.05_allGTmasks.config).

TL;DR:
The best combination of filters we found was VQSR PASS + GTmasks + per-site
genotype missingness less than or equal to 5%. We refer to these as the
"minimal filters". Phasing was performed on the "minimal filters" dataset.

The "BEDs" included:
- CpG islands from []()
- Large genomic duplications in the reference from UCSC []()
- Heng Li's 35-mer mappability mask []()
- Any site not in the intersection of the accessibility masks from the four high-depth archaic hominin genomes

The "Sampledep" criteria included:
- Excess heterozygosity p-value ≤ 1e-4 using `bcftools +fill-tags -- -t ExcHet`
- Deviation from Hardy-Weinberg Equilibrium (p-value ≤ 1e-4) using `bcftools +fill-tags -- -t HWE`
- Per-site missingness of genotypes less than or equal to 5%
- Per-site missingness of genotypes equal to 0%

The "GTmasks" criteria included:
- Minimum filtered per-sample read depth (`FORMAT/DP`) of 10 reads
- Maximum filtered per-sample read depth (`FORMAT/DP`) of approximately the 99.5th percentile of the per-sample depth distribution
- Minimum PHRED-scaled genotype quality (`FORMAT/GQ`) of 30

Genotypes failing these criteria were masked using `bcftools +setGT`
Because the "GTmasks" masking changes genotypes, the "Sampledep" criteria had to
be recalculated after performing this masking.

### "MAFLD"

After annotating and applying these minimal filters to the VCFs, we generated
a derived "MAFLD" dataset by applying a minimum minor allele frequency (MAF)
threshold of 0.01 (i.e. 1%) using bcftools commit 1eba45c and then LD pruning
with `--indep-pairwise` in 50 SNP windows with 5 SNP step size and an r<super>2</super>
threshold of 0.5 using PLINK 2.00 alpha 3.6. These LD pruning parameters are
consistent with those used in [Choin et al. 2021](https://doi.org/10.1038/s41586-021-03236-5).
These two filters were implemented in the Nextflow pipeline [`MAFLD_VCF.nf`](/Analysis_Pipelines/MAFLD_VCF.nf).

The config used for this pipeline run is available [here](/Configs/Used/MAFLD/PIBv1_noVanuatuContamRelated_globalMAFge0.01_LDpruneChoin.config).

### "Mappability50"

Some downstream analyses (i.e. PCA projection, D statistics,
f<sub>4</sub> ratio statistics) required that our dataset of modern humans
be merged with the high-depth archaic hominins as well as the chimpanzee
genome. The high-depth archaic hominin genomes were masked using Heng Li's
35-mer mappability mask, so we generated a version of our modern human dataset
also including this filter ("Mappability50" in our test of variant filters).
This "Mappability50" dataset consists of the minimal filters VCFs further
filtered to exclude the "Mappability50" `FILTER` column value, followed by
filtering for minimum MAF of 0.01, and then the same LD pruning parameters
as for the above "MAFLD" dataset. All of these tasks were performed using
the [`MAFLD_VCF.nf`](/Analysis_Pipelines/MAFLD_VCF.nf) pipeline.

The config used for this pipeline run is available [here](/Configs/Used/MAFLD/PIBv1_noVanuatuContamRelated_Mappability50_globalMAFge0.01_LDpruneChoin.config).

### "ROH"

Some downstream analyses (i.e. ROH) required slightly different LD pruning, so
we prepared a dataset using the Nextflow pipeline [`MAFLD_VCF.nf`](/Analysis_Pipelines/MAFLD_VCF.nf)
with different LD pruning parameters: `--indep` in 50 SNP windows with 5 SNP step
size and a variance inflation factor (VIF) threshold of 2. These parameters are
recommended for ROH analysis by [Howrigan et al. 2011](https://doi.org/10.1186/1471-2164-12-460)
and were used by [Tucci et al. 2018](https://doi.org/10.1126/science.aar8486).

The config used for this pipeline run is available [here](/Configs/Used/MAFLD/PIBv1_noVanuatuContamRelated_globalMAFge0.01_LDpruneROH.config).

## Geographic distribution of genetic variation

### Novel variants



### Population-specific variants



## Population structure analysis

### Principal components analysis (PCA)

This analysis was performed by Patrick F. Reilly and corroborated using a
different pipeline by Chang Liu.
We ran a principal components analysis (PCA) on the MAFLD dataset using
the `snpgdsPCA()` function of the [SNPRelate](https://doi.org/10.18129/B9.bioc.SNPRelate) R package
after converting from PLINK format to GDS format using `snpgdsBED2GDS()`.
The PCA was computed using the exact algorithm and returning all eigenvectors.
We generated a scree plot of the percent of variance explained by each principal
component, and biplots of the projections of each sample on PCs 1 and 2,
3 and 4, and 5 and 6, respectively. This process was automated in an R script
[`SNPRelate_analyses.R`](/Analysis_Pipelines/SNPRelate_analyses.R).

We ran this script on two different datasets: first on the full worldwide
dataset, and second on the subset consisting only of Taiwan, Island Southeast
Asia, and Oceania ("TaiwanISEAOCN"). The command lines for doing so were as
follows.

Worldwide dataset:

```bash
[path to script]/SNPRelate_analyses.R "pruned_PLINK/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin" \
   "PCA/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_test" \
   PIBv1_noVanuatu_metadata_filtered.tsv \
   '' \
   '' \
   '' \
   "AMR,EUR,AFR,MDE,CSA,EAS,ISEA,OCN" \
   '' \
   '' \
   Population \
   "0,1,2,3,4,0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,0,1,2,3,0,1,2,3,4,5,6,7,8,0,1,2,0,1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,0,1,2,0,1,2,0,1,3,4,2,5,6,7,8,9,10,3,11,4,5,6,7,12,8,9"
```

TaiwanISEAOCN dataset:

```bash
[path to script]/SNPRelate_analyses.R "pruned_PLINK/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin" \
   "PCA/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin" \
   PIBv1_noVanuatu_metadata_filtered.tsv \
   '' \
   Island \
   "Taiwan,Philippines,Flores,New Guinea,New Britain,Mussau,New Hanover,New Ireland,Bougainville,Solomon Islands Taiwan,Philippines,Flores,New Guinea,New Britain,Mussau,New Hanover,New Ireland,Bougainville,Solomon Islands" \
   Island \
   "Taiwan,Philippines,Flores,New Guinea,New Britain,Mussau,New Hanover,New Ireland,Bougainville,Solomon Islands" \
   Population \
   "1,2,8,11,2,0,1,2,0,1,5,3,4,7,6,5,4,0,1,2,8,11,2,6,0,1,11,9"
```

(Yes, I know, those command lines are ugly af... We did a ton of customization
of plot points, labels, and colours, so this was the fastest way I could add
it in without hard-coding anything...)

The biplots were later remade for publication by Audrey Tjahjadi using a
modified version of this R code and adjusted for improved visual contrast
with Adobe Illustrator.

Dependencies used:

- R version 4.1.0
- SNPRelate version 1.28.0
- tidyverse version 1.3.1
- ggrepel
- RColorBrewer
- viridis

### ADMIXTURE

This analysis was performed by Daniela Tejada Martinez and Patrick F. Reilly.
We ran ADMIXTURE with 20 replicates for a range of K from 2 to 12 and integrated
the results across replicates and Ks using CLUMPAK as implemented in
[`ADMIXPIPE.nf`](/Analsis_Pipelines/ADMIXPIPE.nf). [`ADMIXTURE.nf`](/Analysis_Pipelines/ADMIXTURE.nf)
is a drop-in replacement written in DSL2. This analysis was repeated for both
the complete dataset and the subset only including Taiwan, Island Southeast
Asia, and Oceania (`TaiwanISEAOCN`). The configs used for these two runs can
be found [here for the worldwide dataset](/Configs/Used/ADMIXTURE/) and
[here for the TaiwanISEAOCN dataset](/Configs/Used/ADMIXTURE/).

Instructions for running this pipeline are provided [here](/Instructions/ADMIXTURE_README.md).

The MAFLD dataset's PLINK files were used as input for the pipeline, along with
a population map file (a TSV consisting of two columns: Sample ID, and
population label), a file of sample IDs to exclude, and an optional
`drawparams` file serving as a custom configuration for `distruct`.

Dependencies used:

- PLINK version 1.90 beta 7 devel from 2023/01/16
- [ADMIXTURE version 1.3.0](https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz)
- [CLUMPAK version March 26, 2015](https://tau.evolseq.net/clumpak/download/CLUMPAK.zip)
- scripts from YourePrettyGood/HumanPopGenScripts

See the [instructions README](/Instructions/ADMIXTURE_README.md) for detailed
instructions on installing CLUMPAK dependencies, as that was a pain to sort
out, and existing online instructions we found were incorrect or incomplete.

### Outgroup f<sub>3</sub>

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

### F<sub>ST</sub> NJ tree

This analysis was performed by Daniela Tejada Martinez and Patrick F. Reilly.
We calculated pairwise F<sub>ST</sub> for all analysis groups in the dataset
using the estimator of [Hudson, Slatkin, and Maddison 1992](https://doi.org/10.1093/genetics/132.2.583)
as implemented in [PLINK 2.00 alpha 3.6](https://www.cog-genomics.org/plink/2.0/).
We used the "MAFLD" dataset as input. We then constructed a PxP matrix of these
pairwise F<sub>ST</sub> values (where P is the number of populations) to serve
as a distance matrix and used the `nj()` function in the [ape](https://cran.r-project.org/web/packages/ape/index.html)
to construct a neighbour-joining tree from this distance matrix using the
method of Saitou and Nei (1987) MBE. The resulting tree (output in Newick
format using `write.tree()`) was then imported into [iTOL](https://itol.embl.de/)
for visualization.

The command line used for calculating pairwise F<sub>ST</sub> was:
```bash
plink2_linux_avx2_intel --threads 1 --memory 3000 \
  --bfile pruned_PLINK/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin \
  --fst AnalysisGroup method=hudson blocksize=10000 \
  --pheno iid-only PIBv1_Region_AnalysisGroup_phenos.tsv \
  --pheno-name Region,AnalysisGroup \
  --out FST/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_byAnalysisGroup_Hudson_blocksize10000
```

### Local ancestry inference with RFMix

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

### f<sub>4</sub> statistics

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

## Demographic inference

### Runs of Homozygosity (ROH) analysis

Documentation for this step is a work in progress, these analyses were performed
by Audrey Tjahjadi.

### SMC++

This analysis was performed by Patrick F. Reilly.
We inferred single-population trajectories of effective population size over
time (N<sub>e</sub>(t)) from whole-genome sequences of at least 8 samples
per population using [SMC++](https://github.com/popgenmethods/smcpp). For each
of 58 populations, we performed composite likelihood inference based on 10
randomly chosen distinguished lineages per population (or all possible
distinguished lineages for those with less than 10 samples). To distinguish
between long tracts of reference homozygosity and regions of missing data,
we generated a mask of sites deemed uncallable based on read depth. A site
was uncallable if the read depth of any individual in the population was
below 10 or above the 99.5th quantile genome-wide of depth for that individual.
Said a different way, this is equivalent to generating per-sample masks based
on individual read depth thresholds and then taking the union of these masks.
Per-site read depth was calculated using [mosdepth](https://github.com/brentp/mosdepth)
with a minimum mapping quality (MAPQ) threshold of 30. Generation of this mask
as well as the subsequent composite likelihood runs of SMC++ were implemented
in the Nextflow pipeline [`smcpp.nf`](/Analysis_Pipelines/smcpp.nf). The config
used can be found [here](/Configs/Used/SMCpp/PIBv1_smcpp.config).

One note: The config was for a previous version of the `smcpp.nf` pipeline that
simultaneously ran MSMC2 and CHIMP on the same data (though with bugs in the
MSMC2 pipeline). The `smcpp.nf` pipeline found in this repository contains the
same SMC++ steps as the previous pipeline, but the MSMC2 and CHIMP sections
have been refactored out. The `msmc2.nf` pipeline found here is a result of
that refactoring, though with bugs fixed and extended to two-population runs
for cross-coalescence rate analysis as well as MSMC-IM. Neither the MSMC2 nor
CHIMP results are included in the current draft of the manuscript.

Dependencies used:

- mosdepth 0.3.2
- bedtools commit cc714eb
- SMC++ commit 8bdecdf with Python 3.8.6 and GCC 10.2.0
- bcftools commit 1eba45c (htslib commit a1dec95)
- scripts from YourePrettyGood/HumanPopGenScripts

## Archaic introgression

### PCA projection

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

### D statistics

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

### f<sub>4</sub> ratio statistics

Documentation for this step is a work in progress, these analyses were performed
by Chang Liu.

### Sprime



#### Phased projections



#### Archaic coverage



#### Match rate mixture models



#### Core haplotypes



## Selection and adaptive introgression

Documentation for this step is a work in progress, these analyses were performed
by Daniela Tejada Martinez.

## Functional annotation

Documentation for this step is a work in progress, these analyses were performed
by Stephen Rong.

## Massively Parallel Reporter Assay



### MPRA oligo design



### MPRA analysis

Documentation for this step is a work in progress, these analyses were performed
by Stephen Rong.
