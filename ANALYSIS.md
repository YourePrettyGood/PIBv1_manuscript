# Analyses and analysis pipelines for the manuscript titled "Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania"

## Analysis summary

- [Variant filtering](#variant-filtering-and-datasets)
  - [Minimal filters](#minimal-filters)
  - [MAFLD](#mafld)
  - [Mappability50](#mappability50)
  - [ROH](#roh)
- [Geographic distribution of genetic variation](#geographic-distribution-of-genetic-variation)
  - [Novel vs. known variation using dbSNP, 1kGP, and gnomAD](#novel-variants)
  - [Region- or population-specific genetic variation](#population-specific-variants)
- [Population structure analysis](#population-structure-analysis)
  - [Principal components analysis (PCA)](#principal-components-analysis-pca)
  - [ADMIXTURE](#admixture)
  - Outgroup f<sub>3</sub>
  - Neighbour-Joining tree from pairwise F<sub>ST</sub>
  - Local ancestry inference (RFMix)
  - f<sub>4</sub> statistics
- [Demographic inference](#demographic-inference)
  - Runs of Homozygosity (ROH) analysis
  - [SMC++](#smc)
- [Archaic introgression](#archaic-introgression)
  - [PCA projection](#pca-projection)
  - [D statistics](#d-statistics)
  - f<sub>4</sub> ratio statistics
  - [Sprime](#sprime)
    - [Phased projections](#phased-projections)
    - [Archaic coverage](#archaic-coverage)
    - [Match rate mixture models](#match-rate-mixture-models)
    - [Core haplotypes](#core-haplotypes)
    - [Archaic deserts](#archaic-deserts)
- [Selection and adaptive introgression](#selection-and-adaptive-introgression)
- [Functional annotation](#functional-annotation)
- [Massively Parallel Reporter Assay](#massively-parallel-reporter-assay)
  - [MPRA oligo design](#mpra-oligo-design)
  - [MPRA analysis](#mpra-analysis)

Most of these pipelines involve awk scripts from my YourePrettyGood/HumanPopGenScripts
repository, so I omit listing GNU awk as a dependency below. I also don't list
any dependencies for the building/installation of listed dependencies (e.g.
GCC 10.2.0, zlib, etc.).

## Variant filtering and datasets

### Minimal filters

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
- CpG mask from [Vernot et al. 2016](https://doi.org/10.1126/science.aad9416) which cites [Prüfer et al. 2014](https://doi.org/10.1038/nature12886)
- Large genomic duplications in the reference from UCSC [genomicSuperDups.txt.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz)
- Heng Li's 35-mer mappability mask (generated with [SNPable](https://lh3lh3.users.sourceforge.net/snpable.shtml))
- Any site not in the intersection of the accessibility masks from the four high-depth archaic hominin genomes

(Quick note: I couldn't find an original source link for the 35-mer mask, as the
mask linked on the SNPable page is for hs36, not hs37, and all the documentation
we have points to a link on the cdna.eva.mpg.de FTP server, but I can't find it
anywhere on there...)

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

### MAFLD

After annotating and applying these minimal filters to the VCFs, we generated
a derived "MAFLD" dataset by applying a minimum minor allele frequency (MAF)
threshold of 0.01 (i.e. 1%) using bcftools commit 1eba45c and then LD pruning
with `--indep-pairwise` in 50 SNP windows with 5 SNP step size and an r<super>2</super>
threshold of 0.5 using PLINK 2.00 alpha 3.6. These LD pruning parameters are
consistent with those used in [Choin et al. 2021](https://doi.org/10.1038/s41586-021-03236-5).
These two filters were implemented in the Nextflow pipeline [`MAFLD_VCF.nf`](/Analysis_Pipelines/MAFLD_VCF.nf).

The config used for this pipeline run is available [here](/Configs/Used/MAFLD/PIBv1_noVanuatuContamRelated_globalMAFge0.01_LDpruneChoin.config).

### Mappability50

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

The high-depth archaic hominin genome VCFs were originallly in GRCh37 coordinate
space, so we converted the VCFs into hs37d5 coordinate space (and added missing
VCF header lines) using the Nextflow pipeline [`prepare_archaics_dbsnp.nf`](/Processing_Pipelines/gvcf_to_annotated_vcf_nextflow/prepare_archaics_dbsnp.nf).
This pipeline used the scaffold maps from NCBI/GRCh37.p13 coordinates to
UCSC/hg19 coordinates ([GCF_000001405.25_GRCh37.p13_assembly_report.txt](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt))
as well as from UCSC/hg19 coordinates to 1000 Genomes/hs37d5 coordinates
([g1kToUcsc.txt](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/chromAlias/g1kToUcsc.txt)) with a slight modification for the mtDNA contig name
to construct a map of equivalent scaffolds between GRCh37.p13 and hs37d5.
The pipeline also adjusted the REF alleles at a small number of sites that
differ in encoding [between hg19 and hs37d5](https://gatk.broadinstitute.org/hc/en-us/articles/360035890711-GRCh37-hg19-b37-humanG1Kv37-Human-Reference-Discrepancies).

The chimpanzee allele calls are derived from the EPO6 primate alignment from
Ensembl release 71, which includes both the CHIMP2.1.4 and GRCh37 assemblies.
This alignment is provided in EMF format, so we extracted the CHIMP2.1.4
allele as a homozygous diploid genotype with the GRCh37 allele as REF in
VCF format using the `EMFtoVCFheader.awk` script, and then sorted the resulting
VCF with `bcftools sort`. This VCF consisted of sites in regions of the EPO6
alignment where the CHIMP2.1.4 and GRCh37 assemblies had one-to-one alignments.

### ROH

Some downstream analyses (i.e. ROH) required slightly different LD pruning, so
we prepared a dataset using the Nextflow pipeline [`MAFLD_VCF.nf`](/Analysis_Pipelines/MAFLD_VCF.nf)
with different LD pruning parameters: `--indep` in 50 SNP windows with 5 SNP step
size and a variance inflation factor (VIF) threshold of 2. These parameters are
recommended for ROH analysis by [Howrigan et al. 2011](https://doi.org/10.1186/1471-2164-12-460)
and were used by [Tucci et al. 2018](https://doi.org/10.1126/science.aar8486).

The config used for this pipeline run is available [here](/Configs/Used/MAFLD/PIBv1_noVanuatuContamRelated_globalMAFge0.01_LDpruneROH.config).

## Geographic distribution of genetic variation

### Novel variants

This analysis was performed by Patrick F. Reilly.
We counted the number of polymorphic variants found in our minimal filters
dataset ("PIBv1") that matched polymorphic variants found in three major
resources for human polymorphic variants:

- 1000 Genomes Project phase 3 variant calls from 2,504 individuals
- dbSNP build 154
- gnomAD release 2.1.1 (exomes and genomes)

These aren't necessarily the latest versions of each dataset, but are
relatively recent. In the case of the 1000 Genomes Project, the phase 3
variant calls are widely used compared to the NYGC 30x ("phase 4") calls
from [Byrska-Bishop et al. 2022](https://doi.org/10.1016/j.cell.2022.08.004),
although that is starting to change. For dbSNP, this build was selected
in early 2021 while I was developing the variant calling Nextflow pipelines,
and is what was used in the `prepare_archaics_dbsnp.nf` run for this manuscript.
dbSNP build 155 was released around the same time, but I missed out on updating
it. Finally, the gnomAD release is constrained to 2.1.1, as this was the last
release to be in GRCh37 coordinate space.

As noted by [Choin et al. 2021](https://doi.org/10.1038/s41586-021-03236-5),
humans in Oceania harbour a significant amount of genetic variation not
found in existing human whole-genome datasets. Similarly, we found a
substantial proportion of novel genetic variation as compared to the
above-mentioned databases, much of which was polymorphic in the Oceanic
genomes in our dataset.

We assessed this by normalizing the variant representation across our dataset
and the databases using `bcftools norm -m -any -f [ref]` and then matching
variants using `bcftools isec`. Although it's not compiled into a Nextflow
pipeline, these are the commands that were used for this analysis and the
related analysis of region-specific/-private genetic variation:

```bash
#Load the appropriate modules:
module load parallel/20210222-GCCcore-10.2.0
module load bcftools/1eba45c
module load htslib/a1dec95

#Make the subdirectories:
mkdir databases
mkdir noVanuatu_ACrecalc_noGT_norm all_variant_counts
mkdir region_private region_all region_private_counts region_all_counts

#Normalize variant representations for dbSNP 154 in hs37d5 space:
PERL5LIB="" parallel -j16 --eta 'bcftools view -G -r {1} -Ou [path to prepare_archaics_dbsnp.nf run]/final_VCFs/dbSNP_154_hs37d5.vcf.gz | bcftools norm -m -any -f [path to refs]/1kGP/hs37d5/hs37d5.fa -Oz -o databases/dbSNP_154_chr{1}_sites_norm.vcf.gz -; tabix -f databases/dbSNP_154_chr{1}_sites_norm.vcf.gz' ::: {1..22}

#Retain polymorphic variants and normalize variant representations for
# 1000 Genomes Project phase 3 calls:
PERL5LIB="" parallel -j16 --eta 'bcftools view -G -Ou [path to 1kGP resources]/phase3_release20130502_dl20210224/ALL.chr{1}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools norm -m -any -f [path to refs]/1kGP/hs37d5/hs37d5.fa -Oz -o databases/1kGP_phase3_chr{1}_sites_norm.vcf.gz -; tabix -f databases/1kGP_phase3_chr{1}_sites_norm.vcf.gz' ::: {1..22}

#Retain polymorphic variants and normalize variant representations for
# gnomAD 2.1.1 genomes:
PERL5LIB="" parallel -j16 --eta 'bcftools view -G -i "AC_raw>0" -Ou [path to databases]/gnomAD/v2.1.1/gnomad.genomes.r2.1.1.sites.{1}.vcf.bgz | bcftools norm -m -any -f [path to refs]/1kGP/hs37d5/hs37d5.fa -Oz -o databases/gnomAD_genomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz -; tabix -f databases/gnomAD_genomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz' ::: {1..22}

#Retain polymorphic variants and normalize variant representations for
# gnomAD 2.1.1 exomes:
PERL5LIB="" parallel -j16 --eta 'bcftools view -G -i "AC_raw>0" -Ou [path to databases]/gnomAD/v2.1.1/gnomad.exomes.r2.1.1.sites.{1}.vcf.bgz | bcftools norm -m -any -f [path to refs]/1kGP/hs37d5/hs37d5.fa -Oz -o databases/gnomAD_exomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz -; tabix -f databases/gnomAD_exomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz' ::: {1..22}

#Recalculate ACs, remove genotypes, normalize variant representations,
# and remove the now-invalid MAF_* tags from the PIBv1 VCFs:
PERL5LIB="" parallel -j16 --eta 'bcftools +fill-tags -Ou [path to MAFLD_VCF.nf direcory]/MAF_annotated_VCFs/PIBv1_noVanuatu_chr{1}_MAFannotated.vcf.gz -- -t AC,AF -S PIBv1_noVanuatu_samplemap.tsv | bcftools view -G -Ou | bcftools norm -m -any -f [path to refs]/1kGP/hs37d5/hs37d5.fa -Ou | bcftools annotate -x "INFO/MAF_AFR,INFO/MAF_AMR,INFO/MAF_CSA,INFO/MAF_EAS,INFO/MAF_EUR,INFO/MAF_ISEA,INFO/MAF_MDE,INFO/MAF_OCN" -Oz -o noVanuatu_ACrecalc_noGT_norm/PIBv1_noVanuatu_chr{1}_sites_norm.vcf.gz; tabix -f noVanuatu_ACrecalc_noGT_norm/PIBv1_noVanuatu_chr{1}_sites_norm.vcf.gz' ::: {1..22}

#Subset the variants polymorphic in each major geographic region of PIBv1:
PERL5LIB="" parallel -j16 --eta 'inclexpr="MAX(AC_{1})>0";bcftools view -i "${inclexpr}" -Oz -o region_all/PIBv1_noVanuatu_chr{2}_sites_norm_{1}_all.vcf.gz noVanuatu_ACrecalc_noGT_norm/PIBv1_noVanuatu_chr{2}_sites_norm.vcf.gz; tabix -f region_all/PIBv1_noVanuatu_chr{2}_sites_norm_{1}_all.vcf.gz' ::: AFR AMR CSA EAS EUR ISEA MDE OCN ::: {1..22}
#Subset the variants private to each major geographic region of PIBv1:
PERL5LIB="" parallel -j16 --eta 'inclstr="MAX(AC_AFR)=0&&MAX(AC_AMR)=0&&MAX(AC_CSA)=0&&MAX(AC_EAS)=0&&MAX(AC_EUR)=0&&MAX(AC_ISEA)=0&&MAX(AC_MDE)=0&&MAX(AC_OCN)=0";inclexpr=${inclstr//AC_{1})=/AC_{1})>};bcftools view -i "${inclexpr}" -Oz -o region_private/PIBv1_noVanuatu_chr{2}_sites_norm_{1}_private.vcf.gz noVanuatu_ACrecalc_noGT_norm/PIBv1_noVanuatu_chr{2}_sites_norm.vcf.gz; tabix -f region_private/PIBv1_noVanuatu_chr{2}_sites_norm_{1}_private.vcf.gz' ::: AFR AMR CSA EAS EUR ISEA MDE OCN ::: {1..22}

#Match all variants from the PIBv1 dataset to the databases:
PERL5LIB="" parallel -j16 --eta 'bcftools isec -c none -n+1 noVanuatu_ACrecalc_noGT_norm/PIBv1_noVanuatu_chr{1}_sites_norm.vcf.gz databases/dbSNP_154_chr{1}_sites_norm.vcf.gz databases/1kGP_phase3_chr{1}_sites_norm.vcf.gz databases/gnomAD_exomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz databases/gnomAD_genomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz | [path to scripts]/HumanPopGenScripts/PrivateAlleles/isecAlleleCounts.awk -v "chrom={1}" -v "query=ALL" > all_variant_counts/PIBv1_noVanuatu_chr{1}_ALL_isecAlleleCounts.tsv' ::: {1..22}
#Match variants private to major geographic regions of PIBv1 to the databases:
PERL5LIB="" parallel -j16 --eta 'bcftools isec -c none -n+1 region_private/PIBv1_noVanuatu_chr{1}_sites_norm_{2}_private.vcf.gz databases/dbSNP_154_chr{1}_sites_norm.vcf.gz databases/1kGP_phase3_chr{1}_sites_norm.vcf.gz databases/gnomAD_exomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz databases/gnomAD_genomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz | [path to scripts]/HumanPopGenScripts/PrivateAlleles/isecAlleleCounts.awk -v "chrom={1}" -v "query={2}" > region_private_counts/PIBv1_noVanuatu_chr{1}_{2}_isecAlleleCounts.tsv' ::: {1..22} ::: AFR AMR CSA EAS EUR ISEA MDE OCN
#Match variants found in major geographic regions of PIBv1 to the databases:
PERL5LIB="" parallel -j16 --eta 'bcftools isec -c none -n+1 region_all/PIBv1_noVanuatu_chr{1}_sites_norm_{2}_all.vcf.gz databases/dbSNP_154_chr{1}_sites_norm.vcf.gz databases/1kGP_phase3_chr{1}_sites_norm.vcf.gz databases/gnomAD_exomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz databases/gnomAD_genomes_r2.1.1_ACrawgt0_chr{1}_sites_norm.vcf.gz | [path to scripts]/HumanPopGenScripts/PrivateAlleles/isecAlleleCounts.awk -v "chrom={1}" -v "query={2}_all" > region_all_counts/PIBv1_noVanuatu_chr{1}_{2}_isecAlleleCounts.tsv' ::: {1..22} ::: AFR AMR CSA EAS EUR ISEA MDE OCN

#Combine the counts from the above three matching runs into a summary file:
[path to scripts]/HumanPopGenScripts/PrivateAlleles/catIsecCounts.awk -v "numfiles=5" all_variant_counts/PIBv1_noVanuatu_chr{1..22}_ALL_isecAlleleCounts.tsv region_all_counts/PIBv1_noVanuatu_chr{1..22}_*_isecAlleleCounts.tsv region_private_counts/PIBv1_noVanuatu_chr{1..22}_*_isecAlleleCounts.tsv > PIBv1_noVanuatu_combined_isecAlleleCounts.tsv
```

### Population-specific variants

This analysis was performed by Patrick F. Reilly.
The code for this analysis is listed in the [novel variants](#novel-variants)
section. We identified genetic variants in our dataset specific to each major
geographic region by (re-)calculating allele counts using `bcftools +fill-tags -- -t AC,AF`
and filtering the resulting VCF based on an include expression with the `AC_*`
tag for the major geographic region being > 0 and all other `AC_*` tags being
0.

As noted previously in [Bergström et al. 2020](https://doi.org/10.1126/science.aay5012),
humans in Oceania harbour a substantial amount of unique genetic variation. Our
analyses concur with this result, and integrating the results with the novel
variant analysis we find that this remains true even when comparing against
the deep sampling of gnomAD. It is possible that a newer release of gnomAD might
contain some of these putatively Oceania-specific variants, but they would be
at such incredibly low frequency compared to that found in Oceania that the
change in result would be more semantic than real. For instance, such variants
are equally problematic from the standpoint of GWAS, as the power to detect
associations would be much higher in an Oceanic cohort than in current European-
biased GWAS cohorts regardless of whether the frequency in Europeans is truly 0
versus 1 in 100,000. The same logic would apply to the prevalence of any
phenotype associated with these variants.

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

> [!NOTE]
> The config was for a previous version of the `smcpp.nf` pipeline that
> simultaneously ran MSMC2 and CHIMP on the same data (though with bugs in the
> MSMC2 pipeline). The `smcpp.nf` pipeline found in this repository contains the
> same SMC++ steps as the previous pipeline, but the MSMC2 and CHIMP sections
> have been refactored out. The `msmc2.nf` pipeline found here is a result of
> that refactoring, though with bugs fixed and extended to two-population runs
> for cross-coalescence rate analysis as well as MSMC-IM. Neither the MSMC2 nor
> CHIMP results are included in the current draft of the manuscript.

Dependencies used:

- mosdepth version 0.3.2
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

These analyses were performed by Patrick F. Reilly.
We used the reference-free approach implemented in [Sprime](https://github.com/browning-lab/sprime)
to identify genomic regions putatively introgressed from archaic hominins.
"Reference-free" here means that the core algorithm for introgressed tract
detection does not use information about matching to archaic hominin genomes.
We do, however, use the references in a post-processing step to calculate
the match rates of these Sprime tracts to Neanderthals and Denisovans, and
this matching information is used to filter out tracts that are likely false
positives, so only the detection step is truly reference-free. This is similar
to how the [S*](https://github.com/bvernot/freezing-archer) method is also
reference-free during detection, though not during classification of tracts
as Neanderthal or Denisovan (step 3 [here](https://github.com/bvernot/freezing-archer#general-pipeline-details-below)).

Samples from each target population were extracted along with samples from an
outgroup population that is expected to have little to no archaic introgression.
For the manuscript's analyses, we used HGDP Yorubans as an outgroup, as
Yorubans show minimal ([though not quite zero](https://doi.org/10.1016/j.cell.2020.01.012))
signals of archaic introgression, and Yorubans from the 1000 Genomes Project
[were previously used as outgroup with Sprime](https://doi.org/10.1016/j.cell.2018.02.031).

After removing any sites with missing genotypes from this subset VCF and
concatenating the per-chromosome VCFs into an autosomal VCF, we ran Sprime
with the [HapMap phase II genetic map](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip)
with one run per chromosome but inputting the autosomal VCF. Sprime uses
the autosomal variant calls to calculate a normalization for variation in
local mutation rate, hence requiring the autosomal VCF as input despite
processing each chromosome in independent runs. Sprime produces a `.score` file
as output which consists of one line per putative archaic introgressed allele
found in the target population. Each allele is associated with a `SEGMENT`,
which is an ID of a population-level putative introgressed tract (which we
refer to as an "Sprime tract"). These Sprime tracts are effectively a tiling
path across multiple introgressed haplotypes amongst the individuals of the
target population, so the haplotype represented by the full sequence of Sprime
putative introgressed alleles ("Sprime alleles") with the same `SEGMENT` in
the `.score` file might not actually exist as an individual haplotype.

To facilitate downstream analyses, we perform some post-processing of these
`.score` files to provide more information and filtering. These post-processing
steps include:

- Matching Sprime alleles to alleles found amongst the high-depth archaic hominins
- Calculating a "match rate" to Neanderthals and Denisovans based on these allele matches
- Projecting Sprime tracts back onto the genotypes of individuals to identify introgressed regions in each diploid individual
- Calculating an approximate frequency of each Sprime tract in each target population
- Identifying genes that overlap with each Sprime tract

These steps are automated in the Nextflow pipeline [`sprime.nf`](/Analysis_Pipelines/sprime.nf).

There are certain limitations to the projection and tract frequency steps that
we attempt to address with separate pipelines. For the genotypic projection
step, one problem is that the two haplotypes of a diploid individual need not
have introgressed tracts of equal length, which leads to ambiguity, as the
projection must be expressed as an interval in a diploid individual (i.e.
one projection per individual, not per haplotype). This implies some loss
of information in the process, although we attempt to retain some of this
information by annotating the genotypic projection (which is taken as the
longest projected interval that contains homozygous or heterozygous archaic
alleles from a given Sprime tract) with the counts of genotypes as polarized
by the putative archaic alleles. As an example, consider the following Sprime
tract (`ARC`) and modern human genotypes (`TEST`):

```
CHROM:  1   1   1   1   1   1   1   1   1   1   1   1 
POS:   101 105 112 122 146 150 163 181 218 248 257 259
ARC:    A   C   T   C   C   T   C   A   T   G   G   G 
TEST:  C/C G/G A/T A/C C/C T/T C/C A/A C/T T/T C/C T/T
MATCH: M/M M/M M/A M/A A/A A/A A/A A/A A/M M/M M/M M/M
PRE:    P   P   -   -   -   -   -   -   -   -   -   - 
SUF:    -   -   -   -   -   -   -   -   -   S   S   S 
```

(`PRE` and `SUF` are the non-matching prefix and suffix of the projection,
and are ignored from the site counts. `MATCH` indicates whether the alleles
of the genotype match the Sprime allele, though the genotype is not phased,
so the matching is unordered.)

The corresponding genotype projection would be encoded in BED format as:

```
1	111	218	TractID=POP_1;Individual=TEST;State=mixed;HomSprimeSites=4;HetSprimeSites=3;HomModernSites=0	.	.
```

Multiple different combinations of haplotypes could have produced this same
match pattern of genotypes, such as:

```
POS:   101 105 112 122 146 150 163 181 218 248 257 259
H1:     M   M   M   M   A   A   A   A   M   M   M   M 
H2:     M   M   A   A   A   A   A   A   A   M   M   M 
```

or

```
POS:   101 105 112 122 146 150 163 181 218 248 257 259
H1:     M   M   M   M   A   A   A   A   A   M   M   M 
H2:     M   M   A   A   A   A   A   A   M   M   M   M 
```

or

```
POS:   101 105 112 122 146 150 163 181 218 248 257 259
H1:     M   M   A   A   A   A   A   A   M   M   M   M 
H2:     M   M   M   M   A   A   A   A   A   M   M   M 
```

We can't tell the difference between these cases without phased genotypes,
and the sum of haplotype lengths is not equal to two times the genotypic
projection length, so this can affect interpretation if you use the sum of
lengths of genotypic projections as a measure of the total amount of
introgression per individual. Since the total amount of introgression per
individual is frequently a quantity of interest, this is a notable problem.
To address this limitation, we perform haplotype/phased projections with a
[separate pipeline below](#phased-projections).

The other limitation comes from the method for roughly approximating the
frequency of an Sprime tract. The very simple idea we use in this pipeline
is to calculate the frequency of the Sprime allele at each site in the
Sprime tract for a given population, and use the median of these frequencies
as the frequency of the Sprime tract in the given population. In principle,
introgressed tracts should consist of alleles in strong LD, so the frequencies
of Sprime alleles across these sites should be highly correlated, meaning that
the median is a decent way to estimate tract frequency. The problem is that
Sprime tracts are tiling paths, so we are no longer measuring a single
"frequency" for the whole tract, but rather some function of the frequencies
of the multiple distinct haplotypes forming the tiling path. To address this
problem, we follow the approach used by [Gittelman et al. 2016](https://doi.org/10.1016/j.cub.2016.10.041)
and identify a set of "core haplotypes" that consist of Sprime alleles from
the same archaic origin that are in strong LD. This process is performed in a
[separate pipeline below](#core-haplotypes).

> [!IMPORTANT]
> We performed a bit of post-processing to filter and summarize the tracts
> and [core haplotypes](#core-haplotypes) in BED format using the script
> [`PIBv1_Sprime_BEDgeneration_20230309.R`](/Analysis_Pipelines/PIBv1_Sprime_BEDgeneration_20230309.R).
> The resulting BEDs were used as inputs for other analyses including the
> [archaic coverage](#archaic-coverage) analysis below as well as the
> [adaptive introgression](#selection-and-adaptive-introgression) section.

Dependencies used:

- bcftools commit 1eba45c
- htslib commit a1dec95
- Sprime commit 5ff0e15 (i.e. version 20May22.855)
- bedtools commit cc714eb

#### Phased projections

As noted above in the [Sprime section](#Sprime), we address the limitations
of genotypic projections of Sprime tracts by using phased genotypes to perform
haplotype/phased projections. This is performed using the Nextflow pipeline
[`sprime_phased_projection.nf`](/Analysis_Pipelines/sprime_phased_projection.nf).

I use the term "projection" here in a manner analogous to the mathematical
concept of projection, especially visualized like the projection operation in
geometry. The idea is that an Sprime tract is a tiling path, so we take the
set of Sprime alleles that comprise an Sprime tract and map them onto the
sequence of alleles or genotypes (which are just unordered tuples of alleles)
of a given diploid individual or haplotype. It should be rather intuitive
that this mapping is idempotent, as projecting a projection onto the same
individual's genotypes or the same haplotype will result in the same
projection. Technically, though, genotypic projection doesn't quite fit
the definition of a mathematical projection, as the domain and codomain are
not the same set.

The phased projection pipeline is quite simple, and uses the same awk scripts
with slightly different arguments to handle phased genotypes. There are slight
differences in the output BED relative to the genotypic projection BED:

1. The `Individual` tag in column 4 is replaced with a `Haplotype` tag, which is the individual ID followed by `_1` or `_2` depending on the haplotype
2. The `State` tag in column 4 takes on values of `archaic` or `nonarchaic` (or, rarely, `unphased`) rather than `homozygous`, `heterozygous`, or `mixed`
3. Column 4 has tags `ArchaicSprimeSites` and `ModernSprimeSites` instead of `HomSprimeSites`, `HetSprimeSites`, and `HomModernSites` to keep counts of matches along the projection

One important note for the phased projection pipeline:

Results can be influenced in subtle ways by the value of `params.tract_max_gap`.
Genotyping errors and errors in Sprime identifying putative archaic alleles
can lead to projections where putative modern alleles are interspersed among
putative archaic alleles. If everything were perfect, projections would only
consist of a consecutive sequence of putative archaic alleles surrounded by
putative modern alleles. As a heuristic to address these imperfections, we
simply bridge these interspersed "gaps" of putative modern alleles between
segments of putative archaic alleles if the gaps are sufficiently short.
Thus, `params.tract_max_gap` is the maximum number of consecutive modern
alleles allowed without splitting a projection. For the analyses in the
manuscript such as the total amount of introgression (or classified
introgression) per individual and the [archaic coverage analysis](#archaic-coverage)
below, we set this parameter to 0. We originally set this parameter to a
large value (1000000) under the assumption that switches from archaic to
modern and back to archaic ancestry in a short span should be rare, and
had run with both parameter values to compare. Varying the parameter has
a small impact on the [archaic coverage analysis](#archaic-coverage), but
a noticeable impact on per individual introgression values and the pattern
of introgression levels with New Guinean ancestry.

This parameter had a larger effect on some unpublished preliminary analyses
when we were looking at inferring introgression waves and the times of those
waves using MultiWaver3.1 as well as naïve Bayesian estimators on the genetic
map lengths of phased projections. Values of `params.tract_max_gap` larger than
0 produced a strong downward bias on the age estimates of introgression waves
relative to those found from a parallel ArchaicSeeker2.0 + MultiWaver3.1
analysis. However, we aren't certain if this is an indication of the projections
at `params.tract_max_gap = 0` being more accurate or if the MultiWaver3.1
algorithm being biased toward older time estimates by an excess of short
tracts. To a first approximation, MultiWaver3.1 fits a mixture model of
exponentials, but this is a very difficult problem that runs into issues of
identifiability quickly, so we did not proceed further with this analysis.

However, the results of this and a much later comparative analysis with
[hmmix](https://github.com/LauritsSkov/Introgression-detection) indicated
that `params.tract_max_gap=0` should be preferred due to greater consistency
between methods in the inferred relationship between archaic (and especially
Denisovan) introgression and New Guinean ancestry. This does raise an
interesting question about why there is population heterogeneity in the
effect of this max gap parameter (e.g. why are ancestry switches more closely
spaced in some populations than others?)

Dependencies used:

- bcftools commit 1eba45c

#### Archaic coverage

While the amount of archaic introgressed DNA found in any given individual in a
population can be informative, one can also look at the total amount of the
genome that is covered by archaic introgressed DNA across individuals in a
population. In a sense, this can be interpreted as the amount of the archaic
hominin genome that can be [reconstructed](https://doi.org/10.1126/science.aad9416)
from modern human genomes. If you're familiar with genomics, this can also be
understood as the difference between measuring the "depth" of archaic
introgression (in the sense of the average per-individual probability of that
region being of archaic origin) versus the "coverage" of archaic introgression
(in the sense of the amount of the human genome "covered" by archaic
introgressed sequence). The contrast between these two measures and the impact
on the distribution of archaic genetic variants among modern human populations
is nicely laid out in [Witt et al. 2022](https://doi.org/10.1098/rstb.2020.0411).
Comparing these two measures can also lead to inferences about the evolutionary
and population histories of different populations.

Although it can be informative to use a rarefaction/subsampling approach to
calculating archaic coverage (as demonstrated by Witt et al. 2022), we only
calculated archaic coverage for the full sample for each population. We took
the [haplotype/phased projections](#phased-projections) from the previous
section, partitioned them by population or major geographic region, filtered
for only projections from Sprime tracts passing match rate thresholds, and then
took the union of introgressed intervals with `bedtools merge`. The resulting
tiling paths then had lengths summed to come up with the final results in
supplementary table S7 (rounded to the nearest 0.1 Mbp). The code used to
perform this is found below. We also have a secondary estimate using the Sprime
tracts instead of phased projections. These estimates only slightly differ,
as the phased projection-based estimates are slightly conservative.

```bash
module load bedtools/b891a0b

#Summarize archaic coverage by archaic origin for all samples together:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   cat {MDE,EUR,CSA,AMR,EAS,ISEA,OCN}_pops.txt | while read p;
      do
      cat ../Sprime/${p}_chr{1..22}_Sprime_phased_tracts_perSample_maxgap0.bed;
   done | \
      [path to scripts]/HumanPopGenScripts/Sprime/filterTractProjections.awk [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed - | \
      sort -k1,1V -k2,2n -k3,3n | \
      bedtools merge -i - | \
      [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "region=All";
   h="";
done >> PIBv1_total_Sprime_path_length_fromProjections_maxgap0.tsv
#Now by region:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   for r in "MDE" "EUR" "CSA" "AMR" "EAS" "ISEA" "OCN" "nonOCN";
      do
      cat ${r}_pops.txt | while read p;
         do
         cat ../Sprime/${p}_chr{1..22}_Sprime_phased_tracts_perSample_maxgap0.bed;
      done | \
         [path to scripts]/HumanPopGenScripts/Sprime/filterTractProjections.awk [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed - | \
         sort -k1,1V -k2,2n -k3,3n | \
         bedtools merge -i - | \
         [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "region=${r}";
      h="";
   done;
done > PIBv1_perRegion_Sprime_path_length_fromProjections_maxgap0.tsv
#And per population:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   cat ../../PIBv1_Sprime_target_populations.txt | while read p;
      do
      cat ../Sprime/${p}_chr{1..22}_Sprime_phased_tracts_perSample_maxgap0.bed | \
         [path to scripts]/HumanPopGenScripts/Sprime/filterTractProjections.awk [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed - | \
         sort -k1,1V -k2,2n -k3,3n | \
         bedtools merge -i - | \
         [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "pop=${p}";
      h="";
   done;
   h="";
done > PIBv1_perPop_Sprime_path_length_fromProjections_maxgap0.tsv

#If you want the summary based on Sprime tracts instead of phased projections:
#All samples together:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   cat {MDE,EUR,CSA,AMR,EAS,ISEA,OCN}_pops.txt | while read p;
      do
      [path to scripts]/HumanPopGenScripts/Sprime/selectSprimeTracts.awk -v "keycol=SprimePopulation" -v "key=${p}" [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed;
   done | \
      sort -k1,1V -k2,2n -k3,3n | \
      bedtools merge -i - | \
      [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "region=All";
   h="";
done > PIBv1_total_Sprime_path_length_fromTracts.tsv
#Per major geographic region:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   for r in "MDE" "EUR" "CSA" "AMR" "EAS" "ISEA" "OCN" "nonOCN";
      do
      cat ${r}_pops.txt | while read p;
         do
         [path to scripts]/HumanPopGenScripts/Sprime/selectSprimeTracts.awk -v "keycol=SprimePopulation" -v "key=${p}" [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed;
      done | \
         sort -k1,1V -k2,2n -k3,3n | \
         bedtools merge -i - | \
         [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "region=${r}";
      h="";
   done;
done > PIBv1_perRegion_Sprime_path_length_fromTracts.tsv
#Per population:
h="1"
for o in "Ambiguous" "Denisovan" "Neandertal";
   do
   cat ../../PIBv1_Sprime_target_populations.txt | while read p;
      do
      [path to scripts/HumanPopGenScripts/Sprime/selectSprimeTracts.awk -v "keycol=SprimePopulation" -v "key=${p}" [path to final Sprime BEDs]/PIBv1_${o}_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed | \
         sort -k1,1V -k2,2n -k3,3n | \
         bedtools merge -i - | \
         [path to scripts]/HumanPopGenScripts/Sprime/tractSummary.awk -v "header=${h}" -v "origin=${o}" -v "pop=${p}";
      h="";
   done;
   h="";
done > PIBv1_perPop_Sprime_path_length_fromTracts.tsv
```

As a simple metric for detecting reduced archaic coverage relative to
expectation, we also computed a very naïve estimator of expected archaic
coverage from the sum of phased projections across individuals. By analogy
to the simple estimator of sequencing read coverage from sequencing read
depth-of-coverage:

```math
E[C]=1-\exp(-\frac{NL}{G})
```

where $C$ is the coverage, $NL$ is the product of read count and length,
and $G$ is the genome size, we estimate expected archaic coverage as:

```math
E[C_{arc}]=1-\exp(-\frac{\sum_{i,j} L_{i,j}}{G})
```

where $C_{arc}$ is the archaic coverage and $L_{i,j}$ is the length of
projections of tract $j$ in individual $i$ from the population of interest.

Observed archaic coverage will almost certainly be less than this estimate
of expected archaic coverage, but the degree to which archaic coverage is
reduced depends on a variety of factors including demographic factors like
population bottlenecks and even the efficacy of purifying selection on
archaic introgressed sequence. We observed that the five Oceanic populations
with the strongest bottleneck signals had significantly reduced archaic
coverage relative to expectation.

#### Match rate mixture models

Based on the Denisovan and Neanderthal match rates provided by the
postprocessing steps of the [`sprime.nf`](#Sprime) pipeline, we fit Gaussian
and Beta mixture models to the match rate distributions with varying
numbers of components and selected the best fitting model (i.e. the optimal
number of mixture components along with their parameter estimates) for each
population. We provide an R script for automating this process for Denisovan
match rates ([`SprimeGMMBMM.R`](/Analysis_Pipelines/SprimeGMMBMM.R)) as a
reference implementation, though the final plots and Neanderthal match rate
mixture model fits included in the supplement were generated with a different
though similar script.

Mixture models with between 1 and 10 components were fit, and the optimal model
was chosen using BIC. For Gaussian mixture models, we were able to retain all
fit models, but for Beta mixture models the `betareg` interface only allowed
us to retain the best fit model. We defined Denisovan tracts for the mixture
models as those with at least 30 Sprime sites with non-missing alleles in both
archaic hominin species, Neanderthal match rate below 30%, and Denisovan match
rate above 30%. Due to limitations in the betareg model, we encoded match rates
of 100% as 99.99999% (i.e. 1.0 as 0.9999999). Though we used the logit link
function for `betamix()`, this limitation persisted regardless of choice of
link function.

Previous work ([Browning et al. 2018](https://doi.org/10.1016/j.cell.2018.02.031),
[Jacobs et al. 2019](https://doi.org/10.1016/j.cell.2019.02.035),
[Choin et al. 2021](https://doi.org/10.1038/s41586-021-03236-5))
used Gaussian mixture models to identify these multiple Denisovan-like
introgression waves, though limited their analyses to a maximum of 2 components
per population. Our initial analyses thus performed Gaussian mixture model
fits, though it quickly became apparent that this approach (in addition to the
use of likelihood ratio tests or the Akaike Information Criterion (AIC))
resulted in poor model fits and overestimates of the optimal number of
components. We are grateful to John Novembre for suggesting the use of
Beta mixture models for this analysis.

Dependencies used:

- R 4.1.0
- flexmix 2.3-18
- betareg 3.1-4
- tidyverse 1.3.1
- doParallel
- foreach

#### Core haplotypes

As mentioned in the [Sprime](#Sprime) section earlier, in principle
introgressed tracts should consist of alleles in strong LD, so the frequencies
of Sprime alleles across these sites should be highly correlated. However,
Sprime tracts are also tiling paths of introgressed haplotypes. Thus, we
needed a way to identify segments of these Sprime tracts that are in strong
LD. To address this problem, we followed the approach used by [Gittelman et al. 2016](https://doi.org/10.1016/j.cub.2016.10.041)
to identify a set of "core haplotypes" that consist of Sprime alleles from
the same archaic origin that are in strong LD.

First, we classified the full set of Sprime sites from a population based
on the Denisovan and Neanderthal match rates of the Sprime tracts containing
these sites. Sprime tracts were first filtered for a minimum of 30 Sprime
sites with non-missing alleles in both archaic hominin species. Tracts also
needed to have Neanderthal or Denisovan match rate (or both) greater than or
equal to 30% (i.e. 0.3). Tracts passing these criteria were then classified
into "Neandertal", "Denisovan", or "Ambiguous" based on:

- "Neandertal" if Neanderthal match rate >= 0.3 and Denisovan match rate <= 0.3
- "Denisovan" if Neanderthal match rate <= 0.3 and Denisovan match rate >= 0.3
- "Ambiguous" if Neanderthal match rate >= 0.3 and Denisovan match rate >= 0.3

We evaluate each archaic origin independently, and further filter sites based
on whether the Sprime allele directly matches the relevant archaic hominin
species (or either of the two species, in the case of "Ambiguous" origin).
This extra filter is enabled by the `only_matches=1` parameter of
`extract_Sprime_arcmatch_sites.awk`, and this is the case in the Nextflow
pipeline [`adaptiveintrogression.nf`](/Analysis_Pipelines/adaptiveintrogression.nf).

This step generates a list of Sprime sites from each archaic origin, roughly
equivalent to the "tag SNPs" used in the Gittelman et al. 2016 approach.
However, any correspondence between an Sprime site in this list and it's
original Sprime tract is lost. We do expect Sprime sites from the same tract
to be in strong LD though.

Thus, the next step is to calculate pairwise LD amongst all Sprime sites of a
given archaic origin on a particular chromosome. In the `adaptiveintrogression.nf`
pipeline, we do this with [vcftools](https://github.com/vcftools/vcftools),
calculating r<super>2</super> from either unphased genotypes (`params.phased = false`)
or phased genotypes (`params.phased = true`), which triggers either the
`--geno-r2` or `--hap-r2` flag to `vcftools`. Only pairs with r<super>2</super> >= `params.min_r2` (default=0.3)
are retained in the output, as these sites are in sufficiently strong LD to
be considered candidates for a core haplotype.

As a slight digression, r<super>2</super> can be calculated at a different
sample size scale than the original Sprime target group if desired. For example,
let's say the Sprime target group was HGDP Dai Chinese where the sample size
is n=8. Calculating r<super>2</super> from 8 diploid individuals will produce
a fairly high variance estimate of r<super>2</super>, which isn't great if you
want to identify sites in reliably strong LD for a core haplotype. Instead, you
could calculate LD amongst those sites in East Asians where your sample size
would then be in the 100s and thus much lower variance. This, of course, assumes
that there aren't big differences in LD among the constituent populations, but
if that seems like a reasonable assumption, then it should work. We chose to
follow this approach in the manuscript, as Sprime was run with target groups
being analysis groups/populations, so sample sizes were well below 100. If you
were to run this pipeline on the 1000 Genomes Project dataset, you could just
calculate r<super>2</super> on the target populations themselves, since those
sample sizes are in the hundreds per population.

Once we have the pairs of Sprime sites with r<super>2</super> >= 0.3, we
construct a simple undirected graph where Sprime sites are nodes and site pairs
with r<super>2</super> >= 0.3 are edges. Core haplotypes are then identified
as the connected components of this graph, and the R script outputs sites
labeled by the connected component in which they are found (a simple integer
label). This simple label is then translated into a core haplotype ID by
prepending with the Sprime target group, chromosome number, and archaic origin.

We then perform some post-processing of these core haplotypes:

1. The core haplotypes are mapped back onto the original Sprime outputs to annotate with source tract and archaic hominin match information
2. We perform (genotypic) projection of the core haplotypes onto the genotypes of individuals
3. We estimate the frequency of each core haplotype by taking the median frequency amongst the core haplotype Sprime sites
4. We identify genes overlapping each core haplotype
5. We construct a map file that indicates the correspondence between Sprime tract IDs and core haplotype IDs

The frequency estimate here is probably the most useful, as we have now assured
that the sites in a core haplotype are in strong LD, thus the allele frequencies
should be strongly correlated and not from a tiling path of haplotypes. One
obvious use is to filter core haplotypes by frequency to identify high-frequency
introgressed genomic regions in each target group. These high-frequency core
haplotypes are decent candidates for putative sites of adaptive introgression.
In the manuscript supplement, we show that these high-frequency core haplotypes
reproduce previously-identified candidates of adaptive introgression in
Europeans and East Asians (notably, many of these previous studies used 1000
Genomes Project data, so we are able to reproduce these signals from a different
dataset using a different approach).

Dependencies used:

- bcftools commit 1eba45c
- vcftools 0.1.16
- R 4.1.0
- tidyverse 1.3.1
- igraph 1.3.5
- bedtools commit cc714eb

#### Archaic deserts

Documentation for this step is a work in progress, most of these analyses were
performed by Stephen Rong.
Please see Stephen's [MPRA and functional analyses Github repository](https://github.com/stephenrong/PIBv1_MPRA_analyses)
for details and scripts, or see [the git submodule here](/Analysis_Pipelines/PIBv1_MPRA_analyses/).

Identification of archaic deserts was performed using a combination of
BEDtools commands and awk scripts. Specifically, we first merged all phased
projections from all non-African individuals with `bedtools merge`, then
took the autosome-wide complement using `bedtools complement` with the `-g`
flag pointing to a `.genome` file of only the autosomes of hs37d5. We then
subtracted out the masked regions of hs37d5 with `bedtools subtract`, and
finally thresholded the resulting intervals absent of introgression to
extract only those intervals of at least 1 Mbp using a custom awk script.
Precise commands used can be found below:

```bash
#Merge all phased projections to get the set of introgressed regions:
for m in "0" "1000000";
   do
   while read p;
      do
      cat projections/${p}_Sprime_{Ambiguous,Denisovan,Neandertal}_phased_tracts_perSample_maxgap${m}.bed | \
         egrep -v "^#";
   done < <(cat ../PIBv1_Sprime/phased_projections/tiling_paths/{AMR,CSA,EAS,EUR,ISEA,MDE,OCN}_pops.txt) | \
      sort -k1,1V -k2,2n -k3,3n | \
      bedtools merge -i - > projections/PIBv1_total_tiling_paths_maxgap${m}.bed;
done
#Take the autosome-wide complement and subtract the masked regions:
for r in "hg19" "hs37d5";
   do
   for m in "0" "1000000";
      do
      bedtools complement -i projections/PIBv1_total_tiling_paths_maxgap${m}.bed -g <(head -n22 [path to ref]/${r}.genome) | \
         bedtools subtract -a - -b [path to ref]/${r}_assembly_gaps.bed > projections/PIBv1_total_tiling_paths_complement_${r}unmasked_maxgap${m}.bed;
   done
done

#Now extract only no-introgression intervals of at least 1 Mbp:
#Ran this using both hg19 and hs37d5 masked regions and both maxgap
# values to compare. There are no differences whatsoever in the results.
for r in "hg19" "hs37d5";
   do
   for m in "0" "1000000";
      do
      awk 'BEGIN{FS="\t";OFS=FS;}{if ($3-$2 >= 1000000) {print;};}' projections/PIBv1_total_tiling_paths_complement_${r}unmasked_maxgap${m}.bed > deserts/PIBv1_archaic_deserts_${r}unmasked_maxgap${m}.bed;
   done
done
```

Dependencies used:

- bedtools commit b891a0b

## Selection and adaptive introgression

Documentation for this step is a work in progress, these analyses were performed
by Daniela Tejada Martinez.

In response to a reviewer comment, we compared the adaptive introgression
candidate regions identified using our method with windows identified as
adaptive introgression candidates using the U and Q95 statistics of
[Racimo et al. 2017 MBE](https://doi.org/10.1093/molbev/msw216). Since we
weren't able to find a reference implementation of these statistics, we
re-implemented them as a Nextflow pipeline [RacimoUQ95.nf](/Analysis_Pipelines/RacimoUQ95.nf).
In particular, we evaluated `U_{A,B,C,D}` and `Q95_{A,B,C,D}` for `A` as all
HGDP African populations, `B` being one of 18 target Oceanic groups (we ran
the Sepik and Goroka of New Guinea separately for simplicity), `C` being the
Vindija Neanderthal, and `D` being the Altai Denisovan. We ran each group
three times for the different combinations of `y` and `z` corresponding to
different archaic origins: `(y,z)=(1,0)` for Neanderthal, `(y,z)=(0,1)` for
Denisovan, and `(y,z)=(1,1)` for Ambiguous. We evaluated U and Q95 with
`w=0.01` and U with `x=0.5`, following the main analyses of Racimo et al. 2017.
We tested this pipeline using the 1000 Genomes Project phase 3 variant calls
along with the Vindija33.19 and Denisova-Phalanx calls from the merged
archaics callset used elsewhere in this manuscript, and found hits for EUR
and EAS in many of the same regions (e.g. BNC2, CHMP1A, POU2F3). All runs
used non-overlapping windows of 40 kbp. For our own dataset, we identified
windows with U and Q95 greater than the 95th percentile genome-wide of each
statistic, and checked these windows for overlap with adaptive introgression
candidate regions from our method, finding an average of 70.1% (s.d. 11.7%)
of AI candidate regions from our methods per population overlapping with a
significant U and Q95 window per population, whereas only an average of
24.8% (s.d. 7.3%) of significant U and Q95 windows per population
overlapping an AI candidate region from our method. This indicates that our
method is both accurate and conservative.

## Functional annotation

Documentation for this step is a work in progress, these analyses were performed
by Stephen Rong.
Please see Stephen's [MPRA and functional analyses Github repository](https://github.com/stephenrong/PIBv1_MPRA_analyses)
for details and scripts, or see [the git submodule here](/Analysis_Pipelines/PIBv1_MPRA_analyses/).

## Massively Parallel Reporter Assay

The massively parallel reporter assay (MPRA) portion of the manuscript was
performed in collaboration with the lab of Steven K. Reilly in the Department
of Genetics at Yale School of Medicine. The Reilly lab also performed the
functional analyses described in [the previous section](#functional-annotation).

### MPRA oligo design

MPRA oligo design was performed by Samantha L. Miller with code written by
Patrick F. Reilly and design guidance by Steven K. Reilly.
These oligos were designed to test the effects of single archaic introgressed
variants identified from a pilot dataset of 92 samples from 5 Near Oceanic
populations. Variants from Sprime tracts with median frequency estimates of
at least 30% (i.e. 0.3) in any of the 5 Near Oceanic populations were initially
selected before further filtration. A much more detailed explanation of this
process and the filters considered can be found [here](/Analysis_Pipelines/MPRA/README.md).

Variants were only retained if the Sprime allele was an ALT allele (it is
possible that Sprime emits a REF allele as the putative archaic allele in a
tract), and oligos were designed based on the hs37d5 reference sequence
including 85 bp of 5' flanking sequence and 84 bp of 3' flanking sequence.
Two oligos were designed for each variant: one with the REF allele and one with
the putative archaic ALT allele. 15 bp of adapter sequence was then added to
each end of the target oligo sequence, and these oligos (along with control
oligos) were submitted to Twist Biosciences for synthesis.

### MPRA analysis

Documentation for this step is a work in progress, these analyses were performed
by Stephen Rong.
Please see Stephen's [MPRA and functional analyses Github repository](https://github.com/stephenrong/PIBv1_MPRA_analyses)
for details and scripts, or see [the git submodule here](/Analysis_Pipelines/PIBv1_MPRA_analyses/).
