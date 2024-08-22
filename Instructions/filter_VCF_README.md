## Variant filtering for the Pacific Islander Biobank v1 dataset

### Quick start:

VCFs are available for each nuclear chromosome (i.e. 1-22, X, and Y).
There are two VCFs for each chromosome, one with only the
category 1 and 2 filters annotated, and one with the category 3
masks applied in addition. Filters have only been annotated,
not applied.

To generate a VCF with category 1 and 2 filters applied:

```
bcftools view -Oz -o [output VCF name] -i 'FILTER=="PASS"' filtered_VCFs/PIBv1_chr[chromosome]_sampledep_and_bed_filtering.vcf.gz
```

To generate a VCF with category 1, 2, and 3 filters/masks applied:

```
bcftools view -Oz -o [output VCF name] -i 'FILTER=="PASS"' filtered_VCFs/PIBv1_chr[chromosome]_allmasks.vcf.gz
```

To generate a VCF with category 1 and 2 filters except
for `ArchaicMasked` applied:

```
bcftools view -Oz -o [output VCF name] -i 'FILTER=="ArchaicMasked"||FILTER=="PASS"' filtered_VCFs/PIBv1_chr[chromosome]_sampledep_and_bed_filtering.vcf.gz
```

To generate a VCF with category 1 and 2 filters except
for `HWE0.0001` and `ExcHet0.0001` applied:

```
bcftools view -Oz -o [output VCF name] -i 'FILTER=="HWE0.0001"||FILTER=="ExcHet0.0001"||FILTER=="HWE0.0001;ExcHet0.0001"||FILTER=="PASS"' filtered_VCFs/PIBv1_chr[chromosome]_sampledep_and_bed_filtering.vcf.gz
```

(Note that you need to include all possible combinations
of the filters to ignore, and the combinations have
filter IDs delimited by semicolons `;`.)

### General terminology:

A filter in the context of VCFs is a label for a given site that
can be used for inclusion or (more commonly) exclusion of a site.

However, filtering can be used more generally to include a threshold
that determines whether or not a genotype should be included
(or alternatively, set to missing). Sometimes this process of
setting genotypes that fail a threshold to missing is called
"masking" a genotype.

Similarly, sometimes a set of annotations indicating to exclude a
site is called a "mask".

Here, we will use two layers of meaning for filter and mask. We
distinguish between the two by the context of whether they describe
an external annotation file, or describe an operation on a VCF.

When describing external annotation files, "filter" should be taken
to mean a set of annotations which indicate *inclusion* of a site
in the output, whereas a "mask" means a set of annotations which
indicate *exclusion* of a site from the output.

When describing operations on a VCF, "filter" should be taken as
a value in the `FILTER` column of the VCF, indicating exclusion
of a site from the output. "Applying a filter" means the process
of excluding such sites based on that value in the `FILTER` column.
"Applying a mask" in the context of VCF operations should be taken
to mean changing specific genotypes to missing based on a threshold.

Some examples may help clarify:

A BED file containing the coordinates of the pseudo-autosomal regions
(PAR) of the sex chromosomes which we want to exclude -> "mask" file

The value `PAR` in the `FILTER` column of the VCF indicating sites
in the VCF that overlap with the PARs -> "filter" in the VCF

Excluding sites with the value `PAR` in the `FILTER` column of the VCF
from the output -> "filtering" the VCF

Setting all genotypes for males in the PAR regions to missing in the
output VCF -> "masking" the VCF

A BED file containing the coordinates of sites in the Altai Neanderthal
genome to retain as reliable -> "filter" file

The value `ArchaicMasked` in the `FILTER` column column of the VCF
indicating sites where not all of the archaic "filter" files
overlap (i.e. not in the intersection of all of these BED files)
and thus to be excluded from the output -> "filter" in the VCF

Etc.

### Types of filters/masks:

We annotated the jointly-genotyped VCFs for this dataset with a variety
of filters and masks, but they generally fall into three categories:

1. BED-type per-site masks

These masks are generally genotype-agnostic mask files that simply
indicate whether or not a site should be considered reliable given
external features. These features include:

- CpG islands (estimated by GC content within a window, plus a dinucleotide enrichment)
- Segmental duplications (aka genomicSuperDups)
- Heng Li's 35-mer mappability mask (aka hs37m_mask35_50)
- PAR (pseudoautosomal regions of the sex chromosomes)
- Archaic hominin minimal filters (Altai, Vindija, and Chagyrskaya Neanderthals plus the Altai Denisovan)

The archaic hominin minimal filters are actually composite filters,
so they somewhat incorporate the other masks, plus some archaic-specific
depth filters.

Each of these mask files is then translated into a single VCF `FILTER`
value:

- `CpG`
- `SegDup`
- `Mappability50`
- `PAR`
- `ArchaicMasked`

Note that the `ArchaicMasked` filter is the genomic complement of the
intersection of all the archaic minimal filter BEDs, so it is quite
stringent. The size of that intersection is about 1.7 Gbp.

2. Sample-dependent per-site filters

These filters rely on data from the samples within the VCF, but
are applied to the entire site. Thus, they represent dataset-dependent
featuress of the variant calls. These features include:

- Hardy-Weinberg Equilibrium
- Excess heterozygosity (similar to HWE, but one-tailed in the direction of excess heterozygotes)
- Proximity of SNPs to indels
- Proximity of indels to indels
- Per-site genotype missingness

Each of these features is thresholded and then translated into a
VCF `FILTER` value:

- `HWE0.0001` (p-value <= 0.0001 for any allele pair)
- `ExcHet0.0001` (p-value <= 0.0001 for any allele pair)
- `SnpGap` (SNP within 5 bp of any called indel)
- `IndelGap` (indel within 5 bp of any called indel)
- `Missingness0.05` (more than 5% of samples have missing genotypes)
- `Missingness0` (any samples have missing genotypes)

There are some additional `FILTER` values that are included in the
VCF and qualify under this category. These are generated by either
GATK GenotypeGVCFs or GATK VQSR, so are more complicated to explain.
Please refer to the GATK documentation for more details. These
`FILTER`s include:

- LowQual (GenotypeGVCFs annotation)
- VQSRTrancheSNP99.90to100.00+ (VQSR annotation, SNP mode)
- VQSRTrancheSNP99.90to100.00 (VQSR annotation, SNP mode)
- VQSRTrancheINDEL99.90to100.00+ (VQSR annotation, INDEL mode)
- VQSRTrancheINDEL99.90to100.00 (VQSR annotation, INDEL mode)
- VQSRTrancheINDEL99.50to99.90 (VQSR annotation, INDEL mode)
- VQSRTrancheINDEL99.00to99.50 (VQSR annotation, INDEL mode)

If you would like to apply VQSR filters, you must exclude sites with
any of these `FILTER` values.

3. Sample-dependent genotype masks

These are thresholds applied at the genotype level, so they are
not visible in the `FILTER` column, they are only noticeable by
comparing the VCF with categories 1 and 2 applied to the VCF with
all three categories applied. The thresholds include:

- Minimum `FMT/DP` >= 10
- Minimum `FMT/GQ` >= 30
- Maximum `FMT/DP` >= a sample-dependent threshold

The sample-dependent maximum `FMT/DP` threshold is based on
optimal k-medians clustering of the 99.5th quantile of autosomal
depth for each sample. k-medians clustering is performed for
all k between 1 and 20, and the best k is chosen by BIC for
a Gaussian mixture model based on the k-medians clustering result.
See the paper for the Ckmeans.1d.df R package for details.

Any samples failing any of these thresholds has its genotype
masked.

~~Note that these category 3 masks are applied after/downstream
of the category 1 and 2 filters, meaning that HWE, ExcHet, and
Missingness are all calculated without accounting for this
masking. We may add a step that recalculates and re-annotates
these filters after applying the category 3 masks in the future.~~

As of 2022/09/21, the HWE, ExcHet, Missingness0, and Missingness0.05
filters are re-calculated in the `*_allmasks.vcf.gz` files *after*
the category 3 masks are applied. So there shouldn't be any issues.
