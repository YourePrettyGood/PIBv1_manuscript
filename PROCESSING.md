# Data processing pipelines for the manuscript titled "Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania"

## Data processing summary

- [Sample screening](#sample-screening)
- [Mapping and per-sample processing](#per-sample-processing)
- [Joint genotyping](#joint-genotyping)
- [Phasing](#phasing)
- [Dataset QC](#dataset-qc)

## Sample screening

After baseline screening of samples based on DNA concentration and volume,
initial PCR or PCR-free libraries were prepared. We performed further screening
of these libraries by sequencing ~1% of a NovaSeq 6000 S4 lane per sample,
resulting in between 1 and 4x depth per sample. This low-depth data was used
for a variety of QC and screening steps compiled in a `bash` pipeline:
[low depth QC pipeline](/Processing_Pipelines/low_depth_QC)

These QC steps include:

- QC of the sequencing run and libraries with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Adapter trimming, read pair merging, and insert size distribution with [SeqPrep](https://github.com/jstjohn/SeqPrep)
- Contamination screening with [FastQ-Screen](https://github.com/StevenWingett/FastQ-Screen) and [kraken2](https://github.com/DerrickWood/kraken2)+[Bracken](https://github.com/jenniferlu717/Bracken) using the [PlusPFP](https://benlangmead.github.io/aws-indexes/k2) database
- Read mapping statistics with [BWA-MEM](https://github.com/lh3/bwa) and [SAMtools](https://github.com/samtools/samtools)
- mtDNA haplogroup classification with [mutserve](https://github.com/seppinho/mutserve)+[haplogrep](https://github.com/seppinho/haplogrep-cmd)
- Y haplogroup classification with [Yleaf](https://github.com/genid/Yleaf)
- Genetic sex estimation based on read mapping rates and read depth including [Skoglund et al.'s R_Y](https://doi.org/10.1016/j.jas.2013.07.004)

The results of the latter three steps were compared with previous results
on mtDNA haplogroups, Y haplogroups, and self-reported sex.

This pipeline was developed very early on, so I never got around to rewriting
it as a Nextflow pipeline.

Furthermore, we performed multi-sample joint QC at multiple stages (92 samples,
154 samples, 170 samples, 184 samples) using [ANGSD](https://github.com/ANGSD/angsd).
We ran ANGSD with arguments
`-P ${cores} -doMaf 1 -doMajorMinor 1 -doGlf 2 -GL 2 -SNP_pval 1e-6`
followed by [PCAngsd](https://github.com/Rosemeis/pcangsd) and
[ngsRelate v2](https://github.com/ANGSD/NgsRelate).

## Per-sample processing

Once samples were sequenced to full depth, we processed them and published WGS
FASTQs from [Vernot et al. 2016](https://doi.org/10.1126/science.aad9416),
[Tucci et al. 2018](https://doi.org/10.1126/science.aar8486),
[Bergström et al. 2020](https://doi.org/10.1126/science.aay5012), and
[Choin et al. 2021](https://doi.org/10.1038/s41586-021-03236-5) from FASTQ
to gVCF following the GATK Best Practices [Van der Auwera et al. 2013](https://10.1002/0471250953.bi1110s43)
with some modifications.

1. Adapters were trimmed from FASTQs using [AdapterRemoval v2](https://github.com/MikkelSchubert/adapterremoval)
2. Adapter-trimmed reads were mapped to the human reference genome [hs37d5](https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/) using [BWA-MEM](https://github.com/lh3/bwa)
3. Alignments were sorted and BAMs from the same sample were merged using [SAMtools](https://github.com/samtools/samtools)
4. Per-sample BAMs were scattered into 16 groups of intervals of approximately
200 Mbp each, only ever splitting at gaps in the reference
5. Scattered BAMs had PCR and optical duplicates marked using [Picard](https://broadinstitute.github.io/picard/) MarkDuplicates,
realigned around indels using GATK `IndelRealigner` v3.8-1-0-gf15c1c3ef, and
base quality score recalibration (BQSR) statistics were collected using GATK
`BaseRecalibrator` v4.1.8.1
6. BAMs were gathered back into one BAM per sample with SAMtools, and BQSR
statistics were also gathered using GATK `GatherBQSRReports`
7. BQSR was applied to these per-sample BAMs using GATK `ApplyBQSR` v4.1.8.1 to
produce a single "analysis-ready" BAM per sample
8. Due to a small issue with regions splitting at a single base gap on chr6,
these "analysis-ready" BAMs were then reprocessed with [dedup_merged_BAMs](https://github.com/YourePrettyGood/dedup_merged_BAMs),
which deduplicated alignments that spanned the gap
9. Call variants for each sample using GATK `HaplotypeCaller` across 29 groups
of intervals of approximately 100 Mbp each, splitting only at gaps in the ref
as previously, and then these scattered gVCFs are merged with Picard `MergeVcfs`

Steps 1 through 7 and 9 were performed with the Nextflow pipeline [`fastq_to_gvcf.nf`](/Processing_Pipelines/fastq_to_gvcf_nextflow/fastq_to_gvcf.nf),
while steps 8 and 9 were performed with the Nextflow pipeline [`recall_bams.nf`](/Processing_Pipelines/fastq_to_gvcf_nextflow/recall_bams.nf).
The latest version of `fastq_to_gvcf.nf` directly incorporates step 8, but the
results should be equivalent. Step 9 is not affected directly by the gap
spanning issue.

Two extra options of note were set:

1. During step 5, the `OPTICAL_DUPLICATE_PIXEL_DISTANCE` option of Picard
`MarkDuplicates` was set to `2500` for all samples except those from SGDP
(i.e. all samples except those from unpatterned flowcells), whereas SGDP
libraries were set to `100`, as recommended by the `MarkDuplicates` documentation
2. During step 9, the `--pcr-indel-model` argument was set to `NONE` for
PCR-free libraries, and was left at default value for PCR libraries

This pipeline also has QC steps:

- Adapter trimming statistics from AdapterRemoval
- mtDNA haplogroup estimation with mutserve and haplogrep
- Y haplogroup estimation with Yleaf
- Genetic sex estimation and genome-wide read depth distribution
- BAM metrics with Picard `CollectWgsMetrics` and `CollectGcBiasMetrics`
- The final gVCF is validated using GATK `ValidateVariants` with dbSNP 138 and `--validation-type-to-exclude ALLELES`

## Joint genotyping

Once all per-sample gVCFs were prepared, we performed joint genotyping,
annotated known variant IDs with dbSNP 154, and performed variant quality score
recalibration (VQSR) as an initial filter. These steps were performed using the
Nextflow pipeline [`gvcf_to_annotated_vcf.nf`](/Processing_Pipelines/gvcf_to_annotated_vcf_nextflow/gvcfs_to_annotated_vcf.nf).
The joint genotyping and VQSR steps generally follow the GATK Best Practices
workflow for joint genotyping of germline variants (Joint Genotyping WDL version
1.6.3) with slight modifications.

This involved a few main steps:

1. Import of gVCFs using GATK `GenomicsDBImport` in 55 different groups of intervals of approximately 50 Mbp each across the major nuclear chromosomes
2. Joint genotyping using GATK `GenotypeGVCFs` v4.1.8.1 for each group of genomic intervals
3. Merging of scattered jointly-genotyped VCFs and VQSR on SNPs followed by INDELs, as well as recovery of `FILTER` and `INFO` annotations lost by VQSR
4. Splitting of the single VQSR-annotated VCF into per-chromosome VCFs using `bcftools +scatter`
5. Annotation of `ID` column with dbSNP 154 rsids using `bcftools annotate`

The dbSNP 154 VCF used was in GRCh37 coordinate space, and the high-depth
archaic hominin VCFs were in hg19 coordinate space, so both were mapped over to
hs37d5 space using the Nextflow pipeline [`prepare_archaics_dbsnp.nf`](/Processing_Pipelines/gvcf_to_annotated_vcf_nextflow/prepare_archaics_dbsnp.nf).

Although the `gvcfs_to_annotated_vcf.nf` pipeline is capable of merging the
genotypes of the high-depth archaic hominins into the final VCF, we opted to
skip this step for the main VCFs and perform merging as necessary later for
downstream analyses.

## Phasing

The unphased "minimal filters" per-chromosome VCFs generated by filtering the
output of [Joint genotyping](#joint-genotyping) above (see
[this page](/ANALYSIS.md#) for details on that filtering) was then used along
with the BAMs in order to phase the VCFs using a read-informed statistical
phasing approach. Specifically, we performed initial read-backed phasing with
[WhatsHap](https://github.com/whatshap/whatshap) 1.6 and then used these
read-backed phased genotypes to inform statistical phasing with
[SHAPEIT 4.2.2](https://odelaneau.github.io/shapeit4/). This phasing was
performed with the Nextflow pipeline [`phase_final_vcf.nf`](/Processing_Pipelines/phase_final_vcf_nextflow/phase_final_vcf.nf).

Due to differences in BAM validation between GATK `HaplotypeCaller` and WhatsHap,
we performed an extra step of alignment deduplication on the BAMs for chromosome
6 to handle supplementary alignments with different CIGAR strings. This
deduplication is handled by the Nextflow pipeline [`dedup_bams_bychrom.nf`](/Processing_Pipelines/phase_final_vcf_nextflow/dedup_bams_bychrom.nf).

For computational efficiency, the read-backed phasing was performed on VCFs
consisting of all samples from a single analysis group, which was generally
tens of samples. This subsetting does not affect the results of read-backed
phasing, as this is done serially per sample. We included the `--indels`
argument to WhatsHap to ensure that indels were included for phasing.

Read-backed phased VCFs were then merged together and passed into SHAPEIT 4.2.2
with the b37 genetic map provided by SHAPEIT4, the settings preset for phasing
of whole genome sequencing data (`--sequencing`), an argument for incorporating
the WhatsHap phased genotypes as a prior with error rate `--use-PS 0.0001`, and
arguments recommended for higher accuracy (`--mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m`
and `--pbwt-depth 8`). Each chromosome was phased separately using 16 threads.

During this process SHAPEIT4 drops all multiallelic sites as well as all site
and genotype metadata (e.g. `INFO`, `FILTER`, and `FORMAT` tags), so we
transferred the phased genotypes back into the original unphased VCFs
using `bcftools annotate -c FORMAT/GT --pair-logic exact --single-overlaps`.

Initial QC of this phasing was performed using `whatshap stats` for overall
phasing and phase block statistic as well as `bcftools +trio-stats` for
trio-based phasing switch error rate estimates. We note that the switch
error rate calculated by `bcftools +trio-stats` on the final phased VCF
is calculated assuming that all sites along a chromosome are phased, so
artifactual switch errors may be counted at unphased positions such as
multiallelic sites excluded by SHAPEIT4.

We also performed non-pipelined checks of phasing quality using a select
subset of 8 HGDP individuals with 10x Genomics linked read data generated
by [Bergström et al. 2020](https://doi.org/10.1126/science.aay5012). Since
these 10x Genomics VCFs were generated based on GRCh38, we lifted these VCFs
over to hs37d5 using Picard `LiftoverVcf` using the hg38 to hg19 chain file
from UCSC. Phasing switch error was calculated using `whatshap compare` with
the 10x Genomics VCF as ground truth. We further validated that the liftover
process did not significantly affect error rate estimates by comparing the
switch error rates of the Bergström et al. 2020 statistically phased VCFs
for these same samples both in their original GRCh38 space and lifted over
to hs37d5 space. Switch error rate estimates were quite similar between
the lifted and unlifted call sets.

## Dataset QC

Pipelines and code used for QC of the full dataset are generally hosted
in the [Analysis_Pipelines](/Analysis_Pipelines) folder. However, we describe
them here as well.

Beyond the initial sample screening described above, we also screened samples
after variant calling. This included:

- Checking for contamination by other human samples with [VerifyBamID2](https://github.com/Griffan/VerifyBamID) using the 100k SNP HGDP panel and the FREEMIX model
- Heteroplasmy checks using [haplocheck](https://github.com/genepi/haplocheck)
- Relatedness checks using [somalier](https://github.com/brentp/somalier)
- Checking for heterozygosity outliers using `bcftools +smpl-stats`
- Calculating the KING and IBS0 coefficients for relatedness checks using PLINK 2.00 alpha 3.6 and cross-checking with [SNPRelate](https://doi.org/10.18129/B9.bioc.SNPRelate) v1.28.0

The first three steps are implemented in a Nextflow pipeline [`BAM_contam_QC.nf`](/Analysis_Pipelines/BAM_contam_QC.nf).
Heterozygosity was calculated in the Nextflow pipeline [`VCF_stats.nf`](/Analysis_Pipelines/VCF_stats.nf),
and outliers were identified using an R script.
The PLINK and SNPRelate relatedness analyses were performed without a Nextflow
pipeline. In particular, this QC analysis required the generation of relatedness
graphs using an R script as well as manual curation of these graphs to identify
the optimal nodes to exclude in order to minimize connectivity in the graph.
Relatedness pruning algorithms such as those implemented in KING and PLINK can
produce suboptimal results in this respect, hence the manual curation.
