# Selection of variants for the PIB92 pilot MPRA

## Initial variant set

Putative archaic variants were identified from a jointly genotyped and VQSR-
filtered variant set comprised of 92 samples from the Pacific Islander Biobank,
279 worldwide samples from the SGDP, and other samples not used for analyses.
These samples were jointly genotyped on all 22 autosomes.

Each of the 5 Bismarck Archipelago populations sufficiently represented in
the 92 PIB samples was run through Sprime (Browning et al. 2018 Cell) with
the SGDP Africans as outgroup. For each resulting Sprime tract, we calculated
the match rate to the Altai Denisovan and to the set of alleles from the
Altai, Chagyrskaya, and Vindija Neandertals, as well as the frequency of each
putative archaic allele in each of the 5 PIB populations. Sprime tracts were
selected if there were at least 30 sites with non-missing data in the Altai
Denisovan as well as at least 30 sites with at least one non-missing genotype
amongst the three Neandertals, had greater than 30% match rate to the Altai
Denisovan or the allele sets combined across all three Neandertals, and if
the median archaic allele frequency along the tract was at least 30% in at
least one Bismarck Archipelago population.

The first two of these filters ensure that the tract identified by Sprime
is likely to be of archaic hominin origin (and not a false positive), whereas
the third filter serves to enrich for candidates of adaptive introgression.

This initial variant set consisted of 116,324 putatively archaic alleles.

## Variant filtration based on genomic features

To further narrow the list of candidates and minimize the risk of false
positive variant calls and poor quality oligo design, we filtered the above
variant set by excluding any variants that overlapped with the following
masks, and including any variants that overlapped with the following filters:

### Masks

1. UCSC CpG islands for hg19 (cpgIslandExt)
1. UCSC genomic segmental duplications for hg19 (genomicSuperDups)
1. UCSC simple repeats/low complexity sequence mask for hg19 (hg19.TRF.bed.gz)
1. Heng Li's 35-mer mappability mask (hs37m_mask35_50pct) generated with SNPable

### Filters

1. ENCODE version 4 candidate *cis* regulatory elements (cCREs) lifted over from GRCh38 to hs37d5 with UCSC liftOver

### Mask/Filter preprocessing

Any masks/filters intended for hg19 were preprocessed to be compatible with
the hs37d5 reference genome by removing any "chr" prefixed to scaffold names,
and dropping any records with scaffolds that contained an underscore (i.e.
all the unplaced scaffolds and patches in hg19). Thus, only the autosomes,
X, and Y generally remained. BED files were further sorted with the scaffold
column version-sorted rather than lexicographically sorted. This made the BEDs
compatible with the order of scaffolds in the hs37d5 reference.

### Filtration Results:

| Mask 1 | Mask 2 | Mask 3 | Mask 4 | Filter 1 | Variant Count |
|:------:|:------:|:------:|:------:|:--------:| -------------:|
|  |  |  |  |  | 116,324 |
<!-- cpg2_windows	90035 -->
| X |  |  |  |  | 90,035 |
<!-- genomicSuperDups	114014 -->
|  | X |  |  |  | 114,014 |
<!-- hg19_TRF_nochrprefix_nounplaced	114832 -->
|  |  | X |  |  | 114,832 |
<!-- hs37m_mask35_50	95822 -->
|  |  |  | X |  | 95,822 |
<!-- hg19lifted-cCREs_V4_nochrprefix_nounplaced	27808 -->
|  |  |  |  | X | 27,808 |
<!-- cpg2_windows	genomicSuperDups	88237 -->
| X | X |  |  |  | 88,237 |
<!-- cpg2_windows	hg19_TRF_nochrprefix_nounplaced	88669 -->
| X |  | X |  |  | 88,669 |
<!-- cpg2_windows	hs37m_mask35_50	75612 -->
| X |  |  | X |  | 75,612 |
<!-- cpg2_windows	hg19lifted-cCREs_V4_nochrprefix_nounplaced	20451 -->
| X |  |  |  | X | 20,451 |
<!-- genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	112547 -->
|  | X | X |  |  | 112,547 |
<!-- genomicSuperDups	hs37m_mask35_50	94742 -->
|  | X |  | X |  | 94,742 |
<!-- genomicSuperDups	hg19lifted-cCREs_V4_nochrprefix_nounplaced	27291 -->
|  | X |  |  | X | 27,291 |
<!-- hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	95095 -->
|  |  | X | X |  | 95,095 |
<!-- hg19_TRF_nochrprefix_nounplaced	hg19lifted-cCREs_V4_nochrprefix_nounplaced	27458 -->
|  |  | X |  | X | 27,458 |
<!-- hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	25202 -->
|  |  |  | X | X | 25,202 |
<!-- cpg2_windows	genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	86894 -->
| X | X | X |  |  | 86,894 |
<!-- cpg2_windows	genomicSuperDups	hs37m_mask35_50	74796 -->
| X | X |  | X |  | 74,796 |
<!-- cpg2_windows	genomicSuperDups	hg19lifted-cCREs_V4_nochrprefix_nounplaced	20074 -->
| X | X |  |  | X | 20,074 |
<!-- cpg2_windows	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	74971 -->
| X |  | X | X |  | 74,971 |
<!-- cpg2_windows	hg19_TRF_nochrprefix_nounplaced	hg19lifted-cCREs_V4_nochrprefix_nounplaced	20145 -->
| X |  | X |  | X | 20,145 |
<!-- cpg2_windows	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	18654 -->
| X |  |  | X | X | 18,654 |
<!-- genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	94027 -->
|  | X | X | X |  | 94,027 |
<!-- genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hg19lifted-cCREs_V4_nochrprefix_nounplaced	26945 -->
|  | X | X |  | X | 26,945 |
<!-- genomicSuperDups	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	24904 -->
|  | X |  | X | X | 24,904 |
<!-- hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	25017 -->
|  |  | X | X | X | 25,017 |
<!-- cpg2_windows	genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	74166 -->
| X | X | X | X |  | 74,166 |
<!-- cpg2_windows	genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hg19lifted-cCREs_V4_nochrprefix_nounplaced	19771 -->
| X | X | X |  | X | 19,771 |
<!-- cpg2_windows	genomicSuperDups	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	18446 -->
| X | X |  | X | X | 18,446 |
<!-- cpg2_windows	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	18499 -->
| X |  | X | X | X | 18,499 |
<!-- genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	24722 -->
|  | X | X | X | X | 24,722 |
<!-- cpg2_windows	genomicSuperDups	hg19_TRF_nochrprefix_nounplaced	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	18293 -->
| X | X | X | X | X | 18,293 |

Since our target is 20-30k variants and some of these masks are less
applicable (ex. CpG islands may matter less, TRF may be redundant with
hs37m_mask35_50), it seems that we should go with a minimum of
genomicSuperDups, hs37m_mask35_50, and the ENCODE cCREs filter.

One thing to note is that our mappability/TRF mask was applied directly
to the variants, rather than to flanking sequence. We may want to screen
for low complexity sequence within the entire oligo region (85+ bp on
each flank):

| Mask 2 | Mask 4 | Filter 1 | Mask 3 slop | Variant Count |
|:------:|:------:|:--------:|:-----------:| -------------:|
<!-- genomicSuperDups	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	hg19_TRF_nochrprefix_nounplaced(slop0)	24722 -->
| X | X | X | 0 | 24,722 |
<!-- genomicSuperDups	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	hg19_TRF_nochrprefix_nounplaced(slop85)	24255 -->
| X | X | X | 85 | 24,255 |
<!-- genomicSuperDups	hs37m_mask35_50	hg19lifted-cCREs_V4_nochrprefix_nounplaced	hg19_TRF_nochrprefix_nounplaced(slop110)	24143 -->
| X | X | X | 110 | 24,143 |

We also should screen for variants in close proximity to each other,
especially if they are in the same Sprime tract.

## Final oligo set results:

### Masks/filters used:

1. UCSC genomicSuperDups mask
1. Heng Li's 35-mer mappability mask
1. ENCODE candidate *cis* regulatory elements v4 filter (lifted from GRCh38 to hs37d5)
1. Dropped indel variants
1. Dropped variants where the Sprime allele was the REF allele
1. Dropped variants where the Sprime allele was a spanned deletion (* in VCF notation) -- this was added post-hoc

### Final oligo count:

| Filtering stage | Variant count |
|:---------------:| -------------:|
| First 3 masks | 24,904 |
| First 5 masks | 22,379 |
| All 6 masks | 22,313 |

Thus the total number of oligos ordered was 22,313 variants * 2 oligos/variant + 709 controls = 45,335 oligos.

## Code run for filtering:

```
PATH="/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/MPRA/PIB92_pilot:${PATH}" ./prepare_mpra_oligos.sh site_selection/ match_rates/ designed_oligos/PIB92_pilot_MPRA_final_variants

PATH="/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/MPRA/PIB92_pilot:${PATH}" ./add_mpra_adapters_controls.sh designed_oligos/PIB92_pilot_MPRA_final_variants ../ReillyLab_MPRA_Controls.txt > designed_oligos/PIB92_pilot_MPRA_final_oligos_inclControls_toOrder.txt
```

### Notes about output:

The output of `add_mpra_adapters_controls.sh` wasn't quite the same as
the format of the controls because I misunderstood a tab as a newline,
so I thought it was supposed to be a FASTA-like format rather than a
simple TSV.

Somewhat ironically, the way the interleaving gets done puts the file
in the correct format, but I undo this with the subsequent AWK script...

In addition, I made a mistake in the filtering code and failed to filter
out cases where the Sprime allele is a spanning deletion allele (* in
VCF notation), so Steve had to filter 66 of those out before ordering.

Steve fixed this (and some issues with the controls) for final ordering,
and checked the oligos via BLAT.

### Notes for future oligo design:

Adjust scripts to output TSV rather than FASTA-like, and drop variants
with non-ACGT alleles.
