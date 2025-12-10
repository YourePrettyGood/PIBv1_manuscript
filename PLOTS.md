# Plotting code for the manuscript titled "Long-term isolation and archaic introgression shape functional genetic variation in Near Oceania"

Many of these plots were made using:
- R 4.2.0 or 4.3.2
- tidyverse 1.3.1 or 2.0.0
- cowplot 1.1.1 or 1.1.3
- ggrepel
- viridisLite

## Dataset QC

[Supplementary figures S1-S2](/Plotting/PIBv1_lowdepth_QC_20230819.R)

[Supplementary figures S3-S4](/Plotting/PIBv1_QC_20220923.R)

[Supplementary figures S5-S6](/Plotting/PIBv1_HeterozygosityVarPerGenome_20230523.R)

[Supplementary figures S7-S11](/Plotting/PIBv1_QC_20220923.R)

## Geographic distribution of genetic variation

[Figure 1 panel B and supplementary figures S12-S13](/Plotting/PIBv1_SharedPrivateNovelKnown_20240708.R)

## Population structure analysis

### Principal components analysis (PCA)

[Figure 1 panel C and supplementary figures S14-S16,S20](/Plotting/PIBv1_final_PCAs_20231017.R)

## Demographic inference

### SMC++

[Figure 1 panel F and supplementary figures S32-S34](/Plotting/PIBv1_smcpp_20231006.R)

## Archaic introgression

### Match rate mixture models

[Figure 2 panels A,B,C,E and Supplementary figures S45-S48,S52-S55](/Plotting/PIBv1_Sprime_FiguresTables_maxgap0redo_20250203.R)

[Supplementary figure S56](/Plotting/PIBv1_Sprime_NEAGMMBMM_20230604.R)

[Supplementary figures S57-S59](/Plotting/PIBv1_Sprime_DENGMMBMMsupplement_20230605.R)

[Supplementary figure S60](/Plotting/Browning_Sprime_DENGMMBMM_20230919.R)

[Supplementary figure S61-S62](/Plotting/PIBv1_Sprime_DENGMMBMMsupplement_20230605.R)

## Selection and adaptive introgression

[Supplementary figures S64-S65](/Plotting/PIBv1_Sprime_KnownAIManhattanPlots_20230606.R)

[Supplementary figure S66](/Analysis_Pipelines/ReviewerResponses/Dani_results/PIBv1_pFCS_intersection_effect_20251210.R)

[Supplementary figure S76](/Analysis_Pipelines/ReviewerResponses/haplotype_divergence_plot.R) (see also [the description here](/ANALYSIS.md#haplotype-divergence-at-trps1))

Supplementary figures S77, S88, S90, S92, and S94 were made using [haplotype_plot.R](/Analysis_Pipelines/haplotype_plot.R) as run in the [haplotype_plot.nf](/Analysis_Pipelines/haplotype_plot.nf) pipeline.
Necessary input files that are not subject to controlled access can be found in
the [haplotype_plots/ subdirectory](/Plotting/haplotype_plots/), and the
Nextflow pipeline config file can be found in [Configs/Used/haplotype_plots/](/Configs/Used/haplotype_plots/PIBv1_haplotype_plots.config).

Figure 3C and supplementary figures S89, S91, and S93 were made using [PopArt](https://popart.maths.otago.ac.nz)
based on Nexus files output by the [haplotype_plot.nf](/Analysis_Pipelines/haplotype_plot.nf) pipeline
and then manually curated.
