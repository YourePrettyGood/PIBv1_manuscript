#PIBv1 lowdepth QC:

#Load the libraries:
library(tidyverse)
library(ggrepel)
library(viridis)

#Helper functions for navigating directories:
pushd <- function(path) {
   assign("wd_stack", c(get0("wd_stack",
                             envir=globalenv(),
                             ifnotfound=character(0)),
                        getwd()),
          envir=globalenv());
   setwd(path)
}
popd <- function() {
   path_stack <- get("wd_stack",
                     envir=globalenv());
   setwd(path_stack[length(path_stack)]);
   assign("wd_stack",
          path_stack[-length(path_stack)],
          envir=globalenv())
}

#Get to the parent directory:
pushd('[path redacted]')

#Helper functions:
pca_screeplot <- function(pve, plot_prefix) {
   pve_plot <- pve %>%
      mutate(PC=as.numeric(PC)) %>%
      ggplot(aes(x=PC, y=PVE)) +
         geom_point() +
         theme_bw() +
         labs(y="Percent of Variance Explained") +
         ylim(0, NA)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
}
pca_biplot <- function(loadings, pve, PCs, colour, shape, palette, shape_palette, title, plot_prefix, pointsize=1, label=NULL, labelsize=1, labellinewidth=0.5) {
   if (!is.null(label)) {
      library(ggrepel)
   }
   pve_precision <- 3
   pc_one <- paste0("PC", PCs[1])
   pc_one_pve <- signif(pve[pve$PC == PCs[1], "PVE"], pve_precision)
   pc_two <- paste0("PC", PCs[2])
   pc_two_pve <- signif(pve[pve$PC == PCs[2], "PVE"], pve_precision)
   biplot <- loadings %>%
      filter(PC %in% PCs) %>%
      pivot_wider(names_from="PC",
                  names_prefix="PC",
                  values_from="Loading") %>%
      ggplot(data=., aes(x=.data[[pc_one]], y=.data[[pc_two]], colour=.data[[colour]], shape=.data[[shape]])) +
         geom_point(size=pointsize, alpha=0.5) +
         {
            if (!is.null(label)) {
               geom_text_repel(aes(label=.data[[label]]), size=labelsize, segment.size=labellinewidth)
            }
         } +
         theme_bw() +
         scale_colour_manual(values=palette) +
         scale_shape_manual(values=shape_palette) +
         labs(x=paste0(pc_one, " (", pc_one_pve, "%)"),
              y=paste0(pc_two, " (", pc_two_pve, "%)"),
              subtitle=title) +
         guides(shape=guide_legend(direction="vertical",
                                   ncol=3,
                                   byrow=TRUE,
                                   title.theme=element_text(size=7),
                                   label.theme=element_text(size=6),
                                   keywidth=0.5,
                                   keyheight=0.5,
                                   title=""),
                colour=guide_legend(direction="vertical",
                                    ncol=3,
                                    byrow=TRUE,
                                    title.theme=element_text(size=7),
                                    label.theme=element_text(size=6),
                                    keywidth=0.5,
                                    keyheight=0.5,
                                    title="")) +
         theme(legend.spacing.y=unit(0, "cm"))
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
}

#Load the data:
#Sample metadata:
#Better geographic region labels:
region_names <- c("America", "Europe",
                  "Africa", "Middle_East",
                  "Central_South_Asia", "East_Asia",
                  "Island_Southeast_Asia", "Oceania")
region_codes <- c("AMR", "EUR",
                  "AFR", "MDE",
                  "CSA", "EAS",
                  "ISEA", "OCN")
region_pretty <- c("America", "Europe",
                   "Africa", "Middle East",
                   "Central/South Asia", "East Asia",
                   "Island Southeast Asia", "Oceania")
PIBv1_metadata <- read_tsv('[path redacted]/Metadata/PIBv1_metadata_v0.3.tsv',
                           col_types='ccccclcccnncccclc',
                           na=c("NA", "")) %>%
   mutate(Region=factor(Region,
                        levels=region_names,
                        labels=region_codes))
#PCA results:
#Sample IDs:
ids <- read_tsv('QC_results/PIB184_HC3_lowdepth/ANGSD_sample_order.tsv',
                col_names=c("ID"))
#Sample covariance matrix:
cov_mat <- as.matrix(read.table('QC_results/PIB184_HC3_lowdepth/PIB184_HC3_autosomes.cov'))
#Per-sample QC summaries:
#Note that we add in metadata lines for the original three contaminated samples
# and adjust their batch IDs accordingly.
#The inner join also only retains samples used in the final set of 177 30x genomes.
qc_df <- read_tsv('QC_results/Summaries/PIB184_HC3_QC_summaries.tsv') %>%
   inner_join(bind_rows(PIBv1_metadata,
                        PIBv1_metadata %>%
                           filter(SampleID %in% c("UV0057HC", "UV0058HC", "UV0841HC")) %>%
                           mutate(SampleID=str_replace(SampleID, "HC", ""),
                                  Batch=case_when(str_detect(SampleID, "UV005[78]") ~ "RQ13381",
                                                  str_detect(SampleID, "UV0841") ~ "RQ15941"))),
              by=c("Sample"="SampleID"))

#Set the colour and shape palettes for the PCA biplots:
test_colour_palette <- viridis(12, end=0.9, option="C", direction=-1)
test_shape_palette <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

#Perform the PCA by taking the eigenvalue decomposition of the covariance matrix:
autosomes_eigen <- eigen(cov_mat)

#Generate data.frames for plotting the PCA:
#Percent of Variance Explained by each PC:
pve_df <- autosomes_eigen$values %>%
   as_tibble() %>%
   transmute(PC=as.character(row_number()),
             Eigenvalue=value) %>%
   mutate(PVE=abs(Eigenvalue)*100/sum(abs(Eigenvalue)))
#PC scores:
pca_df <- autosomes_eigen$vectors %>%
   as.data.frame() %>%
   as_tibble(.name_repair="check_unique") %>%
   rename_with(~ gsub("V", "PC", .x, fixed=TRUE)) %>%
   mutate(ID=ids$ID) %>%
   pivot_longer(cols=-c(ID),
                names_prefix="PC",
                names_to="PC",
                values_to="Loading") %>%
   inner_join(bind_rows(PIBv1_metadata,
                        PIBv1_metadata %>%
                           filter(SampleID %in% c("UV0057HC", "UV0058HC", "UV0841HC")) %>%
                           mutate(SampleID=str_replace(SampleID, "HC", ""),
                                  Batch=case_when(str_detect(SampleID, "UV005[78]") ~ "RQ13381",
                                                  str_detect(SampleID, "UV0841") ~ "RQ15941"))),
              by=c("ID"="SampleID"))

#Generate the screeplot:
pca_screeplot(pve_df, 'QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth')

#Supplementary figure S2 in the PIBv1 manuscript:
#Generate the biplot:
pca_biplot(loadings=pca_df,
           pve=pve_df,
           PCs=c(1, 2),
           colour="Population",
           shape="Batch",
           palette=test_colour_palette,
           shape_palette=test_shape_palette,
           title="",
           plot_prefix="QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_IDlabelled",
           label="ID",
           labellinewidth=0.01)

#Supplementary figure S1 in the PIBv1 manuscript:
#QC summary plot of Kraken/Bracken results:
bracken_plot <- qc_df %>%
   select(Sample, Region, Population,
          PctClassifiedMerged, PctClassifiedUnmerged,
          HumanAbundanceMerged, HumanAbundanceUnmerged) %>%
   pivot_longer(cols=c(PctClassifiedMerged, PctClassifiedUnmerged,
                       HumanAbundanceMerged, HumanAbundanceUnmerged),
                names_pattern="(.*)(Unmerged|Merged)",
                names_to=c("PercentType", "ReadType"),
                values_to="Percent") %>%
   pivot_wider(id_cols=c(Sample, Region, Population, ReadType),
               names_from="PercentType",
               values_from="Percent") %>%
   ggplot(aes(x=PctClassified, y=HumanAbundance, colour=Population)) +
      geom_point() +
      theme_bw() +
      xlim(NA, 100) +
      ylim(NA, 100) +
      labs(x="Percent of reads classified by Bracken (%)",
           y="Percent of reads classified as Human by Bracken (%)") +
      scale_colour_viridis_d(end=0.9) +
      facet_wrap(~ ReadType, ncol=1)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_BrackenHumanAbundance_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=bracken_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_BrackenHumanAbundance_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=bracken_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_BrackenHumanAbundance_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=bracken_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped this plot in the supplement
#QC summary plot of estimated genetic sex vs. self-reported sex:
sex_mismatch_plot <- qc_df %>%
   transmute(Sample=Sample,
             Sex=Sex,
             Region=Region,
             Population=Population,
             Batch=Batch,
             Y_normdp_merged=YNormalizedDepthMerged,
             Y_normdp_unmerged=YNormalizedDepthUnmerged,
             X_normdp_merged=XNormalizedDepthMerged,
             X_normdp_unmerged=XNormalizedDepthUnmerged,
             Sex_mismatch=case_when(YEstimatedSexMerged != Sex | YEstimatedSexUnmerged != Sex | XEstimatedSexMerged != Sex | XEstimatedSexUnmerged != Sex ~ Sample,
                                    TRUE ~ NA_character_)) %>%
   pivot_longer(cols=-c(Sample, Sex, Region, Population, Batch, Sex_mismatch),
                names_pattern="([XY])_normdp_(unmerged|merged)",
                names_to=c("Chromosome", "ReadType"),
                values_to="NormalizedDepth") %>%
   pivot_wider(id_cols=-c(Chromosome, NormalizedDepth),
               names_from="Chromosome",
               values_from="NormalizedDepth") %>%
   ggplot(aes(x=2*X, y=2*Y, shape=Sex, colour=Sex, label=Sex_mismatch)) +
      geom_point(alpha=0.7) +
      geom_text_repel(segment.size=0.3, size=4) +
      theme_bw() +
      labs(x="X dosage",
           y="Y dosage") +
      scale_x_continuous(limits=c(1.0, 2.2), expand=c(0, 0)) +
      scale_y_continuous(limits=c(0.0, 1.7), expand=c(0, 0)) +
      scale_colour_viridis_d(end=0.7) +
      facet_wrap(~ ReadType, ncol=1)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_sexmismatch_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=sex_mismatch_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_sexmismatch_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=sex_mismatch_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('QC_results/PIB184_HC3_lowdepth/PIBv1_lowdepth_sexmismatch_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=sex_mismatch_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Clean up:
rm(region_names, region_codes, region_pretty, PIBv1_metadata)
rm(ids, cov_mat, autosomes_eigen, pve_df, pca_df)
rm(qc_df)
rm(bracken_plot, sex_mismatch_plot)
rm(test_colour_palette, test_shape_palette)
rm(pca_screeplot, pca_biplot)
