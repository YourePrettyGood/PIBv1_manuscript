#PIBv1 sample and site filter/mask QC:

#Load the libraries:
library(tidyverse)

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

#Get to the directory:
pushd('[path redacted]/PIBv1_QC/')

#Load the data:
#Sample metadata file:
PIBv1_metadata <- read_tsv('[path redacted]/Metadata/PIBv1_metadata_v0.2.tsv',
                           col_types='ccccclcccnncccclc',
                           na=c("NA", "")) %>%
   mutate(Region=factor(Region,
                        levels=c("America", "Europe",
                                 "Africa", "Middle_East",
                                 "Central_South_Asia", "East_Asia",
                                 "Island_Southeast_Asia", "Oceania"),
                        labels=c("America", "Europe",
                                 "Africa", "Middle East",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania")))

#VerifyBamID2 contamination results:
VBID2_contam <- read_tsv('PIBv1_VBID2.tsv')

#haplocheck contamination results:
haplocheck <- read_tsv('PIBv1_haplocheck.tsv')

#Set up the level order of the filter set for plotting:
#filterset_levels <- c("nofilters", "VQSRpass",
#                      "noArchaicMaskedNoGTmasks", "noArchaicMaskedAllGTmasks",
#                      "allfilters", "allmasks")
filterset_levels <- c("nofilters", "VQSRpass",
                      "noArchaicMaskedNoGTmasks", "noArchaicMaskedAllGTmasks",
                      "allfilters", "allmasks",
                      "VQSRpassAllGTmasks", "VQSRpassCpGAllGTmasks",
                      "VQSRpassSegDupAllGTmasks", "VQSRpassMappability50AllGTmasks",
                      "VQSRpassArchaicMaskedAllGTmasks", "VQSRpassExcHet0.0001AllGTmasks",
                      "VQSRpassHWE0.0001AllGTmasks", "VQSRpassIndelProxAllGTmasks",
                      "VQSRpassMissingness0.05AllGTmasks", "VQSRpassMissingness0AllGTmasks")
filterset_extrema <- c("nofilters", "VQSRpass",
                       "allfilters", "allmasks")
filterset_vqsrbed <- c("VQSRpassAllGTmasks", "VQSRpassCpGAllGTmasks",
                        "VQSRpassSegDupAllGTmasks", "VQSRpassMappability50AllGTmasks",
                        "VQSRpassArchaicMaskedAllGTmasks")
filterset_vqsrsampledep <- c("VQSRpassAllGTmasks", "VQSRpassExcHet0.0001AllGTmasks",
                             "VQSRpassHWE0.0001AllGTmasks", "VQSRpassIndelProxAllGTmasks",
                             "VQSRpassMissingness0.05AllGTmasks", "VQSRpassMissingness0AllGTmasks")

#Note: We only use the autosomal estimates, since the X and Y are funky
# for trio and heterozygosity stats.
#Genotype discordance isn't available for the X or Y anyways due to our
# ground truth dataset (autosomal Human Origins array data).
#Basic sample stats from VCF:
smplstats <- read_tsv('PIBv1_smplstats_allfiltersets.tsv.gz') %>%
   mutate(FilterSet=factor(FilterSet,
                           levels=filterset_levels)) %>%
   filter(Chromosomes == "autosomes")

#Genotype discordance stats from VCF:
gtcheck <- read_tsv('PIBv1_gtcheck_allfiltersets.tsv.gz') %>%
   mutate(FilterSet=factor(FilterSet, levels=filterset_levels)) %>%
   filter(Chromosomes == "autosomes")

#Trio-based Mendelian error stats from VCF:
triostats <- read_tsv('PIBv1_triostats_allfiltersets.tsv.gz') %>%
   mutate(FilterSet=factor(FilterSet, levels=filterset_levels)) %>%
   filter(Chromosomes == "autosomes")

#Some basic QC plots:

#Skipped in supplement
#General regional and population breakdown of samples:
PIBv1_metadata %>%
   ggplot(aes(x=Region, colour=Region, fill=Region)) +
      geom_bar() +
      geom_text(stat="count", aes(label=..count..), colour="black", vjust=1.5) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      guides(colour=FALSE, fill=FALSE) +
      labs(y="Sample size") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_samplecount_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
PIBv1_metadata %>%
   ggplot(aes(x=Data_from, colour=Data_from, fill=Data_from)) +
      geom_bar() +
      geom_text(stat="count", aes(label=..count..), colour="black", vjust=1.5) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      guides(colour=FALSE, fill=FALSE) +
      labs(x="Data source",
           y="Sample size") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_samplecount_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Silly little plot of sample size of analysis groups split by region:
PIBv1_metadata %>%
   ggplot(aes(x=Region, group=AnalysisGroup)) +
      geom_bar(position=position_dodge2(preserve="total")) +
      theme_bw() +
      labs(y="Sample size") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_samplecount_byAnalysisGroup_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Sample QC:
#Supplementary figure S7 in PIBv1 manuscript:
#Per-sample missingness:
#We only retain the filtersets that don't filter on missingness
# since those that do filter on missingness result in all 0s.
smplstats %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   filter(FilterSet %in% c("nofilters", "VQSRpass")) %>%
   ggplot(aes(x=Region, y=MissingRate*100, colour=Region)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      ylim(0, NA) +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      facet_wrap(~ FilterSet, ncol=2) +
      labs(y="Per-sample missingness (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('SampleQC/PIBv1_persample_missingness_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('SampleQC/PIBv1_persample_missingness_byRegion.png',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
smplstats %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   filter(FilterSet %in% c("nofilters", "VQSRpass")) %>%
   ggplot(aes(x=Data_from, y=MissingRate*100, colour=Data_from)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      ylim(0, NA) +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      facet_wrap(~ FilterSet, ncol=2) +
      labs(x="Data source",
           y="Per-sample missingness (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('SampleQC/PIBv1_persample_missingness_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('SampleQC/PIBv1_persample_missingness_byDataSource.png',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Corrected version plotted in PIBv1_HeterozygosityVarPerGenome_20230523.R
#Per-sample heterozygosity:
smplstats %>%
   filter(FilterSet %in% filterset_extrema) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, y=HetRate*100, colour=Region)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      facet_wrap(~ FilterSet, ncol=2) +
      labs(y="Heterozygosity (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_persample_heterozygosity_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
smplstats %>%
   filter(FilterSet %in% filterset_extrema) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Data_from, y=HetRate*100, colour=Data_from)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      facet_wrap(~ FilterSet, ncol=2) +
      labs(x="Data source",
           y="Heterozygosity (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_persample_heterozygosity_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Mitochondrial heteroplasmy contamination status:
haplocheck %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, colour=`Contamination Status`, fill=`Contamination Status`)) +
      geom_bar(position="fill") +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      labs(y="% samples with heteroplasmy",
           colour="mtDNA contaminated?",
           fill="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_heteroplasmy_proportion_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
haplocheck %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Data_from, colour=`Contamination Status`, fill=`Contamination Status`)) +
      geom_bar(position="fill") +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      labs(x="Data source",
           y="% samples with heteroplasmy",
           colour="mtDNA contaminated?",
           fill="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_heteroplasmy_proportion_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Now plot counts rather than proportion:
haplocheck %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, colour=`Contamination Status`, fill=`Contamination Status`)) +
      geom_bar(position="stack") +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      labs(y="# samples with heteroplasmy",
           colour="mtDNA contaminated?",
           fill="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_heteroplasmy_count_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
haplocheck %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Data_from, colour=`Contamination Status`, fill=`Contamination Status`)) +
      geom_bar(position="stack") +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      labs(x="Data source",
           y="# samples with heteroplasmy",
           colour="mtDNA contaminated?",
           fill="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_heteroplasmy_count_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S10 in PIBv1 manuscript:
#VerifyBamID2 contamination estimates:
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, y=FREEMIX*100, colour=Region)) +
      geom_boxplot() +
      geom_jitter(size=0.75, alpha=0.5) +
      geom_hline(yintercept=2, col="red", linetype=2) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      labs(y="Estimated contamination (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('Contamination/PIBv1_VBID2_contamination_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('Contamination/PIBv1_VBID2_contamination_byRegion.png',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Data_from, y=FREEMIX*100, colour=Data_from)) +
      geom_boxplot() +
      geom_jitter(size=0.75, alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      guides(colour=FALSE) +
      labs(x="Data source",
           y="Estimated contamination (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave('PIBv1_VBID2_contamination_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Stratify VBID2 contamination estimates by heteroplasmy contamination status:
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   inner_join(haplocheck %>%
                 select(Sample, `Contamination Status`),
              by=c("Sample"="Sample")) %>%
   ggplot(aes(x=Region, y=FREEMIX*100, colour=`Contamination Status`)) +
      geom_jitter(size=1, alpha=0.7) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      labs(y="Estimated contamination (%)",
           colour="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_VBID2_contamination_byHeteroplasmy_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement:
#Repeat the plot, but by data source:
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   inner_join(haplocheck %>%
                 select(Sample, `Contamination Status`),
              by=c("Sample"="Sample")) %>%
   ggplot(aes(x=Data_from, y=FREEMIX*100, colour=`Contamination Status`)) +
      geom_jitter(size=1, alpha=0.7) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      labs(x="Data source",
           y="Estimated contamination (%)",
           colour="mtDNA contaminated?") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_VBID2_contamination_byHeteroplasmy_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Plot heterozygosity outliers per-region:
smplstats %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks") %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region, AnalysisGroup),
             by=c("Sample"="SampleID")) %>%
   group_by(Region) %>%
   mutate(HetOutlier=case_when(HetRate*100>=mean(HetRate*100)+3*sd(HetRate*100) ~ "High",
                               HetRate*100<=mean(HetRate*100)-3*sd(HetRate*100) ~ "Low",
                               TRUE ~ "No"),
          SampleSize=n()) %>%
   mutate(AnalysisGroup=str_c(as.character(AnalysisGroup), " (", as.character(SampleSize), ")"),
          HetOutlier=factor(HetOutlier,
                            levels=c("No", "Low", "High"))) %>%
   ggplot(aes(x=AnalysisGroup, y=HetRate*100, colour=HetOutlier)) +
      geom_jitter(alpha=0.5) +
      ylim(0, NA) +
      theme_bw() +
      scale_colour_manual(values=c("black", "purple", "red")) +
      labs(x="Population",
           y="Heterozygosity (%)",
           colour="Heterozygosity Outlier?") +
      facet_wrap(~ Region, ncol=2, scales="free") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            legend.position="top")
ggsave('SampleQC/PIBv1_heterozygosity_outliers_byRegion.pdf',
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Supplementary figure S11 in PIBv1 manuscript:
#Plot heterozygosity outliers per-analysis group:
smplstats %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks") %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region, AnalysisGroup),
             by=c("Sample"="SampleID")) %>%
   group_by(AnalysisGroup) %>%
   mutate(HetOutlier=case_when(HetRate*100>=mean(HetRate*100)+3*sd(HetRate*100) ~ "High",
                               HetRate*100<=mean(HetRate*100)-3*sd(HetRate*100) ~ "Low",
                               TRUE ~ "No"),
          SampleSize=n()) %>%
   mutate(AnalysisGroup=str_c(as.character(AnalysisGroup), " (", as.character(SampleSize), ")"),
          HetOutlier=factor(HetOutlier,
                            levels=c("No", "Low", "High"))) %>%
   ggplot(aes(x=AnalysisGroup, y=HetRate*100, colour=HetOutlier)) +
      geom_jitter(alpha=0.5) +
      ylim(0, NA) +
      theme_bw() +
      scale_colour_manual(values=c("black", "purple", "red")) +
      labs(x="Population",
           y="Heterozygosity (%)",
           colour="Heterozygosity Outlier?") +
      facet_wrap(~ Region, ncol=2, scales="free") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            legend.position="top")
ggsave('SampleQC/PIBv1_heterozygosity_outliers_byAnalysisGroup.pdf',
       width=24.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave('SampleQC/PIBv1_heterozygosity_outliers_byAnalysisGroup.png',
       width=24.0,
       height=24.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat but plotting points by heterozygosity scaled by a fixed denominator:
smplstats %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks") %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region, AnalysisGroup),
             by=c("Sample"="SampleID")) %>%
   group_by(AnalysisGroup) %>%
   mutate(HetOutlier=case_when(HetRate*100>=mean(HetRate*100)+3*sd(HetRate*100) ~ "High",
                               HetRate*100<=mean(HetRate*100)-3*sd(HetRate*100) ~ "Low",
                               TRUE ~ "No"),
          SampleSize=n()) %>%
   mutate(AnalysisGroup=str_c(as.character(AnalysisGroup), " (", as.character(SampleSize), ")"),
          HetOutlier=factor(HetOutlier,
                            levels=c("No", "Low", "High"))) %>%
   ggplot(aes(x=AnalysisGroup, y=Het*100/2684573005, colour=HetOutlier)) +
      geom_jitter(alpha=0.5) +
      ylim(0, NA) +
      theme_bw() +
      scale_colour_manual(values=c("black", "purple", "red")) +
      labs(x="Population",
           y="Heterozygosity (%)",
           colour="Heterozygosity Outlier?") +
      facet_wrap(~ Region, ncol=2, scales="free_x") +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            legend.position="top")
ggsave('SampleQC/PIBv1_heterozygosity_outliers_correctScale_byAnalysisGroup.pdf',
       width=24.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave('SampleQC/PIBv1_heterozygosity_outliers_correctScale_byAnalysisGroup.png',
       width=24.0,
       height=24.0,
       units="cm",
       dpi=500)

#Generate a table of samples to exclude due to contamination or being
# heterozygosity outliers:
#From per-region outlier calculations:
smplstats %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks") %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region, AnalysisGroup, Outlier),
             by=c("Sample"="SampleID")) %>%
   left_join(VBID2_contam %>%
                mutate(Sample=`#SEQ_ID`) %>%
                select(Sample, FREEMIX),
             by=c("Sample"="Sample")) %>%
   group_by(Region) %>%
   mutate(HetOutlier=case_when(HetRate*100>=mean(HetRate*100)+3*sd(HetRate*100) ~ "High",
                               HetRate*100<=mean(HetRate*100)-3*sd(HetRate*100) ~ "Low",
                               TRUE ~ "No"),
          Contaminated=FREEMIX>=0.02) %>%
   rename(ChoinOutlier=Outlier) %>%
   filter(Contaminated | HetOutlier != "No" | ChoinOutlier != "FALSE") %>%
   transmute(Sample=Sample, Region=Region, AnalysisGroup=AnalysisGroup,
             Contaminated_VBID2=Contaminated, Outlier_Heterozygosity=HetOutlier,
             Outlier_Choin=ChoinOutlier, HeterozygosityRate=HetRate, VBID2_ContaminationRate=FREEMIX) %>%
   arrange(Sample) %>%
   write_tsv('SampleQC/PIBv1_ContaminatedOrPerRegionHetOutlier_samples.tsv')

#From per-analysis group outlier calculations:
smplstats %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks") %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region, AnalysisGroup, Outlier),
             by=c("Sample"="SampleID")) %>%
   left_join(VBID2_contam %>%
                mutate(Sample=`#SEQ_ID`) %>%
                select(Sample, FREEMIX),
             by=c("Sample"="Sample")) %>%
   group_by(AnalysisGroup) %>%
   mutate(HetOutlier=case_when(HetRate*100>=mean(HetRate*100)+3*sd(HetRate*100) ~ "High",
                               HetRate*100<=mean(HetRate*100)-3*sd(HetRate*100) ~ "Low",
                               TRUE ~ "No"),
          Contaminated=FREEMIX>=0.02) %>%
   rename(ChoinOutlier=Outlier) %>%
   filter(Contaminated | HetOutlier != "No" | ChoinOutlier != "FALSE") %>%
   transmute(Sample=Sample, Region=Region, AnalysisGroup=AnalysisGroup,
             Contaminated_VBID2=Contaminated, Outlier_Heterozygosity=HetOutlier,
             Outlier_Choin=ChoinOutlier, HeterozygosityRate=HetRate, VBID2_ContaminationRate=FREEMIX) %>%
   arrange(Sample) %>%
   write_tsv('SampleQC/PIBv1_ContaminatedOrPerAnalysisGroupHetOutlier_samples.tsv')


#Variant QC:
#Skipped in supplement
#Trio-based Mendelian error rates by filterset:
#MErate1 is the Mendelian error rate when normalized by the number
# of sites where there is at least one ALT allele in the trio.
#This error rate basically excludes sites where all members of
# the trio are homozygous REF, so the rate is higher than MErate2.
#First plot only the extremes of the filter sets:
triostats %>%
   filter(FilterSet %in% filterset_extrema) %>%
   ggplot(aes(x=ValidSitesWithAlt, y=MErate1*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid sites with >=1 ALT allele in trio",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_trioerrors_ValidSitesWithAltDenom_byFiltersetExtremes.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat including BED-based filters:
triostats %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrbed)) %>%
   ggplot(aes(x=ValidSitesWithAlt, y=MErate1*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid sites with >=1 ALT allele in trio",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_trioerrors_ValidSitesWithAltDenom_byFiltersetExtremesBEDs.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat including sample-dependent filters instead:
triostats %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrsampledep)) %>%
   ggplot(aes(x=ValidSitesWithAlt, y=MErate1*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid sites with >=1 ALT allele in trio",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_trioerrors_ValidSitesWithAltDenom_byFiltersetExtremesSampledep.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#MErate2 is the Mendelian error rate when normalized by the number
# of sites where all trio genotypes pass filters.
#Skipped in supplement
#First plot only the extremes of the filter sets:
triostats %>%
   filter(FilterSet %in% filterset_extrema) %>%
   ggplot(aes(x=ValidSites, y=MErate2*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid trio sites",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_trioerrors_ValidSitesDenom_byFiltersetExtremes.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S4 in PIBv1 manuscript:
#Repeat including BED-based filters:
triostats %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrbed)) %>%
   ggplot(aes(x=ValidSites, y=MErate2*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid trio sites",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_trioerrors_ValidSitesDenom_byFiltersetExtremesBEDs.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S3 in PIBv1 manuscript:
#Repeat including sample-dependent filters instead:
triostats %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrsampledep)) %>%
   ggplot(aes(x=ValidSites, y=MErate2*100, colour=FilterSet)) +
      geom_point() +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE)) +
      labs(x="Valid trio sites",
           y="Mendelian Error Rate (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="top")
ggsave('PIBv1_trioerrors_ValidSitesDenom_byFiltersetExtremesSampledep.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Genotype discordance by region, faceted by filterset extrema:
gtcheck %>%
   filter(FilterSet %in% filterset_extrema) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=Region, y=GenotypeDiscordance*100, colour=Region)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=FALSE) +
      labs(y="Genotype Discordance (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      facet_wrap(~ FilterSet, ncol=2)
ggsave('PIBv1_genotypediscordance_byFiltersetExtrema_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
gtcheck %>%
   filter(FilterSet %in% filterset_extrema) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=Data_from, y=GenotypeDiscordance*100, colour=Data_from)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=FALSE) +
      labs(x="Data source",
           y="Genotype Discordance (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      facet_wrap(~ FilterSet, ncol=2)
ggsave('PIBv1_genotypediscordance_byFiltersetExtrema_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#And now only jitter so we can see the sample size differences:
gtcheck %>%
   filter(FilterSet %in% filterset_extrema) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=Data_from, y=GenotypeDiscordance*100, colour=Data_from)) +
      geom_jitter(size=0.5, alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=FALSE) +
      labs(x="Data source",
           y="Genotype Discordance (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      facet_wrap(~ FilterSet, ncol=2)
ggsave('PIBv1_genotypediscordance_jitter_byFiltersetExtrema_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S8 in PIBv1 manuscript:
#Let's plot a different way to see the effects of the filters.
#Facet by Region, and have the x-axis and colour be FilterSet:
#We include the BED filters in this one.
gtcheck %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrbed)) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=FilterSet, y=GenotypeDiscordance*100, colour=FilterSet)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=FALSE) +
      labs(y="Genotype Discordance (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      facet_wrap(~ Region, ncol=2)
ggsave('PIBv1_genotypediscordance_byRegion_byFiltersetExtremaBEDs.pdf',
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Supplementary figure S9 in PIBv1 manuscript:
#Include the sample-dependent filters instead:
gtcheck %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrsampledep)) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=FilterSet, y=GenotypeDiscordance*100, colour=FilterSet)) +
      geom_boxplot() +
      geom_jitter(size=0.25, alpha=0.3) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=FALSE) +
      labs(y="Genotype Discordance (%)") +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      facet_wrap(~ Region, ncol=2)
ggsave('PIBv1_genotypediscordance_byRegion_byFiltersetExtremaSampledep.pdf',
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Now plot the genotype discordance of each filterset against the denominator:
#BED filters:
gtcheck %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrbed)) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=SitesCompared, y=GenotypeDiscordance*100, colour=FilterSet)) +
      geom_point(size=0.5, alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE, override.aes=list(size=2))) +
      labs(x="# Sites Compared",
           y="Genotype Discordance (%)") +
      theme(legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_genotypediscordance_vsSitesCompared_byFiltersetExtremaBEDs.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Sample-dependent per-site filters:
gtcheck %>%
   filter(FilterSet %in% c(filterset_extrema, filterset_vqsrsampledep)) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("QuerySample"="SampleID")) %>%
   ggplot(aes(x=SitesCompared, y=GenotypeDiscordance*100, colour=FilterSet)) +
      geom_point(size=0.5, alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=2, byrow=TRUE, override.aes=list(size=2))) +
      labs(x="# Sites Compared",
           y="Genotype Discordance (%)") +
      theme(legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_genotypediscordance_vsSitesCompared_byFiltersetExtremaSampledep.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)


#Skipped in supplement
#We might suspect that contamination would lead to genotype discordance,
# especially if contamination happened in between the array being assayed
# and the WGS data being generated.
#Genotype discordance against estimated WGS contamination (from VerifyBamID2):
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Region),
             by=c("Sample"="SampleID")) %>%
   inner_join(gtcheck,
              by=c("Sample"="QuerySample")) %>%
   filter(FilterSet %in% filterset_extrema) %>%
   ggplot(aes(x=FREEMIX*100, y=GenotypeDiscordance*100, colour=Region)) +
      geom_point(size=0.5, alpha=0.8) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(ncol=4, byrow=TRUE, override.aes=list(size=3))) +
      labs(x="Estimated contamination (%)",
           y="Genotype discordance (%)") +
      facet_wrap(~ FilterSet, ncol=2) +
      theme(legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_genotypediscordance_byContamination_byFiltersetExtrema_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Repeat the plot, but by data source:
VBID2_contam %>%
   mutate(Sample=`#SEQ_ID`) %>%
   select(Sample, FREEMIX) %>%
   left_join(PIBv1_metadata %>%
                select(SampleID, Data_from),
             by=c("Sample"="SampleID")) %>%
   inner_join(gtcheck,
              by=c("Sample"="QuerySample")) %>%
   filter(FilterSet %in% filterset_extrema) %>%
   ggplot(aes(x=FREEMIX*100, y=GenotypeDiscordance*100, colour=Data_from)) +
      geom_point(size=0.5, alpha=0.8) +
      theme_bw() +
      scale_colour_brewer(palette="Paired") +
      guides(colour=guide_legend(nrow=2, byrow=TRUE, override.aes=list(size=3))) +
      labs(x="Estimated contamination (%)",
           y="Genotype discordance (%)",
           colour="Data source") +
      facet_wrap(~ FilterSet, ncol=2) +
      theme(legend.position="top",
            legend.spacing.y=unit(0.02, "cm"))
ggsave('PIBv1_genotypediscordance_byContamination_byFiltersetExtrema_byDataSource.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#All PCA plotting code found in a different script, so ignore below:
#PCA:
#Function to make the PCA biplot:
pca_biplot <- function(loadings, pve, PCs, colour, palette, title, pointsize=1) {
   pve_precision <- 3
   pc_one <- paste0("PC", PCs[1])
   pc_one_pve <- signif(pve[pve$PC == PCs[1], "PVE"], pve_precision)
   pc_two <- paste0("PC", PCs[2])
   pc_two_pve <- signif(pve[pve$PC == PCs[2], "PVE"], pve_precision)
   loadings %>%
      filter(PC %in% PCs) %>%
      pivot_wider(names_from="PC",
                  names_prefix="PC",
                  values_from="Loading") %>%
      ggplot(aes_string(x=pc_one, y=pc_two, colour=colour)) +
         geom_point(size=pointsize, alpha=0.7) +
         theme_bw() +
         {if(length(palette) > 1) {scale_colour_manual(values=palette)} else {scale_colour_brewer(palette=palette)}} +
         labs(x=paste0(pc_one, " (", pc_one_pve, "%)"),
              y=paste0(pc_two, " (", pc_two_pve, "%)"),
              title=title,
              subtitle=paste0(pc_one, " vs. ", pc_two))
}

#Load the PCA eigenvalues and eigenvectors from SNPRelate:
pca_eigval <- read.table('PIBv1_autosomes_VQSRpassMissingness0.05AllGTmasks_allMAFge0.05_bSNPs_LDpruned_SNPRelate_PCA.eigenval',
                         header=FALSE,
                         colClasses=c("numeric"),
                         col.names=c("Eigenvalue"),
                         row.names=NULL)
pca_eigvec <- read.table('PIBv1_autosomes_VQSRpassMissingness0.05AllGTmasks_allMAFge0.05_bSNPs_LDpruned_SNPRelate_PCA.eigenvec',
                         colClasses=c("character", rep("numeric", nrow(pca_eigval))),
                         col.names=c("ID", paste0("PC", seq(1, nrow(pca_eigval)))),
                         skip=1)

#Reformat a bit:
pca_eigval <- pca_eigval %>%
   mutate(PC=as.character(row_number()))
pca_eigvec <- pca_eigvec %>%
   pivot_longer(cols=-ID,
                names_prefix="PC",
                names_to="PC",
                values_to="Loading")
#Create a data.frame with the eigenvalues and percent of variance explained by each:
pve_df <- pca_eigval %>%
   mutate(PVE=abs(Eigenvalue)*100/sum(abs(Eigenvalue)))

#Create a data.frame for the PC loadings (with metadata):
pca_df <- pca_eigvec %>%
   inner_join(PIBv1_metadata, by=c("ID"="SampleID"))

#Now plot PC 1 vs. 2, labeled by region:
pca_biplot(pca_df, pve_df, c(1, 2), "Region", "Paired", "")
ggsave('PIBv1_autosomes_minimalFilter_allMAFge0.05_LDprunedChoin_PCA_1v2_byRegion.pdf',
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Quick check on the visual outliers:
pca_df %>%
   filter(PC %in% c(1, 2)) %>%
   pivot_wider(names_from="PC",
               names_prefix="PC",
               values_from="Loading") %>%
   filter(PC1 < 0.025, PC1 > 0.000, PC2 > 0.00, PC2 < 0.02) %>%
   View()


#Clean up:
rm(filterset_levels, filterset_extrema, filterset_vqsrbed, filterset_vqsrsampledep)
rm(PIBv1_metadata)
rm(VBID2_contam, haplocheck)
rm(smplstats, gtcheck, triostats)
