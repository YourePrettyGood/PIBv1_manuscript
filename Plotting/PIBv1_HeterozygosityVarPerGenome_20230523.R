#PIBv1 heterozygosity and variants per genome comparison

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
#Sample metadata:
PIBv1_metadata <- read_tsv('[path redacted]/Metadata/PIBv1_metadata_v0.3.tsv',
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

#Sample variant/genotype statistics:
#Including all types of variants:
smplstats_allfiltersets <- read_tsv('PIBv1_smplstats_allfiltersets.tsv.gz')
#Including only the minimal filters and only SNPs:
smplstats <- read_tsv('smplstats_autosomes_minimalFilters_snps.tsv')

#Supplementary figure S5 in PIBv1 manuscript:
#Plot SNP-only per-sample autosomal heterozygosity:
smplstats %>%
   left_join(PIBv1_metadata,
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, y=Het*100/2684573005, colour=Region)) +
      geom_jitter(alpha=0.5) +
      theme_bw() +
      scale_y_continuous(breaks=seq(0, 0.1, by=0.01),
                         limits=c(0, NA)) +
      scale_colour_brewer(palette="Set2", direction=-1) +
      guides(colour="none") +
      labs(x="Region",
           y="Heterozygosity (%)") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(paste0('PIBv1_AutosomalHeterozygosity_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_AutosomalHeterozygosity_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_AutosomalHeterozygosity_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S6 in PIBv1 manuscript:
#Plot number of variants per genome including all variant types:
#This version is using only VQSR as filter to be comparable to
# 1000 Genomes Project 30x (Byrska-Bishop et al. 2022 Cell).
smplstats_allfiltersets %>%
   filter(FilterSet == "VQSRpass",
          Chromosomes == "autosomes") %>%
   left_join(PIBv1_metadata,
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, y=(Het+HomAlt)/1e6, colour=Region)) +
      geom_jitter(alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Set2", direction=-1) +
      guides(colour="none") +
      labs(x="Region",
           y="Variants per genome (M)") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(paste0('PIBv1_VarPerGenome_byRegion_VQSRonly_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_VarPerGenome_byRegion_VQSRonly_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_VarPerGenome_byRegion_VQSRonly_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Now with minimal filters:
smplstats_allfiltersets %>%
   filter(FilterSet == "VQSRpassMissingness0.05AllGTmasks",
          Chromosomes == "autosomes") %>%
   left_join(PIBv1_metadata,
             by=c("Sample"="SampleID")) %>%
   ggplot(aes(x=Region, y=(Het+HomAlt)/1e6, colour=Region)) +
      geom_jitter(alpha=0.5) +
      theme_bw() +
      scale_colour_brewer(palette="Set2", direction=-1) +
      guides(colour="none") +
      labs(x="Region",
           y="Variants per genome (M)") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(paste0('PIBv1_VarPerGenome_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_VarPerGenome_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_VarPerGenome_byRegion_minimalFilters_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Clean up:
rm(PIBv1_metadata, smplstats, smplstats_allfiltersets)
