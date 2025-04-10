#PIBv1 shared, private, novel, and known variant plots

#Load the libraries:
library(tidyverse)
library(cowplot)

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
pushd('[path redacted]/PIBv1_results/VariantMatching/')

#Load the data:
isec <- read_tsv('PIBv1_noVanuatu_combined_isecAlleleCounts.tsv',
                 col_types=cols(CHROM=col_character(),
                                Subset=col_character(),
                                VariantType=col_character(),
                                .default=col_integer()))

#Aggregate allele counts into novel, known by 1 database, and well-known:
#Aggregation happens across chromosomes, and alleles specific to gnomAD
# exomes or genomes are lumped together.
allele_counts <- isec %>%
   pivot_longer(cols=-c(CHROM, Subset, VariantType),
                names_to="GroupBitstring",
                values_to="Count") %>%
   filter(str_detect(GroupBitstring, "^1")) %>%
   mutate(Group=case_when(str_count(GroupBitstring, "1") > 1 ~ "Known",
                          TRUE ~ "Novel")) %>%
   group_by(Subset, VariantType, Group) %>%
   summarize(Count=sum(Count))

#Calculate the partition between shared vs. private alleles rather than all vs. private:
#This helps with plotting.
allele_counts <- allele_counts %>%
   filter(Subset != "ALL") %>%
   mutate(Subset=case_when(str_detect(Subset, "_all$", negate=TRUE) ~ str_c(Subset, "_private"),
                           TRUE ~ Subset)) %>%
   separate(col=Subset, into=c("Region", "VariantSubset"), sep="_")
shared_private_counts <- allele_counts %>%
   pivot_wider(id_cols=c(Region, VariantType, Group),
               names_from="VariantSubset",
               values_from="Count") %>%
   mutate(shared=all - private) %>%
   pivot_longer(cols=c(all, shared, private),
                names_to="VariantSubset",
                values_to="Count")

#Figure 1 panel B of the PIBv1 manuscript:
#Note that the outline of the bars for OCN were tweaked in Adobe Illustrator
# for the final figure.
#Main text figure of shared and private x novel and known:
shared_private_novel_known_plot <- shared_private_counts %>%
   filter(VariantSubset != "all",
          VariantType == "SNP") %>%
   unite(col="VariantGroup",
         c(Group, VariantSubset),
         sep="_") %>%
   pivot_wider(id_cols=c(Region, VariantType),
               names_from="VariantGroup",
               values_from="Count") %>%
   mutate(AllVariants=Known_shared+Known_private+Novel_shared+Novel_private) %>%
   mutate(frac_shared=(Known_shared+Novel_shared)/AllVariants,
          frac_knownprivate=Known_private/AllVariants,
          frac_novelprivate=Novel_private/AllVariants) %>%
   select(-c(Known_shared, Known_private, Novel_shared, Novel_private, AllVariants)) %>%
   pivot_longer(cols=c(frac_shared, frac_knownprivate, frac_novelprivate),
                names_to="VariantGroup",
                values_to="Fraction",
                names_prefix="frac_") %>%
   mutate(VariantGroup=factor(VariantGroup,
                              levels=c("novelprivate",
                                       "knownprivate",
                                       "shared"),
                              labels=c("Population-specific and novel",
                                       "Population-specific and known",
                                       "Shared among populations")),
          Percentage=str_c(signif(Fraction*100, digits=3), "%"),
          Region=factor(Region,
                        levels=c("AFR", "OCN", "EAS", "CSA",
                                 "EUR", "MDE", "AMR", "ISEA"))) %>%
   ggplot(aes(x=Region, y=Fraction, fill=VariantGroup)) +
      geom_col(position="stack") +
      geom_text(aes(label=Percentage), position=position_fill(vjust=0.5), col="white") +
      theme_classic() +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(labels=function(x) {paste0(format(x*100), "%")},
                         expand=c(0,0)) +
      scale_fill_manual(values=c("#bf3b36", "#f28d8f", "#92c4e8")) +
      guides(fill=guide_legend(nrow=2)) +
      labs(x="",
           y="",
           fill="") +
      theme(axis.text.x=element_text(size=12, face="bold"),
            axis.text.y=element_text(face="bold"),
            legend.position="bottom",
            panel.border=element_rect(fill=NA, colour="black"))
ggsave(paste0('PIBv1_Fig1B_SharedPrivateNovelKnown_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=shared_private_novel_known_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig1B_SharedPrivateNovelKnown_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=shared_private_novel_known_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig1B_SharedPrivateNovelKnown_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=shared_private_novel_known_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#This plot was skipped in the supplement due to partial redundancy with Fig. 1B
#Plot proportion of shared vs. private variants per region:
shared_private <- shared_private_counts %>%
   group_by(Region, VariantType, VariantSubset) %>%
   summarize(Count=sum(Count)) %>%
   pivot_wider(id_cols=c(Region, VariantType),
               names_from="VariantSubset",
               values_from="Count") %>%
   mutate(frac_shared=shared/all,
          frac_private=private/all) %>%
   select(-c(shared, private, all)) %>%
   pivot_longer(cols=c(frac_shared, frac_private),
                names_to="VariantSubset",
                values_to="Fraction",
                names_prefix="frac_") %>%
   mutate(VariantSubset=factor(VariantSubset,
                               levels=c("private", "shared"),
                               labels=c("Population-specific", "Shared among populations")),
          Percentage=str_c(signif(Fraction*100, digits=3), "%"),
          Region=factor(Region,
                        levels=c("AFR", "OCN", "EAS", "CSA",
                                 "EUR", "MDE", "AMR", "ISEA"))) %>%
   filter(VariantType == "SNP") %>%
   ggplot(aes(x=Region, y=Fraction, fill=VariantSubset)) +
      geom_col(position="stack") +
      geom_text(aes(label=Percentage), position=position_fill(vjust=0.5), col="white") +
      theme_bw() +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(labels=function(x) {paste0(format(x*100), "%")},
                         expand=c(0,0)) +
      scale_fill_manual(values=c("#9d1e25", "#92c4e8")) +
      guides(fill=guide_legend(nrow=2)) +
      labs(x="",
           y="",
           fill="") +
      theme(axis.text.x=element_text(size=12, face="bold"),
            axis.text.y=element_text(face="bold"),
            legend.position="bottom")
ggsave(paste0('PIBv1_SharedPrivateVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=shared_private,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_SharedPrivateVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=shared_private,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_SharedPrivateVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=shared_private,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#This was an old plot we used in a grant proposal, but was skipped
#Plot proportion novel vs. known region-private variants as pie charts:
private_novel_known <- shared_private_counts %>%
   filter(VariantSubset == "private") %>%
   group_by(Region, VariantType, Group) %>%
   summarize(Count=sum(Count)) %>%
   pivot_wider(id_cols=c(Region, VariantType),
               names_from="Group",
               values_from="Count") %>%
   mutate(frac_novel=Novel/(Novel+Known),
          frac_known=Known/(Novel+Known)) %>%
   select(-c(Novel, Known)) %>%
   pivot_longer(cols=c(frac_novel, frac_known),
                names_to="Group",
                values_to="Fraction",
                names_prefix="frac_") %>%
   mutate(Group=factor(Group,
                       levels=c("known", "novel"),
                       labels=c("Population-specific known", "Population-specific novel")),
          Percentage=str_c(signif(Fraction*100, digits=3), "%"),
          Region=factor(Region,
                        levels=c("AFR", "OCN", "EAS", "CSA",
                                 "EUR", "MDE", "AMR", "ISEA"))) %>%
   filter(VariantType == "SNP") %>%
   ggplot(aes(x=factor(1), y=Fraction, fill=Group)) +
      geom_col(position="stack") +
      theme_void() +
      coord_polar(theta="y") +
      facet_wrap(~ Region, nrow=1) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_manual(values=c("#f28d8f", "#bf3b36")) +
      guides(fill=guide_legend(nrow=2)) +
      labs(x="",
           y="",
           fill="") +
      theme(legend.position="bottom",
            axis.text=element_blank())
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=private_novel_known,
       width=16.0,
       height=6.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=private_novel_known,
       width=16.0,
       height=6.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=private_novel_known,
       width=16.0,
       height=6.0,
       units="cm",
       dpi=500)

#Supplementary figure S12 in the PIBv1 manuscript:
#Plot proportion novel vs. known variants per region, faceted by variant type:
novel_known <- allele_counts %>%
   filter(VariantSubset == "all",
          VariantType %in% c("SNP", "INS", "DEL")) %>%
   pivot_wider(id_cols=c(Region, VariantType),
               names_from="Group",
               values_from="Count") %>%
   mutate(frac_Novel=Novel/(Novel+Known),
          frac_Known=Known/(Novel+Known)) %>%
   select(-c(Novel, Known)) %>%
   pivot_longer(cols=c(frac_Known, frac_Novel),
                names_to="Found",
                values_to="Fraction",
                names_prefix="frac_") %>%
   mutate(Found=factor(Found,
                       levels=c("Novel", "Known")),
          Percentage=str_c(signif(Fraction*100, digits=3), "%"),
          Region=factor(Region,
                        levels=c("OCN", "CSA", "EAS", "MDE",
                                 "EUR", "ISEA", "AMR", "AFR")),
          VariantType=factor(VariantType,
                             levels=c("SNP", "INS", "DEL"))) %>%
   ggplot(aes(x=Region, y=Fraction, fill=Found)) +
      geom_col(position="stack") +
      geom_text(aes(label=Percentage), position=position_fill(vjust=0.25), col="black", fontface="bold") +
      theme_bw() +
      facet_wrap(~ VariantType, ncol=1) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_manual(values=c("#bf3b36", "#f28d8f")) +
      guides(fill=guide_legend(nrow=2)) +
      labs(x="",
           y="",
           fill="") +
      theme(axis.text.x=element_text(angle=90, size=12, face="bold"),
            axis.text.y=element_text(face="bold"),
            legend.position="bottom")
ggsave(paste0('PIBv1_NovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=novel_known,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_NovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=novel_known,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_NovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=novel_known,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S13 in the PIBv1 manuscript:
#Plot proportion novel vs. known region-private variants, faceted by variant type:
private_novel_known_bar <- shared_private_counts %>%
   filter(VariantSubset == "private") %>%
   group_by(Region, VariantType, Group) %>%
   summarize(Count=sum(Count)) %>%
   pivot_wider(id_cols=c(Region, VariantType),
               names_from="Group",
               values_from="Count") %>%
   mutate(frac_novel=Novel/(Novel+Known),
          frac_known=Known/(Novel+Known)) %>%
   select(-c(Novel, Known)) %>%
   pivot_longer(cols=c(frac_novel, frac_known),
                names_to="Group",
                values_to="Fraction",
                names_prefix="frac_") %>%
   mutate(Group=factor(Group,
                       levels=c("novel", "known"),
                       labels=c("Population-specific novel", "Population-specific known")),
          Percentage=str_c(signif(Fraction*100, digits=3), "%"),
          Region=factor(Region,
                        levels=c("OCN", "ISEA", "CSA", "MDE",
                                 "AMR", "EUR", "EAS", "AFR"))) %>%
   filter(VariantType %in% c("SNP", "INS", "DEL")) %>%
   mutate(VariantType=factor(VariantType,
                             levels=c("SNP", "INS", "DEL"))) %>%
   ggplot(aes(x=Region, y=Fraction, fill=Group)) +
      geom_col(position="stack") +
      geom_text(aes(label=Percentage), position=position_fill(vjust=0.25), col="black", fontface="bold") +
      theme_bw() +
      facet_wrap(~ VariantType, ncol=1) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_manual(values=c("#bf3b36", "#f28d8f")) +
      guides(fill=guide_legend(nrow=2)) +
      labs(x="",
           y="",
           fill="") +
      theme(axis.text.x=element_text(angle=90, size=12, face="bold"),
            axis.text.y=element_text(face="bold"),
            legend.position="bottom")
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=private_novel_known_bar,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=private_novel_known_bar,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_PrivateNovelKnownVariants_byRegion_byType_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=private_novel_known_bar,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Clean up:
rm(isec, allele_counts, shared_private_counts)
rm(shared_private_novel_known_plot, shared_private, private_novel_known, novel_known, private_novel_known_bar)
