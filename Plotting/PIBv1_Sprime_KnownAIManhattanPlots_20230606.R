#PIBv1 Sprime adaptive introgression plotting

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

#Load the libraries:
library(tidyverse)
library(ggrepel)
library(cowplot)

#Get to the directory:
pushd('[path redacted]/PIBv1_results/Sprime')

#Set the region codes:
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

#Load the metadata file:
PIBv1_metadata <- read_tsv('[path redacted]/Metadata/PIBv1_metadata_v0.3.tsv',
                           col_types='ccccclcccnncccclc',
                           na=c("NA", "")) %>%
   mutate(Region=factor(Region,
                        levels=region_names,
                        labels=region_codes))

lumped_pops <- data.frame(AnalysisGroup=c("Ambae,Maewo", "Bellona,Rennell",
                                          "Goroka,Sepik", "Nakanai,Mangseng"),
                          Region=c("OCN", "OCN",
                                   "OCN", "OCN"),
                          IslandGroup=c("Vanuatu", "Solomon Islands",
                                        "New Guinea", "New Britain"))

pop_region_map <- bind_rows(PIBv1_metadata %>%
                               group_by(AnalysisGroup) %>%
                               summarize(Region=first(Region),
                                         IslandGroup=first(Island)),
                            lumped_pops)

#Set up the data.frame for the background of the Manhattan plot and
# the offsets for tract positions:
hs37d5_fai <- read_tsv('[path redacted]/hs37d5.fa.fai',
                       col_types='cnnii',
                       col_names=c("Chromosome", "Length", "Offset", "Wrap", "WrapLF"))
manhattan_df <- hs37d5_fai %>%
   transmute(Chromosome=Chromosome,
             ChromStart=c(0, cumsum(Length))[-(n()+1)]+1,
             ChromEnd=cumsum(Length),
             Background=case_when(row_number() %% 2 == 1 ~ "black",
                                  TRUE ~ "grey30")) %>%
   mutate(Center=(ChromStart+ChromEnd)/2)

#Load in the core haplotype frequencies in the Sprime target populations:
target_corehap_freqs <- read_tsv('corehaps/PIBv1_Sprime_targetpop_corehap_freqs.tsv.gz',
                                 col_types='ciiccnnni')
#Load the set of protein-coding genes for trimming the gene lists:
coding_genes <- read_tsv('corehaps/gencode.v38lift37.annotation.codingGENES_sorted.bed',
                         col_types='cnncc',
                         col_names=c("CHROM", "BEDStart", "BEDEnd", "ENSG", "GeneName"),
                         skip=0)
#Load in the genes overlapping these core haplotypes:
overlapping_corehap_genes <- read_tsv('corehaps/PIBv1_Sprime_corehaps_gene_name_lists.tcsv.gz',
                                      col_types='ciicc') %>%
   mutate(PCGs=map_chr(OverlappingGenes,
                       .f=function(x, keepgenes) {
                          paste(intersect(str_split(x, ",")[[1]], keepgenes),
                                collapse=",")
                       },
                       keepgenes=coding_genes$GeneName))

#Annotate core haplotype frequencies with metadata and the overlapping genes:
annotated_corehap_freqs <- target_corehap_freqs %>%
   left_join(pop_region_map,
             by=c("QueryPop"="AnalysisGroup")) %>%
   left_join(overlapping_corehap_genes,
             by=c("Chromosome", "Start", "End", "TractID"))

#Set up the plotting data.frame for Sprime core haplotypes:
corehap_plotting_df <- annotated_corehap_freqs %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   group_by(QueryPop) %>%
   mutate(freq95=quantile(MedianAF, 0.95)) %>%
   mutate(Significant=MedianAF>=freq95,
          GeneList=case_when(MedianAF>=freq95 ~ str_wrap(str_replace_all(PCGs, ",", ", "), width=40),
                             TRUE ~ NA_character_)) %>%
   left_join(manhattan_df %>%
                select(Chromosome, ChromStart),
             by="Chromosome") %>%
   transmute(QueryPop=QueryPop,
             Chromosome=Chromosome,
             ManhattanPos=ChromStart+Start, TractID=TractID, Region=Region,
             MedianAF=MedianAF,
             Significant=Significant,
             GeneList=GeneList,
             freq95=freq95) %>%
   separate(col=TractID,
            into=c(NA, NA, "Origin"),
            sep="_",
            remove=FALSE) %>%
   mutate(Origin=str_replace(Origin, "[.][0-9]+$", ""))

#Supplementary figure S64 in the PIBv1 manuscript:
#Show AI candidates from the literature with hits in our analysis in Europeans:
EUR_known_AI <- corehap_plotting_df %>%
   filter(Region == "EUR") %>%
   group_by(GeneList) %>%
   mutate(PointSize=case_when(Significant & str_detect(GeneList, "(BNC2|OAS[123]|WDR88|GPATCH1|CNTNAP4|EYS|SGCZ|PRKCQ)") ~ 8,
                              TRUE ~ 1)) %>%
   mutate(GeneList=case_when(MedianAF == max(MedianAF) & Significant & str_detect(GeneList, "(BNC2|OAS[123]|WDR88|GPATCH1|CNTNAP4|EYS|SGCZ|PRKCQ)") ~ GeneList,
                             TRUE ~ NA_character_)) %>%
   ungroup() %>%
   mutate(SignificantPop=case_when(Significant ~ QueryPop,
                                TRUE ~ NA_character_)) %>%
   ggplot(aes(x=ManhattanPos, y=MedianAF)) +
      geom_rect(aes(xmin=ChromStart, xmax=ChromEnd, ymin=0, ymax=1, fill=Background),
                data=manhattan_df %>%
                   filter(Chromosome %in% as.character(1:22)),
                inherit.aes=FALSE,
                alpha=0.1,
                colour=NA) +
      geom_point(aes(colour=SignificantPop, size=PointSize, alpha=Significant)) +
      geom_text_repel(aes(label=GeneList),
                      size=3,
                      segment.size=0.2,
                      min.segment.length=0,
                      nudge_x=-50000000,
                      nudge_y=0.2,
                      force=10,
                      max.iter=1e6) +
      theme_classic() +
      scale_x_continuous(labels=1:22,
                         breaks=(manhattan_df %>%
                                    filter(Chromosome %in% as.character(1:22)))$Center,
                         expand=c(0,0)) +
      scale_y_continuous(limits=c(0, 1),
                         expand=c(0, 0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_colour_viridis_d(na.value="grey40", begin=0.1, end=0.9, name="") +
      scale_fill_manual(values=c("black", "grey30")) +
      scale_size_area(max_size=2) +
      scale_alpha_discrete(na.value=0.1, breaks=NULL) +
      guides(colour=guide_legend(override.aes=list(size=3),
                                 byrow=TRUE),
             fill="none",
             alpha="none",
             size="none") +
      labs(x="Chromosome",
           y="S' core haplotype frequency",
           title="EUR S' AI core haplotypes (known)") +
      theme(legend.position="top",
            legend.spacing.x=unit(0, "cm"),
            legend.spacing.y=unit(0, "cm"),
            legend.margin=margin(t=-0.25, r=0, b=-0.25, l=0, unit="cm"))
ggsave(paste0('PIBv1_FigS_EUR_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=EUR_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS_EUR_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=EUR_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS_EUR_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=EUR_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Supplementary figure S65 in the PIBv1 manuscript:
#Show AI candidates from the literature with hits in our analysis in East Asians:
EAS_known_AI <- corehap_plotting_df %>%
   filter(Region == "EAS") %>%
   group_by(GeneList) %>%
   mutate(PointSize=case_when(Significant & str_detect(GeneList, "(POU2F3|SIPA1L2|TLR[16][0]?|ZNF169|TXN(,|$)|TBC1D1|CHMP1A|SLC35F3)") ~ 8,
                              TRUE ~ 1)) %>%
   mutate(GeneList=case_when(MedianAF == max(MedianAF) & Significant & str_detect(GeneList, "(POU2F3|SIPA1L2|TLR[16][0]?|ZNF169|TXN(,|$)|TBC1D1|CHMP1A|SLC35F3)") ~ GeneList,
                             TRUE ~ NA_character_)) %>%
   ungroup() %>%
   mutate(SignificantPop=case_when(Significant ~ QueryPop,
                                   TRUE ~ NA_character_)) %>%
   ggplot(aes(x=ManhattanPos, y=MedianAF)) +
      geom_rect(aes(xmin=ChromStart, xmax=ChromEnd, ymin=0, ymax=1, fill=Background),
                data=manhattan_df %>%
                   filter(Chromosome %in% as.character(1:22)),
                inherit.aes=FALSE,
                alpha=0.1,
                colour=NA) +
      geom_point(aes(colour=SignificantPop, size=PointSize, alpha=Significant)) +
      geom_text_repel(aes(label=GeneList),
                      size=2,
                      segment.size=0.2,
                      min.segment.length=0,
                      nudge_x=-50000000,
                      nudge_y=0.2,
                      force=10,
                      max.iter=1e6) +
      theme_classic() +
      scale_x_continuous(labels=1:22,
                         breaks=(manhattan_df %>%
                                    filter(Chromosome %in% as.character(1:22)))$Center,
                         expand=c(0,0)) +
      scale_y_continuous(limits=c(0, 1),
                         expand=c(0, 0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_colour_viridis_d(na.value="grey40", begin=0.1, end=0.9, name="") +
      scale_fill_manual(values=c("black", "grey30")) +
      scale_size_area(max_size=2) +
      scale_alpha_discrete(na.value=0.1, breaks=NULL) +
      guides(colour=guide_legend(override.aes=list(size=1.5),
                                 byrow=TRUE,
                                 nrow=3),
             fill="none",
             alpha="none",
             size="none") +
      labs(x="Chromosome",
           y="S' core haplotype frequency",
           title="EAS S' AI core haplotypes (known)") +
      theme(legend.position="top",
            legend.spacing.x=unit(0, "cm"),
            legend.spacing.y=unit(0, "cm"),
            legend.key.size=unit(0.5, "cm"),
            legend.margin=margin(t=-0.25, r=0, b=-0.25, l=0, unit="cm"))
ggsave(paste0('PIBv1_FigS_EAS_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=EAS_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS_EAS_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=EAS_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS_EAS_AI_knowncandidates_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=EAS_known_AI,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Table of the effect of trimming very short and very long core haplotypes:
annotated_corehap_freqs %>%
   separate(col=TractID, into=c(NA, NA, "Origin"), sep="_", remove=FALSE) %>%
   mutate(Origin=str_replace(Origin, "[.][0-9]+$", "")) %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   mutate(TractPhysLen=End-Start) %>%
   group_by(QueryPop) %>%
   mutate(freq95=quantile(MedianAF, 0.95),
          NumTracts=n()) %>%
   filter(TractPhysLen >= 100, TractPhysLen <= 5000000) %>%
   mutate(trimmedfreq95=quantile(MedianAF, 0.95),
          NumTrimmedTracts=n()) %>%
   summarize(freq95=first(freq95),
             NumTracts=last(NumTracts),
             trimmedfreq95=first(trimmedfreq95),
             NumTrimmedTracts=last(NumTrimmedTracts)) %>%
   mutate(threshdiff=trimmedfreq95-freq95,
          tractdiff=NumTracts-NumTrimmedTracts) %>%
   arrange(desc(threshdiff)) %>%
   View()

#Plotting high-frequency core haplotypes per population stratified by archaic origin:
annotated_corehap_freqs %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   group_by(QueryPop) %>%
   mutate(freq95=quantile(MedianAF, 0.95)) %>%
   mutate(Significant=MedianAF>=freq95) %>%
   filter(Significant) %>%
   group_by(Region, IslandGroup, QueryPop, Origin) %>%
   summarize(NumAIHits=n()) %>%
   mutate(QueryPop=factor(QueryPop,
                          levels=unique(arrange(., Region, IslandGroup, QueryPop)$QueryPop))) %>%
   ggplot(aes(x=QueryPop, y=NumAIHits, colour=Origin, fill=Origin)) +
      geom_col(position="stack") +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      scale_fill_brewer(palette="Set2") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
            legend.position="top")

#Clean up:
rm(region_names, region_codes, region_pretty, PIBv1_metadata, lumped_pops, pop_region_map)
rm(hs37d5_fai, manhattan_df)
rm(target_corehap_freqs, coding_genes, overlapping_corehap_genes)
rm(annotated_corehap_freqs, corehap_plotting_df)
rm(EUR_known_AI, EAS_known_AI)
