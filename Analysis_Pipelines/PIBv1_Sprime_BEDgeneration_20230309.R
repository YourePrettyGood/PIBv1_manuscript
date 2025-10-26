#PIBv1 Sprime tract and core haplotype BED generation (including AI-specific BEDs)

#Helper functions:

#Load the libraries:
library(tidyverse)

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

#Load the data:
PIBv1_metadata <- read_tsv('../../Metadata/PIBv1_metadata_v0.3.tsv',
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

#Generate a map from population to geographic region:
pop_region_map <- bind_rows(PIBv1_metadata %>%
                               group_by(AnalysisGroup) %>%
                               summarize(Region=first(Region),
                                         IslandGroup=first(Island)),
                            lumped_pops)

#Tract frequencies, overlapping genes, and match rates:
target_tract_freqs <- read_tsv('allpops/PIBv1_Sprime_TargetPops_tract_freqs.tsv.gz',
                               col_types='ciiccnnni')
overlapping_tract_genes <- read_tsv('allpops/PIBv1_Sprime_gene_name_lists.tcsv.gz',
                                    col_types='ciicc')
target_corehap_freqs <- read_tsv('corehaps/PIBv1_Sprime_targetpop_corehap_freqs.tsv.gz',
                                 col_types='ciiccnnni')
overlapping_corehap_genes <- read_tsv('corehaps/PIBv1_Sprime_corehaps_gene_name_lists.tcsv.gz',
                                      col_types='ciicc')
match_rates <- read_tsv('allpops/PIBv1_perPop_Sprime_autosomal_match_rates.tsv.gz',
                        col_types='ciiciininininini')

#Annotate the tracts with geographic region and overlapping gene lists:
annotated_tract_freqs <- target_tract_freqs %>%
   left_join(pop_region_map,
             by=c("QueryPop"="AnalysisGroup")) %>%
   left_join(overlapping_tract_genes,
             by=c("Chromosome", "Start", "End", "TractID"))
annotated_corehap_freqs <- target_corehap_freqs %>%
   left_join(pop_region_map,
             by=c("QueryPop"="AnalysisGroup")) %>%
   left_join(overlapping_corehap_genes,
             by=c("Chromosome", "Start", "End", "TractID"))

#Identify the adaptive introgression threshold frequencies:
#For S' tracts, we apply a match rate filter:
# Ngood >= 30, Dgood >= 30, Nmatchrate >= 0.3 | Dmatchrate >= 0.3
AI_tracts <- annotated_tract_freqs %>%
   left_join(match_rates, by=c("Chromosome"="CHROM", "TractID"="TractID")) %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate >= 0.3 | Dmatchrate >= 0.3) %>%
   group_by(QueryPop) %>%
   mutate(freq95=quantile(MedianAF, 0.95)) %>%
   mutate(Significant=MedianAF>=freq95)

AI_corehaps <- annotated_corehap_freqs %>%
   group_by(QueryPop) %>%
   mutate(freq95=quantile(MedianAF, 0.95)) %>%
   mutate(Significant=MedianAF>=freq95)

#Generate the BED files of Sprime tracts and core haplotypes:
#All populations, S' tracts:
AI_tracts %>%
   ungroup() %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             Denisovan=Dmatchrate,
             Neandertal=Nmatchrate,
             AltaiNea=Amatchrate,
             VindijaNea=Vmatchrate,
             ChagyrskayaNea=Cmatchrate,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_tracts_MRfilter_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Oceania excluding Vanuatu, S' tracts:
AI_tracts %>%
   ungroup() %>%
   filter(Region == "OCN", IslandGroup != "Vanuatu") %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             Denisovan=Dmatchrate,
             Neandertal=Nmatchrate,
             AltaiNea=Amatchrate,
             VindijaNea=Vmatchrate,
             ChagyrskayaNea=Cmatchrate,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_tracts_OCNnoVanuatu_MRfilter_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#All populations, core haplotypes:
AI_corehaps %>%
   ungroup() %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_corehaps_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Oceania excluding Vanuatu, core haplotypes:
AI_corehaps %>%
   ungroup() %>%
   filter(Region == "OCN", IslandGroup != "Vanuatu") %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_corehaps_OCNnoVanuatu_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Generate the BED files of adaptive introgression candidates:
#All populations, S' tracts:
AI_tracts %>%
   ungroup() %>%
   filter(Significant) %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             Denisovan=Dmatchrate,
             Neandertal=Nmatchrate,
             AltaiNea=Amatchrate,
             VindijaNea=Vmatchrate,
             ChagyrskayaNea=Cmatchrate,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_putative_adaptive_introgressed_tracts_MRfilterBeforeTFfilter_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Oceania excluding Vanuatu, S' tracts:
AI_tracts %>%
   ungroup() %>%
   filter(Region == "OCN", IslandGroup != "Vanuatu") %>%
   filter(Significant) %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             Denisovan=Dmatchrate,
             Neandertal=Nmatchrate,
             AltaiNea=Amatchrate,
             VindijaNea=Vmatchrate,
             ChagyrskayaNea=Cmatchrate,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_putative_adaptive_introgressed_tracts_OCNnoVanuatu_MRfilterBeforeTFfilter_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#All populations, core haplotypes:
AI_corehaps %>%
   ungroup() %>%
   filter(Significant) %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_putative_adaptive_introgressed_corehaps_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Oceania excluding Vanuatu, core haplotypes:
AI_corehaps %>%
   ungroup() %>%
   filter(Region == "OCN", IslandGroup != "Vanuatu") %>%
   filter(Significant) %>%
   transmute(`#Chromosome`=Chromosome,
             Start=Start,
             End=End,
             SprimePopulation=QueryPop,
             TractID=TractID,
             NumSprimeSites=NumSites,
             TractFreq=MedianAF,
             MinAF=MinAF,
             MaxAF=MaxAF,
             OverlappingGenes=OverlappingGenes) %>%
   write_tsv('PIBv1_Sprime_putative_adaptive_introgressed_corehaps_OCNnoVanuatu_wGeneList.bed',
             col_names=TRUE,
             quote_escape="none")

#Clean up:
rm(region_names, region_codes, region_pretty)
rm(PIBv1_metadata, lumped_pops, pop_region_map)
rm(target_tract_freqs, overlapping_tract_genes, target_corehap_freqs, overlapping_corehap_genes)
rm(match_rates)
rm(annotated_tract_freqs, annotated_corehap_freqs)
rm(AI_tracts, AI_corehaps)
