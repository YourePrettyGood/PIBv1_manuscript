#PIBv1 Sprime per-individual and per-haplotype tract length plots:

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
pushd('[path redacted]/PIBv1_results/Sprime/')

#Helper functions:
#Wrapper for overplotting the Gaussian mixture model density curves
# on a 1-D histogram:
geom_gmm_density <- function(fits_df, binwidth) {
   mapply(function(pop, alpha, mean, sd, n, binwidth) {
      stat_function(data=fits_df %>%
                       filter(Population == pop, mu == mean),
                    fun=function(x) {
                       alpha*dnorm(x, mean=mean, sd=sd)*n*binwidth
                    })
   },
   pop=fits_df$Population,
   alpha=fits_df$alpha,
   mean=fits_df$mu,
   sd=fits_df$stdev,
   n=fits_df$TractCount,
   binwidth=binwidth)
}
#Wrapper for overplotting the beta mixture model density curves
# on a 1-D histogram:
geom_bmm_density <- function(fits_df, binwidth) {
   mapply(function(pop, alpha, shape1, shape2, n, binwidth) {
      stat_function(data=fits_df %>%
                       filter(Population == pop, shape1 == shape1),
                    fun=function(x) {
                       alpha*dbeta(x, shape1=shape1, shape2=shape2)*n*binwidth
                    })
   },
   pop=fits_df$Population,
   alpha=fits_df$alpha,
   shape1=fits_df$shape1,
   shape2=fits_df$shape2,
   n=fits_df$TractCount,
   binwidth=binwidth)
}

#Plotting choices:
archaic_origins <- c("Ambiguous", "Denisovan", "Neanderthal")
archaic_origin_colours <- c("#B588AF", "#41BBEC", "#E84B42")
archaic_origin_colours_browning <- c(archaic_origin_colours, "grey70")
#Match rate boundaries for example plot:
archaic_boxes_pfr <- data.frame(xmin=c(0.31, 0, 0.3),
                                xmax=c(1.0, 0.29, 1.0),
                                ymin=c(0.31, 0.3, 0),
                                ymax=c(1.0, 1.0, 0.29),
                                TractOrigin=archaic_origins) %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins))
archaic_boxes_browning <- data.frame(xmin=c(0.31, 0, 0.6, 0.30),
                                     xmax=c(1.0, 0.29, 1.0, 0.57),
                                     ymin=c(0.41, 0.4, 0, 0),
                                     ymax=c(1.0, 1.0, 0.38, 0.38),
                                     TractOrigin=c(archaic_origins, "NA")) %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=c(archaic_origins, "NA")))
#Region colours:
region_colours <- c("#BE1E2D", #AFR
                    "#DB7A24", #MDE
                    "#81B758", #EUR
                    "#F9B817", #AMR
                    "#E292A8", #CSA
                    "#842B6E", #EAS
                    "#5050FF", #ISEA
                    "#3CC0C4") #OCN

#Colours and shapes for Sprime vs. RFMix plot from Chang and Audrey:
Fig2C_shapes_colours <- data.frame(Population=c("Agta", "Cebuano", "Rampasasa", #ISEA
																"Kove", "Nakanai-Mangseng", #West New Britain
																"Mamusi", "Ata", "Melamela", #West New Britain
																"Baining-Kagat", "Baining-Mali", #East New Britain
																"Lavongai-Mussau", #Mussau and New Hanover
																"Nailik-Notsi-Tigak", #New Ireland
																"Nasioi", "Saposa", #Bougainville
																"Vella Lavella", "Malaita", #non-Polynesian outliers
																"Bellona-Rennell", "Tikopia", #Polynesian outliers
																"Santa Cruz"), #non-Polynesian outliers
											  PopLevels=c("Agta", "Cebuano", "Rampasasa", #ISEA
											  				  "Kove", "Nakanai-Mangseng", #West New Britain
											  				  "Mamusi", "Ata", "Melamela", #West New Britain
											  				  "Baining-Kagat", "Baining-Mali", #East New Britain
											  				  "Lavongai-Mussau", #Mussau and New Hanover
											  				  "Nailik-Notsi-Tigak", #New Ireland
											  				  "Nasioi", "Saposa", #Bougainville
											  				  "Vella Lavella", "Malaita", #non-Polynesian outliers
											  				  "Bellona-Rennell", "Tikopia", #Polynesian outliers
											  				  "Santa Cruz"), #non-Polynesian outliers
											  colours=c(rep("#0C775E", 3), #ISEA
															rep("#FFE59D", 2), #West New Britain
															rep("#FFB531", 3), #West New Britain
															rep("#FF7575", 2), #East New Britain
															"#BA9EF1", #Mussau and New Hanover
															"#5897D2", #New Ireland
															rep("#7DDFEC", 2), #Bougainville
															rep("#83CF8E", 2), #non-Polynesian outliers
															rep("#C3557A", 2), #Polynesian outliers
															"#83CF8E"), #non-Polynesian outliers
											  shapes=c(16, 15, 18, #ISEA
											  			  16, 17, 15, 18, 16, #West New Britain
											  			  15, 16, #East New Britain
											  			  15, #Mussau and New Hanover
											  			  16, #New Ireland
											  			  15, 16, #Bougainville
											  			  16, 17, 15, 16, 18)) #Solomon Islands

#Load the data:
#Sample metadata file:
PIBv1_metadata <- read_tsv('../../Metadata/PIBv1_metadata_v0.3.tsv',
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

#Generate the population-to-region map:
pop_region_map <- bind_rows(PIBv1_metadata %>%
                               group_by(AnalysisGroup) %>%
                               summarize(Region=first(Region),
                                         IslandGroup=first(Island)),
                            data.frame(AnalysisGroup=c("Ambae,Maewo", "Bellona,Rennell",
                                                       "Goroka,Sepik", "Nakanai,Mangseng"),
                                       Region=c("Oceania", "Oceania",
                                                "Oceania", "Oceania"),
                                       IslandGroup=c("Vanuatu", "Solomon Islands",
                                                     "New Guinea", "New Britain"))) %>%
   mutate(IslandGroup=factor(IslandGroup,
                             levels=c("Taiwan", "Philippines", "Flores",
                                      "New Guinea", "New Britain", "New Hanover",
                                      "New Ireland", "Bougainville", "Solomon Islands",
                                      "Vanuatu")))

#Set up a data.frame of the populations to select for the figure:
curated_pops <- data.frame(name=c("French", "Maya",
                                  "Hazara", "Dai",
                                  "Rampasasa", "Lavongai-Mussau",
                                  "Ata", "Mamusi",
                                  "Baining-Kagat", "Baining-Mali",
                                  "Goroka", "Sepik"),
                           label=c("EUR-French", "AMR-Maya",
                                   "CSA-Hazara", "EAS-Dai",
                                   "ISEA-Rampasasa", "OCN-Lavongai-Mussau",
                                   "OCN-Ata", "OCN-Mamusi",
                                   "OCN-Baining-Kagat", "OCN-Baining-Mali",
                                   "OCN-Goroka", "OCN-Sepik"))

#Load the match rates and projected tract lengths:
match_rates <- read_tsv('PIBv1_perPop_Sprime_autosomal_match_rates.tsv.gz',
                        col_types='ciiciininininini')
arc_tract_lengths_maxgap0 <- read_tsv('PIBv1_Sprime_perChrom_perHaplotype_perPop_tract_lengths_maxgap0_sumtractsfix.tsv.gz',
                                        col_types='cccciii')
#Load the genotype-level tract projections:
arc_tract_lengths <- read_tsv('PIBv1_Sprime_perChrom_perIndiv_perPop_tract_lengths.tsv.gz',
                              col_types='cccciii')

#Load the BMM and GMM results:
load('PIBv1_Denisovan_gmms_bmms_k1thru10_noinputs_Rv2.Rdata')

#Match rate plots:
#Annotate tracts and apply a baseline filter on ascertainabile sites:
filtered_match_rates <- match_rates %>%
   separate(col=TractID, into=c("Population"), sep="_", remove=FALSE, extra="drop") %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 10 | Dgood >= 10)

#Tract length plots:
#Annotate and filter tracts with baseline filter on ascertainable sites and
# minimum match rates:
filtered_tracts_lengths <- arc_tract_lengths_maxgap0 %>%
   left_join(match_rates, by=c("CHROM", "TractID")) %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(RegionCode=factor(Region,
                            levels=c("Africa", "Middle East",
                                     "Europe", "America",
                                     "Central/South Asia", "East Asia",
                                     "Island Southeast Asia", "Oceania"),
                            labels=c("AFR", "MDE",
                                     "EUR", "AMR",
                                     "CSA", "EAS",
                                     "ISEA", "OCN")),
   		 Population=str_replace(Population, ",", "-")) %>%
	mutate(Population=case_when(Population %in% c("Santa-Cruz", "Vella-Lavella") ~ str_replace(Population, "-", " "),
										 TRUE ~ Population)) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., RegionCode, Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 30, Dgood >= 30, Nmatchrate >= 0.3 | Dmatchrate >= 0.3)

#Archaic tract homozygosity plots:
#Get total number of Sprime sites per tract from match rates and filter
# tracts:
filtered_tract_homozygosity <- arc_tract_lengths %>%
   left_join(match_rates, by=c("CHROM", "TractID")) %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(RegionCode=factor(Region,
                            levels=c("Africa", "Middle East",
                                     "Europe", "America",
                                     "Central/South Asia", "East Asia",
                                     "Island Southeast Asia", "Oceania"),
                            labels=c("AFR", "MDE",
                                     "EUR", "AMR",
                                     "CSA", "EAS",
                                     "ISEA", "OCN"))) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., RegionCode, Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 30, Dgood >= 30, Nmatchrate >= 0.3 | Dmatchrate >= 0.3) %>%
#   mutate(TractOrigin=case_when(Ngood >= 30 & Dgood >= 30 & Nmatchrate >= 0.3 & Dmatchrate <= 0.3 ~ "Neanderthal",
#                                Ngood >= 30 & Dgood >= 30 & Nmatchrate <= 0.3 & Dmatchrate >= 0.3 ~ "Denisovan",
#                                Ngood >= 30 & Dgood >= 30 & Nmatchrate > 0.3 & Dmatchrate > 0.3 ~ "Ambiguous",
#                                TRUE ~ NA_character_)) %>%
#   group_by(RegionCode, Region, IslandGroup, Population, TractOrigin, Sample) %>%
   group_by(RegionCode, Region, IslandGroup, Population, Sample) %>%
   summarize(HetSites=sum(ARCMOD),
             HomSites=sum(ARCARC),
             TotalSites=sum(SNPLEN)) %>%
   mutate(Heterozygous=HetSites*100/TotalSites,
          HomozygousModern=100-(HetSites+HomSites)*100/TotalSites,
          HomozygousArchaic=HomSites*100/TotalSites)

#Determine total tract length per individual regardless of origin:
total_perindiv_lengths <- filtered_tracts_lengths %>%
   mutate(TractOrigin="Total") %>%
   group_by(RegionCode, Region, IslandGroup, Population, Haplotype, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(TRACTLEN)) %>%
   separate(Haplotype, into=c("Sample"), remove=FALSE, extra="drop") %>%
   group_by(RegionCode, Region, IslandGroup, Population, Sample, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(ArchaicTractTotal))

#Also determine total tract length per individual partitioned by tract origin:
#Using rough criteria here.
rough_classified_perindiv_lengths <- filtered_tracts_lengths %>%
   mutate(TractOrigin=case_when(Ngood >= 30 & Dgood >= 30 & Nmatchrate >= 0.3 & Dmatchrate <= 0.3 ~ "Neanderthal",
                                Ngood >= 30 & Dgood >= 30 & Nmatchrate <= 0.3 & Dmatchrate >= 0.3 ~ "Denisovan",
                                Ngood >= 30 & Dgood >= 30 & Nmatchrate > 0.3 & Dmatchrate > 0.3 ~ "Ambiguous",
                                TRUE ~ NA_character_)) %>%
   group_by(RegionCode, Region, IslandGroup, Population, Haplotype, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(TRACTLEN)) %>%
   separate(Haplotype, into=c("Sample"), remove=FALSE, extra="drop") %>%
   group_by(RegionCode, Region, IslandGroup, Population, Sample, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(ArchaicTractTotal))

#Using Browning et al. 2018 criteria here:
browning_classified_perindiv_lengths <- filtered_tracts_lengths %>%
   mutate(TractOrigin=case_when(Agood >= 30 & Dgood >= 30 & Amatchrate >= 0.6 & Dmatchrate <= 0.4 ~ "Neanderthal",
                                Agood >= 30 & Dgood >= 30 & Amatchrate <= 0.3 & Dmatchrate >= 0.4 ~ "Denisovan",
                                Agood >= 30 & Dgood >= 30 & Amatchrate > 0.3 & Dmatchrate > 0.4 ~ "Ambiguous",
                                TRUE ~ NA_character_)) %>%
   group_by(RegionCode, Region, IslandGroup, Population, Haplotype, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(TRACTLEN)) %>%
   separate(Haplotype, into=c("Sample"), remove=FALSE, extra="drop") %>%
   group_by(RegionCode, Region, IslandGroup, Population, Sample, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(ArchaicTractTotal))

#Figure 2A:
tractlen_byorigin_byregion_stackedbar <- rough_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins)) %>%
   filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
   filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
   group_by(Region, RegionCode, TractOrigin) %>%
   summarize(MinLen=min(ArchaicTractTotal),
             Q1Len=quantile(ArchaicTractTotal, probs=c(0.25)),
             MedianLen=median(ArchaicTractTotal),
             Q3Len=quantile(ArchaicTractTotal, probs=c(0.75)),
             MaxLen=max(ArchaicTractTotal),
             MeanLen=mean(ArchaicTractTotal),
             SdLen=sd(ArchaicTractTotal)) %>%
   ggplot(aes(x=RegionCode, y=MedianLen/1e6, ymin=MaxLen/1e6, ymax=MaxLen/1e6, fill=TractOrigin)) +
      geom_col(position="stack", width=0.7, col="black", size=0.25) +
      theme_classic(base_line_size=0.25,
      				  base_rect_size=0.5) +
	   scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(0, 200)) +
      scale_fill_manual(values=archaic_origin_colours) +
      guides(fill="none") +
      labs(x="",
           y="",
           fill="") +
      theme(axis.line=element_line(linewidth=0.05),
      		axis.text.x=element_text(angle=90, hjust=1,
      										 vjust=0.5, face="bold",
      										 size=12))
tractlen_byorigin_OCN_stackedbar <- rough_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins)) %>%
   filter(!is.na(IslandGroup),
          IslandGroup != "Vanuatu",
   		 RegionCode == "OCN") %>%
   filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
   group_by(Region, RegionCode, IslandGroup, Population, TractOrigin) %>%
   summarize(MinLen=min(ArchaicTractTotal),
             Q1Len=quantile(ArchaicTractTotal, probs=c(0.25)),
             MedianLen=median(ArchaicTractTotal),
             Q3Len=quantile(ArchaicTractTotal, probs=c(0.75)),
             MaxLen=max(ArchaicTractTotal),
             MeanLen=mean(ArchaicTractTotal),
             SdLen=sd(ArchaicTractTotal)) %>%
   ggplot(aes(x=Population, y=MedianLen/1e6, ymin=MaxLen/1e6, ymax=MaxLen/1e6, fill=TractOrigin)) +
      geom_col(position="stack", width=0.7, col="black", size=0.25) +
      theme_classic(base_line_size=0.25,
      				  base_rect_size=0.5) +
	   scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0), limits=c(0, 200)) +
      scale_fill_manual(values=archaic_origin_colours) +
      labs(x="",
           y="",
           fill="") +
      theme(axis.line=element_line(linewidth=0.05),
      		axis.text.x=element_text(angle=90, hjust=1,
      										 vjust=0.5, face="bold",
      										 size=12),
      		panel.border=element_blank())
#Fig2A_init <- plot_grid(tractlen_byorigin_byregion_stackedbar, tractlen_byorigin_OCN_stackedbar,
#                   align='hv',
#                   axis='tblr',
#                   ncol=2,
#                   rel_widths=c(7, 18))
Fig2A_legend <- get_legend(tractlen_byorigin_OCN_stackedbar + theme(legend.position="top",
																						  legend.direction="vertical"))
Fig2A_vertalign <- align_plots(tractlen_byorigin_byregion_stackedbar,
										 tractlen_byorigin_OCN_stackedbar + theme(legend.position="none"),
										 align="v",
										 axis="l")
Fig2A_top <- plot_grid(Fig2A_vertalign[[1]], Fig2A_legend,
							  ncol=2,
							  rel_widths=c(1,1))
#Coordinates for drawing the axis label and zoom lines between panels:
arclen_label_coords <- list(x=0.02,
									 y=0.62)
ocnzoom <- list(x=c(0.085, 0.446, 0.487, 0.990),
					 y=c(0.625, 0.67, 0.75))
#Compose the plot and add the common y-axis label and the zoom-in lines:
Fig2A <- plot_grid(Fig2A_top, Fig2A_vertalign[[2]],
						 nrow=2,
						 rel_heights=c(1, 1.75)) +
	draw_text("Median Per-Individual Archaic Sequence (Mbp)",
				 angle=90,
				 x=arclen_label_coords$x,
				 y=arclen_label_coords$y) +
	draw_line(x=c(ocnzoom$x[2], ocnzoom$x[2], ocnzoom$x[1]),
				 y=c(ocnzoom$y[3], ocnzoom$y[2], ocnzoom$y[1]),
				 colour="black",
				 size=0.25) +
	draw_line(x=c(ocnzoom$x[3], ocnzoom$x[3], ocnzoom$x[4]),
				 y=c(ocnzoom$y[3], ocnzoom$y[2], ocnzoom$y[1]),
				 colour="black",
				 size=0.25) +
	draw_line(x=c(ocnzoom$x[2], ocnzoom$x[3]),
				 y=c(ocnzoom$y[2], ocnzoom$y[2]),
				 colour="black",
				 size=0.25)
ggsave(paste0('PIBv1_Fig2A_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=Fig2A,
       width=16.0,
       height=16.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig2A_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=Fig2A,
       width=16.0,
       height=16.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig2A_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=Fig2A,
       width=16.0,
       height=16.0,
       units="cm",
       dpi=500)

#Figure 2B was made by Stephen with separate code

#Actual Fig. 2C:
#Load Chang's original per-individual dataset:
sprime_vs_rfmix <- read_tsv('Chang_Sprime_vs_RFMix_Fig2C_perindiv_ignoreSprime_20250124.txt')
#Drop the Sprime proportions from Chang's dataset and replace with the
# updated maxgap=0 results:
#Remember to drop the redundant individuals (because I ran Sprime with
# and without lumping some populations).
sprime_vs_rfmix_adjusted <- rough_classified_perindiv_lengths %>%
	filter(!(Population %in% c("Nakanai", "Bellona"))) %>%
	mutate(ArchaicTractProportion=ArchaicTractTotal/(2*2684673005)) %>%
	inner_join(sprime_vs_rfmix %>%
				  	filter(!(population %in% c("Nakanai", "Bellona"))) %>%
				  	select(ID, EA, Papuan),
				  by=c("Sample"="ID"))
#Now make the plot:
Fig2C <- sprime_vs_rfmix_adjusted %>%
	mutate(Population=factor(Population,
									 levels=Fig2C_shapes_colours$PopLevels,
									 labels=Fig2C_shapes_colours$Population)) %>%
	filter(TractOrigin == "Denisovan") %>%
	group_by(Population) %>%
	summarize(DenisovanMean=mean(ArchaicTractProportion),
				 DenisovanSD=sd(ArchaicTractTotal)/(2*2684673005),
				 PapuanMean=mean(Papuan),
				 PapuanSD=sd(Papuan)) %>%
	ggplot(aes(x=PapuanMean, xmin=PapuanMean-PapuanSD, xmax=PapuanMean+PapuanSD,
				  y=DenisovanMean, ymin=DenisovanMean-DenisovanSD, ymax=DenisovanMean+DenisovanSD,
				  colour=Population, shape=Population)) +
	   geom_point(size=4) +
	   geom_errorbar() +
	   geom_errorbarh() +
	   geom_smooth(aes(x=PapuanMean, y=DenisovanMean),
	   				method="lm",
	   				inherit.aes=FALSE,
	   				se=FALSE,
	   				col="grey70",
	   				linetype=2) +
	   scale_x_continuous(limits=c(0, 1),
	   						 expand=c(0, 0.02),
	   						 labels=scales::percent_format(suffix="")) +
	   scale_y_continuous(limits=c(0, NA),
	   						 expand=c(0, 0),
	   						 labels=scales::percent_format(suffix="")) +
	   scale_colour_manual(values=Fig2C_shapes_colours$colours,
	   						  labels=Fig2C_shapes_colours$Population) +
	   scale_shape_manual(values=Fig2C_shapes_colours$shapes,
	   						 labels=Fig2C_shapes_colours$Population) +
	   theme_classic() +
	   labs(x="New Guinean ancestry (%)",
	   	  y="Denisovan ancestry (%)") +
	   theme(panel.border=element_rect(fill=NA, colour="black"))
ggsave(paste0('PIBv1_Fig2C_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
		 plot=Fig2C + theme(legend.position="none"),
		 width=14.0,
		 height=14.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_Fig2C_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.png'),
		 plot=Fig2C + theme(legend.position="none"),
		 width=14.0,
		 height=14.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_Fig2C_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
		 plot=Fig2C + theme(legend.position="none"),
		 width=14.0,
		 height=14.0,
		 units="cm",
		 dpi=500)
Fig2C_legend <- get_legend(Fig2C)

#Quick regressions for archaic ancestries against New Guinean ancestry:
sprime_vs_rfmix_adjusted %>%
	filter(TractOrigin == "Neanderthal") %>%
	lm(data=., ArchaicTractProportion ~ Papuan) %>%
	summary()
#Call:
#lm(formula = ArchaicTractProportion ~ Papuan, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0047014 -0.0005831  0.0001983  0.0008747  0.0029505
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.0092051  0.0002045   45.00   <2e-16 ***
#Papuan      -0.0002446  0.0003097   -0.79     0.43 
#Residual standard error: 0.001302 on 306 degrees of freedom
#Multiple R-squared:  0.002034,    Adjusted R-squared:  -0.001227
#F-statistic: 0.6237 on 1 and 306 DF,  p-value: 0.4303
sprime_vs_rfmix_adjusted %>%
	filter(TractOrigin == "Denisovan") %>%
	lm(data=., ArchaicTractProportion ~ Papuan) %>%
	summary()
#Call:
#lm(formula = ArchaicTractProportion ~ Papuan, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0037487 -0.0008045 -0.0001845  0.0004872  0.0050194
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.0017513  0.0002090   8.379 1.96e-15 ***
#Papuan       0.0070056  0.0003165  22.136  < 2e-16 ***
#Residual standard error: 0.00133 on 306 degrees of freedom
#Multiple R-squared:  0.6156,    Adjusted R-squared:  0.6143
#F-statistic: 490 on 1 and 306 DF,  p-value: < 2.2e-16
sprime_vs_rfmix_adjusted %>%
	filter(TractOrigin == "Ambiguous") %>%
	lm(data=., ArchaicTractProportion ~ Papuan) %>%
	summary()
#Call:
#lm(formula = ArchaicTractProportion ~ Papuan, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0058501 -0.0010643 -0.0001588  0.0010254  0.0044529
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.0062618  0.0002523   24.82   <2e-16 ***
#Papuan       0.0051924  0.0003820   13.59   <2e-16 ***
#Residual standard error: 0.001606 on 306 degrees of freedom
#Multiple R-squared:  0.3764,    Adjusted R-squared:  0.3744
#F-statistic: 184.7 on 1 and 306 DF,  p-value: < 2.2e-16

#Remove Agta and recompute regressions for DEN and AMB:
sprime_vs_rfmix_adjusted %>%
	filter(TractOrigin == "Denisovan",
			 Population != "Agta") %>%
	lm(data=., ArchaicTractProportion ~ Papuan) %>%
	summary()
#Call:
#lm(formula = ArchaicTractProportion ~ Papuan, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0036957 -0.0004261  0.0000018  0.0005260  0.0031805
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.0003708  0.0001680   2.207   0.0281 *
#Papuan       0.0088147  0.0002473  35.646   <2e-16 ***
#Residual standard error: 0.0009425 on 287 degrees of freedom
#Multiple R-squared:  0.8157,    Adjusted R-squared:  0.8151
#F-statistic: 1271 on 1 and 287 DF,  p-value: < 2.2e-16
sprime_vs_rfmix_adjusted %>%
	filter(TractOrigin == "Ambiguous",
			 Population != "Agta") %>%
	lm(data=., ArchaicTractProportion ~ Papuan) %>%
	summary()
#Call:
#lm(formula = ArchaicTractProportion ~ Papuan, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0058033 -0.0008459 -0.0001727  0.0008477  0.0031914
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.0050397  0.0002494   20.21   <2e-16 ***
#Papuan       0.0067940  0.0003671   18.51   <2e-16 ***
#Residual standard error: 0.001399 on 287 degrees of freedom
#Multiple R-squared:  0.5441,    Adjusted R-squared:  0.5425
#F-statistic: 342.5 on 1 and 287 DF,  p-value: < 2.2e-16

#Figure S52:
sprime_rfmix_plot <- function(df, arcorigin, modernanc, shapes_colours) {
	df %>%
		pivot_longer(cols=c(EA, Papuan),
						 names_to="ModernAnc",
						 values_to="AncProportion") %>%
		mutate(Population=factor(Population,
										 levels=shapes_colours$PopLevels,
										 labels=shapes_colours$Population),
				 ModernAnc=case_when(ModernAnc == "EA" ~ "East Asian",
				 						   ModernAnc == "Papuan" ~ "New Guinean",
				 						   TRUE ~ "NA")) %>%
		filter(TractOrigin == arcorigin,
				 ModernAnc == modernanc) %>%
		group_by(Population) %>%
		summarize(OriginMean=mean(ArchaicTractProportion),
					 OriginSD=sd(ArchaicTractTotal)/(2*2684673005),
					 ModernMean=mean(AncProportion),
					 ModernSD=sd(AncProportion)) %>%
		ggplot(aes(x=ModernMean, xmin=ModernMean-ModernSD, xmax=ModernMean+ModernSD,
					  y=OriginMean, ymin=OriginMean-OriginSD, ymax=OriginMean+OriginSD,
					  colour=Population, shape=Population)) +
		   geom_point(size=4) +
		   geom_errorbar() +
		   geom_errorbarh() +
		   scale_x_continuous(limits=c(0, 1),
		   						 expand=c(0, 0.02),
		   						 labels=scales::percent_format(suffix="")) +
		   scale_y_continuous(limits=c(0, NA),
		   						 expand=c(0, 0),
		   						 labels=scales::percent_format(suffix="")) +
		   scale_colour_manual(values=shapes_colours$colours,
		   						  labels=shapes_colours$Population) +
		   scale_shape_manual(values=shapes_colours$shapes,
		   						 labels=shapes_colours$Population) +
		   guides(colour=guide_legend(ncol=5),
		   		 shape=guide_legend(ncol=5)) +
		   theme_classic() +
		   labs(x=paste(modernanc, "ancestry (%)"),
		   	  y=paste(arcorigin, "ancestry (%)"),
		   	  colour="",
		   	  shape="") +
		   theme(legend.position="none")
}
FigS52A <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Neanderthal",
									  "East Asian",
									  Fig2C_shapes_colours)
FigS52B <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Neanderthal",
									  "New Guinean",
									  Fig2C_shapes_colours)
FigS52C <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Denisovan",
									  "East Asian",
									  Fig2C_shapes_colours)
FigS52D <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Denisovan",
									  "New Guinean",
									  Fig2C_shapes_colours)
FigS52E <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Ambiguous",
									  "East Asian",
									  Fig2C_shapes_colours)
FigS52F <- sprime_rfmix_plot(sprime_vs_rfmix_adjusted,
									  "Ambiguous",
									  "New Guinean",
									  Fig2C_shapes_colours)
FigS52_legend <- get_legend(FigS52A + theme(legend.position="bottom"))
FigS52_plots <- plot_grid(FigS52A, FigS52B, FigS52C, FigS52D, FigS52E, FigS52F,
								  align="v",
								  axis="lr",
								  ncol=2,
								  labels="AUTO")
FigS52 <- plot_grid(FigS52_plots, FigS52_legend,
						  ncol=1,
						  rel_heights=c(9, 1))
ggsave(paste0('PIBv1_FigS52_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
		 plot=FigS52,
		 width=24.0,
		 height=32.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS52_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.png'),
		 plot=FigS52,
		 width=24.0,
		 height=32.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS52_SprimeVsRFMix_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
		 plot=FigS52,
		 width=24.0,
		 height=32.0,
		 units="cm",
		 dpi=500)

#Figure 2D was made by Audrey, so see her code.
#It's based on the homozygosity data read in here though.

#Skipped in supplement
#Figure SX:
#homozygosity_curated <- filtered_tract_homozygosity %>%
#   filter(Population %in% curated_pops$name) %>%
#   mutate(Population=factor(Population,
#                            levels=curated_pops$name,
#                            labels=curated_pops$label)) %>%
#   ggplot(aes(x=Population, y=HomozygousArchaic, colour=Region)) +
#      geom_violin(draw_quantiles=c(0.5)) +
#      geom_jitter(alpha=0.75, size=0.25) +
#      theme_bw() +
#      ylim(0, NA) +
#      scale_colour_viridis_d(end=0.9) +
#      labs(x="",
#           y="Homozygous archaic Sprime sites (%)") +
#      guides(colour="none") +
#      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"))
#ggsave(paste0('PIBv1_FigSX_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
#       plot=homozygosity_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_FigSX_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.png'),
#       plot=homozygosity_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_FigSX_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
#       plot=homozygosity_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)

#Formerly figure 2B, now skipped
##Example contour plot with highlighted tract origin criteria rectangles:
#matchrate_example <- filtered_match_rates %>%
#   filter(Population %in% c("Baining-Kagat")) %>%
#   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
#      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
#      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
#      geom_rect(data=archaic_boxes_pfr,
#                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour=TractOrigin),
#                fill=NA,
#                size=2,
#                inherit.aes=FALSE) +
#      scale_x_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_y_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_fill_distiller(palette="Spectral") +
#      scale_colour_manual(values=archaic_origin_colours) +
#      theme_bw() +
#      guides(colour=guide_legend(direction="vertical"),
#             fill=guide_colourbar(title.position="top")) +
#      labs(x="Neanderthal Match Rate",
#           y="Denisovan Match Rate",
#           fill="Normalized\ndensity",
#           colour="Tract Origin") +
#      theme(strip.text=element_text(face="bold"),
#            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#            legend.position="none")
#matchrate_example_legend <- get_legend(matchrate_example + theme(legend.position="top",
#                                                                 legend.box="vertical"))
##Plot one curated example population per world region:
#matchrate_curated <- filtered_match_rates %>%
#   filter(Population %in% curated_pops$name) %>%
#   mutate(Population=factor(Population,
#                            levels=curated_pops$name,
#                            labels=curated_pops$label)) %>%
#   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
#      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
#      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
#      scale_x_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_y_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_fill_distiller(palette="Spectral") +
#      theme_bw() +
#      facet_wrap(~ Population, ncol=3) +
#      guides(fill="none") +
#      labs(x="Neanderthal Match Rate",
#           y="Denisovan Match Rate",
#           fill="Normalized\ndensity") +
#      theme(strip.text=element_text(size=12, face="bold"),
#            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#            axis.title=element_text(size=12))
#Fig2B_bottom_alignment <- align_plots(matchrate_example, matchrate_curated,
#                                      align="h",
#                                      axis="b")
#Fig2B_left <- plot_grid(matchrate_example_legend, Fig2B_bottom_alignment[[1]],
#                        nrow=2,
#                        rel_heights=c(2,1))
#Fig2B <- plot_grid(Fig2B_left, Fig2B_bottom_alignment[[2]],
#                   align="h",
#                   axis="b",
#                   ncol=2,
#                   rel_widths=c(1, 3))
#ggsave(paste0('PIBv1_Fig2B_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
#       plot=Fig2B,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_Fig2B_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.png'),
#       plot=Fig2B,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_Fig2B_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
#       plot=Fig2B,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)

#Actual Figure 2E (sorry about the variable names, it used to be panel D):
#Generate the plot of BMM fits overlaid on histograms of DEN match rate:
#Uses scale_fill_viridis_d() here, and was matched to
# Fig. 2C and 2D in Adobe Illustrator by Dani, though
# it would be possible to use Fig2C_shapes_colours with
# scale_fill_manual() instead.
den_bmmfits_curated <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% curated_pops$name) %>%
   mutate(Population=factor(Population,
                            levels=curated_pops$name,
                            labels=curated_pops$label)) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_bmm_density(Denisovan_bmm_bestfits %>%
                          filter(Population %in% curated_pops$name) %>%
                          mutate(Population=factor(Population,
                                                   levels=curated_pops$name,
                                                   labels=curated_pops$label)),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_viridis_d() +
      facet_wrap(~ Population,
                 ncol=3) +
      guides(fill="none") +
      labs(x="Denisovan match rate",
           y="Number of tracts") +
      theme(strip.text=element_text(face="bold"))
ggsave(paste0('PIBv1_Fig2D_Sprime_BMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=den_bmmfits_curated,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig2D_Sprime_BMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=den_bmmfits_curated,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_Fig2D_Sprime_BMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=den_bmmfits_curated,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#This plot got replaced by plots in PIBv1_Sprime_DENGMMBMMsupplement_20230605.R
#Generate the plot of GMM fits overlaid on histograms of DEN match rate:
#den_gmmfits_curated <- filtered_match_rates %>%
#   filter(Ngood >= 30, Dgood >= 30,
#          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
#   filter(Population %in% curated_pops$name) %>%
#   mutate(Population=factor(Population,
#                            levels=curated_pops$name,
#                            labels=curated_pops$label)) %>%
#   ggplot(aes(x=Dmatchrate, fill=Population)) +
#      geom_histogram(binwidth=0.02) +
#      geom_gmm_density(Denisovan_gmm_bestfits %>%
#                          filter(Population %in% curated_pops$name) %>%
#                          mutate(Population=factor(Population,
#                                                   levels=curated_pops$name,
#                                                   labels=curated_pops$label)),
#                       binwidth=0.02) +
#      theme_bw() +
#      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_fill_viridis_d() +
#      facet_wrap(~ Population,
#                 ncol=3) +
#      guides(fill="none") +
#      labs(x="Denisovan match rate",
#           y="Number of tracts") +
#      theme(strip.text=element_text(face="bold"))
#ggsave(paste0('PIBv1_Fig2D_Sprime_GMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
#       plot=den_gmmfits_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_Fig2D_Sprime_GMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.png'),
#       plot=den_gmmfits_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_Fig2D_Sprime_GMMfits_curated_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
#       plot=den_gmmfits_curated,
#       width=16.0,
#       height=12.0,
#       units="cm",
#       dpi=500)

#Other extra figures:
#Figure S45: Contour plot composite for supplement:
#Example contour plot with highlighted tract origin criteria rectangles:
matchrate_example <- filtered_match_rates %>%
   filter(Population %in% c("Baining-Kagat")) %>%
   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      #Rough criteria:
      geom_rect(data=archaic_boxes_pfr,
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour=TractOrigin),
                fill=NA,
                size=3,
                inherit.aes=FALSE) +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      scale_colour_manual(values=archaic_origin_colours) +
      theme_bw() +
      guides(colour=guide_legend(direction="vertical"),
             fill=guide_colourbar(title.position="top")) +
      labs(x="Neanderthal Match Rate",
           y="Denisovan Match Rate",
           fill="Normalized\ndensity",
           colour="Tract Origin") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            legend.position="top")
#Plot one curated example population per world region:
matchrate_regions <- filtered_match_rates %>%
   filter(Population %in% c("French", "Maya",
                            "Hazara", "Dai",
                            "Rampasasa", "Baining-Kagat")) %>%
   mutate(Population=factor(Population,
                            levels=c("French", "Maya",
                                     "Hazara", "Dai",
                                     "Rampasasa", "Baining-Kagat"),
                            labels=c("EUR-French", "AMR-Maya",
                                     "CSA-Hazara", "EAS-Dai",
                                     "ISEA-Rampasasa", "OCN-Baining-Kagat"))) %>%
   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      theme_bw() +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal Match Rate",
           y="Denisovan Match Rate",
           fill="Normalized\ndensity") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#Now plot only Oceania:
matchrate_OCN <- filtered_match_rates %>%
   filter(Region == "Oceania",
          IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      theme_bw() +
      facet_wrap(~ Population, ncol=3) +
      guides(fill="none") +
      labs(x="Neanderthal Match Rate",
           y="Denisovan Match Rate",
           fill="Normalized\ndensity") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
FigS1_left <- plot_grid(matchrate_example, matchrate_regions,
                        align="v",
                        axis="l",
                        ncol=1,
                        rel_heights=c(1, 1),
                        labels=c("A", "B"))
FigS1 <- plot_grid(FigS1_left, matchrate_OCN,
                   ncol=2,
                   rel_widths=c(2, 3),
                   labels=c("", "C"))
ggsave(paste0('PIBv1_FigS1_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=FigS1,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS1_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=FigS1,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS1_SprimeMatchRates_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=FigS1,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)

#Figure S46: Contour plots for each Neanderthal
#Plot one curated example population per world region for Altai NEA:
matchrate_regions_Altai <- filtered_match_rates %>%
   filter(Population %in% c("French", "Maya",
                            "Hazara", "Dai",
                            "Rampasasa", "Baining-Kagat")) %>%
   mutate(Population=factor(Population,
                            levels=c("French", "Maya",
                                     "Hazara", "Dai",
                                     "Rampasasa", "Baining-Kagat"),
                            labels=c("EUR-French", "AMR-Maya",
                                     "CSA-Hazara", "EAS-Dai",
                                     "ISEA-Rampasasa", "OCN-Baining-Kagat"))) %>%
   ggplot(aes(x=Amatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      theme_bw() +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Altai Match Rate",
           y="Denisovan Match Rate",
           fill="Normalized\ndensity") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#Plot one curated example population per world region for Chagyrskaya NEA:
matchrate_regions_Chagyrskaya <- filtered_match_rates %>%
   filter(Population %in% c("French", "Maya",
                            "Hazara", "Dai",
                            "Rampasasa", "Baining-Kagat")) %>%
   mutate(Population=factor(Population,
                            levels=c("French", "Maya",
                                     "Hazara", "Dai",
                                     "Rampasasa", "Baining-Kagat"),
                            labels=c("EUR-French", "AMR-Maya",
                                     "CSA-Hazara", "EAS-Dai",
                                     "ISEA-Rampasasa", "OCN-Baining-Kagat"))) %>%
   ggplot(aes(x=Cmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      theme_bw() +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Chagyrskaya Match Rate",
           y="",
           fill="Normalized\ndensity") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
#Plot one curated example population per world region for Vindija NEA:
matchrate_regions_Vindija <- filtered_match_rates %>%
   filter(Population %in% c("French", "Maya",
                            "Hazara", "Dai",
                            "Rampasasa", "Baining-Kagat")) %>%
   mutate(Population=factor(Population,
                            levels=c("French", "Maya",
                                     "Hazara", "Dai",
                                     "Rampasasa", "Baining-Kagat"),
                            labels=c("EUR-French", "AMR-Maya",
                                     "CSA-Hazara", "EAS-Dai",
                                     "ISEA-Rampasasa", "OCN-Baining-Kagat"))) %>%
   ggplot(aes(x=Vmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      theme_bw() +
      facet_wrap(~ Population, ncol=2) +
      labs(x="Vindija Match Rate",
           y="",
           fill="Normalized\ndensity") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="none")
matchrate_regions_legend <- get_legend(matchrate_regions_Vindija + theme(legend.position="top"))
FigS46_plots <- plot_grid(matchrate_regions_Altai, matchrate_regions_Chagyrskaya, matchrate_regions_Vindija,
                          align="h",
                          axis="tb",
                          ncol=3,
                          rel_widths=c(1.1, 1, 1),
                          labels="AUTO")
FigS46 <- plot_grid(matchrate_regions_legend, FigS46_plots,
                   nrow=2,
                   rel_heights=c(0.1, 1))
ggsave(paste0('PIBv1_FigS46_SprimeMatchRates_perNEAindiv_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=FigS46,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS46_SprimeMatchRates_perNEAindiv_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=FigS46,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS46_SprimeMatchRates_perNEAindiv_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=FigS46,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)

#Figure S47: Results using Browning et al. 2018 classification criteria:
#Example contour plot with highlighted tract origin criteria rectangles:
matchrate_example_Browning <- filtered_match_rates %>%
   filter(Population %in% c("Baining-Kagat")) %>%
   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
      #Browning criteria:
      geom_rect(data=archaic_boxes_browning,
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour=TractOrigin),
                fill=NA,
                size=2,
                inherit.aes=FALSE) +
      scale_x_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_y_continuous(expand=c(0,0),
                         labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_distiller(palette="Spectral") +
      scale_colour_manual(values=archaic_origin_colours_browning) +
      theme_bw() +
      guides(colour=guide_legend(direction="vertical",
                                 keywidth=0.5,
                                 keyheight=0.5),
             fill="none") +
      labs(x="NEA",
           y="DEN",
           fill="",
           colour="Tract Origin",
           subtitle="Browning 2018 criteria") +
      theme(strip.text=element_text(face="bold"),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            legend.position="top")
#Classified archaic sequence per individual using Browning criteria
tractlen_byorigin_byregion_Browning <- browning_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins)) %>%
   filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
   filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
   ggplot(aes(x=RegionCode, y=ArchaicTractTotal/1e6, colour=RegionCode)) +
      geom_violin() +
      geom_jitter(alpha=0.125, size=0.25) +
      ylim(0, 80) +
      theme_bw() +
      scale_colour_viridis_d(end=0.9) +
      facet_wrap(~ TractOrigin, ncol=1) +
      labs(x="",
           y="Per-Individual Total Archaic Sequence (Mbp)",
           subtitle="") +
      guides(colour="none") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
            strip.text=element_text(face="bold"))
#Plot origin-classified total archaic sequence boxplots for populations in Oceania:
tractlen_byorigin_OCN_Browning <- browning_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins)) %>%
   filter(RegionCode == "OCN",
          IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
   ggplot(aes(x=Population, y=ArchaicTractTotal/1e6, colour=IslandGroup)) +
      geom_violin() +
      geom_jitter(alpha=0.375, size=0.25) +
      ylim(0, 80) +
      theme_bw() +
      scale_colour_brewer(palette="Set2") +
      facet_wrap(~ TractOrigin, ncol=1) +
      labs(x="",
           y="",
           colour="") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
            strip.text=element_text(face="bold"),
            legend.position="top",
            legend.text=element_text(face="bold", size=9),
            legend.key.size=unit(0.25, units="cm"))
FigS3_matchrate_plot <- plot_grid(NULL, matchrate_example_Browning,
                                  ncol=1,
                                  rel_heights=c(0.5, 1),
                                  labels=c("", "A"))
FigS3_tractlen_plots <- align_plots(tractlen_byorigin_byregion_Browning, tractlen_byorigin_OCN_Browning,
                                    align="h",
                                    axis="tb")
FigS3 <- plot_grid(FigS3_matchrate_plot,
                   FigS3_tractlen_plots[[1]], FigS3_tractlen_plots[[2]],
                   align="h",
                   axis="b",
                   ncol=3,
                   rel_widths=c(0.5, 0.5, 1),
                   labels=c("", "B", "C"))
ggsave(paste0('PIBv1_FigS47_SprimeByRegionBrowning_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=FigS3,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS47_SprimeByRegionBrowning_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=FigS3,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS47_SprimeByRegionBrowning_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=FigS3,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)

#Skipped in supplement
#Figure S3b:
#Example contour plot with highlighted tract origin criteria rectangles:
#matchrate_example <- filtered_match_rates %>%
#   filter(Population %in% c("Baining-Kagat")) %>%
#   ggplot(aes(x=Nmatchrate, y=Dmatchrate)) +
#      stat_density_2d(geom="raster", aes(fill=after_stat(ndensity)), contour=FALSE) +
#      geom_density_2d(binwidth=0.1, contour_var="ndensity") +
#      #Rough criteria:
#      geom_rect(data=archaic_boxes_pfr,
#                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour=TractOrigin),
#                fill=NA,
#                size=2,
#                inherit.aes=FALSE) +
#      scale_x_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_y_continuous(expand=c(0,0),
#                         labels=function(x) {paste0(format(x*100), "%")}) +
#      scale_fill_distiller(palette="Spectral") +
#      scale_colour_manual(values=archaic_origin_colours) +
#      theme_bw() +
#      guides(colour=guide_legend(direction="vertical",
#                                 keywidth=0.5,
#                                 keyheight=0.5),
#             fill="none") +
#      labs(x="NEA",
#           y="DEN",
#           fill="",
#           colour="Tract Origin") +
#      theme(strip.text=element_text(face="bold"),
#            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#            legend.position="top")
##Classified archaic sequence per individual using our criteria
#tractlen_byorigin_byregion <- rough_classified_perindiv_lengths %>%
#   mutate(TractOrigin=factor(TractOrigin,
#                             levels=archaic_origins)) %>%
#   filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
#   filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
#   ggplot(aes(x=RegionCode, y=ArchaicTractTotal/1e6, colour=RegionCode)) +
#      geom_violin() +
#      geom_jitter(alpha=0.125, size=0.25) +
#      ylim(0, 150) +
#      theme_bw() +
#      scale_colour_viridis_d(end=0.9) +
#      facet_wrap(~ TractOrigin, ncol=1) +
#      labs(x="",
#           y="Per-Individual Total Archaic Sequence (Mbp)") +
#      guides(colour=FALSE) +
#      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
#            strip.text=element_text(face="bold"))
##Plot origin-classified total archaic sequence boxplots for populations in Oceania:
#tractlen_byorigin_OCN <- rough_classified_perindiv_lengths %>%
#   mutate(TractOrigin=factor(TractOrigin,
#                             levels=archaic_origins)) %>%
#   filter(RegionCode == "OCN",
#          IslandGroup != "Vanuatu") %>%
#   filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
#   ggplot(aes(x=Population, y=ArchaicTractTotal/1e6, colour=IslandGroup)) +
#      geom_violin() +
#      geom_jitter(alpha=0.375, size=0.25) +
#      ylim(0, 150) +
#      theme_bw() +
#      scale_colour_brewer(palette="Set2") +
#      facet_wrap(~ TractOrigin, ncol=1) +
#      labs(x="",
#           y="",
#           colour="") +
#      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
#            strip.text=element_text(face="bold"),
#            legend.position="top",
#            legend.text=element_text(face="bold", size=9),
#            legend.key.size=unit(0.25, units="cm"))
#FigS3b_matchrate_plot <- plot_grid(NULL, matchrate_example,
#                                   ncol=1,
#                                   rel_heights=c(0.5, 1),
#                                   labels=c("", "A"))
#FigS3b_tractlen_plots <- align_plots(tractlen_byorigin_byregion, tractlen_byorigin_OCN,
#                                     align="h",
#                                     axis="tb")
#FigS3b <- plot_grid(FigS3b_matchrate_plot,
#                    FigS3b_tractlen_plots[[1]], FigS3b_tractlen_plots[[2]],
#                    align="h",
#                    axis="b",
#                    ncol=3,
#                    rel_widths=c(0.5, 0.5, 1),
#                    labels=c("", "B", "C"))
#ggsave(paste0('PIBv1_FigS3b_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
#       plot=FigS3b,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_FigS3b_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.png'),
#       plot=FigS3b,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)
#ggsave(paste0('PIBv1_FigS3b_SprimeByRegion_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
#       plot=FigS3b,
#       width=24.0,
#       height=18.0,
#       units="cm",
#       dpi=500)

#Figure S48: Comparison of archaic sequence per individual with subsampled Sprime runs
#Load the subsampled match rates and projection tract lengths:
match_rates_subsample <- read_tsv('subset_10samples/PIBv1_perPop_Sprime_autosomal_match_rates.tsv.gz',
                                  col_types='ciiciininininini')
arc_tract_lengths_maxgap0_subsample <- read_tsv('subset_10samples/PIBv1_Sprime_subsample_perChrom_perHaplotype_perPop_tract_lengths_maxgap0_sumtractsfix.tsv.gz',
                                                  col_types='cccciii')
#Filter the projections:
filtered_tracts_lengths_subsample <- arc_tract_lengths_maxgap0_subsample %>%
   left_join(match_rates_subsample, by=c("CHROM", "TractID")) %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(RegionCode=factor(Region,
                            levels=c("Africa", "Middle East",
                                     "Europe", "America",
                                     "Central/South Asia", "East Asia",
                                     "Island Southeast Asia", "Oceania"),
                            labels=c("AFR", "MDE",
                                     "EUR", "AMR",
                                     "CSA", "EAS",
                                     "ISEA", "OCN")),
   		 Population=str_replace(Population, ",", "-")) %>%
	mutate(Population=case_when(Population %in% c("Santa-Cruz", "Vella-Lavella") ~ str_replace(Population, "-", " "),
										 TRUE ~ Population)) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., RegionCode, Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 30, Dgood >= 30, Nmatchrate >= 0.3 | Dmatchrate >= 0.3)
#Classify projections:
rough_classified_perindiv_lengths_subsample <- filtered_tracts_lengths_subsample %>%
   mutate(TractOrigin=case_when(Ngood >= 30 & Dgood >= 30 & Nmatchrate >= 0.3 & Dmatchrate <= 0.3 ~ "Neanderthal",
                                Ngood >= 30 & Dgood >= 30 & Nmatchrate <= 0.3 & Dmatchrate >= 0.3 ~ "Denisovan",
                                Ngood >= 30 & Dgood >= 30 & Nmatchrate > 0.3 & Dmatchrate > 0.3 ~ "Ambiguous",
                                TRUE ~ NA_character_)) %>%
   group_by(RegionCode, Region, IslandGroup, Population, Haplotype, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(TRACTLEN)) %>%
   separate(Haplotype, into=c("Sample"), remove=FALSE, extra="drop") %>%
   group_by(RegionCode, Region, IslandGroup, Population, Sample, TractOrigin) %>%
   summarize(ArchaicTractTotal=sum(ArchaicTractTotal))
#Get median per-individual archaic sequence values per population from all-sample vs. subsample:
median_perindiv_classified_seq <- rough_classified_perindiv_lengths %>%
   filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
   filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
   group_by(Region, RegionCode, IslandGroup, Population, TractOrigin) %>%
   summarize(MedianLen=median(ArchaicTractTotal),
             MeanLen=mean(ArchaicTractTotal),
             Analysis="full sample")
median_perindiv_classified_seq_subsample <- rough_classified_perindiv_lengths_subsample %>%
   filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
   filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
   group_by(Region, RegionCode, IslandGroup, Population, TractOrigin) %>%
   summarize(MedianLen=median(ArchaicTractTotal),
             MeanLen=mean(ArchaicTractTotal),
             Analysis="n=10 subsample")
#Combine and plot:
FigS48 <- bind_rows(median_perindiv_classified_seq,
                   median_perindiv_classified_seq_subsample) %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=archaic_origins)) %>%
   ggplot(aes(x=Population, y=MedianLen/1e6, fill=TractOrigin)) +
      geom_col(position="stack", width=0.5) +
      theme_bw() +
      scale_y_continuous(expand=c(0,0), limits=c(0, 200)) +
      scale_fill_manual(values=archaic_origin_colours) +
      guides(fill="none") +
      labs(x="",
           y="Median Per-Individual Archaic Sequence (Mbp)",
           fill="") +
      facet_grid(Analysis ~ RegionCode, scales="free_x", space="free_x") +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
            legend.position="top",
            legend.direction="horizontal")
ggsave(paste0('PIBv1_FigS48_SprimeSubsample_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
       plot=FigS48,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS48_SprimeSubsample_', format(Sys.Date(), format="%Y%m%d"), '.png'),
       plot=FigS48,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
ggsave(paste0('PIBv1_FigS48_SprimeSubsample_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
       plot=FigS48,
       width=24.0,
       height=18.0,
       units="cm",
       dpi=500)
#To examine the per-population differences subsampling makes:
bind_rows(rough_classified_perindiv_lengths %>%
             filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
             filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
             mutate(Analysis="full"),
          rough_classified_perindiv_lengths_subsample %>%
             filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
             filter(!(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
             mutate(Analysis="subsample")) %>%
   pivot_wider(id_cols=c(RegionCode, Region, IslandGroup, Population, Sample, TractOrigin),
               names_from=Analysis,
               values_from=ArchaicTractTotal) %>%
   mutate(diff=full-subsample) %>%
   group_by(Region, RegionCode, IslandGroup, Population, Sample) %>%
   summarize(diff=sum(diff)) %>%
   group_by(Region, RegionCode, IslandGroup, Population) %>%
   summarize(mindiff=min(diff),
             mediandiff=median(diff),
             maxdiff=max(diff)) %>%
   arrange(mediandiff) %>%
   View()

#Figures S53 and S54 (Sprime-ROH overlap):
library(ggstar)
sprime_roh_overlap <- read_tsv('ROH/PIBv1_Sprime_projection_ROH_overlap_maxgap0.tsv')
sprime_vs_roh <- rough_classified_perindiv_lengths %>%
	filter((is.na(IslandGroup) | IslandGroup != "Vanuatu"),
			 !(Population %in% c("Goroka-Sepik", "Nakanai", "Bellona"))) %>%
	left_join(sprime_roh_overlap %>%
				 	 mutate(Sample=SampleID,
				 			  TractOrigin=case_when(ArchaicOrigin == "Neandertal" ~ "Neanderthal",
				 			 							   TRUE ~ ArchaicOrigin),
				 			  ROHOverlap=case_when(!is.na(OverlapSum) ~ OverlapSum,
				 			 							  TRUE ~ 0),
				 			  ProjOverlap=case_when(!is.na(ASum) ~ ASum,
				 			 							   TRUE ~ 0),
				 			  .keep="none"),
				 by=c("Sample", "TractOrigin")) %>%
	pivot_wider(names_from=TractOrigin,
					values_from=c(ArchaicTractTotal, ROHOverlap, ProjOverlap)) %>%
	mutate(ArchaicTractTotal_Archaic=ArchaicTractTotal_Ambiguous+ArchaicTractTotal_Denisovan+ArchaicTractTotal_Neanderthal,
			 ROHOverlap_Archaic=ROHOverlap_Ambiguous+ROHOverlap_Denisovan+ROHOverlap_Neanderthal,
			 ProjOverlap_Archaic=ProjOverlap_Ambiguous+ProjOverlap_Denisovan+ProjOverlap_Neanderthal) %>%
	pivot_longer(cols=-c(RegionCode, Region, IslandGroup, Population, Sample),
					 names_to=c(".value", "TractOrigin"),
					 names_sep="_") %>%
	mutate(TractOrigin=factor(TractOrigin,
									  levels=c("Archaic", "Neanderthal", "Denisovan", "Ambiguous"),
									  labels=c("All Archaic", "Neanderthal", "Denisovan", "Ambiguous")))

chang_ww_colours <- c(rep("#8AC626", 4), rep("#C7A57A", 8),
							 rep("#FC2C07", 5), rep("#CF53CD", 9),
							 rep("#A0101A", 3), "#773410", rep("#A0101A", 4), "#773410", rep("#A0101A", 2), rep("#773410", 2), "#A0101A", rep("#773410", 2), "#A0101A", "#773410", rep("#A0101A", 2),
							 rep("#5D6D7D", 3),
							 rep("#05006C", 2), "#1C79FE", rep("#3350D9", 2), rep("#1C79FE", 6), rep("#76C3E6", 2), "#16B1FE", rep("#76C3E6", 2), "#16B1FE", "#76C3E6")
chang_ww_shapes <- c(23, 15, 29, 13,
							22, 11, 15, 5, 13, 29, 23, 12,
							15, 29, 13, 5, 23,
							23, 15, 22, 13, 11, 29, 12, 5, 6,
							15, 23, 13, 13, 18, 5, 28, 29, 29, 6, 17, 12, 23, 1, 5, 11, 11, 15, 12, 22,
							15, 23, 13,
							15, 23, 15, 23, 13, 29, 12, 22, 1, 5, 6,
							23, 29, 13, 15, 13, 15, 5)
chang_ocn_colours <- c(rep("#8AC626", 2),
							  rep("#C7A57A", 7),
							  "#5D6D7D", "#1F7C1A",
							  rep("#C55418", 2),
							  rep("#CF53CD", 5))
chang_ocn_shapes <- c(23, 15,
							 5, 23, 15, 22, 29, 12, 13,
							 15, 15,
							 15, 23,
							 13, 5, 29, 15, 12)

sprime_vs_roh_plot <- function(arcorigin, df, colours, shapes) {
	df %>%
		group_by(TractOrigin, Population) %>%
		summarize(ArchaicTractTotal=mean(ArchaicTractTotal)/1e6,
					 ROHOverlap=mean(ROHOverlap)/1e6,
					 ProjOverlap=mean(ProjOverlap)/1e6) %>%
		filter(TractOrigin == arcorigin) %>%
		ggplot(aes(x=ROHOverlap, y=ArchaicTractTotal, colour=Population, fill=Population, starshape=Population)) +
		   geom_star() +
	      theme_classic() +
	      labs(x=paste("Overlap of", arcorigin, "sequence with ROH (Mbp)"),
	      	  y=paste(arcorigin, "per individual (Mbp)"),
	      	  colour="",
	      	  fill="",
	      	  starshape="") +
		   scale_x_continuous(limits=c(0, NA)) +
	      scale_y_continuous(limits=c(0, NA)) +
		   scale_colour_manual(values=colours) +
		   scale_fill_manual(values=colours) +
		   scale_starshape_manual(values=shapes) +
		   theme(axis.title=element_text(size=8),
		   		legend.position="none",
		   		legend.key.size=unit(0.3, "cm"))
}

FigS5354_arc_origins <- c("All Archaic", "Neanderthal", "Denisovan", "Ambiguous")

FigS53_plots <- sprime_vs_roh %>%
	lapply(FigS5354_arc_origins, sprime_vs_roh_plot, ., chang_ww_colours, chang_ww_shapes)
FigS53_legend <- get_legend(FigS53_plots[[1]] + theme(legend.position="bottom", legend.key.size=unit(0.3, "cm")))
FigS53_nolegend <- plot_grid(plotlist=FigS53_plots,
									  align="hv",
									  axis="trbl",
									  ncol=2,
									  labels="AUTO")
FigS53 <- plot_grid(FigS53_nolegend, FigS53_legend,
						  ncol=1,
						  rel_heights=c(3, 1))
ggsave(paste0('PIBv1_FigS53_SprimeROH_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
		 plot=FigS53,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS53_SprimeROH_', format(Sys.Date(), format="%Y%m%d"), '.png'),
		 plot=FigS53,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS53_SprimeROH_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
		 plot=FigS53,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)

FigS54_plots <- sprime_vs_roh %>%
	filter(RegionCode == "OCN") %>%
	lapply(FigS5354_arc_origins, sprime_vs_roh_plot, ., chang_ocn_colours, chang_ocn_shapes)
FigS54_legend <- get_legend(FigS54_plots[[1]] + theme(legend.position="bottom", legend.key.size=unit(0.3, "cm")))
FigS54_nolegend <- plot_grid(plotlist=FigS54_plots,
									  align="hv",
									  axis="trbl",
									  ncol=2,
									  labels="AUTO")
FigS54 <- plot_grid(FigS54_nolegend, FigS54_legend,
						  ncol=1,
						  rel_heights=c(6, 1))
ggsave(paste0('PIBv1_FigS54_SprimeROHOCN_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
		 plot=FigS54,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS54_SprimeROHOCN_', format(Sys.Date(), format="%Y%m%d"), '.png'),
		 plot=FigS54,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS54_SprimeROHOCN_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
		 plot=FigS54,
		 width=24.0,
		 height=24.0,
		 units="cm",
		 dpi=500)

#Skipped in supplement
#Figure S4: Archaic homozygosity (regional and OCN population distributions):
#homozygosity_byregion <- filtered_tract_homozygosity %>%
#	filter((is.na(IslandGroup) | IslandGroup != "Vanuatu")) %>%
#	filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
#	ggplot(aes(x=RegionCode, y=HomozygousArchaic, colour=RegionCode)) +
#	   geom_violin(draw_quantiles=c(0.5)) +
#	   geom_jitter(alpha=0.75, size=0.25) +
#	   theme_bw() +
#	   ylim(0, NA) +
#	   scale_colour_viridis_d(end=0.9) +
#	   labs(x="",
#	   	  y="",
#	   	  subtitle="Homozygous archaic Sprime sites (%)") +
#	   guides(colour="none") +
#	   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"))
#homozygosity_OCN <- filtered_tract_homozygosity %>%
#	filter(RegionCode == "OCN",
#			 IslandGroup != "Vanuatu") %>%
#	filter(!(Population %in% c("Goroka,Sepik", "Nakanai", "Bellona"))) %>%
#	ggplot(aes(x=Population, y=HomozygousArchaic, colour=IslandGroup)) +
#	   geom_violin(draw_quantiles=c(0.5)) +
#	   geom_jitter(alpha=0.75, size=0.25) +
#	   theme_bw() +
#	   ylim(0, NA) +
#	   scale_colour_brewer(palette="Set2") +
#	   labs(x="",
#	   	  y="",
#	   	  colour="") +
#	   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, face="bold"),
#	   		legend.position="top",
#	   		legend.direction="horizontal",
#	   		legend.text=element_text(size=8))
#FigS4 <- plot_grid(homozygosity_byregion, homozygosity_OCN,
#						 align='v',
#						 axis='lr',
#						 nrow=2,
#						 rel_heights=c(1, 2),
#						 labels="AUTO")
#ggsave(paste0('PIBv1_FigS4_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.pdf'),
#		 plot=FigS4,
#		 width=16.0,
#		 height=12.0,
#		 units="cm",
#		 dpi=500)
#ggsave(paste0('PIBv1_FigS4_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.png'),
#		 plot=FigS4,
#		 width=16.0,
#		 height=12.0,
#		 units="cm",
#		 dpi=500)
#ggsave(paste0('PIBv1_FigS4_SprimeHomozygosity_', format(Sys.Date(), format="%Y%m%d"), '.ps'),
#		 plot=FigS4,
#		 width=16.0,
#		 height=12.0,
#		 units="cm",
#		 dpi=500)

#Skipped in supplement, but useful
#Sprime tract score distributions passing vs. failing the match rate filters:
annotated_match_rates <- match_rates %>%
   separate(col=TractID,
            into=c("Population"),
            sep="_",
            remove=FALSE,
            extra="drop") %>%
   left_join(pop_region_map,
             by=c("Population"="AnalysisGroup")) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, IslandGroup)$Population)),
          MRpass=case_when(Ngood >= 30 & Dgood >= 30 & (Nmatchrate >= 0.3 | Dmatchrate >= 0.3) ~ TRUE,
                           TRUE ~ FALSE))

OCN_score_dist_filter <- annotated_match_rates %>%
   filter(Region == "Oceania",
          is.na(IslandGroup) | IslandGroup != "Vanuatu",
          !(Population %in% c("Goroka", "Sepik", "Nakanai", "Bellona"))) %>%
   ggplot(aes(x=SCORE, colour=MRpass)) +
      geom_step(stat="ecdf") +
      theme_bw() +
      scale_x_continuous(name="Sprime tract score",
                         limits=c(1e5, 1e6),
                         labels=function(x) {format(x, scientific=TRUE)}) +
      scale_colour_brewer(palette="Set2") +
      facet_wrap(~ Population, ncol=2) +
      labs(y="eCDF proportion",
           colour="Passing match rate filters") +
      theme(legend.position="top",
            legend.direction="horizontal",
            axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave(paste0("PIBv1_Sprime_OCN_score_distributions_matchratefilter_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=OCN_score_dist_filter,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_Sprime_OCN_score_distributions_matchratefilter_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=OCN_score_dist_filter,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_Sprime_OCN_score_distributions_matchratefilter_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=OCN_score_dist_filter,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Supplementary table S_TotalSprime aka S5:
#by population:
total_perindiv_lengths %>%
   mutate(Region=factor(Region,
                        levels=c("Middle East", "Europe", "America",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"))) %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Goroka-Sepik", "Bellona", "Nakanai"))) %>%
   group_by(Region, Population) %>%
   summarize(minlen=min(ArchaicTractTotal)/1e6,
             Q1len=quantile(ArchaicTractTotal, prob=0.25)/1e6,
             medianlen=median(ArchaicTractTotal)/1e6,
             Q3len=quantile(ArchaicTractTotal, prob=0.75)/1e6,
             maxlen=max(ArchaicTractTotal)/1e6) %>%
   View()
#by region, removing redundant individuals by retaining lumped populations:
total_perindiv_lengths %>%
   mutate(Region=factor(Region,
                        levels=c("Middle East", "Europe", "America",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"))) %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Bellona", "Nakanai", "Goroka-Sepik"))) %>%
   group_by(Region) %>%
   summarize(minlen=min(ArchaicTractTotal)/1e6,
             Q1len=quantile(ArchaicTractTotal, prob=0.25)/1e6,
             medianlen=median(ArchaicTractTotal)/1e6,
             Q3len=quantile(ArchaicTractTotal, prob=0.75)/1e6,
             maxlen=max(ArchaicTractTotal)/1e6) %>%
   View()

#Supplementary table S_ClassifiedArchaic aka S6:
#by population:
rough_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=c("Neanderthal", "Denisovan", "Ambiguous")),
          Region=factor(Region,
                        levels=c("Middle East", "Europe", "America",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"))) %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Goroka-Sepik", "Bellona", "Nakanai"))) %>%
   group_by(Region, Population, TractOrigin) %>%
   summarize(medianlen=median(ArchaicTractTotal)/1e6,
             minlen=min(ArchaicTractTotal)/1e6,
             maxlen=max(ArchaicTractTotal)/1e6) %>%
   View()
#by region, removing redundant individuals by retaining lumped populations:
rough_classified_perindiv_lengths %>%
   mutate(TractOrigin=factor(TractOrigin,
                             levels=c("Neanderthal", "Denisovan", "Ambiguous")),
          Region=factor(Region,
                        levels=c("Middle East", "Europe", "America",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"))) %>%
   filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
   filter(!(Population %in% c("Bellona", "Nakanai", "Goroka-Sepik"))) %>%
   group_by(Region, TractOrigin) %>%
   summarize(medianlen=median(ArchaicTractTotal)/1e6,
             minlen=min(ArchaicTractTotal)/1e6,
             maxlen=max(ArchaicTractTotal)/1e6) %>%
   View()

#Archaic coverage analyses:
#Load the archaic coverage data:
arc_cov_maxgap0 <- read_tsv('PIBv1_Sprime_path_length_fromProjections_maxgap0.tsv')

#Correlations among archaic origins with populations as data points:
arc_cov_maxgap0 %>%
	filter(!(Region %in% c("All", "OCN", "ISEA", "EAS",
								  "CSA", "AMR", "EUR", "MDE",
								  "nonOCN"))) %>%
	select(TractOrigin, Region, TotalPathLength) %>%
	pivot_wider(id_cols=Region,
					names_from=TractOrigin,
					values_from=TotalPathLength) %>%
	mutate(RegionGroup=case_when(Region %in% c("Ata", "Baining-Kagat",
															 "Baining-Mali", "Bellona,Rennell",
															 "Goroka,Sepik", "Kove",
															 "Lavongai-Mussau", "Malaita",
															 "Mamusi", "Melamela",
															 "Nailik-Notsi-Tigak", "Nakanai,Mangseng",
															 "Nasioi", "Santa-Cruz",
															 "Saposa", "Tikopia",
															 "Vella-Lavella", "Agta",
															 "Cebuano", "Rampasasa") ~ "ISEAOCN",
										  TRUE ~ "nonISEAOCN")) %>%
	group_by(RegionGroup) %>%
	summarize(AMBDEN_r=cor(Ambiguous, Denisovan, method="pearson"),
				 AMBDEN_rho=cor(Ambiguous, Denisovan, method="spearman"),
				 AMBDEN_n=n(),
				 AMBNEA_r=cor(Ambiguous, Neandertal, method="pearson"),
				 AMBNEA_rho=cor(Ambiguous, Neandertal, method="spearman"),
				 AMBNEA_n=n()) %>%
	print(n=Inf, width=Inf)
#RegionGroup AMBDEN_r AMBDEN_rho AMBDEN_n AMBNEA_r AMBNEA_rho AMBNEA_n
#ISEAOCN        0.831      0.863       20    0.545      0.592       20
#nonISEAOCN     0.457      0.386       46    0.914      0.859       46

#Figure S55: Observed vs. Expected archaic coverage:
#Based on a naïve estimator of expected archaic coverage derived from
# the naïve sequencing coverage estimator E[C]=1-e^{-N*L/G}
#We estimate E[C] as E[C]=1-e^{-(1/G)*sum_i{L_i}}.
arc_total_perpop_lengths <- rough_classified_perindiv_lengths %>%
	filter(is.na(IslandGroup) | IslandGroup != "Vanuatu") %>%
	filter(!(Population %in% c("Goroka", "Sepik", "Nakanai", "Bellona"))) %>%
	group_by(RegionCode, Population, TractOrigin) %>%
	summarize(TotalArchaicLength=sum(ArchaicTractTotal))
expected_arc_cov <- bind_rows(arc_total_perpop_lengths %>%
											ungroup() %>%
											select(Population, TractOrigin, TotalArchaicLength) %>%
											rename(Region=Population),
										arc_total_perpop_lengths %>%
											group_by(RegionCode, TractOrigin) %>%
											summarize(TotalArchaicLength=sum(TotalArchaicLength)) %>%
											ungroup() %>%
											mutate(Region=as.character(RegionCode)) %>%
											select(-RegionCode)) %>%
	mutate(EArcCov=1-exp(-TotalArchaicLength/2684673005))
OE_arc_cov <- arc_cov_maxgap0 %>%
	filter(!(Region %in% c("All", "nonOCN"))) %>%
	mutate(Region=str_replace(Region, ",", "-")) %>%
	mutate(Region=case_when(Region %in% c("Santa-Cruz", "Vella-Lavella") ~ str_replace(Region, "-", " "),
									TRUE ~ Region)) %>%
	transmute(Region=Region,
				 TractOrigin=case_when(TractOrigin == "Neandertal" ~ "Neanderthal",
				 							 TRUE ~ TractOrigin),
				 TotalPathLength=TotalPathLength,
				 OArcCov=TotalPathLength/2684673005) %>%
	inner_join(expected_arc_cov,
				  by=c("Region"="Region", "TractOrigin"="TractOrigin")) %>%
	mutate(ArcCovRatio=OArcCov/EArcCov) %>%
	filter(!(Region %in% c("OCN", "ISEA", "EAS", "CSA", "MDE", "EUR", "AMR"))) %>%
	mutate(Bottlenecked=case_when(Region %in% c("Ata", "Baining-Kagat", "Baining-Mali", "Bellona-Rennell", "Mamusi") ~ TRUE,
											TRUE ~ FALSE)) %>%
	rename(Population=Region) %>%
	left_join(pop_region_map %>%
				 	 mutate(AnalysisGroup=str_replace(AnalysisGroup, ",", "-")) %>%
				 	 mutate(AnalysisGroup=case_when(AnalysisGroup %in% c("Santa-Cruz", "Vella-Lavella") ~ str_replace(AnalysisGroup, "-", " "),
				 	 										  TRUE ~ AnalysisGroup)),
				 by=c("Population"="AnalysisGroup"))
arc_cov_bottlenecked_plot <- OE_arc_cov %>%
	ggplot(aes(x=EArcCov, y=OArcCov, colour=Region, shape=Bottlenecked)) +
	   geom_smooth(aes(x=EArcCov, y=OArcCov, linetype=Bottlenecked), inherit.aes=FALSE,
	   				method="lm", se=FALSE, col="grey50") +
	   geom_point() +
	   geom_abline(slope=1, intercept=0, linetype=2, col="grey80") +
	   theme_bw() +
	   scale_x_continuous(limits=c(0, NA),
	   						 expand=c(0.01, 0.01)) +
	   scale_y_continuous(limits=c(0, NA),
	   						 expand=c(0.01, 0.01)) +
	   scale_colour_brewer(palette="Set2") +
	   scale_linetype_manual(values=c(4, 5)) +
	   labs(x="Expected Archaic Coverage",
	   	  y="Observed Archaic Coverage") +
	   facet_wrap(~ TractOrigin, ncol=3)
ggsave(paste0('PIBv1_FigS55new_Sprime_archaic_coverage_bottleneck_effect_', format(Sys.Date(), format="%Y%m%d"), ".pdf"),
		 plot=arc_cov_bottlenecked_plot,
		 width=24.0,
		 height=18.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS55new_Sprime_archaic_coverage_bottleneck_effect_', format(Sys.Date(), format="%Y%m%d"), ".png"),
		 plot=arc_cov_bottlenecked_plot,
		 width=24.0,
		 height=18.0,
		 units="cm",
		 dpi=500)
ggsave(paste0('PIBv1_FigS55new_Sprime_archaic_coverage_bottleneck_effect_', format(Sys.Date(), format="%Y%m%d"), ".ps"),
		 plot=arc_cov_bottlenecked_plot,
		 width=24.0,
		 height=18.0,
		 units="cm",
		 dpi=500)
#Regression for Denisovan coverage:
OE_arc_cov %>%
	filter(TractOrigin == "Denisovan") %>%
	lm(data=., OArcCov ~ EArcCov + Bottlenecked) %>%
	summary()
#Call:
#lm(formula = OArcCov ~ EArcCov + Bottlenecked, data = .)
#Residuals:
#       Min         1Q     Median         3Q        Max
#-0.0280676 -0.0024331 -0.0008751  0.0004878  0.0213591
#Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.0030929  0.0009837   3.144  0.00254 **
#EArcCov           0.2788429  0.0107736  25.882  < 2e-16 ***
#BottleneckedTRUE -0.0249718  0.0036450  -6.851 3.63e-09 ***
#Residual standard error: 0.006584 on 63 degrees of freedom
#Multiple R-squared:  0.9218,    Adjusted R-squared:  0.9193
#F-statistic: 371.5 on 2 and 63 DF,  p-value: < 2.2e-16
#
#So residuals are symmetric and approximately 0-centered, R^2 of fit is
# decent, intercept is minimal as expected, unit effect size of expected
# archaic coverage is about 28%, and effect of being a bottlenecked group
# is negative and significant (~-2.5%).
OE_arc_cov %>%
	filter(TractOrigin == "Neanderthal") %>%
	lm(data=., OArcCov ~ EArcCov + Bottlenecked) %>%
	summary()
#Call:
#lm(formula = OArcCov ~ EArcCov + Bottlenecked, data = .)
#Residuals:
#      Min        1Q    Median        3Q       Max
#-0.032424 -0.008140  0.000503  0.008241  0.032119
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.046553   0.006726   6.922 2.74e-09 ***
#EArcCov           0.120196   0.028805   4.173 9.39e-05 ***
#BottleneckedTRUE -0.034356   0.007311  -4.699 1.46e-05 ***
#Residual standard error: 0.01514 on 63 degrees of freedom
#Multiple R-squared:  0.3313,    Adjusted R-squared:  0.3101
#F-statistic: 15.61 on 2 and 63 DF,  p-value: 3.119e-06
#
#So residuals are symmetric and approximately 0-centered, R^2 of fit is
# not great, intercept is minimal as expected, unit effect size of expected
# archaic coverage is about 12%, and effect of being a bottlenecked group
# is negative and significant (~-3.4%).
OE_arc_cov %>%
	filter(TractOrigin == "Ambiguous") %>%
	lm(data=., OArcCov ~ EArcCov + Bottlenecked) %>%
	summary()
#Call:
#lm(formula = OArcCov ~ EArcCov + Bottlenecked, data = .)
#Residuals:
#      Min        1Q    Median        3Q       Max
#-0.022287 -0.005441  0.001544  0.006398  0.015577
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.015191   0.002196   6.918 2.77e-09 ***
#EArcCov           0.192720   0.013510  14.265  < 2e-16 ***
#BottleneckedTRUE -0.022934   0.004589  -4.998 4.90e-06 ***
#Residual standard error: 0.008624 on 63 degrees of freedom
#Multiple R-squared:  0.7678,    Adjusted R-squared:  0.7605
#F-statistic: 104.2 on 2 and 63 DF,  p-value: < 2.2e-16
#
#So residuals are symmetric and approximately 0-centered, R^2 of fit is
# decent, intercept is minimal as expected, unit effect size of expected
# archaic coverage is about 19%, and effect of being a bottlenecked group
# is negative and significant (~-2.3%).


#Cleanup:
rm(geom_bmm_density, geom_gmm_density)
rm(PIBv1_metadata, pop_region_map, curated_pops)
rm(archaic_origins, archaic_origin_colours, archaic_origin_colours_browning, archaic_boxes_pfr, archaic_boxes_browning)
rm(region_colours)
rm(arc_tract_lengths_maxgap0, filtered_tracts_lengths)
rm(total_perindiv_lengths, rough_classified_perindiv_lengths, browning_classified_perindiv_lengths)
rm(tractlen_byorigin_byregion_stackedbar, tractlen_byorigin_OCN_stackedbar, Fig2A_legend, Fig2A_top, Fig2A_vertalign, Fig2A)
rm(arclen_label_coords, ocnzoom)
rm(match_rates, filtered_match_rates)
rm(matchrate_example, matchrate_example_legend, matchrate_curated, Fig2B_bottom_alignment, Fig2B_left, Fig2B)
rm(sprime_vs_rfmix, sprime_vs_rfmix_adjusted)
rm(Fig2C_shapes_colours, Fig2C, Fig2C_legend)
rm(FigS52A, FigS52B, FigS52C, FigS52D, FigS52E, FigS52F, FigS52_legend)
rm(FigS52_plots, FigS52, sprime_rfmix_plot)
rm(arc_tract_lengths, filtered_tract_homozygosity)
rm(homozygosity_curated)
rm(matchrate_regions, matchrate_OCN, FigS1_left, FigS1)
rm(matchrate_regions_Altai, matchrate_regions_Chagyrskaya, matchrate_regions_Vindija, matchrate_regions_legend, FigS46_plots, FigS46)
rm(matchrate_example_Browning, tractlen_byorigin_byregion_Browning, tractlen_byorigin_OCN_Browning, FigS3_matchrate_plot, FigS3_tractlen_plots, FigS3)
rm(matchrate_example, tractlen_byorigin_byregion, tractlen_byorigin_OCN, FigS3b_matchrate_plot, FigS3b_tractlen_plots, FigS3b)
rm(homozygosity_byregion, homozygosity_OCN, FigS4)
rm(match_rates_subsample, arc_tract_lengths_maxgap0_subsample, filtered_tracts_lengths_subsample, rough_classified_perindiv_lengths_subsample)
rm(median_perindiv_classified_seq, median_perindiv_classified_seq_subsample, FigS48)
rm(annotated_match_rates, OCN_score_dist_filter)
rm(chang_ww_colours, chang_ww_shapes, chang_ocn_colours, chang_ocn_shapes)
rm(sprime_roh_overlap, sprime_vs_roh, sprime_vs_roh_plot, FigS5354_arc_origins)
rm(FigS53_plots, FigS53_nolegend, FigS53_legend, FigS53)
rm(FigS54_plots, FigS54_nolegend, FigS54_legend, FigS54)
rm(total_arc_seq_plot)
rm(Denisovan_gmms, Denisovan_gmm_bestfits, Denisovan_bmms, Denisovan_bmm_bestfits)
rm(den_bmmfits_curated, den_gmmfits_curated)
rm(arc_total_perpop_lengths)
rm(arc_cov_maxgap0, expected_arc_cov, OE_arc_cov, arc_cov_bottlenecked_plot)
