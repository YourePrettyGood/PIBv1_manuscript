#PIBv1 Sprime DEN GMM and BMM supplementary plots

#Load the libraries:
library(tidyverse)
library(cowplot)
library(doParallel)
library(foreach)

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
                                      "New Guinea", "New Britain", "Lavongai",
                                      "New Ireland", "Bougainville", "Solomon Islands",
                                      "Vanuatu")))

#Load the match rates and projected tract lengths:
match_rates <- read_tsv('PIBv1_perPop_Sprime_autosomal_match_rates.tsv.gz',
                        col_types='ciiciininininini')

#Load the Denisovan BMM and GMM results:
load('PIBv1_Denisovan_gmms_bmms_k1thru10_noinputs_Rv2.Rdata')

#Annotate tracts and apply a baseline filter on ascertainabile sites:
filtered_match_rates <- match_rates %>%
   separate(col=TractID, into=c("Population"), sep="_", remove=FALSE, extra="drop") %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 10 | Dgood >= 10)

gmm_bmm_plot <- function(matchrates, gmm_fits, bmm_fits, pops, origin=c("DEN", "NEA"), prefix) {
   binwidth <- 0.02
   originrate <- case_when(origin == "DEN" ~ "Dmatchrate",
                           origin == "NEA" ~ "Nmatchrate")
   originlong <- case_when(origin == "DEN" ~ "Denisovan",
                           origin == "NEA" ~ "Neanderthal")
   if (origin == "DEN") {
      filtered_matchrates <- matchrates %>%
         filter(Ngood >= 30, Dgood >= 30,
                Dmatchrate > 0.3, Nmatchrate < 0.3) %>%
         filter(Population %in% pops)
   } else if (origin == "NEA") {
      filtered_matchrates <- matchrates %>%
         filter(Ngood >= 30, Dgood >= 30,
                Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
         filter(Population %in% pops)
   }
   filtered_gmm_fits <- gmm_fits %>%
      filter(Population %in% pops)
   filtered_bmm_fits <- bmm_fits %>%
      filter(Population %in% pops)
   gmm_plot <- filtered_matchrates %>%
      ggplot(aes(x=.data[[originrate]], fill=Population)) +
         geom_histogram(binwidth=binwidth) +
         geom_gmm_density(filtered_gmm_fits,
                          binwidth=binwidth) +
         theme_bw() +
         scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
         scale_fill_brewer(palette="Spectral") +
         facet_wrap(~ Population, ncol=2) +
         guides(fill="none") +
         labs(x=paste0(originlong, " match rate"),
              y="Number of tracts",
              subtitle="GMM") +
         theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
   bmm_plot <- filtered_matchrates %>%
      ggplot(aes(x=.data[[originrate]], fill=Population)) +
         geom_histogram(binwidth=binwidth) +
         geom_bmm_density(filtered_bmm_fits,
                          binwidth=binwidth) +
         theme_bw() +
         scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
         scale_fill_brewer(palette="Spectral") +
         facet_wrap(~ Population, ncol=2) +
         guides(fill="none") +
         labs(x=paste0(originlong, " match rate"),
              y="",
              subtitle="BMM") +
         theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
   full_plot <- plot_grid(gmm_plot, bmm_plot,
                          align="h",
                          axis="tb",
                          ncol=2,
                          rel_widths=c(1,1),
                          labels="AUTO")
   ggsave(paste0(prefix, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=full_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(prefix, "_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=full_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(prefix, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=full_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
}

#Supplementary figure S57 in PIBv1 manuscript:
#Generate the Central/South Asian population plot:
gmm_bmm_plot(matchrates=filtered_match_rates,
             gmm_fits=Denisovan_gmm_bestfits,
             bmm_fits=Denisovan_bmm_bestfits,
             pops=c("Balochi", "Brahui", "Burusho",
                    "Hazara", "Kalash", "Makrani",
                    "Pathan", "Sindhi", "Uygur"),
             origin="DEN",
             prefix="PIBv1_FigS_DENMMCSA")

#Supplementary figure S58 in PIBv1 manuscript:
#Generate the Northern East Asian population plot:
gmm_bmm_plot(matchrates=filtered_match_rates,
             gmm_fits=Denisovan_gmm_bestfits,
             bmm_fits=Denisovan_bmm_bestfits,
             pops=c("Yakut", "Oroqen", "Mongolian", "Daur", "Hezhen",
                    "Xibo", "Japanese", "Tu", "NorthernHan", "Han"),
             origin="DEN",
             prefix="PIBv1_FigS_DENMMNEAS")

#Supplementary figure S59 in PIBv1 manuscript:
#Generate the Southern East Asian population plot:
gmm_bmm_plot(matchrates=filtered_match_rates,
             gmm_fits=Denisovan_gmm_bestfits,
             bmm_fits=Denisovan_bmm_bestfits,
             pops=c("Tujia", "Yi", "Miao", "She", "Naxi",
                    "Atayal", "Paiwan", "Lahu", "Dai", "Cambodian"),
             origin="DEN",
             prefix="PIBv1_FigS_DENMMSEAS")

#Supplementary figure S61 in PIBv1 manuscript:
#Generate the Oceanic PIB-specific population plot:
gmm_bmm_plot(matchrates=filtered_match_rates,
             gmm_fits=Denisovan_gmm_bestfits,
             bmm_fits=Denisovan_bmm_bestfits,
             pops=c("Ata", "Baining-Kagat", "Baining-Mali", "Kove", "Lavongai-Mussau",
                    "Mamusi", "Melamela", "Nailik-Notsi-Tigak", "Nakanai,Mangseng", "Saposa"),
             origin="DEN",
             prefix="PIBv1_FigS_DENMMOCNPIB")

#Supplementary figure S62 in PIBv1 manuscript:
#Generate the Oceanic non-PIB population plot:
gmm_bmm_plot(matchrates=filtered_match_rates,
             gmm_fits=Denisovan_gmm_bestfits,
             bmm_fits=Denisovan_bmm_bestfits,
             pops=c("Goroka", "Sepik", "Nasioi", "Bellona,Rennell",
                    "Tikopia", "Malaita", "Santa-Cruz", "Vella-Lavella"),
             origin="DEN",
             prefix="PIBv1_FigS_DENMMOCNnonPIB")

#Clean up:
rm(PIBv1_metadata, pop_region_map)
rm(match_rates, filtered_match_rates)
rm(Denisovan_gmms, Denisovan_bmms, Denisovan_gmm_bestfits, Denisovan_bmm_bestfits)
rm(geom_gmm_density, geom_bmm_density, gmm_bmm_plot)
