#PIBv1 Sprime NEA GMM and BMM fits and plots

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
#load('PIBv1_Denisovan_gmms_bmms_k1thru10_noinputs_Rv2.Rdata')

#Annotate tracts and apply a baseline filter on ascertainabile sites:
filtered_match_rates <- match_rates %>%
   separate(col=TractID, into=c("Population"), sep="_", remove=FALSE, extra="drop") %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, IslandGroup)$Population))) %>%
   filter(Ngood >= 10 | Dgood >= 10)

#Exclude Vanuatu populations:
pops <- setdiff((filtered_match_rates %>%
                    summarize(Pops=unique(Population)))$Pops,
                c("Ambae", "Ambae,Maewo", "Ambrym",
                  "Efate", "Emae","Espiritu-Santo",
                  "Maewo", "Malakula", "Pentecost",
                  "Tanna", "Ureparapara"))

#Set the PRNG seed for the mixture model fits:
set.seed(42)
#Set the number of parallel threads for the model fits:
registerDoParallel(4)

#Run the GMM fits in parallel:
Neanderthal_gmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   cat(paste0("Running GMM for ", p, "\n"));
   cat(paste0(filtered_match_rates %>%
                 filter(Ngood >= 30, Dgood >= 30,
                        Nmatchrate > 0.3, Dmatchrate < 0.3,
                        Population == p) %>%
                 nrow(),
              " tracts retained for ", p, "\n"));
   gmm <- filtered_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate > 0.3, Dmatchrate < 0.3,
             Population == p) %>%
      mutate(Nmatchrate=case_when(Nmatchrate == 1.0 ~ 0.999999,
                                  TRUE ~ Nmatchrate)) %>%
      stepFlexmix(Nmatchrate ~ 1,
                  data=.,
                  model=FLXMRglm(family="gaussian"),
                  k=seq(1, 10, by=1));
   l <- list(gmm);
   names(l) <- p;
   l
}
#Run the BMM fits in parallel:
Neanderthal_bmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   cat(paste0("Running BMM for ", p, "\n"));
   cat(paste0(filtered_match_rates %>%
                 filter(Ngood >= 30, Dgood >= 30,
                        Nmatchrate > 0.3, Dmatchrate < 0.3,
                        Population == p) %>%
                 nrow(),
              " tracts retained for ", p, "\n"));
   bmm <- filtered_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate > 0.3, Dmatchrate < 0.3,
             Population == p) %>%
      mutate(Nmatchrate=case_when(Nmatchrate == 1.0 ~ 0.999999,
                                  TRUE ~ Nmatchrate)) %>%
      possibly(.f=betamix, otherwise=NULL)(Nmatchrate ~ 1,
                                           data=.,
                                           link="logit",
                                           k=seq(1, 10, by=1),
                                           which="BIC",
                                           verbose=TRUE);
   l <- list(bmm);
   names(l) <- p;
   l
}

#Summarize the best fits by their parameter values:
#These parameter values will be used for the density overlays.
Neanderthal_gmm_bestfits <- data.frame(Population=c(),
                                       component=c(),
                                       alpha=c(),
                                       mu=c(),
                                       stdev=c())
Neanderthal_bmm_bestfits <- data.frame(Population=c(),
                                       component=c(),
                                       alpha=c(),
                                       shape1=c(),
                                       shape2=c())
for (p in names(Neanderthal_gmms)) {
   best_fit <- getModel(Neanderthal_gmms[[p]], "BIC");
   Neanderthal_gmm_bestfits <- bind_rows(Neanderthal_gmm_bestfits,
                                         data.frame(Population=p,
                                                    component=seq(1, best_fit@k),
                                                    alpha=best_fit@prior,
                                                    mu=unname(t(parameters(best_fit)))[,1],
                                                    stdev=unname(t(parameters(best_fit)))[,2]));
}
for (p in names(Neanderthal_bmms)) {
   if (!is.null(Neanderthal_bmms[[p]])) {
      k <- Neanderthal_bmms[[p]]$flexmix@k;
      if (k < 2) {
         mu <- plogis(coef(Neanderthal_bmms[[p]])[1]);
         nu <- exp(coef(Neanderthal_bmms[[p]])[2]);
      } else {
         mu <- plogis(coef(Neanderthal_bmms[[p]])[,1]);
         nu <- exp(coef(Neanderthal_bmms[[p]])[,2]);
      };
      alpha <- Neanderthal_bmms[[p]]$flexmix@prior;
      Neanderthal_bmm_bestfits <- bind_rows(Neanderthal_bmm_bestfits,
                                            data.frame(Population=p,
                                                       component=seq(1, k),
                                                       alpha=alpha,
                                                       shape1=unname(mu*nu),
                                                       shape2=unname((1-mu)*nu)));
   };
}
#Add the total tract counts per population to the fits so that
# the densities can be adjusted into counts:
Neanderthal_gmm_bestfits <- Neanderthal_gmm_bestfits %>%
   left_join(filtered_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Nmatchrate=0.0)
Neanderthal_bmm_bestfits <- Neanderthal_bmm_bestfits %>%
   left_join(filtered_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Nmatchrate=0.0)

#Save the inputs and results in case we need them in the future:
save(filtered_match_rates, Neanderthal_gmms, Neanderthal_bmms,
     Neanderthal_gmm_bestfits, Neanderthal_bmm_bestfits,
     geom_gmm_density, geom_bmm_density,
     file='PIBv1_Neanderthal_gmms_bmms_k1thru10_Rv2.Rdata',
     compress=TRUE,
     compression_level=9,
     version=2)

#If you want to load from the .Rdata file, start here:
load('PIBv1_Neanderthal_gmms_bmms_k1thru10_Rv2.Rdata')

#Supplementary figure S56 in the PIBv1 manuscript:
#Generate the subplot for European population GMM fits:
EUR_NEA_GMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Adygei", "Basque", "French", "Italian-Bergamo",
                            "Orcadian", "Russian", "Sardinian", "Tuscan")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_gmm_density(Neanderthal_gmm_bestfits %>%
                          filter(Population %in% c("Adygei", "Basque", "French", "Italian-Bergamo",
                                                   "Orcadian", "Russian", "Sardinian", "Tuscan")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="Number of tracts",
           subtitle="GMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#Now the BMM subplot:
EUR_NEA_BMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Adygei", "Basque", "French", "Italian-Bergamo",
                            "Orcadian", "Russian", "Sardinian", "Tuscan")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_bmm_density(Neanderthal_bmm_bestfits %>%
                          filter(Population %in% c("Adygei", "Basque", "French", "Italian-Bergamo",
                                                   "Orcadian", "Russian", "Sardinian", "Tuscan")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="",
           subtitle="BMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
EUR_NEA_plot <- plot_grid(EUR_NEA_GMM_plot, EUR_NEA_BMM_plot,
                          align="h",
                          axis="tb",
                          ncol=2,
                          rel_widths=c(1,1),
                          labels="AUTO")
ggsave('PIBv1_FigS_NEAMMEUR_20230604.pdf',
       plot=EUR_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMEUR_20230604.png',
       plot=EUR_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMEUR_20230604.ps',
       plot=EUR_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in the supplement
#Generate the subplot for Central/South Asian population GMM fits:
CSA_NEA_GMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Balochi", "Brahui", "Burusho",
                            "Hazara", "Kalash", "Makrani",
                            "Pathan", "Sindhi", "Uygur")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_gmm_density(Neanderthal_gmm_bestfits %>%
                          filter(Population %in% c("Balochi", "Brahui", "Burusho",
                                                   "Hazara", "Kalash", "Makrani",
                                                   "Pathan", "Sindhi", "Uygur")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="Number of tracts",
           subtitle="GMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#Now the BMM subplot:
CSA_NEA_BMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Balochi", "Brahui", "Burusho",
                            "Hazara", "Kalash", "Makrani",
                            "Pathan", "Sindhi", "Uygur")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_bmm_density(Neanderthal_bmm_bestfits %>%
                          filter(Population %in% c("Balochi", "Brahui", "Burusho",
                                                   "Hazara", "Kalash", "Makrani",
                                                   "Pathan", "Sindhi", "Uygur")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="",
           subtitle="BMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
CSA_NEA_plot <- plot_grid(CSA_NEA_GMM_plot, CSA_NEA_BMM_plot,
                          align="h",
                          axis="tb",
                          ncol=2,
                          rel_widths=c(1,1),
                          labels="AUTO")
ggsave('PIBv1_FigS_NEAMMCSA_20230604.pdf',
       plot=CSA_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMCSA_20230604.png',
       plot=CSA_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMCSA_20230604.ps',
       plot=CSA_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in the supplement
#Generate the subplot for Northern East Asian population GMM fits:
NEAS_NEA_GMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Yakut", "Oroqen", "Mongolian", "Daur",
                            "Hezhen", "Xibo", "Japanese", "Tu",
                            "NorthernHan", "Han")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_gmm_density(Neanderthal_gmm_bestfits %>%
                          filter(Population %in% c("Yakut", "Oroqen", "Mongolian", "Daur",
                                                   "Hezhen", "Xibo", "Japanese", "Tu",
                                                   "NorthernHan", "Han")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="Number of tracts",
           subtitle="GMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#Now the BMM subplot:
NEAS_NEA_BMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Yakut", "Oroqen", "Mongolian", "Daur",
                            "Hezhen", "Xibo", "Japanese", "Tu",
                            "NorthernHan", "Han")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_bmm_density(Neanderthal_bmm_bestfits %>%
                          filter(Population %in% c("Yakut", "Oroqen", "Mongolian", "Daur",
                                                   "Hezhen", "Xibo", "Japanese", "Tu",
                                                   "NorthernHan", "Han")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="",
           subtitle="BMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
NEAS_NEA_plot <- plot_grid(NEAS_NEA_GMM_plot, NEAS_NEA_BMM_plot,
                           align="h",
                           axis="tb",
                           ncol=2,
                           rel_widths=c(1,1),
                           labels="AUTO")
ggsave('PIBv1_FigS_NEAMMNEAS_20230604.pdf',
       plot=NEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMNEAS_20230604.png',
       plot=NEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMNEAS_20230604.ps',
       plot=NEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Skipped in the supplement
#Generate the subplot for Southern East Asian population GMM fits:
SEAS_NEA_GMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Tujia", "Yi", "Miao", "She",
                            "Naxi", "Atayal", "Paiwan", "Lahu",
                            "Dai", "Cambodian")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_gmm_density(Neanderthal_gmm_bestfits %>%
                          filter(Population %in% c("Tujia", "Yi", "Miao", "She",
                                                   "Naxi", "Atayal", "Paiwan", "Lahu",
                                                   "Dai", "Cambodian")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="Number of tracts",
           subtitle="GMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
#Now the BMM subplot:
SEAS_NEA_BMM_plot <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate > 0.3, Dmatchrate < 0.3) %>%
   filter(Population %in% c("Tujia", "Yi", "Miao", "She",
                            "Naxi", "Atayal", "Paiwan", "Lahu",
                            "Dai", "Cambodian")) %>%
   ggplot(aes(x=Nmatchrate, fill=Population)) +
      geom_histogram(binwidth=0.02) +
      geom_bmm_density(Neanderthal_bmm_bestfits %>%
                          filter(Population %in% c("Tujia", "Yi", "Miao", "She",
                                                   "Naxi", "Atayal", "Paiwan", "Lahu",
                                                   "Dai", "Cambodian")),
                       binwidth=0.02) +
      theme_bw() +
      scale_x_continuous(labels=function(x) {paste0(format(x*100), "%")}) +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population, ncol=2) +
      guides(fill="none") +
      labs(x="Neanderthal match rate",
           y="",
           subtitle="BMM") +
      theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
SEAS_NEA_plot <- plot_grid(SEAS_NEA_GMM_plot, SEAS_NEA_BMM_plot,
                           align="h",
                           axis="tb",
                           ncol=2,
                           rel_widths=c(1,1),
                           labels="AUTO")
ggsave('PIBv1_FigS_NEAMMSEAS_20230604.pdf',
       plot=SEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMSEAS_20230604.png',
       plot=SEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave('PIBv1_FigS_NEAMMSEAS_20230604.ps',
       plot=SEAS_NEA_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Clean up:
rm(PIBv1_metadata, pop_region_map, pops)
rm(match_rates, filtered_match_rates)
rm(Neanderthal_gmms, Neanderthal_bmms, Neanderthal_gmm_bestfits, Neanderthal_bmm_bestfits)
rm(EUR_NEA_GMM_plot, EUR_NEA_BMM_plot, EUR_NEA_plot)
rm(CSA_NEA_GMM_plot, CSA_NEA_BMM_plot, CSA_NEA_plot)
rm(NEAS_NEA_GMM_plot, NEAS_NEA_BMM_plot, NEAS_NEA_plot)
rm(SEAS_NEA_GMM_plot, SEAS_NEA_BMM_plot, SEAS_NEA_plot)
rm(geom_gmm_density, geom_bmm_density)
