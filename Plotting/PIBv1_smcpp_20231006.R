#PIBv1 SMC++ results:

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
library(cowplot)

#Get to the directory:
pushd('[path redacted]/PIBv1_results/Demography/')

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
excluded_samples <- read_tsv('[path redacted]/PIBv1_samples_to_exclude.txt',
                             col_names="SampleID",
                             skip=0)
PIBv1_metadata <- read_tsv('[path redacted]/PIBv1_metadata_v0.3.tsv',
                           col_types='ccccclcccnncccclc',
                           na=c("NA", "")) %>%
   mutate(Region=factor(Region,
                        levels=region_names,
                        labels=region_codes)) %>%
   anti_join(excluded_samples,
             by="SampleID")

pop_region_map <- PIBv1_metadata %>%
   group_by(AnalysisGroup) %>%
   summarize(Region=first(Region),
             IslandGroup=first(Island),
             SampleSize=n())

#Scalars for normalizing estimates:
generation_time <- 29
mu <- 1.25e-8
#xlim for these plots:
time_limits <- c(1e3, 2e6)
#ylim for these plots:
Ne_limits <- c(3e2, 7e5)
plot_params <- list("generation_time"=generation_time,
                    "Ne_limits"=Ne_limits,
                    "time_limits"=time_limits)

#Load the SMC++ fits:
smcpp <- read_csv('smcpp/PIBv1_58pops_smcpp_fits.csv.gz') %>%
   rename(Population=label,
          Time=x,
          Ne=y,
          Iteration=plot_num) %>%
   mutate(Iteration=case_when(Iteration > 0 ~ as.character(Iteration),
                              TRUE ~ "final")) %>%
   left_join(pop_region_map, by=c("Population"="AnalysisGroup")) %>%
   mutate(PopWithSize=str_c(Population, " (n=", SampleSize, ")")) %>%
   rename(PopWithoutSize=Population) %>%
   rename(Population=PopWithSize)

#Plotting function:
smcpp_plot <- function(df, params, title, palette, theme="bw", section="supplement") {
   generation_time <- params[["generation_time"]]
   Ne_limits <- params[["Ne_limits"]]
   time_limits <- params[["time_limits"]]
   theme_func <- match.fun(paste0("theme_", theme))
   df %>%
      ggplot(aes(x=Time*generation_time, y=Ne, colour=Population)) +
         geom_step(direction="hv") +
         theme_func() +
         labs(x=paste0("Years (G=", generation_time, ")"),
              y=expression(N[e](t)),
              colour=title) +
         scale_x_log10(limits=time_limits,
                       expand=c(0,0)) +
         scale_y_log10(limits=Ne_limits,
                       expand=c(0,0)) +
         annotation_logticks(sides="bl") +
         {
            if (length(palette) == 1 && str_length(palette) == 1) {
               scale_colour_viridis_d(option=palette)
            } else if (length(palette) == 1) {
               if (substr(palette, 1, 1) == "-") {
                  if (str_length(palette) == 2) {
                     scale_colour_viridis_d(option=substr(palette, 2, str_length(palette)),
                                            direction=-1)
                  } else {
                     scale_colour_brewer(palette=substr(palette, 2, str_length(palette)),
                                         direction=-1)
                  }
               } else {
                  scale_colour_brewer(palette=palette)
               }
            } else {
               scale_colour_manual(values=palette)
            }
         } +
         guides(colour=guide_legend(byrow=TRUE)) +
         {
            if (section == "main") {
               guides(colour=guide_legend(nrow=4,
                                          keywidth=0.7,
                                          keyheight=0.7))
            } else {
               guides(colour=guide_legend(byrow=TRUE))
            }
         } +
         theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
               legend.text=element_text(size=10),
               legend.spacing.y=unit(0, "cm"),
               legend.key.height=unit(0.4, "cm")) +
         {
            if (section == "main") {
               if (theme == "classic") {
                  theme(panel.border=element_rect(fill=NA, colour="black"),
                        legend.position="top",
                        legend.title=element_blank())
               } else {
                  theme(legend.position="top",
                        legend.title=element_blank())
               }
            } else {
               theme(legend.position="right",
                     legend.title=element_text(face="bold"))
            }
         }
}

#Figure 1 panel F was based on this code with some colour tweaks:
#Main text figure of curated populations:
curated_pops <- c("Yoruba", "French", "Han",
                  "Bellona", "Tikopia", "Goroka",
                  "Nakanai", "Ata", "Baining-Kagat",
                  "Baining-Mali")
curated_smcpp_plot <- smcpp %>%
   filter(PopWithoutSize %in% curated_pops,
          Iteration == "final") %>%
   mutate(PopWithoutSize=factor(PopWithoutSize,
                                levels=curated_pops)) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., PopWithoutSize)$Population))) %>%
   smcpp_plot(params=plot_params,
              title="",
              palette="-Paired",
              theme="classic",
              section="main")
ggsave(paste0("PIBv1_Fig1F_smcpp_curated_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=curated_smcpp_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_Fig1F_smcpp_curated_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=curated_smcpp_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_Fig1F_smcpp_curated_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=curated_smcpp_plot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Generate each subplot by major geographic region:
#We split OCN in two, so skip those for now:
smcpp_byregion <- list()
for (r in setdiff(region_codes, "OCN")) {
   smcpp_byregion[[r]] <- smcpp %>%
      filter(Region == r,
             Iteration == "final") %>%
      smcpp_plot(params=plot_params,
                 title=r,
                 palette="C")
}
#OCN is split into populations in the Bismarck archipelago vs. not:
smcpp_byregion[["Bismarck"]] <- smcpp %>%
   filter(Region == "OCN",
          IslandGroup %in% c("New Britain", "New Hanover", "New Ireland"),
          Iteration == "final") %>%
   smcpp_plot(params=plot_params,
              title="Bismarck",
              palette="C")
smcpp_byregion[["non-Bismarck"]] <- smcpp %>%
   filter(Region == "OCN",
          IslandGroup %in% c("New Guinea", "Bougainville", "Solomon Islands"),
          Iteration == "final") %>%
   smcpp_plot(params=plot_params,
              title="non-Bismarck",
              palette="C")

#Combined plots:
#Supplementary figure S32 in the PIBv1 manuscript:
AFR_MDE_CSA_AMR <- plot_grid(smcpp_byregion[["AFR"]],
                             smcpp_byregion[["MDE"]],
                             smcpp_byregion[["CSA"]],
                             smcpp_byregion[["AMR"]],
                             ncol=1,
                             align="v",
                             axis="lr",
                             labels="AUTO")
ggsave(paste0("PIBv1_smcpp_AFR_MDE_CSA_AMR_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=AFR_MDE_CSA_AMR,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_AFR_MDE_CSA_AMR_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=AFR_MDE_CSA_AMR,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_AFR_MDE_CSA_AMR_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=AFR_MDE_CSA_AMR,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Supplementary figure S33 in the PIBv1 manuscript:
AFR_EUR_EAS_ISEA <- plot_grid(smcpp_byregion[["AFR"]],
                              smcpp_byregion[["EUR"]],
                              smcpp_byregion[["EAS"]],
                              smcpp_byregion[["ISEA"]],
                              ncol=1,
                              align="v",
                              axis="lr",
                              labels="AUTO")
ggsave(paste0("PIBv1_smcpp_AFR_EUR_EAS_ISEA_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=AFR_EUR_EAS_ISEA,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_AFR_EUR_EAS_ISEA_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=AFR_EUR_EAS_ISEA,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_AFR_EUR_EAS_ISEA_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=AFR_EUR_EAS_ISEA,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#Supplementary figure S34 in the PIBv1 manuscript:
EAS_ISEA_OCN <- plot_grid(smcpp_byregion[["EAS"]],
                          smcpp_byregion[["ISEA"]],
                          smcpp_byregion[["Bismarck"]],
                          smcpp_byregion[["non-Bismarck"]],
                          ncol=1,
                          align="v",
                          axis="lr",
                          labels="AUTO")
ggsave(paste0("PIBv1_smcpp_EAS_ISEA_OCN_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=EAS_ISEA_OCN,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_EAS_ISEA_OCN_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=EAS_ISEA_OCN,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)
ggsave(paste0("PIBv1_smcpp_EAS_ISEA_OCN_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=EAS_ISEA_OCN,
       width=16.0,
       height=24.0,
       units="cm",
       dpi=500)

#

#Clean up:
rm(r, region_codes, region_names, region_pretty)
rm(generation_time, Ne_limits, time_limits, plot_params, mu, curated_pops)
rm(PIBv1_metadata, excluded_samples, pop_region_map, smcpp, smcpp_byregion)
rm(AFR_MDE_CSA_AMR, AFR_EUR_EAS_ISEA, EAS_ISEA_OCN, curated_smcpp_plot)
rm(smcpp_plot)
