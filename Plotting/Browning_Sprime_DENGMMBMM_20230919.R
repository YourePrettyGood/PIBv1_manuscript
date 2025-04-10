#Browning et al. 2018 1kGP GMM and BMM fits for PIBv1 manuscript:

#Load the libraries:
library(tidyverse)
library(cowplot)
library(doParallel)
library(foreach)
library(flexmix)
library(betareg)

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

#Load the data:
#Load the match file and calculate match rates, plus some data cleaning:
browning_match_rates <- read_tsv('[path redacted]/Browning_etal_2018_autosomes_Sprime_byPop_match_rates.tsv.gz') %>%
   unite("Segment_ID", CHROM, SEGMENT) %>%
   group_by(Group, Segment_ID) %>%
   summarize(snplen=n_distinct(POS),
             posstart=min(POS),
             posend=max(POS),
             ONmatch=sum(NMATCH == "match"),
             ONgood=sum(NMATCH %in% c("match", "mismatch")),
             Nmatch=sum(NeandertalMatch == "match"),
             Ngood=sum(NeandertalMatch %in% c("match", "mismatch")),
             ODmatch=sum(DMATCH == "match"),
             ODgood=sum(DMATCH %in% c("match", "mismatch")),
             Dmatch=sum(DenisovanMatch == "match"),
             Dgood=sum(DenisovanMatch %in% c("match", "mismatch")),
             Amatch=sum(AltaiMatch == "match"),
             Agood=sum(AltaiMatch %in% c("match", "mismatch")),
             Vmatch=sum(VindijaMatch == "match"),
             Vgood=sum(VindijaMatch %in% c("match", "mismatch")),
             Cmatch=sum(ChagyrskayaMatch == "match"),
             Cgood=sum(ChagyrskayaMatch %in% c("match", "mismatch"))) %>%
   mutate(poslen=posend-posstart+1,
          ONmatchrate=ONmatch/ONgood,
          Nmatchrate=Nmatch/Ngood,
          ODmatchrate=ODmatch/ODgood,
          Dmatchrate=Dmatch/Dgood,
          Amatchrate=Amatch/Agood,
          Vmatchrate=Vmatch/Vgood,
          Cmatchrate=Cmatch/Cgood) %>%
   mutate(Group=str_replace(Group, fixed("Sprime/Browning_etal_2018_results/"), "")) %>%
   mutate(Group=factor(Group,
                       levels=c("CEU", "FIN", "GBR", "IBS", "TSI",
                                "CDX", "CHB", "CHS", "JPT", "KHV",
                                "BEB", "GIH", "ITU", "PJL", "STU",
                                "MXL", "PUR", "CLM", "PEL", "Papuans")))

#Now set up for the parallelized mixture models:
set.seed(42)
registerDoParallel(4)
pops <- unique(browning_match_rates$Group)

#Run the parallelized GMMs with PFR match rates:
Browning_DEN_gmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   gmm <- browning_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate < 0.3, Dmatchrate > 0.3,
             Group == p) %>%
      mutate(Dmatchrate=case_when(Dmatchrate == 1.0 ~ 0.999999,
                                  TRUE ~ Dmatchrate)) %>%
      stepFlexmix(Dmatchrate ~ 1,
                  data=.,
                  model=FLXMRglm(family="gaussian"),
                  k=seq(1, 10, by=1));
   l <- list(gmm);
   names(l) <- p;
   l
}

#Now the parallelized BMMs with PFR match rates:
Browning_DEN_bmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   bmm <- browning_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate < 0.3, Dmatchrate > 0.3,
             Group == p) %>%
      mutate(Dmatchrate=case_when(Dmatchrate == 1.0 ~ 0.999999,
                                  TRUE ~ Dmatchrate)) %>%
      possibly(.f=betamix, otherwise=NULL)(Dmatchrate ~ 1,
                                           data=.,
                                           link="logit",
                                           k=seq(1, 10, by=1),
                                           which="BIC",
                                           verbose=TRUE);
   l <- list(bmm);
   names(l) <- p;
   l
}

#Get the best fit models for each population:
Browning_DEN_gmm_bestfits <- data.frame(Population=c(),
                                        component=c(),
                                        alpha=c(),
                                        mu=c(),
                                        stdev=c())
Browning_DEN_bmm_bestfits <- data.frame(Population=c(),
                                        component=c(),
                                        alpha=c(),
                                        shape1=c(),
                                        shape2=c())
for (p in names(Browning_DEN_gmms)) {
   best_fit <- getModel(Browning_DEN_gmms[[p]], "BIC");
   Browning_DEN_gmm_bestfits <- bind_rows(Browning_DEN_gmm_bestfits,
                                          data.frame(Population=p,
                                                     component=seq(1, best_fit@k),
                                                     alpha=best_fit@prior,
                                                     mu=unname(t(parameters(best_fit)))[,1],
                                                     stdev=unname(t(parameters(best_fit)))[,2]));
}
for (p in names(Browning_DEN_bmms)) {
   if (!is.null(Browning_DEN_bmms[[p]])) {
      k <- Browning_DEN_bmms[[p]]$flexmix@k;
      if (k < 2) {
         mu <- plogis(coef(Browning_DEN_bmms[[p]])[1]);
         nu <- exp(coef(Browning_DEN_bmms[[p]])[2]);
      } else {
         mu <- plogis(coef(Browning_DEN_bmms[[p]])[,1]);
         nu <- exp(coef(Browning_DEN_bmms[[p]])[,2]);
      };
      alpha <- Browning_DEN_bmms[[p]]$flexmix@prior;
      Browning_DEN_bmm_bestfits <- bind_rows(Browning_DEN_bmm_bestfits,
                                             data.frame(Population=p,
                                                        component=seq(1, k),
                                                        alpha=alpha,
                                                        shape1=unname(mu*nu),
                                                        shape2=unname((1-mu)*nu)));
   };
}

#Add the number of tracts used in each model fit, and a placeholder
# match rate column to make ggplot happy:
Browning_DEN_gmm_bestfits <- Browning_DEN_gmm_bestfits %>%
   left_join(browning_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
                rename(Population=Group) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)
Browning_DEN_bmm_bestfits <- Browning_DEN_bmm_bestfits %>%
   left_join(browning_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
                rename(Population=Group) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)

#Now plot the best fits against histograms of tract match rates:
#This plot is used as Supplementary figure S60 in the PIBv1 manuscript.
gmm_bmm_plot(matchrates=browning_match_rates %>%
                rename(Population=Group),
             gmm_fits=Browning_DEN_gmm_bestfits,
             bmm_fits=Browning_DEN_bmm_bestfits,
             pops=c("CDX", "CHB", "CHS", "JPT", "KHV",
                    "BEB", "GIH", "ITU", "PJL", "STU"),
             origin="DEN",
             prefix="Browning2018_DENMM")
#Save the runs in case we want them later:
save(Browning_DEN_gmms, Browning_DEN_bmms,
     Browning_DEN_gmm_bestfits, Browning_DEN_bmm_bestfits,
     geom_gmm_density, geom_bmm_density, gmm_bmm_plot,
     file='Browning2018_DENMM_PFRmatchrates_20230919.Rdata',
     compression_level=9,
     compress=TRUE)

#Repeat the process for the Browning match rates:
#Run the parallelized GMMs with Browning match rates:
Browning_ODEN_gmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   gmm <- browning_match_rates %>%
      filter(ONgood >= 30, ODgood >= 30,
             ONmatchrate < 0.3, ODmatchrate > 0.3,
             Group == p) %>%
      mutate(ODmatchrate=case_when(ODmatchrate == 1.0 ~ 0.999999,
                                   TRUE ~ ODmatchrate)) %>%
      stepFlexmix(ODmatchrate ~ 1,
                  data=.,
                  model=FLXMRglm(family="gaussian"),
                  k=seq(1, 10, by=1));
   l <- list(gmm);
   names(l) <- p;
   l
}

#Now the parallelized BMMs with Browning match rates:
Browning_ODEN_bmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   bmm <- browning_match_rates %>%
      filter(ONgood >= 30, ODgood >= 30,
             ONmatchrate < 0.3, ODmatchrate > 0.3,
             Group == p) %>%
      mutate(ODmatchrate=case_when(ODmatchrate == 1.0 ~ 0.999999,
                                   TRUE ~ ODmatchrate)) %>%
      possibly(.f=betamix, otherwise=NULL)(ODmatchrate ~ 1,
                                           data=.,
                                           link="logit",
                                           k=seq(1, 10, by=1),
                                           which="BIC",
                                           verbose=TRUE);
   l <- list(bmm);
   names(l) <- p;
   l
}

#Get the best fit models for each population:
Browning_ODEN_gmm_bestfits <- data.frame(Population=c(),
                                         component=c(),
                                         alpha=c(),
                                         mu=c(),
                                         stdev=c())
Browning_ODEN_bmm_bestfits <- data.frame(Population=c(),
                                         component=c(),
                                         alpha=c(),
                                         shape1=c(),
                                         shape2=c())
for (p in names(Browning_ODEN_gmms)) {
   best_fit <- getModel(Browning_ODEN_gmms[[p]], "BIC");
   Browning_ODEN_gmm_bestfits <- bind_rows(Browning_ODEN_gmm_bestfits,
                                           data.frame(Population=p,
                                                      component=seq(1, best_fit@k),
                                                      alpha=best_fit@prior,
                                                      mu=unname(t(parameters(best_fit)))[,1],
                                                      stdev=unname(t(parameters(best_fit)))[,2]));
}
for (p in names(Browning_ODEN_bmms)) {
   if (!is.null(Browning_ODEN_bmms[[p]])) {
      k <- Browning_ODEN_bmms[[p]]$flexmix@k;
      if (k < 2) {
         mu <- plogis(coef(Browning_ODEN_bmms[[p]])[1]);
         nu <- exp(coef(Browning_ODEN_bmms[[p]])[2]);
      } else {
         mu <- plogis(coef(Browning_ODEN_bmms[[p]])[,1]);
         nu <- exp(coef(Browning_ODEN_bmms[[p]])[,2]);
      };
      alpha <- Browning_ODEN_bmms[[p]]$flexmix@prior;
      Browning_ODEN_bmm_bestfits <- bind_rows(Browning_ODEN_bmm_bestfits,
                                              data.frame(Population=p,
                                                         component=seq(1, k),
                                                         alpha=alpha,
                                                         shape1=unname(mu*nu),
                                                         shape2=unname((1-mu)*nu)));
   };
}

#Add the number of tracts used in each model fit, and a placeholder
# match rate column to make ggplot happy:
Browning_ODEN_gmm_bestfits <- Browning_ODEN_gmm_bestfits %>%
   left_join(browning_match_rates %>%
                filter(ONgood >= 30, ODgood >= 30,
                       ONmatchrate < 0.3, ODmatchrate > 0.3) %>%
                rename(Population=Group) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)
Browning_ODEN_bmm_bestfits <- Browning_ODEN_bmm_bestfits %>%
   left_join(browning_match_rates %>%
                filter(ONgood >= 30, ODgood >= 30,
                       ONmatchrate < 0.3, ODmatchrate > 0.3) %>%
                rename(Population=Group) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)

#Now plot the best fits against histograms of tract match rates:
gmm_bmm_plot(matchrates=browning_match_rates %>%
                select(-c(Nmatchrate, Ngood, Dmatchrate, Dgood)) %>%
                rename(Population=Group,
                       Nmatchrate=ONmatchrate,
                       Ngood=ONgood,
                       Dmatchrate=ODmatchrate,
                       Dgood=ODgood),
             gmm_fits=Browning_ODEN_gmm_bestfits,
             bmm_fits=Browning_ODEN_bmm_bestfits,
             pops=c("CDX", "CHB", "CHS", "JPT", "KHV",
                    "BEB", "GIH", "ITU", "PJL", "STU"),
             origin="DEN",
             prefix="Browning2018_ODENMM")
#Save the runs in case we want them later:
save(Browning_ODEN_gmms, Browning_ODEN_bmms,
     Browning_ODEN_gmm_bestfits, Browning_ODEN_bmm_bestfits,
     geom_gmm_density, geom_bmm_density, gmm_bmm_plot,
     file='Browning2018_ODENMM_PFRmatchrates_20230919.Rdata',
     compression_level=9,
     compress=TRUE)

#Clean up:
rm(browning_match_rates)
rm(Browning_DEN_gmms, Browning_DEN_bmms, Browning_ODEN_gmms, Browning_ODEN_bmms)
rm(Browning_DEN_gmm_bestfits, Browning_DEN_bmm_bestfits, Browning_ODEN_gmm_bestfits, Browning_ODEN_bmm_bestfits)
rm(alpha, best_fit, k, mu, nu, p, pops)
rm(geom_gmm_density, geom_bmm_density, gmm_bmm_plot)
