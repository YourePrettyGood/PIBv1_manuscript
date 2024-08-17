#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

#Load the libraries:
check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, repos="https://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}
check_package("tidyverse")
check_package("flexmix")
check_package("betareg")
check_package("doParallel")
check_package("foreach")

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

#Read in arguments:
metadata <- options[1]
match_rates <- options[2]
output_prefix <- options[3]
prng_seed <- options[4]
max_K <- options[5]
num_threads <- options[6]

#Set the PRNG seed:
if (is.null(prng_seed)) {
   set.seed(42)
} else {
   set.seed(as.integer(prng_seed))
}

#Default max K is 10 if not provided:
if (is.null(max_K) || as.integer(max_K) < 1) {
   max_K <- 10
} else {
   max_K <- as.integer(max_K)
}

#Default number of threads is 1 if not provided:
if (is.null(num_threads) || as.integer(num_threads) < 1) {
   num_threads <- 1
} else {
   num_threads <- as.integer(num_threads)
}

#For now, we're hard-coding model selection to use BIC:
criterion <- "BIC"
#And the binwidth for the histogram at 0.02:
binwidth <- 0.02

#Load in the data:
cat("Generating pop->region map\n")
#Sample metadata file:
PIBv1_metadata <- read_tsv(metadata,
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
#Include the lumped populations:
lumped_pops <- data.frame(AnalysisGroup=c("Ambae,Maewo", "Bellona,Rennell",
                                          "Goroka,Sepik", "Nakanai,Mangseng"),
                          Region=c("Oceania", "Oceania",
                                   "Oceania", "Oceania"),
                          IslandGroup=c("Vanuatu", "Solomon Islands",
                                        "New Guinea", "New Britain"))
#Generate the population-to-region map:
pop_region_map <- bind_rows(PIBv1_metadata %>%
                               group_by(AnalysisGroup) %>%
                               summarize(Region=first(Region),
                                         IslandGroup=first(Island)),
                            lumped_pops)

#Load the Sprime tract match rates:
cat("Loading Sprime match rate file\n")
match_rates <- read_tsv(match_rates,
                        col_types='ciiciininininini') %>%
   separate(TractID,
            into=c("Population"),
            sep="_",
            remove=FALSE,
            extra="drop")
cat(paste0("Read in ", nrow(match_rates), " match rates from ", n_distinct(match_rates$Population), " populations\n"))

#Filter the Sprime tract match rates for contour plots:
#Also add population metadata
#The filters from Browning et al. 2018 for contour plots were:
# Ngood >= 10 | Dgood >= 10
cat("Adding population metadata and minimally filtering Sprime tracts (incl. removal of tracts with NA Denisovan match rate)\n")
filtered_match_rates <- match_rates %>%
   left_join(pop_region_map,
             by=c("Population"="AnalysisGroup")) %>%
   filter(Ngood >= 10 | Dgood >= 10) %>%
   filter(!is.na(Dmatchrate))
cat(paste0(nrow(filtered_match_rates), " Sprime tracts retained\n"))

#Set up the parallel processing backend:
registerDoParallel(num_threads)
#Run the GMMs and BMMs for each population:
#Make sure to adjust match rates that are 1.0 to very slightly less than 1 to avoid
# bounds errors in the BMM due to the link function.
pops <- (filtered_match_rates %>%
   summarize(Pops=unique(Population)))$Pops
#Denisovan_gmms <- list()
#Denisovan_bmms <- list()
Denisovan_gmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   cat(paste0(filtered_match_rates %>%
                 filter(Ngood >= 30, Dgood >= 30,
                        Nmatchrate < 0.3, Dmatchrate > 0.3,
                        Population == p) %>%
                 nrow(),
              " tracts retained for ", p, "\n"));
   cat(paste0("Running GMM for ", p, "\n"));
   gmm <- filtered_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate < 0.3, Dmatchrate > 0.3,
             Population == p) %>%
      mutate(Dmatchrate=case_when(Dmatchrate == 1.0 ~ 0.9999999,
                                  TRUE ~ Dmatchrate)) %>%
      stepFlexmix(Dmatchrate ~ 1,
                  data=.,
                  model=FLXMRglm(family="gaussian"),
                  k=seq(1, max_K, by=1));
   l <- list(gmm);
   names(l) <- p;
   l
}
Denisovan_bmms <- foreach(p=pops, .inorder=FALSE, .combine=c) %dopar% {
   cat(paste0(filtered_match_rates %>%
                 filter(Ngood >= 30, Dgood >= 30,
                        Nmatchrate < 0.3, Dmatchrate > 0.3,
                        Population == p) %>%
                 nrow(),
              " tracts retained for ", p, "\n"));
   cat(paste0("Running BMM for ", p, "\n"));
   bmm <- filtered_match_rates %>%
      filter(Ngood >= 30, Dgood >= 30,
             Nmatchrate < 0.3, Dmatchrate > 0.3,
             Population == p) %>%
      mutate(Dmatchrate=case_when(Dmatchrate == 1.0 ~ 0.9999999,
                                  TRUE ~ Dmatchrate)) %>%
      possibly(.f=betamix,
               otherwise=NULL)(Dmatchrate ~ 1,
                               data=.,
                               link="logit",
                               k=seq(1, max_K, by=1),
                               which=criterion,
                               verbose=TRUE);
   l <- list(bmm);
   names(l) <- p;
   l
}

#Determine the best fit K and reformat the model for use
# in plotting the mixture components:
cat("Identifying GMM best fit K and prepping mixture component data.frame\n")
Denisovan_gmm_bestfits <- data.frame(Population=c(),
                                     component=c(),
                                     alpha=c(),
                                     mu=c(),
                                     stdev=c())
#The Gaussian mixture model estimates coefficients that are the
# mean and standard deviation of the mixture components, so we
# store those per-component.
#"alpha" in this data.frame corresponds to the vector of prior
# probabilities of assignment to each cluster/component.
for (p in names(Denisovan_gmms)) {
   best_fit <- getModel(Denisovan_gmms[[p]], criterion);
   Denisovan_gmm_bestfits <- bind_rows(Denisovan_gmm_bestfits,
                                       data.frame(Population=p,
                                                  component=seq(1, best_fit@k),
                                                  alpha=best_fit@prior,
                                                  mu=unname(t(parameters(best_fit)))[,1],
                                                  stdev=unname(t(parameters(best_fit)))[,2]));
}
cat("Identifying BMM best fit K and prepping mixture component data.frame\n")
Denisovan_bmm_bestfits <- data.frame(Population=c(),
                                     component=c(),
                                     alpha=c(),
                                     shape1=c(),
                                     shape2=c())
#Note that for the beta mixture model, we have to convert between the
# beta parameterizations, as the mixture model uses the mean (mu) and
# ~precision (nu, aka phi) parameterization while the R dbeta function
# uses the shape1 and shape2 (alpha and beta) parameterization.
#Furthermore, the beta mixture model is a bit like a GLM where it has
# link functions between the coefficient estimates and the beta
# parameter values. We specify a logit link for the mean and use the
# default log link for the ~precision (phi), so the corresponding
# inverse link functions are plogis() and exp().
#Thus, the value of shape1 is plogis(b1)*exp(b2) and the value of
# shape2 is (1-plogis(b1))*exp(b2)
# where b1 and b2 are the two coefficients inferred by the mixture
# model.
#Also, for clarity, "alpha" in the data.frame is the vector of
# prior probabilities of assignment to a given cluster/component.
#It does *not* correspond to alpha in the R parameterization of
# the beta distribution.
for (p in names(Denisovan_bmms)) {
   if (!is.null(Denisovan_bmms[[p]])) {
      k <- Denisovan_bmms[[p]]$flexmix@k;
      if (k < 2) {
         mu <- plogis(coef(Denisovan_bmms[[p]])[1]);
         nu <- exp(coef(Denisovan_bmms[[p]])[2]);
      } else {
         mu <- plogis(coef(Denisovan_bmms[[p]])[,1]);
         nu <- exp(coef(Denisovan_bmms[[p]])[,2]);
      };
      alpha <- Denisovan_bmms[[p]]$flexmix@prior;
      Denisovan_bmm_bestfits <- bind_rows(Denisovan_bmm_bestfits,
                                          data.frame(Population=p,
                                                     component=seq(1, k),
                                                     alpha=alpha,
                                                     shape1=unname(mu*nu),
                                                     shape2=unname((1-mu)*nu)));
   };
}

#Add the tract counts and a dummy match rate column to the bestfits
# data.frames so they're compatible with the ggplot geom wrappers:
cat("Adding tract counts to best fit data.frames\n")
Denisovan_gmm_bestfits <- Denisovan_gmm_bestfits %>%
   left_join(filtered_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)

Denisovan_bmm_bestfits <- Denisovan_bmm_bestfits %>%
   left_join(filtered_match_rates %>%
                filter(Ngood >= 30, Dgood >= 30,
                       Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
                group_by(Population) %>%
                summarize(TractCount=n()),
             by="Population") %>%
   mutate(Dmatchrate=0.0)

#Now make the mixture model fit plots:
#One population selected per geographic region:
worldwide_pops <- c("French", "Maya",
                    "Hazara", "Dai",
                    "Rampasasa", "Baining-Kagat")
worldwide_pops_annot <- c("EUR-French", "AMR-Maya",
                          "CSA-Hazara", "EAS-Dai",
                          "ISEA-Rampasasa", "OCN-Baining-Kagat")
#GMM:
cat("Plotting worldwide GMM fits\n")
worldwide_den_gmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% worldwide_pops) %>%
   mutate(Population=factor(Population,
                            levels=worldwide_pops,
                            labels=worldwide_pops_annot)) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_gmm_density(Denisovan_gmm_bestfits %>%
                          filter(Population %in% worldwide_pops) %>%
                          mutate(Population=factor(Population,
                                                   levels=worldwide_pops,
                                                   labels=worldwide_pops_annot)),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
           y="Number of tracts")
ggsave(paste0(output_prefix, '_byRegionCuratedPops_Denisovan_GMMfits.pdf'),
       plot=worldwide_den_gmmfits,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
#BMM:
cat("Plotting worldwide BMM fits\n")
worldwide_den_bmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% worldwide_pops) %>%
   mutate(Population=factor(Population,
                            levels=worldwide_pops,
                            labels=worldwide_pops_annot)) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_bmm_density(Denisovan_bmm_bestfits %>%
                          filter(Population %in% worldwide_pops) %>%
                          mutate(Population=factor(Population,
                                                   levels=worldwide_pops,
                                                   labels=worldwide_pops_annot)),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
           y="Number of tracts")
ggsave(paste0(output_prefix, '_byRegionCuratedPops_Denisovan_BMMfits.pdf'),
       plot=worldwide_den_bmmfits,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Only PIB populations::
PIB_pops <- c("Ata", "Baining-Kagat",
              "Baining-Mali", "Kove",
              "Lavongai-Mussau", "Mamusi",
              "Mangseng", "Melamela",
              "Nailik-Notsi-Tigak", "Nakanai",
              "Saposa")
#GMM:
cat("Plotting PIB GMM fits\n")
PIB_den_gmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% PIB_pops) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_gmm_density(Denisovan_gmm_bestfits %>%
                          filter(Population %in% PIB_pops),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
      y="Number of tracts")
ggsave(paste0(output_prefix, '_byAnalysisGroupPIB_Denisovan_GMMfits.pdf'),
       plot=PIB_den_gmmfits,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
#BMM:
cat("Plotting PIB BMM fits\n")
PIB_den_bmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% PIB_pops) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_bmm_density(Denisovan_bmm_bestfits %>%
                          filter(Population %in% PIB_pops),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
      y="Number of tracts")
ggsave(paste0(output_prefix, '_byAnalysisGroupPIB_Denisovan_BMMfits.pdf'),
       plot=PIB_den_bmmfits,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)


#Oceania without Solomons (or Vanuatu):
OCNnoSolomons_pops <- c("Ata", "Baining-Kagat",
                        "Baining-Mali", "Goroka",
                        "Kove", "Lavongai-Mussau",
                        "Mamusi", "Mangseng",
                        "Melamela", "Nailik-Notsi-Tigak",
                        "Nakanai", "Nasioi",
                        "Saposa", "Sepik"))
#GMM:
cat("Plotting Oceania (no Solomons or Vanuatu) GMM fits\n")
OCNnoSolomons_den_gmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% OCNnoSolomons_pops) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_gmm_density(Denisovan_gmm_bestfits %>%
                          filter(Population %in% OCNnoSolomons_pops),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
      y="Number of tracts")
ggsave(paste0(output_prefix, '_byAnalysisGroupOCNnoSolomons_Denisovan_GMMfits.pdf'),
       plot=OCNnoSolomons_den_gmmfits,
       width=16.0,
       height=18.0,
       units="cm",
       dpi=500)
#BMM:
cat("Plotting Oceania (no Solomons or Vanuatu) BMM fits\n")
OCNnoSolomons_den_bmmfits <- filtered_match_rates %>%
   filter(Ngood >= 30, Dgood >= 30,
          Nmatchrate < 0.3, Dmatchrate > 0.3) %>%
   filter(Population %in% OCNnoSolomons_pops) %>%
   ggplot(aes(x=Dmatchrate, fill=Population)) +
      geom_histogram(binwidth=binwidth) +
      geom_bmm_density(Denisovan_bmm_bestfits %>%
                          filter(Population %in% OCNnoSolomons_pops),
                       binwidth=binwidth) +
      theme_bw() +
      scale_fill_brewer(palette="Spectral") +
      facet_wrap(~ Population,
                 ncol=2) +
      guides(fill=FALSE) +
      labs(x="Denisovan match rate",
      y="Number of tracts")
ggsave(paste0(output_prefix, '_byAnalysisGroupOCNnoSolomons_Denisovan_BMMfits.pdf'),
       plot=OCNnoSolomons_den_bmmfits,
       width=16.0,
       height=18.0,
       units="cm",
       dpi=500)

#Save the inputs and fits:
cat("Saving inputs and results to Rdata file\n")
save(PIBv1_metadata, lumped_pops, pop_region_map,
     match_rates, filtered_match_rates,
     Denisovan_gmms, Denisovan_gmm_bestfits,
     Denisovan_bmms, Denisovan_bmm_bestfits,
     geom_gmm_density, geom_bmm_density,
     file=paste0(output_prefix, "_Denisovan_gmms_bmms_k1thru", max_K, ".Rdata"),
     compress=TRUE,
     compression_level=9)

