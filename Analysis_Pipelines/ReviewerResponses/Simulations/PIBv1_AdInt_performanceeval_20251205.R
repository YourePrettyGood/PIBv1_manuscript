#!/usr/bin/env Rscript
#PIBv1 adaptive introgression method performance evaluation 2025/12/05:
#This script compares one simulation each of the original
# PapuansOutOfAfrica_10J19 model with variants where Nb_Papua is
# either 500 (Sim1b) or 1000 (Sim1c).

#Load the libraries:
library(RcppRoll)
library(tidyverse)
library(tidymodels)

#Helper functions:
#Generate multiple binary factor columns based on thresholds:
cut_binary <- function(y, thresholds, thresh_labels, ...) {
   binary_list <- list()
   for (i in seq_along(thresholds)) {
      binary_factor <- as.character(i)
      break_set <- c(-Inf, thresholds[i], Inf)
      label_set <- thresh_labels[c(1,i+1)]
      binary_list[[binary_factor]] <- cut(y,
                                          breaks=break_set,
                                          labels=label_set,
                                          ...)
   }
   bind_cols(binary_list)
}

#Read in observed values of PBSn1 and windowed_max(XP-EHH) for a given AG:
prep_obs_predictors <- function(AG, w_size, w_step) {
   #Read in the windowed PBSn1 scores:
   pbs_df <- read_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_pbs.tsv.gz'),
                      col_select=c("#CHROM", "STARTPOS", "ENDPOS", "PBSn1")) %>%
      rename(chr=`#CHROM`,
             start=STARTPOS,
             end=ENDPOS,
             pbs=PBSn1) %>%
      mutate(clippbs=case_when(pbs < 0 ~ 0,
                               TRUE ~ pbs)) %>%
      mutate(pbs_rank=rank(pbs, ties.method="average", na.last="keep"),
             clippbs_rank=rank(clippbs, ties.method="average", na.last="keep")) %>%
      mutate(pbs_quantile_rank=pbs_rank/max(pbs_rank, na.rm=TRUE),
             clippbs_quantile_rank=clippbs_rank/max(clippbs_rank, na.rm=TRUE)) %>%
      dplyr::select(-c(pbs_rank, clippbs_rank))
   #Read in the per-SNP XP-EHH scores:
   xpehh_df <- read_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_xpehh.tsv.gz'),
                        col_select=c("#CHROM", "POS", "rawXPEHH", "XPEHH")) %>%
      rename(chr=`#CHROM`,
             pos=POS,
             xpehh=rawXPEHH,
             z_xpehh=XPEHH) %>%
      mutate(xpehh_rank=rank(xpehh, ties.method="average", na.last="keep"),
             z_xpehh_rank=rank(z_xpehh, ties.method="average", na.last="keep")) %>%
      mutate(quantile_rank_xpehh=xpehh_rank/max(xpehh_rank, na.rm=TRUE),
             quantile_rank_z_xpehh=z_xpehh_rank/max(z_xpehh_rank, na.rm=TRUE)) %>%
      dplyr::select(-c(xpehh_rank, z_xpehh_rank))
   #Window the XP-EHH scores in the same way as PBSn1:
   xpehh_windowed <- xpehh_df %>%
      group_by(chr) %>%
      reframe(start=roll_min(pos, n=w_size, by=w_step, normalize=FALSE),
              end=roll_max(pos, n=w_size, by=w_step, normalize=FALSE),
              max_xpehh=roll_max(xpehh, n=w_size, by=w_step, normalize=FALSE),
              max_z_xpehh=roll_max(z_xpehh, n=w_size, by=w_step, normalize=FALSE),
              max_quantile_rank_xpehh=roll_max(quantile_rank_xpehh, n=w_size, by=w_step, normalize=FALSE),
              max_quantile_rank_z_xpehh=roll_max(quantile_rank_z_xpehh, n=w_size, by=w_step, normalize=FALSE)) %>%
      ungroup() %>%
      mutate(max_xpehh_rank=rank(max_xpehh, ties.method="average", na.last="keep"),
             max_z_xpehh_rank=rank(max_z_xpehh, ties.method="average", na.last="keep")) %>%
      mutate(max_xpehh_quantile_rank=max_xpehh_rank/max(max_xpehh_rank, na.rm=TRUE),
             max_z_xpehh_quantile_rank=max_z_xpehh_rank/max(max_z_xpehh_rank, na.rm=TRUE)) %>%
      dplyr::select(-c(max_xpehh_rank, max_z_xpehh_rank))
   #Now join PBSn1 and windowed_max(XP-EHH) by window and output:
   pbs_xpehh_df <- pbs_df %>%
      inner_join(xpehh_windowed, by=c("chr", "start", "end"))
   pbs_xpehh_df
}

#Calculate p_{FCS} and -log10(p_{FCS}) values for the observed data and add
# to the data.frame of windowed PBSn1 and XP-EHH:
eval_FCS <- function(predictors_df) {
   predictors_df %>%
      mutate(p_FCS=pchisq(-2*(log(1-pbs_quantile_rank)+log(1-max_xpehh_quantile_rank)), df=2*2, lower.tail=FALSE),
             p_FCS_c=pchisq(-2*(log(1-clippbs_quantile_rank)+log(1-max_xpehh_quantile_rank)), df=2*2, lower.tail=FALSE),
             p_FCS_z=pchisq(-2*(log(1-pbs_quantile_rank)+log(1-max_z_xpehh_quantile_rank)), df=2*2, lower.tail=FALSE),
             p_FCS_cz=pchisq(-2*(log(1-clippbs_quantile_rank)+log(1-max_z_xpehh_quantile_rank)), df=2*2, lower.tail=FALSE),
             p_FCS_m=pchisq(-2*(log(1-pbs_quantile_rank)+log(1-max_quantile_rank_xpehh)), df=2*2, lower.tail=FALSE),
             p_FCS_cm=pchisq(-2*(log(1-clippbs_quantile_rank)+log(1-max_quantile_rank_xpehh)), df=2*2, lower.tail=FALSE),
             p_FCS_zm=pchisq(-2*(log(1-pbs_quantile_rank)+log(1-max_quantile_rank_z_xpehh)), df=2*2, lower.tail=FALSE),
             p_FCS_czm=pchisq(-2*(log(1-clippbs_quantile_rank)+log(1-max_quantile_rank_z_xpehh)), df=2*2, lower.tail=FALSE)) %>%
      mutate(neglog10_pFCS=case_when(p_FCS == 0 ~ 8,
                                     TRUE ~ -log10(p_FCS)),
             neglog10_pFCSc=case_when(p_FCS_c == 0 ~ 8,
                                      TRUE ~ -log10(p_FCS_c)),
             neglog10_pFCSz=case_when(p_FCS_z == 0 ~ 8,
                                      TRUE ~ -log10(p_FCS_z)),
             neglog10_pFCScz=case_when(p_FCS_cz == 0 ~ 8,
                                       TRUE ~ -log10(p_FCS_cz)),
             neglog10_pFCSm=case_when(p_FCS_m == 0 ~ 8,
                                      TRUE ~ -log10(p_FCS_m)),
             neglog10_pFCScm=case_when(p_FCS_cm == 0 ~ 8,
                                       TRUE ~ -log10(p_FCS_cm)),
             neglog10_pFCSzm=case_when(p_FCS_zm == 0 ~ 8,
                                       TRUE ~ -log10(p_FCS_zm)),
             neglog10_pFCSczm=case_when(p_FCS_czm == 0 ~ 8,
                                        TRUE ~ -log10(p_FCS_czm)))
}

#Read in simulated archaic core haplotypes for a given simulation:
read_sim_corehaps <- function(simnum) {
   read_tsv(paste0('Sim', simnum, '/Sprime/PIBv1_Sprime_targetpop_corehap_freqs.tsv.gz'),
            col_select=c("Chromosome", "Start", "End", "MedianAF")) %>%
      rename(chr=Chromosome, start=Start, end=End, TractFreq=MedianAF) %>%
      mutate(Simulation=simnum)
}

#Read in simulated values of PBSn1 and windowed_max(XP-EHH) for a given simulation:
prep_sim_predictors <- function(simnum, w_size, w_step) {
   #Read in the windowed PBSn1 scores:
   sim_pbs_df <- read_tsv(paste0('Sim', simnum, '/Sprime/PIBv1_perPop_pbs_w', w_size, '_step', w_step, '.tsv.gz'),
                          col_select=c("#CHROM", "STARTPOS", "ENDPOS", "PBSn1")) %>%
      rename(chr=`#CHROM`,
             start=STARTPOS,
             end=ENDPOS,
             pbs=PBSn1) %>%
      mutate(clippbs=case_when(pbs < 0 ~ 0,
                               TRUE ~ pbs))
   #Read in the per-SNP XP-EHH scores:
   sim_xpehh_df <- read_tsv(paste0('Sim', simnum, '/Sprime/PIBv1_perPop_xpehh.tsv.gz'),
                            col_select=c("#CHROM", "POS", "rawXPEHH", "XPEHH")) %>%
      rename(chr=`#CHROM`,
             pos=POS,
             xpehh=rawXPEHH,
             z_xpehh=XPEHH)
   #Window the XP-EHH scores in the same way as PBSn1:
   sim_xpehh_windowed <- sim_xpehh_df %>%
      group_by(chr) %>%
      reframe(start=roll_min(pos, n=w_size, by=w_step, normalize=FALSE),
              end=roll_max(pos, n=w_size, by=w_step, normalize=FALSE),
              max_xpehh=roll_max(xpehh, n=w_size, by=w_step, normalize=FALSE),
              max_z_xpehh=roll_max(z_xpehh, n=w_size, by=w_step, normalize=FALSE))
   #Now join PBSn1 and windowed_max(XP-EHH) by window and output:
   sim_pbs_xpehh_df <- sim_pbs_df %>%
      inner_join(sim_xpehh_windowed, by=c("chr", "start", "end")) %>%
      mutate(Simulation=simnum)
   sim_pbs_xpehh_df
}

#Helper function taking a row representing a selection scan window (i.e.
# (chr_w, start_w, end_w, group_w)) and the data.frame of high frequency
# core haplotypes as arguments, and counting the number of 
count_overlapping_windows <- function(chr_a, start_a, end_a, group_a, b_df) {
   b_df %>%
      filter(group == group_a, chr == chr_a, start <= end_a, end >= start_a) %>%
      nrow()
}

#Training function for a simple linear interpolation model:
linterp.train <- function(s, xvar, yvar, distinct.only=FALSE) {
   if (inherits(s, "rsplit")) {
      train_df <- analysis(s)
   } else if (inherits(s, "data.frame")) {
      train_df <- s
   } else {
      cli::cli_abort("{.arg s} must be either an rsplit object or a data.frame (or derivative class like tbl).")
   }
   if (distinct.only) {
      train_df <- train_df %>%
         distinct(.data[[xvar]], .data[[yvar]])
   }
   approxfun(x=as.data.frame(train_df)[,xvar],
             y=as.data.frame(train_df)[,yvar],
             rule=2)
}

#Prediction function for simple linear interpolation model:
predict.linterp <- function(f, newdata) {
   f(newdata)
}

#Pipeline and evaluation parameters:
#Threshold(s) for -log10(p_{FCS}) to use for classification in confusion matrix:
neglog10FCS_thresholds <- c(2, 3, 4)
#Labels for threshold levels for classification:
neglog10FCS_thresh_labels <- c("no_hit", "weak_hit", "hit", "strong_hit")
#Archaic core haplotype frequency quantile threshold to use:
corehap_quantile_threshold <- 0.95
#Number of threads to use:
num_threads <- 4
#Window size:
w_size <- 20
#Window step size:
w_step <- 5
#PRNG seed for determining evaluation set:
prng_seed <- 42
#Number of folds for v-fold CV:
cv_v <- 5
#Number of bins to use for stratified v-fold CV:
n_bins <- 5
#Proportion of the inputs to use as evaluation set:
eval_set_size <- 0.2

#Metric sets to use for model evaluation:
regression_metrics <- metric_set(ccc, rsq, rsq_trad, mae, rmse, mape)

#Analysis groups in the observed data:
AGs <- c("Ata", "Baining-Kagat", "Baining-Mali", "Bellona-Rennell", "Kove", "Lavongai-Mussau", "Malaita", "Mamusi", "Melamela", "Nailik-Notsi-Tigak", "Nakanai-Mangseng", "Nasioi", "Santa-Cruz", "Saposa", "Sepik-Goroka", "Tikopia", "Vella-Lavella")
corehap_AGs <- c("Ata"="Ata",
                 "Baining-Kagat"="Baining-Kagat",
                 "Baining-Mali"="Baining-Mali",
                 "Bellona-Rennell"="Bellona,Rennell",
                 "Kove"="Kove",
                 "Lavongai-Mussau"="Lavongai-Mussau",
                 "Malaita"="Malaita",
                 "Mamusi"="Mamusi",
                 "Melamela"="Melamela",
                 "Nailik-Notsi-Tigak"="Nailik-Notsi-Tigak",
                 "Nakanai-Mangseng"="Nakanai,Mangseng",
                 "Nasioi"="Nasioi",
                 "Santa-Cruz"="Santa-Cruz",
                 "Saposa"="Saposa",
                 "Sepik-Goroka"="Goroka,Sepik",
                 "Tikopia"="Tikopia",
                 "Vella-Lavella"="Vella-Lavella")

#Simulations to include:
sim_nums <- c("1", "1b", "1c")
Nb_Papua <- data.frame(Simulation=sim_nums,
                       NbPapua=c(243L, 500L, 1000L))

AG_sim_stats <- list()

#Read in the full set of observed archaic core haplotypes:
cat("Reading in observed archaic core haplotype data\n")
obs_corehaps <- read_tsv('/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_Sprime/PIBv1_Sprime_corehaps_wGeneList.bed',
                         col_select=c("#Chromosome", "Start", "End", "SprimePopulation", "TractID", "TractFreq")) %>%
   mutate(Start=Start+1) %>%
   rename(chr=`#Chromosome`,
          start=Start,
          end=End)
#Read in the full set of simulated archaic core haplotypes:
cat("Reading in simulated archaic core haplotype data\n")
sim_corehaps <- bind_rows(lapply(sim_nums, read_sim_corehaps))
n_sim_corehaps <- sim_corehaps %>%
   group_by(Simulation) %>%
   summarize(sim_all_corehaps=n())
#Read in the full set of simulated selection scan scores:
cat("Reading in simulated selection scan windows\n")
sim_pbs_xpehh_df <- bind_rows(lapply(sim_nums, prep_sim_predictors, w_size=w_size, w_step=w_step))
n_sim_selscan_windows <- sim_pbs_xpehh_df %>%
   group_by(Simulation) %>%
   summarize(sim_all_windows=n())

for (AG in AGs) {
   gc()
   corehap_AG <- corehap_AGs[AG]
   #Extract the observed selection scan data:
   cat("Reading in observed", AG, "selection scan data\n")
   pbs_xpehh_fcs_df <- prep_obs_predictors(AG, w_size, w_step) %>%
      eval_FCS()
   #Model with clipping of negative PBSn1 to 0 and windowing on max rawXPEHH quantile rank:
   #Prep the observed selection scan data for model fitting and evaluation:
   cat("Preparing observed", AG, "selection scan data for model fitting and evaluation\n")
   set.seed(prng_seed)
   training_eval_split <- initial_split(pbs_xpehh_fcs_df %>%
                                           dplyr::select(clippbs, max_xpehh, clippbs_quantile_rank, max_quantile_rank_xpehh, neglog10_pFCScm),
                                        prop=1-eval_set_size,
                                        strata=neglog10_pFCScm)
   training_cv_set <- vfold_cv(training(training_eval_split),
                               strata=neglog10_pFCScm,
                               v=cv_v,
                               breaks=n_bins)
   eval_set <- testing(training_eval_split)

   #Set up and fit linterp models for basic metrics:
   cat("Fitting", cv_v, "fold cross-validation linear interpolation models for", AG, "clipped PBSn1 and windowed_max(quantile_rank(XP-EHH))\n")
   training_cv_models <- training_cv_set %>%
      rowwise() %>%
      mutate(.pbsn1_model=list(linterp.train(splits,
                                             xvar="clippbs",
                                             yvar="clippbs_quantile_rank",
                                             distinct.only=TRUE)),
             .max_xpehh_model=list(linterp.train(splits,
                                                 xvar="max_xpehh",
                                                 yvar="max_quantile_rank_xpehh",
                                                 distinct.only=TRUE)))

   training_cv_metrics <- training_cv_models %>%
      rowwise() %>%
      mutate(.pbsn1_metrics=list(splits %>%
         assessment() %>%
         mutate(.pred=predict.linterp(.pbsn1_model, clippbs)) %>%
         regression_metrics(truth=clippbs_quantile_rank,
                            estimate=.pred)),
             .max_xpehh_metrics=list(splits %>%
         assessment() %>%
         mutate(.pred=predict.linterp(.max_xpehh_model, max_xpehh)) %>%
         regression_metrics(truth=max_quantile_rank_xpehh,
                            estimate=.pred)),
             .neglog10FCS_metrics=list(splits %>%
         assessment() %>%
         mutate(.pbsn1_pred=predict.linterp(.pbsn1_model, clippbs),
                .max_xpehh_pred=predict.linterp(.max_xpehh_model, max_xpehh)) %>%
         mutate(p_FCS_cm=pchisq(-2*(log(1-.pbsn1_pred)+log(1-.max_xpehh_pred)), df=2*2, lower.tail=FALSE)) %>%
         mutate(.pred=case_when(p_FCS_cm == 0 ~ 8,
                                TRUE ~ -log10(p_FCS_cm))) %>%
         regression_metrics(truth=neglog10_pFCScm,
                            estimate=.pred)))
   cat(cv_v, "fold cross-validation results for", AG, "clipped PBSn1 linear interpolation model\n")
   training_cv_metrics %>%
      dplyr::pull(.pbsn1_metrics) %>%
      print()
   cat(cv_v, "fold cross-validation results for", AG, "windowed_max(quantile_rank(XP-EHH)) linear interpolation model\n")
   training_cv_metrics %>%
      dplyr::pull(.max_xpehh_metrics) %>%
      print()
   cat(cv_v, "fold cross-validation results for", AG, "-log10(p_FCS) derived from the above linear interpolation models\n")
   training_cv_metrics %>%
      dplyr::pull(.neglog10FCS_metrics) %>%
      print()

   cat("Fitting final linear interpolation models for", AG, "clipped PBSn1 and windowed_max(quantile_rank(XP-EHH))\n")
   pbsn1_final_model <- training_eval_split %>%
      linterp.train(xvar="clippbs",
                    yvar="clippbs_quantile_rank",
                    distinct.only=TRUE)
   max_xpehh_final_model <- training_eval_split %>%
      linterp.train(xvar="max_xpehh",
                    yvar="max_quantile_rank_xpehh",
                    distinct.only=TRUE)
   final_model_preds <- training_eval_split %>%
      assessment() %>%
      mutate(.pbsn1_pred=predict.linterp(pbsn1_final_model, clippbs),
             .max_xpehh_pred=predict.linterp(max_xpehh_final_model, max_xpehh)) %>%
      mutate(p_FCS_cm=pchisq(-2*(log(1-.pbsn1_pred)+log(1-.max_xpehh_pred)), df=2*2, lower.tail=FALSE)) %>%
      mutate(.neglog10FCS_pred=case_when(p_FCS_cm == 0 ~ 8,
                                         TRUE ~ -log10(p_FCS_cm))) %>%
      mutate(FCS_hits_thresh=cut_binary(neglog10_pFCScm,
                                        thresholds=neglog10FCS_thresholds,
                                        thresh_labels=neglog10FCS_thresh_labels,
                                        right=FALSE,
                                        include.lowest=TRUE,
                                        ordered_results=TRUE),
             .pred_class=cut_binary(.neglog10FCS_pred,
                                    thresholds=neglog10FCS_thresholds,
                                    thresh_labels=neglog10FCS_thresh_labels,
                                    right=FALSE,
                                    include.lowest=TRUE,
                                    ordered_results=TRUE)) %>%
      unpack(c(FCS_hits_thresh, .pred_class), names_sep="_")
   cat("Final results for", AG, "clipped PBSn1 linear interpolation model\n")
   final_model_preds %>%
      regression_metrics(truth=clippbs_quantile_rank,
                         estimate=.pbsn1_pred) %>%
      print(n=Inf, width=Inf)
   cat("Final results for", AG, "windowed_max(quantile_rank(XP-EHH)) linear interpolation model\n")
   final_model_preds %>%
      regression_metrics(truth=max_quantile_rank_xpehh,
                         estimate=.max_xpehh_pred) %>%
      print(n=Inf, width=Inf)
   cat("Final results for", AG, "-log10(p_FCS) derived from the above linear interpolation models\n")
   final_model_preds %>%
      regression_metrics(truth=neglog10_pFCScm,
                         estimate=.neglog10FCS_pred) %>%
      print(n=Inf, width=Inf)
   for (i in seq_along(neglog10FCS_thresholds)) {
      truth_col_name <- paste0("FCS_hits_thresh_", i)
      estimate_col_name <- paste0(".pred_class_", i)
      cat("Confusion matrix for evaluation set of", AG, "linear interpolation model combo using -log10(p_FCS) threshold of", neglog10FCS_thresholds[i], "\n")
      final_model_preds %>%
         conf_mat(truth=all_of(truth_col_name),
                  estimate=all_of(estimate_col_name)) %>%
         print(n=Inf, width=Inf)
   }

#   cat("Fitting full linear interpolation models for", AG, "PBSn1 and windowed_max(XP-EHH)\n")
#   pbsn1_full_model <- pbs_xpehh_fcs_df %>%
#      linterp.train(xvar="pbs",
#                    yvar="pbs_quantile_rank",
#                    distinct.only=TRUE)
#   max_xpehh_full_model <- pbs_xpehh_fcs_df %>%
#      linterp.train(xvar="max_xpehh",
#                    yvar="max_xpehh_quantile_rank",
#                    distinct.only=TRUE)
#   #No point in evaluating the fit for the full model, since it would be guaranteed perfect.

   #Extract the observed corehaps for this AG and filter:
   cat("Filtering and counting archaic core haplotypes for", AG, "\n")
   AG_corehaps <- obs_corehaps %>%
      filter(SprimePopulation == corehap_AG)
   n_AG_corehaps <- nrow(AG_corehaps)
   corehap_freq_threshold <- quantile(AG_corehaps$TractFreq, corehap_quantile_threshold)
   AG_highfreq_corehaps <- AG_corehaps %>%
      filter(TractFreq >= corehap_freq_threshold) %>%
      mutate(group=AG) %>%
      select(chr, start, end, group)
   n_AG_highfreq_corehaps <- nrow(AG_highfreq_corehaps)

   #Evaluate FP count for corehap_freq_threshold from simulations:
   sim_highfreq_corehaps <- sim_corehaps %>%
      filter(TractFreq >= corehap_freq_threshold) %>%
      rename(group=Simulation) %>%
      select(chr, start, end, group)
   n_sim_highfreq_corehaps <- sim_highfreq_corehaps %>%
      rename(Simulation=group) %>%
      group_by(Simulation) %>%
      summarize(sim_signif_corehaps=n())

   #Predict -log10(p_FCS) for all the simulated windows:
   cat("Predicting -log10(p_FCS) for simulated windows using the", AG, "linear interpolation model combo\n")
   sim_predictions_df <- sim_pbs_xpehh_df %>%
      mutate(clippbs_quantile_rank=predict.linterp(pbsn1_final_model, clippbs),
             max_quantile_rank_xpehh=predict.linterp(max_xpehh_final_model, max_xpehh)) %>%
      mutate(p_FCS_cm=pchisq(-2*(log(1-clippbs_quantile_rank)+log(1-max_quantile_rank_xpehh)), df=2*2, lower.tail=FALSE)) %>%
      mutate(neglog10_pFCScm=case_when(p_FCS_cm == 0 ~ 8,
                                       TRUE ~ -log10(p_FCS_cm)))

   for (thresh in neglog10FCS_thresholds) {
      #Evaluate FP count for selection scans with given threshold:
      n_AG_selscan_windows <- nrow(pbs_xpehh_fcs_df)
      AG_signif_selscan_windows <- pbs_xpehh_fcs_df %>%
         filter(neglog10_pFCScm >= thresh) %>%
         mutate(group=AG) %>%
         select(chr, start, end, group)
      n_AG_signif_selscan_windows <- nrow(AG_signif_selscan_windows)

      #Now do the overlap for the combined method:
      #First predict the -log10(p_{FCS}) of the simulations and identify significant
      # simulated windows:
      cat("Identifying simulated selection scan windows with -log10(p_FCS) >=", thresh, "using the regression model trained on", AG, "\n")
      AG_sim_signif_windows <- sim_predictions_df %>%
         filter(neglog10_pFCScm >= thresh) %>%
         rename(group=Simulation) %>%
         select(chr, start, end, group)
      n_sim_signif_selscan_windows <- AG_sim_signif_windows %>%
         rename(Simulation=group) %>%
         group_by(Simulation) %>%
         summarize(sim_signif_windows=n())
      #We already identified high frequency core haplotypes in the simulated data above
      # (sim_highfreq_corehaps), so now we can count overlaps between the two sets of
      # intervals:
      cat("Counting overlaps between simulated >=", round(corehap_freq_threshold*100, 2), "% freq. core haplotypes and simulated selection scan windows with -log10(p_FCS) >=", thresh, "(", AG, "thresholds )\n")
      n_sim_intersect_windows <- sim_highfreq_corehaps %>%
         rowwise() %>%
         mutate(n_overlap=count_overlapping_windows(chr, start, end, group, AG_sim_signif_windows)) %>%
         rename(Simulation=group) %>%
         group_by(Simulation) %>%
         summarize(sim_signif_overlaps=sum(n_overlap))
      #Do the same for the observed windows and high frequency core haplotypes:
      cat("Counting overlaps between observed >=", round(corehap_freq_threshold*100, 2), "% freq. core haplotypes and observed selection scan windows with -log10(p_FCS) >=", thresh, "(", AG, "thresholds )\n")
      n_AG_intersect_windows <- AG_highfreq_corehaps %>%
         pmap_int(count_overlapping_windows, AG_signif_selscan_windows) %>%
         sum()

      #Compile all the performance statistics for this AG:
      hit_counts <- data.frame(Simulation=sim_nums,
                               obs_all_corehaps=n_AG_corehaps,
                               obs_all_windows=n_AG_selscan_windows,
                               obs_signif_corehaps=n_AG_highfreq_corehaps,
                               obs_signif_windows=n_AG_signif_selscan_windows,
                               obs_signif_overlaps=n_AG_intersect_windows) %>%
         left_join(n_sim_corehaps, by="Simulation") %>%
         left_join(n_sim_selscan_windows, by="Simulation") %>%
         left_join(n_sim_highfreq_corehaps, by="Simulation") %>%
         left_join(n_sim_signif_selscan_windows, by="Simulation") %>%
         left_join(n_sim_intersect_windows, by="Simulation") %>%
         mutate(obs_all_overlaps=obs_all_windows,
                sim_all_overlaps=sim_all_windows) %>%
         pivot_longer(cols=-c(Simulation),
                      names_sep="_",
                      names_to=c("Source", "Signif", "Stage"),
                      values_to="Count") %>%
         pivot_wider(id_cols=c(Simulation, Stage),
                     names_from=c(Source, Signif),
                     names_sep="_",
                     values_from=Count) %>%
         mutate(across(!c(Simulation, Stage), ~ replace_na(.x, 0)))
      hit_counts %>%
         tibble() %>%
         print(n=Inf, width=Inf)
      AG_sim_stats[[paste(AG, thresh)]] <- hit_counts %>%
         mutate(FPR=sim_signif/sim_all,
                P_O=obs_signif/obs_all) %>%
         mutate(FDR=FPR/P_O) %>%
         dplyr::select(Simulation, Stage, FPR, FDR) %>%
         pivot_longer(cols=c(FPR, FDR),
                      names_to="Statistic",
                      values_to="Value") %>%
         mutate(Stage=case_when(Stage == "corehaps" ~ "HighFreqCorehap",
                                Stage == "windows" ~ "SelScan",
                                Stage == "overlaps" ~ "Combined",
                                TRUE ~ Stage)) %>%
         pivot_wider(id_cols=c(Simulation, Statistic),
                     names_from=Stage,
                     values_from=Value) %>%
         mutate(AnalysisGroup=AG,
                SelScanThresh=thresh) %>%
         left_join(Nb_Papua, by="Simulation") %>%
         dplyr::select(AnalysisGroup, NbPapua, SelScanThresh, Statistic,
                       HighFreqCorehap, SelScan, Combined)
#      AG_sim_stats[[paste(AG, thresh)]] <- data.frame(AnalysisGroup=AG,
#                                                      SelScanThresh=thresh,
#                                                      Statistic=c("FPR",
#                                                                  "FDR"),
#                                                      HighFreqCorehap=c(n_sim_highfreq_corehaps/n_sim_corehaps,
#                                                                        (n_sim_highfreq_corehaps/n_sim_corehaps)/(n_AG_highfreq_corehaps/n_AG_corehaps)),
#                                                      SelScan=c(n_sim_signif_selscan_windows/n_sim_selscan_windows,
#                                                                (n_sim_signif_selscan_windows/n_sim_selscan_windows)/(n_AG_signif_selscan_windows/n_AG_selscan_windows)),
#                                                      Combined=c(n_sim_intersect_windows/n_sim_selscan_windows,
#                                                                 (n_sim_intersect_windows/n_sim_selscan_windows)/(n_AG_intersect_windows/n_AG_selscan_windows)))
      AG_sim_stats[[paste(AG, thresh)]] %>%
         tibble() %>%
         print(n=Inf, width=Inf)
   }
   #Clean up:
   rm(pbs_xpehh_fcs_df, training_eval_split, training_cv_set, eval_set)
   rm(training_cv_models, training_cv_metrics)
   rm(pbsn1_final_model, max_xpehh_final_model, final_model_preds)
#   rm(pbsn1_full_model, max_xpehh_full_model)
   rm(sim_predictions_df, AG_sim_signif_windows)
}

bind_rows(AG_sim_stats) %>%
   write_tsv(paste0('PIBv1_adaptive_introgression_performance_eval_', format(Sys.Date(), format="%Y%m%d"), '.tsv'))

bind_rows(AG_sim_stats) %>%
   as_tibble() %>%
   print(n=Inf, width=Inf)

