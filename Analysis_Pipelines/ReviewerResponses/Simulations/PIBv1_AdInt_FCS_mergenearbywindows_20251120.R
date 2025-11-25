#!/usr/bin/env Rscript
#PIBv1 adaptive introgression method FCS and merging of nearby windows steps 2025/11/20:

#Load the libraries:
library(RcppRoll)
library(tidyverse)

#Helper functions:
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

#Output the intersection between an archaic core haplotype (*_a) and a set
# of selection scan windows (b_df):
#First filter for overlaps (and add archaic origin, then calculate intersection
# start and end conditional on overlap if requested (i.e. overlap=FALSE).
#Retains any extra annotation columns in b_df such as -log10(p_FCS)
intersect_windows <- function(chr_a, start_a, end_a, group_a, origin_a, b_df, overlap=FALSE) {
   overlap_df <- b_df %>%
      filter(group == group_a, chr == chr_a, start <= end_a, end >= start_a) %>%
      mutate(origin=origin_a)
   if (overlap) {
      overlap_df
   } else {
      overlap_df %>%
         mutate(start=pmax(start, start_a),
                end=pmin(end, end_a))
   }
}

#Merge windows separated by gaps of a maximum size:
merge_nearby_windows <- function(windows, max_gap=10000) {
   #The following is derived from Dani's merge_windows_within_10kb code:
   windows %>%
      group_by(group, origin, chr) %>%
      arrange(start) %>%
      mutate(group_id=cumsum(start > lag(end, default=first(start)) + max_gap)) %>%
      ungroup() %>%
      group_by(group, origin, chr, group_id) %>%
      summarize(start=min(start),
                end=max(end),
                p_FCS_cm=min(p_FCS_cm),
                neglog10_pFCScm=max(neglog10_pFCScm)) %>%
      ungroup() %>%
      dplyr::select(chr, start, end, group, origin, p_FCS_cm, neglog10_pFCScm)
}

#Pipeline parameters:
#Archaic core haplotype frequency quantile threshold to use:
corehap_quantile_threshold <- 0.95
#Number of threads to use:
#This isn't really used anywhere though
num_threads <- 4
#Window size:
w_size <- 20
#Window step size:
w_step <- 5
#Selection scan -log10(p_FCS) threshold:
selscan_thresh <- 2
#Maximum gap size for merging:
max_gap <- 10000

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

#Read in the full set of observed archaic core haplotypes:
cat("Reading in observed archaic core haplotype data\n")
obs_corehaps <- read_tsv('/gpfs/gibbs/pi/tucci/pfr8/Friedlaender/AnalysisFreezes/PIBv1_Sprime/PIBv1_Sprime_corehaps_wGeneList.bed',
                         col_select=c("#Chromosome", "Start", "End", "SprimePopulation", "TractID", "TractFreq")) %>%
   mutate(Start=Start+1,
          Origin=str_extract(TractID, "[A-Za-z0-9,-]+_\\d+_([A-Za-z]+)[.]\\d+", group=1)) %>%
   rename(chr=`#Chromosome`,
          start=Start,
          end=End)

for (AG in AGs) {
   gc()
   corehap_AG <- corehap_AGs[AG]
   #Extract the observed selection scan data:
   cat("Reading in observed", AG, "selection scan data\n")
   pbs_xpehh_fcs_df <- prep_obs_predictors(AG, w_size, w_step) %>%
      eval_FCS()

   #Extract windows under selection for this AG:
   cat("Filtering selection scan windows for", AG, "\n")
   AG_signif_selscan_windows <- pbs_xpehh_fcs_df %>%
      filter(neglog10_pFCScm >= selscan_thresh) %>%
      mutate(group=AG) %>%
      dplyr::select(chr, start, end, group, p_FCS_cm, neglog10_pFCScm)
   #Extract the observed corehaps for this AG and filter:
   cat("Filtering archaic core haplotypes for", AG, "\n")
   AG_corehaps <- obs_corehaps %>% filter(SprimePopulation == corehap_AG)
   corehap_freq_threshold <- quantile(AG_corehaps$TractFreq, corehap_quantile_threshold)
   AG_highfreq_corehaps <- AG_corehaps %>%
      filter(TractFreq >= corehap_freq_threshold) %>%
      mutate(group=AG) %>%
      dplyr::select(chr, start, end, group, Origin)

   #Take the intersection of the significant selscan windows and high frequency
   # core haplotypes:
   cat("Intersecting selection scan windows with -log10(p_FCS) >=", selscan_thresh, "and core haplotypes with frequency >=", corehap_freq_threshold, "for", AG, "\n")
   AG_selscan_corehap_overlap <- AG_highfreq_corehaps %>%
      pmap(\(chr, start, end, group, Origin) intersect_windows(chr, start, end, group, Origin, b_df=AG_signif_selscan_windows, overlap=TRUE)) %>%
      bind_rows()
   AG_selscan_corehap_intersect <- AG_highfreq_corehaps %>%
      pmap(\(chr, start, end, group, Origin) intersect_windows(chr, start, end, group, Origin, b_df=AG_signif_selscan_windows, overlap=FALSE)) %>%
      bind_rows()

   cat("Outputting the overlapping windows as a BED-like file\n")
   AG_selscan_corehap_overlap %>%
      arrange(chr, start, end) %>%
      transmute(`#chr`=chr,
                start=start-1,
                end=end,
                pop=group,
                FCS=p_FCS_cm,
                log10_FCS=neglog10_pFCScm) %>%
      write_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_AdInt_overlapping_hits_neglog10pFCSge2.bed'),
                quote="none",
                escape="none")

   cat("Outputting the intersections as a BED-like file\n")
   AG_selscan_corehap_intersect %>%
      arrange(chr, start, end) %>%
      transmute(`#chr`=chr,
                start=start-1,
                end=end,
                pop=group,
                FCS=p_FCS_cm,
                log10_FCS=neglog10_pFCScm) %>%
      write_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_AdInt_intersected_hits_neglog10pFCSge2.bed'),
                quote="none",
                escape="none")

   cat("Outputting gap-merged overlapping windows as a BED-like file\n")
   AG_selscan_corehap_overlap %>%
      merge_nearby_windows(max_gap=max_gap) %>%
      arrange(chr, start, end) %>%
      transmute(`#chr`=chr,
                start=start-1,
                end=end,
                pop=group,
                origin=origin,
                FCS=p_FCS_cm,
                log10_FCS=neglog10_pFCScm) %>%
      write_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_AdInt_overlapping_hits_neglog10pFCSge2_merged10k.bed'),
                quote="none",
                escape="none")

   cat("Outputting selection scan windows and high frequency core haplotypes as separate files for bedtools intersect cross-check\n")
   AG_highfreq_corehaps %>%
      arrange(chr, start, end) %>%
      transmute(`#chr`=chr,
                start=start-1,
                end=end,
                pop=group,
                origin=Origin) %>%
      write_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_highfreq_corehaps.bed'),
                quote="none",
                escape="none")
   AG_signif_selscan_windows %>%
      arrange(chr, start, end) %>%
      transmute(`#chr`=chr,
                start=start-1,
                end=end,
                pop=group,
                FCS=p_FCS_cm,
                log10_FCS=neglog10_pFCScm) %>%
      write_tsv(paste0('PIBv1_OCN_scans/PIBv1_', AG, '_selscan_hits_neglog10pFCSge2.bed'),
                quote="none",
                escape="none")
   }
   #Clean up:
   rm(pbs_xpehh_fcs_df, AG_highfreq_corehaps, AG_signif_selscan_windows, AG_selscan_corehap_intersect, AG_selscan_corehap_overlap)
}
