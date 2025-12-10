#!/usr/bin/env Rscript
#PIBv1 adaptive introgression method p_FCS histograms before and after
# overlapping with high frequency core haplotypes

#Load the libraries:
library(tidyverse)

#Read in the p_FCS values for each population before and after
# overlapping with high frequency core haplotypes:
selscan_candidate_lofn <- Sys.glob("*_PIBv1_selection_pvalues.tsv.gz")
AdInt_candidate_lofn <- Sys.glob("*_PIBv1_selection_highfreqcorehapoverlap_pvalues.tsv.gz")
histo_inputs <- bind_rows(lapply(selscan_candidate_lofn,
                                 function(fn) {read_tsv(fn,
                                                        col_select=c("population", "p_FCS")) %>%
                                                  mutate(Step="FCS")}),
                          lapply(AdInt_candidate_lofn,
                                 function(fn) {read_tsv(fn,
                                                        col_select=c("population", "p_FCS")) %>%
                                                  mutate(Step="intersect(FCS, high freq. core haps)")}))

#Now generate the plot:
histo_plot <- histo_inputs %>%
   ggplot(aes(x=p_FCS, fill=Step)) +
      geom_histogram(aes(y=after_stat(density)),
                     binwidth=0.01,
                     alpha=0.5,
                     position="identity") +
      theme_bw() +
      facet_wrap(~ population,
                 ncol=5,
                 scales="free_y") +
      labs(x=expression(p[F[CS]]),
           fill="") +
      theme(legend.position="bottom")

ggsave('PIBv1_FigS66_AdIntscans_FCS_intersection_pvalue_histograms_20251210.pdf',
       plot=histo_plot,
       width=32.0,
       height=24.0,
       units="cm",
       dpi=500)
