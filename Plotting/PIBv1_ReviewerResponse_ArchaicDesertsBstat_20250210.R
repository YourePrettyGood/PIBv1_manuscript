#PIBv1 archaic deserts B statistic distribution (reviewer response) 2025/01/03:

#Load the libraries:
library(tidyverse)
library(viridisLite)

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

#Get to the right directory:
pushd('[path redacted]/PIBv1_results/Sprime/')

#Load the data:
deserts_Bstat <- read_tsv('deserts/bStatistic_MurphyAndMcVicker_deserts_introgression_genomewide_hist_maxgap0.tsv')

#eCDF comparison:
Bstat_eCDFs_plot <- deserts_Bstat %>%
	mutate(source=factor(source,
								levels=c("genomewide", "nondeserts",
											"introgressed", "nonintrogressed",
											"deserts"),
								labels=c("genomewide", "non-deserts",
											"introgressed", "non-introgressed",
											"deserts"))) %>%
	ggplot(aes(x=Bstat, y=eCDF, colour=source)) +
	   geom_step() +
	   theme_bw() +
	   scale_colour_viridis_d(option="H", end=0.9) +
	   facet_wrap(~ BstatType, ncol=1)
ggsave(paste0('PIBv1_Bstat_deserts_introgressed_eCDFs_', format(Sys.Date(), format='%Y%m%d'), '.pdf'),
		 plot=Bstat_eCDFs_plot,
		 width=16.0,
		 height=12.0,
		 units="cm",
		 dpi=500)

#Cliff's delta from ePMFs:
cliff_delta <- function(nXgtY, nYgtX, nX, nY) {
	n1 <- sum(nX)
	n2 <- sum(nY)
	x1gty2 <- sum(nXgtY*nY)
	x1lty2 <- sum(nYgtX*nX)
	delta <- (x1gty2 - x1lty2)/(n1*n2)
	list(delta=delta, n1=n1, n2=n2)
}

#Mann-Whitney U test with cumulative counts:
#mwu <- function(nXgtY, nX, nY) {
#	n <- sum(nX)
#	m <- sum(nY)
#	ties <- nX*nY
#	U <- sum(nXgtY*nY) + sum(ties)/2
#	n_tot <- n+m
#	mu_U <- n*m/2
#	coef_var <- n*m/12
#	var_noties <- coef_var*(n_tot+1)
#	var_ties <- coef_var*sum(ties**3-ties)/(n_tot*(n_tot-1))
#	sd_U <- sqrt(var_noties - var_ties)
#	z_U <- (U - mu_U)/sd_U
#	list(U=U, z_U=z_U, p=pnorm(z_U, lower.tail=FALSE), n=n, m=m)
#}
deserts_Bstat_cumulcounts <- deserts_Bstat %>%
	group_by(BstatType, source) %>%
	pivot_wider(id_cols=c(BstatType, Bstat),
					names_from=source,
					values_from=c(count, eCDF)) %>%
	pivot_longer(cols=-c(BstatType, Bstat),
					 names_to=c("Statistic", "source"),
					 names_sep="_",
					 values_to="Value") %>%
	filter(Statistic != "eCDF") %>%
	transmute(BstatType=BstatType,
				 source=source,
				 Bstat=as.integer(Bstat),
				 count=replace_na(Value, 0)) %>%
	group_by(BstatType, source) %>%
	mutate(countgt=sum(count)-cumsum(count)) %>%
	pivot_wider(id_cols=c(BstatType, Bstat),
					names_from=source,
					values_from=c(count, countgt))
deserts_Bstat_cumulcounts %>%
	group_by(BstatType) %>%
	summarize(desert_nondesert_delta=cliff_delta(countgt_nondeserts, countgt_deserts, count_nondeserts, count_deserts)$delta,
				 desert_nondesert_n1=cliff_delta(countgt_nondeserts, countgt_deserts, count_nondeserts, count_deserts)$n1,
				 desert_nondesert_n2=cliff_delta(countgt_nondeserts, countgt_deserts, count_nondeserts, count_deserts)$n2,
				 nonintro_intro_delta=cliff_delta(countgt_introgressed, countgt_nonintrogressed, count_introgressed, count_nonintrogressed)$delta,
				 nonintro_intro_n1=cliff_delta(countgt_introgressed, countgt_nonintrogressed, count_introgressed, count_nonintrogressed)$n1,
				 nonintro_intro_n2=cliff_delta(countgt_introgressed, countgt_nonintrogressed, count_introgressed, count_nonintrogressed)$n2) %>%
	print(n=Inf, width=Inf)
#BstatType        desert_nondesert_delta                   nonintro_intro_delta
#McVicker         delta=0.288, n1=2383311214, n2=292527918 delta=0.0969, n1=1894710545, n2=781128590
#Murphy_CADD      delta=0.369, n1=2391832032, n2=292740973 delta=0.188, n1=1896877250, n2=787707897
#Murphy_phastCons delta=0.363, n1=2391832032, n2=292740973 delta=0.193, n1=1896877250, n2=787707897
#deserts_Bstat_cumulcounts %>%
#	group_by(BstatType) %>%
#	summarize(desert_nondesert_MWU_U=mwu(countgt_nondeserts, count_nondeserts, count_deserts)$U,
#				 desert_nondesert_MWU_n=mwu(countgt_nondeserts, count_nondeserts, count_deserts)$n,
#				 desert_nondesert_MWU_m=mwu(countgt_nondeserts, count_nondeserts, count_deserts)$m,
#				 desert_nondesert_MWU_zU=mwu(countgt_nondeserts, count_nondeserts, count_deserts)$z_U,
#				 desert_nondesert_MWU_p=mwu(countgt_nondeserts, count_nondeserts, count_deserts)$p,
#				 nonintro_intro_MWU_U=mwu(countgt_introgressed, count_introgressed, count_nonintrogressed)$U,
#				 nonintro_intro_MWU_n=mwu(countgt_introgressed, count_introgressed, count_nonintrogressed)$n,
#				 nonintro_intro_MWU_m=mwu(countgt_introgressed, count_introgressed, count_nonintrogressed)$m,
#				 nonintro_intro_MWU_zU=mwu(countgt_introgressed, count_introgressed, count_nonintrogressed)$z_U,
#				 nonintro_intro_MWU_p=mwu(countgt_introgressed, count_introgressed, count_nonintrogressed)$p) %>%
#	print(n=Inf, width=Inf)
#Mann-Whitney U test code with tie correction gives warnings about
# NaNs produced by sqrt().
#Results below came from code without tie correction.
#BstatType         desert_nondesert_MWU        nonintro_intro_MWU
#McVicker          U=4.49e17, zU=8056,  p=~0   U=7.91e17, zU=4173,  p=~0     
#Murphy_CADD       U=4.79e17, zU=10332, p=~0   U=8.69e17, zU=8084,  p=~0
#Murphy_phastCons  U=4.77e17, zU=10161, p=~0   U=8.73e17, zU=8306,  p=~0

#Clean up:
rm(deserts_Bstat, deserts_Bstat_cumulcounts)
rm(Bstat_eCDFs_plot)
rm(cliff_delta, mwu)
