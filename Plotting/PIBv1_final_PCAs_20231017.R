#PIBv1 PCA inspection

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

#Get to the directory:
pushd('[path redacted]/PIBv1_results/PCA/')

#Helper functions:
pca_screeplot <- function(pve, plot_prefix) {
   pve_plot <- pve %>%
      mutate(PC=as.numeric(PC)) %>%
      ggplot(aes(x=PC, y=PVE)) +
         geom_point() +
         theme_bw() +
         labs(y="Percent of Variance Explained") +
         ylim(0, NA)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=pve_plot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
}
pca_biplot <- function(loadings, pve, PCs, colour, shape, palette, shape_palette, title, plot_prefix, pointsize=1, label=NULL, labelsize=1, labellinewidth=0.5, theme="bw", section="supplement", plot_guide=TRUE) {
   if (!is.null(label)) {
      library(ggrepel)
   }
   pve_precision <- 3
   pc_one <- paste0("PC", PCs[1])
   pc_one_pve <- signif(pve[pve$PC == PCs[1], "PVE"], pve_precision)
   pc_two <- paste0("PC", PCs[2])
   pc_two_pve <- signif(pve[pve$PC == PCs[2], "PVE"], pve_precision)
   theme_func <- match.fun(paste0("theme_", theme))
   biplot <- loadings %>%
      filter(PC %in% PCs) %>%
      pivot_wider(names_from="PC",
                  names_prefix="PC",
                  values_from="Loading") %>%
      ggplot(data=., aes(x=.data[[pc_one]], y=.data[[pc_two]], colour=.data[[colour]], fill=.data[[colour]], shape=.data[[shape]])) +
         {
            if (section == "main") {
               geom_point(size=pointsize*1.2, alpha=1, colour="black", fill="black")
            }
         } +
         {
            if (section == "main") {
               geom_point(size=pointsize, alpha=1)
            } else {
               geom_point(size=pointsize, alpha=0.5)
            }
         } +
         {
            if (!is.null(label)) {
               if (is.data.frame(label)) {
                  geom_text_repel(data=label,
                                  aes(x=.data[[pc_one]], y=.data[[pc_two]], label=.data[[colour]]),
                                  size=labelsize,
                                  segment.size=labellinewidth,
                                  min.segment.length=0,
                                  force=10,
                                  inherit.aes=FALSE)
               } else {
                  geom_text_repel(aes(label=.data[[label]]), size=labelsize, segment.size=labellinewidth)
               }
            }
         } +
         theme_func() +
         scale_colour_manual(values=palette) +
         scale_fill_manual(values=palette) +
         scale_shape_manual(values=shape_palette) +
         labs(x=paste0(pc_one, " (", pc_one_pve, "%)"),
              y=paste0(pc_two, " (", pc_two_pve, "%)"),
              subtitle=title) +
         {
            if (plot_guide) {
               guides(shape=guide_legend(direction="vertical",
                                         ncol=3,
                                         byrow=TRUE,
                                         title.theme=element_text(size=7),
                                         label.theme=element_text(size=6),
                                         keywidth=0.5,
                                         keyheight=0.5,
                                         title=""),
                      colour=guide_legend(direction="vertical",
                                          ncol=3,
                                          byrow=TRUE,
                                          title.theme=element_text(size=7),
                                          label.theme=element_text(size=6),
                                          keywidth=0.5,
                                          keyheight=0.5,
                                          title=""),
                      fill=guide_legend(direction="vertical",
                                        ncol=3,
                                        byrow=TRUE,
                                        title.theme=element_text(size=7),
                                        label.theme=element_text(size=6),
                                        keywidth=0.5,
                                        keyheight=0.5,
                                        title=""))
            } else {
               guides(shape="none",
                      colour="none",
                      fill="none")
            }
         } +
         theme(legend.spacing.y=unit(0, "cm")) +
         {
            if (section == "main" && theme == "classic") {
               theme(panel.border=element_rect(fill=NA, colour="black"))
            } else {
               theme()
            }
         }
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
   ggsave(paste0(plot_prefix, "_", pc_one, "vs", pc_two, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=biplot,
          width=16.0,
          height=12.0,
          units="cm",
          dpi=500)
}

#Load the data:
#Sample metadata file:
PIBv1_metadata <- read_tsv('[path redacted]/Metadata/PIBv1_metadata_v0.3.tsv',
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

#Load the worldwide PCA results:
ww_pca_eigval <- read.table('worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA.eigenval',
                            header=FALSE,
                            colClasses=c("numeric"),
                            col.names=c("Eigenvalue"),
                            row.names=NULL) %>%
   mutate(PC=as.character(row_number()))
ww_pca_eigvec <- read.table('worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA.eigenvec',
                            colClasses=c("character", rep("numeric", nrow(ww_pca_eigval))),
                            col.names=c("ID", paste0("PC", seq(1, nrow(ww_pca_eigval)))),
                            skip=1) %>%
   pivot_longer(cols=-ID,
                names_prefix="PC",
                names_to="PC",
                values_to="Loading")

#Generate the data.frame of PVE per PC:
ww_pve_df <- ww_pca_eigval %>%
   mutate(PVE=abs(Eigenvalue)*100/sum(abs(Eigenvalue)))

#Generate the data.frame for plotting:
#Also do some reordering of factors to make the legends consistent.
ww_pca_df <- ww_pca_eigvec %>%
   inner_join(PIBv1_metadata, by=c("ID"="SampleID")) %>%
   mutate(Region=factor(Region,
                        levels=c("America", "Europe",
                                 "Africa", "Middle East",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"),
                        labels=c("AMR", "EUR",
                                 "AFR", "MDE",
                                 "CSA", "EAS",
                                 "ISEA", "OCN")),
          Island=factor(Island,
                        levels=c("Taiwan", "Philippines",
                                 "Flores", "New Guinea",
                                 "New Britain", "Mussau",
                                 "New Hanover", "New Ireland",
                                 "Bougainville", "Solomon Islands"))) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, Island, AnalysisGroup, Population)$Population)))

#Define the shapes and colours for the plots:
#First column is the shape, second column is the sample size for the region.
ww_shapes_suppl <- unlist(apply(cbind(c(8, 3,
                                        4, 6,
                                        5, 2,
                                        0, 1),
                                      c(5, 8,
                                        7, 4,
                                        9, 20,
                                        3, 23)),
                                1,
                                function(x) {rep(x[1], x[2])}))
#First column is sample size for the region, second is viridis palette, third is direction of palette.
ww_shape_colours_suppl <- unlist(apply(cbind(times=c(5, 8,
                                                     7, 4,
                                                     9, 20,
                                                     3, 23),
                                             palette=c("D", "D",
                                                       "D", "D",
                                                       "D", "D",
                                                       "D", "C"),
                                             c(1, 1,
                                               1, 1,
                                               1, 1,
                                               1, -1)),
                                       1,
                                       function(x) {viridis(as.integer(x[1]), end=0.9, option=x[2], direction=as.numeric(x[3]))}))
#Audrey's palettes:
ww_shapes <- c(21, 22, 23, 24, 25, #Americas
               21, 22, 23, 24, 25, 3, 4, 7, #Europe
               21, 22, 23, 24, 25, 3, 4, #Africa
               21, 22, 23, 24, #Middle East
               21, 22, 23, 24, 25, 3, 4, 7, 8, #Central/South Asia
               21, 8, 22, 23, 21, 24, 4, 24, 3, 4, 23, 7, 22, 25, 9, 3, 25, 7, 8, 10, #East Asia
               21, 22, 23, #Island Southeast Asia
               17, 16, #New Guinea
               15, 15, 16, 16, 17, 4, 3, 18, #New Britain
               15, 15, #Mussau and New Hanover
               15, 15, 15, #New Ireland
               15, 16, #Bougainville
               15, 17, 4, 18, 16, 3) #Solomon Islands
#Audrey's palettes with only open shapes:
#ww_shapes <- c(1, 0, 5, 2, 6, #Americas
#               1, 0, 5, 2, 6, 3, 4, 7, #Europe
#               1, 0, 5, 2, 6, 3, 4, #Africa
#               1, 0, 5, 2, #Middle East
#               1, 0, 5, 2, 6, 3, 4, 7, 8, #Central/South Asia
#               1, 8, 0, 5, 1, 2, 4, 2, 3, 4, 5, 7, 0, 6, 9, 3, 6, 7, 8, 10, #East Asia
#               1, 0, 5, #Island Southeast Asia
#               2, 1, #New Guinea
#               0, 0, 1, 1, 2, 4, 3, 5, #New Britain
#               0, 0, #Mussau and New Hanover
#               0, 0, 0, #New Ireland
#               0, 1, #Bougainville
#               0, 2, 4, 5, 1, 3) #Solomon Islands
ww_shape_colours <- c(rep("#F9B817", 5), #Americas
                      rep("#81B758", 8), #Europe
                      rep("#BE1E2D", 7), #Africa
                      rep("#DB7A24", 4), #Middle East
                      rep("#E292A8", 9), #Central/South Asia
                      case_when(c("S", "S", "S", "S", "N",
                                  "S", "N", "N", "S", "S",
                                  "N", "S", "N", "N", "S",
                                  "N", "S", "N", "N", "S") == "N" ~ "#842B6E",
                                TRUE ~ "#93678D"), #East Asia
                      rep("#5050FF", 3), #Island Southeast Asia
                      rep("#506CB3", 2), #New Guinea
                      case_when(c("W", "E", "E", "W", "W", "W", "W", "W") == "W" ~ "#41BBEC",
                                TRUE ~ "#3CC0C4"), #New Britain
                      rep("#90CBEF", 2), #Mussau and New Hanover
                      rep("#5897D2", 3), #New Ireland
                      rep("#7DDFEC", 2), #Bougainville
                      case_when(c("P", "S", "P", "S", "P", "S") == "P" ~ "#015D5D",
                                TRUE ~ "#018C8C")) #Solomon Islands

#Worldwide scree plot:
pca_screeplot(pve=ww_pve_df,
              plot_prefix="worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Supplementary figure S14 in the PIBv1 manuscript:
#PC1 vs. PC2 worldwide:
pca_biplot(loadings=ww_pca_df,
           pve=ww_pve_df,
           PCs=c(1,2),
           colour="Population",
           shape="Population",
           palette=ww_shape_colours_suppl,
           shape_palette=ww_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Supplementary figure S15 in the PIBv1 manuscript:
#PC3 vs. PC4 worldwide:
pca_biplot(loadings=ww_pca_df,
           pve=ww_pve_df,
           PCs=c(3,4),
           colour="Population",
           shape="Population",
           palette=ww_shape_colours_suppl,
           shape_palette=ww_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Supplementary figure S16 in the PIBv1 manuscript:
#PC5 vs. PC6 worldwide:
pca_biplot(loadings=ww_pca_df,
           pve=ww_pve_df,
           PCs=c(5,6),
           colour="Population",
           shape="Population",
           palette=ww_shape_colours_suppl,
           shape_palette=ww_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="worldwide/PIBv1_noVanuatu_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Clean up worldwide PCA:
rm(ww_pca_eigval, ww_pca_eigvec, ww_pve_df, ww_pca_df)
rm(ww_shapes, ww_shape_colours, ww_shapes_suppl, ww_shape_colours_suppl)

#Load the TaiwanISEAOCN PCA results:
ocn_pca_eigval <- read.table('TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA.eigenval',
                             header=FALSE,
                             colClasses=c("numeric"),
                             col.names=c("Eigenvalue"),
                             row.names=NULL) %>%
   mutate(PC=as.character(row_number()))
ocn_pca_eigvec <- read.table('TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA.eigenvec',
                             colClasses=c("character", rep("numeric", nrow(ocn_pca_eigval))),
                             col.names=c("ID", paste0("PC", seq(1, nrow(ocn_pca_eigval)))),
                             skip=1) %>%
   pivot_longer(cols=-ID,
                names_prefix="PC",
                names_to="PC",
                values_to="Loading")

#Generate the data.frame of PVE per PC:
ocn_pve_df <- ocn_pca_eigval %>%
   mutate(PVE=abs(Eigenvalue)*100/sum(abs(Eigenvalue)))

#Generate the data.frame for plotting:
#Also do some reordering of factors to make the legends consistent.
ocn_pca_df <- ocn_pca_eigvec %>%
   inner_join(PIBv1_metadata, by=c("ID"="SampleID")) %>%
   mutate(Region=factor(Region,
                        levels=c("America", "Europe",
                                 "Africa", "Middle East",
                                 "Central/South Asia", "East Asia",
                                 "Island Southeast Asia", "Oceania"),
                        labels=c("AMR", "EUR",
                                 "AFR", "MDE",
                                 "CSA", "EAS",
                                 "ISEA", "OCN")),
          Island=factor(Island,
                        levels=c("Taiwan", "Philippines",
                                 "Flores", "New Guinea",
                                 "New Britain", "Mussau",
                                 "New Hanover", "New Ireland",
                                 "Bougainville", "Solomon Islands"),
                        labels=c("Taiwan", "Philippines",
                                 "Flores", "New Guinea",
                                 "New Britain", "Mussau",
                                 "New Hanover", "New Ireland",
                                 "Bougainville", "Solomon Islands"))) %>%
   mutate(Population=factor(Population,
                            levels=unique(arrange(., Region, Island, AnalysisGroup, Population)$Population)))

#Define the shapes and colours for the plots:
#First column is the shape, second column is the sample size for the island(s).
ocn_shapes_suppl <- unlist(apply(cbind(c(8, 6, 7, 0, 1, 3, 2, 4, 5),
                                       c(2, 2, 1, 2, 8, 1+1, 3, 2, 6)),
                                 1,
                                 function(x) {rep(x[1], x[2])}))
#First column is sample size for the island(s), second is viridis palette, third is direction of palette.
ocn_shape_colours_suppl <- unlist(apply(cbind(times=c(2+2+1, 2, 8, 1+1+3, 2, 6),
                                              palette=c("D", "B", "C", "D", "D", "D"),
                                              c(1, 1, -1, 1, 1, 1)),
                                        1,
                                        function(x) {viridis(as.integer(x[1]), end=0.9, option=x[2], direction=as.numeric(x[3]))}))
#Audrey's palettes:
ocn_shapes <- c(21, 8, #Taiwan
                21, 22, 23, #Island Southeast Asia
                17, 16, #New Guinea
                15, 15, 16, 16, 17, 4, 3, 18, #New Britain
                15, 15, #Mussau and New Hanover
                15, 15, 15, #New Ireland
                15, 16, #Bougainville
                15, 17, 4, 18, 16, 3) #Solomon Islands
#Audrey's palettes with only open shapes:
#ocn_shapes <- c(1, 8, #Taiwan
#                1, 0, 5, #Island Southeast Asia
#                2, 1, #New Guinea
#                0, 0, 1, 1, 2, 4, 3, 5, #New Britain
#                0, 0, #Mussau and New Hanover
#                0, 0, 0, #New Ireland
#                0, 1, #Bougainville
#                0, 2, 4, 5, 1, 3) #Solomon Islands
ocn_shape_colours <- c(case_when(c("S", "S") == "N" ~ "#842B6E",
                                 TRUE ~ "#93678D"), #Taiwan
                       rep("#5050FF", 3), #Island Southeast Asia
                       rep("#506CB3", 2), #New Guinea
                       case_when(c("W", "E", "E", "W", "W", "W", "W", "W") == "W" ~ "#41BBEC",
                                 TRUE ~ "#3CC0C4"), #New Britain
                       rep("#90CBEF", 2), #Mussau and New Hanover
                       rep("#5897D2", 3), #New Ireland
                       rep("#7DDFEC", 2), #Bougainville
                       case_when(c("P", "S", "P", "S", "P", "S") == "P" ~ "#015D5D",
                                 TRUE ~ "#018C8C")) #Solomon Islands

#TaiwanISEAOCN scree plot:
pca_screeplot(pve=ocn_pve_df,
              plot_prefix="TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Figure 1 panel C of the PIBv1 manuscript:
#Note that font sizes and labels were tweaked in Adobe Illustrator
# for clarity of presentation.
#PC1 vs. PC2 TaiwanISEAOCN for main text:
pop_labels_df <- ocn_pca_df %>%
   filter(PC %in% c(1,2)) %>%
   pivot_wider(names_from="PC",
               names_prefix="PC",
               values_from="Loading") %>%
   group_by(Population) %>%
   summarize(PC1=mean(PC1),
             PC2=mean(PC2))
pca_biplot(loadings=ocn_pca_df,
           pve=ocn_pve_df,
           PCs=c(1,2),
           colour="Population",
           shape="Population",
           palette=ocn_shape_colours,
           shape_palette=ocn_shapes,
           title="",
           pointsize=2,
           label=pop_labels_df,
           labelsize=3,
           labellinewidth=0.25,
           theme="classic",
           section="main",
           plot_guide=FALSE,
           plot_prefix="TaiwanISEAOCN/PIBv1_Fig1C_TaiwanISEAOCN_PC1vsPC2")

#This plot actually got skipped because it's redundant with Fig. 1C:
#PC1 vs. PC2 TaiwanISEAOCN for supplement:
pca_biplot(loadings=ocn_pca_df,
           pve=ocn_pve_df,
           PCs=c(1,2),
           colour="Population",
           shape="Population",
           palette=ocn_shape_colours_suppl,
           shape_palette=ocn_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Supplementary figure S20 in the PIBv1 manuscript:
#PC3 vs. PC4 TaiwanISEAOCN:
pca_biplot(loadings=ocn_pca_df,
           pve=ocn_pve_df,
           PCs=c(3,4),
           colour="Population",
           shape="Population",
           palette=ocn_shape_colours_suppl,
           shape_palette=ocn_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#This plot also got skipped from the supplement:
#PC5 vs. PC6 TaiwanISEAOCN:
pca_biplot(loadings=ocn_pca_df,
           pve=ocn_pve_df,
           PCs=c(5,6),
           colour="Population",
           shape="Population",
           palette=ocn_shape_colours_suppl,
           shape_palette=ocn_shapes_suppl,
           title="",
           pointsize=2,
           plot_prefix="TaiwanISEAOCN/PIBv1_noVanuatu_TaiwanISEAOCN_autosomes_globalMAFge0.01_bSNPs_Choin_PCA")

#Clean up TaiwanISEAOCN PCA:
rm(ocn_pca_eigval, ocn_pca_eigvec, ocn_pve_df, ocn_pca_df)
rm(ocn_shapes, ocn_shape_colours, ocn_shapes_suppl, ocn_shape_colours_suppl)
rm(pop_labels_df)
