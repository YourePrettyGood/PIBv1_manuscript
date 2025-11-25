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
check_package("cowplot")
check_package("viridis")
check_package("ape")

#Helper functions:
cluster_haps <- function(haps, dist_method="manhattan", cluster_method="single") {
   hap_dists <- dist(haps, method=dist_method)
   if (any(is.na(hap_dists))) {
      stop("Pairwise haplotype distances include some NAs, please check your input for missingness.")
   }
   hap_clustering <- hclust(hap_dists, method=cluster_method)
   hap_clustering
}
sort_haps <- function(haps, ref_hap, clustering=NULL) {
   ref_dists <- rowSums(abs(sweep(haps, 2, ref_hap)), na.rm=TRUE)
   if (is.null(clustering)) {
      order(ref_dists)
   } else {
      library(ape)
      order.dendrogram(reorder(as.dendrogram(clustering), ref_dists, agglo.FUN=min))
   }
}
order_haps <- function(hap_matrix, pop_map, regions, ref_hap, do_clustering=FALSE, dist_method="manhattan", cluster_method="single") {
   hap_order <- c()
   for (r in regions) {
      if (sum(pop_map$Region == r) > 0) {
         region_haps <- hap_matrix[pop_map$Region == r, , drop=FALSE]
         region_haps_clustered <- NULL
         if (do_clustering) {
            region_haps_clustered <- cluster_haps(region_haps,
                                                  dist_method=dist_method,
                                                  cluster_method=cluster_method)
         }
         region_haps_order <- sort_haps(region_haps,
                                        haplotype_matrix[ref_hap,],
                                        region_haps_clustered)
         hap_order <- c(hap_order, rownames(region_haps)[region_haps_order])
      }
   }
   hap_order
}

#Read in arguments:
haps <- options[1]
site_list <- options[2]
modern_metadata <- options[3]
archaic_metadata <- options[4]
locus <- options[5]
ref_hap <- options[6]
output_prefix <- options[7]

#Set some of the plotting parameters:
plot_width <- 16.0
plot_height <- 24.0
plot_units <- "cm"
plot_dpi <- 500
#Relative panel heights for archaic and modern:
rel_panel_heights <- c(1, 4)

#Settings for clustering of haplotypes:
do_clustering <- TRUE
dist_method <- "euclidean"
cluster_method <- "average"

#Region names and short codes:
#For now we hard-code this.
regions <- data.frame(name=c("DenisovaDenisovan", "AltaiNeandertal",
                             "VindijaNeandertal", "ChagyrskayaNeandertal",
                             "Oceania", "Island_Southeast_Asia",
                             "East_Asia", "Central_South_Asia",
                             "America", "Europe",
                             "Middle_East", "Africa"),
                      code=c("DEN", "Altai",
                             "Vindija", "Chagyrskaya",
                             "OCN", "ISEA",
                             "EAS", "CSA",
                             "AMR", "EUR",
                             "MDE", "AFR"))

#Load the metadata files and combine:
archaic <- read_tsv(archaic_metadata) %>%
   transmute(SampleID=Sample,
             Region=str_c(Subpopulation, Region),
             Population=Population)
modern <- read_tsv(modern_metadata,
                   col_types='ccccclcccnncccclc',
                   na=c("NA", "")) %>%
   select(SampleID, Region, Population)
metadata <- bind_rows(archaic, modern) %>%
   mutate(Region=factor(Region,
                        levels=regions$name,
                        labels=regions$code))

#Generate a plot for each locus:
plot_prefix <- paste0(output_prefix, "_", locus)
cat(paste0("Loading data for ", locus, " locus\n"))
#Load the locus-specific data:
#Note: These are all sorted by CHROM:POS which should work for most loci,
# but any loci near a change in POS order of magnitude will be out of order.
haplotype_df <- read_tsv(haps,
                         col_types=cols(`CHROM:POS`=col_character(),
                                        .default=col_integer())) %>%
   pivot_longer(cols=-`CHROM:POS`,
                names_to="HaplotypeID",
                values_to="Haplotype")
all_sites <- read_tsv(site_list,
                      col_types='cc',
                      col_names=c("CHROM", "POS"),
                      skip=0) %>%
   unite("CHROM:POS", CHROM, POS, sep=":")

#Set up a map from haplotype ID to region:
pop_map <- bind_rows(metadata %>%
                        mutate(HaplotypeID=str_c(SampleID, "_1")),
                     metadata %>%
                        mutate(HaplotypeID=str_c(SampleID, "_2"))) %>%
   filter(HaplotypeID %in% unique(haplotype_df$HaplotypeID))

#Determine the plotting order by sorting (and possibly clustering) by match
# to the reference haplotype (e.g. Denisova_1):
cat(paste0("Sorting and ordering haplotypes for ", locus, " based on reference haplotype ", ref_hap, "\n"))
haplotype_matrix <- haplotype_df %>%
   filter(HaplotypeID != "Ancestral") %>%
   pivot_wider(id_cols=HaplotypeID,
               names_from="CHROM:POS",
               values_from="Haplotype") %>%
   column_to_rownames(var="HaplotypeID") %>%
   as.matrix()

#Reorder the pop_map based on the haplotype matrix:
pop_map <- pop_map[order(match(pop_map$HaplotypeID, rownames(haplotype_matrix))), , drop=FALSE]

haplotype_order_sorted <- order_haps(hap_matrix=haplotype_matrix,
                                     pop_map=pop_map,
                                     regions=regions$code,
                                     ref_hap=ref_hap,
                                     do_clustering=do_clustering,
                                     dist_method=dist_method,
                                     cluster_method=cluster_method)

#Now (re-)calculate the divergences from the reference haplotype:
cat(paste0("(Re-)calculating divergences from reference haplotype ", ref_hap, " for ", locus, " locus\n"))
ref_dists <- rowSums(abs(sweep(haplotype_matrix, 2, haplotype_matrix[ref_hap,])), na.rm=TRUE)
dist_min <- min(ref_dists)
dist_max <- max(ref_dists)

#Make the ref divergence plot:
cat(paste0("Generating haplotype divergence plot for ", locus, " locus\n"))
haplotype_dists <- data.frame(HaplotypeID=rownames(haplotype_matrix),
                              Divergence=ref_dists) %>%
   left_join(pop_map, by="HaplotypeID") %>%
   mutate(HaplotypeID=factor(HaplotypeID,
                             levels=rev(haplotype_order_sorted)))
arc_div_plot <- haplotype_dists %>%
   filter(Region %in% c("DEN", "Altai", "Vindija", "Chagyrskaya")) %>%
   ggplot(aes(x=Divergence, y=HaplotypeID, colour=Divergence)) +
      geom_point(size=2) +
      theme_minimal() +
      scale_x_continuous(limits=c(dist_min, dist_max)) +
      scale_colour_viridis_c(option="C",
                             end=0.8,
                             direction=-1,
                             limits=c(dist_min, dist_max)) +
      facet_grid(rows=vars(Region),
                 scales="free_y",
                 space="free_y",
                 switch="both") +
      labs(x="",
           y="") +
      guides(colour="none") +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            strip.text.y.left=element_text(angle=0, hjust=1, vjust=0.5, face="bold", size=12),
            plot.margin=unit(c(0.25, 0, -0.5, 0), "cm"),
            panel.border=element_rect(fill=NA, colour="black"),
            panel.spacing.y=unit(0, "cm"),
            panel.grid=element_blank())
modern_div_plot <- haplotype_dists %>%
   filter(!(Region %in% c("DEN", "Altai", "Vindija", "Chagyrskaya",
                          "AMR", "CSA", "MDE"))) %>%
   ggplot(aes(x=Divergence, y=HaplotypeID, colour=Divergence)) +
      geom_point(size=0.5) +
      theme_minimal() +
      scale_x_continuous(limits=c(dist_min, dist_max)) +
      scale_colour_viridis_c(option="C",
                             end=0.8,
                             direction=-1,
                             limits=c(dist_min, dist_max)) +
      facet_grid(rows=vars(Region),
                 scales="free_y",
                 space="free_y",
                 switch="both") +
      labs(x=paste0("Divergence from ", ref_hap),
           y="") +
      theme(axis.text.y=element_blank(),
            strip.text.y.left=element_text(angle=0, hjust=0, vjust=0.5, face="bold", size=12),
            plot.margin=unit(c(0, 0, 0, 0), "cm"),
            panel.border=element_rect(fill=NA, colour="black"),
            panel.spacing.y=unit(0, "cm"),
            panel.grid=element_blank())
haplotype_div_plot <- plot_grid(arc_div_plot, modern_div_plot,
                                align="v",
                                axis="lr",
                                nrow=length(rel_panel_heights),
                                rel_heights=rel_panel_heights)
ggsave(paste0(plot_prefix, "_haplotype_divergence_plot_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=haplotype_div_plot,
       width=plot_width,
       height=plot_height,
       units=plot_units,
       dpi=plot_dpi)
ggsave(paste0(plot_prefix, "_haplotype_divergence_plot_", format(Sys.Date(), format="%Y%m%d"), ".png"),
       plot=haplotype_div_plot,
       width=plot_width,
       height=plot_height,
       units=plot_units,
       dpi=plot_dpi)
ggsave(paste0(plot_prefix, "_haplotype_divergence_plot_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=haplotype_div_plot,
       width=plot_width,
       height=plot_height,
       units=plot_units,
       dpi=plot_dpi)
#Clean up within loop variables:
rm(haplotype_df, all_sites)
rm(pop_map)
rm(haplotype_matrix, haplotype_order_sorted)
rm(ref_dists, dist_min, dist_max)
rm(arc_div_plot, modern_div_plot)
rm(haplotype_div_plot)

#Clean up:
rm(cluster_haps, sort_haps, order_haps)
rm(do_clustering, dist_method, cluster_method, regions)
rm(plot_width, plot_height, plot_units, plot_dpi)
rm(metadata)
