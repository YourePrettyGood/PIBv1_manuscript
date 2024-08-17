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
locus_data <- options[1]
modern_metadata <- options[2]
archaic_metadata <- options[3]
input_prefix <- options[4]
minmaf <- options[5]
output_prefix <- options[6]

input_suffix <- paste0("_regionspan_MAFge", minmaf, "_wAA_")

#Set some of the plotting parameters:
plot_width <- 32.0
plot_height <- 24.0
plot_units <- "cm"
plot_dpi <- 500

#Settings for clustering of haplotypes:
do_clustering <- TRUE
dist_method <- "euclidean"
cluster_method <- "average"

#Set up the data.frame of loci we're plotting:
#
loci <- read_tsv(locus_data,
                 col_types=cols(.default=col_character()))

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

#Set the colours of the haplotype plot:
na_colour <- "gray80"
#Colours to use for VEP subplot:
vep_colours <- c("white", "maroon")
#Colours to use for protein-coding genes track:
pcg_colours <- c("white", "blue")
#Colours to use for variant source subplot:
annot_colours <- c("white", "black")
#Colours to use when ancestral allele is found:
anc_colours <- c("white", "purple", "blue", "red", "lavender")
#Colours to use when no ancestral allele is found:
noanc_colours <- c("gray60", "gray40")

#Generate a plot for each locus:
for (i in seq.int(nrow(loci))) {
   locus <- loci$locus[i]
   ref_hap <- loci$refhap[i]
   variant_source_levels <- str_split(loci$variant_source_levels[i], ",")[[1]]
   variant_source_labels <- str_split(loci$variant_source_labels[i], ",")[[1]]
   rel_panel_heights <- as.numeric(str_split(loci$rel_panel_heights[i], ",")[[1]])
   plot_prefix <- paste0(output_prefix, "_", locus)
   cat(paste0("Loading data for ", locus, " locus\n"))
   #Load the locus-specific data:
   #Note: These are all sorted by CHROM:POS which should work for most loci,
   # but any loci near a change in POS order of magnitude will be out of order.
   haplotype_df <- read_tsv(paste0("haplotypes/", input_prefix, "_merged_", locus, input_suffix, "haps_forR.tsv.gz"),
                            col_types=cols(`CHROM:POS`=col_character(),
                                           .default=col_integer())) %>%
      pivot_longer(cols=-`CHROM:POS`,
                   names_to="HaplotypeID",
                   values_to="Haplotype")
   all_sites <- read_tsv(paste0("annotations/", input_prefix, "_", locus, "_modernMAFfiltered_sites.tsv"),
                         col_types='cc',
                         col_names=c("CHROM", "POS"),
                         skip=0) %>%
      unite("CHROM:POS", CHROM, POS, sep=":")
   cadd <- read_tsv(paste0("annotations/", input_prefix, "_merged_", locus, input_suffix, "CADD.tsv.gz"),
                    col_types='cn') %>%
      right_join(all_sites, by=c("CHROM:POS"))
   vep <- read_tsv(paste0("annotations/", input_prefix, "_merged_", locus, input_suffix, "VEP_worst.tsv.gz"),
                   col_types=cols(.default=col_character()),
                   col_names=paste0("V", 1:31),
                   skip=0) %>%
      transmute("CHROM:POS"=V1, Allele=V2, VariantType=V3, Severity=V4, GeneSymbol=V5) %>%
      mutate(Missense=str_detect(VariantType, fixed("missense_variant", ignore_case=TRUE)),
             UTR=str_detect(VariantType, fixed("UTR_variant", ignore_case=TRUE)),
             CodingLoF=str_detect(VariantType, regex("frameshift_variant|stop_gained|stop_lost|start_lost", ignore_case=TRUE)),
             Regulatory=str_detect(VariantType, regex("regulatory|TF", ignore_case=TRUE)),
             Splice=str_detect(VariantType, fixed("splice", ignore_case=TRUE))) %>%
      right_join(all_sites, by=c("CHROM:POS"))
   annot <- read_tsv(paste0("annotations/", input_prefix, "_", locus, "_site_presence_matrix.tsv.gz"),
                     col_types=cols(.default=col_character())) %>%
      right_join(all_sites, by=c("CHROM:POS"))
   pcgs <- read_tsv(paste0("annotations/", input_prefix, "_", locus, "_modernMAFfiltered_PCGs.tsv.gz"),
                    col_types=cols(.default=col_character())) %>%
      right_join(all_sites, by=c("CHROM:POS")) %>%
      mutate(Exon=!is.na(Strand),
             Strand=case_when(Strand == "+" ~ ">",
                              Strand == "-" ~ "<",
                              TRUE ~ NA_character_)) %>%
      group_by(Feature) %>%
      arrange(`CHROM:POS`) %>%
      mutate(GeneLabel=case_when(`CHROM:POS` == nth(`CHROM:POS`, 1) ~ str_c(Feature, Strand),
                                 TRUE ~ NA_character_)) %>%
      ungroup()

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

   #Prep a plotting dataframe that is sorted by region and match to ref_hap:
   #We polarize alleles based on ancestral vs. derived, and derived alleles
   # are further polarized based on match to any Denisovan, Neandertal,
   # both, or neither.
   #If we only wanted to polarize based on ancestral vs. derived, use the
   # following mutate:
   #   mutate(across(-`CHROM:POS`,
   #                 ~ case_when(is.na(.data[["Ancestral"]]) ~ .x + 2L,
   #                             TRUE ~ as.integer(.x != .data[["Ancestral"]])))) %>%
   cat(paste0("Polarizing haplotypes for ", locus, " locus\n"))
   haplotype_plotting_df <- haplotype_df %>%
      pivot_wider(id_cols=`CHROM:POS`,
                  names_from="HaplotypeID",
                  values_from="Haplotype") %>%
      mutate(Shared=case_when(!is.na(.data[["Denisova_1"]]) & !is.na(.data[["Vindija33.19_1"]]) & .data[["Denisova_1"]] == .data[["Vindija33.19_1"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["Vindija33.19_1"]]) & .data[["Denisova_2"]] == .data[["Vindija33.19_1"]] ~ .data[["Denisova_2"]],
                              !is.na(.data[["Denisova_1"]]) & !is.na(.data[["Vindija33.19_2"]]) & .data[["Denisova_1"]] == .data[["Vindija33.19_2"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["Vindija33.19_2"]]) & .data[["Denisova_2"]] == .data[["Vindija33.19_2"]] ~ .data[["Denisova_2"]],
                              !is.na(.data[["Denisova_1"]]) & !is.na(.data[["AltaiNeandertal_1"]]) & .data[["Denisova_1"]] == .data[["AltaiNeandertal_1"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["AltaiNeandertal_1"]]) & .data[["Denisova_2"]] == .data[["AltaiNeandertal_1"]] ~ .data[["Denisova_2"]],
                              !is.na(.data[["Denisova_1"]]) & !is.na(.data[["AltaiNeandertal_2"]]) & .data[["Denisova_1"]] == .data[["AltaiNeandertal_2"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["AltaiNeandertal_2"]]) & .data[["Denisova_2"]] == .data[["AltaiNeandertal_2"]] ~ .data[["Denisova_2"]],
                              !is.na(.data[["Denisova_1"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) & .data[["Denisova_1"]] == .data[["Chagyrskaya-Phalanx_1"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) & .data[["Denisova_2"]] == .data[["Chagyrskaya-Phalanx_1"]] ~ .data[["Denisova_2"]],
                              !is.na(.data[["Denisova_1"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) & .data[["Denisova_1"]] == .data[["Chagyrskaya-Phalanx_2"]] ~ .data[["Denisova_1"]],
                              !is.na(.data[["Denisova_2"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) & .data[["Denisova_2"]] == .data[["Chagyrskaya-Phalanx_2"]] ~ .data[["Denisova_2"]],
                              TRUE ~ NA_integer_)) %>%
      mutate(DEN=case_when(is.na(.data[["Shared"]]) & !is.na(.data[["Denisova_1"]]) ~ .data[["Denisova_1"]],
                           !is.na(.data[["Shared"]]) & !is.na(.data[["Denisova_1"]]) & .data[["Shared"]] != .data[["Denisova_1"]] ~ .data[["Denisova_1"]],
                           is.na(.data[["Shared"]]) & !is.na(.data[["Denisova_2"]]) ~ .data[["Denisova_2"]],
                           !is.na(.data[["Shared"]]) & !is.na(.data[["Denisova_2"]]) & .data[["Shared"]] != .data[["Denisova_2"]] ~ .data[["Denisova_2"]],
                           TRUE ~ NA_integer_),
             NEAv=case_when(is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) ~ .data[["Vindija33.19_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) & .data[["Shared"]] != .data[["Vindija33.19_1"]] ~ .data[["Vindija33.19_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) ~ .data[["Vindija33.19_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) & .data[["Shared"]] != .data[["Vindija33.19_2"]] ~ .data[["Vindija33.19_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) ~ .data[["AltaiNeandertal_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_1"]] ~ .data[["AltaiNeandertal_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) ~ .data[["AltaiNeandertal_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_2"]] ~ .data[["AltaiNeandertal_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) ~ .data[["Chagyrskaya-Phalanx_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_1"]] ~ .data[["Chagyrskaya-Phalanx_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) ~ .data[["Chagyrskaya-Phalanx_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_2"]] ~ .data[["Chagyrskaya-Phalanx_2"]],
                            TRUE ~ NA_integer_),
             NEAc=case_when(is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) ~ .data[["Chagyrskaya-Phalanx_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_1"]] ~ .data[["Chagyrskaya-Phalanx_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) ~ .data[["Chagyrskaya-Phalanx_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_2"]] ~ .data[["Chagyrskaya-Phalanx_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) ~ .data[["Vindija33.19_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) & .data[["Shared"]] != .data[["Vindija33.19_1"]] ~ .data[["Vindija33.19_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) ~ .data[["Vindija33.19_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) & .data[["Shared"]] != .data[["Vindija33.19_2"]] ~ .data[["Vindija33.19_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) ~ .data[["AltaiNeandertal_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_1"]] ~ .data[["AltaiNeandertal_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) ~ .data[["AltaiNeandertal_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_2"]] ~ .data[["AltaiNeandertal_2"]],
                            TRUE ~ NA_integer_),
             NEAa=case_when(is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) ~ .data[["AltaiNeandertal_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_1"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_1"]] ~ .data[["AltaiNeandertal_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) ~ .data[["AltaiNeandertal_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["AltaiNeandertal_2"]]) & .data[["Shared"]] != .data[["AltaiNeandertal_2"]] ~ .data[["AltaiNeandertal_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) ~ .data[["Chagyrskaya-Phalanx_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_1"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_1"]] ~ .data[["Chagyrskaya-Phalanx_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) ~ .data[["Chagyrskaya-Phalanx_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Chagyrskaya-Phalanx_2"]]) & .data[["Shared"]] != .data[["Chagyrskaya-Phalanx_2"]] ~ .data[["Chagyrskaya-Phalanx_2"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) ~ .data[["Vindija33.19_1"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_1"]]) & .data[["Shared"]] != .data[["Vindija33.19_1"]] ~ .data[["Vindija33.19_1"]],
                            is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) ~ .data[["Vindija33.19_2"]],
                            !is.na(.data[["Shared"]]) & !is.na(.data[["Vindija33.19_2"]]) & .data[["Shared"]] != .data[["Vindija33.19_2"]] ~ .data[["Vindija33.19_2"]],
                            TRUE ~ NA_integer_)) %>%
      mutate(across(-c(`CHROM:POS`, `Ancestral`, `Shared`, `DEN`, `NEAv`, `NEAc`, `NEAa`),
                    ~ factor(case_when(is.na(.data[["Ancestral"]]) & !is.na(.x) & .x == 0L ~ 5L,
                                       is.na(.data[["Ancestral"]]) & !is.na(.x) ~ 6L,
                                       is.na(.x) ~ NA_integer_,
                                       .x == .data[["Ancestral"]] ~ 0L,
                                       !is.na(.data[["Shared"]]) & .x == .data[["Shared"]] ~ 1L,
                                       !is.na(.data[["DEN"]]) & .x == .data[["DEN"]] ~ 2L,
                                       !is.na(.data[["NEAv"]]) & .x == .data[["NEAv"]] ~ 3L,
                                       !is.na(.data[["NEAc"]]) & .x == .data[["NEAc"]] ~ 3L,
                                       !is.na(.data[["NEAa"]]) & .x == .data[["NEAa"]] ~ 3L,
                                       TRUE ~ 4L),
                             levels=0L:6L,
                             labels=c("Ancestral", "Shared", "DEN", "NEA", "NonArcDerived",
                                      "UnpolarizedREF", "UnpolarizedALT")))) %>%
      select(-c(Ancestral, Shared, DEN, NEAv, NEAc, NEAa)) %>%
      pivot_longer(cols=-`CHROM:POS`,
                   names_to="HaplotypeID",
                   values_to="Haplotype") %>%
      left_join(pop_map, by="HaplotypeID") %>%
      mutate(HaplotypeID=factor(HaplotypeID,
                                levels=rev(haplotype_order_sorted)))

   cat(paste0("Generating plot for ", locus, " locus\n"))
   #Subplot for the protein-coding genes track:
   pcg_track <- pcgs %>%
      right_join(annot %>% select(`CHROM:POS`), by="CHROM:POS") %>%
      ggplot(aes(x=`CHROM:POS`, y="Genes", fill=Exon, label=GeneLabel)) +
         geom_raster(alpha=0.5) +
         geom_text() +
         theme_minimal() +
         labs(x="",
              y="") +
         scale_fill_manual(values=pcg_colours) +
         guides(alpha="none",
                fill="none") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_text(face="bold", size=12),
               legend.position="none",
               plot.margin=unit(c(0, 0.5, -0.25, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Subplot for the CADD track:
   cadd_track <- cadd %>%
      right_join(annot %>% select(`CHROM:POS`), by="CHROM:POS") %>%
      ggplot(aes(x=`CHROM:POS`, y="CADD (PHRED)", fill=CADD)) +
         geom_raster() +
         theme_minimal() +
         labs(x="",
              y="") +
         scale_fill_viridis_c(option="magma",
                              na.value=na_colour) +
         guides(fill="none") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_text(face="bold", size=12),
               legend.position="top",
               plot.margin=unit(c(0, 0.5, -0.25, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Subplot for the variant effects tracks:
   variant_effects <- vep %>%
      select(`CHROM:POS`, Missense, UTR, CodingLoF, Splice, Regulatory) %>%
      pivot_longer(cols=-c(`CHROM:POS`),
                   names_to="VariantType",
                   values_to="Indicator") %>%
      ggplot(aes(x=`CHROM:POS`, y=VariantType, fill=Indicator)) +
         geom_raster() +
         theme_minimal() +
         labs(x="",
              y="") +
         scale_fill_manual(values=vep_colours,
                           na.value=na_colour,
                           drop=FALSE) +
         guides(fill="none") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_text(face="bold", size=12),
               plot.margin=unit(c(0, 0.5, -0.25, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Subplot for the variant sources:
   variant_sources <- annot %>%
      pivot_longer(cols=-c(`CHROM:POS`),
                   names_to="VariantSource",
                   values_to="Indicator") %>%
      mutate(VariantSource=factor(VariantSource,
                                  levels=variant_source_levels,
                                  labels=variant_source_labels),
             Indicator=case_when(Indicator == "1" ~ TRUE,
                                 Indicator == "0" ~ FALSE,
                                 TRUE ~ NA)) %>%
      ggplot(aes(x=`CHROM:POS`, y=VariantSource, fill=Indicator)) +
         geom_raster() +
         theme_minimal() +
         labs(x="",
              y="") +
         scale_fill_manual(values=annot_colours,
                           na.value=na_colour,
                           drop=FALSE) +
         guides(fill="none") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_text(face="bold", size=12),
               plot.margin=unit(c(0, 0.5, -0.25, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Subplot for just the archaic haplotypes:
   archaic_haplotypes <- haplotype_plotting_df %>%
      filter(Region %in% c("DEN", "Altai", "Vindija", "Chagyrskaya")) %>%
      ggplot(aes(x=`CHROM:POS`, y=HaplotypeID, fill=Haplotype)) +
         geom_raster() +
         theme_minimal() +
         scale_fill_manual(values=c(anc_colours, noanc_colours),
                           na.value=na_colour,
                           drop=FALSE) +
         facet_grid(rows=vars(Region),
                    scales="free_y",
                    space="free_y",
                    switch="both") +
         labs(x="",
              y="") +
         guides(fill="none") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               strip.text.y.left=element_text(angle=0, hjust=1, vjust=0.5, face="bold", size=12),
               plot.margin=unit(c(0, 0.5, -0.25, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Subplot for the modern haplotypes excluding x-axis labels:
   modern_haplotypes <- haplotype_plotting_df %>%
      filter(!(Region %in% c("DEN", "Altai", "Vindija", "Chagyrskaya",
                             "AMR", "CSA", "MDE"))) %>%
      ggplot(aes(x=`CHROM:POS`, y=HaplotypeID, fill=Haplotype)) +
         geom_raster() +
         theme_minimal() +
         scale_fill_manual(values=c(anc_colours, noanc_colours),
                           na.value=na_colour,
                           drop=FALSE) +
         facet_grid(rows=vars(Region),
                    scales="free_y",
                    space="free_y",
                    switch="both") +
         labs(x="",
              y="") +
         theme(axis.text.x=element_blank(),
               axis.text.y=element_blank(),
               strip.text.y.left=element_text(angle=0, hjust=0, vjust=0.5, face="bold", size=12),
               plot.margin=unit(c(0, 0.5, -0.5, 0), "cm"),
               panel.border=element_rect(fill=NA, colour="black"),
               panel.spacing.y=unit(0, "cm"))
   #Construct the final plot from subplots:
   haplotype_plot <- plot_grid(variant_effects, pcg_track, cadd_track, variant_sources,
                               archaic_haplotypes, modern_haplotypes,
                               align="v",
                               axis="lr",
                               nrow=length(rel_panel_heights),
                               rel_heights=rel_panel_heights)
   ggsave(paste0(plot_prefix, "_haplotype_plot_polarized_wVEP_wAnnot_wCADD_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
          plot=haplotype_plot,
          width=plot_width,
          height=plot_height,
          units=plot_units,
          dpi=plot_dpi)
   ggsave(paste0(plot_prefix, "_haplotype_plot_polarized_wVEP_wAnnot_wCADD_", format(Sys.Date(), format="%Y%m%d"), ".png"),
          plot=haplotype_plot,
          width=plot_width,
          height=plot_height,
          units=plot_units,
          dpi=plot_dpi)
   ggsave(paste0(plot_prefix, "_haplotype_plot_polarized_wVEP_wAnnot_wCADD_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
          plot=haplotype_plot,
          width=plot_width,
          height=plot_height,
          units=plot_units,
          dpi=plot_dpi)
   #Clean up within loop variables:
   rm(locus, ref_hap, plot_prefix, variant_source_levels, variant_source_labels, rel_panel_heights)
   rm(haplotype_df, all_sites, cadd, vep, annot, pcgs)
   rm(pop_map)
   rm(haplotype_matrix, haplotype_order_sorted, haplotype_plotting_df)
   rm(pcg_track, cadd_track, variant_effects, variant_sources, archaic_haplotypes, modern_haplotypes)
   rm(haplotype_plot)
}

#Clean up:
rm(cluster_haps, sort_haps, order_haps)
rm(do_clustering, dist_method, cluster_method, loci, regions)
rm(na_colour, vep_colours, pcg_colours, annot_colours, anc_colours, noanc_colours)
rm(plot_width, plot_height, plot_units, plot_dpi)
rm(archaic, modern, metadata)
rm(i)
