#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, repos="https://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("SNPRelate")
check_package("tidyverse")
check_package("ggrepel")
check_package("RColorBrewer")
check_package("viridis")

#Read in arguments:
plink_prefix <- options[1]
output_prefix <- options[2]
#Optional arguments for PCA plotting details:
if (length(options) < 3 || is.na(options[3]) || !file.exists(options[3])) {
   stop("Unable to proceed with PCA, the sample group map doesn't exist or wasn't provided")
} else {
   sample_group_map <- options[3]
}
if (length(options) < 4 || is.na(options[4]) || options[4] == '') {
   sample_group_map_layout <- 'ccccclcccnncccclc'
} else {
   sample_group_map_layout <- options[4]
}
if (length(options) < 5 || is.na(options[5]) || options[5] == '') {
   group_colname <- 'Region'
} else {
   group_colname <- options[5]
}
if (length(options) < 6 || is.na(options[6]) || options[6] == '') {
   group_levels <- c("America", "Europe", "Africa", "Middle_East",
                     "Central_South_Asia", "East_Asia", "Island_Southeast_Asia", "Oceania")
} else {
   group_levels <- unlist(strsplit(options[6], ",", fixed=TRUE))
   group_levels[group_levels == "NA"] <- NA_character_
}
if (length(options) < 7 || is.na(options[7]) || options[7] == '') {
   group_labels <- c("America", "Europe", "Africa", "Middle East",
                     "Central/South Asia", "East Asia", "Island Southeast Asia", "Oceania")
} else {
   group_labels <- unlist(strsplit(options[7], ",", fixed=TRUE))
}
if (length(options) < 9 || is.na(options[8]) || is.na(options[9]) || options[8] == '' || options[9] == '') {
   make_subset <- FALSE
} else {
   make_subset <- TRUE
   subset_colname <- options[8]
   subset_groupnames <- unlist(strsplit(options[9], ",", fixed=TRUE))
   print(paste("Keeping ", subset_colname, subset_groupnames, sep="", collapse="\n"))
}
if (length(options) < 11 || is.na(options[10]) || is.na(options[11]) || options[10] == '' || options[11] == '') {
   use_subgroups <- FALSE
} else {
   use_subgroups <- TRUE
   subgroup_colname <- options[10]
   subgroup_shapes <- as.integer(unlist(strsplit(options[11], ",", fixed=TRUE)))
#   print(paste0("Read in list of Subgroup shapes of length ", length(subgroup_shapes), ": ", paste(subgroup_shapes, sep=",", collapse=";")))
}

#Load the sample group map:
sample_groups <- read_tsv(sample_group_map,
                          col_types=sample_group_map_layout,
                          na=c("NA", ""))

if (make_subset) {
   sample_groups <- sample_groups %>%
      filter(.data[[subset_colname]] %in% subset_groupnames)
   keep_ids <- sample_groups$SampleID
   print(paste0("Subsetting ", length(keep_ids), " samples from ", length(subset_groupnames), " groups defined by column ", subset_colname))
} else {
   keep_ids <- sample_groups$SampleID
}

sample_groups <- sample_groups %>%
   mutate(Group=factor(.data[[group_colname]],
                       levels=group_levels,
                       labels=group_labels,
                       exclude=NULL))
print(paste0("After adding Group column (", group_colname, "), sample_groups has ", nrow(sample_groups), " rows and ", ncol(sample_groups), " columns"))
#print(paste(sample_groups$Group, sep=",", collapse=";"))

if (use_subgroups) {
   sample_groups <- sample_groups %>%
      mutate(Subgroup=factor(.data[[subgroup_colname]],
                             levels=unique(arrange(., Group, .data[[subgroup_colname]])[[subgroup_colname]])))
   print(paste0("After adding Subgroup column (", subgroup_colname, "), sample_groups has ", nrow(sample_groups), " rows and ", ncol(sample_groups), " columns"))
#   print(paste(sample_groups$Subgroup, sep=",", collapse=";"))
   subgroup_colour_indices <- (sample_groups %>%
      arrange(Group, Subgroup) %>%
      group_by(Group, Subgroup) %>%
      summarize(Colour=as.integer(first(Group))))$Colour
   palette_size <- length(unique(subgroup_colour_indices))
   if (palette_size <= 8) {
      subgroup_palette_vec <- rev(brewer.pal(palette_size, "Set2"))
   } else  {
      subgroup_palette_vec <- viridis_pal(option="turbo", end=0.9)(palette_size)
   }
   subgroup_colours <- subgroup_palette_vec[subgroup_colour_indices]
}

if (!file.exists(paste0(output_prefix, ".gds"))) {
   #Convert PLINK to GDS:
   snpgdsBED2GDS(paste0(plink_prefix, ".bed"),
                 paste0(plink_prefix, ".fam"),
                 paste0(plink_prefix, ".bim"),
                 paste0(output_prefix, ".gds"))
} else {
   print("Skipping PLINK to GDS conversion, since the GDS file already exists.")
}

if (!file.exists(paste0(output_prefix, "_PCA.Rdata"))) {
   #Open the file handle for the GDS:
   gds_fh <- snpgdsOpen(paste0(output_prefix, ".gds"))

   #Run a basic PCA:
   #eigen.cnt=0 returns all eigenvectors
   pca <- snpgdsPCA(gds_fh,
                    algorithm="exact",
                    eigen.cnt=0,
                    num.thread=1,
                    sample.id=keep_ids)

   #Postprocess the results into something easier to plot:
   eigvec <- pca$eigenvect
   rownames(eigvec) <- pca$sample.id
   colnames(eigvec) <- paste0("PC", seq(1, ncol(pca$eigenvect)))
   #And write them to files for later re-use:
   write.table(pca$eigenval,
               paste0(output_prefix, "_PCA.eigenval"),
               row.names=FALSE,
               col.names=FALSE,
               sep="\t",
               eol="\n",
               quote=FALSE)
   write.table(eigvec,
               paste0(output_prefix, "_PCA.eigenvec"),
               row.names=TRUE,
               col.names=TRUE,
               sep="\t",
               eol="\n",
               quote=FALSE)
   save(pca,
        file=paste0(output_prefix, "_PCA.Rdata"),
        compress=TRUE,
        compression_level=9)

   #Clean up and move on to the next SNPRelate analysis:
   rm(eigvec)

   rm(pca)
   #Close out the GDS file handle:
   snpgdsClose(gds_fh)
} else {
   print("Skipping PCA since the .Rdata file already exists")
}

#Reload the PCA results in more plottable format:
pca_eigval <- read.table(paste0(output_prefix, "_PCA.eigenval"),
                         header=FALSE,
                         colClasses=c("numeric"),
                         col.names=c("Eigenvalue"),
                         row.names=NULL)
pca_eigvec <- read.table(paste0(output_prefix, "_PCA.eigenvec"),
                         colClasses=c("character", rep("numeric", nrow(pca_eigval))),
                         col.names=c("ID", paste0("PC", seq(1, nrow(pca_eigval)))),
                         skip=1)
pca_eigval <- pca_eigval %>%
   mutate(PC=as.character(row_number()))
pca_eigvec <- pca_eigvec %>%
   pivot_longer(cols=-ID,
                names_prefix="PC",
                names_to="PC",
                values_to="Loading")

#Generate a data.frame of the percent of variance explained by each PC:
pve_df <- pca_eigval %>%
   mutate(PVE=abs(Eigenvalue)*100/sum(abs(Eigenvalue)))

#Add sample grouping labels:
#We assume that the sample group map has a sample ID column named "SampleID".
pca_df <- pca_eigvec %>%
   inner_join(sample_groups, by=c("ID"="SampleID"))

#Define the pca screeplot and biplot functions:
pca_screeplot <- function(pve) {
   pve %>%
      mutate(PC=as.numeric(PC)) %>%
      ggplot(aes(x=PC, y=PVE)) +
         geom_point() +
         theme_bw() +
         labs(y="Percent of Variance Explained") +
         ylim(0, NA)
}

pca_biplot <- function(loadings, pve, PCs, colour, shape=NULL, label=NULL, palette, shape_palette=NULL, shape_colours=NULL,
                       title, pointsize=1, labelsize=1, labellinewidth=0.5) {
   pve_precision <- 3
   pc_one <- paste0("PC", PCs[1])
   pc_one_pve <- signif(pve[pve$PC == PCs[1], "PVE"], pve_precision)
   pc_two <- paste0("PC", PCs[2])
   pc_two_pve <- signif(pve[pve$PC == PCs[2], "PVE"], pve_precision)
   n_groups <- loadings %>% count(.data[[colour]]) %>% nrow()
   loadings %>%
      filter(PC %in% PCs) %>%
      pivot_wider(names_from="PC",
                  names_prefix="PC",
                  values_from="Loading") %>%
      ggplot(data=., aes_string(x=pc_one, y=pc_two, colour=colour)) +
         {
            if (!is.null(shape)) {
               geom_point(aes_string(shape=shape), size=pointsize, alpha=0.5)
            } else {
               geom_point(size=pointsize, alpha=0.5)
            }
         } +
         {
            if (!is.null(label)) {
               geom_text_repel(aes_string(label=label), size=labelsize, segment.size=labellinewidth)
            }
         } +
         theme_bw() +
         {
            if (length(palette) > 1) {
               scale_colour_manual(values=palette)
            } else if (length(palette) == 1) {
               scale_colour_brewer(palette=palette, direction=-1)
            } else {
               if (n_groups <= 8) {
                  scale_colour_brewer(palette="Set2", direction=-1)
#               } else if (n_groups <= 12) {
#                  scale_colour_brewer(palette="Paired", direction=-1)
#               } else if (n_groups <= 13) {
#                  scale_colour_manual(values=c("black", "grey70", "deeppink1", "brown",
#                                               "darkred", "orange", "tomato", "blue",
#                                               "purple", "deepskyblue", "cornflowerblue", "darkblue",
#                                               "cyan"))
               } else {
                  scale_colour_viridis_d(option="turbo", end=0.9)
#                  scale_colour_discrete()
               }
            }
         } +
         {
            if (!is.null(shape_palette)) {
               scale_shape_manual(values=shape_palette)
            }
         } +
         labs(x=paste0(pc_one, " (", pc_one_pve, "%)"),
              y=paste0(pc_two, " (", pc_two_pve, "%)"),
              title=title) +
         guides(colour=guide_legend(direction="vertical",
                                    ncol=2,
                                    keywidth=0.7,
                                    keyheight=0.7,
                                    title=""),
                shape=guide_legend(direction="vertical",
                                   ncol=3,
                                   byrow=TRUE,
                                   override.aes=list(colour=shape_colours),
                                   title.theme=element_text(size=7),
                                   label.theme=element_text(size=6),
                                   keywidth=0.5,
                                   keyheight=0.5,
                                   title="")) +
         theme(legend.spacing.y=unit(0, "cm"))
}

#Generate the scree plot:
screeplot <- pca_screeplot(pve_df)
ggsave(paste0(output_prefix, "_PCA_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=screeplot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0(output_prefix, "_PCA_screeplot_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=screeplot,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

#Generate the PCA biplots:
biplot_1v2 <- pca_biplot(loadings=pca_df,
                         pve=pve_df,
                         PCs=c(1, 2),
                         colour="Group",
                         shape=ifelse(use_subgroups, "Subgroup", NULL),
                         label=NULL,
                         palette=NULL,
                         shape_palette=if(use_subgroups) subgroup_shapes else NULL,
                         shape_colours=if(use_subgroups) subgroup_colours else NULL,
                         title="")
ggsave(paste0(output_prefix, "_PCA_PC1vs2_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=biplot_1v2,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0(output_prefix, "_PCA_PC1vs2_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=biplot_1v2,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

biplot_3v4 <- pca_biplot(loadings=pca_df,
                         pve=pve_df,
                         PCs=c(3, 4),
                         colour="Group",
                         shape=ifelse(use_subgroups, "Subgroup", NULL),
                         label=NULL,
                         palette=NULL,
                         shape_palette=if(use_subgroups) subgroup_shapes else NULL,
                         shape_colours=if(use_subgroups) subgroup_colours else NULL,
                         title="")
ggsave(paste0(output_prefix, "_PCA_PC3vs4_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=biplot_3v4,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0(output_prefix, "_PCA_PC3vs4_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=biplot_3v4,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)

biplot_5v6 <- pca_biplot(loadings=pca_df,
                         pve=pve_df,
                         PCs=c(5, 6),
                         colour="Group",
                         shape=ifelse(use_subgroups, "Subgroup", NULL),
                         label=NULL,
                         palette=NULL,
                         shape_palette=if(use_subgroups) subgroup_shapes else NULL,
                         shape_colours=if(use_subgroups) subgroup_colours else NULL,
                         title="")
ggsave(paste0(output_prefix, "_PCA_PC5vs6_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".pdf"),
       plot=biplot_5v6,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
ggsave(paste0(output_prefix, "_PCA_PC5vs6_by", group_colname, "_", format(Sys.Date(), format="%Y%m%d"), ".ps"),
       plot=biplot_5v6,
       width=16.0,
       height=12.0,
       units="cm",
       dpi=500)
