#!/usr/bin/env Rscript

options <- commandArgs(trailingOnly=TRUE)

check_package <- function(pkg_name) {
   if (!require(pkg_name, character.only=TRUE)) {
      install.packages(pkg_name, repos="https://cran.us.r-project.org")
   }
   library(pkg_name, character.only=TRUE)
}

check_package("Ckmeans.1d.dp")

#Read in arguments:
maxdp_file <- options[1]
output_file <- options[2]

#Read in the input file:
maxdp <- read.table(maxdp_file,
                    colClasses=c("character", "integer"),
                    col.names=c("Sample", "maxDP"),
                    header=FALSE,
                    sep="\t")

#Determine the optimal number of clusters and extents of the clusters:
maxdp_clustered <- Ckmedian.1d.dp(maxdp$maxDP,
                                  k=c(1,20),
                                  y=1,
                                  method="loglinear",
                                  estimate.k="BIC")

#We need to output one file per cluster containing sample IDs:
maxdp_num_clusters <- max(maxdp_clustered$cluster)
for (i in 1:maxdp_num_clusters) {
   write.table(maxdp[maxdp_clustered$cluster == i, "Sample", drop=FALSE],
               paste0("maxdp", i, ".txt"),
               sep="\t",
               col.names=FALSE,
               row.names=FALSE,
               quote=FALSE)
}

#We also need to output a file that has one line per cluster with
# the sample ID filename for that cluster and the maxDP threshold
# to use (which we choose to be the median of the cluster):
maxdp_clusters_df <- data.frame(ClusterFiles=paste0("maxdp", 1:maxdp_num_clusters, ".txt"),
                                ClusterSamples=1:maxdp_num_clusters,
                                ClusterCenters=maxdp_clustered$centers)
write.table(maxdp_clusters_df[,c("ClusterFiles", "ClusterCenters")],
            output_file,
            sep="\t",
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE)

#Also just make a summary plot of the clustering:
pdf("maxDP_clustering.pdf",
    width=6.3,
    height=4.7,
    title=paste0("maxDP clustering with k=", maxdp_num_clusters))
plot(maxdp_clustered,
     main=paste0("k-medians clustering of maxDP for k=", maxdp_num_clusters))
dev.off()
