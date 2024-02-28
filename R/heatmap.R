#! /usr/bin/env Rscript

library(pheatmap)
library(peRReo)

args <- commandArgs(trailingOnly = TRUE)

csv <- "/media/fabian/Thot/Nextcloud/Documents/PhD/WP1/heatmap.csv"
out_path <- args[[2]]

data <- read.csv(csv, sep="\t")
rownames(data) <- data[, 1]
data$X <- NULL
heatmap <- pheatmap(data)

svg(out_path)
heatmap
dev.off()
