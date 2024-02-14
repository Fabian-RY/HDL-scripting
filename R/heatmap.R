#! /usr/bin/env Rscript

library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

csv <- args[[1]]
out_path <- args[[2]]

data <- read.csv(csv, sep="\t")
rownames(data) <- data[, 1]
data$X <- NULL
heatmap <- pheatmap(data)

svg(out_path)
heatmap
dev.off()
