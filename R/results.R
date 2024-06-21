#! /usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(svglite)


args <- commandArgs(trailingOnly = TRUE)

#csv <- "/home/fabian/Escritorio/PhD/data/WP1/analisis/head-result.txt"
csv <- args[1]
print(csv)
correlation.data <- fread(csv, sep=args[2])
outfile_prefix <- args[3]
print(outfile_prefix)
width <- as.numeric(args[4])
height <- as.numeric(args[5])

heatmap_data = correlation.data[,1:3]

data_to_plot = as.data.frame(pivot_wider(heatmap_data,
                                         values_from = rg, names_from = p2))
rownames(data_to_plot) = lapply(data_to_plot$p1, basename)
rownames(data_to_plot) = sub("munged-", "", rownames(data_to_plot))
#rownames(data_to_plot) = sub("GCST[0-9]*_", "", rownames(data_to_plot))
rownames(data_to_plot) = sub(".sumstats.gz", "", rownames(data_to_plot))
data_to_plot$p1 = NULL
colnames(data_to_plot) = lapply(colnames(data_to_plot), basename)
colnames(data_to_plot) = sub("munged-", "", colnames(data_to_plot))
#colnames(data_to_plot) = sub("GCST[0-9]*_", "", colnames(data_to_plot))
colnames(data_to_plot) = sub(".sumstats.gz", "", colnames(data_to_plot))
#sum(!is.na(data_to_plot & data_to_plot > 1 | data_to_plot < -1))
discarded = sum(!is.na(data_to_plot & data_to_plot > 1 | data_to_plot < -1))
data_to_plot[, colSums(is.na(data_to_plot)) == nrow(data_to_plot)] = NULL
data_to_plot[data_to_plot > 1 | data_to_plot < -1] = NA
#discarded / nrow(data_to_plot)^2
heatmap <- pheatmap(data_to_plot, cluster_rows = T, cluster_cols = T,
                    show_rownames = T, show_colnames = T)
pval <- correlation.data[correlation.data$p < 0.05 & correlation.data$p > 0, ]
write.csv(pval, file=paste0(outfile_prefix, "_0.05.tsv"), sep = ",",
          col.names = T, row.names = F, quote = Finstall.packages("MendelianRandomization")
)

ggsave(paste0(outfile_prefix,"ggsave.png"), plot=heatmap)
ggsave(paste0(outfile_prefix,"ggsave.svg"), plot=heatmap)
ggsave(paste0(outfile_prefix,"ggsave.pdf"), plot=heatmap)
