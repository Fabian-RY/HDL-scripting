#! /usr/bin/env Rscript

## Mendelian Randomization
##
##

library(MendelianRandomization)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

sums1 <- args[1] #"/home/fabian/Escritorio/PhD/data/WP1/analisis/rawdata/GCST90027158_Alzheimer-disease.tsv.gz" #  #
sums2 <- args[2] #"/home/fabian/Escritorio/PhD/data/WP1/analisis/rawdata/GCST90027164_Amyotrophic-lateral-sclerosis.tsv.gz" #
variantid1 <- args[3]
variantid2 <- args[4]

sumstats_1 <- fread(sums1)
sumstats_2 <- fread(sums2)

col1 <- sumstats_1[[variantid1]]
col2 <- sumstats_2[[variantid2]]

common <- intersect(col1, col2)

sumstats_2 <- sumstats_2[ col2 %in% common, ]
sumstats_1 <- sumstats_1[ col1 %in% common, ]

sumstats_1 <- sumstats_1[!duplicated(sumstats_1[[variantid1]]), ]
sumstats_2 <- sumstats_2[!duplicated(sumstats_2[[variantid2]]), ]


MR_input_1 <- mr_input(bx = sumstats_1$beta,
                       bxse = sumstats_1$standard_error,
                       by = sumstats_2$beta,
                       byse = sumstats_2$standard_error)

IVWObject <- mr_ivw(MR_input_1,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)

IVWObject
