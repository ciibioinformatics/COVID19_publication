rm(list = ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
#BiocManager::install("DESeq2")
#BiocManager::install("SummarizedExperiment")

library(readr)
library(dplyr)
library(preprocessCore)

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

args <- vector()
args[1] <- "IgG_aggregate_file_covid_pub.csv"
args[2] <- "10000"

#Load aggregate files
aggregate_file <- read.csv(args[1])
#Filter aggregate file
aggregate_filtered <- aggregate_file[apply(aggregate_file[,18:(ncol(aggregate_file))]>=as.numeric(args[2]), 1, any),]

#For differential analysis between groups, functions from package "limma" have been used
#To be able to load values from aggregate file into DGEList object, make dataset containing probe sequences and signal intensities
#from filtered aggregate file
analysis_df <- aggregate_filtered[,c(13,18:ncol(aggregate_filtered))]
#Assign probe ids as row names since they are unique identifier
rownames(analysis_df) <- analysis_df$PROBE_ID
analysis_df$PROBE_ID <- NULL

#Calculate threshold from intensity values for random probe
sig_mat <- as.matrix(analysis_df)
test_nm <- normalize.quantiles(sig_mat,copy=TRUE)
test_bc <- rma.background.correct(test_nm,copy=TRUE)
test_bc <- as.data.frame(test_bc)

rownames(test_bc) <- rownames(analysis_df)
colnames(test_bc) <- colnames(analysis_df)
analysis_df <- test_bc
analysis_df <- ceiling(analysis_df)
length(unique(rownames(analysis_df)))
#Add probe sequences to analysis_df from sndf file
sndf_file <- read_delim("200129_Columbia_HCOV_serochip_PEP.sndf", delim = "\t")
analysis_df$PROBE_SEQUENCE <- sndf_file$PROBE_SEQUENCE[match(rownames(analysis_df), sndf_file$PROBE_ID)]

#Laod Wuhan epitopes from significant peptides
wuhan_epitopes <- read.csv("wuhan_epitopes_from_sig_pep_IgG_covid_pub.csv")
#Load 29 shortlisted epitopes
#epitopes_29 <- read.csv("epitopes_29.csv")

#wuhan_epitopes_29 <- wuhan_epitopes[wuhan_epitopes$EPITOPE_SEQ %in% epitopes_29$EPITOPE_SEQ,]
wuhan_epitopes <- wuhan_epitopes[,2:9]
wuhan_epitopes <- left_join(wuhan_epitopes, analysis_df, by = "PROBE_SEQUENCE")
write.csv(wuhan_epitopes, "heat_map_data_all_epitopes_covid_pub.csv")


