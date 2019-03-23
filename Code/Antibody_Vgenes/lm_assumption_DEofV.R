######### DE of V genes ##################
library(edgeR)
library(tidyverse)
library(readxl)
############### rnaseq
setwd("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/Vgenes/")
cnts_raw <- read_excel("HIV_Infected_vs_Uninfected_relaxed_mapping_geneCounts.xlsx")
head(cnts_raw)
############ v genes lists ###################
vgenes <- read_excel("/home/guanshim/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/Vgenes_list.xlsx")
head(vgenes)
######clinical and microbiome data #########
clin_micro <- read_excel("/home/guanshim/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/clinical_Microbiome_pre.xlsx")
head(clin_micro)
sample_size <- length(clin_micro$pid)
sample_size
## test the ids between rna-seq and clin microbiome data
sum( colnames(cnts_raw)[-c(1:3)] != clin_micro$pid )
colnames(cnts_raw)[-c(1:3)]
## get cnts for norm and DE 
## filter to guarantee cnts larger than 1 
cnts_fsym <- cnts_raw[rowSums(as.matrix(cnts_raw[, -c(1:3)]) )>=(3*ncol(cnts_raw[, 1:32])), ]
## processing of rnaseq data
cnts <- cnts_fsym %>% dplyr::select(-c(Symbol, Length ) ) %>% tibble::column_to_rownames("Gene_ID")
cnts <- as.matrix(cnts)
dim(cnts)
######## edgeR TMM ###############
## edgeR object
group <- c(rep("Control", 13), rep("HIV", 19))
y <- DGEList(counts= cnts, group=group)
## normalization 
y <- calcNormFactors(y, method = "TMM")
y$samples
#  TMM counts 
# If you run the cpm function on a DGEList object which contains TMM normalisation factors 
# then you will get TMM normalised counts. 
# Here is a snippet of the source code for the cpm function:
cnts.edger <- edgeR::cpm(y)
dim(cnts.edger) == dim(cnts)
# as data frame

cnts_fsym <- data.frame(cnts_fsym)
cnts.edger <- data.frame(cnts.edger)
cnts.edger$Symbol <- cnts_fsym$Symbol
cnts.edger[]
# adjusting for age gender


########################### check the lm assumption of single mutation  ############################