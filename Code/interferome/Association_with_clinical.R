######################
## Set up workspace
######################
rm(list = ls())
library(edgeR)
library(EDASeq)
library(DESeq2)
library(knitr)
library(tidyverse)
library(magrittr)
library(stats)
library(BiocParallel)
library(readxl)
library(openxlsx)
library(limma)
library(scater)
library(rgl)
library(pca3d)
library(EnhancedVolcano)
library(RColorBrewer)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
getwd()
## not in function
'%nin%' <- Negate('%in%')

# ######## clean memory ######################
# rm(list = ls())
# gc()
# is(dds)
# slotNames(dds)

######## set wd #############
## ubuntu 
## setwd("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/")
## Asus
setwd("C:/Users/hithr/Documents/Stats/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/")
getwd()

#########################################################
### Common Preprocesssing for this project #############
########################################################

########### IMPORT dataset #########
# import unnormalized counts table
cnts.raw <- read.delim("All_Sample_geneCounts_raw_counts.txt", header = TRUE, sep = "\t")

#### another dataset, treated HIV
cnts.treated.raw <- read.xlsx("HIV-1_infected_HAART_treated_Raw_Counts.xlsx")

## numeric counts table 
cnts <- cnts.raw %>% dplyr::select(-c(Symbol, Length ) ) %>% tibble::column_to_rownames("Gene_ID")
cnts <- as.matrix(cnts)

## samples and pheno
rna.pid <- colnames(cnts)
# now we have the common counts table 
## pheno
ctrl.id <-  colnames(cnts)[1:13]
ctrl.id
hiv.id <- colnames(cnts)[14:32]
hiv.id
## from dim() we know there are 32 samples
pheno <- data.frame(pid = rna.pid, txt = as.factor(c(rep("Control", 13), 
                                                     rep("HIV", 19) )) )
pheno$txt %<>% relevel("Control")

######## filter genes ###########
## add Symbol to the data  
cnts <- as.data.frame(cnts)
cnts$Symbol <- cnts.raw$Symbol
cnts$Length <- cnts.raw$Length
head(cnts)
##filter the raw data and check dim
cnts_fsym <- cnts[rowSums(cnts[, 1:32])>=(5*ncol(cnts[, 1:32])), ]
## should end up around 15 - 20K genes 
ngenes <- nrow(cnts_fsym)
paste("The number of remaining genes: ", ngenes, sep = '')
# all integer in cnts_f
######## filtered counts table ############## 
cnts_f <- as.matrix(cnts_fsym[, 1:32])

############ setup for EDA ##########
# using the function from EDASeq
set <- newSeqExpressionSet(as.matrix(cnts_f),
                           phenoData = data.frame(condition=as.factor(pheno$txt),
                                                  row.names=colnames(cnts_f)))

######### setup for scater plots #########
## cnts. list
cnts.list <- list(cnts_f, cnts.edger, cnts.deseq2, cnts.tpm)
cnts.names <- c("Filtered Counts", "TMM Normalized Counts",
                "DESeq2 Normalized Counts", "TPM Normalized Counts")
for (i in 1:4){
  j  = cnts.list[[i]]
  ## using scater, start with an object
  example_sce <- SingleCellExperiment(
    assays = list(counts = j),
    colData = rna_info)
  ## RLE plots
  print( scater::plotRLE(example_sce, exprs_values = "counts", exprs_logged = FALSE,
                         colour_by = "Group", style = "minimal") + scale_x_discrete("Samples", labels = rna.pid) + labs(caption = paste("RLE Plot:" ,cnts.names[i]) ) )
  ## pca plots
  example_sce <- SingleCellExperiment(
    assays = list(logcounts = j),
    colData = rna_info)
  print( scater::plotPCA(example_sce, colour_by = "Group", size_by = "Group") + labs(caption = paste(  "PCA Plot:",cnts.names[i]) ))
}
