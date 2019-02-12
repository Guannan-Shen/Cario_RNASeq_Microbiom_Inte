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

####################### package info #############
sessionInfo()
citation("edgeR")
packageVersion("edgeR")

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

#############################################################
######### PCA analysis by prcomp is preferred #######
### prcomp is preferred for numerical accuracy #########
###########################################################


## 'princomp' can only be used with more units than variables
cnts.f.pca <- prcomp(t(cnts_f), center = TRUE, scale. = TRUE, retx = TRUE)
plot(cnts.f.pca)
##################### cnts.f.pca$x  ###############
### contains sample level pc components, a matrix which is sample by principal components
### from here, a new set of variable formed by the linear combination of original variables, such as pc1, pc2
## the pca plot is the value of each sample on the axis of different pcs 


############## cnts.f.pca$sdev, ###########
## a vector contains s.d. for each principal components,
## from here, calculate the percentage of variance holded by each pcs
cnts.f.sd <- ( cnts.f.pca$sdev/sum(cnts.f.pca$sdev) )* 100

########## pc names ##############
pc <- colnames(cnts.f.pca$x)

## ggplot2 to check the % of variance pca analysis 
pc.f.sd <- data.frame(PCs = pc,
                      Variance = cnts.f.sd)

p <- ggplot(data=pc.f.sd, aes(x=1:length(pc), y=Variance)) +
   geom_bar(stat="identity", fill = "Black") +
   ylab("% Variance") +
   scale_x_continuous(name = "Principal components", breaks = 1:)
   coord_flip() +
   theme_bw()
p

# Horizontal bar plot



