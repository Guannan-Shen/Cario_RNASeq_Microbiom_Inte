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
setwd("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/")
## Asus
## setwd("C:/Users/hithr/Documents/Stats/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/")
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

############## cnts.f cnts.f.pca$rotation #######
## cnts.f.pca$rotation is the decomposition of ecah pc components at the gene level
# add the gene symbol
f.pca.symbol <- cnts.f.pca$rotation
row.names(f.pca.symbol) <- cnts_fsym$Symbol
f.pca.symbol <- as.matrix(f.pca.symbol)
############## cnts.f.pca$sdev, ###########
## a vector contains s.d. for each principal components,
## from here, calculate the percentage of variance holded by each pcs
cnts.f.sd <- ( cnts.f.pca$sdev/sum(cnts.f.pca$sdev) )* 100

########## pc names ##############
pc <- colnames(cnts.f.pca$x)

## ggplot2 to check the % of variance pca analysis 
# Horizontal bar plot
pc.f.sd <- data.frame(PCs = pc,
                      Variance = cnts.f.sd)

p.var.f <- ggplot(data=pc.f.sd, aes(x=1:length(pc), y=Variance)) +
   geom_bar(stat="identity", fill = "Black") +
   ylab("% Variance") +
   scale_x_continuous(name = "Principal components", breaks = 1:length(pc),
                      labels = pc)  +
  # Horizontal bar plot
   coord_flip() +
   theme_bw()

p.var.f
pc.f.sd$Variance[1:2]
## pc1 vs pc2
c("Healthy Controls", "HIV-Infected")

rownames(cnts.f.pca$x) == pheno$pid

pca12.f <- ggplot(data = data.frame(cnts.f.pca$x), aes(x = PC1, y = PC2, color = pheno$txt)) +
        geom_point(shape = 19) +
        xlab(paste("PC1 (", round(pc.f.sd$Variance[1],2) ,"%", ")", 
                   sep = "") ) +
        ylab(paste("PC2 (", round(pc.f.sd$Variance[2],2) ,"%", ")", 
             sep = "") ) +
       scale_color_manual(values=c("steelblue", "tomato"), 
                    name="Filtered Raw Counts\n\nGroups",
                    labels=c("Healthy Control",  "HIV-Infected")
                    ) +
       stat_ellipse(type = "t", linetype = "solid") +
      # annotate("text", x = max(df_male[,3])*0.85,
      #      y = quantile(df_male[,2], 0.95),
      #      label = paste(cpgname, gender, "M values",sep = "\n")) +
       theme_bw()
pca12.f
  
## 
## 'princomp' can only be used with more units than variables
cnts.f.pca <- prcomp(t(cnts_f), center = TRUE, scale. = TRUE, retx = TRUE)
plot(cnts.f.pca)
##################### cnts.f.pca$x  ###############
### contains sample level pc components, a matrix which is sample by principal components
### from here, a new set of variable formed by the linear combination of original variables, such as pc1, pc2
## the pca plot is the value of each sample on the axis of different pcs 

############## cnts.f cnts.f.pca$rotation #######
## cnts.f.pca$rotation is the decomposition of ecah pc components at the gene level
# add the gene symbol
f.pca.symbol <- cnts.f.pca$rotation
row.names(f.pca.symbol) <- cnts_fsym$Symbol
f.pca.symbol <- as.matrix(f.pca.symbol)
############## cnts.f.pca$sdev, ###########
## a vector contains s.d. for each principal components,
## from here, calculate the percentage of variance holded by each pcs
cnts.f.sd <- ( cnts.f.pca$sdev/sum(cnts.f.pca$sdev) )* 100

########## pc names ##############
pc <- colnames(cnts.f.pca$x)

## ggplot2 to check the % of variance pca analysis 
# Horizontal bar plot
pc.f.sd <- data.frame(PCs = pc,
                      Variance = cnts.f.sd)

p.var.f <- ggplot(data=pc.f.sd, aes(x=1:length(pc), y=Variance)) +
  geom_bar(stat="identity", fill = "Black") +
  ylab("% Variance") +
  scale_x_continuous(name = "Principal components", breaks = 1:length(pc),
                     labels = pc)  +
  # Horizontal bar plot
  coord_flip() +
  theme_bw()

p.var.f
pc.f.sd$Variance[1:2]
## pc1 vs pc2
c("Healthy Controls", "HIV-Infected")

rownames(cnts.f.pca$x) == pheno$pid

pca12.f <- ggplot(data = data.frame(cnts.f.pca$x), aes(x = PC1, y = PC2, color = pheno$txt)) +
  geom_point(shape = 19) +
  xlab(paste("PC1 (", round(pc.f.sd$Variance[1],2) ,"%", ")", 
             sep = "") ) +
  ylab(paste("PC2 (", round(pc.f.sd$Variance[2],2) ,"%", ")", 
             sep = "") ) +
  scale_color_manual(values=c("steelblue", "tomato"), 
                     name="Filtered Raw Counts\n\nGroups",
                     labels=c("Healthy Control",  "HIV-Infected")
  ) +
  stat_ellipse(type = "t", linetype = "solid") +
  # annotate("text", x = max(df_male[,3])*0.85,
  #      y = quantile(df_male[,2], 0.95),
  #      label = paste(cpgname, gender, "M values",sep = "\n")) +
  theme_bw()
pca12.f

################ deseq2
## data from another r script
cnts.deseq2.pca <- prcomp(t(cnts.deseq2), center = TRUE, scale. = TRUE, retx = TRUE)

deseq2.pca.symbol <- cnts.deseq2.pca$rotation
row.names(deseq2.pca.symbol) <- cnts_fsym$Symbol
deseq2.pca.symbol <- as.matrix(deseq2.pca.symbol)
############## cnts.deseq2.pca$sdev, ###########
## a vector contains s.d. for each principal components,
## from here, calculate the percentage of variance holded by each pcs
cnts.deseq2.sd <- ( cnts.deseq2.pca$sdev/sum(cnts.deseq2.pca$sdev) )* 100

########## pc names ##############
pc <- colnames(cnts.deseq2.pca$x)

## ggplot2 to check the % of variance pca analysis 
# Horizontal bar plot
pc.deseq2.sd <- data.frame(PCs = pc,
                      Variance = cnts.deseq2.sd)

## pc1 vs pc2
rownames(cnts.deseq2.pca$x) == pheno$pid

pca12.deseq2 <- ggplot(data = data.frame(cnts.deseq2.pca$x), aes(x = PC1, y = PC2, color = pheno$txt)) +
  geom_point(shape = 19) +
  xlab(paste("PC1 (", round(pc.deseq2.sd$Variance[1],2) ,"%", ")", 
             sep = "") ) +
  ylab(paste("PC2 (", round(pc.deseq2.sd$Variance[2],2) ,"%", ")", 
             sep = "") ) +
  scale_color_manual(values=c("steelblue", "tomato"), 
                     name="DESeq2 Normalized Counts\n\nGroups",
                     labels=c("Healthy Control",  "HIV-Infected")
  ) +
  stat_ellipse(type = "t", linetype = "solid") +
  # annotate("text", x = max(df_male[,3])*0.85,
  #      y = quantile(df_male[,2], 0.95),
  #      label = paste(cpgname, gender, "M values",sep = "\n")) +
  theme_bw()
pca12.deseq2

################ edger
## data from another r script
cnts.edger.pca <- prcomp(t(cnts.edger), center = TRUE, scale. = TRUE, retx = TRUE)

edger.pca.symbol <- cnts.edger.pca$rotation
row.names(edger.pca.symbol) <- cnts_fsym$Symbol
edger.pca.symbol <- as.matrix(edger.pca.symbol)
############## cnts.edger.pca$sdev, ###########

cnts.edger.sd <- ( cnts.edger.pca$sdev/sum(cnts.edger.pca$sdev) )* 100


## ggplot2 to check the % of variance pca analysis 
# Horizontal bar plot
pc.edger.sd <- data.frame(PCs = pc,
                           Variance = cnts.edger.sd)

## pc1 vs pc2
rownames(cnts.edger.pca$x) == pheno$pid

pca12.edger <- ggplot(data = data.frame(cnts.edger.pca$x), aes(x = PC1, y = PC2, color = pheno$txt)) +
  geom_point(shape = 19) +
  xlab(paste("PC1 (", round(pc.edger.sd$Variance[1],2) ,"%", ")", 
             sep = "") ) +
  ylab(paste("PC2 (", round(pc.edger.sd$Variance[2],2) ,"%", ")", 
             sep = "") ) +
  scale_color_manual(values=c("steelblue", "tomato"), 
                     name="edgeR Normalized Counts\n\nGroups",
                     labels=c("Healthy Control",  "HIV-Infected")
  ) +
  stat_ellipse(type = "t", linetype = "solid") +
  # annotate("text", x = max(df_male[,3])*0.85,
  #      y = quantile(df_male[,2], 0.95),
  #      label = paste(cpgname, gender, "M values",sep = "\n")) +
  theme_bw()
pca12.edger


################ tpm
## data from another r script
cnts.tpm.pca <- prcomp(t(cnts.tpm), center = TRUE, scale. = TRUE, retx = TRUE)

tpm.pca.symbol <- cnts.tpm.pca$rotation
row.names(tpm.pca.symbol) <- cnts_fsym$Symbol
tpm.pca.symbol <- as.matrix(tpm.pca.symbol)
############## cnts.tpm.pca$sdev, ###########

cnts.tpm.sd <- ( cnts.tpm.pca$sdev/sum(cnts.tpm.pca$sdev) )* 100


## ggplot2 to check the % of variance pca analysis 
# Horizontal bar plot
pc.tpm.sd <- data.frame(PCs = pc,
                          Variance = cnts.tpm.sd)

## pc1 vs pc2
rownames(cnts.tpm.pca$x) == pheno$pid

pca12.tpm <- ggplot(data = data.frame(cnts.tpm.pca$x), aes(x = PC1, y = PC2, color = pheno$txt)) +
  geom_point(shape = 19) +
  xlab(paste("PC1 (", round(pc.tpm.sd$Variance[1],2) ,"%", ")", 
             sep = "") ) +
  ylab(paste("PC2 (", round(pc.tpm.sd$Variance[2],2) ,"%", ")", 
             sep = "") ) +
  scale_color_manual(values=c("steelblue", "tomato"), 
                     name="TPM Normalized Counts\n\nGroups",
                     labels=c("Healthy Control",  "HIV-Infected")
  ) +
  stat_ellipse(type = "t", linetype = "solid") +
  # annotate("text", x = max(df_male[,3])*0.85,
  #      y = quantile(df_male[,2], 0.95),
  #      label = paste(cpgname, gender, "M values",sep = "\n")) +
  theme_bw()
pca12.tpm

##################### sig altered genes pca by gene lists ##########################





