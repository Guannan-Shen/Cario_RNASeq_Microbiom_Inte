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
