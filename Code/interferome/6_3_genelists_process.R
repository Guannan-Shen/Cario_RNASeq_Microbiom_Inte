## load in new 252 ISGs and 424 Beta
# using copy from excel in ubuntu and paste to Text Editor in Ubuntu 
# then read delim to get gene lists
# we have all the gene lists in the DataRaw dir: old: coreISG genesbeta and coreISG252 and genesbeta424
dir <- "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/"
isgs <- as.data.frame(read.delim( paste(dir, "coreISG252", sep = "")  ))
dim(isgs)
######### Beta specific ##############
genesbeta <- as.data.frame(  read.delim(  paste(dir,"genesbeta424", sep = "") )  )
dim(genesbeta)
