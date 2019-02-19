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

########################## get cnts.edger deseq2 tpm ##################
################# deseq2 ###########################
# counts from EDASeq (DESeq2)
# pData is phenoData from Biobase
countData <- counts(set) #Matrix with transcripts IDs as rows and sample IDs as columns
colData <- pData(set) #Vector of type list in which the condition column is the treat/control identfier, and the rownames are sample IDs

#Run DESeq function using above objects

## now using deseq2
dds <- DESeqDataSetFromMatrix(countData = counts(set), colData = pData(set),design = ~ condition)
dds <- estimateSizeFactors(dds)
## normalization factors
sizeFactors(dds)
cnts.deseq2 <- counts(dds, normalized=TRUE)

######## edgeR TMM ###############
## edgeR object
group <- c(rep(1, 13), rep(2, 19))
y <- DGEList(counts= as.matrix(cnts_f), group=group)
## normalization 
y <- calcNormFactors(y, method = "TMM")
cnts.edger <- edgeR::cpm(y)
################### DE analysis by edger ################
y <- estimateDisp(y)
et <- exactTest(y, pair = 1:2)
res.edger <- et$table %>% dplyr::mutate(Symbol = cnts_fsym$Symbol,
                                        FDR = p.adjust(et$table$PValue, method = "BH") )
rownames(res.edger) <- rownames(et$table)

########## TPM #########
cnts.tpm <- calculateTPM(cnts_f, effective_length = cnts_fsym$Length)

######### setup for scater plots #########
## cnts. list
cnts.list <- list(cnts_f, cnts.edger, cnts.deseq2, cnts.tpm)
cnts.names <- c("Filtered Counts", "TMM Normalized Counts",
                "DESeq2 Normalized Counts", "TPM Normalized Counts")
rna_info <- data.frame(pid = rna.pid, Group = as.factor(c(rep("Control", 13), 
                                                          rep("HIV Infected", 19) )) )
rna_info$Group %<>% relevel("Control")
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
  print( scater::plotPCA(example_sce, colour_by = "Group", shape_by = "Group") + labs(caption = paste(  "PCA Plot:",cnts.names[i]) ))
}


############# read in all clinical data
clinical_ready <- as.data.frame( read_excel("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/clinical_ready.xlsx") )
clinical_names <- c("Blood CD4 T Cell Counts (cells/ul)", "Plasma Viral Load", "Tissue HIV RNA (per CD4 T cell)",
                    "Tissue CD4 T Cell Counts (number/g)", "IL-6 (pg/ml)", "CRP (ug/ml)", "iFABP (pg/ml)",
                    "sCD27 (U/ml)", "CD14 (ng/ml)", "LPS (pg/ml)", "LTA (OD)", 
                    base::paste("IFN", '\u03b1', sep = "" ),  base::paste("IFN", '\u03b2', sep = "" ),
                    ##### Adding this clinical parameter ###########
                    "CD4 T cells (% viable CD45+ cells)")
dim(clinical_ready)
n_clinical <- length(clinical_names)
clinical_sum <- matrix(NA, 14, 8)

########### read in altered gene edger data and result ######
genesbeta.edger.05 <- read.csv("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/genesbeta.edger.05.csv")
isgs.edger.05 <- read.csv("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/isgs.edger.05.csv")
cnts.isgs.edger.05 <- cnts.edger[ rownames(cnts.edger) %in% isgs.edger.05$Gene_ID, ]
dim(cnts.isgs.edger.05)
write.csv(cnts.isgs.edger.05,
          "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/cnts.isgs.edger.05.csv")
cnts.genesbeta.edger.05 <- cnts.edger[ rownames(cnts.edger) %in% genesbeta.edger.05$Gene_ID, ]
dim(cnts.genesbeta.edger.05)
write.csv(cnts.genesbeta.edger.05,
          "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/cnts.genesbeta.edger.05.csv")
##### vocalno plots for all gene list DE compare ###

############# read in clinical_order and gene rlog data
clinical_order <- read.csv("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/clinical_order.csv")
isgs.05.lin <- read.csv("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/isgs.05.lin.csv")
genesbeta.05.lin <- read.csv("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/corr/genesbeta.05.lin.csv")
clinical_order <- clinical_order[,-1]
isgs.05.lin <- isgs.05.lin[,-1]
genesbeta.05.lin <- genesbeta.05.lin[,-1]

########################## linear regression ##############################
## equal length of outcomes and covariates
gene_IFNReg <- function(gene_matrix, clinical_variable, clin_var_name){
  # get names ready
  genelistname = base::colnames(gene_matrix)
  ## number of gene to test, also the number of multiple test
  n_gene = ncol(gene_matrix)
  ## outcome lm
  outcome_lm = lapply(1:n_gene, function(i){
    lm = lm(gene_matrix[,i]~ clinical_variable + clinical_order$age + clinical_order$sex )
    coef = summary(lm)$coefficients[2, ]
    return(coef)
  })
  outcome_lm = data.frame(matrix(unlist(outcome_lm), ncol = 4, byrow = TRUE,
                                 dimnames = list(
                                   c(colnames(gene_matrix)),
                                   c("Estimate", "Std.Error", "t.statistic", "p.value"))))
  
  # adjusted p-value
  outcome_lm =  outcome_lm %>% 
    dplyr::mutate(FDR = p.adjust(p.value, "BH", n_gene ),
                  names = colnames(gene_matrix)) %>% 
    dplyr::mutate(Estimate = round(Estimate, 10),
                  Std.Error = round(Std.Error, 10),
                  t.statistic = round(t.statistic,4)
    )%>% 
    select(names, everything())
  # sort by p.value
  outcome_lm = outcome_lm[order(outcome_lm$p.value), ]
  
  ## sample size 
  size = sum(!is.na(clinical_variable))
  
  ## summary table 
  return(list(results = data.frame(outcome_lm), size = size, clinical = clin_var_name ))
}

################ get the sig altered gene list #########################


######## finally ############
for(i in 1:n_clinical) {
  ## number of clinical virable
  j = c(6:19)[i]
  cliname = base::colnames(clinical_order)[j]
  ## linear regression
  lin_res_isgs = gene_IFNReg(isgs.05.lin, clinical_order[,j], clinical_names[i])
  lin_res_genesbeta = gene_IFNReg(genesbeta.05.lin, clinical_order[,j], clinical_names[i])
  ## save data
  
  write.xlsx(lin_res_isgs$results,
             paste("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/isgs_", cliname,".xlsx", sep = ""),
             sheetName= paste("ISGs_", cliname, sep = "") )
  write.xlsx(lin_res_genesbeta$results,
             paste("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/ifnbeta_", cliname,".xlsx", sep = ""),
             sheetName= paste("ifnbeta_", cliname, sep = ""))
  # }
  ## check basic
  if(lin_res_isgs$size == lin_res_genesbeta$size){
    print("Good")
  }else stop("sample size wrong")
  if(lin_res_isgs$clinical == lin_res_genesbeta$clinical){
    print("Good")
  }else stop("clinical parameter wrong")
  ## number of genes
  nisgs = nrow(lin_res_isgs$results)
  ngenesbeta = nrow(lin_res_genesbeta$results)
  ################## summary table ###################
  ###  sig prop
  isgs.sig.no =  sum(lin_res_isgs$results$FDR <= 0.05)
  genesbeta.sig.no =  sum(lin_res_genesbeta$results$FDR <= 0.05)
  if( (isgs.sig.no >= 5) & (genesbeta.sig.no >= 5) ){
    prop.test.sig = prop.test(x = c(isgs.sig.no, genesbeta.sig.no), n = c(nisgs, ngenesbeta), correct = FALSE)
    sig.p = prop.test.sig$p.value
  }else{
    prop.test.sig = prop.test(x = c(isgs.sig.no, genesbeta.sig.no), n = c(nisgs, ngenesbeta), correct = TRUE)
    sig.p = prop.test.sig$p.value
  }
  ########### volcano plots ################
  if(min(lin_res_isgs$results$FDR) <= 0.05){
    rownames(lin_res_isgs$results) <- lin_res_isgs$results$names
    p1 <- EnhancedVolcano(lin_res_isgs$results ,
                          lab = rownames(lin_res_isgs$results),
                          x = "Estimate",
                          y = "p.value",
                          pCutoff = ifelse( min(lin_res_isgs$results$FDR) <= 0.05,
                                            lin_res_isgs$results$p.value[max(which(lin_res_isgs$results$FDR <= 0.05))] + 0.005, 
                                            0.05),
                          FCcutoff = round(max(abs(quantile(lin_res_isgs$results$Estimate, c(0.25, 0.75) ))),1),
                          title = paste("Core ISGs: association with ", lin_res_genesbeta$clinical , sep = ""),
                          xlim = range(lin_res_isgs$results$Estimate)*1.3,
                          ylim = c(0, -log10(min(lin_res_isgs$results$p.value)/5000)),
                          legend=c("NS","Slope","FDR Significant",
                                   "FDR Significant & Slope"),
                          ## select labels to show
                          # selectLab = c("cg18587484","cg00803922", " cg19425295"),
                          ## point and label size 
                          transcriptPointSize = 2,
                          transcriptLabSize = 3,
                          xlab = bquote(~Log[2]~ "Slope"),
                          ylab = bquote(~-Log[10]~italic(P)),
                          
                          #Modify border and remove gridlines
                          gridlines.major = FALSE,
                          gridlines.minor = FALSE,
                          border = "full",
                          borderWidth = 1.0,
                          borderColour = "black",
                          # the transparence of the dots
                          colAlpha = 0.8,
                          
                          # adjust the legend
                          legendPosition = "bottom",
                          legendLabSize = 9,
                          legendIconSize = 3,
                          # connectors
                          DrawConnectors = TRUE,
                          # 
                          widthConnectors = 0.2,
                          # 
                          colConnectors = "grey40",
                          col = c("grey30", "forestgreen", "royalblue", "tomato")
    )
    ggsave(paste("Volcanol_Plot_ISGs_", cliname, sep = ""), width = 6, height = 9, device = "tiff", path = "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/Reports/plots", dpi = 300)
  } else{print("No Sig.")}
  
  if(min(lin_res_genesbeta$results$FDR) <= 0.05){
    rownames(lin_res_genesbeta$results) <- lin_res_genesbeta$results$names
    p1 <- EnhancedVolcano(lin_res_genesbeta$results ,
                          lab = rownames(lin_res_genesbeta$results),
                          x = "Estimate",
                          y = "p.value",
                          pCutoff = ifelse( min(lin_res_genesbeta$results$FDR) <= 0.05,
                                            lin_res_genesbeta$results$p.value[max(which(lin_res_genesbeta$results$FDR <= 0.05))] + 0.005, 
                                            0.005),
                          FCcutoff = max(abs(quantile(lin_res_genesbeta$results$Estimate, c(0.25, 0.75) ))) ,
                          title = paste("IFNBeta Genes: association with ", lin_res_genesbeta$clinical , sep = ""),
                          xlim = range(lin_res_genesbeta$results$Estimate)*1.3,
                          ylim = c(0, -log10(min(lin_res_genesbeta$results$p.value)/5000)),
                          legend=c("NS","Slope","FDR Significant",
                                   "FDR Significant & Slope"),
                          ## select labels to show
                          # selectLab = c("cg18587484","cg00803922", " cg19425295"),
                          ## point and label size 
                          transcriptPointSize = 2,
                          transcriptLabSize = 3,
                          xlab = bquote("Slope"),
                          ylab = bquote(~-Log[10]~italic(P)),
                          
                          #Modify border and remove gridlines
                          gridlines.major = FALSE,
                          gridlines.minor = FALSE,
                          border = "full",
                          borderWidth = 1.0,
                          borderColour = "black",
                          # the transparence of the dots
                          colAlpha = 0.8,
                          
                          # adjust the legend
                          legendPosition = "bottom",
                          legendLabSize = 9,
                          legendIconSize = 3,
                          # connectors
                          DrawConnectors = TRUE,
                          # 
                          widthConnectors = 0.2,
                          # 
                          colConnectors = "grey40",
                          col = c("grey30", "forestgreen", "royalblue", "tomato")
    )
    ggsave(paste("Volcanol_Plot_IFNBeta_", cliname, sep = ""), width = 6, height = 9, device = "tiff", path = "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/Reports/plots", dpi = 300)
  } else{print("No Sig.")}
  
  
  
  ## positive prop
  isgs.pos.no = sum(lin_res_isgs$results$Estimate > 0)
  genesbeta.pos.no = sum(lin_res_genesbeta$results$Estimate > 0)
  if( (isgs.pos.no >= 5) & (genesbeta.pos.no >= 5) ){
    prop.test.pos = prop.test(x = c(isgs.pos.no, genesbeta.pos.no), n = c(nisgs, ngenesbeta), correct = FALSE)
    pos.p = prop.test.pos$p.value
  }else{
    prop.test.pos = prop.test(x = c(isgs.pos.no, genesbeta.pos.no), n = c(nisgs, ngenesbeta), correct = TRUE)
    pos.p = prop.test.pos$p.value
  }
  ## If you check prop.test, chisq.test and z-test on your data then they all give you the same p-value
  ## continuity correction if anyone less than 5
  ## the clinical summary table 
  clinical_sum[i, ] <- c(lin_res_genesbeta$clinical, lin_res_genesbeta$size,
                         isgs.sig.no/nisgs, genesbeta.sig.no/ngenesbeta , sig.p, 
                         isgs.pos.no/nisgs,
                         genesbeta.pos.no/ngenesbeta , pos.p)
  
  colnames(clinical_sum) <- c("Clinical Parameter", "Sample Size", "ISGs: Sig. Prop.",
                              "IFNbeta Genes: Sig.", "p value",
                              "ISGs: Positive Corr. Prop.",
                              "IFNbeta Genes: Pos.", "p value (Pos)")
}


## clinical_sum
write.xlsx(clinical_sum, "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/clinical_summary.xlsx")
clinical_sum[,c(3,4,6,7)] <- round(as.numeric( as.character( clinical_sum[,c(3,4,6,7) ]) ), 4)
clinical_sum[,c(5,8)] <- format(as.numeric( as.character( clinical_sum[,c(5,8) ]) ), digits = 3, scientific = TRUE)
clinical_sum[12:13,1] <- c("IFNalpha", "IFNbeta")
kable(data.frame(clinical_sum[,c(1:5)]), digits = c(2,2,4,4,40))
kable(data.frame(clinical_sum[,c(1,2,6:8)]), digits = c(2,2,4,4,40))

## filtered clinical_sum by non-zero
clinical_sum_sig <- clinical_sum[ ,]


