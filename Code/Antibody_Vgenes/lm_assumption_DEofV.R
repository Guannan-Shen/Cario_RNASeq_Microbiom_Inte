######### DE of V genes ##################
library(edgeR)
library(tidyverse)
library(readxl)
library(openxlsx)
############### rnaseq
setwd("~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataRaw/Vgenes/")
cnts_raw <- read_excel("HIV_Infected_vs_Uninfected_relaxed_mapping_geneCounts.xlsx")
head(cnts_raw)
############ v genes lists ###################
vgenes <- read_excel("/home/guanshim/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/Vgenes_list.xlsx")
head(vgenes)
######clinical and microbiome data #########
clin_micro <- read_excel("/home/guanshim/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/clin_microbiome.xlsx")
head(clin_micro)
sample_size <- length(clin_micro$ID )
sample_size
## test the ids between rna-seq and clin microbiome data
# sum( colnames(cnts_raw)[-c(1:3)] != clin_micro$pid )
colnames(cnts_raw)[-c(1:3)]
## get cnts for norm and DE 
## filter to guarantee cnts larger than 1 
cnts_fsym <- cnts_raw[rowSums(as.matrix(cnts_raw[, -c(1:3)]) )>=(3*ncol(cnts_raw[, 1:32])), ]
## processing of rnaseq data
######################## targeted v genes of y ########################
cnts <- cnts_fsym %>% filter(Symbol %in% vgenes$Vsymbols ) %>%
   dplyr::select(-c(Gene_ID, Length ) ) %>% tibble::column_to_rownames("Symbol")
cnts <- as.matrix(cnts)
dim(cnts)
######## edgeR TMM Normalization ###############
## edgeR object
group <- c(rep("Control", 13), rep("HIV", 19))
y <- DGEList(counts= cnts, group=group)
## normalization 
y <- calcNormFactors(y, method = "TMM")
y$samples
#  TMM counts 
# If you run the cpm function on a DGEList object which contains TMM normalisation factors 
# then you will get TMM normalised counts. 
# Here is a snippet of the source code for the cpm function
cnts.edger <- edgeR::cpm(y)
nrow(cnts.edger) == 45
## 

# write.xlsx(cnts.edger ,
#            "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/Data_generated/vgenes_tmm.xlsx",
#            sheetName="rlog counts_correlation test")
#################### adjusting for age gender DE analysis ################### 
# get the design matrix
# the group
gender <- clin_micro$Gender
group 
age <- clin_micro$Age
design <- model.matrix(~ age + gender + group)
design
# general glm 
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
##################### global norm then filtered DE #######################
# the same design, different y
cnts <- cnts_fsym %>% dplyr::select(-c(Symbol, Length ) ) %>% tibble::column_to_rownames("Gene_ID")
cnts <- as.matrix(cnts)
dim(cnts)
##### 
y <- DGEList(counts= cnts, group=group)
## normalization 
y <- calcNormFactors(y, method = "TMM")
y$samples

###
head(y$counts )
sum( rownames(y$counts) !=  cnts_fsym$Gene_ID )
norm_filtered <- data.frame(edgeR::cpm(y))
norm_filtered$Symbol <- cnts_fsym$Symbol
norm_v <- norm_filtered %>% filter(Symbol %in% vgenes$Vsymbols ) %>% tibble::column_to_rownames("Symbol")
dim(norm_v)
head(norm_v)
norm_v <- as.matrix(norm_v)

##### 
y <- DGEList(counts= norm_v, group=group)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)

########################### whole global ###############################
##################### global norm then filtered DE #######################
# the same design, different y
cnts <- cnts_fsym %>% dplyr::select(-c(Symbol, Length ) ) %>% tibble::column_to_rownames("Gene_ID")
cnts <- as.matrix(cnts)
dim(cnts)
##### 
y <- DGEList(counts= cnts, group=group)
## normalization 
y <- calcNormFactors(y, method = "TMM")
y$samples
## get the counts 
cnts.edger <- edgeR::cpm(y, normalized.lib.sizes = TRUE,
                         log = FALSE)
dim(cnts.edger)
sum(rownames(cnts.edger) != cnts_fsym$Gene_ID)
cnts.edger.vgenes <- data.frame(cnts.edger) %>% mutate(Symbol = cnts_fsym$Symbol) %>%
                     filter(Symbol %in% vgenes$Vsymbols) %>% select(Symbol, everything())
dim(cnts.edger.vgenes)
tail(cnts.edger.vgenes)

write.xlsx(cnts.edger.vgenes ,
           "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/Data_generated/vgenes_tmm.xlsx",
           sheetName="TMM counts, unlogged")

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
sum( rownames(lrt$table) != cnts_fsym$Gene_ID )
lrt_v <- data.frame(lrt$table) %>% mutate(Symbol = cnts_fsym$Symbol ,
                                          FDR = p.adjust(PValue, method = "BH" ) ) %>% 
                                     filter( Symbol %in% vgenes$Vsymbols ) %>% 
                                    select( Symbol, everything()) %>% arrange(FDR)

write.xlsx(lrt_v ,
          "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/Data_generated/DE_vgenes_hiv.xlsx",
           sheetName="Global_DE_age_gender_hiv_Vgenes")
########################### check the lm assumption of single mutation  ############################
library(readxl)
library(boot)
setwd("/home/guanshim/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/DataProcessed/Vgenes/")
miseq <- read_excel("miseq.xlsx")
head(miseq)
hist(miseq$iga_smu )
hist(miseq$igk_smu )
hist(miseq$igl_smu )
# diagnostic plots 
# Xs
Clin <- read_excel("clin_microbiome.xlsx")
########## iga
par(mfrow = c(2,2))
lm = lm(miseq$iga_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender )
plot(lm, main = "lm: iga smu")
par(mfrow = c(2,2))
glmGamma <- glm(miseq$iga_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender, 
                family = Gamma(link = "log") )
glm.diag.plots(glmGamma )

########## igk
par(mfrow = c(2,2))
lm = lm(miseq$igk_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender )
plot(lm , main = "lm: igk smu")
par(mfrow = c(2,2))
glmGamma <- glm(miseq$igk_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender, 
                family = Gamma(link = "log") )
glm.diag.plots(glmGamma )
############ igl
par(mfrow = c(2,2))
lm = lm(miseq$igl_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender )
plot(lm , main = "lm: igl smu")
par(mfrow = c(2,2))
glmGamma <- glm(miseq$igl_smu ~ Clin$`CD4 T cells (% viable)` + Clin$Age + Clin$Gender, 
                family = Gamma(link = "log") )
glm.diag.plots(glmGamma )






