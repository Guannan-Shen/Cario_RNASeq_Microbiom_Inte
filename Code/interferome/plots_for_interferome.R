##### validation plots ##########
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
  ggsave(paste("RLE_Plot_" ,cnts.names[i], sep = ""), width = 4, height = 4, device = "tiff", path = "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/Reports/plots", dpi = 300)
  ## pca plots
  example_sce <- SingleCellExperiment(
    assays = list(logcounts = j),
    colData = rna_info)
  print( scater::plotPCA(example_sce, colour_by = "Group", shape_by = "Group") + labs(caption = paste(  "PCA Plot:",cnts.names[i]) ))
  ggsave(paste("PCA_Plot_" ,cnts.names[i], sep = ""), width = 4, height = 4,device = "tiff", path = "~/Documents/gitlab/Cario_RNASeq_Microbiom_Inte/Reports/plots", dpi = 300)
}
