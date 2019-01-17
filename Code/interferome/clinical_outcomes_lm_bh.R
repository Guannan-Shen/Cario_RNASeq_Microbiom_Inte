






############ Volcano Plots #####
## tempora
# EnhancedVolcano(genesbeta.de,
#                 lab = rownames(genesbeta.de),
#                 x = "Log2FoldChange",
#                 y = "padj",
#                 pCutoff = 5e-2,
#                 FCcutoff = 1,
#                 pLabellingCutoff = 5e-2,
#                 ## select labels to show
#                 # selectLab = c("cg18587484","cg00803922", " cg19425295"),
#                 ## point and label size 
#                 transcriptPointSize = 2.0,
#                 transcriptLabSize = 3.5,
#                 xlab = bquote(~Log[2]~ "Fold Change"),
#                 ylab = bquote(~-Log[10]~adjusted~italic(P)),
#                 title = "IFN-beta Genes: HIV Infected vs Health Control",
#                 #Modify border and remove gridlines
#                 gridlines.major = FALSE,
#                 gridlines.minor = FALSE,
#                 border = "full",
#                 borderWidth = 1.0,
#                 borderColour = "black",
#                 # the transparence of the dots
#                 colAlpha = 0.8,
#                 xlim = c(-3, 3),
#                 ylim = c(0, -log10(10e-25)),
#                 # adjust the legend
#                 legend=c("NS","log2 Fold Change >= 1","adjusted p-value <= 0.05",
#                          "adjusted p-value <= 0.05 & log2 Fold Change >= 1"),
#                 legendPosition = "bottom",
#                 legendLabSize = 9,
#                 legendIconSize = 3,
#                 # connectors
#                 DrawConnectors = TRUE,
#                 # 
#                 widthConnectors = 0.3,
#                 # 
#                 colConnectors = "grey40",
#                 col = c("grey30", "forestgreen", "royalblue", "tomato")
# )
