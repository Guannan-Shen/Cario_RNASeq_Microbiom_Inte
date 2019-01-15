## plot counts
group <- c(1,2,3)

is.matrix(cnts.edger) & is.vector(group.plot)
gs_plotcounts <- function(group, counts, deresults, top, topn, twocutoffs, fc, lists){
## group and lists vector, top, fc: logical, fc numeric, counts, deresults can be matrix 
`%nin%` <- Negate(`%in%`)
if("Symbol" %nin% colnames(deresults))
  stop("The DE analysis results data.frame should contain Symbol")
if("FDR" %nin% colnames(deresults))
  stop("The differential analysis results data.frame should contain FDR(BH)")
if("logFC" %nin% colnames(deresults))
  stop("The differential analysis results data.frame should contain logFC")
if(!is.logical(top) && is.logical(twocutoffs))
  stop("TRUE or FALSE for top and TRUE or FALSE for twocutoffs")
if(!nrow(counts) != nrow(deresults))
  stop("Normalized counts data should have same number of genes with DE analysis results data.frame")
if(!rownames(counts) == rownames(deresults))
  stop("Gene sets in counts and DE results are not the same")
if(!top == twocutoffs)
  stop("The logical values should be different")
  if (top){ 
  ## the top genes mode
    plot(group.plot + runif(ncol(cnts.edger), -0.1, 0.1), cnts.edger, 
     xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n", 
     xlab = "Groups", ylab = ylab, main = main)
  } else if (twocutoffs){
     
   }
  else
  {
  ## the genelists mode, the genes by fold change and pvalue cutoffs
   }
  
}


!is.null(rownames(cnts.edger)) 
top <- TRUE
twocutoffs <- FALSE
!is.logical(top) && is.logical(twocutoffs)
