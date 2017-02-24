require(dtplyr)
require(data.table)
require(tibble)
require(dplyr)

load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/resultsCompressed.rda")
MDCKDETGFb <- data.table::as.data.table(resultsCompressed[["MDCK"]][["sleuth_results.gene"]]$conditionMDCKTGFb)
MDCKDETGFb[order(MDCKDETGFb$b, decreasing = T),]

random500Up <- dplyr::sample_n(MDCKDETGFb[MDCKDETGFb$b > 0], size = 500, replace = F)$target_id
random500Down <- dplyr::sample_n(MDCKDETGFb[MDCKDETGFb$b < 0], size = 500, replace = F)$target_id


l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[X4 %in% random500Down,c(start:end)]
    return(st)
  }
})
l1 <- l1[1:4]
names(l1) <- allGenes_normal.matrix.runDef$sample_labels

l2 <- lapply(names(l1), function(x){
  print(x)
  m1 <- data.frame(bin = c(1:300), 
                   value = apply(l1[[x]], 2, mean), 
                   sample = x,
                   group = paste(unlist(strsplit(x, "-"))[1:2], collapse = "-"))
  return(m1)
})
m1 <- as.data.frame(do.call("rbind", l2))
