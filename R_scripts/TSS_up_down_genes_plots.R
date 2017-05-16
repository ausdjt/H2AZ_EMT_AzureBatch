require(dtplyr)
require(data.table)
require(tibble)
require(dplyr)

load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/resultsCompressed.rda")
MDCKDETGFb <- data.table::as.data.table(resultsCompressed[["MDCK"]][["sleuth_results.gene"]]$conditionMDCKTGFb)
MDCKDETGFb[order(MDCKDETGFb$b, decreasing = T),]

random100Up <- dplyr::sample_n(MDCKDETGFb[MDCKDETGFb$b > 0], size = 100, replace = F)$target_id
random100Down <- dplyr::sample_n(MDCKDETGFb[MDCKDETGFb$b < 0], size = 100, replace = F)$target_id


# random 100 down ---------------------------------------------------------
dataDir <- "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/deepTools/computeMatrix/reference-point/duplicates_marked/TSS"
f1 <- list.files(dataDir, pattern = "random100down_normal", full.names = T)
random100Down <- readr::read_tsv(gzfile(f1),
                                 comment = "@",
                                 col_names = FALSE)

l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[allGenes_normal.matrix$X4 %in% random100Down,c(start:end)]
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
                   group = "random100Down")
  return(m1)
})
m1 <- as.data.frame(do.call("rbind", l2))


# random 100 up -----------------------------------------------------------
l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[allGenes_normal.matrix$X4 %in% random100Up,c(start:end)]
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
                   group = "random100Up")
  return(m1)
})
m2 <- as.data.frame(do.call("rbind", l2))


# merge for plotting ------------------------------------------------------
m1 <- rbind(m1, m2)

pdf("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/TSS_coverage_random.pdf")
ggplot(m1, aes(bin, value)) + geom_line() + facet_wrap(c("sample", "group"), ncol = 2)
dev.off()

# only H2AZ signal, with difference map
m1 <- m1[grep("H2AZ", m1$sample),]
m1$sample <- as.character(m1$sample)
m1$group <- as.character(m1$group)
d1 <- m1[m1$sample == "H2AZ-TGFb_normal_RPKM", "value"] - m1[m1$sample == "H2AZ-WT_normal_RPKM", "value"]
d1 <- data.frame(bin = rep(1:300, 2), value = d1, sample = "diff", group = c(rep("random100Down", 300), rep("random100Up", 300)))
m1 <- rbind(m1, d1)
m1$sample <- as.factor(m1$sample)
m1$group <- as.factor(m1$group)
pdf("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/TSS_H2AZ_coverage_random_w_diff.pdf")
ggplot(m1, aes(bin, value)) + geom_line() + facet_wrap(c("sample", "group"), ncol = 2)
dev.off()

# EMT genes plotting ------------------------------------------------------

# EMT genes
EMTup <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_up_genes.bed")
EMTdown <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_down_genes.bed")

l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[allGenes_normal.matrix$X4 %in% EMTup$name,c(start:end)]
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
                   group = "EMTup")
  return(m1)
})
EMTupData <- as.data.frame(do.call("rbind", l2))

l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[allGenes_normal.matrix$X4 %in% EMTdown$name,c(start:end)]
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
                   group = "EMTdown")
  return(m1)
})
EMTdownData <- as.data.frame(do.call("rbind", l2))

EMTData <- rbind(EMTupData, EMTdownData)

pdf("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/TSS_coverage_DE_EMT_genes.pdf")
ggplot(EMTData, aes(bin, value)) + geom_line() + facet_wrap(c("sample", "group"), ncol = 2)
dev.off()


# look at nucleosome-associated genes -------------------------------------
nucleosomeAssociatedGenes <- read.csv("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/GO_nucleosome_associated_genes.csv", header = T, row.names = 2)
nucleosomeAssociatedGenes <- as.data.table(nucleosomeAssociatedGenes)
setkey(nucleosomeAssociatedGenes, "Gene.ID")
t1 <- MDCKDETGFb[MDCKDETGFb$target_id %in% nucleosomeAssociatedGenes$Gene.ID]
setkey(t1, "target_id")
plot(t1$b, -log10(t1$pval))
t1 <- t1[order(t1$b, decreasing = T),]
t1 <- merge(t1, nucleosomeAssociatedGenes[, c("Gene.ID", "Associated.Gene.Name", "Description")], by.x = "target_id", by.y = "Gene.ID", all.x = T, all.y = F)
t1 <- t1[!duplicated(t1$target_id)]
t1[t1$b > 0 & t1$qval < 0.01]
write_csv(t1, "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/GO_nucleosome_associated_DE_genes.csv")
