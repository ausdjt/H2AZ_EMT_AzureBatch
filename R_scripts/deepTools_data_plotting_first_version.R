
#****************************************
# load the deepTools computeMatrix output
allGenes_normal.matrix <- readr::read_tsv(gzfile("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/deepTools/computeMatrix/reference-point/duplicates_marked/TSS/allSamples_allGenes_normal.matrix.gz"),
                                          comment = "@",
                                          col_names = FALSE)

# run definition is encoded as JSON in the first line of file
allGenes_normal.matrix.runDef <- readLines(gzfile("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/deepTools/computeMatrix/reference-point/duplicates_marked/TSS/allSamples_allGenes_normal.matrix.gz"),
                                           n = 1)
allGenes_normal.matrix.runDef <- jsonlite::fromJSON(gsub("@", "", allGenes_normal.matrix.runDef))

gr.allGenes_normal <- GenomicRanges::GRanges(seqnames = allGenes_normal.matrix$X1,
                                             IRanges(start = allGenes_normal.matrix$X2,
                                                     end = allGenes_normal.matrix$X3,
                                                     names = allGenes_normal.matrix$X4),
                                             strand = allGenes_normal.matrix$X6)
allGenes_normal.matrix <- allGenes_normal.matrix[,-c(1:3,5:6)]

l1 <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x != 1200){
    st <- allGenes_normal.matrix[,c(start:end)]
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

ggplot(m1, aes(bin, value)) + geom_smooth(method = "lm", formula = y ~ splines::bs(x, 200), se = TRUE) + facet_wrap(~group)
ggplot(m1, aes(bin, value)) + geom_line() + facet_wrap(~group)

# EMT genes
EMTup <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_up_genes.bed")
EMTdown <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_down_genes.bed")

EMTup_genes.matrix <- allGenes_normal.matrix[which(allGenes_normal.matrix$X4 %in% EMTup$name),]
EMTdown_genes.matrix <- allGenes_normal.matrix[which(allGenes_normal.matrix$X4 %in% EMTdown$name),]

# EMT up genes
EMTup_genes <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x < max(allGenes_normal.matrix.runDef$sample_boundaries)){
    st <- EMTup_genes.matrix[,c(start:end)]
    return(st)
  }
})
EMTup_genes <- EMTup_genes[1:length(allGenes_normal.matrix.runDef$sample_labels)]
names(EMTup_genes) <- allGenes_normal.matrix.runDef$sample_labels
EMTup_genes <- lapply(names(EMTup_genes), function(x){
  print(x)
  m1 <- data.frame(bin = c(1:300), 
                   value = apply(l1[[x]], 2, mean), 
                   sample = x,
                   group = paste(unlist(strsplit(x, "-"))[1:2], collapse = "-"),
                   geneset = "EMTup")
  return(m1)
})
EMTup_genes <- as.data.frame(do.call("rbind", EMTup_genes))
checkPlot <- function(x){
  ggplot(x, aes(bin, value)) + geom_line() + facet_wrap(c("sample", "geneset"), ncol = 2) + ylim(-20,80)
}
checkPlot(EMTup_genes)

# EMT down genes
EMTdown_genes <- lapply(allGenes_normal.matrix.runDef$sample_boundaries, function(x){
  start <- x + 2
  end <- start + 299
  if (x < max(allGenes_normal.matrix.runDef$sample_boundaries)){
    st <- EMTdown_genes.matrix[,c(start:end)]
    return(st)
  }
})
EMTdown_genes <- EMTdown_genes[1:length(allGenes_normal.matrix.runDef$sample_labels)]
names(EMTdown_genes) <- allGenes_normal.matrix.runDef$sample_labels
EMTdown_genes <- lapply(names(EMTdown_genes), function(x){
  print(x)
  m1 <- data.frame(bin = c(1:300), 
                   value = apply(l1[[x]], 2, mean), 
                   sample = x,
                   group = paste(unlist(strsplit(x, "-"))[1:2], collapse = "-"),
                   geneset = "EMTdown")
  return(m1)
})
EMTdown_genes <- as.data.frame(do.call("rbind", EMTdown_genes))
checkPlot(EMTdown_genes)

EMTgenes <- rbind(EMTup_genes, EMTdown_genes)
checkPlot(EMTgenes)

EMTgenes <- EMTgenes[grep("H2AZ", EMTgenes$sample),]
EMTgenes$sample <- as.character(EMTgenes$sample)
EMTgenes$group <- as.character(EMTgenes$group)
d1 <- EMTgenes[EMTgenes$sample == "H2AZ-TGFb_normal_RPKM", "value"] - EMTgenes[EMTgenes$sample == "H2AZ-WT_normal_RPKM", "value"]
d1 <- data.frame(bin = rep(1:300, 2), 
                 value = d1, 
                 sample = "diff", 
                 group = "diff", 
                 geneset = c(rep("EMTup", table(EMTgenes$geneset)["EMTup"]), 
                             rep("EMTdown", table(EMTgenes$geneset)["EMTdown"])))
EMTgenes <- rbind(EMTgenes, d1)
EMTgenes$sample <- as.factor(EMTgenes$sample)
EMTgenes$group <- as.factor(EMTgenes$group)
pdf("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/TSS_H2AZ_coverage_Tan_et_al_EMT_w_diff.pdf")
ggplot(EMTgenes, aes(bin, value)) + geom_line() + facet_wrap(c("sample", "geneset"), ncol = 2)
dev.off()

