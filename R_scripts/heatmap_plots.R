require(rtracklayer)

setwd("~/Data/Tremethick/EMT/GenomeWide/")

source("~/Development/GeneralPurpose/R/heatmap.3.R")

# boxplot of all genes ----------------------------------------------------
# load RPKM normalised data from bigwig files
load("~/Data/Tremethick/EMT/GenomeWide/Cfam3.genes.rda")
gr.which <- promoters(Cfam3.genes, upstream = 1000, downstream = 1000)

bw.dir <- c("~/Data/Tremethick/EMT/GenomeWide/bigwig/")
bw.TGFb_H2AZ <- import(paste(bw.dir, "H2AZ_TGFb_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.TGFb_Input <- import(paste(bw.dir, "Input_TGFb_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.WT_H2AZ <- import(paste(bw.dir, "H2AZ_WT_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.WT_Input <- import(paste(bw.dir, "Input_WT_10bp_RPKM.bw", sep = ""), which = gr.which)

bw.TGFb_H2AZ_bkgrd_sub <- import(paste(bw.dir, "H2AZ_TGFb_Input_subtracted.bw", sep = ""), which = gr.which)
bw.WT_H2AZ_bkgrd_sub <- import(paste(bw.dir, "H2AZ_WT_Input_subtracted.bw", sep = ""), which = gr.which)

gl1 <- promoters(Cfam3.genes, upstream = 1000, downstream = 1000)
l1 <- lapply(gl1, function(x){
  v1 <- c(h2az_tgfb = mean(subsetByOverlaps(bw.TGFb_H2AZ, x)$score))
  v2 <- c(h2az_wt = mean(subsetByOverlaps(bw.WT_H2AZ, x)$score))
  v3 <- mean(subsetByOverlaps(bw.TGFb_H2AZ_bkgrd_sub, x)$score)
  v4 <- mean(subsetByOverlaps(bw.WT_H2AZ_bkgrd_sub, x)$score)
  v5 <- c(v1, v2, v3, v4)
})
df1 <- data.frame(matrix(unlist(l1), nrow = length(l1), byrow=T))
rownames(df1) <- names(l1)
colnames(df1) <- c("h2az_tgfb", "h2az_wt", "h2az_tgfb_bkgr_sub", "h2az_wt_bkgr_sub")

# heatmap/boxplot of qPCR data --------------------------------------------
# load RPKM normalised data from bigwig files
gr.which <- promoters(Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)], upstream = 1000, downstream = 1000)
bw.dir <- c("~/Data/Tremethick/EMT/GenomeWide/bigwig/")
bw.TGFb_H2AZ <- import(paste(bw.dir, "H2AZ_TGFb_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.TGFb_Input <- import(paste(bw.dir, "Input_TGFb_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.WT_H2AZ <- import(paste(bw.dir, "H2AZ_WT_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.WT_Input <- import(paste(bw.dir, "Input_WT_10bp_RPKM.bw", sep = ""), which = gr.which)
bw.TGFb_H2AZ_bkgrd_sub <- import(paste(bw.dir, "H2AZ_TGFb_Input_subtracted.bw", sep = ""), which = gr.which)
bw.WT_H2AZ_bkgrd_sub <- import(paste(bw.dir, "H2AZ_WT_Input_subtracted.bw", sep = ""), which = gr.which)

gl1 <- promoters(Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)], upstream = 1000, downstream = 1000)

l1 <- lapply(gl1, function(x){
  v1 <- c(h2az_tgfb = mean(subsetByOverlaps(bw.TGFb_H2AZ, x)$score))
  v2 <- c(h2az_wt = mean(subsetByOverlaps(bw.WT_H2AZ, x)$score))
  v3 <- mean(subsetByOverlaps(bw.TGFb_H2AZ_bkgrd_sub, x)$score)
  v4 <- mean(subsetByOverlaps(bw.WT_H2AZ_bkgrd_sub, x)$score)
  v5 <- c(v1, v2, v3, v4)
})
df1 <- data.frame(matrix(unlist(l1), nrow = length(l1), byrow=T))
rownames(df1) <- names(l1)
colnames(df1) <- c("h2az_tgfb", "h2az_wt", "h2az_tgfb_bkgr_sub", "h2az_wt_bkgr_sub")
# df1[which(df1$h2az_tgfb_bkgr_sub < 0), ]$h2az_tgfb_bkgr_sub <- 0
# df1[which(df1$h2az_wt_bkgr_sub < 0), ]$h2az_wt_bkgr_sub <- 0
df1$hgnc_symbol <- Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)]$hgnc_symbol
heatmap2 <- heatmap.3(as.matrix((df1[, c(4,3)])), 
                      trace = "none", 
                      cexCol = 0.6,
                      cexRow = 0.5,
                      main = "qPCR genes\nH2A.Z coverage TSS +/-1Kb [mean RPKM]", 
                      hclustfun=function(x) hclust(x,method="ward.D"),
                      labRow = df1$hgnc_symbol,
                      col = redgreen(15),
                      symm = F)
                      #labCol = c("Epithelial", "Mesenchymal"))

# create a histogram of the data (here log2 transformed)
df2 <- as((df1[, c("h2az_wt")]), "matrix")
df2 <- rbind(df2, as((df1[, c("h2az_tgfb")]), "matrix"))
df2 <- data.frame(df2)
df2$var <- c(rep(paste("MDCK - Untreated [N = ", length(gl1), "]", sep = ""), length(gl1)), rep(paste("MDCK - TGFb-treated [N = ", length(gl1), "]", sep = ""), length(gl1)))
df2$var <- as.factor(df2$var)
histo1 <- ggplot(df2,aes(x=df2, group=var))
#pdf("Histogram_H2AZ_nucleosome_500TSS0_all_genes_FDR0.01.pdf", height = 10, width = 10)
histo1 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
dev.off()
pdf("Boxplot_H2AZ_nucleosome_500TSS0_aPCR_genes_FDR0.01.pdf", height = 10, width = 10)
histo1 + geom_boxplot(position = "identity", aes(y = df2, x= var)) + labs(title = "qPCR genes\nH2A.Z coverage TSS +/-1Kb", x = "Sample", y = "H2A.Z ChIP [mean RPKM]") 
dev.off()
