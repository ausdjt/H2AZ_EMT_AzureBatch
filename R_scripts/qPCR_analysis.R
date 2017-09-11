#qPCR data analysis
library(HTqPCR)
library(ggplot2)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
require(deepToolsUtils)

# modified function for limma moderated t-test on Ct data
source("~/Development/GeneralPurpose/R/limmaCtDataMod.R")

# set working directory
setwd("~/Data/Tremethick/EMT/qPCR_data_analysis")


# set up biomaRt connection -----------------------------------------------
dog <- useEnsembl(biomart = "ensembl", dataset = "cfamiliaris_gene_ensembl", host = "asia.ensembl.org")
attribs <- listAttributes(dog)


load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/R_Analysis/cfam.qPCRGenesTab.rda")
qPCRGeneList <- readLines("~/Data/Tremethick/EMT/MDCK qPCR data/genelist.txt")
qPCRGenesTab <- cfam.qPCRGenesTab
qPCRGeneList[which(!qPCRGeneList %in% qPCRGenesTab$hgnc_symbol)]
# qPCRGeneList.missing <- c("KRT7" = "ENSCAFG00000007307",
#                           "LOC488207" = NULL,
#                           "OCLN" = "ENSCAFG00000007805",
#                           "SIP1" = "ENSCAFG00000013859",
#                           "TCF4" = "ENSCAFG00000000140",
#                           "TGFB1" = "ENSCAFG00000005014",
#                           "TMEFF1" = "ENSCAFG00000002577",
#                           "TWIST1" = "ENSCAFG00000012469", #using TWIST2 - TWIST does not seem to exist in dog genome
#                           "LOC478215/H2AZ" = "ENSCAFG00000010615",
#                           "HPRT1" = "ENSCAFG00000018870",
#                           "LDHAL6B" = "ENSCAFG00000009211")
qPCRGeneList.missing <- c("SIP1" = "ENSCAFG00000013859", "LOC488207" = NULL)
t1 <- cbind(rownames(as.data.frame(qPCRGeneList.missing)), qPCRGeneList.missing)
colnames(t1) <- c("hgnc_symbol", "ensembl_gene_id")
t2 <- rbind(qPCRGenesTab[,c("ensembl_gene_id", "hgnc_symbol")], 
            t1[,c("ensembl_gene_id", "hgnc_symbol")])
rownames(t2) <- t2$hgnc_symbol
qPCRGenesTab <- t2
rm(list = c("t1", "t2"))
qPCRGenesTab$type <- "GOI"
qPCRGenes.control <- c("B2M", "GUSB", "HPRT1", "GAPDH", "LDHAL6B")
qPCRGenesTab[qPCRGenes.control,]$type <- "CONTROL"
save(qPCRGenesTab, file = "qPCRGenesTab.rda")

qPCRGenesPositions <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), filters = "ensembl_gene_id", values = qPCRGenesTab$ensembl_gene_id, dog)
rownames(qPCRGenesPositions) <- qPCRGenesPositions$ensembl_gene_id
gr.qPCRGenesPositions <- GRanges(seqnames = qPCRGenesPositions$chromosome_name, 
                                 IRanges(qPCRGenesPositions$start_position, qPCRGenesPositions$end_position, names = qPCRGenesPositions$ensembl_gene_id),
                                 strand = c("+", "-")[match(qPCRGenesPositions$strand, c(1, -1))],
                                 hgnc_symbol = qPCRGenesPositions$hgnc_symbol)


# Analysis of TGFb-treatment experiment -----------------------------------
qPCRData.TGFb_exp1 <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/TGFb-treatment/20120620 PCR - TGFb-treatment Experiment 1.csv", header = T, row.names = 1, as.is = T)
qPCRData.TGFb_exp2 <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/TGFb-treatment/20120625 PCR - TGFb-treatment Experiment 2.csv", header = T, row.names = 1, as.is = T)
qPCRData.TGFb <- cbind(qPCRData.TGFb_exp1[, c(1:2)], qPCRData.TGFb_exp2[,c(1:2)], qPCRData.TGFb_exp1[, c(3:4)], qPCRData.TGFb_exp2[,c(3:4)])
colnames(qPCRData.TGFb) <- c(paste("Control_", c(1:4), sep= ""), paste("TGFb_treated_", c(1:4), sep= ""))

lapply(seq_along(1:ncol(qPCRData.TGFb)), function(x) {
  Ct <- qPCRData.TGFb[,x]
  genes <- rownames(qPCRData.TGFb)
  pos <- c(1:nrow(qPCRData.TGFb))
  flag <- rep("OK", nrow(qPCRData.TGFb))
  type <- c(rep("GOI", 84), rep("HKG", 5))
  fn <- paste("~/Data/Tremethick/EMT/MDCK qPCR data/TGFb-treatment/", colnames(qPCRData.TGFb)[x], ".csv", sep = "")
  print(fn)
  tab <- data.frame(cbind(genes, Ct, pos, flag, type ))
  write.table(tab, file = fn, row.names = F, col.names = F, sep = "\t")
})

files <- data.frame(ID = 1:ncol(qPCRData.TGFb))
files$Files <- unlist(lapply(seq_along(1:ncol(qPCRData.TGFb)), function(x) { paste(colnames(qPCRData.TGFb)[x], ".csv", sep = "") }))
files$Treatment <- c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated")

column.info <- list(flag = "flag", Ct = "Ct", feature = "genes", position = "pos", type = "type")

qPCRdata.raw <- readCtData(files$Files, 
                           format = "plain", 
                           path = "~/Data/Tremethick/EMT/MDCK qPCR data/TGFb-treatment/", 
                           n.features = 89, 
                           column.info = list(flag = 4, Ct = 2, feature = 1, position = 3, type = 5), 
                           n.data = 1,
                           header = F)

pData(qPCRdata.raw) <- data.frame(Sample = c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"), Replicate = rep(1:4,2 ))

g <- featureNames(qPCRdata.raw)[1:10]
g <- c("TGFB1", "MMP9", "CDH1", "TGFB2", "TGFB3")
plotCtOverview(qPCRdata.raw, calibrator = "Control", genes = g, conf.int = T, ylim = c(0, 2), groups = c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"))

sr.norm <- normalizeCtData(qPCRdata.raw, norm = "scale.rank")
q.norm <- normalizeCtData(qPCRdata.raw, norm = "quantile")
dCT.norm <- normalizeCtData(qPCRdata.raw, norm = "deltaCt", deltaCt.genes = rownames(qPCRdata.raw)[85:89], verbose = T)

plotCtDensity(dCT.norm)
plotCtDensity(qPCRdata.raw)
plotCtDensity(sr.norm)
plotCtDensity(q.norm)
plotCtBoxes(sr.norm, stratify = "type")
plotCtBoxes(dCT.norm, stratify = "type")
plotCtBoxes(q.norm, stratify = "type")

plotCtScatter(sr.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(sr.norm, cards = c(3,4), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(3,4), col = "type", diag = T)

plotCtPairs(sr.norm, col = "type", diag = TRUE)
plotCtPairs(dCT.norm, col = "type", diag = TRUE)

plotCVBoxes(qPCRdata.raw, stratify = "type")
plotCVBoxes(sr.norm, stratify = "type")
plotCVBoxes(dCT.norm, stratify = "type")

plotCtHeatmap(qPCRdata.raw, gene.names = "", dist = "euclidean")
plotCtHeatmap(dCT.norm, dist = "euclidean")

clusterCt(dCT.norm, type = "samples")
clusterCt(dCT.norm, type = "genes")
cluster.list <- clusterCt(dCT.norm, type = "genes", n.cluster = 4, cex = 0.2)
cluster.list[[3]]

plotCtPCA(qPCRdata.raw, features = F)
groups <- as.factor(c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"))
qDE.ttest <- ttestCtData(dCT.norm, groups = groups)
qDE.ttest$genes <- as(qDE.ttest$genes, "character")
qDE.ttest$ensembl_gene_id <- qPCRGenesTab[qDE.ttest$genes, "ensembl_gene_id"]
qDE.ttest <- qDE.ttest[- which(is.na(qDE.ttest$ensembl_gene_id)), ]
rownames(qDE.ttest) <- qDE.ttest$ensembl_gene_id
rownames(qDE.ttest) <- qDE.ttest$genes

design <- model.matrix(~0 + files$Treatment)
colnames(design) <- c("Control","Treatment")
contrasts <- makeContrasts(Control - Treatment, levels = design)
dCT.norm2 <- dCT.norm[order(featureNames(dCT.norm)),]
# qDE.limma <- limmaCtData(dCT.norm2, design = design, contrasts = contrasts, ndups = 1, spacing = 1)
# qDE.limma.fit <- limmaCtDataMod(dCT.norm2, design = design, contrasts = contrasts, ndups = 1, spacing = 1, ret.fit = FALSE, topTableOut = TRUE)
# 
# qDE.limma.tab <- qDE.limma.fit$`topTable`
# qDE.limma.tab$genes <- rownames(qDE.limma.tab)
# 
# save(qDE.limma.tab, file = "~/Data/Tremethick/EMT/MDCK qPCR data/qDE.limma.tab.rda")
# write.csv(qDE.limma.tab, "~/Data/Tremethick/EMT/MDCK qPCR data/qDE.limma.tab.TGFb_treated.csv")

qDE.TGFb.limma.tab <- read.csv("~/Data/Tremethick/EMT/qPCR_data_analysis/qDE.limma.tab.TGFb_treated.csv", row.names = 1)
write.csv(qDE.TGFb.limma.tab[which(qDE.TGFb.limma.tab$logFC > 0), ], file = "~/Data/Tremethick/EMT/qPCR_data_analysis/qDE.limma.tab.up.TGFb_treated.csv")
write.csv(qDE.TGFb.limma.tab[which(qDE.TGFb.limma.tab$logFC < 0), ], file = "~/Data/Tremethick/EMT/qPCR_data_analysis/qDE.limma.tab.down.TGFb_treated.csv")

deepToolsUtils::WriteGRangesToBED(gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% qDE.TGFb.limma.tab[which(qDE.TGFb.limma.tab$logFC >= 0), ]$genes)], 
                                  out_file = "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/qPCR_genes_up.bed")
deepToolsUtils::WriteGRangesToBED(gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% qDE.TGFb.limma.tab[which(qDE.TGFb.limma.tab$logFC < 0), ]$genes)], 
                                  out_file = "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/qPCR_genes_down.bed")

gr.top5EMTUp <- gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% qDE.TGFb.limma.tab[order(qDE.TGFb.limma.tab$logFC, decreasing = T),]$genes[1:5])]
save(gr.top5EMTUp, file = "~/Data/Tremethick/EMT/GenomeWide/gr.top5EMTUp.rda")

gr.top5EMTDown <- gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% qDE.TGFb.limma.tab[order(qDE.TGFb.limma.tab$logFC, decreasing = F),]$genes[1:5])]
save(gr.top5EMTDown, file = "~/Data/Tremethick/EMT/GenomeWide/gr.top5EMTDown.rda")

# plot log2 FC of mesenchymal/epithelial markers
ptab <- qDE.TGFb.limma.tab[c("FN1", "ZEB1", "TGFB1", "TGFB2", "TGFB3", "SPARC", "TWIST1", "CDH1", "SPP1", "FGFBP1", "MMP9"), c("genes", "logFC", "CI.L", "CI.R")]
ptab$marker <- c(rep("mesenchymal", 7), rep("epithelial", 4))
ptab$marker <- factor(ptab$marker)
ptab$geneOrd <- reorder(ptab$gene, ptab$logFC)
colnames(ptab)[c(3,4)] <- c("CI_L", "CI_R")

marker <- factor(levels(ptab$marker), levels = levels(ptab$marker))
limits <- aes(ymax = logFC + CI_R, ymin = logFC - CI_L)
p1 <- ggplot(data = ptab, aes(y = logFC, fill = marker, group = marker, x = geneOrd)) 
pdf("Barplot_EMT_markers_TGFb_treated.pdf")
p1 + geom_bar(stat = "identity") + ggtitle("TGFb-treated vs Control MDCK cells") + scale_y_continuous(limits = c(-5, 8))
#geom_errorbar(aes(ymin = logFC-CI_L, ymax = logFC+CI_R), width = 0.2, position = position_dodge(0.9))
dev.off()

# Volcano plot of limma processed qPCR data -------------------------------
pdf("Volcano_plot.pdf", height = 4, width = 4)
plot(qDE.TGFb.limma.tab$logFC, -log10(qDE.TGFb.limma.tab$adj.P.Val), axes = F, xlab = "", ylab = "", frame = F, xlim = c(-10,10), cex = 0.3, pch = 16)
points(qDE.TGFb.limma.tab[which(-log10(qDE.TGFb.limma.tab$adj.P.Val) >= 1), "logFC"], -log10(qDE.TGFb.limma.tab[which(-log10(qDE.TGFb.limma.tab$adj.P.Val) >= 1), "adj.P.Val"]), col = "red", pch = 16, cex = 1.2)
axis(2, pos = 0, lwd = 3, at = c(seq(0,8,2)), labels = c("", seq(2,8,2)))
axis(1, pos = 0, lwd = 3)
mtext("-log10(adjusted p-value)", side = 2)
mtext("log2(FC)", side = 1, line = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)
# only adding labels to genes with adjuste p-value <= 0.01
# text(qDE.TGFb.limma.tab[which(-log10(qDE.TGFb.limma.tab$adj.P.Val) >= 2), "logFC"], -log10(qDE.TGFb.limma.tab[which(-log10(qDE.TGFb.limma.tab$adj.P.Val) >= 2), "adj.P.Val"]), 
#     labels = qDE.TGFb.limma.tab[which(-log10(qDE.TGFb.limma.tab$adj.P.Val) >= 2), "genes" ],
#     cex = 0.7,
#     pos = 4, offset = 0.3)
#text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
dev.off()

#load("danpos2.anno.rda")
danpos2.anno[which(danpos2.anno$feature %in% qPCRGenesTab$ensembl_gene_id)]

#-----------Analysis of H2AZ-knockdown experiment--------------------
qPCRData.H2AZ_KD <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/H2AZ-knockdown-qPCR-data.csv", header = T, row.names = 1, as.is = T)
colnames(qPCRData.H2AZ_KD) <- c(paste("Control_", c(1:4), sep= ""), paste("H2AZ_Knockdown_", c(1:4), sep= ""))

lapply(seq_along(1:ncol(qPCRData.H2AZ_KD)), function(x) {
  Ct <- qPCRData.H2AZ_KD[,x]
  genes <- rownames(qPCRData.H2AZ_KD)
  pos <- c(1:nrow(qPCRData.H2AZ_KD))
  flag <- rep("OK", nrow(qPCRData.H2AZ_KD))
  type <- c(rep("GOI", 84), rep("HKG", 5))
  fn <- paste("~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/", colnames(qPCRData.H2AZ_KD)[x], ".csv", sep = "")
  print(fn)
  tab <- data.frame(cbind(genes, Ct, pos, flag, type ))
  write.table(tab, file = fn, row.names = F, col.names = F, sep = "\t")
})

files <- data.frame(ID = 1:ncol(qPCRData.H2AZ_KD))
files$Files <- unlist(lapply(seq_along(1:ncol(qPCRData.H2AZ_KD)), function(x) { paste(colnames(qPCRData.H2AZ_KD)[x], ".csv", sep = "") }))
files$Treatment <- c("Control", "Control", "Control", "Control", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown")

column.info <- list(flag = "flag", Ct = "Ct", feature = "genes", position = "pos", type = "type")

qPCRdata.raw <- readCtData(files$Files, 
                           format = "plain", 
                           path = "~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/", 
                           n.features = 89, 
                           column.info = list(flag = 4, Ct = 2, feature = 1, position = 3, type = 5), 
                           n.data = 1,
                           header = F)

pData(qPCRdata.raw) <- data.frame(Sample = c("Control", "Control", "Control", "Control", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown"), Replicate = rep(1:4,2 ))

g <- featureNames(qPCRdata.raw)[1:10]
g <- c("TGFB1", "MMP9", "CDH1", "H2AFZ")
plotCtOverview(qPCRdata.raw, calibrator = "Control", genes = g, conf.int = T, ylim = c(0, 2), groups = c("Control", "Control", "Control", "Control", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown"))

sr.norm <- normalizeCtData(qPCRdata.raw, norm = "scale.rank")
q.norm <- normalizeCtData(qPCRdata.raw, norm = "quantile")
dCT.norm <- normalizeCtData(qPCRdata.raw, norm = "deltaCt", deltaCt.genes = rownames(qPCRdata.raw)[85:89], verbose = T)

plotCtDensity(dCT.norm)
plotCtDensity(qPCRdata.raw)
plotCtDensity(sr.norm)
plotCtDensity(q.norm)
plotCtBoxes(sr.norm, stratify = "type")
plotCtBoxes(dCT.norm, stratify = "type")
plotCtBoxes(q.norm, stratify = "type")

plotCtScatter(sr.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(sr.norm, cards = c(3,4), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(3,4), col = "type", diag = T)

plotCtPairs(sr.norm, col = "type", diag = TRUE)
plotCtPairs(dCT.norm, col = "type", diag = TRUE)

plotCVBoxes(qPCRdata.raw, stratify = "type")

plotCtHeatmap(qPCRdata.raw, gene.names = "", dist = "euclidean")
plotCtHeatmap(dCT.norm, gene.names = "", dist = "euclidean")
clusterCt(dCT.norm, type = "samples")
clusterCt(dCT.norm, type = "genes")
cluster.list <- clusterCt(dCT.norm, type = "genes", n.cluster = 4, cex = 0.2)
cluster.list[[3]]

plotCtPCA(qPCRdata.raw, features = F)
groups <- as.factor(c("Control", "Control", "Control", "Control", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown", "H2AZ_Knockdown"))
qDE.ttest <- ttestCtData(dCT.norm, groups = groups)
qDE.ttest$genes <- as(qDE.ttest$genes, "character")
qDE.ttest$ensembl_gene_id <- qPCRGenesTab[qDE.ttest$genes, "ensembl_gene_id"]
qDE.ttest <- qDE.ttest[- which(is.na(qDE.ttest$ensembl_gene_id)), ]
rownames(qDE.ttest) <- qDE.ttest$ensembl_gene_id

design <- model.matrix(~0 + files$Treatment)
colnames(design) <- c("Control","Treatment")
contrasts <- makeContrasts(Control - Treatment, levels = design)
dCT.norm2 <- dCT.norm[order(featureNames(dCT.norm)),]
qDE.limma <- limmaCtData(dCT.norm2, design = design, contrasts = contrasts, ndups = 1, spacing = 1)
qDE.limma.fit <- limmaCtDataMod(dCT.norm2, design = design, contrasts = contrasts, ndups = 1, spacing = 1, ret.fit = FALSE, topTableOut = TRUE)

qDE.limma.tab <- qDE.limma.fit$`topTable`
qDE.limma.tab$genes <- rownames(qDE.limma.tab)

qDE.shZKD.limma.tab <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/qDE.limma.tab.H2AZ_KD.csv")
write.csv(qDE.shZKD.limma.tab[which(qDE.shZKD.limma.tab$logFC >= 0), ], file = "~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/qDE.limma.tab.up.H2AZ_KD.csv")
write.csv(qDE.shZKD.limma.tab[which(qDE.shZKD.limma.tab$logFC < 0), ], file = "~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/qDE.limma.tab.down.H2AZ_KD.csv")

# plot log2 FC of mesenchymal/epithelial markers
ptab <- qDE.shZKD.limma.tab[c("FN1", "ZEB1", "TGFB1", "TGFB2", "TGFB3", "SPARC", "TWIST1", "CDH1", "SPP1", "FGFBP1", "MMP9"), c("genes", "logFC", "CI.L", "CI.R")]
ptab$marker <- c(rep("mesenchymal", 7), rep("epithelial", 4))
ptab$marker <- factor(ptab$marker)
ptab$geneOrd <- reorder(ptab$gene, ptab$logFC)
colnames(ptab)[c(3,4)] <- c("CI_L", "CI_R")

marker <- factor(levels(ptab$marker), levels = levels(ptab$marker))
limits <- aes(ymax = logFC + CI_R, ymin = logFC - CI_L)
p1 <- ggplot(data = ptab, aes(y = logFC, fill = marker, group = marker, x = geneOrd)) 
pdf("Barplot_EMT_markers_H2AZ_Knockdown.pdf") 
p1 + geom_bar(stat = "identity") + ggtitle("H2A.Z Knockdown vs scramble siRNA MDCK cells") + scale_y_continuous(limits = c(-5, 8))
#   geom_errorbar(aes(ymin = logFC-CI_L, ymax = logFC+CI_R), 
#                 width = 0.2, 
#                 position = position_dodge(0.9))
dev.off()

# Volcano plot of qPCR data -----------------------------------------------
pdf("Volcano_plot_H2AZ_KD.pdf", height = 8, width = 8)
plot(qDE.shZKD.limma.tab$logFC, -log10(qDE.shZKD.limma.tab$adj.P.Val), axes = F, xlab = "", ylab = "", frame = F, xlim = c(-10,10), cex = 0.3, pch = 16)
points(qDE.shZKD.limma.tab[which(-log10(qDE.shZKD.limma.tab$adj.P.Val) >= 1), "logFC"], -log10(qDE.shZKD.limma.tab[which(-log10(qDE.shZKD.limma.tab$adj.P.Val) >= 1), "adj.P.Val"]), col = "red", pch = 16, cex = 1.2)
axis(2, pos = 0, lwd = 3, at = c(seq(0,8,2)), labels = c("", seq(2,8,2)))
axis(1, pos = 0, lwd = 3)
mtext("-log10(adjusted p-value)", side = 2)
mtext("log2(FC)", side = 1, line = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)
# only adding labels to genes with adjuste p-value <= 0.01
 text(qDE.shZKD.limma.tab[which(-log10(qDE.shZKD.limma.tab$adj.P.Val) >= 2), "logFC"], -log10(qDE.shZKD.limma.tab[which(-log10(qDE.shZKD.limma.tab$adj.P.Val) >= 2), "adj.P.Val"]), 
     labels = qDE.shZKD.limma.tab[which(-log10(qDE.shZKD.limma.tab$adj.P.Val) >= 2), "genes" ],
     cex = 0.7,
     pos = 4, offset = 0.3)
text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
dev.off()



#load("danpos2.anno.rda")
danpos2.anno[which(danpos2.anno$feature %in% qPCRGenesTab$ensembl_gene_id)]


# scratchpad --------------------------------------------------------------
# Volcano plot of ttest data
plot(log2(qDE.ttest$FC), -log10(qDE.ttest$adj.p.value), axes = F, xlab = "", ylab = "", frame = F, xlim = c(-10,10), cex = 0.3, pch = 16)
points(log2(qDE.ttest[which(-log10(qDE.ttest$adj.p.value) >= 1), "FC"]), 
       -log10(qDE.ttest[which(-log10(qDE.ttest$adj.p.value) >= 1), "adj.p.value"]), col = "red", pch = 16, cex = 1.2)
axis(2, pos = 0, lwd = 3, at = c(seq(0,8,2)), labels = c("", seq(2,8,2)))
axis(1, pos = 0, lwd = 3)
mtext("-log10(adjusted p-value)", side = 2)
mtext("log2(FC)", side = 1, line = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)

plot(qDE.shZKD.limma.tab$logFC, log2(qDE.ttest$FC))


