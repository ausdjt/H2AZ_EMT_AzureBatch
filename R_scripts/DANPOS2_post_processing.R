# DANPOS2 post-processing

#------------load libraries------------------
library("GenomicFeatures")
library("ChIPpeakAnno")
library("ggplot2")
library("Gviz")
library("GenomicRanges")
library("rtracklayer")
library("GenomicFeatures")
require(org.Cf.eg.db)

setwd("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/")

#------------import external functions------------------
source("~/Development/GeneralPurpose/R/heatmap.3.R")


# create TxDb object of the Canis familiaris 3.1 genome annotation --------
# prepare annotation data &
# create Canis familiaris TXDB object for peak annotation
chromInfo <- read.table("~/mount/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/chromInfo.txt", header = F, as.is = T, sep = "\t")
colnames(chromInfo) <- c("chrom", "length")
TxDb.Cfam3.Ensembl <- makeTxDbFromGFF("~/Data/Annotations/CanFam3/Ensembl/Canis_familiaris.CanFam3.1.83.chr.gtf.gz", 
                                      organism = "Canis familiaris", 
                                      chrominfo = chromInfo)
# TxDb.Cfam3.RefSeq <- makeTxDbFromGFF("~/Data/Annotations/CanFam3/RefSeq/ref_CanFam3.1_top_level.gff3", 
#                                       organism = "Canis familiaris", 
#                                       chrominfo = chromInfo)

Cfam3.genes <- genes(TxDb.Cfam3.Ensembl)
Cfam3.genes.hgnc <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = Cfam3.genes$gene_id, mart = dog)
rownames(Cfam3.genes.hgnc) <- Cfam3.genes.hgnc$ensembl_gene_id
Cfam3.genes$hgnc_symbol <- Cfam3.genes.hgnc[Cfam3.genes$gene_id, ]$hgnc_symbol
save(Cfam3.genes, file = "Cfam3.genes.rda")

load("~/Data/Tremethick/EMT/GenomeWide/Cfam3.genes.rda")


# load repeat regions to exclude from analysis of nucleosome posit --------
Cfam3.repeats <- import("~/Data/Annotations/CanFam3/canFam3_repeat_regions.bed")
seqlevels(Cfam3.repeats) <- gsub("chr", "", seqlevels(Cfam3.repeats))
seqlevels(Cfam3.repeats, force = T) <- seqlevels(Cfam3.genes)
seqinfo(Cfam3.repeats, force = T) <- seqinfo(Cfam3.genes)

# load DANPOS2 results ----------------------------------------------------
# changed 2015-12-14:
# now using the results from DANPOS2 analysis with default settings
# (~/Development/JCSMR-Tremethick-Lab/shell_scripts/danpos2_command_lines.sh)
# previous run seemed to create too large sliding windows, i.e. twice nucleosome size
danpos2.results <- read.table("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT_147bp/result/TGFb_H2AZ-WT_H2AZ.positions.integrative.xls",
                              header = T, 
                              as.is = T,
                              sep = "\t")
gr.danpos2.results <- GRanges(danpos2.results$chr, IRanges(danpos2.results$start, danpos2.results$end), strand = "*", danpos2.results[, c(4:23)])

# Using Gviz for visualization of some of the data
i <- 2
gr1 <- promoters(gr.mesenchymalMarkers.genes[i], upstream = 10000, downstream = 10000)
dT.WT <- DataTrack(subsetByOverlaps(gr.WT_H2AZ_ChIP_bgsub_Fnor, gr1), type = "h", col = "black", name = "Control")
dT.TGFb <- DataTrack(subsetByOverlaps(gr.TGFb_H2AZ_ChIP_bgsub_Fnor, gr1), type = "h", col = "grey", name = "TGFb")
dT.Diff <- DataTrack(subsetByOverlaps(gr.TGFb_vs_WT_diff, gr1), type = "h", col = "blue", name = "Difference [+/- log10p-val]")

biomTrack.ensembl <- BiomartGeneRegionTrack(genome = "canFam3", 
                                    chromosome = as(seqlevels(gr1)[i], "character"),
                                    start = as(start(gr1), "integer"),
                                    end = as(end(gr1), "integer"),
                                    name = paste(mcols(gr1)$hgnc_symbol),
                                    mart = dog)
gat <- GenomeAxisTrack()
plotTracks(list(gat, biomTrack, dT.WT, dT.TGFb, dT.Diff))


# trying to directly visualize the DANPOS2 resuls from the integrative presentation of data
DT1 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$control_smt_val, type = "l")
DT2 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$treat_smt_val, type = "l")
DT3 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$smt_log2FC,type = "l")
DT4 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = -1 * log10(mcols(subsetByOverlaps(gr.danpos2.results, gr1))$smt_diff_FDR), type = c("p", "g"))

max.y <- max(max(values(DT1)), max(values(DT2)))
displayPars(DT1) <- list(ylim = c(0,max.y))
displayPars(DT2) <- list(ylim = c(0,max.y))

  
plotTracks(list(gat, biomTrack, DT1, DT2, DT3, DT4), from = start(gr1), to = end(gr1))

df1 <- as(mcols(subsetByOverlaps(gr.danpos2.results, gr1[1]))$control_smt_val, "matrix")
df1 <- rbind(df1, as(mcols(subsetByOverlaps(gr.danpos2.results, gr1[1]))$treat_smt_val, "matrix"))
df1 <- data.frame(df1)
df1$var <- c(rep("ctrl", length(subsetByOverlaps(gr.danpos2.results, gr1[1]))), 
             rep("treat", length(subsetByOverlaps(gr.danpos2.results, gr1[1]))))
df1$pos <- rep(c(1:length(subsetByOverlaps(gr.danpos2.results, gr1[1]))), 2)
p <- ggplot(df1, aes(x = pos, y = df1 , group = var, colour = var))

# annotate DANPOS2 peaks --------------------------------------------------
danpos2.anno <- annotatePeakInBatch(gr.danpos2.results, AnnotationData = Cfam3.genes)
danpos2.anno <- addGeneIDs(annotatedPeak = danpos2.anno, orgAnn = org.Cf.eg.db, feature_id_type = "ensembl_gene_id", IDs2Add = c("symbol", "entrez_id"))

# load pre-computed data
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/danpos2.anno.rda")
seqlevels(danpos2.anno) <- gsub("chr", "", seqlevels(danpos2.anno))
seqlevels(danpos2.anno)[grep("M", seqlevels(danpos2.anno))] <- "MT"
seqlevels(danpos2.anno, force = T) <- seqlevels(Cfam3.genes)
seqinfo(danpos2.anno, force = T) <- seqinfo(Cfam3.genes)

# filter out peak positions that do overlap with annotated repeats
danpos2.anno <- danpos2.anno[!overlapsAny(danpos2.anno, Cfam3.repeats)]

# filter out nucleosomes with low coverage
danpos2.anno.minCov <- danpos2.anno[which(danpos2.anno$control_smt_val > 10 & danpos2.anno$treat_smt_val > 10)]

# only consider peaks/nucleosome position upstream or on the annotated TSS
upTSS <- which(mcols(danpos2.anno)$insideFeature %in% c("overlapStart", "upstream"))
gr.danpos2.upTSS <- danpos2.anno[upTSS]
hist(mcols(gr.danpos2.upTSS)$smt_diff_FDR)
# gr.danpos2.upTSS <- gr.danpos2.upTSS[which(mcols(gr.danpos2.upTSS)$smt_diff_FDR < 0.001)]
gr.danpos2.upTSS.5kb <- gr.danpos2.upTSS[which(mcols(gr.danpos2.upTSS)$distancetoFeature <= 5000)]

df3 <- as(log2(mcols(gr.danpos2.upTSS)$"control_smt_val" + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(gr.danpos2.upTSS)$"treat_smt_val" + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep("ctrl", length(gr.danpos2.upTSS)), rep("treat", length(gr.danpos2.upTSS)))
hm2 <- ggplot(df3,aes(x=df3, group=var))
hm2 + geom_histogram(alpha = 0.6, position = "identity", group= var)
hm2 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))

df3.1 <- mcols(gr.danpos2.upTSS)[c("treat_smt_val", "control_smt_val")]
df3.1 <- as.matrix(df3.1)
heatmap.3(log2(df3.1 + 1), trace = "none", cexCol = 0.7)

#------------check GO enrichment of peaks------------------
# look at control and treatment enriched separately - here for all nucleosome positions [with minimum log2 2-fold change]
# control results
GO.danpos2.anno_ctrl_enriched <- getEnrichedGO(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC < -2)], 
                                               orgAnn = "org.Cf.eg.db", 
                                               maxP=0.05, 
                                               multiAdj = F, 
                                               minGOterm = 10, 
                                               multiAdjMethod = "BH")
# collapse results at GO term level
go.id1 <- unique(as.character(GO.danpos2.anno_ctrl_enriched[["bp"]]$go.id))
go.term1 <- unique(as.character(GO.danpos2.anno_ctrl_enriched[["bp"]]$go.term))
d1 <- data.frame(GO.danpos2.anno_ctrl_enriched[["bp"]])
d1$go.id <- as.character(d1$go.id)
p1 <- sapply(go.id1, function(x) {w1 <- which(d1$go.id == x); p1 <- unique(d1[w1, "pvalue"])})
GO.danpos2.anno_ctrl_enriched.BP <- data.frame(cbind("go.id" = as.character(go.id1), "go.term" = as.character(go.term1), pvalue = as.vector(p1)))
GO.danpos2.anno_ctrl_enriched.BP$pvalue <- as.numeric(as.character(GO.danpos2.anno_ctrl_enriched.BP$pvalue))
GO.danpos2.anno_ctrl_enriched.BP$FDR <- p.adjust(GO.danpos2.anno_ctrl_enriched.BP$pvalue, method = "fdr")

# treatment results
GO.danpos2.anno_treat_enriched <- getEnrichedGO(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC > 2)],
                                                orgAnn = "org.Cf.eg.db", 
                                                maxP=0.05, 
                                                multiAdj = F, 
                                                minGOterm = 10, 
                                                multiAdjMethod = "BH")
# collapse results at GO term level
go.id1 <- unique(as.character(GO.danpos2.anno_treat_enriched[["bp"]]$go.id))
go.term1 <- unique(as.character(GO.danpos2.anno_treat_enriched[["bp"]]$go.term))
d1 <- data.frame(GO.danpos2.anno_treat_enriched[["bp"]])
d1$go.id <- as.character(d1$go.id)
p1 <- sapply(go.id1, function(x) {w1 <- which(d1$go.id == x); p1 <- unique(d1[w1, "pvalue"])})
GO.danpos2.anno_treat_enriched.BP <- data.frame(cbind("go.id" = as.character(go.id1), "go.term" = as.character(go.term1), pvalue = as.vector(p1)))
GO.danpos2.anno_treat_enriched.BP$pvalue <- as.numeric(as.character(GO.danpos2.anno_treat_enriched.BP$pvalue))
GO.danpos2.anno_treat_enriched.BP$FDR <- p.adjust(GO.danpos2.anno_treat_enriched.BP$pvalue, method = "fdr")

# check enrichment for nucleosome position up to 5kb upstream of TSS
# control results
GO.danpos2.ctrl.upTSS.5kb <- getEnrichedGO(gr.danpos2.upTSS.5kb[which(mcols(gr.danpos2.upTSS.5kb)$smt_log2FC < -2)],
                                           orgAnn = "org.Cf.eg.db", 
                                           maxP=0.05, 
                                           multiAdj = F, 
                                           minGOterm = 10, 
                                           multiAdjMethod = "BH")
# collapse results at GO term level
go.id1 <- unique(as.character(GO.danpos2.ctrl.upTSS.5kb[["bp"]]$go.id))
go.term1 <- unique(as.character(GO.danpos2.ctrl.upTSS.5kb[["bp"]]$go.term))
d1 <- data.frame(GO.danpos2.ctrl.upTSS.5kb[["bp"]])
d1$go.id <- as.character(d1$go.id)
p1 <- sapply(go.id1, function(x) {w1 <- which(d1$go.id == x); p1 <- unique(d1[w1, "pvalue"])})
GO.danpos2.ctrl.upTSS.5kb.BP <- data.frame(cbind("go.id" = as.character(go.id1), "go.term" = as.character(go.term1), pvalue = as.vector(p1)))
GO.danpos2.ctrl.upTSS.5kb.BP$pvalue <- as.numeric(as.character(GO.danpos2.ctrl.upTSS.5kb.BP$pvalue))
GO.danpos2.ctrl.upTSS.5kb.BP$FDR <- p.adjust(GO.danpos2.ctrl.upTSS.5kb.BP$pvalue, method = "fdr")

# treatment results
GO.danpos2.treat.upTSS.5kb <- getEnrichedGO(gr.danpos2.upTSS.5kb[which(mcols(gr.danpos2.upTSS.5kb)$smt_log2FC > 2)],
                                           orgAnn = "org.Cf.eg.db", 
                                           maxP=0.05, 
                                           multiAdj = F, 
                                           minGOterm = 10, 
                                           multiAdjMethod = "BH")
# collapse results at GO term level
go.id1 <- unique(as.character(GO.danpos2.treat.upTSS.5kb[["bp"]]$go.id))
go.term1 <- unique(as.character(GO.danpos2.treat.upTSS.5kb[["bp"]]$go.term))
d1 <- data.frame(GO.danpos2.treat.upTSS.5kb[["bp"]])
d1$go.id <- as.character(d1$go.id)
p1 <- sapply(go.id1, function(x) {w1 <- which(d1$go.id == x); p1 <- unique(d1[w1, "pvalue"])})
GO.danpos2.treat.upTSS.5kb.BP <- data.frame(cbind("go.id" = as.character(go.id1), "go.term" = as.character(go.term1), pvalue = as.vector(p1)))
GO.danpos2.treat.upTSS.5kb.BP$pvalue <- as.numeric(as.character(GO.danpos2.treat.upTSS.5kb.BP$pvalue))
GO.danpos2.treat.upTSS.5kb.BP$FDR <- p.adjust(GO.danpos2.treat.upTSS.5kb.BP$pvalue, method = "fdr")


# check for enrichment of REACTOME pathways
reactome.danpos2.anno_ctrl_enriched <- getEnrichedPATH(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC < -2)], 
                                                       orgAnn = "org.Cf.eg.db", pathAnn = "reactome.db", 
                                                       maxP = 0.1, 
                                                       minPATHterm = 10, 
                                                       feature_id_type = "ensembl_gene_id")
path1 <- path1 <- unique(as.character(reactome.danpos2.anno_ctrl_enriched$PATH))
p1 <- sapply(path1, function(x) {w1 <- which(reactome.danpos2.anno_ctrl_enriched$PATH == x); p1 <- unique(reactome.danpos2.anno_ctrl_enriched[w1, "pvalue"])})
reactome.danpos2.anno_ctrl_enriched.pathLevel <- data.frame(pathway = cbind(as.character(path1), pvalue = as.vector(p1)))
reactome.danpos2.anno_ctrl_enriched.pathLevel$pvalue <- as.numeric(as.character(reactome.danpos2.anno_ctrl_enriched.pathLevel$pathway.pvalue))
reactome.danpos2.anno_ctrl_enriched.pathLevel$FDR <- p.adjust(reactome.danpos2.anno_ctrl_enriched.pathLevel$pvalue, method = "fdr")

reactome.danpos2.anno_treat_enriched <- getEnrichedPATH(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC > 2)], 
                                                        orgAnn = "org.Cf.eg.db", 
                                                        pathAnn = "reactome.db", 
                                                        maxP = 0.1, 
                                                        minPATHterm = 10, 
                                                        feature_id_type = "ensembl_gene_id")

path1 <- unique(as.character(reactome.danpos2.anno_treat_enriched$PATH))
p1 <- sapply(path1, function(x) {w1 <- which(reactome.danpos2.anno_treat_enriched$PATH == x); p1 <- unique(reactome.danpos2.anno_treat_enriched[w1, "pvalue"])})
reactome.danpos2.anno_treat_enriched.pathLevel <- data.frame(pathway = cbind(as.character(path1), pvalue = as.vector(p1)))
reactome.danpos2.anno_treat_enriched.pathLevel$pvalue <- as.numeric(as.character(reactome.danpos2.anno_treat_enriched.pathLevel$pathway.pvalue))
reactome.danpos2.anno_treat_enriched.pathLevel$FDR <- p.adjust(reactome.danpos2.anno_treat_enriched.pathLevel$pvalue, method = "fdr")

#------------MSigDB analysis------------------
seqlevels(danpos2.anno) <- gsub("chr", "", seqlevels(danpos2.anno))
gr1 <- subsetByOverlaps(danpos2.anno, promoters(gr.MSigDB.EMT_associated.cfam, upstream = 500, downstream = 0))
gr1 <- gr1[which(mcols(gr1)$smt_diff_FDR <= 0.01)]
gr1.ctrl <- gr1[which(mcols(gr1)$smt_log2FC < 0)]
gr1.treat <- gr1[which(mcols(gr1)$smt_log2FC > 0)]

# create a histogram of the data (here log2 transformed)
df1 <- as(log2(mcols(gr1.ctrl)[, c("control_smt_val")] + 1), "matrix")
df1 <- rbind(df1, as(log2(mcols(gr1.treat)[, c("treat_smt_val")] + 1), "matrix"))
df1 <- data.frame(df1)
df1$var <- c(rep(paste("MDCK - Untreated [N = ", length(gr1.ctrl), "]", sep = ""), length(gr1.ctrl)), rep(paste("MDCK - TGFb-treated [N = ", length(gr1.treat), "]", sep = ""), length(gr1.treat)))
histo1 <- ggplot(df1,aes(x=df1, group=var))
pdf("Histogram_H2AZ_nucleosome_500TSS0_EMT_associated_genes_FDR0.01.pdf", height = 8, width = 8)
histo1 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
dev.off()
pdf("Boxplot_H2AZ_nucleosome_500TSS0_EMT_associated_genes_FDR0.01.pdf", height = 8, width = 8)
histo1 + geom_boxplot(position = "identity", aes(y = df1, x= var)) + labs(title = "EMT-associated genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))
dev.off()

# display data in heatmap
df1 <- as(mcols(gr1)[,c("control_smt_val", "treat_smt_val")], "data.frame")
pdf("Heatmap_H2AZ_nucleosome_500TSS0_EMT_associated_genes_FDR0.01.pdf", height = 8, width = 8)
heatmap1 <- heatmap.3(as.matrix(log2(df1 + 1)), trace = "none", cexCol = 0.6, main = "EMT-associated genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", hclustfun=function(x) hclust(x,method="ward.D"))
dev.off()


# as a "background" data set we look at the totallity of genes
gr2 <- subsetByOverlaps(danpos2.anno, promoters(Cfam3.genes, upstream = 250, downstream = 0))
gr2 <- gr2[which(mcols(gr2)$smt_diff_FDR <= 0.01)]
gr2.ctrl <- gr2[which(mcols(gr2)$smt_log2FC < 0)]
gr2.treat <- gr2[which(mcols(gr2)$smt_log2FC > 0)]
# create a histogram of the data (here log2 transformed)
df2 <- as(log2(mcols(gr2.ctrl)[, c("control_smt_val")] + 1), "matrix")
df2 <- rbind(df2, as(log2(mcols(gr2.treat)[, c("treat_smt_val")] + 1), "matrix"))
df2 <- data.frame(df2)
df2$var <- c(rep(paste("MDCK - Untreated [N = ", length(gr2.ctrl), "]", sep = ""), length(gr2.ctrl)), rep(paste("MDCK - TGFb-treated [N = ", length(gr2.treat), "]", sep = ""), length(gr2.treat)))
histo2 <- ggplot(df2,aes(x=df2, group=var))
pdf("Histogram_H2AZ_nucleosome_500TSS0_all_genes_FDR0.01.pdf", height = 10, width = 10)
histo2 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
dev.off()
pdf("Boxplot_H2AZ_nucleosome_500TSS0_all_genes_FDR0.01.pdf", height = 10, width = 10)
histo2 + geom_boxplot(position = "identity", aes(y = df2, x= var)) + labs(title = "All genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))
dev.off()

df2 <- as(mcols(gr2)[,c("control_smt_val", "treat_smt_val")], "data.frame")
pdf("Heatmap_H2AZ_nucleosome_500TSS0_all_genes_FDR0.01.pdf", height = 10, width = 10)
heatmap2 <- heatmap.3(as.matrix(log2(df2 + 1)), trace = "none", cexCol = 0.6, main = "All genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", hclustfun=function(x) hclust(x,method="ward.D"))
dev.off()

# trying to determine the different clusters, so that we can identify genes with strong differential H2A.Z occupation
hcrow2 <- as.hclust(heatmap2$rowDendrogram)
silhouette2 <- list()
for(i in 1:9){
  silhouette2[[i]] <- silhouette(cutree(hcrow2, k = i+1), daisy(as.matrix(log2(df2 + 1))))
}

pdf("shilhouetteHeatAutosome.pdf",pointsize=10, height=8,width=8)
par(mfrow=c(3,3))
for(i in 1:9)
  plot(silhouette2[[i]])
dev.off()

c1 <- as.data.frame(cbind(autok2=cutree(hccol, k = 2),autok3=cutree(hccol, k = 3),autok4=cutree(hccol, k = 4)))

# quick look at qPCR genes
Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)]

# changed region to TSS +/-1Kb 
gr3 <- subsetByOverlaps(danpos2.anno, promoters(Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)], upstream = 1000, downstream = 1000))
gr3 <- gr3[which(mcols(gr3)$smt_diff_FDR <= 0.01)]
# create a histogram of the data (here log2 transformed)
df3 <- as(log2(mcols(gr3)[, c("control_smt_val")] + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(gr3)[, c("treat_smt_val")] + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep(paste("MDCK - Untreated [N = ", length(gr3), "]", sep = ""), length(gr3)), rep(paste("MDCK - TGFb-treated [N = ", length(gr3), "]", sep = ""), length(gr3)))
histo3 <- ggplot(df3,aes(x=df3, group=var))
pdf("Histogram_H2AZ_nucleosome_500TSS0_all_genes_FDR0.01.pdf", height = 10, width = 10)
histo3 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
dev.off()
pdf("Boxplot_H2AZ_nucleosome_500TSS0_aPCR_genes_FDR0.01.pdf", height = 10, width = 10)
histo3 + geom_boxplot(position = "identity", aes(y = df3, x= var)) + labs(title = "qPCR genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))
dev.off()

df3 <- as(mcols(gr3)[,c("control_smt_val", "treat_smt_val", "symbol")], "data.frame")
pdf("Heatmap_H2AZ_nucleosome_500TSS0_qPCR_genes_FDR0.01.pdf", height = 10, width = 10)
heatmap2 <- heatmap.3(as.matrix(log2(df3[,c(1,2)] + 1)), 
                      trace = "none", 
                      cexCol = 0.6, 
                      main = "qPCR genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", 
                      hclustfun=function(x) hclust(x,method="ward.D"),
                      labRow = df3[,3])
dev.off()

#---------------plotting coverage maps based on DANPOS2 normalized data------------
gr.TGFb <- subsetByOverlaps(gr.TGFb_H2AZ_ChIP_bgsub_Fnor, gr.which)
gr.WT <- subsetByOverlaps(gr.WT_H2AZ_ChIP_bgsub_Fnor, gr.which)

dT.TGFb <- DataTrack(subsetByOverlaps(gr.TGFb_H2AZ_ChIP_bgsub_Fnor, gr.which), type = "histogram")
dT.WT <- DataTrack(subsetByOverlaps(gr.WT_H2AZ_ChIP_bgsub_Fnor, gr.which), type = "histogram")
i <- 2
biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", 
                                    chromosome = as(seqnames(gr.mesenchymalMarkers), "character")[i],
                                    start = as.integer(start(gr.mesenchymalMarkers[i]), "integer"),
                                    end = as.integer(end(gr.mesenchymalMarkers[i]), "integer"),
                                    name = paste(mcols(gr.mesenchymalMarkers[i])$hgnc_symbol, "transcript",  mcols(gr.mesenchymalMarkers[i])$ensembl_transcript_id, sep = " "),
                                    mart = dog)

aT1 <- AnnotationTrack(gr.mesenchymalMarkers[i])

plotTracks(list(aT1, dT.WT, dT.TGFb))


#-----------extracting regions of EMT-associated genes (FDR <= 0.01), treatment > control------------
gr1 <- subsetByOverlaps(danpos2.anno, promoters(gr.MSigDB.EMT_associated.cfam, upstream = 500, downstream = 0))
gr1[which(gr1$treat_smt_val > 5 & gr1$control_smt_val <= 1)]

# nucleosomes -1000/500+ TSS -----------------------------------------------
gr1 <- subsetByOverlaps(danpos2.anno, promoters(Cfam3.genes, upstream = 1500, downstream = 500))
# Volcanon plot
plot(gr1$smt_log2FC, -log(gr1$smt_diff_FDR))
#gr1 <- gr1[gr1$smt_diff_FDR <= 0.001]

df <- as(mcols(gr1)[,c("control_smt_val", "treat_smt_val")], "data.frame")
hm <- heatmap.3(as.matrix(log2(df + 1)), 
                trace = "none", 
                cexCol = 0.6,
                col = redgreen(19),
                symkey = F,
                main = "All Genes [-1000/500+ bp TSS]\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.001", 
                hclustfun=function(x) hclust(x,method="ward.D"))

df <- as(mcols(gr1)[,c("control_fuzziness_score", "treat_fuzziness_score")], "data.frame")
hm <- heatmap.3(as.matrix(df), 
                trace = "none", 
                cexCol = 0.6,
                col = redgreen(19),
                symkey = F,
                main = "All Genes [-1000/500+ bp TSS]\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.001", 
                hclustfun=function(x) hclust(x,method="ward.D"))

# TODO:
# convert .wig files to .bw and import
#gr.TGFb_H2AZ_ChIP_bgsub_Fnor <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT_147bp/result/pooled//pooled/TGFb_ChIP.bgsub.Fnor.smooth.bw") 
#gr.WT_H2AZ_ChIP_bgsub_Fnor <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/pooled/WT_ChIP.bgsub.Fnor.smooth.bw")
#gr.TGFb_vs_WT_diff <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/diff/TGFb_vs_WT.pois_diff.bw")

heatmap.3(mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.EMT_associated.cfam.tss1500))[,c("control_smt_val", "treat_smt_val")])

# create a histogram of the complete data (here log2 transformed)
df1 <- as(log2(danpos2.results[, c("control_smt_val")] + 0.0001), "matrix")
#df1 <- as((danpos2.results[, c("control_smt_val")]), "matrix")
df1 <- rbind(df1, as(log2(danpos2.results[, c("treat_smt_val")] + 0.0001), "matrix"))
# df1 <- rbind(df1, as((danpos2.results[, c("treat_smt_val")]), "matrix"))

df1 <- data.frame(df1)
df1$var <- c(rep("ctrl", nrow(danpos2.results)), rep("treat", nrow(danpos2.results)))
hm1 <- ggplot(df1,aes(x=df1, group=var))
hm1 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))

# histogram of summit counts of nucleosomes located in TSS+/-1500 of TGFb-induced EMT genes
df2 <- as(mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.TGFb_induced_EMT.cfam.tss1500))[,c("control_smt_val", "treat_smt_val")], "data.frame")
rownames(df2) <- mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.TGFb_induced_EMT.cfam.tss1500))$row_id
heatmap.3(as.matrix(log2(df2 + 1)), trace = "none")

gr.which <- promoters(Cfam3.genes, upstream = 40000, downstream = 40000)
gr.which <- trim(gr.which)
at.danpos2 <- AnnotationTrack(subsetByOverlaps(gr.danpos2.results, gr.which), shape = "box")

