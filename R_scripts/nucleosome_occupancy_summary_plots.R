require(ggplot2)
source("~/Development/JCSMR-Tremethick-Lab/General_Purpose/multiplot.R")

df3 <- as(log2(mcols(danpos2.anno.minCov )$"control_smt_val" + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(danpos2.anno.minCov )$"treat_smt_val" + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep("ctrl", length(danpos2.anno.minCov )), rep("treat", length(danpos2.anno.minCov)))
hm2 <- ggplot(df3,aes(x=df3, group=var))
hm2 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))


# all promoters - minimum coverage
gr3 <- subsetByOverlaps(danpos2.anno.minCov, promoters(Cfam3.genes, upstream = 1000, downstream = 1000))
gr3 <- gr3[which(mcols(gr3)$smt_diff_FDR <= 0.01)]
# create a histogram of the data (here log2 transformed)
df3 <- as(log2(mcols(gr3)[, c("control_smt_val")] + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(gr3)[, c("treat_smt_val")] + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep(paste("MDCK - Untreated [N = ", length(gr3), "]", sep = ""), length(gr3)), rep(paste("MDCK - TGFb-treated [N = ", length(gr3), "]", sep = ""), length(gr3)))
histo3 <- ggplot(df3,aes(x=df3, group=var))
histo3 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
p1 <- histo3 + geom_boxplot(position = "identity", aes(y = df3, x= var)) + labs(title = "All genes +/- 1Kb TSS\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))

# all promoters - no minimum coverage
gr3 <- subsetByOverlaps(danpos2.anno, promoters(Cfam3.genes, upstream = 1000, downstream = 1000))
gr3 <- gr3[which(mcols(gr3)$smt_diff_FDR <= 0.01)]
df3 <- as(log2(mcols(gr3)[, c("control_smt_val")] + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(gr3)[, c("treat_smt_val")] + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep(paste("MDCK - Untreated [N = ", length(gr3), "]", sep = ""), length(gr3)), rep(paste("MDCK - TGFb-treated [N = ", length(gr3), "]", sep = ""), length(gr3)))
histo3 <- ggplot(df3,aes(x=df3, group=var))
histo3 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))
p2 <- histo3 + geom_boxplot(position = "identity", aes(y = df3, x= var)) + labs(title = "All genes +/- 1Kb TSS\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))
multiplot(p1, p2, cols = 2)


# qPCR genes
gr3 <- subsetByOverlaps(danpos2.anno.minCov, promoters(Cfam3.genes[which(Cfam3.genes$gene_id %in% qPCRGenesTab$ensembl_gene_id)], upstream = 1000, downstream = 1000))
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
histo3 + geom_boxplot(position = "identity", aes(y = df3, x= var)) + labs(title = "qPCR genes\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))
dev.off()

