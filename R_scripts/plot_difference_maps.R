require(ggrepel)
library(gridExtra)
require(zoo)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

cmTGFbUp <- dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix
cmWTUp <- dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix 
diffUp <- cmTGFbUp[, 2:301] - cmWTUp[, 2:301]
rownames(diffUp) <- cmTGFbUp$X4


cmTGFbDown <- dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix
cmWTDown <- dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix 
diffDown <- cmTGFbDown[, 2:301] - cmWTDown[, 2:301]
rownames(diffDown) <- diffDown$X4
cmWT <- rbind(cmWTDown, cmWTUp)
cmTGFb <- rbind(cmTGFbDown, cmTGFbUp)
  
diffUp <- as.data.table(diffUp)
diffUp$geneID <- cmTGFbUp$X4
datUp <- melt(diffUp, id.vars = "geneID", measure.vars <- c(1:300))
datUp <- datUp[order(datUp$geneID),]
setkey(datUp, "geneID")

diffDown <- data.table(diffDown)
diffDown$geneID <- cmTGFbDown$X4
datDown <- melt(diffDown, id.vars = "geneID", measure.vars <- c(1:300))
datDown <- datDown[order(datDown$geneID),]
setkey(datDown, "geneID")

cols <- c('#ca0020','#f7f7f7','#0571b0')
cols <- c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')

rollWidth <- 10
rollBy <- 5

datDownSum <- lapply(unique(datDown$geneID), function(x){
  n <- rollapply(datDown[x]$value, width = rollWidth, by = rollBy, FUN = mean)
  df <- data.frame(geneID = x,
                   bin = 1:length(n),
                   value = n)
  return(df)
})
datDownSum <- as.data.table(as.data.frame(do.call("rbind", datDownSum)))
datDownSum$cat <- quantcut(datDownSum$value, q = 3)
levels(datDownSum$cat) <- c("loss", "no change", "gain")
#levels(datDownSum$cat) <- c("high loss", "loss", "no change", "gain", "high gain")

datUpSum <- lapply(unique(datUp$geneID), function(x){
  n <- rollapply(datUp[x]$value, width = rollWidth, by = rollBy, FUN = mean)
  df <- data.frame(geneID = x,
                   bin = 1:length(n),
                   value = n)
  return(df)
})
datUpSum <- as.data.table(as.data.frame(do.call("rbind", datUpSum)))
datUpSum$cat <- quantcut(datUpSum$value, q = 3)
levels(datUpSum$cat) <- c("loss", "no change", "gain")

ggplot(datUpSum, aes(bin, geneID)) + geom_tile(aes(fill = cat)) + scale_fill_manual(values = cols, levels(datUpSum$cat))
ggplot(datDownSum, aes(bin, geneID)) + geom_tile(aes(fill = cat)) + scale_fill_manual(values = cols, levels(datDownSum$cat))

# compare categories
setkey(datUpSum, geneID)
setkey(datDownSum, geneID)
catCountUp <- table(datUpSum[TGFbUpb2]$bin, datUpSum[TGFbUpb2]$cat)
catCountDown <- table(datDownSum[TGFbDownb2]$bin, datDownSum[TGFbDownb2]$cat)

# test for significance
setkey(datDownSum, bin)
setkey(datUpSum, bin)
lt <- lapply(1:max(datDownSum$bin), function(x){
  t.test(datUpSum[bin == x]$value, datDownSum[bin == x]$value)
})
unlist(lapply(lt, function(x) x$p.value))
plot(1:length(unlist(lapply(lt, function(x) x$p.value))), qvalue::qvalue(unlist(lapply(lt, function(x) x$p.value)))$qvalues, type = "l")

l1 <- lapply(1:nrow(catCountDown), function(i){
 t1 <- rbind(catCountDown[i,], catCountUp[i,])
 chisq.test(t1, simulate.p.value = T, B = 1000)
})
l2 <- lapply(1:nrow(catCountDown), function(i){
  t1 <- rbind(catCountDown[i,], catCountUp[i,])
  fisher.test(t1)
})
unlist(lapply(l1, function(x) x$p.value))
unlist(lapply(l1, function(x) x$statistic))
unlist(lapply(l2, function(x) x$p.value))

# multiple testing adjustment comparison
pdf("Compare_multiple_testing_adjustment.pdf")
par(mfrow = c(2,2))
hist(p.adjust(unlist(lapply(l1, function(x) x$p.value))))
hist(p.adjust(unlist(lapply(l1, function(x) x$p.value)), method = "fdr"))
hist(qvalue::qvalue(unlist(lapply(l1, function(x) x$p.value)))$qvalues)
hist(qvalue::lfdr(unlist(lapply(l1, function(x) x$p.value))))
dev.off()

pvals <- unlist(lapply(l1, function(x) x$p.value))
compAdjust <- data.frame(bh = p.adjust(pvals),
                   fdr = p.adjust(pvals, method = "fdr"),
                   qval = qvalue::qvalue(pvals)$qvalues,
                   lfdr = qvalue::lfdr(pvals))

compAdjust <- reshape::melt(compAdjust)
compAdjust$bin <- rep(1:length(pvals), length(unique(compAdjust$variable)))
compAdjust <- as.data.table(compAdjust)
setkey(compAdjust, variable)
ggplot(compAdjust[c("qval")], aes(bin, value)) + facet_wrap("variable") + geom_line()

# correlation tests -------------------------------------------------------
TGFbUpb2 <- deTGFbTab[cmTGFbUp$X4][deTGFbTab[cmTGFbUp$X4]$b > 0]$target_id
TGFbDownb2 <- deTGFbTab[cmTGFbDown$X4][deTGFbTab[cmTGFbDown$X4]$b < 0]$target_id
minDiff <- 10

datSum <- rbind(datDownSum, datUpSum)
colnames(datSum)[1] <- "ensembl_gene_id"
setkey(datSum, ensembl_gene_id, bin)
datSum <- datSum[order(datSum$ensembl_gene_id, datSum$bin)]
diffExp <- subset(deTGFbTab, target_id %in% c(TGFbDownb2,TGFbUpb2))
diffExp <- diffExp[order(diffExp$target_id)]
diffExp$logFC <- log2(exp(diffExp$b))
diffH2AZDEcor <- lapply(1:max(datDownSum$bin), function(x){
  allDat <- merge(subset(datSum, bin == x), diffExp, by.x = "ensembl_gene_id", by.y = "target_id")
  allDat <- merge(allDat, cfamEnsGenesSigEMTCells[,c("ensembl_gene_id", "external_gene_name", "epi_mes")], all.x = T)
  abLine <- coef(lm(logFC ~ value, data = allDat))
  allPlot <- ggplot(allDat, aes(x = value, y = logFC)) + 
                    geom_point() +
                    geom_abline(slope = abLine[2], intercept = abLine[1], colour = "darkgrey", size = 1, linetype = "longdash")
  allCor1 <- cor.test(allDat$value, allDat$logFC, method = "spearman")
  allCor2 <- cor.test(allDat[which(abs(allDat$value) > minDiff),]$value, allDat[which(abs(allDat$value) > minDiff),]$logFC, method = "spearman")
  return(list(all1 = allCor1,
              all2 = allCor2,
              allPlot = allPlot,
              allDat = allDat))
})

allCor <- lapply(diffH2AZDEcor, function(x){ data.frame(x$all1$p.value, 
                                                        x$all1$estimate,
                                                        x$all2$p.value,
                                                        x$all2$estimate)})
allCor <- do.call("rbind", allCor)
allCor$bin <- 1:nrow(allCor)
minCor <- which.min(allCor$x.all1.estimate) ; allCor[which.min(allCor$x.all1.estimate),]
maxCor <- which.max(allCor$x.all1.estimate)

# plot of Spearman correlation between differential expression and --------
corPlot <- ggplot(allCor, aes(bin, x.all1.estimate)) + 
  geom_line(col = "black", size = lineSize) +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  geom_text_repel(data = allCor[c(maxCor, minCor),], aes(x = bin, y = x.all1.estimate, label = c("Panel B", "Panel C")), label.padding = 1.25) +
  geom_point(data = allCor[c(which.max(allCor$x.all1.estimate), which.min(allCor$x.all1.estimate)),], aes(x = bin, y = x.all1.estimate), colour = "red") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  geom_vline(xintercept = 30, colour = vlineCol, linetype = "longdash")

corPlot

# plot X-Y plot of diffChIP vs DiffExp ------------------------------------

datMax <- diffH2AZDEcor[[maxCor]]$allDat
datMin <- diffH2AZDEcor[[minCor]]$allDat
datMin <- diffH2AZDEcor[[41]]$allDat
abLineMax <- coef(lm(logFC ~ value, data = datMax))
abLineMin <- coef(lm(logFC ~ value, data = datMin))

plotMax <- ggplot(datMax, aes(x = value, y = logFC)) + 
                  geom_point(aes(colour = epi_mes)) +
                  geom_abline(slope = coef(lm(logFC ~ value, data = datMax))[2], intercept = coef(lm(logFC ~ value, data = datMax))[1], colour = "darkgrey", size = 0.6, linetype = "longdash") +
                  geom_text_repel(data = subset(datMax, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, label = external_gene_name, colour = epi_mes), show.legend = F) +
                  ylab(NULL) +
                  xlab(NULL) +
                  labs(color = "Gene/Marker type") + scale_colour_hue(labels = c("Epithelial", "H2A.Z", "Mesenchymal", "TGFB1")) +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_rect(fill = NA), 
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())

plotMin <- ggplot(datMin, aes(x = value, y = logFC)) + 
                  geom_point(aes(colour = epi_mes)) +
                  geom_abline(slope = coef(lm(logFC ~ value, data = datMin))[2], intercept = coef(lm(logFC ~ value, data = datMin))[1], colour = "darkgrey", size = 0.6, linetype = "longdash") +
                  geom_text_repel(data = subset(datMin, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, label = external_gene_name, colour = epi_mes), show.legend = F) +
                  ylab(NULL) +
                  xlab(NULL) +
                  labs(color = "Gene/Marker type") + scale_colour_hue(labels = c("Epithelial", "H2A.Z", "Mesenchymal", "TGFB1")) +
                  theme(legend.position = "none",
                        panel.background = element_blank(),
                        panel.border = element_rect(fill = NA), 
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        legend.background = element_rect(fill = NA))

# plotting ----------------------------------------------------------------
gpFontSize <- gpar(fontsize = 8)
xAx1 <- grid::textGrob("Distance from TSS [bp]", gp = gpFontSize)
xAx2 <- grid::textGrob("Difference ChIP TFGb-WT [RPKM]", gp = gpFontSize)
yAx1 <- grid::textGrob("Correlation\n[Spearman rho-statistic]", rot = 90, gp = gpFontSize)
yAx2 <- grid::textGrob("Log2 fold-change\n[estimated]", rot = 90, gp = gpFontSize)
emptyBox <- grid::rectGrob(gp=grid::gpar(fill="white", lty = 0))
legendPlot <- g_legend(plotMin)

lay <- rbind(c(1,2,2),
             c(3,4,4),
             c(5,6,7),
             c(8,9,9))
gs <- lapply(1:max(lay), function(ii) {
  grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=ii, alpha=0.5)), grid::textGrob(ii))
})
gridExtra::grid.arrange(grobs = gs, layout_matrix = lay, widths = c(1,10,10), heights = c(4,1,6,1))

# remove legend as we can recycle legend from Figure 2 (assemble in Illustrator
gs <- list(yAx1,
           corPlot,
           emptyBox,
           xAx1,
           yAx2,
           plotMax, 
           plotMin, 
           emptyBox,
           xAx2)

grob2 <- arrangeGrob(grobs = gs,
                     layout_matrix = lay,
                     widths = c(1,1,8,8), heights = c(1,4,4,4,1))
ggsave("Figure3_panel_AB_input_subtracted_with_sem_updated.pdf", grob2, height = 98, width = 196, units = "mm")

dev.off()

