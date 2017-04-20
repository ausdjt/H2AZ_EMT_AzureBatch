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

TGFbUpb2 <- deTGFbTab[cmTGFbUp$X4][deTGFbTab[cmTGFbUp$X4]$b > 0]$target_id
TGFbDownb2 <- deTGFbTab[cmTGFbDown$X4][deTGFbTab[cmTGFbDown$X4]$b < 0]$target_id
minDiff <- 15

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
  abLine <- coef(lm(b ~ value, data = allDat))
  allPlot <- ggplot(allDat, aes(x = value, y = b)) + 
    geom_point() +
    geom_abline(slope = abLine[2], intercept = abLine[1], colour = "darkgrey", size = 1, linetype = "longdash")
  allCor1 <- cor.test(allDat$value, allDat$b, method = "spearman")
  allCor2 <- cor.test(allDat[which(abs(allDat$value) > minDiff),]$value, allDat[which(abs(allDat$value) > minDiff),]$b, method = "spearman")
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
minCor1 <- which.min(allCor$x.all1.estimate) ; allCor[which.min(allCor$x.all1.estimate),]
maxCor1 <- which.max(allCor$x.all1.estimate) ; allCor[which.max(allCor$x.all1.estimate),]
minCor2 <- which.min(allCor$x.all2.estimate) ; allCor[which.min(allCor$x.all2.estimate),]
maxCor2 <- which.max(allCor$x.all2.estimate) ; allCor[which.max(allCor$x.all2.estimate),]

# plot of Spearman correlation between differential expression and --------
corPlot1 <- ggplot(allCor, aes(bin, x.all1.estimate)) + 
  geom_line(col = "darkgrey", size = lineSize) +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  geom_text_repel(data = allCor[c(maxCor1, minCor1),], aes(x = bin, y = x.all2.estimate, label =c("-2 Nucleosome", "+2 Nucleosome")), nudge_y = c(0.05, 0.1)) +
  geom_point(data = allCor[c(maxCor1, minCor1),], aes(x = bin, y = x.all1.estimate), colour = "red") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  geom_vline(xintercept = 30, colour = vlineCol, linetype = "longdash")
corPlot1

# with genes that have min. +/- 10 H2A.Z difference in a given bin
corPlot2 <- ggplot(allCor, aes(bin, x.all2.estimate)) + 
  geom_line(col = "darkgrey", size = lineSize) +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  geom_text_repel(data = allCor[c(maxCor2, minCor2),], aes(x = bin, y = x.all2.estimate, label =c("-2", "+1")), nudge_y = c(-0.15, +0.1), size = 3) +
  geom_point(data = allCor[c(maxCor2, minCor2),], aes(x = bin, y = x.all2.estimate), colour = "red") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", axisLineSize)) +
  xlab(NULL) +
  ylab(NULL) +
  geom_vline(xintercept = 30, colour = vlineCol, linetype = "longdash")

corPlot2

# plot X-Y plot of diffChIP vs DiffExp ------------------------------------

datMax <- diffH2AZDEcor[[maxCor2]]$allDat
datMin <- diffH2AZDEcor[[minCor2]]$allDat
abLineMax <- coef(lm(b ~ value, data = datMax))
abLineMin <- coef(lm(b ~ value, data = datMin))

plotMax <- ggplot(datMax, aes(x = value, y = logFC)) + 
  geom_point(aes(colour = epi_mes), size = 0.05) +
  geom_abline(slope = coef(lm(logFC ~ value, data = datMax))[2], intercept = coef(lm(logFC ~ value, data = datMax))[1], colour = "darkgrey", size = 0.6, linetype = "longdash") +
  geom_point(data = subset(datMax, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, colour = epi_mes), size = 1) +
  geom_text_repel(data = subset(datMax, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, label = external_gene_name, colour = epi_mes), show.legend = F, size = 2) +
  ylab(NULL) +
  xlab(NULL) +
  labs(color = "Gene/Marker type") + scale_colour_hue(labels = c("Epithelial", "H2A.Z", "Mesenchymal", "TGFB1")) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", axisLineSize),
        axis.text.y = element_text(margin = margin(0,3,0,0, "mm")))
plotMax

plotMin <- ggplot(datMin, aes(x = value, y = logFC)) + 
  geom_point(aes(colour = epi_mes), size = 0.05) +
  geom_abline(slope = coef(lm(logFC ~ value, data = datMin))[2], intercept = coef(lm(logFC ~ value, data = datMin))[1], colour = "darkgrey", size = 0.6, linetype = "longdash") +
  geom_point(data = subset(datMin, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, colour = epi_mes), size = 1) +
  geom_text_repel(data = subset(datMin, value < -200 | abs(logFC) > 3), aes(x = value, y = logFC, label = external_gene_name, colour = epi_mes), show.legend = F, size = 2) +
  ylab(NULL) +
  xlab(NULL) +
  labs(color = "Gene/Marker type") + scale_colour_hue(labels = c("Epithelial", "H2A.Z", "Mesenchymal", "TGFB1")) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = NA),
        axis.line = element_line(colour = "black", axisLineSize))
plotMin
# plotting ----------------------------------------------------------------
gpFontSize <- gpar(fontsize = 8)
xAx1 <- grid::textGrob("Distance from TSS [bp]", gp = gpFontSize)
xAx2 <- grid::textGrob("TFGb-WT [RPKM]", gp = gpFontSize)
yAx1 <- grid::textGrob("Spearman\n[rho-statistic]", rot = 90, gp = gpFontSize)
yAx2 <- grid::textGrob("Log2 fold-change\n[estimated]", rot = 90, gp = gpFontSize)
emptyBox <- grid::rectGrob(gp=grid::gpar(fill="white", lty = 0))
legendPlot <- g_legend(plotMin)

lay <- rbind(c(1,2,2),
             c(3,4,4),
             c(5,6,7),
             c(8,9,9))
# remove legend as we can recycle legend from Figure 2 (assemble in Illustrator
gs <- list(yAx1,
           corPlot2,
           emptyBox,
           xAx1,
           yAx2,
           plotMax, 
           plotMin, 
           emptyBox,
           xAx2)
grob2 <- arrangeGrob(grobs = gs,
                     layout_matrix = lay,
                     widths = c(1,10,10), heights = c(2.5,0.5,7.5,0.5))
ggsave("Figure3_panel_CDE.pdf", grob2, height = 108, width = 196, units = "mm", useDingbats = F)


gridExtra::grid.arrange(grobs = gs, layout_matrix = lay, widths = c(1,10,10), heights = c(4,1,6,1))
gs <- lapply(1:max(lay), function(ii) {
  grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=ii, alpha=0.5)), grid::textGrob(ii))
})
gridExtra::grid.arrange(grobs = gs, layout_matrix = lay, widths = c(1,10,10), heights = c(4,1,6,1))

