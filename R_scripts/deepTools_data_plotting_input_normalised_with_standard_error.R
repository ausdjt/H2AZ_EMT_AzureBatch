require(readr)
require(jsonlite)
require(GenomicRanges)
require(data.table)
require(ggplot2)
require(deepToolsUtils)

baseDir <- "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC"
dataDir <- paste(baseDir, 
                 "deepTools/computeMatrix/reference-point/bigwigCompare/duplicates_marked/TSS/",
                 sep = "/")
setwd("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/")
list.files(dataDir)

# plot the TSS of 100 up & 100 down-regulated randomly selected ge --------
files <- list.files(dataDir, pattern = "Tan")
files <- files[grep("normal", files)]
files <- files[grep("MNase", files, invert = T)]
files <- files[grep("readCount.subtract", files)]

dataList <- lapply(paste(dataDir, files, sep = "/"), function(x) computeMatrixLoader(x))
plotData <- lapply(dataList, function(x) deepToolsUtils::makePlottingData(x))
names(plotData) <- unlist(lapply(strsplit(files, "\\."), function(x) paste(x[c(1:3)], collapse = ".")))
plotData <- lapply(names(plotData), function(x){
  tab <- plotData[[x]]
  tab$geneset <- x
  return(tab)
})
names(plotData) <- unlist(lapply(strsplit(files, "\\."), function(x) paste(x[c(1:3)], collapse = ".")))
plotData <- do.call("rbind", plotData)
plotDataWT <- plotData[grep("WT", plotData$sample),]
plotDataTGFb <- plotData[grep("TGFb", plotData$sample),]
table(plotData$geneset)
ylimMax <- max(plotData$value + plotData$sem)
ylimMin <- min(plotData$value - plotData$sem)
# make difference data by subtracting WT from TGFb
bin <- 1:300
diff1 <- data.frame(bin = bin,
                    value = plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$value 
                    - plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$value,
                    sample = "TGFb - WT",
                    group = "down_diff",
                    geneset = "TanEMTdown",
                    sem = sqrt((plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$sem)^2 
                              + (plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$sem)^2))
diff2 <- data.frame(bin = bin,
                    value = plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$value 
                    - plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$value,
                    sample = "TGFb - WT",
                    group = "down_diff",
                    geneset = "TanEMTup",
                    sem = sqrt((plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$sem)^2 
                              + (plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$sem)^2))
diffData <- rbind(diff1, diff2)
diffylimMin <- min(diffData$value - diffData$sem)
diffylimMax <- max(diffData$value + diffData$sem)

# publication plot --------------------------------------------------------
#pdf("TSS_coverage_Tan_EMT_genes_with_diff.pdf")
#dev.off()
lineSize <- 2
vlineCol <- "black"
wtCol <- "darkgreen"
TGFbCol <- "darkblue"
diffCol <-"darkgrey"
nSplines <- 30
downXAxisMargin <- 4.8

# plot 1 - WT down-regulated ----------------------------------------------
wtPlotDown <- ggplot(plotDataWT[grep("down", plotDataWT$geneset),], aes(bin, value)) 
wtPlotDown <- wtPlotDown + geom_smooth(method = "lm",
                               formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                               se = T, 
                               col = wtCol,
                               size = lineSize)
wtPlotDown <- wtPlotDown + geom_smooth(aes(x = bin, y = value - sem),
                               method = "lm",
                               formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                               se = F, 
                               col = alpha(wtCol, 0.2),
                               size = 0.5,
                               fullrange = F)
wtPlotDown <- wtPlotDown + geom_smooth(aes(x = bin, y = value + sem),
                               method = "lm",
                               formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                               se = F, 
                               col = alpha(wtCol, 0.2),
                               size = 0.5,
                               fullrange = F)
wtPlotDown <- wtPlotDown + facet_wrap(c("geneset"), ncol = 2, labeller = label_both)
wtPlotDown <- wtPlotDown + ylab(NULL)
wtPlotDown <- wtPlotDown + xlab(NULL)
wtPlotDown <- wtPlotDown + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
wtPlotDown <- wtPlotDown + ylim(0,ylimMax)
wtPlotDown <- wtPlotDown + theme(axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y= element_text(margin = margin(0,downXAxisMargin,0,0, "mm"), size = 14),
                                 strip.background = element_blank(), 
                                 strip.text = element_blank(),
                                 axis.line.y = element_line(colour = "black", size = 1),
                                 panel.background = element_blank(),
                                 panel.border = element_blank(), 
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
wtPlotDownData <- ggplot_build(wtPlotDown)
df1 <- data.frame(x = wtPlotDownData$data[[2]]$x,
                  ymin = wtPlotDownData$data[[2]]$y,
                  ymax = wtPlotDownData$data[[3]]$y)
wtPlotDown <- wtPlotDown + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(wtCol, 0.3), alpha = 0.1, inherit.aes = F, fill = wtCol)


# plot 2 - WT up-regulated ------------------------------------------------
wtPlotUp <- ggplot(plotDataWT[grep("up", plotDataWT$geneset),], aes(bin, value)) 
wtPlotUp <- wtPlotUp + geom_smooth(method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                                       se = T, 
                                       col = wtCol,
                                       size = lineSize)
wtPlotUp <- wtPlotUp + geom_smooth(aes(x = bin, y = value - sem),
                                       method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                                       se = F, 
                                       col = alpha(wtCol, 0.2),
                                       size = 0.5,
                                       fullrange = F)
wtPlotUp <- wtPlotUp + geom_smooth(aes(x = bin, y = value + sem),
                                       method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataWT$bin)/nSplines), 
                                       se = F, 
                                       col = alpha(wtCol, 0.2),
                                       size = 0.5,
                                       fullrange = F)
wtPlotUp <- wtPlotUp + ylab(NULL)
wtPlotUp <- wtPlotUp + xlab(NULL)
wtPlotUp <- wtPlotUp + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
wtPlotUp <- wtPlotUp + ylim(0,ylimMax)
wtPlotUp <- wtPlotUp + theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.line.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             strip.background = element_blank(), 
                             strip.text = element_blank(),
                             panel.background = element_blank(),
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())
wtPlotUpData <- ggplot_build(wtPlotUp)
df1 <- data.frame(x = wtPlotUpData$data[[2]]$x,
                  ymin = wtPlotUpData$data[[2]]$y,
                  ymax = wtPlotUpData$data[[3]]$y)
wtPlotUp <- wtPlotUp + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(wtCol, 0.3), alpha = 0.1, inherit.aes = F, fill = wtCol)


# plot 3 - TGFb down-regulated --------------------------------------------
TGFbPlotDown <- ggplot(plotDataTGFb[grep("down", plotDataTGFb$geneset),], aes(bin, value)) 
TGFbPlotDown <- TGFbPlotDown + geom_smooth(method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                       se = T, 
                                       col = TGFbCol,
                                       size = lineSize)
TGFbPlotDown <- TGFbPlotDown + geom_smooth(aes(x = bin, y = value - sem),
                                       method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                       se = F, 
                                       col = alpha(TGFbCol, 0.2),
                                       size = 0.5,
                                       fullrange = F)
TGFbPlotDown <- TGFbPlotDown + geom_smooth(aes(x = bin, y = value + sem),
                                       method = "lm",
                                       formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                       se = F, 
                                       col = alpha(TGFbCol, 0.2),
                                       size = 0.5,
                                       fullrange = F)
TGFbPlotDown <- TGFbPlotDown + ylab(NULL)
TGFbPlotDown <- TGFbPlotDown + xlab(NULL)
TGFbPlotDown <- TGFbPlotDown + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
TGFbPlotDown <- TGFbPlotDown + ylim(0,ylimMax)
TGFbPlotDown <- TGFbPlotDown + theme(axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 axis.text.y= element_text(margin = margin(0,downXAxisMargin,0,0, "mm"), size = 14),
                                 strip.background = element_blank(), 
                                 strip.text = element_blank(),
                                 axis.line.y = element_line(colour = "black", size = 1),
                                 panel.background = element_blank(),
                                 panel.border = element_blank(), 
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
TGFbPlotDownData <- ggplot_build(TGFbPlotDown)
df1 <- data.frame(x = TGFbPlotDownData$data[[2]]$x,
                  ymin = TGFbPlotDownData$data[[2]]$y,
                  ymax = TGFbPlotDownData$data[[3]]$y)
TGFbPlotDown <- TGFbPlotDown + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(TGFbCol, 0.3), alpha = 0.1, inherit.aes = F, fill = TGFbCol)


# plot 4 - TGFb up-regulated ----------------------------------------------
TGFbPlotUp <- ggplot(plotDataTGFb[grep("up", plotDataTGFb$geneset),], aes(bin, value)) 
TGFbPlotUp <- TGFbPlotUp + geom_smooth(method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = T, 
                                           col = TGFbCol,
                                           size = lineSize)
TGFbPlotUp <- TGFbPlotUp + geom_smooth(aes(x = bin, y = value - sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(TGFbCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
TGFbPlotUp <- TGFbPlotUp + geom_smooth(aes(x = bin, y = value + sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(TGFbCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
TGFbPlotUp <- TGFbPlotUp + ylab(NULL)
TGFbPlotUp <- TGFbPlotUp + xlab(NULL)
TGFbPlotUp <- TGFbPlotUp + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
TGFbPlotUp <- TGFbPlotUp + ylim(0,ylimMax)
TGFbPlotUp <- TGFbPlotUp + theme(axis.title.x = element_blank(),
                                 axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 axis.text.y= element_blank(),
                                 axis.line.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 strip.background = element_blank(), 
                                 strip.text = element_blank(),
                                 panel.background = element_blank(),
                                 panel.border = element_blank(), 
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
TGFbPlotUpData <- ggplot_build(TGFbPlotUp)
df1 <- data.frame(x = TGFbPlotUpData$data[[2]]$x,
                  ymin = TGFbPlotUpData$data[[2]]$y,
                  ymax = TGFbPlotUpData$data[[3]]$y)
TGFbPlotUp <- TGFbPlotUp + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(TGFbCol, 0.3), alpha = 0.1, inherit.aes = F, fill = TGFbCol)


# plot 5 - diff down-regulated --------------------------------------------
diffPlotDown <- ggplot(diffData[grep("down", diffData$geneset),], aes(bin, value)) 
diffPlotDown <- diffPlotDown + geom_smooth(method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = T, 
                                           col = diffCol,
                                           size = lineSize)
diffPlotDown <- diffPlotDown + geom_smooth(aes(x = bin, y = value - sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(diffCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
diffPlotDown <- diffPlotDown + geom_smooth(aes(x = bin, y = value + sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(diffCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
diffPlotDown <- diffPlotDown + ylab(NULL)
diffPlotDown <- diffPlotDown + xlab(NULL)
diffPlotDown <- diffPlotDown + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
diffPlotDown <- diffPlotDown + ylim(diffylimMin, diffylimMax)
diffPlotDown <- diffPlotDown + scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                                          labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500"))
diffPlotDown <- diffPlotDown + theme(axis.line.x = element_line(colour = "black", size = 1),
                                     axis.text.x = element_text(size = 14),
                                     axis.text.y= element_text(margin = margin(0,2,0,0, "mm"), size = 14),
                                     strip.background = element_blank(), 
                                     strip.text = element_blank(),
                                     axis.line.y = element_line(colour = "black", size = 1),
                                     panel.background = element_blank(),
                                     panel.border = element_blank(), 
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank())
diffPlotDownData <- ggplot_build(diffPlotDown)
df1 <- data.frame(x = diffPlotDownData$data[[2]]$x,
                  ymin = diffPlotDownData$data[[2]]$y,
                  ymax = diffPlotDownData$data[[3]]$y)
diffPlotDown <- diffPlotDown + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(diffCol, 0.3), alpha = 0.1, inherit.aes = F, fill = diffCol)

# plot 6 - diff up-regulated --------------------------------------------
diffPlotUp <- ggplot(diffData[grep("up", diffData$geneset),], aes(bin, value)) 
diffPlotUp <- diffPlotUp + geom_smooth(method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = T, 
                                           col = diffCol,
                                           size = lineSize)
diffPlotUp <- diffPlotUp + geom_smooth(aes(x = bin, y = value - sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(diffCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
diffPlotUp <- diffPlotUp + geom_smooth(aes(x = bin, y = value + sem),
                                           method = "lm",
                                           formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/nSplines), 
                                           se = F, 
                                           col = alpha(diffCol, 0.2),
                                           size = 0.5,
                                           fullrange = F)
diffPlotUp <- diffPlotUp + ylab(NULL)
diffPlotUp <- diffPlotUp + xlab(NULL)
diffPlotUp <- diffPlotUp + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
diffPlotUp <- diffPlotUp + scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                                              labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500"))
diffPlotUp <- diffPlotUp + ylim(diffylimMin, diffylimMax)
diffPlotUp <- diffPlotUp + theme(axis.text.y= element_blank(),
                                 axis.line.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.line.x = element_line(colour = "black", size = 1),
                                 axis.text.x = element_text(size = 14),
                                 strip.background = element_blank(), 
                                 strip.text = element_blank(),
                                 panel.background = element_blank(),
                                 panel.border = element_blank(), 
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
diffPlotUpData <- ggplot_build(diffPlotUp)
df1 <- data.frame(x = diffPlotUpData$data[[2]]$x,
                  ymin = diffPlotUpData$data[[2]]$y,
                  ymax = diffPlotUpData$data[[3]]$y)
diffPlotUp <- diffPlotUp + geom_ribbon(data = df1, aes(x = x, ymin = ymin, ymax = ymax), colour = alpha(diffCol, 0.3), alpha = 0.1, inherit.aes = F, fill = diffCol)

# prepare layout grid for panels
lay <- rbind(c(1,2,3,4),
            c(1,5,6,7),
            c(1,8,9,10),
            c(1,11,12,13),
            c(14,14,15,15))

xAx <- grid::textGrob("Distance from TSS [bp]")
yAx <- grid::textGrob("Mean coverage (Input subtracted) [RPKM]", rot = 90)
topX1 <- grid::textGrob("Down-regulated Genes")
topX2 <- grid::textGrob("Up-regulated Genes")
wtLabel <- grid::textGrob("MDCK WT", rot = 90)
tgfbLabel <- grid::textGrob("TGFb-treated", rot = 90)
diffLabel <- grid::textGrob("Diff [TGFb-WT]", rot = 90)
emptyBox <- grid::rectGrob(gp=grid::gpar(fill="white", lty = 0))

pdf("FigureX_panel_A_input_subtracted_with_sem.pdf", height = 7, width = 14)
gridExtra::grid.arrange(yAx, 
                        emptyBox, 
                        topX1, 
                        topX2, 
                        wtLabel,
                        wtPlotDown,
                        wtPlotUp,
                        tgfbLabel,
                        TGFbPlotDown,
                        TGFbPlotUp,
                        diffLabel,
                        diffPlotDown,
                        diffPlotUp,
                        emptyBox,
                        xAx,
                        layout_matrix = lay, 
                        widths = c(1,1,8,8), heights = c(1,4,4,4,1))
dev.off()
# for checking layout grid
gs <- lapply(1:12, function(ii) {
  grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=ii, alpha=0.5)), grid::textGrob(ii))
})
gridExtra::grid.arrange(grobs = gs, layout_matrix = lay, widths = c(1,1,8,8), heights = c(1,4,4,4,1))
