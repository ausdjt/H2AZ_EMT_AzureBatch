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
names(dataList) <- files
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
ylimMax <- max(plotData$value + (plotData$sem/2))
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

# publication plot --------------------------------------------------------
#pdf("TSS_coverage_Tan_EMT_genes_with_diff.pdf")
#dev.off()
diffylimMin <- min(diffData$value - (diffData$sem/2.3))
diffylimMax <- max(diffData$value + (diffData$sem/1.3))
lineSize <- 2
axisLineSize <- 0.75
vlineCol <- "black"
wtCol <- "#1b9e77"
TGFbCol <- "#7570b3"
diffCol <-"darkgrey"
nSplines <- 20
downXAxisMargin <- 3.6
axisTextSize <- 8

# plot 1 - WT down-regulated ----------------------------------------------
nBins <- nrow(plotDataWT[grep("down", plotDataWT$geneset),])
wtPlotDown <- ggplot(plotDataWT[grep("down", plotDataWT$geneset),], aes(bin, value)) + 
  geom_smooth(method = "lm",
             formula = y ~ splines::ns(x, nBins/nSplines), 
             se = F, 
             col = wtCol,
             size = lineSize) + 
  geom_smooth(aes(x = bin, y = value - sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(wtCol, 0.2),
              size = 0.5,
              fullrange = F) +
  geom_smooth(aes(x = bin, y = value + sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(wtCol, 0.2),
              size = 0.5,
              fullrange = F) +
  facet_wrap(c("geneset"), ncol = 2, labeller = label_both) +
  ylab(NULL) + 
  xlab(NULL) + 
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  coord_cartesian(ylim = c(0,ylimMax)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y= element_text(margin = margin(0,downXAxisMargin,0,0, "mm"), size = axisTextSize),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.line.y = element_line(colour = "black", axisLineSize),
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
wtPlotUp <- ggplot(plotDataWT[grep("up", plotDataWT$geneset),], aes(bin, value)) + 
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = wtCol,
              size = lineSize) + 
  geom_smooth(aes(x = bin, y = value - sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(wtCol, 0.2),
              size = 0.5,
              fullrange = F) + 
  geom_smooth(aes(x = bin, y = value + sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(wtCol, 0.2),
              size = 0.5,
              fullrange = F) +
  ylab(NULL) +
  xlab(NULL) + 
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  coord_cartesian(ylim = c(0,ylimMax)) + 
  theme(axis.title.x = element_blank(),
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
TGFbPlotDown <- ggplot(plotDataTGFb[grep("down", plotDataTGFb$geneset),], aes(bin, value)) + 
  geom_smooth(method = "lm",
             formula = y ~ splines::ns(x, nBins/nSplines), 
             se = F, 
             col = TGFbCol,
             size = lineSize) + 
  geom_smooth(aes(x = bin, y = value - sem),
             method = "lm",
             formula = y ~ splines::ns(x, nBins/nSplines), 
             se = F, 
             col = alpha(TGFbCol, 0.2),
             size = 0.5,
             fullrange = F) +
  geom_smooth(aes(x = bin, y = value + sem),
             method = "lm",
             formula = y ~ splines::ns(x, nBins/nSplines), 
             se = F, 
             col = alpha(TGFbCol, 0.2),
             size = 0.5,
             fullrange = F) +
  ylab(NULL) +
  xlab(NULL) + 
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  coord_cartesian(ylim = c(0,ylimMax)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y= element_text(margin = margin(0,downXAxisMargin,0,0, "mm"), size = axisTextSize),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.line.y = element_line(colour = "black", axisLineSize),
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
TGFbPlotUp <- ggplot(plotDataTGFb[grep("up", plotDataTGFb$geneset),], aes(bin, value))  +
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = TGFbCol,
              size = lineSize) +
 geom_smooth(aes(x = bin, y = value - sem),
             method = "lm",
             formula = y ~ splines::ns(x, nBins/nSplines), 
             se = F, 
             col = alpha(TGFbCol, 0.2),
             size = 0.5,
             fullrange = F) +
  geom_smooth(aes(x = bin, y = value + sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(TGFbCol, 0.2),
              size = 0.5,
              fullrange = F) +
  ylab(NULL) +
  xlab(NULL) + 
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  coord_cartesian(ylim = c(0,ylimMax)) + 
  theme(axis.title.x = element_blank(),
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
annotationSize <- 2
diffPlotDown <- ggplot(diffData[grep("down", diffData$geneset),], aes(bin, value))  +
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = diffCol,
              size = lineSize) +
  geom_smooth(aes(x = bin, y = value - sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(diffCol, 0.2),
              size = 0.5,
              fullrange = F) +
  geom_smooth(aes(x = bin, y = value + sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(diffCol, 0.2),
              size = 0.5,
              fullrange = F) +
  ylab(NULL) +
  xlab(NULL) +
  geom_text_repel(data = subset(diffData, geneset == "TanEMTdown")[118], aes(x = bin, y = value, label = "-1"), size = annotationSize, nudge_y = -1) +
  geom_text_repel(data = subset(diffData, geneset == "TanEMTdown")[68], aes(x = bin, y = value, label = "-2"), size = annotationSize, nudge_y = -1) +
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  coord_cartesian(ylim = c(diffylimMax, diffylimMin)) +
  scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  theme(axis.line.x = element_line(colour = "black", axisLineSize),
        axis.text.x = element_text(size = axisTextSize),
        axis.text.y= element_text(margin = margin(0,2,0,0, "mm"), size = axisTextSize),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.line.y = element_line(colour = "black", axisLineSize),
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
diffPlotUp <- ggplot(diffData[grep("up", diffData$geneset),], aes(bin, value))  +
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = diffCol,
             size = lineSize) +
  geom_smooth(aes(x = bin, y = value - sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(diffCol, 0.2),
              size = 0.5,
              fullrange = F) +
  geom_smooth(aes(x = bin, y = value + sem),
              method = "lm",
              formula = y ~ splines::ns(x, nBins/nSplines), 
              se = F, 
              col = alpha(diffCol, 0.2),
              size = 0.5,
              fullrange = F) +
  ylab(NULL) +
  xlab(NULL) +
  geom_text_repel(data = subset(diffData, geneset == "TanEMTup")[202], aes(x = bin, y = value, label = "+1"), size = annotationSize, nudge_y = -1) +
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  coord_cartesian(ylim = c(diffylimMax, diffylimMin)) + 
  theme(axis.text.y= element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "black", axisLineSize),
        axis.text.x = element_text(size = axisTextSize),
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

gpFontSize <- gpar(fontsize = 8)
xAx <- grid::textGrob("Distance from TSS [bp]", gp = gpFontSize)
yAx <- grid::textGrob("Mean coverage (Input subtracted) [RPKM]", rot = 90, gp = gpFontSize)
topX1 <- grid::textGrob("Down-regulated EMT genes", gp = gpFontSize)
topX2 <- grid::textGrob("Up-regulated EMT genes", gp = gpFontSize)
wtLabel <- grid::textGrob("WT", rot = 90, gp = gpFontSize)
tgfbLabel <- grid::textGrob("TGFb", rot = 90, gp = gpFontSize)
diffLabel <- grid::textGrob("TGFb-WT", rot = 90, gp = gpFontSize)
emptyBox <- grid::rectGrob(gp=grid::gpar(fill="white", lty = 0))

gs <- list(yAx, 
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
           xAx)
grob1 <- arrangeGrob(grobs = gs,
                     layout_matrix = lay,
                     widths = c(0.5,0.5,8,8), heights = c(1,4,4,4,1))
ggsave("Figure3_panel_AB_input_subtracted_with_sem_updated_20170317.pdf", grob1, height = 96, width = 196, units = "mm", useDingbats = F)

gridExtra::grid.arrange(grobs = gs,
                        layout_matrix = lay, 
                        widths = c(0.5,0.5,8,8), heights = c(1,4,4,4,1))
dev.off()
# for checking layout grid
gs <- lapply(1:max(lay), function(ii) {
  grid::grobTree(grid::rectGrob(gp=grid::gpar(fill=ii, alpha=0.5)), grid::textGrob(ii))
})
gridExtra::grid.arrange(grobs = gs, layout_matrix = lay, widths = c(0.5,0.5,8,8), heights = c(1,4,4,4,1))
