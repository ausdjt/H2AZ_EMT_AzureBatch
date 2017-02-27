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

# make difference data by subtracting WT from TGFb
bin <- 1:300
diff1 <- data.frame(bin = bin,
                    value = plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$value 
                    - plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal", ]$value,
                    sample = "TGFb - WT",
                    group = "down_diff",
                    geneset = "TanEMTdown")
diff2 <- data.frame(bin = bin,
                    value = plotDataTGFb[plotDataTGFb$geneset == "H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$value 
                    - plotDataWT[plotDataWT$geneset == "H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal", ]$value,
                    sample = "TGFb - WT",
                    group = "down_diff",
                    geneset = "TanEMTup")
diffData <- rbind(diff1, diff2)

# publication plot --------------------------------------------------------
#pdf("TSS_coverage_Tan_EMT_genes_with_diff.pdf")
#dev.off()
lineSize <- 2
vlineCol <- "black"

wtPlot <- ggplot(plotDataWT[grep("WT", plotDataWT$sample),], aes(bin, value)) 
#wtPlot <- wtPlot + geom_line(size = lineSize, col = "darkgreen")
wtPlot <- wtPlot + geom_smooth(method = "lm",
                               formula = y ~ splines::ns(x, length(plotDataWT$bin)/25), 
                               se = T, 
                               col = "darkgreen",
                               size = lineSize)
wtPlot <- wtPlot + facet_wrap(c("geneset"), ncol = 2, labeller = label_both)
wtPlot <- wtPlot + ylab(NULL)
wtPlot <- wtPlot + xlab(NULL)
wtPlot <- wtPlot + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
wtPlot <- wtPlot + ylim(0,75)
wtPlot <- wtPlot + theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         axis.text.y= element_text(margin = margin(0,4,0,0, "mm"), size = 14),
                         strip.background = element_blank(), 
                         strip.text = element_blank(),
                         axis.line.y = element_line(colour = "black", size = 1),
                         panel.background = element_blank(),
                         panel.border = element_blank(), 
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())

tgfbPlot <- ggplot(plotDataTGFb[grep("TGFb", plotDataTGFb$sample),], aes(bin, value)) 
# tgfbPlot <- tgfbPlot + geom_line(size = lineSize, col = "darkblue")
tgfbPlot <- tgfbPlot + geom_smooth(method = "lm", 
                                   formula = y ~ splines::ns(x, length(plotDataTGFb$bin)/25),
                                   se = T, 
                                   col = "darkblue",
                                   size = lineSize)
tgfbPlot <- tgfbPlot + facet_wrap(c("geneset"), ncol = 2, labeller = label_both)
tgfbPlot <- tgfbPlot + ylab(NULL)
tgfbPlot <- tgfbPlot + xlab(NULL)
tgfbPlot <- tgfbPlot + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
tgfbPlot <- tgfbPlot + ylim(0,75)
tgfbPlot <- tgfbPlot + theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y= element_text(margin = margin(0,4,0,0, "mm"), size = 14),
                             strip.background = element_blank(), 
                             strip.text = element_blank(),
                             axis.line.y = element_line(colour = "black", size = 1),
                             panel.background = element_blank(),
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())

diffPlot <- ggplot(diffData, aes(bin, value)) 
# diffPlot <- diffPlot + geom_line(size = lineSize, col = "darkgrey")
diffPlot <- diffPlot + geom_smooth(method = "lm", 
                                   formula = y ~ splines::ns(x, length(diffData$bin)/25), 
                                   se = T, 
                                   col = "darkgrey", 
                                   size = lineSize)
diffPlot <- diffPlot + facet_wrap(c("geneset"), ncol = 2)
diffPlot <- diffPlot + ylab(NULL)
diffPlot <- diffPlot + xlab(NULL)
diffPlot <- diffPlot + geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash")
diffPlot <- diffPlot + scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                                          labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500"))
diffPlot <- diffPlot + theme(strip.background = element_blank(), strip.text = element_blank(),
                             axis.line.y = element_line(colour = "black", size = 1),
                             axis.line.x = element_line(colour = "black", size = 1),
                             axis.text.y = element_text(margin = margin(0,1.1,0,0, "mm"), size = 14),
                             axis.text.x = element_text(size = 14),
                             panel.background = element_blank(),
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())


# prepare layout grid for panels
lay <- rbind(c(1,2,3,4),
             c(1,5,6,6),
             c(1,7,8,8),
             c(1,9,10,10),
             c(11,11,12,12))

xAx <- grid::textGrob("Distance from TSS [bp]")
yAx <- grid::textGrob("Mean coverage (Input subtracted) [RPKM]", rot = 90)
topX1 <- grid::textGrob("Down-regulated Genes")
topX2 <- grid::textGrob("Up-regulated Genes")
wtLabel <- grid::textGrob("MDCK WT", rot = 90)
tgfbLabel <- grid::textGrob("TGFb-treated", rot = 90)
diffLabel <- grid::textGrob("Diff [TGFb-WT]", rot = 90)
emptyBox <- grid::rectGrob(gp=grid::gpar(fill="white", lty = 0))

pdf("FigureX_panel_A_input_subtracted.pdf", height = 7, width = 14)
gridExtra::grid.arrange(yAx, 
                        emptyBox, 
                        topX1, 
                        topX2, 
                        wtLabel,
                        wtPlot,
                        tgfbLabel,
                        tgfbPlot,
                        diffLabel,
                        diffPlot,
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
