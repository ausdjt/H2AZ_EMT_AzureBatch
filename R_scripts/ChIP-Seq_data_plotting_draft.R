
matrixNames <- list.files(dataDir, pattern = "*EMT*")
matrixNames <- matrixNames[grep("allSamples", matrixNames, invert = T)]

l1 <- lapply(matrixNames, function(x){
  print(x)
  runInfo <- jsonlite::fromJSON(gsub("@",
                                     "", 
                                     readLines(gzfile(paste(dataDir, x, sep = "/")),
                                               n = 1)))
  t1 <- readr::read_tsv(gzfile(paste(dataDir, x, sep = "/")),
                        comment = "@",
                        col_names = FALSE)
  t1 <- t1[, -c(1:3,5:6)]
  maxCol <- max(runInfo$sample_boundaries)
  m1 <- as.matrix(t1[,2:(maxCol+1)])
  m1 <- data.frame(bin = c(1:maxCol), 
                   value = apply(m1, 2, mean))
  m1$data <- unlist(strsplit(x, "\\."))[1]
  return(m1)
})
names(l1) <- unlist(lapply(strsplit(matrixNames, "\\."), function (x) x[1]))

plotData <- as.data.frame(do.call("rbind", l1))

plotDataMNase <- plotData[grep("MNase", plotData$data),]
plotDataNormal <- plotData[grep("normal", plotData$data),]

ggplot(plotDataNormal[grep("H2AZ", plotDataNormal$data),], aes(bin, value)) + geom_line() + facet_wrap(~data, ncol = 2)

# plotting of differences between WT & TGFb treated
pdf("MDCK_Tan_et_al_EMT_genes_TSS.pdf")
par(mfrow = c(3,2))
par(mar = c(1.8,1.8,1,0.1))
plot(1:300, l1[["H2AZ-WT_TanEMTdown_normal"]]$value, main = "H2AZ-WT_TanEMTdown_normal", type = "l") #, ylim = c(0,80))
plot(1:300, l1[["H2AZ-WT_TanEMTup_normal"]]$value, main = "H2AZ-WT_TanEMTup_normal", type = "l", ylim = c(0,80))
plot(1:300, l1[["H2AZ-TGFb_TanEMTdown_normal"]]$value, main = "H2AZ-TGFb_TanEMTdown_normal", type = "l", ylim = c(0,80))
plot(1:300, l1[["H2AZ-TGFb_TanEMTup_normal"]]$value, main = "H2AZ-TGFb_TanEMTup_normal", type = "l", ylim = c(0,80))
plot(1:300, l1[["H2AZ-TGFb_TanEMTdown_normal"]]$value - l1[["H2AZ-WT_TanEMTdown_normal"]]$value, main = "TGFb - WT", type = "l", ylim = c(-40,40))
plot(1:300, l1[["H2AZ-TGFb_TanEMTup_normal"]]$value - l1[["H2AZ-WT_TanEMTup_normal"]]$value, main = "TGFb - WT", type = "l", ylim = c(-40,40))
dev.off()

# for MCF10A Data
dataDir <- "~/Data/Tremethick/Breast/ChIP-Seq/merged/processed_data/GRCh37_hg19_ensembl75_ERCC/deepTools/computeMatrix/reference-point/duplicates_removed/TSS/"
matrixNames <- list.files(dataDir, pattern = "*EMT*")

pdf("MDCK_Tan_et_al_EMT_genes_TSS.pdf")
par(mfrow = c(3,2))
par(mar = c(1.8,1.8,1,0.1))
plot(1:300, l1[["H2AZ_10A_high_TanEMTdown_normal"]]$value, main = "H2AZ-WT_TanEMTdown_normal", type = "l", ylim = c(0,40))
plot(1:300, l1[["H2AZ_10A_high_TanEMTup_normal"]]$value, main = "H2AZ-WT_TanEMTup_normal", type = "l", ylim = c(0,40))
plot(1:300, l1[["H2AZ_TGFb_10A_TanEMTdown_normal"]]$value, main = "H2AZ-TGFb_TanEMTdown_normal", type = "l", ylim = c(0,40))
plot(1:300, l1[["H2AZ_TGFb_10A_TanEMTup_normal"]]$value, main = "H2AZ-TGFb_TanEMTup_normal", type = "l", ylim = c(0,40))
plot(1:300, l1[["H2AZ_TGFb_10A_TanEMTdown_normal"]]$value - l1[["H2AZ_10A_high_TanEMTdown_normal"]]$value, main = "TGFb - WT", type = "l", ylim = c(-40,40))
plot(1:300, l1[["H2AZ_TGFb_10A_TanEMTup_normal"]]$value - l1[["H2AZ_10A_high_TanEMTup_normal"]]$value, main = "TGFb - WT", type = "l", ylim = c(-40,40))
dev.off()


# load matrix with data for all genes for further exploration
matrixNames <- list.files(dataDir, pattern = "*allGenes*")
matrixNames <- matrixNames[grep("allSamples", matrixNames, invert = T)]
matrixNames <- matrixNames[-1]

l1 <- lapply(matrixNames, function(x){
  print(x)
  runInfo <- jsonlite::fromJSON(gsub("@",
                                     "", 
                                     readLines(gzfile(paste(dataDir, x, sep = "/")),
                                               n = 1)))
  maxCol <- max(runInfo$sample_boundaries)
  t1 <- readr::read_tsv(gzfile(paste(dataDir, x, sep = "/")),
                        comment = "@",
                        col_names = FALSE)
  t1 <- t1[, -c(1:3,5:6)]
  return(t1)
})

names(l1) <- unlist(lapply(strsplit(matrixNames, "\\."), function (x) x[1]))


allGenesPlotData <- lapply(names(l1), function(x){
  m1 <- as.matrix(l1[[x]][,2:301])
  m1 <- data.frame(bin = c(1:300), 
                   value = apply(m1, 2, mean))
  m1$data <- unlist(strsplit(x, "\\."))[1]
  return(m1)
})
allGenesPlotData <- as.data.frame(do.call("rbind", allGenesPlotData))

ggplot2::ggplot(allGenesPlotData, aes(bin, value)) + geom_line() + facet_wrap(~data)
