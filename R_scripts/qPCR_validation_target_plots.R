# look at expression levels of genes Renae checked via qPCR
require(reshape)

geneList <- c (H2AZ = "ENSCAFG00000010615", 
               ECadherin = "ENSCAFG00000020305", 
               NCadherin = "ENSCAFG00000018115",
               EpCAM = "ENSCAFG00000002653",
               Fibronectin = "ENSCAFG00000013142",
               TGFb1 = "ENSCAFG00000005014")

geneExp <- txi$abundance
gdata <- melt(geneExp[geneList,])
gdata$value <- log2(gdata$value)
p1 <- ggplot(gdata, aes(x = X1, y = value, fill= X2))
p1 + geom_boxplot()

pdf("qPCR_validation_targets_expression_levels.pdf")
sapply(names(geneList), function(x){
  dat <- geneExp[geneList[x], ]
  dat <- melt(dat)
  dat$sample <- rownames(dat)
  dat$value <- log2(dat$value + 1)
  p1 <- ggplot(dat, aes(x = sample, y = value, fill = sample)) + geom_bar(stat="identity") + ggtitle(paste(x)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot(p1)
})
dev.off()
