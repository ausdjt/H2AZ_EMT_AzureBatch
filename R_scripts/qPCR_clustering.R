require(cluster)
require(rtracklayer)
require(ade4)
 
setwd("/Users/u1001407/Data/Tremethick/EMT/R_analysis")

sd1 <- apply(as.matrix(exprs(dCT.norm2)), 1, sd)
m1 <- as.matrix(exprs(dCT.norm2))[which(sd1 > 0.5), ]
m2 <- as.matrix(exprs(dCT.norm2))[which(sd1 <= 0.5), ]
hm.qpcr <- heatmap.3(m1, 
                      trace = "none", 
                      cexCol = 0.6, 
                      main = "qPCR array [ddCT]", 
                      hclustfun=function(x) hclust(x,method="ward.D"),
                      col = redgreen(7))

hcrow <- as.hclust(hm.qpcr$rowDendrogram)
sh1 <- list()
for(i in 1:9){
   sh1[[i]] <- silhouette(cutree(hcrow, k = i+1), daisy(m1))
}

# pdf("shilhouetteHeatAutosome.pdf",pointsize=10, height=8,width=8)
par(mfrow=c(3,3))
for(i in 1:9)
  plot(sh1[[i]])
dev.off()

# looks like three clusters are most stable
c1 <- as.data.frame(cbind(autok2=cutree(hcrow, k = 2),autok3=cutree(hcrow, k = 3),autok4=cutree(hcrow, k = 4)))

rowSideCol <- as.matrix(data.frame(class1 = c("red", "green", "blue", "yellow")[match(c1[,1], c(1, 2))],
                         class2 = c("red", "green", "blue", "yellow")[match(c1[,2], c(1, 2, 3))],
                         class3 = c("red", "green", "blue", "yellow")[match(c1[,3], c(1, 2, 3, 4))]))
                 
hm.qpcr <- heatmap.3(m1, 
                     trace = "none", 
                     cexCol = 0.6, 
                     main = "qPCR array [ddCT]", 
                     hclustfun=function(x) hclust(x,method="ward.D"),
                     col = redgreen(15),
                     RowSideColors = as.matrix(data.frame(" " = rowSideCol[,2], "Clusters" = rowSideCol[,2])),
                     labCol = c(paste("Epithelial", seq(1,4,1), sep = " "), paste("Mesenchymal", seq(1,4,1), sep = " ")))

# extract the genes for preparing BED files corresponding to their groups
class3Genes <- list(cluster1.red = rownames(c1[which(c1$autok4 == 1),]),
                    cluster2.green = rownames(c1[which(c1$autok4 == 2),]),
                    cluster3.blue = rownames(c1[which(c1$autok4 == 3),]),
                    cluster4.yellow = rownames(c1[which(c1$autok4 == 4),]))

lapply(seq_along(class3Genes), function(i){
  gr1 <- gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% class3Genes[[i]])]
  export(gr1, paste("~/Data/Tremethick/EMT/GenomeWide/deepTools/regionFiles/", paste(names(class3Genes[i])), ".bed", sep = ""), format =  "BED")
})

x <- class3Genes
lapply(seq_along(x), function(i) paste(names(x)[[i]], x[[i]]))

# looks like three clusters are most stable
c1 <- as.data.frame(cbind(autok2=cutree(hcrow, k = 2),autok3=cutree(hcrow, k = 3),autok4=cutree(hcrow, k = 4)))
rowSideCol <- as.matrix(data.frame(class1 = c("red", "green", "blue", "yellow")[match(c1[,1], c(1, 2))],
                                   class2 = c("red", "green", "blue", "yellow")[match(c1[,2], c(1, 2, 3))],
                                   class3 = c("red", "green", "blue", "yellow")[match(c1[,3], c(1, 2, 3, 4))]))

class2Genes <- list(sd05.cluster1.red = rownames(c1[which(c1$autok3 == 1),]),
                    sd05.cluster2.green = rownames(c1[which(c1$autok3 == 2),]),
                    sd05.cluster3.blue = rownames(c1[which(c1$autok3 == 3),]))
lowSDgenes <- list(lowSD = rownames(m2))

lapply(seq_along(class2Genes), function(i){
  gr1 <- gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% class2Genes[[i]])]
  export(gr1, paste("~/Data/Tremethick/EMT/GenomeWide/deepTools/regionFiles/", paste(names(class2Genes[i])), ".bed", sep = ""), format =  "BED")
})

lapply(seq_along(lowSDgenes), function(i){
  gr1 <- gr.qPCRGenesPositions[which(gr.qPCRGenesPositions$hgnc_symbol %in% lowSDgenes[[i]])]
  export(gr1, paste("~/Data/Tremethick/EMT/GenomeWide/deepTools/regionFiles/", paste(names(lowSDgenes[i])), ".bed", sep = ""), format =  "BED")
})


