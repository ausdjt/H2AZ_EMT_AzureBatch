# This is the server logic of a Shiny web application. You can run the 
#
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Data loading and preparation
load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/resultsCompressed.rda", .GlobalEnv)
load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/cfam.qPCRGenesTab.rda", .GlobalEnv)
load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/ensGenes_CanFam3.1_ensembl84_ERCC.rda", .GlobalEnv)

globalData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% ensGenes$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(globalData) <- names(resultsCompressed[[1]]$sleuth_results.gene)

emtMarkerData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% cfam.qPCRGenesTab$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(emtMarkerData) <- names(resultsCompressed[[1]]$sleuth_results.gene)
  
# Actual Server logic
shinyServer(function(input, output) {
  output$volcano_tgfb <- renderPlot({
    y <- "conditionMDCKTGFb"
    dat <- emtMarkerData[[y]]$dataTable
    xAxisMax <- max(abs(dat$b)) + 1
    plot(dat$b,
         -log10(dat$qval), 
         axes = F, 
         xlab = "", 
         ylab = "", 
         frame = F,
         cex = 0.3,
         xlim = c(-round(xAxisMax, 0), round(xAxisMax,0)),
         pch = 16, main = paste("Volcano plot\nCondition: ", y, sep = ""))
    points(dat[which(-log10(dat$qval) >= 1), "b"], 
           -log10(dat[which(-log10(dat$qval) >= 1), "qval"]),
           col = "red", 
           pch = 16, 
           cex = 1.1)
    axis(2,
         pos = 0, 
         lwd = 3)
    #                             at = c(seq(0,yAxisMax,10)))
    axis(1,
         pos = 0, 
         lwd = 3,
         at = c(seq(-round(xAxisMax, 0), round(xAxisMax,0), 2)))
    mtext("-log10(q-value)", 
          side = 2)
    mtext("beta value",
          side = 1, 
          line = 2)
    abline(h = 1, col = "red", lty = 2, lwd = 2)
    # only adding labels to genes with adjuste p-value <= 0.01
    text(dat[which(-log10(dat$qval) >= 2), "b"], -log10(dat[which(-log10(dat$qval) >= 2), "qval"]), 
         labels = dat[which(-log10(dat$qval) >= 2), "external_gene_name" ],
         cex = 0.7,
         pos = 4, offset = 0.3)
    text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
  })
  
  output$volcano_shZ <- renderPlot({
    y <- "conditionMDCKshZ"
    dat <- emtMarkerData[[y]]$dataTable
    xAxisMax <- max(abs(dat$b)) + 1
    plot(dat$b,
         -log10(dat$qval), 
         axes = F, 
         xlab = "", 
         ylab = "", 
         frame = F,
         cex = 0.3,
         xlim = c(-round(xAxisMax, 0), round(xAxisMax,0)),
         pch = 16, main = paste("Volcano plot\nCondition: ", y, sep = ""))
    points(dat[which(-log10(dat$qval) >= 1), "b"], 
           -log10(dat[which(-log10(dat$qval) >= 1), "qval"]),
           col = "red", 
           pch = 16, 
           cex = 1.1)
    axis(2,
         pos = 0, 
         lwd = 3)
    #                             at = c(seq(0,yAxisMax,10)))
    axis(1,
         pos = 0, 
         lwd = 3,
         at = c(seq(-round(xAxisMax, 0), round(xAxisMax,0), 2)))
    mtext("-log10(q-value)", 
          side = 2)
    mtext("beta value",
          side = 1, 
          line = 2)
    abline(h = 1, col = "red", lty = 2, lwd = 2)
    # only adding labels to genes with adjuste p-value <= 0.01
    text(dat[which(-log10(dat$qval) >= 2), "b"], -log10(dat[which(-log10(dat$qval) >= 2), "qval"]), 
         labels = dat[which(-log10(dat$qval) >= 2), "external_gene_name" ],
         cex = 0.7,
         pos = 4, offset = 0.3)
    text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
  })
  
  output$global_de_tgfb = DT::renderDataTable(globalData[["conditionMDCKTGFb"]]$dataTable)
  output$global_de_shZ = DT::renderDataTable(globalData[["conditionMDCKshZ"]]$dataTable)
  output$emtMarker_de_tgfb = DT::renderDataTable(emtMarkerData[["conditionMDCKTGFb"]]$dataTable)
  output$emtMarker_de_shZ = DT::renderDataTable(emtMarkerData[["conditionMDCKshZ"]]$dataTable)
})
