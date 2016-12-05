# This is the server logic of a Shiny web application

library(shiny)
library(DT)


# load external data files ------------------------------------------------
load("data/resultsCompressed.rda", .GlobalEnv)
load("data/cfam.qPCRGenesTab.rda", .GlobalEnv)
load("data/ensGenes_CanFam3.1_ensembl84_ERCC.rda", .GlobalEnv)
load("data/cfamEnsGenesSigEMTCells.rda", .GlobalEnv)


# prepare data sets -----------------------------------------------------
# gene expression data for all samples 
globalGeneExpData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% ensGenes$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(globalGeneExpData) <- names(resultsCompressed[[1]]$sleuth_results.gene)


# global DE data used for tables of all genes
globalDEData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% ensGenes$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(globalDEData) <- names(resultsCompressed[[1]]$sleuth_results.gene)

# this is the data corresponding to the genes present on the EMT qPCR array
emtqPCRData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% cfam.qPCRGenesTab$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(emtqPCRData) <- names(resultsCompressed[[1]]$sleuth_results.gene)

# this is the data corresponding to the EMT cell line signature from Tan et al. 2014
emtSignatureData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% cfamEnsGenesSigEMTCells$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- merge(dat, 
               cfamEnsGenesSigEMTCells[,c("ensembl_gene_id", "epi_mes")],
               by.x = "target_id",
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  dat$log10_qval <- -log10(dat$qval)
  return(list(dataTable = dat))
})
names(emtSignatureData) <- names(resultsCompressed[[1]]$sleuth_results.gene)

#********************End of data preparation*******************************

# Actual Server logic -----------------------------------------------------
shinyServer(function(input, output, session) {
  
  # add an HTML file
  getPage<-function() {
    return(includeHTML("mdck_rna-seq.html"))
  }
  output$inc_html <- renderUI({getPage()})
 
  # volcano plots
  output$volcano_tgfb <- renderPlot({
    y <- "conditionMDCKTGFb"
    dat <- emtqPCRData[[y]]$dataTable
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
    dat <- emtqPCRData[[y]]$dataTable
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
  
  selCols <- c("target_id", "qval", "b", "external_gene_name", "description")
  output$global_de_tgfb = DT::renderDataTable(globalDEData[["conditionMDCKTGFb"]]$dataTable[,selCols])
  output$global_de_shZ = DT::renderDataTable(globalDEData[["conditionMDCKshZ"]]$dataTable[,selCols])
  output$emtMarker_de_tgfb = DT::renderDataTable(emtqPCRData[["conditionMDCKTGFb"]]$dataTable[,selCols])
  output$emtMarker_de_shZ = DT::renderDataTable(emtqPCRData[["conditionMDCKshZ"]]$dataTable[,selCols])
  
  # interactive volcano plot of the EMT signature data
  highlight <- reactive({
    switch(input$highlight,
           epi = "epi",
           mes = "mes",
           none = "none")
  })
  
  output$text <- renderText(paste("you have selected", input$dataset))
  
  emtSigPlotData <- reactive({
    switch(input$dataset, 
           conditionMDCKTGFb = emtSignatureData[["conditionMDCKTGFb"]]$dataTable,
           conditionMDCKshZ = emtSignatureData[["conditionMDCKshZ"]]$dataTable
          )
  }) #emtSigPlotData
  
  # setup UI for interactive plot
  output$plotui <- renderUI({
    plotOutput("plot", height=600,
               click = "plot_click",
               dblclick = dblclickOpts(
                 id = "plot_dblclick"
               ),
               hover = hoverOpts(
                 id = "plot_hover"
               ),
               brush = brushOpts(
                 id = "plot_brush"
               )
    )
  }) # output$plotui
  
  # actual plot for interactive volcano plot
  output$plot <- renderPlot({
    dat <- emtSigPlotData()
    
    xAxisMax <- max(abs(dat$b)) + 1
    plot(dat$b, dat$log10_qval, axes = F, xlab = "", ylab = "", frame = F, 
         cex = 0.3,xlim = c(-round(xAxisMax, 0), round(xAxisMax,0)), pch = 16, 
         main = paste(input$dataset))
    if (input$highlight == "epi") {
      w1 <- dat$epi_mes == "epi"
      colPoint <- c("blue", "black")[match(dat[w1,]$epi_mes, c("epi", "mes"))]
      points(dat[w1, "b"], dat[w1, "log10_qval"], col = colPoint, pch = 16, cex = 1.1)
    } else if (input$highlight == "mes") {
      w1 <- dat$epi_mes == "mes"
      colPoint <- c("black", "green")[match(dat[w1,]$epi_mes, c("epi", "mes"))]
      points(dat[w1, "b"], dat[w1, "log10_qval"], col = colPoint, pch = 16, cex = 1.1)
    } else if (input$highlight == "both") {
      w1 <- grep("epi|mes", dat$epi_mes)
      colPoint <- c("blue", "green")[match(dat[w1,]$epi_mes, c("epi", "mes"))]
      points(dat[w1, "b"], dat[w1, "log10_qval"], col = colPoint, pch = 16, cex = 1.1)
    } else {
      w1 <- grep("epi|mes", dat$epi_mes);
      colPoint <- "black"
      points(dat[w1, "b"], dat[w1, "log10_qval"], pch = 16, cex = 0.4, col = colPoint)
    }
    
    # always highlighting H2AFZ & TGFB1
    w2 <- grep("H2AFZ|TGFB1", dat$external_gene_name)
    points(dat[w2, "b"], dat[w2, "log10_qval"], col = "red", pch = 16, cex = 1.3)
    text(dat[w2, "b"], dat[w2, "log10_qval"], 
         labels = dat[w2, "external_gene_name" ],
         cex = 0.8,
         pos = 4, offset = 0.3)
    
    # Axis 
    axis(2, pos = 0, lwd = 3)
    axis(1, pos = 0, lwd = 3, at = c(seq(-round(xAxisMax, 0), round(xAxisMax,0), 2)))
    mtext("-log10(q-value)", side = 2)
    mtext("beta value", side = 1, line = 2)
    abline(h = 1, col = "red", lty = 2, lwd = 2)
   }) # output$plot
  
  # setup which columns to display
  selColsGlobal <- c(selCols, "epi_mes")
  
  # generate datatable for the selected points from the interactive plot
  selColsGlobal <- c(selCols, "epi_mes")
  
  output$plot_brushed_points <- DT::renderDataTable({
    dat <- emtSigPlotData()
    res <- brushedPoints(dat, input$plot_brush, "b", "log10_qval")
    datatable(res[,selColsGlobal])
  }) #output$plot_brushed_points
  
 })
