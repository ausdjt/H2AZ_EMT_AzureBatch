setwd("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/")
load("aT.nfkbSites.rda")
load("aT.ap1sites.rda")
load("biomTrack.Cfam3.Ensembl.rda")
load("aT.TSS.rda")
load("genomeWide.50kbTSS.DataTracks.rda")
load("gr.danpos2.results.rda")

# Plotting of all epithelial markers on single page
displayPars(biomTrack) <- list(stacking = "dense")
dpList <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "black", frame =  FALSE, "cex.title" = 0.4, rotation.title = 270, cex.axis = 0.3)


gr.plot <- promoters(gr.epithelialMarkers.genes, up = 1700, down = 20000)
gr.plot.detail <- promoters(gr.epithelialMarkers.genes, up = 1500, down = 1500)
gr.danpos2.results.subSet <- subsetByOverlaps(gr.danpos2.results, gr.plot)

plot.rows <- 5
plot.cols <- 2
yAxis <- seq(0.1, 0.9, 0.2)
xAxis <- c(0.375, 0.8)

pdf("EpithelialMarkers_SinglePage.pdf", width = 10, height = 15)

grid.newpage()
lapply(seq_along(1:plot.rows), function(i){
  chromAct <- as.character(seqnames(gr.plot))[i]
  chromosome(dT.cov.input.emt_markers.wt) <- chromAct
  chromosome(dT.cov.input.emt_markers.tgfb) <- chromAct
  chromosome(dT.cov.h2az.emt_markers.wt) <- chromAct
  chromosome(dT.cov.h2az.emt_markers.tgfb) <- chromAct
  
  displayPars(dT.cov.input.emt_markers.wt) <- list(transformation = NULL)
  displayPars(dT.cov.input.emt_markers.tgfb) <- list(transformation = NULL)
  displayPars(dT.cov.h2az.emt_markers.wt) <- list(transformation = NULL)
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(transformation = NULL)
  
  displayPars(dT.cov.input.emt_markers.wt) <- dpList
  displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
  displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList
  
  max.y <- max(max(score(dT.cov.input.emt_markers.wt)), max(score(dT.cov.input.emt_markers.tgfb)), max(score(dT.cov.h2az.emt_markers.wt)), max(score(dT.cov.h2az.emt_markers.tgfb)))
  displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y))
  
  displayPars(dT.cov.input.emt_markers.wt) <- dpList
  displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
  displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList
  
  lapply(seq_along(1:plot.cols), function(j){
    if (j == 1){
      vp1 <- viewport(x = xAxis[j], y = yAxis[i], width = 0.7, height = 0.2)
      pushViewport(vp1)
      dpList <- list(showTitle = TRUE, lwd.title = 0.0001, box.width = 0.5)
      displayPars(dT.cov.input.emt_markers.wt) <- dpList
      displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
      displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
      displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList
      
      plotTracks(list(biomTrack,
                      dT.cov.input.emt_markers.wt, 
                      dT.cov.h2az.emt_markers.wt,
                      dT.cov.input.emt_markers.tgfb, 
                      dT.cov.h2az.emt_markers.tgfb
                      
      ),
      chromosome = chromAct,
      from = as.integer(start(gr.plot[i]), "integer"),
      to = as.integer(end(gr.plot[i]), "integer"),
      extend.right = 0,
      extend.left = 1400,
      main = paste(gr.plot$hgnc_symbol[i], " (", width(gr.plot[i]), "bp)", sep = ""),
      strand = "*",
      cex.main = 0.3,
      sizes = c(0.1, rep(0.225, 4)),
      add = TRUE)
      upViewport()
    } else {
      vp1 <- viewport(x = xAxis[j], y = yAxis[i], width = 0.25, height = 0.2)
      pushViewport(vp1)
      # plot details around the maximal position of h2az in TGFB
      dpList <- list(showTitle = FALSE, lwd.title = 0.1, showAxis = F, col.frame = "black")
      
      displayPars(dT.cov.input.emt_markers.wt) <- dpList
      displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
      displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
      displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList
      pos <- start(gr.danpos2.results.subSet[which(seqnames(gr.danpos2.results.subSet) == chromAct)]
                   [which.max(abs(gr.danpos2.results.subSet[which(seqnames(gr.danpos2.results.subSet) == chromAct)]$control_smt_val))])
      max.y <- max(bA.cov.h2az.emt_markers.wt[which(seqnames(bA.cov.h2az.emt_markers.wt) == chromAct)]$mean)
      displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y))
      displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y))
      displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y))
      displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y))
      
      plotTracks(list(biomTrack,
                      dT.cov.input.emt_markers.wt, 
                      dT.cov.h2az.emt_markers.wt,
                      dT.cov.input.emt_markers.tgfb, 
                      dT.cov.h2az.emt_markers.tgfb
                      
      ),
      chromosome = as(seqnames(gr.plot), "character")[i],
      from = pos - 1000,
      to = pos + 1000,
      extend.right = 100,
      extend.left = 100,
      main = paste(gr.plot$hgnc_symbol[i], " [2000 bp]", sep = ""),
      strand = "*",
      cex.main = 0.3,
      sizes = c(0.1, rep(0.225, 4)),
      add = TRUE)
      upViewport()
    }
  })
})

dev.off()
