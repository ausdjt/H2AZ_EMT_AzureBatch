setwd("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/")
load("aT.nfkbSites.gr")
load("aT.ap1sites.rda")
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/biomTrack.Cfam3.Ensembl.rda")
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/aT.TSS.rda")
load("/Users/u1001407/Data/Tremethick/EMT/genomeWide.50kbTSS.DataTracks.rda")

# TGFb-1 plotting
axisTrack <- GenomeAxisTrack()
TxDb.Cfam3.Ensembl <- makeTxDbFromGFF("/Volumes/gduserv/Data/Annotations/CanFam3/TGFB1.gtf")
biomTrack <- GeneRegionTrack(TxDb.Cfam3.Ensembl, showId = T, geneSymbol = T, showExonId = F, name = "", stacking = "hide")
displayPars(biomTrack) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")


dpList <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white", "cex.title" = 0.6, rotation.title = 270, cex.axis = 0.6)

displayPars(dT.cov.input.emt_markers.wt) <- dpList
displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList
displayPars(aT.danpos2) <- dpList

# actual plotting
pdf("H2AZ_occupancy_TGFb-1 - Version 2015-12-09.pdf", width = 8.5, height = 12)
grid.newpage()
vp <- viewport(width = 1, height = 1)
pushViewport(vp)

vp1 <- viewport(y = 0.75, width = 0.9, height = 0.5)
pushViewport(vp1)

gr.plot <- promoters(gr.mesenchymalMarkers, up = 1700, down = 20000)
chr <- 1
i <- 3
chromosome(dT.cov.input.emt_markers.wt) <- chr
chromosome(dT.cov.input.emt_markers.tgfb) <- chr
chromosome(dT.cov.h2az.emt_markers.wt) <- chr
chromosome(dT.cov.h2az.emt_markers.tgfb) <- chr

displayPars(dT.cov.input.emt_markers.wt) <- list(transformation = NULL)
displayPars(dT.cov.input.emt_markers.tgfb) <- list(transformation = NULL)
displayPars(dT.cov.h2az.emt_markers.wt) <- list(transformation = NULL)
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(transformation = NULL)
displayPars(biomTrack) <- list(stacking = "squish")


max.y.tss <- max(max(score(dT.cov.input.emt_markers.wt)), max(score(dT.cov.input.emt_markers.tgfb)), max(score(dT.cov.h2az.emt_markers.wt)), max(score(dT.cov.h2az.emt_markers.tgfb)))
displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))

plotTracks(list(axisTrack,
                biomTrack,
                aT.TSS,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = as(seqnames(gr.plot), "character")[i],
from = as.integer(start(gr.plot[i]), "integer"),
to = as.integer(end(gr.plot[i]), "integer"),
extend.right = 0,
extend.left = 2500,
main = "",
#main = paste(gr.plot$hgnc_symbol[i], " (", width(gr.plot[i]), "bp)", sep = ""),
strand = "*",
cex.main = 0.5,
sizes = c(0.03, 0.14, 0.04, 0.1975, 0.1975, 0.1975, 0.1975),
scale = 0.5,
add = TRUE)

upViewport()
vp2 <- viewport(x = 0.25, y = 0.25, width = 0.45, height = 0.4)
pushViewport(vp2)
i <- 3
gr.plot <- promoters(gr.mesenchymalMarkers, up = 600, down = 1500)

displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
displayPars(biomTrack) <- list(stacking = "dense")

plotTracks(list(biomTrack,
#               aT.danpos2,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = chr,
from = as.integer(start(gr.plot[i]), "integer"),
to = as.integer(end(gr.plot[i]), "integer"),
extend.right = 0,
extend.left = 0,
main = "",
strand = "*",
cex.main = 0.5,
sizes = c(0.04, 0.24, 0.24, 0.24, 0.24),
scale = 0.5,
add = TRUE)
upViewport()

vp3 <- viewport(x = 0.75, y = 0.25, width = 0.45, height = 0.4)
pushViewport(vp3)
i <- 5
gr.plot <- gr.mesenchymalMarkers[i]
start(gr.plot) <- gr.plot$transcript_start
gr.plot <- promoters(gr.plot, up = 200, down = 2500)
displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,5))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,5))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,5))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,5))
displayPars(biomTrack) <- list(stacking = "dense")
plotTracks(list(biomTrack,
#                axisTrack,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = as(seqnames(gr.plot), "character"),
from = as.integer(start(gr.plot)),
to = as.integer(end(gr.plot)),
extend.right = 0,
extend.left = 0,
main = "",
strand = "*",
cex.main = 0.5,
sizes = c(0.04, 0.24, 0.24, 0.24, 0.24),
scale = 0.5,
add = TRUE)

upViewport()
dev.off()
 


