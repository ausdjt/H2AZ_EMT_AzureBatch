setwd("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/")
load("aT.nfkbSites.gr")
load("aT.ap1sites.rda")
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/biomTrack.Cfam3.Ensembl.rda")
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/aT.TSS.rda")
load("/Users/u1001407/Data/Tremethick/EMT/genomeWide.50kbTSS.DataTracks.rda")

# Plotting of all mesenchymal markers on single page
displayPars(biomTrack) <- list(stacking = "dense")

plot.rows <- 5
plot.cols <- 2
gr.plot <- promoters(gr.mesenchymalMarkers.genes, up = 1700, down = 20000)
# remove second B9D2/TGFb gene (identical)
gr.plot <- gr.plot[-3] # remove

grid.newpage()
vp <- viewport(width = 1, height = 1)
pushViewport(vp)
vp1 <- viewport(x = 0.25, y = 0.9, width = 0.45, height = 0.2)
pushViewport(vp1)
gr.plot <- promoters(gr.mesenchymalMarkers.genes, up = 1700, down = 20000)
i <- 3
chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.plot)[i]

displayPars(dT.cov.input.emt_markers.wt) <- list(transformation = NULL)
displayPars(dT.cov.input.emt_markers.tgfb) <- list(transformation = NULL)
displayPars(dT.cov.h2az.emt_markers.wt) <- list(transformation = NULL)
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(transformation = NULL)

displayPars(dT.cov.input.emt_markers.wt) <- dpList
displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList

max.y.tss <- max(max(score(dT.cov.input.emt_markers.wt)), max(score(dT.cov.input.emt_markers.tgfb)), max(score(dT.cov.h2az.emt_markers.wt)), max(score(dT.cov.h2az.emt_markers.tgfb)))
displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))

dpList <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white", "cex.title" = 0.4, rotation.title = 270, cex.axis = 0.3)

displayPars(dT.cov.input.emt_markers.wt) <- dpList
displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList


plotTracks(list(biomTrack,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = as(seqnames(gr.plot), "character")[i],
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

vp2 <- viewport(x = 0.25, y = 0.25, width = 0.45, height = 0.5)
pushViewport(vp2)
i <- 3
gr.plot <- promoters(gr.mesenchymalMarkers, up = 600, down = 1500)

displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,400))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,400))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,400))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,400))

plotTracks(list(dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb
                
),
chromosome = as(seqnames(gr.plot), "character")[i],
from = as.integer(start(gr.plot[i]), "integer"),
to = as.integer(end(gr.plot[i]), "integer"),
extend.right = 0,
extend.left = 0,
main = "",
strand = "*",
cex.main = 0.5,
sizes = c(0.25, 0.25, 0.25, 0.25),
scale = 0.5,
add = TRUE)
upViewport()

vp3 <- viewport(x = 0.75, y = 0.25, width = 0.45, height = 0.5)
pushViewport(vp3)
i <- 5
gr.plot <- gr.mesenchymalMarkers[i]
start(gr.plot) <- gr.plot$transcript_start
gr.plot <- promoters(gr.plot, up = 200, down = 2500)
displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,200))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,200))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,200))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,200))

plotTracks(list(dT.cov.input.emt_markers.wt, 
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
sizes = c(0.25, 0.25, 0.25, 0.25),
scale = 0.5,
add = TRUE)

upViewport()




