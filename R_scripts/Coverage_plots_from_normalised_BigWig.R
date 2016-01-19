# coverage from the merged BAM file
counts.h2az_tgfb <- countBam("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb.bam", param = SBParam.all)
reads.h2az_tgfb <- readGAlignments("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb.bam", param = SBParam, use.names = F)
reads.h2az_tgfb <- as(reads.h2az_tgfb, "GRangesList")
reads.h2az_tgfb <- unlist(reads.h2az_tgfb)
width(reads.h2az_tgfb) <- 160
scale.h2az_tgfb <-1000000/counts.h2az_tgfb$records
cov.h2az_tgfb <- coverage(reads.h2az_tgfb) * scale.h2az_tgfb
cov.h2az_tgfb <- cov.h2az_tgfb[which(names(cov.h2az_tgfb) %in% seqlevels(gr.which.tiles))]
cov.h2az_tgfb <- cov.h2az_tgfb[names(cov.h2az_tgfb)[order(names(cov.h2az_tgfb))]]
bA.cov.h2az_tgfb <- binnedAverage(gr.which.tiles, cov.h2az_tgfb, "mean")
dT.h2az_tgfb <- DataTrack(bA.cov.h2az_tgfb, type = "h", col = "black", name = "H2AZ TGFb merged [rpm]")

counts.h2az_wt <- countBam("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_WT.bam", param = SBParam.all)
reads.h2az_wt <- readGAlignments("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_WT.bam", param = SBParam, use.names = F)
reads.h2az_wt <- as(reads.h2az_wt, "GRangesList")
reads.h2az_wt <- unlist(reads.h2az_wt)
width(reads.h2az_wt) <- 160
scale.h2az_wt <-1000000/counts.h2az_wt$records
cov.h2az_wt <- coverage(reads.h2az_wt) * scale.h2az_wt
cov.h2az_wt <- cov.h2az_wt[which(names(cov.h2az_wt) %in% seqlevels(gr.which.tiles))]
cov.h2az_wt <- cov.h2az_wt[names(cov.h2az_wt)[order(names(cov.h2az_wt))]]
bA.cov.h2az_wt <- binnedAverage(gr.which.tiles, cov.h2az_wt, "mean")
dT.h2az_wt <- DataTrack(bA.cov.h2az_wt, type = "h", col = "black", name = "H2AZ wt merged [rpm]")

# H2AZ TGFb treated separate samples 
cov.h2az.emt_markers.tgfb_rep1 <- coverage(reads.h2az[["H2AZ_TGFb_rep1_S3"]]) * scale.h2az[["H2AZ_TGFb_rep1_S3"]]
cov.h2az.emt_markers.tgfb_rep2 <- coverage(reads.h2az[["H2AZ_TGFb_rep2_S4"]]) * scale.h2az[["H2AZ_TGFb_rep2_S4"]]

cov.h2az.emt_markers.tgfb_rep1 <- cov.h2az.emt_markers.tgfb_rep1[which(names(cov.h2az.emt_markers.tgfb_rep1) %in% seqlevels(gr.which.tiles))]
cov.h2az.emt_markers.tgfb_rep1 <- cov.h2az.emt_markers.tgfb_rep1[names(cov.h2az.emt_markers.tgfb_rep1)[order(names(cov.h2az.emt_markers.tgfb_rep1))]]

cov.h2az.emt_markers.tgfb_rep2 <- cov.h2az.emt_markers.tgfb_rep2[which(names(cov.h2az.emt_markers.tgfb_rep2) %in% seqlevels(gr.which.tiles))]
cov.h2az.emt_markers.tgfb_rep2 <- cov.h2az.emt_markers.tgfb_rep2[names(cov.h2az.emt_markers.tgfb_rep2)[order(names(cov.h2az.emt_markers.tgfb_rep2))]]


# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.h2az.emt_markers.tgfb_rep1 <- binnedAverage(gr.which.tiles, cov.h2az.emt_markers.tgfb_rep1, "mean")
bA.cov.h2az.emt_markers.tgfb_rep2 <- binnedAverage(gr.which.tiles, cov.h2az.emt_markers.tgfb_rep2, "mean")

dT.cov.h2az.emt_markers.tgfb_rep1 <- DataTrack(bA.cov.h2az.emt_markers.tgfb_rep1, type = "h", col = "red", name = "H2AZ TGFb Rep1 [rpm]")
dT.cov.h2az.emt_markers.tgfb_rep2 <- DataTrack(bA.cov.h2az.emt_markers.tgfb_rep2, type = "h", col = "red", name = "H2AZ TGFb Rep2 [rpm]")

# Plot TGFb
axisTrack <- GenomeAxisTrack()
i <- 3

gr.plot <- gr.mesenchymalMarkers.1500TSS1500
gr.plot <- gr.mesenchymalMarkers.genes.25kbTSS25kb

biomTrack <- GeneRegionTrack(TxDb.Cfam3.Ensembl,
                             chromosome = as(seqnames(gr.plot), "character")[i],
                             start = as.integer(start(gr.plot[i])),
                             end = as.integer(end(gr.plot[i])),
                             name = paste(mcols(gr.plot[i])$hgnc_symbol, "transcript",  mcols(gr.plot[i])$ensembl_transcript_id, sep = " "),
                             genome = "canFam3")
displayPars(biomTrack) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.tgfb_rep1) <- seqnames(gr.plot)[i]
chromosome(dT.cov.h2az.emt_markers.tgfb_rep2) <- seqnames(gr.plot)[i]
chromosome(dT.h2az_tgfb) <- seqnames(gr.plot)[i]
chromosome(dT.h2az_wt) <- seqnames(gr.plot)[i]

max.y.tss <- max(max(values(dT.cov.input.emt_markers.wt)), 
                 max(values(dT.cov.input.emt_markers.tgfb)),
                 max(values(dT.cov.h2az.emt_markers.wt)), 
                 max(values(dT.cov.h2az.emt_markers.tgfb)),
                 max(values(dT.cov.h2az.emt_markers.tgfb_rep1)),
                 max(values(dT.cov.h2az.emt_markers.tgfb_rep2)),
                 max(values(dT.h2az_tgfb)),
                 max(values(dT.h2az_wt)))

displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0, max.y.tss))
displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0, max.y.tss))
displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0, max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0, max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb_rep1) <- list(ylim = c(0, max.y.tss))
displayPars(dT.cov.h2az.emt_markers.tgfb_rep2) <- list(ylim = c(0, max.y.tss))
displayPars(dT.h2az_tgfb) <- list(ylim = c(0, max.y.tss))
displayPars(dT.h2az_wt) <- list(ylim = c(0, max.y.tss))

plotTracks(list(axisTrack,
                biomTrack,
                aT.ap1Sites,
                aT.nfkbSites,
                aT.TSS,
                dT.cov.input.emt_markers.wt, 
                dT.cov.h2az.emt_markers.wt,
                align1,
                dT.cov.input.emt_markers.tgfb, 
                dT.cov.h2az.emt_markers.tgfb,
                dT.cov.h2az.emt_markers.tgfb_rep1,
                dT.cov.h2az.emt_markers.tgfb_rep2,
                align2
),
chromosome = as(seqnames(gr.plot), "character")[i],
from = as.integer(start(gr.plot[i]), "integer"),
to = as.integer(end(gr.plot[i]), "integer"),
extend.right = 1000,
extend.left = 1000,
main = mcols(gr.plot[i])$hgnc_symbol,
strand = "*",
cex.main = 0.5,
#sizes = c(0.01, 0.04, 0.02, 0.02, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05),
scale = 0.5)

align1 <- AlignmentsTrack("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_WT.bam", 
                          isPaired = T, 
                          chromosome = as(seqnames(gr.plot), "character")[i],
                          from = as.integer(start(gr.plot[i]), "integer"),
                          to = as.integer(end(gr.plot[i]), "integer"),
                          type = "coverage"
)

align2 <- AlignmentsTrack("~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/H2AZ_TGFb.bam", 
                          isPaired = T, 
                          chromosome = as(seqnames(gr.plot), "character")[i],
                          from = as.integer(start(gr.plot[i]), "integer"),
                          to = as.integer(end(gr.plot[i]), "integer"),
                          type = "coverage"
                          )

