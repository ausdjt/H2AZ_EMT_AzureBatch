# TFBindingSiteSearch.R

# load libraries
library("JASPAR2014")
library("Biostrings")
library("BSgenome.Cfamiliaris.UCSC.canFam3")
library("TFBSTools")

genome <- BSgenome.Cfamiliaris.UCSC.canFam3
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))

#----------searching for AP-1 sites-----------------------------
# Using pattern matching to find putative AP-1 binding sites
# ap1Motif <- DNAString("ATGAGTCAT")
# m1 <- matchPattern(ap1Motif, genome$"1", max.mismatch = 1)
# gr.ap1Sites <- GRanges(seqnames = "1", as(m1, "IRanges"))
# subsetByOverlaps(gr.which.tss.nostrand, gr.ap1Sites)

# AP-1 is a heterodimeric protein, consisting of c-Jun/c-Fos/ATF/JDP families
# JASPAR contains matrices (all human)
# c-Jun - MA0488.1
# c-Fos - MA0476.1
# Atf1_1 - PB0004.1
# Atf1_2 - PB0108.1
opts <- list()
opts[["ID"]] = c("MA0488.1", "MA0476.1", "PB0004.1", "PB0108.1")
opts[["name"]] = "NFKB1"
opts[["all_versions"]] = TRUE
PFMatrixList <- getMatrixSet(JASPAR2014, opts = opts)
pwmList <- PWMatrixList(MA0488.1 = toPWM(PFMatrixList[["MA0488.1"]]), 
                        MA0476.1 = toPWM(PFMatrixList[["MA0476.1"]]), 
                        PB0004.1 = toPWM(PFMatrixList[["PB0004.1"]]),
                        PB0108.1 = toPWM(PFMatrixList[["PB0108.1"]]), 
                        use.names = TRUE)

grl.ap1Sites <- GRangesList(sapply(seqlevels(gr.which), function(x){
  subject <- genome[[as(x, "character")]]
  sitesetList <- searchSeq(pwmList, subject, seqname = as(x, "character"), min.score = "95%", strand = "*") 
  gr.ap1Sites <- c(as(sitesetList[[1]], "GRanges"),
                   as(sitesetList[[2]], "GRanges"), 
                   as(sitesetList[[3]], "GRanges"),
                   as(sitesetList[[4]], "GRanges"))
  return(gr.ap1Sites)
})
)

gr.ap1Sites <- unlist(grl.ap1Sites)
gr.ap1Sites <- sort(gr.ap1Sites)
gr.ap1Sites <- reduce(gr.ap1Sites)
# extend CDS by 10kb
aT.ap1Sites <- AnnotationTrack(subsetByOverlaps(gr.ap1Sites, gr.which), name = "AP1", col = "blue")

displayPars(aT.ap1Sites) <- list("fontcolor.title" = "black", 
                                 "background.title" = "white", 
                                 "col.axis" = "black", 
                                 "col.frame" = "white",
                                 cex.title = 0.5,
                                 rotation.title = 0)

#----------searching for NFKB sites-----------------------------------
# NFKB1 (4 specie consensus) - MA0105.1
# NFKB1 (human) - MA0105.2
# NFKB1 (human) - MA0105.3
opts <- list()
opts[["ID"]] = c("MA0105.1", "MA0105.2", "MA0105.3")
opts[["name"]] = "NFKB1"
opts[["all_versions"]] = TRUE
PFMatrixList <- getMatrixSet(JASPAR2014, opts = opts)
pwmList <- PWMatrixList(MA0105.1 = toPWM(PFMatrixList[["MA0105.1"]]),
                        MA0105.2 = toPWM(PFMatrixList[["MA0105.2"]]),
                        MA0105.3 = toPWM(PFMatrixList[["MA0105.3"]]), 
                        use.names = TRUE)

grl.nfkbSites <- GRangesList(sapply(seqlevels(gr.which), function(x){
  subject <- genome[[as(x, "character")]]
  sitesetList <- searchSeq(pwmList, subject, seqname = as(x, "character"), min.score = "95%", strand = "*") 
  gr.nfkbSites <- c(as(sitesetList[[1]], "GRanges"),
                   as(sitesetList[[2]], "GRanges"), 
                   as(sitesetList[[3]], "GRanges"))
  return(gr.nfkbSites)
}))

gr.nfkbSites <- unlist(grl.nfkbSites)
gr.nfkbSites <- sort(gr.nfkbSites)
gr.nfkbSites <- reduce(gr.nfkbSites)
aT.nfkbSites <- AnnotationTrack(subsetByOverlaps(gr.nfkbSites, gr.which), name = "NFKB", col = "salmon")

displayPars(aT.nfkbSites) <- list("fontcolor.title" = "black", 
                                  "background.title" = "white", 
                                  "col.axis" = "black", 
                                  "col.frame" = "white",
                                  cex.title = 0.5,
                                  rotation.title = 0)

save(aT.ap1Sites, file = "~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/aT.ap1Sites.rda")
save(aT.nfkbSites, file = "~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/aT.nfkbSites.rda")



