# load libraries
library("gdata")
library("Gviz")
library("rtracklayer")
library("biomaRt")
library("GenomicRanges")
library("GenomicAlignments")

options(ucscChromosomeNames=FALSE)

# additional functions
# function to calculate binned averages from a coverage RleList and GRanges object
source("~/Development/GeneralPurpose/R/binnedAverage.R")
source("~/Development/GeneralPurpose/R/binnedSum.R")
# set WD
setwd("~/Data/Tremethick/MDCK_sure_select/BAMs")


#----------connect to Ensembl biomaRt for annotation data----------------------
dog <- useMart("ensembl", dataset = "cfamiliaris_gene_ensembl")
filters <- listFilters(dog)
att <- listAttributes(dog)
# # EMT markers
# # mesenchymal genes HGNC Symbols
# # FN1
# # ZEB1
# # TGFb1 -> TGFB1I1 -> ENSCAFG00000005014
# # SPARC
# # TWIST2
mesenchymalMarkers <- c("FN1", "ZEB1", "SPARC", "TWIST2")
mesenchymalMarkers.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), filters = "hgnc_symbol", values = mesenchymalMarkers, dog)
mesenchymalMarkers.TSS.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "transcription_start_site", "hgnc_symbol"), filters = "hgnc_symbol", values = mesenchymalMarkers, dog)
mesenchymalMarkers.5UTR.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "5_utr_start"), filters = "hgnc_symbol", values = mesenchymalMarkers, dog)
mesenchymalMarkers.transcripts.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "transcript_start", "hgnc_symbol", "strand"), filters = "hgnc_symbol", values = mesenchymalMarkers, dog)

tgfb.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), filters = "ensembl_gene_id", values = "ENSCAFG00000005014", dog)
tgfb.tab$hgnc_symbol <- "B9D2(TGFB1)"
tgfb.TSS.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol", "transcription_start_site"), filters = "ensembl_gene_id", values = "ENSCAFG00000005014", dog)
tgfb.transcripts.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol", "transcript_start"), filters = "ensembl_gene_id", values = "ENSCAFG00000005014", dog)

mesenchymalMarkers.tab <- rbind(mesenchymalMarkers.tab, tgfb.tab)
mesenchymalMarkers.tab$strand <- c("-", "+")[match(mesenchymalMarkers.tab$strand, c("-1", "1"))]
mesenchymalMarkers.tab$marker <- "mesenchymal"
gr.mesenchymalMarkers <- GRanges(mesenchymalMarkers.tab$chromosome_name, IRanges(mesenchymalMarkers.tab$start_position, mesenchymalMarkers.tab$end_position), strand = mesenchymalMarkers.tab$strand, hgnc_symbol = mesenchymalMarkers.tab$hgnc_symbol, marker = mesenchymalMarkers.tab$marker)
# 
# # epithelial markers
# # CDH1
# # SPP1
# # FGFBP1
# # MMP9
# # EPCAM
epithelialMarkers <- c("CDH1", "SPP1", "FGFBP1", "MMP9", "EPCAM")
epithelialMarkers.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), filters = "hgnc_symbol", values = epithelialMarkers, dog)
epithelialMarkers.TSS.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "transcription_start_site", "hgnc_symbol"), filters = "hgnc_symbol", values = epithelialMarkers, dog)
epithelialMarkers.tab$strand <- c("-", "+")[match(epithelialMarkers.tab$strand, c("-1", "1"))]
epithelialMarkers.tab$marker <- "epithelial"
gr.epithelialMarkers <- GRanges(epithelialMarkers.tab$chromosome_name, IRanges(epithelialMarkers.tab$start_position, epithelialMarkers.tab$end_position), strand = epithelialMarkers.tab$strand, hgnc_symbol = epithelialMarkers.tab$hgnc_symbol, marker = epithelialMarkers.tab$marker)
# 
# create table to define region around TSS (5' most coordinate stored)
which.tss.tab <- rbind(mesenchymalMarkers.tab, epithelialMarkers.tab)
gr.which.tss <- GRanges(which.tss.tab$chromosome_name, IRanges(which.tss.tab$start_position, which.tss.tab$end_position), strand = which.tss.tab$strand, hgnc_symbol = which.tss.tab$hgnc_symbol, marker = which.tss.tab$marker)

gr.which.cds <- gr.which.tss
start(gr.which.cds) <- start(gr.which.cds) - 5000
end(gr.which.cds) <- end(gr.which.cds) + 5000
strand(gr.which.cds) <- "*"

gr.which.tss <- promoters(gr.which.tss, upstream = 1500, downstream = 1500)
gr.which.tss.nostrand <- gr.which.tss
strand(gr.which.tss.nostrand) <- "*"

which.tss.tab$start_position <- start(gr.which.tss)
which.tss.tab$end_position <- end(gr.which.tss)
colnames(which.tss.tab)[c(3,4)] <- c("promoter_start_position", "promoter_end_position")

designTab <- data.frame(chromosome_name = c(5, 32, 3, 24, 10, 37, 2, 1, 4, 25), 
                        start_position = c(80737492, 11334119, 64429408, 33254255, 49465757, 22436951, 15282325, 112608937, 57639351, 49145461),
                        end_position = c(80777492, 11374119, 64469409, 33294255, 49505757, 22476951, 15322325, 112648937, 57679351, 49185460),
                        hgnc_symbol = c("CDH1", "SPP1", "FGFBP1", "MMP9", "EPCAM", "FN1", "ZEB1", "B9D2(TGFB1)", "SPARC", "TWIST2"),
                        marker = c(rep("epithelial", 5), rep("mesenchymal", 5)))
designTab$hgnc_symbol <- as.character(designTab$hgnc_symbol)
designTab$marker <- as.character(designTab$marker)

designTab <- merge(designTab, which.tss.tab[,c("ensembl_gene_id", "promoter_start_position", "promoter_end_position", "strand", "hgnc_symbol")], by.x = "hgnc_symbol", by.y = "hgnc_symbol")
designTab <- designTab[order(designTab$marker),]

gr.which <- GRanges(seqnames = designTab$chromosome_name, 
                    IRanges(start = designTab$start_position - 100000, end = designTab$end_position + 1000000),
                    strand = "*",
                    hgnc_symbol = as.character(designTab$hgnc_symbol),
                    marker = designTab$marker,
                    promoter_start = designTab$promoter_start_position,
                    promoter_end = designTab$promoter_end_position)

#----------load ChIP-qPCR amplicon data----------------------------------------
primers <- read.xls("../ChIP-Primers.xlsx")
gr.primers <- GRanges(seqnames = primers$Amplicon.location...Chr, IRanges(start = primers$Amplicon.location...Start, end = primers$Amplicon.location...Start), strand = "*", primerID = primers$Primer)
seqlevels(gr.primers, force = T) <- gsub("chr", "", seqlevels(gr.primers))
aT.primers <- AnnotationTrack(gr.primers, name = "ChIP-qPCR amplicons", col = "darkgrey")

#----------load annotation of capture probes-----------------------------------
gr.captureProbes <- import("../MDCK_1_Covered.bed")
seqlevels(gr.captureProbes, force = TRUE) <- gsub("chr", "", seqlevels(gr.captureProbes))
aT.captureProbes <- AnnotationTrack(gr.captureProbes, name = "Capture probes", col = "yellow")

#----------H2AZ ChIP data------------------------------------------------------
# load ChIP coverage from BAM files
# if (length(grep("chr", seqlevels(gr.which))) > 0){
#   seqlevels(gr.which, force = T) <- gsub("chr", "", seqlevels(gr.which))
# }
flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = NA)
SBParam <- ScanBamParam(flag = flag, simpleCigar = F, what = scanBamWhat(), which = gr.which)

h2az.wt.rep1 <- readGAlignmentPairs("H2AZ-WT-rep1_S1_L001_sorted_MkDup.bam", param = SBParam)
h2az.wt.rep2 <- readGAlignmentPairs("H2AZ-WT-rep2_S2_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep1 <- readGAlignmentPairs("H2AZ-TGFb-rep1_S3_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep2 <- readGAlignmentPairs("H2AZ-TGFb-rep2_S4_L001_sorted_MkDup.bam", param = SBParam)

# scaling factor for coverage across on-target regions [RPM]
h2az.wt.rep1.scale <- 1000000/length(h2az.wt.rep1)
h2az.wt.rep2.scale <- 1000000/length(h2az.wt.rep2)
h2az.tgfb.rep1.scale <- 1000000/length(h2az.tgfb.rep1)
h2az.tgfb.rep2.scale <- 1000000/length(h2az.tgfb.rep2)

# calculate coverage from both replicates
cov.h2az.wt <- (coverage(h2az.wt.rep1) * h2az.wt.rep1.scale + coverage(h2az.wt.rep2) * h2az.wt.rep2.scale) / 2
cov.h2az.tgfb <- (coverage(h2az.tgfb.rep1) * h2az.tgfb.rep1.scale + coverage(h2az.tgfb.rep2) * h2az.tgfb.rep2.scale) / 2

#----------Input data---read coverage---------------------------------------------------------
input.wt.rep1 <- readGAlignmentPairs("Input-TGFb-rep1_S7_L001_sorted_MkDup.bam", param = SBParam)
input.wt.rep2 <- readGAlignmentPairs("Input-TGFb-rep2_S8_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep1 <- readGAlignmentPairs("Input-WT-rep1_S5_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep2 <- readGAlignmentPairs("Input-WT-rep2_S6_L001_sorted_MkDup.bam", param = SBParam)

# scaling factor for coverage for on-target reads [RPM]
input.wt.rep1.scale <- 1000000/length(input.wt.rep1)
input.wt.rep2.scale <- 1000000/length(input.wt.rep2)
input.tgfb.rep1.scale <- 1000000/length(input.tgfb.rep1)
input.tgfb.rep2.scale <- 1000000/length(input.tgfb.rep2)

# calculate coverage from both replicates
cov.input.wt <- (coverage(input.wt.rep1) * input.wt.rep1.scale + coverage(input.wt.rep2) * input.wt.rep2.scale) / 2
cov.input.tgfb <- (coverage(input.tgfb.rep1) * input.tgfb.rep1.scale + coverage(input.tgfb.rep2) * input.tgfb.rep2.scale) / 2

#----------H2AZ data---dyad coverage---------------------------------------------------------
# convert GAlingmentPairs to GRanges, i.e. turn them into fragments
gr.h2az.wt.rep1 <- as(h2az.wt.rep1, "GRanges")
gr.h2az.wt.rep2 <- as(h2az.wt.rep2, "GRanges")
gr.h2az.tgfb.rep1 <- as(h2az.tgfb.rep1, "GRanges")
gr.h2az.tgfb.rep2 <- as(h2az.tgfb.rep2, "GRanges")

# retain only "genuine" nucleosome, i.e. discard fragments > 200bp - there is considerable variability in fragment length distribution... 
gr.h2az.wt.rep1 <- gr.h2az.wt.rep1[which(width(gr.h2az.wt.rep1) <= 200)]
gr.h2az.wt.rep2 <- gr.h2az.wt.rep2[which(width(gr.h2az.wt.rep2) <= 200)]
gr.h2az.tgfb.rep2 <- gr.h2az.tgfb.rep1[which(width(gr.h2az.tgfb.rep1) <= 200)]
gr.h2az.tgfb.rep2 <- gr.h2az.tgfb.rep2[which(width(gr.h2az.tgfb.rep2) <= 200)]

# resetting fragments to mid-point of nucleosome
gr.h2az.wt.rep1 <- resize(gr.h2az.wt.rep1, width = 1, fix = "center")
gr.h2az.wt.rep2 <- resize(gr.h2az.wt.rep2, width = 1, fix = "center")
gr.h2az.tgfb.rep2 <- resize(gr.h2az.tgfb.rep1, width = 1, fix = "center")
gr.h2az.tgfb.rep2 <- resize(gr.h2az.tgfb.rep2, width = 1, fix = "center")

#----------Input data---dyad coverage---------------------------------------------------------
# convert GAlingmentPairs to GRanges, i.e. turn them into fragments
gr.input.wt.rep1 <- as(input.wt.rep1, "GRanges")
gr.input.wt.rep2 <- as(input.wt.rep2, "GRanges")
gr.input.tgfb.rep1 <- as(input.tgfb.rep1, "GRanges")
gr.input.tgfb.rep2 <- as(input.tgfb.rep2, "GRanges")

# retain only "genuine" nucleosome, i.e. discard fragments > 160bp
gr.input.wt.rep1 <- gr.input.wt.rep1[which(width(gr.input.wt.rep1) <= 200)]
gr.input.wt.rep2 <- gr.input.wt.rep2[which(width(gr.input.wt.rep2) <= 200)]
gr.input.tgfb.rep2 <- gr.input.tgfb.rep1[which(width(gr.input.tgfb.rep1) <= 200)]
gr.input.tgfb.rep2 <- gr.input.tgfb.rep2[which(width(gr.input.tgfb.rep2) <= 200)]

# resetting fragments to mid-point of nucleosome
gr.input.wt.rep1 <- resize(gr.input.wt.rep1, width = 1, fix = "center")
gr.input.wt.rep2 <- resize(gr.input.wt.rep2, width = 1, fix = "center")
gr.input.tgfb.rep2 <- resize(gr.input.tgfb.rep1, width = 1, fix = "center")
gr.input.tgfb.rep2 <- resize(gr.input.tgfb.rep2, width = 1, fix = "center")

#----------preparing data for plotting of sequencing coverage--------------
# creating single-bp level windows for visualization
gr.which.tiles <- tile(gr.which, width = 1) # width = bin size
gr.which.tiles <- unlist(gr.which.tiles)
# if(length(grep("chr", seqlevels(gr.which.tiles))) == 0){
#   seqlevels(gr.which.tiles, force = T) <- paste("chr", seqlevels(gr.which.tiles), sep = "")
# }
seqlevels(gr.which.tiles, force = T) <- seqlevels(gr.which.tiles)[order(seqlevels(gr.which.tiles))]

cov.input.wt <- cov.input.wt[which(names(cov.input.wt) %in% seqlevels(gr.which.tiles))]
cov.input.wt <- cov.input.wt[names(cov.input.wt)[order(names(cov.input.wt))]]

cov.input.tgfb <- cov.input.tgfb[which(names(cov.input.tgfb) %in% seqlevels(gr.which.tiles))]
cov.input.tgfb <- cov.input.tgfb[names(cov.input.tgfb)[order(names(cov.input.tgfb))]]

cov.h2az.wt <- cov.h2az.wt[which(names(cov.h2az.wt) %in% seqlevels(gr.which.tiles))]
cov.h2az.wt <- cov.h2az.wt[names(cov.h2az.wt)[order(names(cov.h2az.wt))]]

cov.h2az.tgfb <- cov.h2az.tgfb[which(names(cov.h2az.tgfb) %in% seqlevels(gr.which.tiles))]
cov.h2az.tgfb <- cov.h2az.tgfb[names(cov.h2az.tgfb)[order(names(cov.h2az.tgfb))]]

# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.input.wt <- binnedAverage(gr.which.tiles, cov.input.wt, "mean")
bA.cov.input.tgfb <- binnedAverage(gr.which.tiles, cov.input.tgfb, "mean")
bA.cov.h2az.wt <- binnedAverage(gr.which.tiles, cov.h2az.wt, "mean")
bA.cov.h2az.tgfb <- binnedAverage(gr.which.tiles, cov.h2az.tgfb, "mean")

#----------Using Gviz for visualization----------------------------------------
# reads across the whole capture region (40kb)
# # sequencing coverage
# seqlevels(bA.cov.input.wt, force = T) <- gsub("chr", "", seqlevels(bA.cov.input.wt))
# seqlevels(bA.cov.input.tgfb, force = T) <- gsub("chr", "", seqlevels(bA.cov.input.tgfb))
# seqlevels(bA.cov.h2az.wt, force = T) <- gsub("chr", "", seqlevels(bA.cov.h2az.wt))
# seqlevels(bA.cov.h2az.tgfb, force = T) <- gsub("chr", "", seqlevels(bA.cov.h2az.tgfb))

dT.cov.input.wt <- DataTrack(bA.cov.input.wt, type = "h", col = "darkgreen", name = "Input WT [rpm]")
dT.cov.input.tgfb <- DataTrack(bA.cov.input.tgfb, type = "h", col = "lightgreen", name = "Input TGFb [rpm]")
dT.cov.h2az.wt <- DataTrack(bA.cov.h2az.wt, type = "h", col = "darkred", name = "H2AZ WT [rpm]")
dT.cov.h2az.tgfb <- DataTrack(bA.cov.h2az.tgfb, type = "h", col = "red", name = "H2AZ TGFb [rpm]")

# dT.cov.input.wt.tss <- DataTrack(bA.cov.input.wt.tss, type = "h", col = "darkgreen", name = "Input WT TSS1500 [rpm]")
# dT.cov.input.tgfb.tss <- DataTrack(bA.cov.input.tgfb.tss, type = "h", col = "lightgreen", name = "Input TGFb TSS1500 [rpm]")
# dT.cov.h2az.wt.tss <- DataTrack(bA.cov.h2az.wt.tss, type = "h", col = "darkred", name = "H2AZ WT TSS1500 [rpm]")
# dT.cov.h2az.tgfb.tss <- DataTrack(bA.cov.h2az.tgfb.tss, type = "h", col = "red", name = "H2AZ TGFb TSS1500 [rpm]")

displayPars(dT.cov.input.wt) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.input.tgfb) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.wt) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.tgfb) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
# displayPars(dT.cov.input.wt.tss) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
# displayPars(dT.cov.input.tgfb.tss) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
# displayPars(dT.cov.h2az.wt.tss) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
# displayPars(dT.cov.h2az.tgfb.tss) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

displayPars(aT.primers) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(aT.captureProbes) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

for(i in 1:nrow(designTab)){
  if (designTab[i,]$hgnc_symbol %in% subsetByOverlaps(gr.which.tss.nostrand, gr.which)$hgnc_symbol){
    chromosome(dT.cov.input.wt) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.input.tgfb) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.wt) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.tgfb) <- as(designTab[i, "chromosome_name"], "integer")
    
    chromosome(dT.cov.input.wt.tss) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.input.tgfb.tss) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.wt.tss) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.tgfb.tss) <- as(designTab[i, "chromosome_name"], "integer")
    
    chromosome(aT.primers) <- as(designTab[i, "chromosome_name"], "integer")
    chromosome(aT.captureProbes) <- as(designTab[i, "chromosome_name"], "integer")
    
    max.y <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
    displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y))
    displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y))
    displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y))
    displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y))
    
    max.y.tss <- max(max(values(dT.cov.input.wt.tss)), max(values(dT.cov.input.tgfb.tss)), max(values(dT.cov.h2az.wt.tss)), max(values(dT.cov.h2az.tgfb.tss)))
    displayPars(dT.cov.input.wt.tss) <- list(ylim = c(0,max.y.tss))
    displayPars(dT.cov.input.tgfb.tss) <- list(ylim = c(0,max.y.tss))
    displayPars(dT.cov.h2az.wt.tss) <- list(ylim = c(0,max.y.tss))
    displayPars(dT.cov.h2az.tgfb.tss) <- list(ylim = c(0,max.y.tss))
  
    # Annotation of the gene region
    biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", chromosome = as(designTab[i, "chromosome_name"], "integer"), start = as.integer(designTab[i, "promoter_start_position"]), end = as.integer(designTab[i, "promoter_end_position"]), name = paste(designTab[i, "hgnc_symbol"], sep = ""), mart = dog)
    displayPars(biomTrack) <- list(showFeatureId = TRUE, showId = TRUE, "fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
    
    # plotting at two different levels of resolution (full capture region & TSS)
    ncols <- 1
    nrows <- 2
    pdf(file = paste(as(designTab[i,"marker"], "character"), "_", as(designTab[i,"hgnc_symbol"], "character"), "combined", ".pdf", sep = ""), width = 12, height = 18)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrows, ncols)))
    
    axisTrack <- GenomeAxisTrack()
    pushViewport(viewport(layout.pos.col = ncols, layout.pos.row = 1))
    plotTracks(main = paste("Coverage ", as(designTab[i,"marker"], "character"), " ", as(designTab[i,"hgnc_symbol"], "character"), " 40kb", sep = ""),
               list(axisTrack, 
                    biomTrack,
                    aT.primers,
                    aT.captureProbes,
                    dT.cov.input.wt,
                    dT.cov.h2az.wt,
                    dT.cov.input.tgfb,
                    dT.cov.h2az.tgfb
                    ), 
               chromosome = designTab[i, "chromosome_name"], 
               from = as(designTab[i, "start_position"], "integer"), 
               to = as(designTab[i, "end_position"], "integer"), 
               extend.left = 2500, 
               extend.right = 2500,
               add = TRUE, 
               littleTicks = TRUE, 
               scale = 0.5)
    popViewport(1)
    pushViewport(viewport(layout.pos.col = ncols, layout.pos.row = 2))
    plotTracks(main = paste("Coverage ", as(designTab[i,"marker"], "character"), " ", as(designTab[i,"hgnc_symbol"], "character"), " TSS1500", sep = ""),
               list(axisTrack, 
                    biomTrack,
                    aT.primers,
                    dT.cov.input.wt.tss,
                    dT.cov.h2az.wt.tss,
                    dT.cov.input.tgfb.tss,
                    dT.cov.h2az.tgfb.tss
                    ), 
               chromosome = designTab[i, "chromosome_name"], 
               from = as(designTab[i, "promoter_start_position"], "integer"), 
               to = as(designTab[i, "promoter_end_position"], "integer"), 
               extend.left = 2500, 
               extend.right = 2500,
               add = TRUE, 
               littleTicks = TRUE, 
               scale = 0.5)
    popViewport(1)
    dev.off()
  } else {
    print("Wrong TSS")
  }
}

#----------H2AZ ChIP data across the complete CDSs-----------------------------------------------
# load ChIP coverage from BAM files
# if (length(grep("chr", seqlevels(gr.which))) > 0){
#   seqlevels(gr.which, force = T) <- gsub("chr", "", seqlevels(gr.which))
# }
flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = NA)
SBParam <- ScanBamParam(flag = flag, simpleCigar = F, what = scanBamWhat(), which = gr.which.cds)

h2az.wt.rep1.cds <- readGAlignmentPairs("H2AZ-WT-rep1_S1_L001_sorted_MkDup.bam", param = SBParam)
h2az.wt.rep2.cds <- readGAlignmentPairs("H2AZ-WT-rep2_S2_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep1.cds <- readGAlignmentPairs("H2AZ-TGFb-rep1_S3_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep2.cds <- readGAlignmentPairs("H2AZ-TGFb-rep2_S4_L001_sorted_MkDup.bam", param = SBParam)

# scaling factor for coverage across on-target regions [RPM]
h2az.wt.rep1.scale.cds <- 1000000/length(h2az.wt.rep1.cds)
h2az.wt.rep2.scale.cds <- 1000000/length(h2az.wt.rep2.cds)
h2az.tgfb.rep1.scale.cds <- 1000000/length(h2az.tgfb.rep1.cds)
h2az.tgfb.rep2.scale.cds <- 1000000/length(h2az.tgfb.rep2.cds)

# calculate coverage from both replicates
cov.h2az.wt.cds <- (coverage(h2az.wt.rep1.cds) * h2az.wt.rep1.scale.cds + coverage(h2az.wt.rep2.cds) * h2az.wt.rep2.scale.cds) / 2
cov.h2az.tgfb.cds <- (coverage(h2az.tgfb.rep1.cds) * h2az.tgfb.rep1.scale.cds + coverage(h2az.tgfb.rep2.cds) * h2az.tgfb.rep2.scale.cds) / 2

#----------Input data---read coverage across the complete CDSs-----------------------------------------
input.wt.rep1.cds <- readGAlignmentPairs("Input-TGFb-rep1_S7_L001_sorted_MkDup.bam", param = SBParam)
input.wt.rep2.cds <- readGAlignmentPairs("Input-TGFb-rep2_S8_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep1.cds <- readGAlignmentPairs("Input-WT-rep1_S5_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep2.cds <- readGAlignmentPairs("Input-WT-rep2_S6_L001_sorted_MkDup.bam", param = SBParam)

# scaling factor for coverage for on-target reads [RPM]
input.wt.rep1.scale.cds <- 1000000/length(input.wt.rep1.cds)
input.wt.rep2.scale.cds <- 1000000/length(input.wt.rep2.cds)
input.tgfb.rep1.scale.cds <- 1000000/length(input.tgfb.rep1.cds)
input.tgfb.rep2.scale.cds <- 1000000/length(input.tgfb.rep2.cds)

# calculate coverage from both replicates
cov.input.wt.cds <- (coverage(input.wt.rep1.cds) * input.wt.rep1.scale.cds + coverage(input.wt.rep2.cds) * input.wt.rep2.scale.cds) / 2
cov.input.tgfb.cds <- (coverage(input.tgfb.rep1.cds) * input.tgfb.rep1.scale.cds + coverage(input.tgfb.rep2.cds) * input.tgfb.rep2.scale.cds) / 2

#----------Table with start/end position of CDS for plotting------------------------------------------
cdsTab <- rbind(mesenchymalMarkers.tab, epithelialMarkers.tab)
cdsTab$promoter_start_position <- start(gr.which.tss)
cdsTab$promoter_end_position <- end(gr.which.tss)

#----------preparing data for plotting of CDS sequencing coverage-------------------------------------
# creating single-bp level windows for visualization
gr.which.tiles <- tile(gr.which, width = 1) # width = bin size
gr.which.tiles <- unlist(gr.which.tiles)
# if(length(grep("chr", seqlevels(gr.which.tiles))) == 0){
#   seqlevels(gr.which.tiles, force = T) <- paste("chr", seqlevels(gr.which.tiles), sep = "")
# }
seqlevels(gr.which.tiles, force = T) <- seqlevels(gr.which.tiles)[order(seqlevels(gr.which.tiles))]

cov.input.wt.cds <- cov.input.wt.cds[which(names(cov.input.wt.cds) %in% seqlevels(gr.which.tiles))]
cov.input.wt.cds <- cov.input.wt.cds[names(cov.input.wt.cds)[order(names(cov.input.wt.cds))]]

cov.input.tgfb.cds <- cov.input.tgfb.cds[which(names(cov.input.tgfb.cds) %in% seqlevels(gr.which.tiles))]
cov.input.tgfb.cds <- cov.input.tgfb.cds[names(cov.input.tgfb.cds)[order(names(cov.input.tgfb.cds))]]

cov.h2az.wt.cds <- cov.h2az.wt.cds[which(names(cov.h2az.wt.cds) %in% seqlevels(gr.which.tiles))]
cov.h2az.wt.cds <- cov.h2az.wt.cds[names(cov.h2az.wt.cds)[order(names(cov.h2az.wt.cds))]]

cov.h2az.tgfb.cds <- cov.h2az.tgfb.cds[which(names(cov.h2az.tgfb.cds) %in% seqlevels(gr.which.tiles))]
cov.h2az.tgfb.cds <- cov.h2az.tgfb.cds[names(cov.h2az.tgfb.cds)[order(names(cov.h2az.tgfb.cds))]]

# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.input.wt.cds <- binnedAverage(gr.which.tiles, cov.input.wt.cds, "mean")
bA.cov.input.tgfb.cds <- binnedAverage(gr.which.tiles, cov.input.tgfb.cds, "mean")
bA.cov.h2az.wt.cds <- binnedAverage(gr.which.tiles, cov.h2az.wt.cds, "mean")
bA.cov.h2az.tgfb.cds <- binnedAverage(gr.which.tiles, cov.h2az.tgfb.cds, "mean")

bA.cov.input.wt.tss.cds <- subsetByOverlaps(bA.cov.input.wt.cds, gr.which.tss.nostrand)
bA.cov.input.tgfb.tss.cds <- subsetByOverlaps(bA.cov.input.tgfb.cds, gr.which.tss.nostrand)
bA.cov.h2az.wt.tss.cds <- subsetByOverlaps(bA.cov.h2az.wt.cds, gr.which.tss.nostrand)
bA.cov.h2az.tgfb.tss.cds <- subsetByOverlaps(bA.cov.h2az.tgfb.cds, gr.which.tss.nostrand)

#----------plotting reads across the complete CDSs-----------------------------------

dT.cov.input.wt.cds <- DataTrack(bA.cov.input.wt.cds, type = "h", col = "darkgreen", name = "Input WT [rpm]")
dT.cov.input.tgfb.cds <- DataTrack(bA.cov.input.tgfb.cds, type = "h", col = "lightgreen", name = "Input TGFb [rpm]")
dT.cov.h2az.wt.cds <- DataTrack(bA.cov.h2az.wt.cds, type = "h", col = "darkred", name = "H2AZ WT [rpm]")
dT.cov.h2az.tgfb.cds <- DataTrack(bA.cov.h2az.tgfb.cds, type = "h", col = "red", name = "H2AZ TGFb [rpm]")

dT.cov.input.wt.tss.cds <- DataTrack(bA.cov.input.wt.tss.cds, type = "h", col = "darkgreen", name = "Input WT TSS1500 [rpm]")
dT.cov.input.tgfb.tss.cds <- DataTrack(bA.cov.input.tgfb.tss.cds, type = "h", col = "lightgreen", name = "Input TGFb TSS1500 [rpm]")
dT.cov.h2az.wt.tss.cds <- DataTrack(bA.cov.h2az.wt.tss.cds, type = "h", col = "darkred", name = "H2AZ WT TSS1500 [rpm]")
dT.cov.h2az.tgfb.tss.cds <- DataTrack(bA.cov.h2az.tgfb.tss.cds, type = "h", col = "red", name = "H2AZ TGFb TSS1500 [rpm]")

displayPars(dT.cov.input.wt.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.input.tgfb.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.wt.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.tgfb.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.input.wt.tss.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.input.tgfb.tss.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.wt.tss.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(dT.cov.h2az.tgfb.tss.cds) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

displayPars(aT.primers) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(aT.captureProbes) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

for(i in 1:nrow(cdsTab)){
  if (cdsTab[i,]$hgnc_symbol %in% subsetByOverlaps(gr.which.tss.nostrand, gr.which)$hgnc_symbol){
    chromosome(dT.cov.input.wt.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.input.tgfb.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.wt.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.tgfb.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    
    chromosome(dT.cov.input.wt.tss.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.input.tgfb.tss.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.wt.tss.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(dT.cov.h2az.tgfb.tss.cds) <- as(cdsTab[i, "chromosome_name"], "integer")
    
    chromosome(aT.primers) <- as(cdsTab[i, "chromosome_name"], "integer")
    chromosome(aT.captureProbes) <- as(cdsTab[i, "chromosome_name"], "integer")
    
    max.y.cds <- max(max(values(dT.cov.input.wt.cds)), max(values(dT.cov.input.tgfb.cds)), max(values(dT.cov.h2az.wt.cds)), max(values(dT.cov.h2az.tgfb.cds)))
    displayPars(dT.cov.input.wt.cds) <- list(ylim = c(0,max.y.cds))
    displayPars(dT.cov.input.tgfb.cds) <- list(ylim = c(0,max.y.cds))
    displayPars(dT.cov.h2az.wt.cds) <- list(ylim = c(0,max.y.cds))
    displayPars(dT.cov.h2az.tgfb.cds) <- list(ylim = c(0,max.y.cds))
    
    max.y.tss.cds <- max(max(values(dT.cov.input.wt.tss.cds)), max(values(dT.cov.input.tgfb.tss.cds)), max(values(dT.cov.h2az.wt.tss.cds)), max(values(dT.cov.h2az.tgfb.tss.cds)))
    displayPars(dT.cov.input.wt.tss.cds) <- list(ylim = c(0,max.y.tss.cds))
    displayPars(dT.cov.input.tgfb.tss.cds) <- list(ylim = c(0,max.y.tss.cds))
    displayPars(dT.cov.h2az.wt.tss.cds) <- list(ylim = c(0,max.y.tss.cds))
    displayPars(dT.cov.h2az.tgfb.tss.cds) <- list(ylim = c(0,max.y.tss.cds))
    
    # Annotation of the gene region
    biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", chromosome = as(cdsTab[i, "chromosome_name"], "integer"), start = as.integer(cdsTab[i, "promoter_start_position"]), end = as.integer(cdsTab[i, "promoter_end_position"]), name = paste(cdsTab[i, "hgnc_symbol"], sep = ""), mart = dog)
    displayPars(biomTrack) <- list(showFeatureId = TRUE, showId = TRUE, "fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
    
    # plotting at two different levels of resolution (full capture region & TSS)
    ncols <- 1
    nrows <- 2
    pdf(file = paste(as(cdsTab[i,"marker"], "character"), "_", as(cdsTab[i,"hgnc_symbol"], "character"), "_combined_CDS", ".pdf", sep = ""), width = 12, height = 18)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrows, ncols)))
    
    axisTrack <- GenomeAxisTrack()
    pushViewport(viewport(layout.pos.col = ncols, layout.pos.row = 1))
    plotTracks(main = paste("Coverage ", as(cdsTab[i,"marker"], "character"), " ", as(cdsTab[i,"hgnc_symbol"], "character"), " CDS", sep = ""),
               list(axisTrack, 
                    biomTrack,
                    aT.primers,
                    aT.captureProbes,
                    dT.cov.input.wt.cds,
                    dT.cov.h2az.wt.cds,
                    dT.cov.input.tgfb.cds,
                    dT.cov.h2az.tgfb.cds
               ), 
               chromosome = cdsTab[i, "chromosome_name"], 
               from = as(cdsTab[i, "start_position"], "integer"), 
               to = as(cdsTab[i, "end_position"], "integer"), 
               extend.left = 2500, 
               extend.right = 2500,
               add = TRUE, 
               littleTicks = TRUE, 
               scale = 0.5)
    popViewport(1)
    pushViewport(viewport(layout.pos.col = ncols, layout.pos.row = 2))
    plotTracks(main = paste("Coverage ", as(cdsTab[i,"marker"], "character"), " ", as(cdsTab[i,"hgnc_symbol"], "character"), " TSS1500", sep = ""),
               list(axisTrack, 
                    biomTrack,
                    aT.primers,
                    dT.cov.input.wt.tss.cds,
                    dT.cov.h2az.wt.tss.cds,
                    dT.cov.input.tgfb.tss.cds,
                    dT.cov.h2az.tgfb.tss.cds
               ), 
               chromosome = cdsTab[i, "chromosome_name"], 
               from = as(cdsTab[i, "promoter_start_position"], "integer"), 
               to = as(cdsTab[i, "promoter_end_position"], "integer"), 
               extend.left = 2500, 
               extend.right = 2500,
               add = TRUE, 
               littleTicks = TRUE, 
               scale = 0.5)
    popViewport(1)
    dev.off()
  } else {
    print("Wrong TSS")
  }
}
