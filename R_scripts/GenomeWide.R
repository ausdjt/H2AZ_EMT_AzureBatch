# load libraries
library("gdata")
library("Gviz")
library("rtracklayer")
library("biomaRt")
library("GenomicRanges")
library("GenomicAlignments")
library(BSgenome.Cfamiliaris.UCSC.canFam3)

options(ucscChromosomeNames=FALSE)

setwd('~/Data/Tremethick/EMT/')

source("~/Development/GeneralPurpose/R/binnedAverage.R")
source("~/Development/GeneralPurpose/R/binnedSum.R")
source("~/Development/JCSMR_Genomics/R/TremethickLab/H2AZ_EMT/calculateCoverage.R")

# define Ensembl IDs prior to lookup at Biomart:
# # EMT markers
# # mesenchymal genes HGNC Symbols
# FN1 -> ENSCAFG00000014345
# ZEB1 -> ENSCAFG00000004023
# TGFb1 -> ENSCAFG00000005014
# SPARC -> ENSCAFG00000017855
# TWIST2 -> ENSCAFG00000012469
# # epithelial markers
# CDH1 -> E-cadherin -> ENSCAFG00000020305 [CDH3 -> ENSCAFG00000020303]
# SPP1 -> ENSCAFG00000009569
# FGFBP1 -> ENSCAFG00000015272
# MMP9 -> ENSCAFG00000009905
# EPCAM -> ENSCAFG00000002653

#----------connect to Ensembl biomaRt for annotation data----------------------
ensembl <- useEnsembl(biomart = "ensembl", host = "asia.ensembl.org")
dog <- useEnsembl(biomart = "ensembl", dataset = "cfamiliaris_gene_ensembl", host = "asia.ensembl.org")
human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "asia.ensembl.org")

filters <- listFilters(dog)
att <- listAttributes(dog)
human.attributes <- listAttributes(human)

#---------------Mesenchymal Markers---------------------------------------------------
mesenchymalMarkers<- c("FN1" = "ENSCAFG00000014345", 
                       "ZEB1" = "ENSCAFG00000004023", 
                       "TGFb1" = "ENSCAFG00000005014",
                       "TGFb3" = "ENSCAFG00000017101",
                       "SPARC" = "ENSCAFG00000017855", 
                       "TWIST2" = "ENSCAFG00000012469"
                       )

mesenchymalMarkers.transcripts.tab <- getBM(attributes = c("ensembl_gene_id",
                                                           "ensembl_transcript_id",
                                                           "chromosome_name",
                                                           "start_position",
                                                           "end_position",
                                                           "transcript_start",
                                                           "hgnc_symbol",
                                                           "external_transcript_name",
                                                           "strand"),
                                            filters = "ensembl_gene_id",
                                            values = mesenchymalMarkers, dog)

mesenchymalMarkers.transcripts.tab$marker <- "mesenchymal"

mesenchymalMarkers.genes.tab <- getBM(attributes = c("ensembl_gene_id",
                                                     "entrezgene",
                                                     "chromosome_name",
                                                     "start_position",
                                                     "end_position",
                                                     "hgnc_symbol",
                                                     "description",
                                                     "strand"),
                                            filters = "ensembl_gene_id",
                                            values = mesenchymalMarkers, dog)
mesenchymalMarkers.genes.tab$marker <- "mesenchymal"

# get human homologs of these dog genes
mesenchymalMarkers.Hsap.homologs <- getBM(attributes = c("ensembl_gene_id",
                                                        "hsapiens_homolog_ensembl_gene"),
                                         filters = "ensembl_gene_id",
                                         values = mesenchymalMarkers, dog)
mesenchymalMarkers.Hsap.homologs.EntrezIDs <- getBM(attributes = c("ensembl_gene_id",
                                                                  "entrezgene",
                                                                  "hgnc_symbol"),
                                                   filters = "ensembl_gene_id",
                                                   values = mesenchymalMarkers.Hsap.homologs$hsapiens_homolog_ensembl_gene, human)
mesenchymalMarkers.Hsap.homologs.EntrezIDs$marker <- "mesenchymal"

gr.mesenchymalMarkers <- GRanges(seqnames = mesenchymalMarkers.transcripts.tab$chromosome_name, 
                                IRanges(mesenchymalMarkers.transcripts.tab$start_position, 
                                        mesenchymalMarkers.transcripts.tab$end_position), 
                                strand = c("-", "+")[match(mesenchymalMarkers.transcripts.tab$strand, c("-1", "1"))],
                                mesenchymalMarkers.transcripts.tab[,c("ensembl_gene_id", 
                                                                      "ensembl_transcript_id", 
                                                                      "hgnc_symbol", 
                                                                      "transcript_start", 
                                                                      "external_transcript_name", 
                                                                      "marker")])

gr.mesenchymalMarkers$hgnc_symbol[3:5] <- paste(gr.mesenchymalMarkers$hgnc_symbol[3:5], "(TGFB-1)", sep = " ")

gr.mesenchymalMarkers.genes <-  GRanges(seqnames = mesenchymalMarkers.genes.tab$chromosome_name, 
                                        IRanges(mesenchymalMarkers.genes.tab$start_position, 
                                                mesenchymalMarkers.genes.tab$end_position), 
                                        strand = c("-", "+")[match(mesenchymalMarkers.genes.tab$strand, c("-1", "1"))],
                                        mesenchymalMarkers.genes.tab[,c("ensembl_gene_id", 
                                                                        "entrezgene", 
                                                                        "hgnc_symbol", 
                                                                        "description", 
                                                                        "marker")])
gr.mesenchymalMarkers.genes$hgnc_symbol[2:3] <- paste(gr.mesenchymalMarkers.genes$hgnc_symbol[2:3], "(TGFB-1)", sep = " ")


#---------------Epithelial Markers---------------------------------------------------
epithelialMarkers <- c("CDH1" = "ENSCAFG00000020305",
                       "SPP1" = "ENSCAFG00000009569",
                       "FGFBP1" = "ENSCAFG00000015272",
                       "MMP9" = "ENSCAFG00000009905",
                       "EPCAM" = "ENSCAFG00000002653")

epithelialMarkers.transcripts.tab <- getBM(attributes = c("ensembl_gene_id",
                                                          "ensembl_transcript_id",
                                                          "entrezgene",
                                                          "chromosome_name",
                                                          "start_position",
                                                          "end_position",
                                                          "transcript_start",
                                                          "hgnc_symbol",
                                                          "strand",
                                                          "external_transcript_name"),
                                           filters = "ensembl_gene_id",
                                           values = epithelialMarkers, dog)
epithelialMarkers.transcripts.tab$marker <- "epithelial"

epithelialMarkers.genes.tab <- getBM(attributes = c("ensembl_gene_id",
                                                    "entrezgene",
                                                    "chromosome_name",
                                                    "start_position",
                                                    "end_position",
                                                    "hgnc_symbol",
                                                    "strand"),
                                     filters = "ensembl_gene_id",
                                     values = epithelialMarkers, dog)
epithelialMarkers.genes.tab$marker <- "epithelial"

# get human homologs of these dog genes
epithelialMarkers.Hsap.homologs <- getBM(attributes = c("ensembl_gene_id",
                                                        "hsapiens_homolog_ensembl_gene"),
                                         filters = "ensembl_gene_id",
                                         values = epithelialMarkers, dog)
epithelialMarkers.Hsap.homologs.EntrezIDs <- getBM(attributes = c("ensembl_gene_id",
                                                                  "entrezgene",
                                                                  "hgnc_symbol"),
                                                   filters = "ensembl_gene_id",
                                                   values = epithelialMarkers.Hsap.homologs$hsapiens_homolog_ensembl_gene, human)
epithelialMarkers.Hsap.homologs.EntrezIDs$marker <- "epithelial"

gr.epithelialMarkers <- GRanges(seqnames = epithelialMarkers.transcripts.tab$chromosome_name, 
                                IRanges(epithelialMarkers.transcripts.tab$start_position, 
                                        epithelialMarkers.transcripts.tab$end_position), 
                                strand = c("-", "+")[match(epithelialMarkers.transcripts.tab$strand, c("-1", "1"))],
                                epithelialMarkers.transcripts.tab[,c("ensembl_gene_id", 
                                                                     "ensembl_transcript_id", 
                                                                     "hgnc_symbol", 
                                                                     "transcript_start", 
                                                                     "marker",
                                                                     "external_transcript_name")])

gr.epithelialMarkers.genes <-  GRanges(seqnames = epithelialMarkers.genes.tab$chromosome_name, 
                                        IRanges(epithelialMarkers.genes.tab$start_position, 
                                                epithelialMarkers.genes.tab$end_position), 
                                        strand = c("-", "+")[match(epithelialMarkers.genes.tab$strand, c("-1", "1"))],
                                        epithelialMarkers.genes.tab[,c("ensembl_gene_id", 
                                                                       "entrezgene",
                                                                       "hgnc_symbol", 
                                                                       "marker")])


#--------------add TGFb-------------------------------------------------------
# TGFb.Hsap <- "ENSG00000105329"
# TGFb.Hsap <- getBM(attributes = c("ensembl_gene_id",
#                                   "entrezgene",
#                                   "hgnc_symbol"),
#                    filters = "ensembl_gene_id",
#                    values = TGFb.Hsap, human)
# 
# getBM(attributes = c("ensembl_gene_id",
#                      "cfamiliaris_homolog_ensembl_gene"),
#       filters = "ensembl_gene_id",
#       values = TGFb.Hsap$ensembl_gene_id, human)

#-------------put together GRanges object for reading BAM files---------------
gr.which <- c(promoters(reduce(gr.mesenchymalMarkers), upstream = 400000, downstream = 400000),
              promoters(reduce(gr.epithelialMarkers), upstream = 400000, downstream = 400000))

# setting chromosome info
seqlevels(gr.which) <- paste("chr", seqlevels(gr.which), sep = "")
seqinfo(gr.which, force = T) <- seqinfo(BSgenome.Cfamiliaris.UCSC.canFam3)[seqlevels(gr.which)]
seqlevels(gr.which) <- gsub("chr", "", seqlevels(gr.which))

# load whole chromosome
# start(gr.which) <- 1
# end(gr.which) <- as.vector(seqlengths(gr.which))

# extending CDH1 region
# start(gr.which[10]) <- start(gr.which[10]) - 150000
# strand(gr.which) <- "*"
# save(gr.which, file = "gr.which.rda")

#-----------------------importing data from BAM files-------------------------------------------
# working with BAM files
# on GDU cluster
#path = "/Volumes/MHS/researchdata/JCSMR/TremethickLab/Illumina_Sequencing/MDCK_ChIPSeq/"

#path.input = "/Volumes/gduserv/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/"
path.input = "~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/"
#path.h2az = "/Volumes/gduserv/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/"
path.h2az = "~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/"

suffix = ".Q10.sorted.MkDup.bam"

files.input = c("Input_TGFb_rep1_S7", "Input_TGFb_rep2_S8", "Input_WT_rep1_S5", "Input_WT_rep2_S6")
files.h2az = c("H2AZ_TGFb_rep1_S3", "H2AZ_TGFb_rep2_S4", "H2AZ_WT_rep1_S1", "H2AZ_WT_rep2_S2")
#files.input <- c("Input_TGFb", "Input_WT")
#files.h2az <- c("H2AZ_TGFb", "H2AZ_WT")

# parameters for reading in BAM files
# flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = F, isSecondaryAlignment = F)#, isFirstMateRead = T, isSecondMateRead = F)
flag <- scanBamFlag(isProperPair = T, isDuplicate = F)
SBParam.all <- ScanBamParam(flag = flag, simpleCigar = T, what = c("rname", "strand", "pos", "qwidth")) #
SBParam <- ScanBamParam(flag = flag, simpleCigar = T, what = c("rname", "strand", "pos", "qwidth")) #, which = gr.which)

counts.input <- lapply(files.input, function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  print(fn)
  countBam(fn, param = SBParam.all)
})
names(counts.input) <- files.input

counts.h2az <- lapply(files.h2az, function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  print(fn)
  countBam(fn, param = SBParam.all)
})
names(counts.h2az) <- files.h2az

# calculate scaling factors per sample (RPM)
scale.input <- lapply(counts.input, function(x) {1000000/x$records})
scale.h2az <- lapply(counts.h2az, function(x) {1000000/x$records})

# import reads
# and calculate coverage
# and scale it according to total library size
# Input
coverage.input <- lapply(files.input, function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  fn.temp <- paste(fn, "temp", sep = "")
  print(fn)
  filterBam(fn, fn.temp , param = ScanBamParam(what = c("rname", "strand", "pos", "qwidth"), which = gr.which))
  rap <- readGAlignmentPairs(fn.temp, use.names = TRUE)
  rap <- as(rap, "GRanges")
  file.remove(fn.temp)
  file.remove(paste(fn.temp, "bai", sep = "."))
  coverage(rap) * scale.input[[x]]
})
names(coverage.input) <- files.input

# H2AZ ChIP
coverage.h2az <- lapply(files.h2az, function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  fn.temp <- paste(fn, "temp", sep = "")
  print(fn)
  filterBam(fn, fn.temp , param = ScanBamParam(what = c("rname", "strand", "pos", "qwidth"), which = gr.which))
  rap <- readGAlignmentPairs(fn.temp, use.names = TRUE)
  rap <- as(rap, "GRanges")
  file.remove(fn.temp)
  file.remove(paste(fn.temp, "bai", sep = "."))
  coverage(rap) * scale.h2az[[x]]
})
names(coverage.h2az) <- files.h2az

# calculate coverage from replicates, applying scaling factor (RPM)
cov.input.emt_markers.wt <- (coverage.input[["Input_WT_rep1_S5"]] + coverage.input[["Input_WT_rep2_S6"]]) /2
cov.input.emt_markers.tgfb <- (coverage.input[["Input_TGFb_rep1_S7"]] + coverage.input[["Input_TGFb_rep2_S8"]]) /2
cov.h2az.emt_markers.wt <- (coverage.h2az[["H2AZ_WT_rep1_S1"]] + coverage.h2az[["H2AZ_WT_rep2_S2"]]) / 2
cov.h2az.emt_markers.tgfb <- (coverage.h2az[["H2AZ_TGFb_rep1_S3"]] + coverage.h2az[["H2AZ_TGFb_rep2_S4"]]) / 2

#----------preparing data for plotting of sequencing coverage--------------
# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.input.emt_markers.wt <- calculateCoverage(step = 1, gr.which, cov.input.emt_markers.wt, func = "mean")
bA.cov.input.emt_markers.tgfb <- calculateCoverage(step = 1, gr.which, cov.input.emt_markers.tgfb, func = "mean")
bA.cov.h2az.emt_markers.wt <- calculateCoverage(step = 1, gr.which, cov.h2az.emt_markers.wt, func = "mean")
bA.cov.h2az.emt_markers.tgfb <- calculateCoverage(step = 1, gr.which, cov.h2az.emt_markers.tgfb, func = "mean")

#----------creating DataTrack objects for visualization using Gviz------------------------------
dT.cov.input.emt_markers.wt <- DataTrack(bA.cov.input.emt_markers.wt, type = "h", col = "darkgreen", name = "Input WT\n[rpm]", strand = "*")
dT.cov.input.emt_markers.tgfb <- DataTrack(bA.cov.input.emt_markers.tgfb, type = "h", col = "lightgreen", name = "Input TGFb\n[rpm]", strand = "*")
dT.cov.h2az.emt_markers.wt <- DataTrack(bA.cov.h2az.emt_markers.wt, type = "h", col = "darkred", name = "H2AZ WT\n[rpm]", strand = "*")
dT.cov.h2az.emt_markers.tgfb <- DataTrack(bA.cov.h2az.emt_markers.tgfb, type = "h", col = "red", name = "H2AZ TGFb\n[rpm]", strand = "*")

dpList <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white", "cex.title" = 0.4, rotation.title = 270, cex.axis = 0.6)

displayPars(dT.cov.input.emt_markers.wt) <- dpList
displayPars(dT.cov.input.emt_markers.tgfb) <- dpList
displayPars(dT.cov.h2az.emt_markers.wt) <- dpList
displayPars(dT.cov.h2az.emt_markers.tgfb) <- dpList

save(file = "genomeWide.50kbTSS.DataTracks.rda", list = c("dT.cov.input.emt_markers.wt", "dT.cov.input.emt_markers.tgfb", "dT.cov.h2az.emt_markers.wt", "dT.cov.h2az.emt_markers.tgfb"))

#-----------plotting coverage across epithelial markers------------------------------------------

dT.bgSubTGFb <- DataTrack("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/pooled/TGFb_ChIP.bgsub.Fnor.smooth.bw", stream = T, name = "TGFb - background subtracted", col = "black", type = "l")
dT.bgSubWT <- DataTrack("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/pooled/WT_ChIP.bgsub.Fnor.smooth.bw", stream = T, name = "WT - background subtracted", col = "grey", type = "l")
dT.Diff <- DataTrack("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/diff/TGFb_vs_WT.pois_diff.bw", stream = T, name = "Difference - DANPOS2", col = "blue", type = "l")

gr.plot <- promoters(gr.mesenchymalMarkers, up = 1500, down = 20000)
gr.plot <- promoters(gr.epithelialMarkers.genes, upstream = 20000, downstream = 20000)

biomTrack <- GeneRegionTrack(TxDb.Cfam3.Ensembl, showId = T, geneSymbol = T, showExonId = F, name = "", stacking = "hide")
displayPars(biomTrack) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
save(biomTrack, file = "biomTrack.Cfam3.Ensembl.rda")

pdf("~/OneDrive/Documents/ANU/Tremethick Lab/Lab Meetings/Lab Meeting 2015-10-14/MDCK_ChIP-Seq_EpitheliaMarkers_coverage_plots_1500TSS1500_incl_DMRs_incl_SureSelect.pdf")
for (i in 1:length(gr.plot)){
  displayPars(biomTrack) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
  
  # Data from SureSelect capture
#   chromosome(dT.cov.input.wt) <- seqnames(gr.plot)[i]
#   chromosome(dT.cov.input.tgfb) <- seqnames(gr.plot)[i]
#   chromosome(dT.cov.h2az.wt) <- seqnames(gr.plot)[i]
#   chromosome(dT.cov.h2az.tgfb) <- seqnames(gr.plot)[i]
#   
#   chromosome(aT.primers) <- seqnames(gr.plot)[i]
#   
#   max.y <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
#   displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y))
#   displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y))
#   displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y))
#   displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y))
#   
  # Data from whole genome ChIP-Seq
  
  chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.plot)[i]
  
  max.y.tss <- max(max(values(dT.cov.input.emt_markers.wt)), max(values(dT.cov.input.emt_markers.tgfb)), max(values(dT.cov.h2az.emt_markers.wt)), max(values(dT.cov.h2az.emt_markers.tgfb)))
  displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
  
  plotTracks(list(axisTrack,
                  biomTrack,
                  aT.ap1Sites,
                  aT.nfkbSites,
                  aT.TSS,
                  dT.cov.input.emt_markers.wt, 
                  dT.cov.h2az.emt_markers.wt,
                  dT.cov.input.emt_markers.tgfb, 
                  dT.cov.h2az.emt_markers.tgfb

  ),
  chromosome = as(seqnames(gr.plot), "character")[i],
  from = as.integer(start(gr.plot[i]), "integer"),
  to = as.integer(end(gr.plot[i]), "integer"),
  extend.right = 1000,
  extend.left = 1000,
  main = paste(gr.plot$hgnc_symbol[i], " (", width(gr.plot[i]), "bp)", sep = ""),
  strand = "*",
  cex.main = 0.5,
  sizes = c(0.02, 0.06, 0.04, 0.04, 0.04, 0.2, 0.2, 0.2, 0.2),
  scale = 0.5)

}
dev.off()

# mesenchymal markers
gr.mesenchymalMarkers.1500TSS1500 <- promoters(gr.mesenchymalMarkers, upstream = 1500, downstream = 1500)

pdf("~/OneDrive/Documents/ANU/Tremethick Lab/Lab Meetings/Lab Meeting 2015-10-14/MDCK_ChIP-Seq_MesenchymalMarkers_coverage_plots_1500TSS1500_incl_DMRs_incl_SureSelect.pdf")
for (i in 1:length(gr.mesenchymalMarkers.1500TSS1500)){
  biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", 
                                      chromosome = as(seqnames(gr.mesenchymalMarkers.1500TSS1500), "character")[i],
                                      start = as.integer(start(gr.mesenchymalMarkers.1500TSS1500[i]), "integer"),
                                      end = as.integer(end(gr.mesenchymalMarkers.1500TSS1500[i]), "integer"),
                                      name = paste(mcols(gr.mesenchymalMarkers.1500TSS1500[i])$hgnc_symbol, "transcript",  mcols(gr.mesenchymalMarkers.1500TSS1500[i])$ensembl_transcript_id, sep = " "),
                                      mart = dog)
  
  displayPars(biomTrack) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
  
  # Data from SureSelect capture
  chromosome(dT.cov.input.wt) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.input.tgfb) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.h2az.wt) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.h2az.tgfb) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  
  chromosome(aT.primers) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  
  max.y <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
  displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y))
  
  chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.mesenchymalMarkers.1500TSS1500)[i]
  
  max.y.tss <- max(max(values(dT.cov.input.emt_markers.wt)), max(values(dT.cov.input.emt_markers.tgfb)), max(values(dT.cov.h2az.emt_markers.wt)), max(values(dT.cov.h2az.emt_markers.tgfb)))
  displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
  
  
  plotTracks(list(axisTrack,
                  biomTrack,
                  aT.ap1Sites,
                  aT.nfkbSites,
                  aT.primers,
                  dT.cov.input.emt_markers.wt, 
                  dT.cov.h2az.emt_markers.wt,
                  dT.cov.input.wt,
                  dT.cov.h2az.wt,
                  dT.cov.input.emt_markers.tgfb, 
                  dT.cov.h2az.emt_markers.tgfb,
                  dT.cov.input.tgfb,
                  dT.cov.h2az.tgfb,
                  dT.dmr
  ),
  chromosome = as(seqnames(gr.mesenchymalMarkers.1500TSS1500), "character")[i],
  from = as.integer(start(gr.mesenchymalMarkers.1500TSS1500[i]), "integer"),
  to = as.integer(end(gr.mesenchymalMarkers.1500TSS1500[i]), "integer"),
  extend.right = 1000,
  extend.left = 1000,
  main = mcols(gr.mesenchymalMarkers.1500TSS1500[i])$hgnc_symbol,
  strand = "*",
  cex.main = 0.5,
  sizes = c(0.01, 0.04, 0.02, 0.02, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05),
  scale = 0.5)
 
}
dev.off()

 
# whole 50kb region
gr.mesenchymalMarkers.genes.25kbTSS25kb <- promoters(gr.mesenchymalMarkers.genes, upstream = 20000, downstream = 20000)

pdf("~/OneDrive/Documents/ANU/Tremethick Lab/Lab Meetings/Lab Meeting 2015-10-14/MDCK_ChIP-Seq_MesenchymalMarkers_coverage_plots_20kbTSS20kb_incl_DMRs_incl_SureSelect.pdf", height = 10, width = 15)
for (i in 1:length(gr.mesenchymalMarkers.genes.25kbTSS25kb)){
  biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", 
                                      chromosome = as(seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb), "character")[i],
                                      start = as.integer(start(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
                                      end = as.integer(end(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
                                      name = mcols(gr.mesenchymalMarkers.genes.25kbTSS25kb[i])$hgnc_symbol,
                                      mart = dog)
  displayPars(biomTrack) <- list(showFeatureId = TRUE, showId = TRUE, "fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
  
  # Data from SureSelect capture
  chromosome(dT.cov.input.wt) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.input.tgfb) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.h2az.wt) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.h2az.tgfb) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  
  chromosome(aT.primers) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  
  max.y <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
  displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y))
  
  chromosome(dT.cov.input.emt_markers.wt) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.input.emt_markers.tgfb) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.h2az.emt_markers.wt) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  chromosome(dT.cov.h2az.emt_markers.tgfb) <- seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb)[i]
  
  max.y.tss <- max(max(values(dT.cov.input.emt_markers.wt)), max(values(dT.cov.input.emt_markers.tgfb)), max(values(dT.cov.h2az.emt_markers.wt)), max(values(dT.cov.h2az.emt_markers.tgfb)))
  displayPars(dT.cov.input.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.input.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.emt_markers.tgfb) <- list(ylim = c(0,max.y.tss))

  plotTracks(list(axisTrack,
                  biomTrack,
                  aT.ap1Sites,
                  aT.nfkbSites,
                  aT.primers,
                  dT.cov.input.emt_markers.wt, 
                  dT.cov.h2az.emt_markers.wt,
                  dT.cov.input.wt,
                  dT.cov.h2az.wt,
                  dT.cov.input.emt_markers.tgfb, 
                  dT.cov.h2az.emt_markers.tgfb,
                  dT.cov.input.tgfb,
                  dT.cov.h2az.tgfb,
                  dT.dmr
                ),
             chromosome = as(seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb), "character")[i],
             from = as.integer(start(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
             to = as.integer(end(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
             extend.right = 1000,
             extend.left = 1000,
             main = mcols(gr.mesenchymalMarkers.genes.25kbTSS25kb[i])$hgnc_symbol,
             strand = "*",
             cex.main = 0.5,
             sizes = c(0.01, 0.04, 0.02, 0.02, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05),
             scale = 0.5)
}
dev.off()

# 
# Input_TGFb_rep1 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep1/NA_control_lambda.bw"), gr.which)
# H2AZ_TGFb_rep1 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep1/NA_treat_pileup.bw"), gr.which)
# H2AZ_TGFb_rep1_FE <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep1/H2AZ_TGFb_rep1_FE.bw"), gr.which)
# H2AZ_TGFb_rep1_logLR <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep1/H2AZ_TGFb_rep1_logLR.bw"), gr.which)
# 
# Input_TGFb_rep2 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep2/NA_control_lambda.bw"), gr.which)
# H2AZ_TGFb_rep2 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep2/NA_treat_pileup.bw"), gr.which)
# H2AZ_TGFb_rep2_FE <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep2/H2AZ_TGFb_rep2_FE.bw"), gr.which)
# H2AZ_TGFb_rep2_logLR <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/TGFb_rep2/H2AZ_TGFb_rep2_logLR.bw"), gr.which)
# 
# Input_WT_rep1 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep1/NA_control_lambda.bw"), gr.which)
# H2AZ_WT_rep1 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep1/NA_treat_pileup.bw"), gr.which)
# H2AZ_WT_rep1_FE <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep1/H2AZ_WT_rep1_FE.bw"), gr.which)
# H2AZ_WT_rep1_logLR <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep1/H2AZ_WT_rep1_logLR.bw"), gr.which)
# 
# Input_WT_rep2 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep2/NA_control_lambda.bw"), gr.which)
# H2AZ_WT_rep2 <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep2/NA_treat_pileup.bw"), gr.which)
# H2AZ_WT_rep2_FE <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep2/H2AZ_WT_rep2_FE.bw"), gr.which)
# H2AZ_WT_rep2_logLR <- subsetByOverlaps(import("~/Data/Tremethick/EMT/GenomeWide/macs2_analysis/WT_rep2/H2AZ_WT_rep2_logLR.bw"), gr.which)
plotTracks(list(biomTrack,
                aT.ap1Sites,
                aT.nfkbSites,
                dT.cov.input.emt_markers.wt
#                 dT.cov.h2az.emt_markers.wt,
#                 dT.cov.input.wt,
#                 dT.cov.h2az.wt,
#                 dT.cov.input.emt_markers.tgfb, 
#                 dT.cov.h2az.emt_markers.tgfb,
#                 dT.cov.input.tgfb,
#                 dT.cov.h2az.tgfb,
#                 dT.dmr
),
chromosome = as(seqnames(gr.mesenchymalMarkers.genes.25kbTSS25kb), "character")[i],
from = as.integer(start(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
to = as.integer(end(gr.mesenchymalMarkers.genes.25kbTSS25kb[i]), "integer"),
extend.right = 1000,
extend.left = 1000)


# plot histograms of insert sizes
par(mfrow = c(2,2))
lapply(seq_along(reads.input), function(x, i, n) {
  w1 <- width(as(x[[i]], "GRanges"))
  med <- median(w1)
  title <- n[[i]]
  hist(w1, 
       main = title,
       xlab = "Insert size [bp]\n (min, median, max)",
       axes = F)
  abline(v = med, col = "red")
  axis(1, at = c(median(w1)), col = c("red"), lty = NULL, col.lab = "red")
  axis(1, at = c(min(w1), max(w1)), col = c("black"))
  axis(2)
}, x = reads.input, n = names(reads.input))


par(mfrow = c(2,2))
lapply(seq_along(reads.h2az), function(x, i, n) {
  w1 <- width(as(x[[i]], "GRanges"))
  med <- median(w1)
  title <- n[[i]]
  hist(w1, 
       main = title,
       xlab = "Insert size [bp]\n (min, median, max)",
       axes = F)
  abline(v = med, col = "red")
  axis(1, at = c(median(w1)), col = c("red"), lty = NULL, col.lab = "red")
  axis(1, at = c(min(w1), max(w1)), col = c("black"))
  axis(2)
}, x = reads.h2az, n = names(reads.h2az))

#-------------loading MACS2 peak calling results----------------------
macs2_path <- "/Volumes/gduserv/Data/Tremethick/EMT/GenomeWide/macs2_analysis"

Peaks_WT <- read.table(paste(macs2_path, "NA_peaks.narrowPeak", sep = "/WT/"), header = F, as.is = T, sep = "\t")
Peaks_WT <- GRanges(Peaks_WT$V1, IRanges(Peaks_WT$V2, Peaks_WT$V3), names = Peaks_WT$V4, mcols = Peaks_WT[, c("V5", "V6", "V7", "V8", "V9", "V10")])
Peaks_WT <- subsetByOverlaps(Peaks_WT, gr.which)
seqlevels(Peaks_WT, force = T) <- seqlevels(gr.which)

Summits_WT <- import(paste(macs2_path, "NA_summits.bed", sep = "/WT/"))
Summits_WT <- subsetByOverlaps(Summits_WT, gr.which)
seqlevels(Summits_WT, force = T) <- seqlevels(gr.which)

Peaks_TGFb <- read.table(paste(macs2_path, "NA_peaks.narrowPeak", sep = "/TGFb/"), header = F, as.is = T, sep = "\t")
Peaks_TGFb <- GRanges(Peaks_TGFb$V1, IRanges(Peaks_TGFb$V2, Peaks_TGFb$V3), names = Peaks_TGFb$V4, mcols = Peaks_TGFb[, c("V5", "V6", "V7", "V8", "V9", "V10")])
Peaks_TGFb <- subsetByOverlaps(Peaks_TGFb, gr.which)
seqlevels(Peaks_TGFb, force = T) <- seqlevels(gr.which)

Summits_TGFb <- import(paste(macs2_path, "NA_summits.bed", sep = "/TGFb/"))
Summits_TGFb <- subsetByOverlaps(Summits_TGFb, gr.which)
seqlevels(Summits_TGFb, force = T) <- seqlevels(gr.which)

DiffPeaks_WT <- import(paste(macs2_path, "diff_WT_vs_TGFb_c3.0_cond1.bed", sep = "/"))
DiffPeaks_WT <- subsetByOverlaps(DiffPeaks_WT, gr.which)
seqlevels(DiffPeaks_WT, force = T) <- seqlevels(gr.which)

DiffPeaks_TGFb <- import(paste(macs2_path, "diff_WT_vs_TGFb_c3.0_cond2.bed", sep = "/"))
DiffPeaks_TGFb <- subsetByOverlaps(DiffPeaks_TGFb, gr.which)
seqlevels(DiffPeaks_TGFb, force = T) <- seqlevels(gr.which)

DiffPeaks_Common <- import(paste(macs2_path, "diff_WT_vs_TGFb_c3.0_common.bed", sep = "/"))
DiffPeaks_Common <- subsetByOverlaps(DiffPeaks_Common, gr.which)
seqlevels(DiffPeaks_Common, force = T) <- seqlevels(gr.which)

atPeaks_WT <- AnnotationTrack(Peaks_WT, name = "H2AZ Peaks, WT", col = "black")
atSummits_WT <-  AnnotationTrack(Summits_WT, name = "H2AZ Summits, WT", col = "blue")
atPeaks_TGFb <-  AnnotationTrack(Peaks_TGFb, name = "H2AZ Peaks, TGFb", col = "black")
atSummits_TGFb <- AnnotationTrack(Summits_TGFb, name = "H2AZ Summits, TGFb", col = "black")

displayPars(atPeaks_WT) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(atSummits_WT) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(atPeaks_TGFb) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
displayPars(atSummits_TGFb) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

atDiffPeaks_WT <- AnnotationTrack(DiffPeaks_WT, name = "Differential Peaks, WT only", col = "black")
atDiffPeaks_TGFb <- AnnotationTrack(DiffPeaks_TGFb, name = "Differential Peaks, TGFb only", col = "black")
atDiffPeaks_Common <- AnnotationTrack(DiffPeaks_Common, name = "Common Peaks", col = "black")
