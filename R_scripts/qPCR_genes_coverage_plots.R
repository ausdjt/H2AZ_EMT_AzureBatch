# load libraries
library("gdata")
library("Gviz")
library("rtracklayer")
library("biomaRt")
library("GenomicRanges")
library("GenomicAlignments")
library(BSgenome.Cfamiliaris.UCSC.canFam3)

options(ucscChromosomeNames=FALSE)

setwd('~/Data/Tremethick/EMT/GenomeWide/')

source("~/Development/GeneralPurpose/R/binnedAverage.R")
source("~/Development/GeneralPurpose/R/binnedSum.R")
source("~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/R_scripts/calculateCoverage.R")


# set upd biomaRt connection ----------------------------------------------
ensembl <- useEnsembl(biomart = "ensembl", host = "asia.ensembl.org")
dog <- useEnsembl(biomart = "ensembl", dataset = "cfamiliaris_gene_ensembl", host = "asia.ensembl.org")

# load qPCR analysis results ----------------------------------------------
load("~/Data/Tremethick/EMT/MDCK qPCR data/qDE.limma.tab.rda")

# plotting H2A.Z coverage across qPCR genes which do not change
i1 <- intersect(rownames(qPCRGenesTab), rownames(m2))
toPlot <- qPCRGenesTab[i1,]$ensembl_gene_id
names(toPlot) <- rownames(qPCRGenesTab[i1,])

gr.which <- promoters(sort(gr.qPCRGenesPositions[toPlot]), upstream = 40000, downstream = 40000)

path.input = "~/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_removed/"
path.h2az = "~/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_removed/"

suffix = ".Q10.sorted.DeDup.bam"

files.input = c("Input_TGFb_rep1_S7", "Input_TGFb_rep2_S8", "Input_WT_rep1_S5", "Input_WT_rep2_S6")
files.h2az = c("H2AZ_TGFb_rep1_S3", "H2AZ_TGFb_rep2_S4", "H2AZ_WT_rep1_S1", "H2AZ_WT_rep2_S2")

# parameters for reading in BAM files
# flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = F, isSecondaryAlignment = F)#, isFirstMateRead = T, isSecondMateRead = F)
flag <- scanBamFlag(isProperPair = T, isDuplicate = F)
SBParam.all <- ScanBamParam(flag = flag, simpleCigar = T, what = c("rname", "strand", "pos", "qwidth")) #
SBParam <- ScanBamParam(flag = flag, simpleCigar = T, what = c("rname", "strand", "pos", "qwidth"), which = gr.which)

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

# calculate coverage from replicates
cov.input.wt <- (coverage.input[["Input_WT_rep1_S5"]] + coverage.input[["Input_WT_rep2_S6"]]) /2
cov.input.tgfb <- (coverage.input[["Input_TGFb_rep1_S7"]] + coverage.input[["Input_TGFb_rep2_S8"]]) /2
cov.h2az.wt <- (coverage.h2az[["H2AZ_WT_rep1_S1"]] + coverage.h2az[["H2AZ_WT_rep2_S2"]]) / 2
cov.h2az.tgfb <- (coverage.h2az[["H2AZ_TGFb_rep1_S3"]] + coverage.h2az[["H2AZ_TGFb_rep2_S4"]]) / 2

#----------preparing data for plotting of sequencing coverage--------------
# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.input.wt <- calculateCoverage(step = 1, gr.which, cov.input.wt, func = "mean")
bA.cov.input.tgfb <- calculateCoverage(step = 1, gr.which, cov.input.tgfb, func = "mean")
bA.cov.h2az.wt <- calculateCoverage(step = 1, gr.which, cov.h2az.wt, func = "mean")
bA.cov.h2az.tgfb <- calculateCoverage(step = 1, gr.which, cov.h2az.tgfb, func = "mean")

# per samples coverage
bA.cov.h2az.tgfb1 <- calculateCoverage(step = 1, gr.which, coverage.h2az[["H2AZ_WT_rep1_S1"]], func = "mean")
bA.cov.h2az.tgfb2 <- calculateCoverage(step = 1, gr.which, coverage.h2az[["H2AZ_WT_rep2_S2"]], func = "mean")
bA.cov.h2az.wt1 <- calculateCoverage(step = 1, gr.which, coverage.h2az[["H2AZ_TGFb_rep1_S3"]], func = "mean")
bA.cov.h2az.wt2 <- calculateCoverage(step = 1, gr.which, coverage.h2az[["H2AZ_TGFb_rep2_S4"]], func = "mean")

#----------creating DataTrack objects for visualization using Gviz------------------------------
# ['#d7191c','#fdae61','#abd9e9','#2c7bb6']
dT.cov.input.wt <- DataTrack(bA.cov.input.wt, type = "h", col = "#abd9e9", name = "Input Epithelial\n[rpm]", strand = "*", cex = 2)
dT.cov.input.tgfb <- DataTrack(bA.cov.input.tgfb, type = "h", col = "#2c7bb6", name = "Input Mesenchymal\n[rpm]", strand = "*", cex = 2)
dT.cov.h2az.wt <- DataTrack(bA.cov.h2az.wt, type = "h", col = "#fdae61", name = "H2AZ Epithelial\n[rpm]", strand = "*", cex = 2)
dT.cov.h2az.tgfb <- DataTrack(bA.cov.h2az.tgfb, type = "h", col = "#d7191c", name = "H2AZ Mesenchymal\n[rpm]", strand = "*", cex = 2)

dpList <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white", "cex.title" = 0.5, rotation.title = 270, cex.axis = 0.6, lwd.title = 1)

displayPars(dT.cov.input.wt) <- dpList
displayPars(dT.cov.input.tgfb) <- dpList
displayPars(dT.cov.h2az.wt) <- dpList
displayPars(dT.cov.h2az.tgfb) <- dpList


# plotting ----------------------------------------------------------------
gr.plot <- promoters(sort(gr.qPCRGenesPositions[toPlot]), upstream = 2000, downstream = 2000)
axisTrack <- GenomeAxisTrack()

suffix = "TSS_2kb"
mainDir = "~/Data/Tremethick/EMT/GenomeWide/Publication_plots/"
subDir = "qPCR Genes - no change"

if (file.exists(paste(mainDir, subDir, sep = ""))){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
}

for (i in 1:length(gr.plot)){
  # Data from whole genome ChIP-Seq
  
  chromosome(dT.cov.input.wt) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.input.tgfb) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.h2az.wt) <- seqnames(gr.plot)[i]
  chromosome(dT.cov.h2az.tgfb) <- seqnames(gr.plot)[i]
  
  max.y.tss <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
  displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y.tss))
  displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y.tss))
  
  biomTrack.ensembl <- BiomartGeneRegionTrack(genome = "canFam3", 
                                              gene = gr.plot[i]$ensembl_gene_id,
                                              mart = dog,
                                              stacking = "squish",
                                              "fontcolor.title" = "black", 
                                              "background.title" = "white", 
                                              "col.axis" = "black", 
                                              "col.frame" = "white", 
                                              "fill" = "black",
                                              "fontcolor.group" = "black",
                                              "protein_coding" = "black",
                                              "utr3" = "black",
                                              "utr5" = "black",
                                              "col.line" = "darkgrey",
                                              name = NULL)
  #displayPars(biomTrack.ensembl) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white", "fill" = "black")
  
  pdf(paste(gr.plot$hgnc_symbol[i], "_", suffix,".pdf", sep = ""))
  plotTracks(list(axisTrack,
                  biomTrack.ensembl,
                  dT.cov.input.wt, 
                  dT.cov.h2az.wt,
                  dT.cov.input.tgfb, 
                  dT.cov.h2az.tgfb
  ),
  chromosome = as(seqnames(gr.plot), "character")[i],
  from = as.integer(start(gr.plot[i]), "integer"),
  to = as.integer(end(gr.plot[i]), "integer"),
  extend.right = 1000,
  extend.left = 1000,
  main = paste(gr.plot$hgnc_symbol[i], " (", width(gr.plot[i]), "bp)", sep = ""),
  strand = "*",
  cex.main = 0.5,
  #sizes = c(0.02, 0.06, 0.04, 0.04, 0.04, 0.2, 0.2, 0.2, 0.2),
  sizes = c(0.02, 0.04, 0.235, 0.235, 0.235, 0.235),
  scale = 0.5)
  dev.off()
}
setwd('~/Data/Tremethick/EMT/GenomeWide/')


