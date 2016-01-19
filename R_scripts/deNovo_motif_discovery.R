# motif discovery
# 1. extract sequences under peaks
# 2. select common (non-diff) peaks
# 3. select condition-specific peaks (differential summit values) 

# load libraries
library(GenomicRanges)
library(rGADEM)
library(BSgenome.Cfamiliaris.UCSC.canFam3)
library(Biostrings)

# go to working directory
setwd('~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/')

# load genome
genome <- BSgenome.Cfamiliaris.UCSC.canFam3
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))

# load GRanges object of DANPOS2 results
load("gr.danpos2.results.rda")
seqlevels(gr.danpos2.results)[grep("MT", seqlevels(gr.danpos2.results))] <- "M"
seqlevels(gr.danpos2.results) <- seqlevels(genome)
gr.danpos2.results <- sort(gr.danpos2.results)

# we start by looking at promoters of all genes
load("Cfam3.genes.rda")
seqlevels(Cfam3.genes)[40] <- "M"

Cfam3.genes.promoters <- promoters(Cfam3.genes, upstream = 400, downstream = 200)
# tile the promoters into 3 200bp regions, to approximate nucleosome spacing
Cfam3.genes.promoters <- unlist(tile(Cfam3.genes.promoters, width = 200))

Cfam3.genes.tss.up2.nucleosome <- Cfam3.genes.promoters[seq(1, 72669, 3)]
Cfam3.genes.tss.up1.nucleosome <- Cfam3.genes.promoters[seq(2, 72669, 3)]
Cfam3.genes.tss.down1.nucleosome <- Cfam3.genes.promoters[seq(3, 72669, 3)]
# subset DANPOS2 nucleosome positions
dp2.subset.up2 <- subsetByOverlaps(gr.danpos2.results, Cfam3.genes.tss.up2.nucleosome)
dp2.subset.up1 <- subsetByOverlaps(gr.danpos2.results, Cfam3.genes.tss.up1.nucleosome)
dp2.subset.down1 <- subsetByOverlaps(gr.danpos2.results, Cfam3.genes.tss.down1.nucleosome)

# test run
ds1 <-  DNAStringSet(BSgenomeViews(genome, dp2.subset.up2), use.names = TRUE)
names(ds1) <- paste(as.character(seqnames(dp2.subset.up2)), start(dp2.subset.up2), end(dp2.subset.up2), sep = "_")
writeXStringSet(ds1, "ds1.fa")

ds2 <-  DNAStringSet(BSgenomeViews(genome, dp2.subset.up1), use.names = TRUE)
names(ds2) <- paste(as.character(seqnames(dp2.subset.up1)), start(dp2.subset.up1), end(dp2.subset.up1), sep = "_")
writeXStringSet(ds2, "ds2.fa")

ds3 <-  DNAStringSet(BSgenomeViews(genome, dp2.subset.down1), use.names = TRUE)
names(ds3) <- paste(as.character(seqnames(dp2.subset.down1)), start(dp2.subset.down1), end(dp2.subset.down1), sep = "_")
writeXStringSet(ds3, "ds3.fa")

gadem <- GADEM(ds1, verbose = 1)
