# spp_post_processing.R

library("GenomicRanges")

# select if running script on cluster or locally
gdu <- F

# load Data

if (gdu == T) {
  path_root <- "/home/skurscheid"
} else {
  path_root <- "/Volumes/gduserv"
}

path.spp <- paste(path_root, "/Data/Tremethick/EMT/GenomeWide/spp_analysis/", sep = "")

# load IDR overlapped peaks
# WT first
WT_peaks <- read.table(paste(path.spp, "WT_rep1_vs_rep2-overlapped-peaks.txt", sep = ""), header = T, as.is = T, sep = " ")
WT_peaks$start <- as(apply(WT_peaks, 1, function(x) min(x["start1"], x["start2"])), "integer")
WT_peaks$stop <- as(apply(WT_peaks, 1, function(x) max(x["stop1"], x["stop2"])), "integer")
# filter for conservative peaks - IDR <= 0.01
WT_peaks <- WT_peaks[WT_peaks$IDR <= 0.01, ]
gr.WT_peaks <- GRanges(WT_peaks$chr1, IRanges(WT_peaks$start, WT_peaks$stop), strand = "*", WT_peaks[, c("sig.value1", "sig.value2", "idr.local", "IDR")])
names(gr.WT_peaks) <- paste("WT_peak", names(gr.WT_peaks), sep = "_")


# TGFb treated
TGFb_peaks <- read.table(paste(path.spp, "TGFB_rep1_vs_rep2-overlapped-peaks.txt", sep = ""), header = T, as.is = T, sep = " ")
TGFb_peaks$start <- as(apply(TGFb_peaks, 1, function(x) min(x["start1"], x["start2"])), "integer")
TGFb_peaks$stop <- as(apply(TGFb_peaks, 1, function(x) max(x["stop1"], x["stop2"])), "integer")
TGFb_peaks <- TGFb_peaks[TGFb_peaks$IDR <= 0.01, ]
gr.TGFb_peaks <- GRanges(TGFb_peaks$chr1, IRanges(TGFb_peaks$start, TGFb_peaks$stop), strand = "*", TGFb_peaks[, c("sig.value1", "sig.value2", "idr.local", "IDR")])
names(gr.TGFb_peaks) <- paste("TGFb_peak", names(gr.TGFb_peaks), sep = "_")


# now we have the set of high confidence peaks for each condition
# 
# create GRanges to obtain reads from all samples
gr.which.allPeaks <- c(gr.TGFb_peaks, gr.WT_peaks)
gr.which.allPeaks <- sort(gr.which.allPeaks)
gr.which.allPeaks <- reduce(gr.which.allPeaks, min.gapwidth = 147)

gr.which <- c(subsetByOverlaps(gr.TGFb_peaks, gr.MSigDB.TGFb_emt_gene_id.cfam.tss), subsetByOverlaps(gr.WT_peaks, gr.MSigDB.TGFb_emt_gene_id.cfam.tss))
gr.which <- sort(gr.which)
reduce(gr.which, min.gapwidth = 147)

path.input = "/Volumes/gduserv/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/"
path.h2az = "/Volumes/gduserv/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/"

suffix = ".Q10.sorted.MkDup.bam.tagAlign.gz"

files.input = c("Input_TGFb_rep1_S7", "Input_TGFb_rep2_S8", "Input_WT_rep1_S5", "Input_WT_rep2_S6")
files.h2az = c("H2AZ_TGFb_rep1_S3", "H2AZ_TGFb_rep2_S4", "H2AZ_WT_rep1_S1", "H2AZ_WT_rep2_S2")

# parameters for reading in tagAlign files
counts.input <- lapply(files.input, function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  print(fn)
  as(system(paste("gzip -dc ", path.input, x, suffix, "|wc -l", sep = ""), intern = T), "integer")
})
names(counts.input) <- files.input

counts.h2az <- lapply(files.h2az, function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  print(fn)
  as(system(paste("gzip -dc ", path.h2az, x, suffix, "|wc -l", sep = ""), intern = T), "integer")
})
names(counts.h2az) <- files.h2az

#--------------calculate scaling factors per sample (RPM)--------------
mean.count <- mean(c(unlist(counts.h2az), unlist(counts.input)))

#--------------total count adjustment: http://bib.oxfordjournals.org/content/14/6/671.long--------------
scale.input <- lapply(counts.input, function(x) {mean.count / x})
scale.h2az <- lapply(counts.h2az, function(x) {mean.count / x})

#--------------Input - only TGFb-induced EMT markers--------------
reads.input <- lapply(files.input, function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  print(fn)
  import.bed(gzfile(fn), which = gr.which)
})
names(reads.input) <- files.input
reads.input <- GRangesList(reads.input)
width(ranges(reads.input)) <- 150

reads.h2az <- lapply(files.h2az, function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  print(fn)
  import.bed(gzfile(fn), which = gr.which)
})
names(reads.h2az) <- files.h2az
reads.h2az <- GRangesList(reads.h2az)
width(ranges(reads.h2az)) <- 150

#--------------summarizeOverlaps for further analysis of TGFb-induced EMT markers--------------
df.h2az <- data.frame(matrix(unlist(lapply(reads.h2az, function(x) data.frame(union = assay(summarizeOverlaps(gr.which, x))))), nrow = 219, byrow = T))
colnames(df.h2az) <- names(reads.h2az)
for (i in colnames(df.h2az)){
  df.h2az[,i] <- df.h2az[,i] * scale.h2az[[i]]
}
colnames(df.h2az) <- unlist(lapply(strsplit(colnames(df.h2az), "_"), function(x) paste(x[2:3], collapse = "_")))

df.input <- data.frame(matrix(unlist(lapply(reads.input, function(x) data.frame(union = assay(summarizeOverlaps(gr.which, x))))), nrow = 219, byrow = T))
colnames(df.input) <- names(reads.input)
for (i in colnames(df.input)){
  df.input[,i] <- df.input[,i] * scale.input[[i]]
}

#--------------Input - all peaks--------------
reads.input.allPeaks <- lapply(files.input, function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  print(fn)
  import.bed(gzfile(fn), which = gr.which.allPeaks)
})
names(reads.input.allPeaks) <- files.input
reads.input.allPeaks <- GRangesList(reads.input.allPeaks)
width(ranges(reads.input.allPeaks)) <- 150

reads.h2az.allPeaks <- lapply(files.h2az, function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  print(fn)
  import.bed(gzfile(fn), which = gr.which.allPeaks)
})
names(reads.h2az.allPeaks) <- files.h2az
reads.h2az.allPeaks <- GRangesList(reads.h2az.allPeaks)
width(ranges(reads.h2az.allPeaks)) <- 150


#--------------summarizeOverlaps for further analysis of all peaks--------------
df.h2az.allPeaks <- data.frame(matrix(unlist(lapply(reads.h2az.allPeaks, function(x) data.frame(union = assay(summarizeOverlaps(gr.which.allPeaks, x))))), nrow = length(gr.which.allPeaks), byrow = T))
colnames(df.h2az.allPeaks) <- names(reads.h2az.allPeaks)
colnames(df.h2az.allPeaks) <- unlist(lapply(strsplit(colnames(df.h2az.allPeaks), "_"), function(x) paste(x[2:3], collapse = "_")))

df.input.allPeaks <- data.frame(matrix(unlist(lapply(reads.input.allPeaks, function(x) data.frame(union = assay(summarizeOverlaps(gr.which.allPeaks, x))))), nrow = length(gr.which.allPeaks), byrow = T))
colnames(df.input.allPeaks) <- names(reads.input)

#--------------un-supervised analysis--------------
source("~/Development/GeneralPurpose/R/heatmap.3.R")


#--------------calculate coverage for plotting...--------------
cov.input.emt_markers.wt.rep1 <- coverage(reads.input[["Input_WT_rep1_S5"]]) * scale.input[["Input_WT_rep1_S5"]]
cov.input.emt_markers.wt.rep2 <- coverage(reads.input[["Input_WT_rep2_S6"]]) * scale.input[["Input_WT_rep2_S6"]]
cov.h2az.emt_markers.wt.rep1 <- coverage(reads.h2az[["H2AZ_WT_rep1_S1"]]) * scale.h2az[["H2AZ_WT_rep1_S1"]]
cov.h2az.emt_markers.wt.rep2 <- coverage(reads.h2az[["H2AZ_WT_rep2_S2"]]) * scale.h2az[["H2AZ_WT_rep2_S2"]]

cov.input.emt_markers.tgfb.rep1 <- coverage(reads.input[["Input_TGFb_rep1_S7"]]) * scale.input[["Input_TGFb_rep1_S7"]]
cov.input.emt_markers.tgfb.rep2 <- coverage(reads.input[["Input_TGFb_rep2_S8"]]) * scale.input[["Input_TGFb_rep2_S8"]] 
cov.h2az.emt_markers.tgfb.rep1 <- coverage(reads.h2az[["H2AZ_TGFb_rep1_S3"]]) * scale.h2az[["H2AZ_TGFb_rep1_S3"]] 
cov.h2az.emt_markers.tgfb.rep2 <- coverage(reads.h2az[["H2AZ_TGFb_rep2_S4"]]) * scale.h2az[["H2AZ_TGFb_rep2_S4"]]


