# nucleR analysis
# load libraries
require("GenomicAlignments")
require("nucleR")

# select if running script on cluster or locally
gdu <- T

# load Data

if (gdu == F) {
  path_root <- "/home/skurscheid"
} else {
  path_root <- "/Volumes/gduserv"
}

path.input <- paste(path_root, "/Data/Tremethick/EMT/GenomeWide/Input/processed_data/duplicates_marked/", sep = "")
path.h2az <- paste(path_root, "/Data/Tremethick/EMT/GenomeWide/H2AZ/processed_data/duplicates_marked/", sep = "")

suffix = ".Q10.sorted.MkDup.bam"

files.input = c("Input_TGFb_rep1_S7", "Input_TGFb_rep2_S8", "Input_WT_rep1_S5", "Input_WT_rep2_S6")
files.h2az = c("H2AZ_TGFb_rep1_S3", "H2AZ_TGFb_rep2_S4", "H2AZ_WT_rep1_S1", "H2AZ_WT_rep2_S2")

# parameters for reading in BAM files
flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = F, isSecondaryAlignment = F)
SBParam.all <- ScanBamParam(flag = flag, simpleCigar = T, what = c("strand", "pos", "qwidth")) #

reads.input <- lapply(files.input[1], function(x){
  fn <- paste(path.input, x, suffix, sep = "")
  print(fn)
  readGAlignmentPairs(fn, param = SBParam.all, use.names = TRUE)
})
names(reads.input) <- files.input[[1]]

rd1 <- as(reads.input[[1]], "GRanges")
rd1 <- as(rd1, "RangedData")
reads_input <- processReads(rd1, type = "paired", fragmentLen = 200)
reads_input_trim <- processReads(rd1, type = "paired", fragmentLen = 200, trim = 40)
coverage

reads.h2az <- lapply(files.h2az[1], function(x){
  fn <- paste(path.h2az, x, suffix, sep = "")
  print(fn)
  readGAlignmentPairs(fn, param = SBParam.all, use.names = F)
})

