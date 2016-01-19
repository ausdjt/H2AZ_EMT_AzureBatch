library(rtracklayer)
cf3.chrom <- read.table("~/Data/RefGenomes/CanFam3.1/CanFam3.chromSizes.txt", header = F, as.is = T, sep = "\t")
sum(as.numeric(cf3.chrom[,2]))
rm1 <- import("canFam3_repeat_regions.bed")
sum(width(reduce(rm1)))

sum(as.numeric(cf3.chrom[,2])) - sum(width(reduce(rm1)))
# [1] 1382640842
