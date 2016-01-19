library("rtracklayer")
# MACS2 post-processing
# results of bdgdiff need to be annotated


setwd("~/Data/Tremethick/EMT/")
common_peaks <- import("diff_WT_vs_TGFb_c3.0_common.bed")
WT_peaks <- import("diff_WT_vs_TGFb_c3.0_cond1.bed")
WT_summits <- import("WT/NA_summits.bed")
WT_H2AZ_coverage <- import.bedGraph("WT/NA_treat_pileup.bdg")

TGFb_peaks <- import("diff_WT_vs_TGFb_c3.0_cond2.bed")
TGFb_summits <- import("TGFb/NA_summits.bed")
TGFb_H2AZ_coverage <- import.bedGraph("TGFb/NA_treat_pileup.bdg")

