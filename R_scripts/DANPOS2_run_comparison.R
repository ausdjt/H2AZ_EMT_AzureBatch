danpos2.rerun.results <- read.table("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT_rerun/result/TGFb_H2AZ-WT_H2AZ.positions.integrative.xls",
                              header = T, 
                              as.is = T,
                              sep = "\t")

gr.danpos2.rerun.results <- GRanges(danpos2.rerun.results$chr, IRanges(danpos2.rerun.results$start, danpos2.rerun.results$end), strand = "*", danpos2.rerun.results[, c(4:23)])
ns <- seqlevels(gr.danpos2.rerun.results)[-which(seqlevels(gr.danpos2.rerun.results) %in% c("X", "MT"))][order(as.integer(seqlevels(gr.danpos2.rerun.results)[-which(seqlevels(gr.danpos2.rerun.results) %in% c("X", "MT"))]))]
seqlevels(gr.danpos2.rerun.results) <- c(ns, "X", "MT")
gr.danpos2.rerun.results <- sort(gr.danpos2.rerun.results)


# compare results of DANPOS2 run 1 (MkDup) & 2 (DeDup + pheight =  --------
# here for TGFb
h1 <- hist(subsetByOverlaps(gr.danpos2.rerun.results, promoters(gr.mesenchymalMarkers.genes[2], upstream = 1500, downstream = 1500))$smt_log2FC)
h2 <- hist(subsetByOverlaps(danpos2.anno, promoters(gr.mesenchymalMarkers.genes[2], upstream = 1500, downstream = 1500))$smt_log2FC)
h3 <- hist(subsetByOverlaps(gr.danpos2.rerun.results, promoters(gr.mesenchymalMarkers.genes[2], upstream = 1500, downstream = 1500))$smt_diff_FDR)
h4 <- hist(subsetByOverlaps(danpos2.anno, promoters(gr.mesenchymalMarkers.genes[2], upstream = 1500, downstream = 1500))$smt_diff_FDR)

pdf("DANPOS2_run_comparison_TGFb.pdf")
par(mfrow = c(2,2))
plot(h1, main = "DANPOS2 re-run TGFb1\nTSS +/- 1.5Kb [log2 FC]")
plot(h2, main = "DANPOS2 original TGFb1\nTSS +/- 1.5Kb [log2 FC]")
plot(h3, main = "DANPOS2 re-run TGFb1\nTSS +/- 1.5Kb [FDR]")
plot(h4, main = "DANPOS2 original TGFb1\nTSS +/- 1.5Kb [FDR]")
dev.off()

# here across all nucleosome positions
h1 <- hist(gr.danpos2.rerun.results$smt_log2FC)
h2 <- hist(danpos2.anno$smt_log2FC)
h3 <- hist(gr.danpos2.rerun.results$smt_diff_FDR)
h4 <- hist(danpos2.anno$smt_diff_FDR)

pdf("DANPOS2_run_comparison_WG.pdf", )
par(mfrow = c(2,2))
plot(h1, main = "DANPOS2 re-run\n whole genome [log2 FC]")
plot(h2, main = "DANPOS2 original TGFb1\n whole genome [log2 FC]")
plot(h3, main = "DANPOS2 re-run TGFb1\n whole genome [FDR]")
plot(h4, main = "DANPOS2 original TGFb1\n whole genome [FDR")
dev.off()

h1 <- hist(gr.danpos2.rerun.results$control_smt_val)
h2 <- hist(danpos2.anno$control_smt_val)
h3 <- hist(gr.danpos2.rerun.results$treat_smt_val)
h4 <- hist(danpos2.anno$treat_smt_val)

summary(gr.danpos2.rerun.results$control_smt_val)
summary(danpos2.anno$control_smt_val)
summary(gr.danpos2.rerun.results$treat_smt_val)
summary(danpos2.anno$treat_smt_val)

