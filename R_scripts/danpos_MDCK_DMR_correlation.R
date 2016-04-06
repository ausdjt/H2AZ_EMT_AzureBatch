
upstream <- 100
downstream <- 100

gr.danpos2.results.treat <- gr.danpos2.results[which(gr.danpos2.results$smt_log2FC > 4)]
gr.danpos2.results.treat.up500 <- promoters(gr.danpos2.results.treat[which(gr.danpos2.results.treat$smt_diff_FDR < 0.01)], 
                                            upstream = upstream, 
                                            downstream = 0)
gr.danpos2.results.treat.down500 <- promoters(gr.danpos2.results.treat[which(gr.danpos2.results.treat$smt_diff_FDR < 0.01)], 
                                              upstream = 0, 
                                              downstream = downstream)

gr.danpos2.results.ctrl <- gr.danpos2.results[which(gr.danpos2.results$smt_log2FC < -1)]
gr.danpos2.results.ctrl.up500 <- promoters(gr.danpos2.results.ctrl[which(gr.danpos2.results.ctrl$smt_diff_FDR < 0.01)], 
                                           upstream = upstream, 
                                           downstream = 0)
gr.danpos2.results.ctrl.down500 <- promoters(gr.danpos2.results.ctrl[which(gr.danpos2.results.ctrl$smt_diff_FDR < 0.01)], 
                                             upstream = 0, 
                                             downstream = downstream)

avgDMR0 <- gr.mdck_dmr$average_DNA_methylation_difference
avgDMR0 <- avgDMR0[!is.na(avgDMR0)]
dens0 <- density(avgDMR0)

avgDMR1 <- subsetByOverlaps(gr.mdck_dmr, gr.danpos2.results.treat[which(gr.danpos2.results.treat$smt_diff_FDR < 0.01)])$average_DNA_methylation_difference
avgDMR1 <- avgDMR1[!is.na(avgDMR1)]
dens1 <- density(avgDMR1)

avgDMR2 <- subsetByOverlaps(gr.mdck_dmr, c(gr.danpos2.results.treat.up500, gr.danpos2.results.treat.up500))$average_DNA_methylation_difference
avgDMR2 <- avgDMR2[!is.na(avgDMR2)]
dens2 <- density(avgDMR2)

avgDMR3 <- subsetByOverlaps(gr.mdck_dmr, c(gr.danpos2.results.treat.down500, gr.danpos2.results.treat.up500))$average_DNA_methylation_difference
avgDMR3 <- avgDMR3[!is.na(avgDMR3)]
dens3 <- density(avgDMR3)


avgDMR4 <- subsetByOverlaps(gr.mdck_dmr, gr.danpos2.results.ctrl[which(gr.danpos2.results.ctrl$smt_diff_FDR < 0.01)])$average_DNA_methylation_difference
avgDMR4 <- avgDMR4[!is.na(avgDMR4)]
dens4 <- density(avgDMR4)

avgDMR5 <- subsetByOverlaps(gr.mdck_dmr, gr.danpos2.results.ctrl.up500)$average_DNA_methylation_difference
avgDMR5 <- avgDMR5[!is.na(avgDMR5)]
dens5 <- density(avgDMR5)

avgDMR6 <- subsetByOverlaps(gr.mdck_dmr, gr.danpos2.results.ctrl.down500)$average_DNA_methylation_difference
avgDMR6 <- avgDMR6[!is.na(avgDMR6)]
dens6 <- density(avgDMR6)

pdf("MDCK_DMR_Density.pdf", height = 7, width = 8)
plot(dens0, 
     ylim = c(0, 2.25), 
     xlim = c(-1, 1), 
     col = "lightgray", 
     lwd = 5, 
     bty = "none", 
     main = "",
     cex = 2,
     axes = F,
     xlab = "",
     ylab = "")
lines(dens1$x, dens1$y, col = "blue", lwd = 5)
lines(dens4$x, dens4$y, col = "brown", lwd = 5)
axis(1, lwd = 3)
axis(2, lwd = 3)
mtext("% Difference in methylation", 1 , cex = 2, padj = 2.2)
mtext("Density", 2, cex = 2, padj = -1.8)
legend("topleft", 
       c("All DMRs", "H2A.Z gain DMRs", "H2A.Z loss DMRs"), 
       col = c("lightgray", "blue", "brown"),
       border = "white",
       bty = "n",
       lwd = 3, 
       cex = 1.4)
dev.off()

ks.test(dens0$y, dens4$y)
