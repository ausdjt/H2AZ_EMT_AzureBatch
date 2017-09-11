# Figure 4 Panels A & B ------------------------------------------------------------
emtUp <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_up_genes.bed")
emtDown <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_down_genes.bed")
emtGenes <- c(emtUp, emtDown)
names(emtGenes) <- emtGenes$name
cfamEnsGenesSigEMTCells[emtUp$name]

# coverage plot for TGFB1
m1 <- melt(dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix[grep("ENSCAFG00000005014", dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix$X4),])
m2 <- melt(dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix[grep("ENSCAFG00000005014", dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTup_normal.matrix.gz"]]$computeMatrix$X4),])
m3 <- melt(dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix[grep("ENSCAFG00000002653", dataList[["H2AZ-TGFb_vs_Input-TGFb_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix$X4),])
m3$set <- "TGFb"
m4 <- melt(dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix[grep("ENSCAFG00000002653", dataList[["H2AZ-WT_vs_Input-WT_normal.readCount.subtract_RPKM_TanEMTdown_normal.matrix.gz"]]$computeMatrix$X4),])
m4$set <- "WT"
m3 <- rbind(m3,m4)
m3$bin <- rep(1:300, 2)
m3$set <- as.factor(m3$set)
m3$set <- relevel(m3$set, ref = "WT")
m3$gene <- "EPCAM"

m1$set <- "TGFb"
m1$bin <- 1:300
m2$set <- "WT"
m2$bin <- 1:300
m1 <- rbind(m1,m2)
m1$gene <- "TGFB1"
m1$set <- as.factor(m1$set)
m1$set <- relevel(m1$set, ref = "WT")

plotTGFb <- ggplot(m1, aes(bin, value, colour = set)) + 
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, 10), 
              se = F, size = lineSize) + 
  facet_wrap(c("set"), nrow = 2) + 
  ylab(NULL) +
  xlab(NULL) +
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  labs(color = "Sample") +
  theme(axis.line = element_line(colour = "black", size = axisLineSize),
        axis.text = element_text(size = 8),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values = c(wtCol, TGFbCol))
ggsave("Figure4_PanelA_TGFB1_coverage.pdf", plotTGFb, height = 64, width = 98, units = "mm")

# coverage plot for EPCAM ENSCAFG00000002653 -------------------------------------------------
plotEPCAM <- ggplot(m3, aes(bin, value, colour = set)) + 
  geom_smooth(method = "lm",
              formula = y ~ splines::ns(x, 10), 
              se = F, size = lineSize) + 
  facet_wrap(c("set"), nrow = 2) + 
  ylab("Mean coverage (input subtracted)\n[RPKM]") +
  xlab(NULL) +
  labs(color = "Sample") +
  geom_vline(xintercept = 150, colour = vlineCol, linetype = "longdash") +
  scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300),
                     labels=c("-1500", "-1000", "-500", "TSS", "+500", "+1000", "+1500")) +
  theme(axis.line = element_line(colour = "black", size = axisLineSize),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual(values = c(wtCol, TGFbCol))
legendPlot <- g_legend(plotEPCAM)
plotEPCAM <- plotEPCAM + theme(legend.position = "none")
ggsave("Figure4_PanelB_EPCAM_coverage.pdf", plotEPCAM, height = 64, width = 98, units = "mm")
ggsave("Figure4_PanelAB_legend.pdf", legendPlot, height = 98, width = 98, units = "mm")
