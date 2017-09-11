# volcano plots
labelSize <- 3

p1 <- ggplot(volcanoData, aes(x = logFC, y = logqval, color = epi_mes)) +
  geom_point(size = 0.8) + 
  geom_point(data = subset(volcanoData, experiment == "TGFb" & logqval > 200), aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = subset(volcanoData, experiment == "TGFb" & logqval > 200), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  geom_point(data = subset(volcanoData, experiment == "TGFb" & epi_mes %in% c("TGFB1", "H2A.Z")), aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = subset(volcanoData, experiment == "TGFb" & epi_mes %in% c("TGFB1", "H2A.Z")), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  geom_point(data = volcanoData[grep("FN1|CDH2", volcanoData$external_gene_name)], aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = volcanoData[grep("FN1|CDH2", volcanoData$external_gene_name)], aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & logqval > 18), aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & logqval > 18), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & epi_mes %in% c("TGFB1")), aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & epi_mes %in% c("TGFB1")), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & target_id == "ENSCAFG00000002653"), aes(x = logFC, y = logqval, color = epi_mes), size = 1.5) +
  geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & target_id == "ENSCAFG00000002653"), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
  facet_wrap("experiment", scales = "free_y", nrow = 2, ncol = 1) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  coord_cartesian(xlim = c(-8,8)) + 
  theme(panel.background = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8), 
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(c(0.5,0.75), units = "cm" )) +
  xlab("log2 fold-change") +
  ylab("-log10(q-value)") + 
  labs(color = "Gene/Marker type") + scale_colour_hue(labels = c("Epithelial", "H2A.Z", "Mesenchymal", "TGFB1"))
p1 <- p1 + theme(legend.position = "none")
p1
ggsave("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/Figure2_DE.pdf", p1, height = 216, width = 196, units = "mm", useDingbats = F)
