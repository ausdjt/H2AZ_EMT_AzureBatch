

hmData <- tidyr::spread(volcanoData[volcanoData$qval < 0.1, c("target_id", "external_gene_name", "expression", "experiment", "epi_mes")], experiment, expression)
hmData <- hmData[order(hmData$TGFb),]
hmData$target_id <- as.factor(hmData$target_id)
hmData$target_id <- ordered(hmData$target_id, levels=levels(hmData$target_id)[unclass(hmData$target_id)])
hmDataLong <- melt(hmData, id.vars = c("target_id", "epi_mes"), measure.vars = c("TGFb", "H2A.Z KD"))
hmDataLong <- hmDataLong[!hmDataLong$epi_mes %in% c("H2A.Z", "TGFB1"),]
hmDataLong[which(is.na(hmDataLong$value)),]$value <- "N.S."
hmDataLong$value <- as.factor(hmDataLong$value)
hmDataLong$epi_mes <- as.factor(hmDataLong$epi_mes)
levels(hmDataLong$epi_mes) <- c("Epithelial", "Mesenchymal")

Fig2A <- ggplot(hmDataLong, aes(target_id, variable, color = value)) + 
  geom_raster(aes(fill=value), interpolate = F, show.legend = F) +
  ylab("Experiment") +
  xlab("Genes") +
  facet_grid(. ~epi_mes, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "lightgrey"), 
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = "white")) +
  scale_fill_manual(values=c("green", "grey", "red"),
                    breaks=c("1","2","3"),
                    labels=c("down","N.S.","up")) +
  guides(fill = guide_legend(title = NULL))
Fig2A
ggsave("Figure_2A_Heatmap.png", Fig2A, height = 55, width = 128, units = "mm", device = "png")
