require(ggplot2)
require(ggrepel)
require(data.table)
require(VennDiagram)

setwd("~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/")

load("~/Data/Tremethick/EMT/RNA-Seq/NB501086_0082_RDomaschenz_JCSMR_mRNAseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/resultsCompressed.rda")
load("~/Development/JCSMR-Tremethick-Lab/EMT_shiny_app/MDCK_EMT_RNA-Seq/data/cfamEnsGenesSigEMTCells.rda")

emtUp <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_up_genes.bed")
emtDown <- rtracklayer::import("~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Tan_et_al_EMT_down_genes.bed")
emtGenes <- c(emtUp, emtDown)
names(emtGenes) <- emtGenes$name

deTGFbTab <- as.data.table(resultsCompressed[["MDCK"]]$sleuth_results.gene[["conditionMDCKTGFb"]])
deshZTab <- as.data.table(resultsCompressed[["MDCK"]]$sleuth_results.gene[["conditionMDCKshZ"]])
setkey(deshZTab, target_id)
setkey(deTGFbTab, target_id)

cfamEnsGenesSigEMTCells <- data.table(cfamEnsGenesSigEMTCells)
setkey(cfamEnsGenesSigEMTCells, ensembl_gene_id)
cfamEnsGenesSigEMTCells["ENSCAFG00000005014"]$epi_mes <- "TGFB1"
cfamEnsGenesSigEMTCells["ENSCAFG00000010615"]$epi_mes <- "H2A.Z"
cfamEnsGenesSigEMTCells$expression <- "NA"
setkey(cfamEnsGenesSigEMTCells, ensembl_gene_id)
cfamEnsGenesSigEMTCells[emtDown$name, "expression"] <- rep("down", length(emtDown))
cfamEnsGenesSigEMTCells[emtUp$name, "expression"] <- rep("up", length(emtUp))

# volcano plots
m1 <- merge(deTGFbTab[,c("target_id", "qval", "b", "pval")], 
            cfamEnsGenesSigEMTCells[,c("ensembl_gene_id", "external_gene_name", "epi_mes", "expression")], 
            by.x = "target_id", 
            by.y = "ensembl_gene_id")
m1$experiment <- "TGFb"
m1$logqval <- -log10(m1$qval)
m1$logFC <- log2(exp(m1$b))
# have to manipulate data in order to deal with Inf when plotting
m1[which(m1$logqval == Inf), "logqval"] <- 310
m2 <- merge(deshZTab[,c("target_id", "qval", "b", "pval")], 
            cfamEnsGenesSigEMTCells[,c("ensembl_gene_id", "external_gene_name", "epi_mes", "expression")], 
            by.x = "target_id", 
            by.y = "ensembl_gene_id")
m2$experiment <- "H2A.Z KD"
m2$logqval <- -log10(m2$qval)
m2$logFC <- log2(exp(m2$b))
table(m1$expression)
setkey(m2, "target_id")
m2[deshZTab[m2$target_id][which(deshZTab[m2$target_id]$b > 0),]$target_id, "expression"] <- "up"
m2[deshZTab[m2$target_id][which(deshZTab[m2$target_id]$b < 0),]$target_id, "expression"] <- "down"

volcanoData <- rbind(m1, m2)
table(volcanoData$experiment, volcanoData$expression)

volcanoData$experiment <- as.factor(volcanoData$experiment)
volcanoData$experiment <- relevel(volcanoData$experiment, "TGFb")

labelSize <- 2

p1 <- ggplot(volcanoData, aes(x = logFC, y = logqval, color = epi_mes)) +
      geom_point(size = 0.05) + 
      geom_point(data = subset(volcanoData, experiment == "TGFb" & logqval > 200), aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = subset(volcanoData, experiment == "TGFb" & logqval > 200), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      geom_point(data = subset(volcanoData, experiment == "TGFb" & epi_mes %in% c("TGFB1", "H2A.Z")), aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = subset(volcanoData, experiment == "TGFb" & epi_mes %in% c("TGFB1", "H2A.Z")), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      geom_point(data = volcanoData[grep("FN1|CDH2", volcanoData$external_gene_name)], aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = volcanoData[grep("FN1|CDH2", volcanoData$external_gene_name)], aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & logqval > 18), aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & logqval > 18), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & epi_mes %in% c("TGFB1")), aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & epi_mes %in% c("TGFB1")), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      geom_point(data = subset(volcanoData, experiment == "H2A.Z KD" & target_id == "ENSCAFG00000002653"), aes(x = logFC, y = logqval, color = epi_mes), size = 1) +
      geom_text_repel(data = subset(volcanoData, experiment == "H2A.Z KD" & target_id == "ENSCAFG00000002653"), aes(x = logFC, y = logqval, label = external_gene_name), show.legend = F, size = labelSize) +
      facet_wrap("experiment", scales = "free") +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
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
l1 <- g_legend(p1)
p1 <- p1 + theme(legend.position = "none")
ggsave("Figure2_CD_Volcano_Plots_sized.pdf", p1, width = 196, height = 98, units = "mm", useDingbats = F)
ggsave("Figure2_CD_Volcano_Plots_legend_only.pdf", l1, width = 196, height = 98, units = "mm", useDingbats = F)

# Venn diagrams (argh) ----------------------------------------------------
TGFbup <- volcanoData[volcanoData$experiment == "TGFb" & volcanoData$expression == "up"]$target_id
TGFbdown <- volcanoData[volcanoData$experiment == "TGFb" & volcanoData$expression == "down"]$target_id
shZup <- volcanoData[volcanoData$experiment == "H2A.Z KD" & volcanoData$expression == "up"]$target_id
shZdown <- volcanoData[volcanoData$experiment == "H2A.Z KD" & volcanoData$expression == "down"]$target_id
upList <- list(TGFb = TGFbup, shZ = shZup)
downList <- list(TGFb = TGFbdown, shZ = shZdown)
vennUp <- gplots::venn(upList, intersection = T)
vennDown <- gplots::venn(downList, intersections = T)

grid.newpage()
vUp <- draw.pairwise.venn(area1 = length(attr(vennUp, "intersections")$TGFb) + length(attr(vennUp, "intersections")$`TGFb:shZ`),
                   area2 = length(attr(vennUp, "intersections")$shZ) + length(attr(vennUp, "intersections")$`TGFb:shZ`),
                   cross.area = length(attr(vennUp, "intersections")$`TGFb:shZ`),
                   category = c("TGFb", "H2A.Z KD"),
                   scaled = F)

grid.newpage()
vDown <- draw.pairwise.venn(area1 = length(attr(vennDown, "intersections")$TGFb) + length(attr(vennDown, "intersections")$`TGFb:shZ`),
                          area2 = length(attr(vennDown, "intersections")$shZ) + length(attr(vennDown, "intersections")$`TGFb:shZ`),
                          cross.area = length(attr(vennDown, "intersections")$`TGFb:shZ`),
                          category = c("TGFb", "H2A.Z KD"),
                          scaled = F)

pdf("Figure2_A_VennUp.pdf", height = 3.85, width = 3.85)
grid.draw(vUp)
dev.off()

pdf("Figure2_B_VennDown.pdf", height = 3.85, width = 3.85)
grid.draw(vDown)
dev.off()


# check H2AFV expression --------------------------------------------------

deshZTab["ENSCAFG00000032723"] # no difference
deTGFbTab["ENSCAFG00000032723"] # slightly down-regulated


