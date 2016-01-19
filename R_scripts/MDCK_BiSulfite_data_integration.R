# integration of MDCK bi-sulfite sequencing data
#----------------Integrate Bi-Sulfite sequencing data--------------------------
library("GEOquery")
GEO_MDCK_Methylation_GSE43697 <- getGEO("GSE43697", GSEMatrix = F) #data comes us supplementary files, downloaded separately
path_gse43697 <- "/Volumes/gduserv/Data/GEO/GSE43697/"

GSM1068566_CpG_context_MDCK_control <- read.table(paste(path_gse43697, "GSM1068566_CpG_context_MDCK_control.pileup.txt", sep = ""), header = F, as.is = T, sep = "\t")
colnames(GSM1068566_CpG_context_MDCK_control) <- c("Chomosome", "position", "methylated_C", "unmethylated_C")
GSM1068566_CpG_context_MDCK_control <- GSM1068566_CpG_context_MDCK_control[-47583803,]
GSM1068566_CpG_context_MDCK_control$perc_methylated <-  GSM1068566_CpG_context_MDCK_control$methylated_C / (GSM1068566_CpG_context_MDCK_control$unmethylated_C + GSM1068566_CpG_context_MDCK_control$methylated_C) * 100
gr.GSM1068566_CpG_context_MDCK_control <- GRanges(GSM1068566_CpG_context_MDCK_control$Chomosome, 
                                                  IRanges(start = as(GSM1068566_CpG_context_MDCK_control$position, "integer"), width = 1),
                                                  perc_methylated = GSM1068566_CpG_context_MDCK_control$perc_methylated)

GSM1068567_CpG_context_MDCK_TGFB <- read.table(paste(path_gse43697, "GSM1068567_CpG_context_MDCK_TGFB.pileup.txt", sep = ""), header = F, as.is = T, sep = "\t")
colnames(GSM1068567_CpG_context_MDCK_TGFB) <- c("Chomosome", "position", "methylated_C", "unmethylated_C")
GSM1068567_CpG_context_MDCK_TGFB <- GSM1068567_CpG_context_MDCK_TGFB[-46385053,]
GSM1068567_CpG_context_MDCK_TGFB$perc_methylated <-  GSM1068567_CpG_context_MDCK_TGFB$methylated_C / (GSM1068567_CpG_context_MDCK_TGFB$unmethylated_C + GSM1068567_CpG_context_MDCK_TGFB$methylated_C) * 100
gr.GSM1068567_CpG_context_MDCK_TGFB <- GRanges(GSM1068567_CpG_context_MDCK_TGFB$Chomosome, 
                                               IRanges(start = as(GSM1068567_CpG_context_MDCK_TGFB$position, "integer"), width = 1),
                                               perc_methylated = GSM1068567_CpG_context_MDCK_TGFB$perc_methylated)
# now we are only looking for methylation across the 10 markers
gr.which.ucsc <- gr.which
seqlevels(gr.which.ucsc) <- paste("chr", seqlevels(gr.which), sep = "")

gr.GSM1068567_CpG_context_MDCK_TGFB <- subsetByOverlaps(gr.GSM1068567_CpG_context_MDCK_TGFB, gr.which.ucsc)
gr.GSM1068566_CpG_context_MDCK_control <- subsetByOverlaps(gr.GSM1068566_CpG_context_MDCK_control, gr.which.ucsc)
gr.intersect <- subsetByOverlaps(gr.GSM1068567_CpG_context_MDCK_TGFB, gr.GSM1068566_CpG_context_MDCK_control)
colnames(mcols(gr.intersect))[1] <- "perc_methylated_TGFb"
mcols(gr.intersect)$perc_methylated_WT <- mcols(subsetByOverlaps(gr.GSM1068566_CpG_context_MDCK_control, gr.intersect))[,1]
dT.methylation <- DataTrack(gr.intersect, type = "p", groups = c("TGFb","WT"), legend = TRUE, name = "DNA methylation")
displayPars(dT.methylation) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")

# perhaps better to use the published Differentially Methylated Regions [Carmona et al 2014]
mdck_dmr <- read.csv("MDCK_WGBS_Carmona_et_al_2014/suppTable2.csv", strip.white = TRUE)
colnames(mdck_dmr) <- c("chrom", "start", "end", "gene_symbol", "ENSEMBL_transcript_id", "No_of_CpG_sites", "average_DNA_methylation_difference")
#mdck_dmr$chrom <- paste("chr", mdck_dmr$chrom, sep = "")
gr.mdck_dmr <- GRanges(mdck_dmr$chrom, IRanges(mdck_dmr$start, mdck_dmr$end), strand = "*", mdck_dmr[, c("gene_symbol", "ENSEMBL_transcript_id", "No_of_CpG_sites", "average_DNA_methylation_difference")])
dT.dmr <- DataTrack(gr.mdck_dmr, type = "histogram", name = "MDCK DMRs", data = mcols(gr.mdck_dmr)$average_DNA_methylation_difference)
displayPars(dT.dmr) <- list("fontcolor.title" = "black", "background.title" = "white", "col.axis" = "black", "col.frame" = "white")
