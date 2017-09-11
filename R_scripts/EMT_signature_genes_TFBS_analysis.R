# load libraries
library(GGally)
library(ggplot2)
library(JASPAR2014)
library(Biostrings)
library(BSgenome.Cfamiliaris.UCSC.canFam3)
library(TFBSTools)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(biomaRt)
library(RSQLite)

# create genome object
genome <- BSgenome.Cfamiliaris.UCSC.canFam3
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))

# setup biomaRt object
hsapEnsemblHost <- "grch37.ensembl.org"
hsapDataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
hsapMart <- biomaRt::useEnsembl(biomart = biomart, dataset = hsapDataset, host = hsapEnsemblHost)
attribs <- biomaRt::listAttributes(hsapMart)
# Cfam mart
cfamEnsemblHost <- "uswest.ensembl.org"
cfamDataset <- "cfamiliaris_gene_ensembl"
cfamMart <- biomaRt::useEnsembl(biomart = biomart, dataset = cfamDataset, host = cfamEnsemblHost)


# load EMT genes list -----------------------------------------------------
sigEMTCells <- readr::read_tsv("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")

hsapEnsGenesSigEMTCells <- data.table::data.table(
                                      biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                                                    "external_gene_name",
                                                                    "chromosome_name",
                                                                    "start_position",
                                                                    "end_position",
                                                                    "strand"), 
                                                    filters = "external_gene_name",
                                                    values = c(sigEMTCells$cellLine_sig, "TGFB1", "H2AFZ"),
                                                    mart = hsapMart)                                      )

hsapEnsGenesSigEMTCells <- merge(hsapEnsGenesSigEMTCells, 
                                 sigEMTCells, 
                                 by.x = "external_gene_name", 
                                 by.y = "cellLine_sig", 
                                 all.x = T)

hsapEnsGenesSigEMTCells[is.na(hsapEnsGenesSigEMTCells$epi_mes),]$epi_mes <- hsapEnsGenesSigEMTCells[is.na(hsapEnsGenesSigEMTCells$epi_mes),]$external_gene_name


# collect dog homologs ----------------------------------------------------
cfamEnsGenesSigEMTCells <- data.table::data.table(biomaRt::getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), 
                                                                 filters = "ensembl_gene_id",
                                                                 values = hsapEnsGenesSigEMTCells$ensembl_gene_id,
                                                                 mart = hsapMart))
cfamEnsGenesSigEMTCells <- merge(cfamEnsGenesSigEMTCells, 
                                 hsapEnsGenesSigEMTCells[,c("ensembl_gene_id", "external_gene_name", "epi_mes")], 
                                 by.x = "ensembl_gene_id",
                                 by.y = "ensembl_gene_id",
                                 all.x = T)
cfamEnsGenesSigEMTCells <- cfamEnsGenesSigEMTCells[!cfamEnsGenesSigEMTCells$cfamiliaris_homolog_ensembl_gene == "",]
cfamEnsGenesSigEMTCells <- cfamEnsGenesSigEMTCells[,-1]
# removing one of the dog H2AFZ genes
cfamEnsGenesSigEMTCells <- cfamEnsGenesSigEMTCells[!cfamiliaris_homolog_ensembl_gene == "ENSCAFG00000015774",] 
# adding TGFB1
cfamEnsGenesSigEMTCells <- rbind(cfamEnsGenesSigEMTCells, data.table::data.table(cfamiliaris_homolog_ensembl_gene = "ENSCAFG00000005014",
                                                                                 external_gene_name = "TGFB1",
                                                                                 epi_mes = "TGFB"))

# adding SPP1, FN1, N-Cadherin <-- this has to be done after loading annotation data!
colnames(cfamEnsGenesSigEMTCells)[1] <- "ensembl_gene_id"
ensGenes <- data.table::data.table(ensGenes)
cfamEnsGenesSigEMTCells <- rbind(cfamEnsGenesSigEMTCells, ensGenes[ensembl_gene_id %in% c("ENSCAFG00000009569", "ENSCAFG00000014345", "ENSCAFG00000018115"), c("ensembl_gene_id", "external_gene_name")], fill = T)

cfamEnsGenesSigEMTCells.pos <- data.table::data.table(
                                getBM(attributes = c("ensembl_gene_id",
                                                     "external_gene_name",
                                                     "chromosome_name",
                                                     "start_position",
                                                     "end_position",
                                                     "strand"),
                                      filters = "ensembl_gene_id",
                                      values = cfamEnsGenesSigEMTCells$ensembl_gene_id,
                                      mart = cfamMart))
gr.cfamEnsGenesSigEMTCells <- GenomicRanges::GRanges(seqnames = cfamEnsGenesSigEMTCells.pos$chromosome_name,
                                                     IRanges(start = cfamEnsGenesSigEMTCells.pos$start_position,
                                                             end = cfamEnsGenesSigEMTCells.pos$end_position,
                                                             names = cfamEnsGenesSigEMTCells.pos$ensembl_gene_id),
                                                     strand = c("+", "-")[match(cfamEnsGenesSigEMTCells.pos$strand, c(1,-1))],
                                                     external_gene_name = cfamEnsGenesSigEMTCells.pos$external_gene_name)

gr.cfamEnsGenesSigEMTCells.promoters <- GenomicRanges::promoters(gr.cfamEnsGenesSigEMTCells,
                                                                 upstream = 1500,
                                                                 downstream = 1500)

gr.cfamEnsGenesSigEMTCells.promoters$external_gene_name <- gr.cfamEnsGenesSigEMTCells$external_gene_name
grl.cfamEnsGenesSigEMTCells.promoters <- as(gr.cfamEnsGenesSigEMTCells.promoters, "GRangesList")

gr.cfamEnsGenesSigEMTCells.promoters.up <- GenomicRanges::promoters(gr.cfamEnsGenesSigEMTCells,
                                                                 upstream = 1500,
                                                                 downstream = 0)

# GRanges of human EMT genes for searching for motifs ---------------------
gr.hsapEnsGenesSigEMTCells <- GenomicRanges::GRanges(seqnames = hsapEnsGenesSigEMTCells$chromosome_name,
                                                     IRanges(start = hsapEnsGenesSigEMTCells$start_position,
                                                             end = hsapEnsGenesSigEMTCells$end_position,
                                                             names = hsapEnsGenesSigEMTCells$ensembl_gene_id),
                                                     strand = c("+", "-")[match(hsapEnsGenesSigEMTCells$strand, c(1,-1))],
                                                     marker = hsapEnsGenesSigEMTCells$epi_mes,
                                                     external_gene_name = hsapEnsGenesSigEMTCells$external_gene_name,
                                                     ensemblStrand = hsapEnsGenesSigEMTCells$strand)

# prune out non-canonical chromosomes
seqlevels(gr.hsapEnsGenesSigEMTCells, pruning.mode = "coarse") <- seqlevels(gr.hsapEnsGenesSigEMTCells)[grep("H", seqlevels(gr.hsapEnsGenesSigEMTCells), invert = T)]
gr.hsapEnsGenesSigEMTCells.promoters <- GenomicRanges::promoters(gr.hsapEnsGenesSigEMTCells,
                                                                 upstream = 1500,
                                                                 downstream = 0)


# get motifs found in the promoter regions of EMT signature genes ---------
ensembl_regulation = biomaRt::useMart(biomart="ENSEMBL_MART_FUNCGEN",
                                      host=ensemblHost,
                                      dataset="hsapiens_motif_feature")
#listDatasets(useMart(biomart="ENSEMBL_MART_FUNCGEN",host=ensemblHost))
#listFilters(ensembl_regulation)
attribs.ensReg <- biomaRt::listAttributes(ensembl_regulation)
searchValues <- paste(seqnames(gr.hsapEnsGenesSigEMTCells.promoters), 
                      start(gr.hsapEnsGenesSigEMTCells.promoters), 
                      end(gr.hsapEnsGenesSigEMTCells.promoters), 
                      gr.hsapEnsGenesSigEMTCells.promoters$ensemblStrand,
                      sep = ":")

hsapEnsRegSigEMTCells.Motifs <- data.table::data.table(
                                  biomaRt::getBM(attributes=c('binding_matrix_id',
                                                              'chromosome_name', 
                                                              'chromosome_start', 
                                                              'chromosome_end',
                                                              'chromosome_strand',
                                                              'display_label',
                                                              'feature_type_name',
                                                              'score',
                                                              'so_name',
                                                              'so_accession'), 
                                                filters='chromosomal_region',
                                                values= searchValues, 
                                                mart=ensembl_regulation))

write.csv(hsapEnsRegSigEMTCells.Motifs, file = "hsapEnsRegSigEMTCells.Motifs.csv")

gr.hsapEnsRegSigEMTCells.Motifs <- GenomicRanges::GRanges(seqnames = hsapEnsRegSigEMTCells.Motifs$chromosome_name,
                                                          IRanges(start = hsapEnsRegSigEMTCells.Motifs$chromosome_start,
                                                                  end = hsapEnsRegSigEMTCells.Motifs$chromosome_end,
                                                                  names = hsapEnsRegSigEMTCells.Motifs$display_label),
                                                          strand = c("+", "-")[match(hsapEnsRegSigEMTCells.Motifs$chromosome_strand, c(1,-1))],
                                                          binding_matrix_id = hsapEnsRegSigEMTCells.Motifs$binding_matrix_id,
                                                          feature_type_name = hsapEnsRegSigEMTCells.Motifs$feature_type_name,
                                                          score = hsapEnsRegSigEMTCells.Motifs$score)


# check expression of TFs in Cfam -----------------------------------------
# now I have the TFs which are actually of relevance in human context
# and have found which ones of them have binding sites in the promoters of EMT genes in dog
# now check their expression so that we can narow down which ones are of interest
# 1) get ENS gene IDs and extract dog homologs
# cfam mart
cfamEnsemblHost <- "uswest.ensembl.org"
cfamDataset <- "cfamiliaris_gene_ensembl"
cfamMart <- biomaRt::useEnsembl(biomart = biomart, dataset = cfamDataset, host = cfamEnsemblHost)

motifGeneNames <- unique(unlist(lapply(strsplit(unique(names(gr.hsapEnsRegSigEMTCells.Motifs)), ":"), function(x) x[1])))

hsapEnsRegSigEMTCells.TFs <- data.table::data.table(
                                biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                              "external_gene_name",
                                                              "cfamiliaris_homolog_ensembl_gene"),
                                               filters = "external_gene_name",
                                               values = unique(c(motifGeneNames, "JUN", "FOS", "MYC")),
                                               mart = hsapMart))


cfamEnsRegSigEMTCells.TFs <- data.table::data.table(
                                getBM(attributes = c("ensembl_gene_id",
                                                     "external_gene_name"),
                                      filters = "ensembl_gene_id",
                                      values = hsapEnsRegSigEMTCells.TFs$cfamiliaris_homolog_ensembl_gene,
                                      mart = cfamMart))

write.csv(cfamEnsRegSigEMTCells.TFs, "cfamEnsRegSigEMTCells.TFs.csv")

plotData <- data.table(results[["MDCK"]]$kallisto_table_genes)
plotData$condition <- factor(plotData$condition,levels(plotData$condition)[c(1,3,2)])

plotDataSubset <- plotData[target_id %in% cfamEnsRegSigEMTCells.TFs$ensembl_gene_id]
plotDataSubset <- merge(plotDataSubset, cfamEnsRegSigEMTCells.TFs, by.x = "target_id", by.y = "ensembl_gene_id", all.x = T, all.y = F)

bp1 <- ggplot(data = plotDataSubset, aes(x = external_gene_name, y = tpm, fill = condition)) + 
              geom_boxplot() +
              theme(axis.text.x = element_text(angle = 45))
bp1

# look at differential expression of the identified TFs
diffExpTGFb <- data.table::data.table(resultsCompressed[["MDCK"]]$sleuth_results.gene$conditionMDCKTGFb)
data.table::setkey(diffExpTGFb, "target_id")
diffExpTGFb.TFs <- diffExpTGFb[cfamEnsRegSigEMTCells.TFs$ensembl_gene_id]
diffExpTGFb.TFs <- merge(diffExpTGFb.TFs, cfamEnsRegSigEMTCells.TFs, by.x = "target_id", by.y = "ensembl_gene_id")
table(diffExpTGFb.TFs$qval < 0.1)
diffExpTGFb.TFs[order(diffExpTGFb.TFs$b),]
upTFs <- diffExpTGFb.TFs[diffExpTGFb.TFs$b > 0 & diffExpTGFb.TFs$qval < 0.1]$external_gene_name
downTFs <- diffExpTGFb.TFs[diffExpTGFb.TFs$b < 0 & diffExpTGFb.TFs$qval < 0.1]$external_gene_name
nsTFs <- diffExpTGFb.TFs[diffExpTGFb.TFs$qval > 0.1]$external_gene_name

# create heatmap of the TPMs, ordered by differential expression
heatmapData <- plotDataSubset[,mean(tpm, na.rm = TRUE), by = c("condition", "external_gene_name")]
heatmapData <- heatmapData[condition != "MDCKshZ"]
colnames(heatmapData)[3] <- "mean_tpm"
heatmapData$external_gene_name <- factor(heatmapData$external_gene_name, levels = rev(diffExpTGFb.TFs[order(diffExpTGFb.TFs$b),]$external_gene_name))
heatmapData$change <- as.character("NA")
setkey(heatmapData, "external_gene_name")
heatmapData[upTFs, ]$change <- "up"
setkey(heatmapData, "external_gene_name")
heatmapData[downTFs, ]$change <- "down"
setkey(heatmapData, "external_gene_name")
heatmapData[nsTFs, ]$change <- "not significant"
setkey(heatmapData, "external_gene_name")
heatmapData <- heatmapData[!"CTCFL"]  
hm1 <- ggplot(data = heatmapData, aes(x = condition, y = external_gene_name)) + geom_tile(aes(fill = mean_tpm), colour = "white") +
              scale_fill_gradient(low = "white", high = "red", name = "Mean TPM") +
              theme(axis.text.x = element_text(angle = 0, vjust = 0.75)) +
              ylab("TF Gene Name") +
              xlab("Condition") +
              facet_grid(change ~ ., scales = "free", shrink = F, space = "free")
hm1
ggsave("TF_heatmap.pdf", hm1, height = 108, width = 196, units = "mm", useDingbats = F)

# perform scanning of Cfam EMT gene promoters for TFBS --------------------
opts <- list()
opts[["ID"]] = unique(gr.hsapEnsRegSigEMTCells.Motifs$binding_matrix_id)
opts[["name"]] = unlist(lapply(unique(hsapEnsRegSigEMTCells.Motifs$binding_matrix_id), function(x) unique(hsapEnsRegSigEMTCells.Motifs[binding_matrix_id == x]$feature_type_name)))
opts[["all_versions"]] = FALSE
PFMatrixList <- getMatrixSet(JASPAR2014, opts = opts)
pwmList <- do.call(PWMatrixList, lapply(PFMatrixList, toPWM))


# create DNAStringSet for the EMT gene promoters --------------------------
dss <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters))
dssEpiDown <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters[epiDown$target_id]))
dssMesDown <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters[mesDown$target_id]))
dssEpiUp <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters[epiUp$target_id]))
dssMesUp <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters[mesUp$target_id]))
TGFb <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters["ENSCAFG00000005014"]))
H2AZ <- DNAStringSet(Views(genome, gr.cfamEnsGenesSigEMTCells.promoters["ENSCAFG00000010615"]))


sitesetList <- TFBSTools::searchSeq(pwmList, dss, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListEpiDown <- TFBSTools::searchSeq(pwmList, dssEpiDown, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListMesDown <- TFBSTools::searchSeq(pwmList, dssMesDown, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListEpiUp <- TFBSTools::searchSeq(pwmList, dssEpiUp, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListMesUp <- TFBSTools::searchSeq(pwmList, dssMesUp, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListTFGb <- TFBSTools::searchSeq(pwmList, TGFb, seqname = as(x, "character"), min.score = "95%", strand = "*") 
sitesetListH2AZ <- TFBSTools::searchSeq(pwmList, H2AZ, seqname = as(x, "character"), min.score = "95%", strand = "*") 

gr.sitesetListEpiDown <- as(sitesetListEpiDown, "GRanges")
gr.sitesetListMesDown <- as(sitesetListMesDown, "GRanges")
gr.sitesetListEpiUp <- as(sitesetListEpiUp, "GRanges")
gr.sitesetListMesUp <- as(sitesetListMesUp, "GRanges")
gr.sitesetListTFGb <- as(sitesetListTFGb, "GRanges")
gr.sitesetListH2AZ <- as(sitesetListH2AZ, "GRanges")

hist(unlist(pvalues(sitesetListTFGb)))
hist(unlist(pvalues(sitesetListH2AZ)))
intersect(as.character(unique(gr.sitesetListTFGb$TF)), as.character(unique(gr.sitesetListH2AZ$TF)))
# try some parallel coordinates plotting
# first summarize TPMs of all three replicats to single mean per gene
pacoData <- plotData[target_id %in% cfamEnsRegSigEMTCells.TFs$ensembl_gene_id, mean(tpm, na.rm = TRUE),by = c("target_id", "condition")]
pacoData <- spread(data = pacoData, key = condition, value = V1)
GGally::ggparcoord(pacoData, columns = c(2:4), scale = "robust")


# get epigenetic regulatory features --------------------------------------
ensembl_regulation = useMart(biomart="ENSEMBL_MART_FUNCGEN",
                             host=ensemblHost,
                             dataset="hsapiens_regulatory_feature")
listAttributes(ensembl_regulation)
listFilters(ensembl_regulation)
searchValues <- paste(seqnames(gr.promoters.hsapEnsGenesSigEMTCells), 
                      start(gr.promoters.hsapEnsGenesSigEMTCells), 
                      end(gr.promoters.hsapEnsGenesSigEMTCells), 
                      sep = ":")
ensRegSigEMTCells.RegFeat <- getBM(attributes=c('activity',
                                                'regulatory_stable_id',
                                                'bound_seq_region_start',
                                                'bound_seq_region_end',
                                                'chromosome_name', 
                                                'chromosome_start', 
                                                'chromosome_end',
                                                'feature_type_name',
                                                'feature_type_description',
                                                'epigenome_name'), 
                                   filters='chromosomal_region',
                                   values= searchValues, 
                                   mart=ensembl_regulation)
