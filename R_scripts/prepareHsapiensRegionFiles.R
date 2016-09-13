library(biomaRt)
library(GenomicRanges)
library(deepToolsUtils)
library(rtracklayer)
require(biovizBase)

setwd(lDir(pathPrefix, "Data/Tremethick/EMT/ChIP-Seq/GenomeWide/danpos_analysis/"))
qPCRGeneList <- readLines(lDir(pathPrefix, "Data/Tremethick/EMT/ChIP-Seq/MDCK qPCR data/genelist.txt"))
ensemblHost <- "mar2016.archive.ensembl.org"
dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"

mart <- human <- useEnsembl(biomart = biomart, host = ensemblHost, dataset = dataset)
human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
listDatasets(mart)
hsap.qPCRGenesTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = qPCRGeneList, mart = mart)

hsap.qPCRGenesPositions <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), 
                                 filters = "ensembl_gene_id", 
                                 values = hsap.qPCRGenesTab$ensembl_gene_id,
                                 mart = human)
rownames(hsap.qPCRGenesPositions) <- hsap.qPCRGenesPositions$ensembl_gene_id
gr1 <- GRanges(seqnames = hsap.qPCRGenesPositions$chromosome_name, 
               IRanges(hsap.qPCRGenesPositions$start_position, hsap.qPCRGenesPositions$end_position, names = hsap.qPCRGenesPositions$ensembl_gene_id),
               strand = c("+", "-")[match(hsap.qPCRGenesPositions$strand, c(1, -1))],
               hgnc_symbol = hsap.qPCRGenesPositions$hgnc_symbol)

# MCF10A ChIP-Seq data was aligned against UCSC reference, need to add "chr"
seqlevels(gr1) <- paste("chr", seqlevels(gr1), sep = "")

# cross-referencing capture regions with EMT genes
promoterCapRegions <- import("~/Data/Tremethick/Breast/PromoterSeqCap/seqCapTargets_hg38.bed")
subsetByOverlaps(promoterCapRegions, gr1)

deepToolsUtils::WriteGRangesToBED(gr1, out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/EMT_markers.bed")
grlTemp <- flatGrl(GRangesList(sapply(subsetByOverlaps(promoterCapRegions ,gr1), function(x) subsetByOverlaps(gr1, x))))
deepToolsUtils::WriteGRangesToBED(grlTemp, out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/seqCap_Targets_EMT_markers.bed")

# use TGFb-treated MDCK differentially expressed EMT markers
mdckDown <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/qDE.limma.tab.down.TGFb_treated.csv")
mdckUp <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/qDE.limma.tab.up.TGFb_treated.csv")
deepToolsUtils::WriteGRangesToBED(gr1[which(gr1$hgnc_symbol %in% mdckDown$genes)], out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/EMT_markers_down.bed")
deepToolsUtils::WriteGRangesToBED(gr1[which(gr1$hgnc_symbol %in% mdckUp$genes)], out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/EMT_markers_up.bed")

grlTemp <- flatGrl(GRangesList(sapply(subsetByOverlaps(promoterCapRegions,gr1[which(gr1$hgnc_symbol %in% mdckDown$genes)]), function(x) subsetByOverlaps(gr1, x))))
deepToolsUtils::WriteGRangesToBED(grlTemp, out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/seqCap_Targets_EMT_markers_down.bed")
grlTemp <- flatGrl(GRangesList(sapply(subsetByOverlaps(promoterCapRegions,gr1[which(gr1$hgnc_symbol %in% mdckUp$genes)]), function(x) subsetByOverlaps(gr1, x))))
deepToolsUtils::WriteGRangesToBED(grlTemp, out_file = "~/Data/Tremethick/Breast/PromoterSeqCap/processed_data/hg38/deepTools/regionFiles/seqCap_Targets_EMT_markers_up.bed")
