require(dtplyr)

qDE.TGFb.limma.tab <- read.csv("~/Data/Tremethick/EMT/qPCR_data_analysis/qDE.limma.tab.TGFb_treated.csv", row.names = 1)
qDE.shZKD.limma.tab <- read.csv("~/Data/Tremethick/EMT/MDCK qPCR data/H2AZ-knockdown/qDE.limma.tab.H2AZ_KD.csv")

RNASeqMDCKshZ <- as.data.table(resultsCompressed[["MDCK"]][["sleuth_results.gene"]]$conditionMDCKshZ)
RNASeqMDCKTGFb <- as.data.table(resultsCompressed[["MDCK"]][["sleuth_results.gene"]]$conditionMDCKTGFb)
qPCRMDCKTGFb <- as.data.table(qDE.TGFb.limma.tab)
qPCRMDCKshZ <- as.data.table(qDE.shZKD.limma.tab)

qPCRGenesTab <- as.data.table(qPCRGenesTab)
length(qPCRMDCKTGFb$genes)

plotENSGIDs <- filter(RNASeqMDCKshZ, target_id %in% filter(qPCRGenesTab, hgnc_symbol %in% qPCRMDCKTGFb$genes)$ensembl_gene_id)$target_id
plotGeneNames <- filter(qPCRGenesTab, ensembl_gene_id %in% plotENSGIDs)$hgnc_symbol
plotENSGIDs <- filter(qPCRGenesTab, ensembl_gene_id %in% plotENSGIDs)$ensembl_gene_id

setkey(RNASeqMDCKTGFb, target_id)
setkey(qPCRMDCKTGFb, genes)
cor(RNASeqMDCKTGFb[plotENSGIDs]$b, qPCRMDCKTGFb[plotGeneNames]$logFC)

setkey(qPCRMDCKshZ, genes)
setkey(RNASeqMDCKshZ, target_id)
cor(RNASeqMDCKshZ[plotENSGIDs]$b, qPCRMDCKshZ[plotGeneNames]$logFC)
