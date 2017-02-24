genome <- BSgenome.Cfamiliaris.UCSC.canFam3
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))
promotersDown <- getSeq(genome, promoters(dataList[["allSamples_TanEMTdown_normal"]][["gr"]], upstream = 500, downstream = 500))
promotersUp <- getSeq(genome, promoters(dataList[["allSamples_TanEMTup_normal"]][["gr"]], upstream = 500, downstream = 500))

writeXStringSet(promotersDown, file = "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/promoters_Tan_down.fa")
writeXStringSet(subseq(promotersDown, start = 1000, end = 1500),
                file = "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/promoters_Tan_down_minus1_nucl.fa")
writeXStringSet(promotersUp, file = "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/promoters_Tan_up.fa")
writeXStringSet(subseq(promotersUp, start = 1500, end = 2500),
                file = "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/promoters_Tan_up_plus1_nucl.fa")

# create set of random promoters as background to test against
cfamEnsemblAnnotation <- "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Canis_familiaris.CanFam3.1.84.gtf"
cfamEnsemblTxDb <- makeTxDbFromGFF(file = cfamEnsemblAnnotation, 
                                   organism = "Canis familiaris",
                                   dataSource = "Ensembl84")

TxDB <- cfamEnsemblTxDb
gr.genes <- genes(TxDB)
seqlevels(gr.genes, force = T) <- seqlevels(genome)
ensGenesIDs <- unique(names(gr.genes))
i1 <- intersect(ensGenesIDs, names(dataList[["allSamples_TanEMTdown_normal"]][["gr"]]))
i1 <- c(i1, intersect(ensGenesIDs, names(dataList[["allSamples_TanEMTup_normal"]][["gr"]])))
table(ensGenesIDs %in% i1)
ensGenesIDs <- ensGenesIDs[!ensGenesIDs %in% i1]
randomPromoters <- sample(ensGenesIDs, size = 2000, replace = F)
randomPromotersSeq <- getSeq(genome, promoters(gr.genes[randomPromoters], upstream = 1500, downstream = 1500))
writeXStringSet(randomPromotersSeq, file = "~/Data/Tremethick/EMT/ChIP-Seq/NB501086_0011_MNekrasov_MDCK_JCSMR_ChIPseq/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/randomPromoters_3kb.fa")
