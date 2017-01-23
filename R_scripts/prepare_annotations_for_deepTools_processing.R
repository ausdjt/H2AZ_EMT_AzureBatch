# prepare annotation files for deepTools processing
require(GenomicFeatures)
require(GenomicRanges)

cfamEnsemblAnnotation <- "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/Canis_familiaris.CanFam3.1.84.gtf"
cfamEnsemblTxDb <- makeTxDbFromGFF(file = cfamEnsemblAnnotation, 
                                   organism = "Canis familiaris",
                                   dataSource = "Ensembl84")

TxDB <- cfamEnsemblTxDb
# create BED file for all genes -------------------------------------------
gr.genes <- genes(TxDB)
gr.promoters <- promoters(gr.genes,
                          upstream = 2000,
                          downstream = 2000)
deepToolsUtils::WriteGRangesToBED(gr.genes, out_file = "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/allGenes.bed")
deepToolsUtils::WriteGRangesToBED(gr.promoters, out_file = "~/Data/References/Annotations/Canis_familiaris/CanFam3.1_ensembl84/allPromoters_1500bp_up_down.bed")
