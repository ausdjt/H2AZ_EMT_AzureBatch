# using Tan et al. 2014 EMT signatures
# supplied as gene symbols :(
sigEMTCells <- readr::read_tsv("~/Data/References/Annotations/Literature/Tan_et_al_2014/Thiery_generic_EMT_sig_cellLine.txt")

# try to map to Ensembl IDs
ensemblHost <- "grch37.ensembl.org"
dataset <- "hsapiens_gene_ensembl"
biomart <- "ensembl"
mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
attribs <- biomaRt::listAttributes(mart)


ensGenesSigEMTCells <- biomaRt::getBM(c("ensembl_gene_id", "external_gene_name"), 
                                         filters = "external_gene_name",
                                         values = sigEMTCells$cellLine_sig,
                                         mart = mart)

ensGenesSigEMTCells <- merge(ensGenesSigEMTCells, 
                              sigEMTCells, 
                              by.x = "external_gene_name", 
                              by.y = "cellLine_sig", 
                              all.x = T)

# get the Cfam homologs
cfamEnsGenesSigEMTCells <- biomaRt::getBM(c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), 
                                            filters = "ensembl_gene_id",
                                            values = c(ensGenesSigEMTCells$ensembl_gene_id),
                                            mart = mart)

cfamEnsGenesSigEMTCells <- cfamEnsGenesSigEMTCells[!cfamEnsGenesSigEMTCells$cfamiliaris_homolog_ensembl_gene == "",]
cfamEnsGenesSigEMTCells <- merge(cfamEnsGenesSigEMTCells, 
                                 ensGenesSigEMTCells, 
                                 by.x = "ensembl_gene_id",
                                 by.y = "ensembl_gene_id",
                                 all.x = T)

cfamEnsGenesSigEMTCells <- cfamEnsGenesSigEMTCells[,-1]
# adding H2AFZ
cfamEnsGenesSigEMTCells <- rbind(cfamEnsGenesSigEMTCells, c("ENSCAFG00000010615", "H2AFZ", NA))
# adding TGFB1
cfamEnsGenesSigEMTCells <- rbind(cfamEnsGenesSigEMTCells, c("ENSCAFG00000005014", "TGFB1", NA))

colnames(cfamEnsGenesSigEMTCells)[1] <- "ensembl_gene_id"
save(cfamEnsGenesSigEMTCells, file = "~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/shiny/MDCK_EMT_RNA-Seq/data/cfamEnsGenesSigEMTCells.rda")
# 

# volcano plot of cell line EMT genes -----------------------------
emtSignatureData <- lapply(names(resultsCompressed[[1]]$sleuth_results.gene), function(y){
  s <- resultsCompressed[[1]]$sleuth_results.gene[[y]]$target_id %in% cfamEnsGenesSigEMTCells$ensembl_gene_id
  dat <- resultsCompressed[[1]]$sleuth_results.gene[[y]][s,]
  dat <- merge(dat,
               ensGenes[,c("ensembl_gene_id", "external_gene_name", "description")], 
               by.x = "target_id", 
               by.y = "ensembl_gene_id")
  dat <- merge(dat, 
               cfamEnsGenesSigEMTCells[,c("ensembl_gene_id", "epi_mes")],
               by.x = "target_id",
               by.y = "ensembl_gene_id")
  dat <- dat[order(dat$qval), ]
  return(list(dataTable = dat))
})
names(emtSignatureData) <- names(resultsCompressed[[1]]$sleuth_results.gene)

y <- "conditionMDCKshZ"
dat <- emtSignatureData[[y]]$dataTable
xAxisMax <- max(abs(dat$b)) + 1
plot(dat$b,
     -log10(dat$qval), 
     axes = F, 
     xlab = "", 
     ylab = "", 
     frame = F,
     cex = 0.3,
     xlim = c(-round(xAxisMax, 0), round(xAxisMax,0)),
     pch = 16, main = paste("Volcano plot\nCondition: ", y, sep = ""))
# points(dat[which(-log10(dat$qval) >= 10), "b"], 
#        -log10(dat[which(-log10(dat$qval) >= 10), "qval"]),
#        col = "red", 
#        pch = 16, 
#        cex = 1.1)
points(dat[dat$epi_mes == "epi","b"],
       -log10(dat[dat$epi_mes == "epi","qval"]),
       col = "blue",
       pch = 16,
       cex = 1)
points(dat[dat$epi_mes == "mes","b"],
       -log10(dat[dat$epi_mes == "mes","qval"]),
       col = "green",
       pch = 16,
       cex = 1)
axis(2,
     pos = 0, 
     lwd = 3)
#                             at = c(seq(0,yAxisMax,10)))
axis(1,
     pos = 0, 
     lwd = 3,
     at = c(seq(-round(xAxisMax, 0), round(xAxisMax,0), 2)))
mtext("-log10(q-value)", 
      side = 2)
mtext("beta value",
      side = 1, 
      line = 2)
# only adding labels to genes with adjuste p-value <= 0.01
text(dat[which(-log10(dat$qval) >= 15), "b"], -log10(dat[which(-log10(dat$qval) >= 15), "qval"]), 
     labels = dat[which(-log10(dat$qval) >= 15), "external_gene_name" ],
     cex = 0.7,
     pos = 4, offset = 0.3)



