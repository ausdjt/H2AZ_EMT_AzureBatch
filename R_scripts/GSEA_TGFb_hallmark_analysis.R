# Analysis of H2A.Z occupation at TSS of TGFb induced EMT hallmark genes
# Geneset from MSigDB

# MSigDB data extraction
library("XML")

# downloaded msigbdb XML file from Broad Institute 2015-09-30
# parse the XML file
doc <- xmlTreeParse("~/Data/Annotations/MSigDB/msigdb_v5.0.xml", useInternalNodes = FALSE)
root <- xmlRoot(doc)
child <- xmlChildren(root)

# get STANDARD_NAME & DESCRIPTION_BRIEF attributes of all child nodes
child.standard_name <- sapply(child, function(x) xmlGetAttr(x, name = "STANDARD_NAME"))
child.description_brief <- sapply(child, function(x) xmlGetAttr(x, name = "DESCRIPTION_BRIEF"))
child.organism <- sapply(child, function(x) xmlGetAttr(x, name = "ORGANISM")) 

names(child.standard_name) <- NULL
names(child.description_brief) <- NULL
names(child.organism) <- NULL

child.df <- data.frame(cbind(as(child.standard_name, "character"), as(child.description_brief, "character"), as(child.organism, "character")))
child.df.emt <- child.df[grep("chymal", child.df[,2]),]
child.df.emt <- child.df.emt[grep("sapiens", child.df.emt[,3]),]

# only the "Cancer Hallmark Gene Set"
MSigDB.hallmark_set_emt_gene_id <- unlist(strsplit(sapply(child[c(10328)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ","))
save(MSigDB.TGFb_emt_gene_id, file = "MSigDB.TGFb_emt_gene_id.rda")
# get Cfam homologs of EMT hallmark genes
MSigDB.TGFb_emt_gene_id.cfam <- getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"),
                                      filters = "entrezgene", 
                                      values = as(MSigDB.hallmark_set_emt_gene_id, "character"), 
                                      mart = human)
MSigDB.TGFb_emt_gene_id.cfam <- getBM(attributes = c("ensembl_gene_id",
                                                     "entrezgene",
                                                     "chromosome_name",
                                                     "start_position",
                                                     "end_position",
                                                     "hgnc_symbol",
                                                     "strand"), 
                                      filters = "ensembl_gene_id", 
                                      values = unique(as(MSigDB.TGFb_emt_gene_id.cfam$cfamiliaris_homolog_ensembl_gene, "character")), mart = dog)
MSigDB.TGFb_emt_gene_id.cfam$set <- "Hallmarks_EMTS"
gr.MSigDB.TGFb_emt_gene_id.cfam <- GRanges(MSigDB.TGFb_emt_gene_id.cfam$chromosome_name, 
                                           IRanges(MSigDB.TGFb_emt_gene_id.cfam$start_position, MSigDB.TGFb_emt_gene_id.cfam$end_position), 
                                           strand = c("-", "+")[match(MSigDB.TGFb_emt_gene_id.cfam$strand, c("-1", "1"))],
                                           MSigDB.TGFb_emt_gene_id.cfam[, c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "set")])

gr.MSigDB.TGFb_emt_gene_id.cfam.tss <- promoters(gr.MSigDB.TGFb_emt_gene_id.cfam, upstream = 1, downstream = 0)
gr.MSigDB.TGFb_emt_gene_id.cfam.tss1500 <- promoters(gr.MSigDB.TGFb_emt_gene_id.cfam, upstream = 1500, downstream = 1500)

# TGFb-induced EMT
MSigDB.TGFb_induced_EMT <- unlist(strsplit(sapply(child[c(4242)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ","))
MSigDB.TGFb_induced_EMT.cfam <- getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), filters = "entrezgene", values = as(MSigDB.TGFb_induced_EMT, "character"), mart = human)
MSigDB.TGFb_induced_EMT.cfam <- getBM(attributes = c("ensembl_gene_id",
                                                     "entrezgene",
                                                     "chromosome_name",
                                                     "start_position",
                                                     "end_position",
                                                     "hgnc_symbol",
                                                     "strand"), 
                                      filters = "ensembl_gene_id", 
                                      values = unique(as(MSigDB.TGFb_induced_EMT.cfam$cfamiliaris_homolog_ensembl_gene, "character")), mart = dog)
MSigDB.TGFb_induced_EMT.cfam$set <- "TGFb_induced_EMT"
gr.MSigDB.TGFb_induced_EMT.cfam <- GRanges(MSigDB.TGFb_induced_EMT.cfam$chromosome_name, 
                                           IRanges(MSigDB.TGFb_induced_EMT.cfam$start_position, MSigDB.TGFb_induced_EMT.cfam$end_position), 
                                           strand = c("-", "+")[match(MSigDB.TGFb_induced_EMT.cfam$strand, c("-1", "1"))],
                                           MSigDB.TGFb_induced_EMT.cfam[, c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "set")])

gr.MSigDB.TGFb_induced_EMT.cfam.tss <- promoters(gr.MSigDB.TGFb_induced_EMT.cfam, upstream = 1, downstream = 0)
gr.MSigDB.TGFb_induced_EMT.cfam.tss1500 <- promoters(gr.MSigDB.TGFb_induced_EMT.cfam, upstream = 1500, downstream = 1500)

# General EMT-associated genes
MSigDB.EMT_associated <- unlist(strsplit(sapply(child[c(4242, 505, 7217, 7218)], function(x) xmlGetAttr(x, name = "MEMBERS_EZID")), ","))
MSigDB.EMT_associated <- unique(MSigDB.EMT_associated)
MSigDB.EMT_associated.cfam <- getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), filters = "entrezgene", values = as(MSigDB.EMT_associated, "character"), mart = human)
MSigDB.EMT_associated.cfam <- getBM(attributes = c("ensembl_gene_id",
                                                     "entrezgene",
                                                     "chromosome_name",
                                                     "start_position",
                                                     "end_position",
                                                     "hgnc_symbol",
                                                     "strand"), 
                                      filters = "ensembl_gene_id", 
                                      values = unique(as(MSigDB.EMT_associated.cfam$cfamiliaris_homolog_ensembl_gene, "character")), mart = dog)
MSigDB.EMT_associated.cfam$set <- "EMT_associated"
gr.MSigDB.EMT_associated.cfam <- GRanges(MSigDB.EMT_associated.cfam$chromosome_name, 
                                           IRanges(MSigDB.EMT_associated.cfam$start_position, MSigDB.EMT_associated.cfam$end_position), 
                                           strand = c("-", "+")[match(MSigDB.EMT_associated.cfam$strand, c("-1", "1"))],
                                           MSigDB.EMT_associated.cfam[, c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "set")])

gr.MSigDB.EMT_associated.cfam.tss <- promoters(gr.MSigDB.EMT_associated.cfam, upstream = 1, downstream = 0)
gr.MSigDB.EMT_associated.cfam.tss1500 <- promoters(gr.MSigDB.EMT_associated.cfam, upstream = 1500, downstream = 1500)

l1 <- lapply(as.integer(rownames(child.df.emt)), function(x) unique(unlist(strsplit(sapply(child[x], function(y) xmlGetAttr(y, name = "MEMBERS_EZID")), ","))))
names(l1) <- child.df.emt$X1
list.cfamId <- lapply(l1, function(x) {
  bm1 <- getBM(attributes = c("ensembl_gene_id", 
                              "cfamiliaris_homolog_ensembl_gene"), 
               filters = "entrezgene", 
               values = as(x, "character"), mart = human)
  bm1 <- getBM(attributes = c("ensembl_gene_id",
                               "entrezgene",
                               "chromosome_name",
                               "start_position",
                               "end_position",
                               "hgnc_symbol",
                               "strand"), 
                filters = "ensembl_gene_id", 
                values = unique(as(bm1$cfamiliaris_homolog_ensembl_gene, "character")), mart = dog)
})
df <- do.call(rbind.data.frame, list.cfamId)
df$set <- unlist(lapply(strsplit(rownames(df), "\\."), function(x) x[1]))
rownames(df) <- c(1:nrow(df))
df <- df[order(df$ensembl_gene_id), ]

df.unique <- df[-which(duplicated(df$ensembl_gene_id)), ]

gr.MSigDB.EMT_associated.cfam <- GRanges(df.unique$chromosome_name, 
                                         IRanges(df.unique$start_position, df.unique$end_position), 
                                         strand = c("-", "+")[match(df.unique$strand, c("-1", "1"))],
                                         df.unique[, c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "set")])




