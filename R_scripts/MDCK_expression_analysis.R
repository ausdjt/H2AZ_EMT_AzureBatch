# TGFb stimulated MDCK cells expression analysis
# GEO dataset http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48525
library("GEOquery")
library("Biobase")

gse48525 <- getGEO(GEO = "GSE48525", GSEMatrix = FALSE)

gsmplatforms <- lapply(GSMList(gse48525),function(x) {Meta(x)$platform})
samples <-lapply(GSMList(gse48525), function(x) {Meta(x)$title})
gsmlist <- Filter(function(gsm) {Meta(gsm)$platform=='GPL17403'},GSMList(gse48525))
probesets <- Table(GPLList(gse48525)[[1]])$ID

#
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) {
  tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
# only retain control and TGFb-treated samples
mdck.expression <- data.matrix[,c(1:4)]
mdck.expression.sd <- apply(mdck.expression, 1, sd)

# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)

# get platform data for probe annotation etc
gpl <- getGEO("GPL17403")
tab <- dataTable(gpl)
gpl17403.annotation <- Table(gpl)[,c("ID", "probeset_id", "seqname", "strand", "start", "stop", "gene_assignment", "mrna_assignment", "category")]

# remove all probes which are lacking genomic locations [not useful in context of H2A.Z analysis]
gpl17403.annotation.genomic <- gpl17403.annotation[which(gpl17403.annotation$category %in% c("main", "rescue", "rrna")), ]
gpl17403.annotation.genomic$start <- as(gpl17403.annotation.genomic$start, "integer")
gpl17403.annotation.genomic$stop <- as(gpl17403.annotation.genomic$stop, "integer")
gpl17403.annotation.genomic$strand <- factor(as.character(gpl17403.annotation.genomic$strand))
gpl17403.annotation.genomic$seqname <- factor(as.character(gpl17403.annotation.genomic$seqname))

gr.gpl17403 <- GRanges(gpl17403.annotation.genomic$seqname, 
                       IRanges(gpl17403.annotation.genomic$start, gpl17403.annotation.genomic$stop), 
                       strand = gpl17403.annotation.genomic$strand, 
                       gpl17403.annotation.genomic[, c("probeset_id", "gene_assignment", "mrna_assignment", "category")])
seqlevels(gr.gpl17403, force = T) <- gsub("chr", "", seqlevels(gr.gpl17403))
gr.gpl17403$gene_assignment <- as(gr.gpl17403$gene_assignment, "character")
gr.gpl17403$mrna_assignment <- as(gr.gpl17403$mrna_assignment, "character")

#
toMatch <- mcols(gr.MSigDB.EMT_associated.cfam)$hgnc_symbol
toMatch <- toMatch[!toMatch == ""]
matches <-  unique(grep(paste(toMatch, collapse = "|"), gpl17403.annotation[,"gene_assignment"], value = F))
MSigDB.EMT_associated.probes <- gpl17403.annotation[matches ,"probeset_id"]

toMatch <- mcols(gr.MSigDB.TGFb_induced_EMT.cfam)$hgnc_symbol
toMatch <- toMatch[!toMatch == ""]
matches <-  unique(grep(paste(toMatch, collapse = "|"), gpl17403.annotation[,"gene_assignment"], value = F))
MSigDB.TGFb_induced_EMT.probes <- gpl17403.annotation[matches, "probeset_id"]

toMatch <- mcols(gr.MSigDB.TGFb_emt_gene_id.cfam)$hgnc_symbol
toMatch <- toMatch[!toMatch == ""]
matches <-  unique(grep(paste(toMatch, collapse = "|"), gpl17403.annotation[,"gene_assignment"], value = F))
MSigDB.EMT.probes <- unique(gpl17403.annotation[matches, "probeset_id"])

