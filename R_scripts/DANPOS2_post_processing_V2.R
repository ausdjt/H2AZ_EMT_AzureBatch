# create TxDb object of the Canis familiaris 3.1 genome annotation --------
chromInfo <- read.table("~/mount/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/chromInfo.txt", header = F, as.is = T, sep = "\t")
colnames(chromInfo) <- c("chrom", "length")
TxDb.Cfam3.Ensembl <- makeTxDbFromGFF("~/Data/Annotations/CanFam3/Ensembl/Canis_familiaris.CanFam3.1.83.chr.gtf.gz", 
                                      organism = "Canis familiaris", 
                                      chrominfo = chromInfo)


# extract genomic features, e.g. genes, intergenic regions, promoters --------
Cfam3.genes <- genes(TxDb.Cfam3.Ensembl)
Cfam3.genes.hgnc <- getBM(c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = Cfam3.genes$gene_id, mart = dog)
rownames(Cfam3.genes.hgnc) <- Cfam3.genes.hgnc$ensembl_gene_id
Cfam3.genes$hgnc_symbol <- Cfam3.genes.hgnc[Cfam3.genes$gene_id, ]$hgnc_symbol
save(Cfam3.genes, file = "Cfam3.genes.rda")

#load("~/Data/Tremethick/EMT/GenomeWide/Cfam3.genes.rda")

Cfam3.cds <- cdsBy(TxDb.Cfam3.Ensembl, by = "gene")

Cfam3.intergenic <- gaps(Cfam3.genes)



# check for overlaps ------------------------------------------------------
# load pre-processed DANPOS2 data
load("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/danpos2.anno.rda")
seqlevels(danpos2.anno) <- gsub("chr", "", seqlevels(danpos2.anno))
seqlevels(danpos2.anno)[grep("M", seqlevels(danpos2.anno))] <- "MT"
seqlevels(danpos2.anno, force = T) <- seqlevels(Cfam3.genes)
seqinfo(danpos2.anno, force = T) <- seqinfo(Cfam3.genes)

gr4 <- subsetByOverlaps(danpos2.anno, Cfam3.intergenic)
gr4 <- gr4[which(mcols(gr2)$smt_diff_FDR <= 0.01)]
# create a histogram of the data (here log2 transformed)
df4 <- as(log2(mcols(gr4)[, c("control_smt_val")] + 1), "matrix")
df4 <- rbind(df4, as(log2(mcols(gr4)[, c("treat_smt_val")] + 1), "matrix"))
df4 <- data.frame(df4)
df4$var <- c(rep(paste("Epithelial [N = ", length(gr4), "]", sep = ""), length(gr4)), rep(paste("Mesenchymal [N = ", length(gr4), "]", sep = ""), length(gr4)))
histo4 <- ggplot(df4,aes(x=df4, group=var))
b <- histo4 + geom_boxplot(position = "identity", aes(y = df4, x= var)) + labs(title = "Intergenic regions\nH2A.Z containing nucleosomes,\nsummit value [log2], FDR <= 0.01", x = "Sample", y = "Summit occupation (BG-subtracted) [log2]") + scale_y_continuous(limits=c(4, 16))

grid.newpage()
pushViewport(viewport(layout = Layout))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(a, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(b, newpage = F)

