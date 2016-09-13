require(deepToolsUtils)
require(biomaRt)
require(tidyr)
require(RColorBrewer)
require(tibble)
require(snowfall)
require(sleuth)
require(tibble)
require(tximport)
require(readr)

# external functions ------------------------------------------------------
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")

# local functions ---------------------------------------------------------
lDir <- function(x, y){
  paste(x, y, sep = "/")
}

# global variables --------------------------------------------------------
ensemblHost <- "uswest.ensembl.org"
dataset <- "cfamiliaris_gene_ensembl"
biomart <- "ensembl"
colors <- RColorBrewer::brewer.pal(3, "Set2")

if (amILocal("JCSMR027564ML")){
  pathPrefix <- "~/mount/gduserv"
  cpus <- 8
} else {
  pathPrefix <- "~"
  cpus <- 16
  options(width = 137)
}
options(mc.cores = cpus)
setwd(lDir(pathPrefix, 
           "Data/Tremethick/EMT/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/R_Analysis/"))
dataPath <- lDir(pathPrefix, 
                 "Data/Tremethick/EMT/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data/CanFam3.1_ensembl84_ERCC/HTSeq/count/")
devPath <- "~/Development"

files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)

# get snakemake run configuration -----------------------------------------
runConfig <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json")
runConfig$references$CanFam3.1$version

# preparing annotation data from Ensembl ----------------------------------
if (!file.exists("ensGenes.rda")){
  mart <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, host = ensemblHost)
  attribs <- biomaRt::listAttributes(mart)
  ensGenes <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                              "external_gene_name",
                              "chromosome_name",
                              "start_position",
                              "end_position",
                              "strand",
                              "band",
                              "description",
                              "percentage_gc_content",
                              "gene_biotype"),
                              mart = mart)
  save(ensGenes, file = "ensGenes.rda")

  # get Ensembl transcripts
  ensTranscripts <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                                  "ensembl_gene_id",
                                                  "transcript_length",
                                                  "version", 
                                                  "transcript_version",
                                                  "external_gene_name"),
                                   mart = mart,
                                   filter = "ensembl_gene_id",
                                   values = ensGenes$ensembl_gene_id)
  save(ensTranscripts, file = "ensTranscripts.rda")
  
  # create t2g object
  t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                            "ensembl_gene_id", 
                            "external_gene_name", 
                            "version", 
                            "transcript_version")]
  t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  save(t2g, file = "t2g.rda")
  
  mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
    y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
    y <- y[which.max(y$transcript_length), ]$transcript_length})
  save(mylength, file = "mylength.rda")
  mygc <- ensGenes$percentage_gc_content
  names(mygc) <- ensGenes$ensembl_gene_id
  save(mygc, file = "mygc.rda")
  mybiotypes <- ensGenes$gene_biotype
  names(mybiotypes) <- ensGenes$ensembl_gene_id
  save(mybiotypes, file = "mybiotypes.rda")
  mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
  save(mychroms, file = "mychroms.rda")
  } else {
  load("ensGenes.rda")
  load("ensTranscripts.rda")
  load("mylength.rda")
  load("mygc.rda")
  load("mybiotypes.rda")
  load("mychroms.rda")
  load("t2g.rda")
}

# load kallisto data with tximport and inspect via PCA -------------------------
base_dir <- paste(pathPrefix, 
                "Data/Tremethick/EMT/RNA-Seq/NB501086_0067_RDomaschenz_JCSMR_RNASeq/processed_data", runConfig$references$CanFam3.1$version, "kallisto", sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) paste(x[1:2], collapse = "_")))
files <- paste(kal_dirs, "abundance.tsv", sep = "/")
names(files) <- sample_id
txi <- tximport::tximport(files, 
                          type = "kallisto",
                          geneIdCol = "ens_gene",
                          txIdCol = "target_id",
                          tx2gene = t2g,
                          reader = read_tsv)

# exploratory analysis of the abundance data ------------------------------
# perform PCA for first inspection of data
pca1 <- ade4::dudi.pca(t(txi$abundance), scannf = F, nf = 6)
pdf("Exploratory_Analysis_sclass_plot.pdf")
ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()
pdf("Exploratory_Analysis_sarrow_plot.pdf")
ade4::s.arrow(pca1$li, clabel = 0.7)
dev.off()

cor(txi$abundance)
pairs(txi$abundance)
sd1 <- apply(txi$abundance, 1, sd)
pdf("Exploratory_Analysis_Heatmap.pdf")
heatmap.3(log2(txi$abundance[sd1 > 15,] + 1),
          trace = "none", 
          distfun = function(x) dist(x, method = "manhattan"), 
          hclustfun = function(y) hclust(y,method="ward.D"),
          cexCol = 0.6,
          labRow = NA,
          main = "Exploratory analysis on abundance estimates")
dev.off()


################################################################################
# IMPORTANT!!!!!
# PCA shows that the knockdown in sample MDCK_shZ_rep1 probably did not work 
# removing it from list
# clustering also shows that MDCK_TGFb_rep1 is actually more similar to
# MDCK_shZ_rep2 than to MDCK_TGFb_rep2
condition <- as.factor(condition)
condition <- relevel(condition, "MDCK_wt")
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c <- s2c[-1, ]
s2c <- s2c[c("5", "6", "3", "4", "2"),]
s2c$sample <- as.character(s2c$sample)
s2c.list <- list(MDCK = s2c)
################################################################################

# actual processing using sleuth------------------------------------------------
if (!file.exists("MDCK_kallisto_analysis_results.rda")){
  results <- lapply(names(s2c.list), function(x){
    print(paste("Processing ", x, sep = ""))
    design <- model.matrix(~ condition, data = s2c.list[[x]])
    #---------------------------------------------------------------------------
    # transcript-level DE
    so <- sleuth::sleuth_prep(s2c.list[[x]], ~ condition, target_mapping = t2g)
    so <- sleuth::sleuth_fit(so, formula = design)
    so <- sleuth::sleuth_fit(so, ~1, "reduced")
    so <- sleuth::sleuth_lrt(so, "reduced", "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so <- sleuth::sleuth_wt(so, i)  
    }
    rt.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt <- sleuth::sleuth_results(so, x)
      rt <- rt[order(rt$qval),]
    })
    names(rt.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    kt <- sleuth::kallisto_table(so, include_covariates = T)
    kt <- sleuth::kallisto_table(so, normalized = T, include_covariates = T)
    kt_wide <- tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm)
    rownames(kt_wide) <- kt_wide[,1]
    kt_wide <- kt_wide[,-1]
    #---------------------------------------------------------------------------
    # gene-level DE  
    so.gene <- sleuth::sleuth_prep(s2c.list[[x]], 
                                   ~ condition, 
                                   target_mapping = t2g, 
                                   aggregation_column = "ens_gene")
    so.gene <- sleuth::sleuth_fit(so.gene, 
                                  formula = design)
    so.gene <- sleuth::sleuth_fit(so.gene,
                                  ~1,
                                  "reduced")
    so.gene <- sleuth::sleuth_lrt(so.gene, 
                                  "reduced",
                                  "full")
    for (i in colnames(design)[grep("Intercept", colnames(design), invert = T)]){
      so.gene <- sleuth::sleuth_wt(so.gene, i)  
    }
    rt.gene.list <- lapply(colnames(design)[grep("Intercept", colnames(design), invert = T)], function(x){
      rt.gene <- sleuth::sleuth_results(so.gene, x)
      rt.gene <- rt.gene[order(rt.gene$qval),]
    })
    names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    # gene-level expression is summed from transcript level data (sum(TPM))
    sfInit(parallel = T, cpus = cpus)
    target_mapping <- so$target_mapping
    sfExport("target_mapping")
    sfExport("kt_wide")
    l1 <- sfLapply(unique(target_mapping$ens_gene), function(x) {
      s <- apply(kt_wide[target_mapping[target_mapping$ens_gene == x, "target_id"], ], 2, sum)
    })
    sfStop()
    names(l1) <- unique(target_mapping$ens_gene)
    kt_genes <- as.data.frame(do.call("rbind", l1))
    kt_genes$ensembl_gene_id <- rownames(kt_genes)
    kt_genes <- kt_genes[, c("ensembl_gene_id", s2c.list[[x]]$sample)]
    kt_genes <- tibble::as_tibble(merge(kt_genes, ensGenes, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id"))
    kt_wide <- tibble::as_tibble(tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm))
    return(list(sleuth_object = so,
                sleuth_results = rt.list,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                sleuth_results.gene = rt.gene.list,
                kallisto_table_genes = kt_genes))
  })
  save(results, file = "MDCK_kallisto_analysis_results.rda")
} else {
  load("MDCK_kallisto_analysis_results.rda")
}

# re-formatting of list object --------------------------------------------
names(results) <- names(s2c.list)
resultsCompressed <- lapply(names(results), function(x){
  results[[x]][grep("sleuth_object", names(results[[x]]), invert = T)]
})
names(resultsCompressed) <- names(results)

resultsCompressed <- lapply(names(resultsCompressed), function(x){
  resultsCompressed[[x]][grep("kallisto_pca", names(resultsCompressed[[x]]), invert = T)]
})
names(resultsCompressed) <- names(results)
save(resultsCompressed, file = "resultsCompressed.rda")

load("resultsCompressed.rda")
resultsCompressedBU <- resultsCompressed

resultsCompressed <- lapply(names(resultsCompressed), function(x) {
  resultsCompressed[[x]]$kallisto_table_wide <- resultsCompressed[[x]]$kallisto_table_wide[, c("target_id", s2c.list[[x]]$sample)]
  return(resultsCompressed[[x]])
})

names(resultsCompressed) <- names(s2c.list)

# get list of EMT genes (from qPCR array) ---------------------------------
if(!file.exists("cfam.qPCRGenesTab.rda")){
  qPCRGeneList <- readLines(lDir(pathPrefix, "Data/Tremethick/EMT/ChIP-Seq/MDCK qPCR data/genelist.txt"))
  mart <- useEnsembl(biomart = biomart, host = ensemblHost, dataset = dataset)
  human <- biomaRt::useEnsembl(biomart = biomart, dataset = "hsapiens_gene_ensembl", host = ensemblHost)
  hsap_ensembl_gene_ids <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
        filter = "hgnc_symbol",
        values = qPCRGeneList,
        mart = human)
  cfam_ensembl_gene_ids <- getBM(attributes = c("ensembl_gene_id", "cfamiliaris_homolog_ensembl_gene"), 
        filter = "ensembl_gene_id",
        values = hsap_ensembl_gene_ids,
        mart = human)
  hsap_ensembl_gene_ids <- merge(hsap_ensembl_gene_ids, 
                                 cfam_ensembl_gene_ids, 
                                 by.x = "ensembl_gene_id",
                                 by.y = "ensembl_gene_id",
                                 all.x = T)
  rownames(hsap_ensembl_gene_ids) <- hsap_ensembl_gene_ids$ensembl_gene_id
  qPCRGeneList.missing <- c("KRT7" = "ENSCAFG00000007307",
                            "LOC488207" = NULL,
                            "OCLN" = "ENSCAFG00000007805",
                            "SIP1" = "ENSCAFG00000013859",
                            "TCF4" = "ENSCAFG00000000140",
                            "TGFB1" = "ENSCAFG00000005014",
                            "TMEFF1" = "ENSCAFG00000002577",
                            "TWIST1" = "ENSCAFG00000012469", #using TWIST2 - TWIST does not seem to exist in dog genome
                            "LOC478215/H2AZ" = "ENSCAFG00000010615",
                            "HPRT1" = "ENSCAFG00000018870",
                            "LDHAL6B" = "ENSCAFG00000009211")
  qPCRGeneList.missing <- data.frame(hgnc_symbol = names(qPCRGeneList.missing),
                                     cfamiliaris_homolog_ensembl_gene = qPCRGeneList.missing)
  rownames(qPCRGeneList.missing) <- qPCRGeneList.missing$hgnc_symbol
  i1 <- intersect(qPCRGeneList.missing$hgnc_symbol, hsap_ensembl_gene_ids$hgnc_symbol)
  hsap_ensembl_gene_ids <- hsap_ensembl_gene_ids[! hsap_ensembl_gene_ids$hgnc_symbol %in% i1, ]
  tab1 <- rbind(hsap_ensembl_gene_ids[, c("hgnc_symbol", "cfamiliaris_homolog_ensembl_gene")],
                qPCRGeneList.missing[i1, c("hgnc_symbol", "cfamiliaris_homolog_ensembl_gene")],
                qPCRGeneList.missing["LOC478215/H2AZ", c("hgnc_symbol", "cfamiliaris_homolog_ensembl_gene")])
  
  cfam.qPCRGenesTab <- tab1[!tab1$cfamiliaris_homolog_ensembl_gene == "",]
  colnames(cfam.qPCRGenesTab)[2] <- "ensembl_gene_id"
  save(cfam.qPCRGenesTab, file = "cfam.qPCRGenesTab.rda")
} else {
  load("cfam.qPCRGenesTab.rda")
}

# make volcano plots for shZ & TGFb vs WT ---------------------------------
if(file.exists(lDir(devPath, "JCSMR-Tremethick-Lab/H2AZ_EMT/R_scripts/MDCK_volcano_plots_RNA-Seq.R"))){
  source(lDir(devPath, "JCSMR-Tremethick-Lab/H2AZ_EMT/R_scripts/MDCK_volcano_plots_RNA-Seq.R"))
  }

# heatmap of samples using MCF10A_wt as reference
df1 <- resultsCompressed[["MDCK"]]$kallisto_table_genes
dim(df1)
head(df1)
df1 <- df1[, c("ensembl_gene_id", s2c.list[["MDCK"]]$sample)]
df1 <- df1[complete.cases(df1),]
df2 <- log2(df1[, c(2:ncol(df1))] + 1)
rownames(df2) <- df1$ensembl_gene_id
filter <- apply(df2, 1, function(y) length(y[y>2])>=0.1)
df2 <- as.matrix(df2[filter, ])
sd1 <- apply(df2, 1, sd)
length(sd1)
table(sd1 > 1)
pca1 <- ade4::dudi.pca(t(df2[sd1 > 1, ]), scannf = F, nf = 5, scale = T, center = T)
pdf("Heatmaps_MDCK_WT_TGFb_shZ.pdf")
heatmap.3(df2[sd1 > 1.5,], 
          trace = "none",
          cexCol = 0.6,
          labRow = NA,
          hclustfun=function(x) hclust(x,method="ward.D"), 
          main = "Heatmap based on log2(sum(tpm) + 1)")
heatmap.3(t(pca1$tab),
          trace = "none",
          cexCol = 0.6,
          labRow = NA,
          hclustfun=function(x) hclust(x,method="ward.D"),
          main = "Heatmap based on scaled, centered data")
dev.off()
head(df2)
table(is.na(sd1))
head(sd1[is.na(sd1)])
head(df2[is.na(sd1),])
