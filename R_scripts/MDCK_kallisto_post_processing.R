library(deepToolsUtils)
library(biomaRt)
library(tidyr)
library(RColorBrewer)
library(tibble)
library(snowfall)
library(sleuth)
library(tibble)
library(tximport)
library(readr)

# external functions ------------------------------------------------------
source("~/Development/GeneralPurpose/R/amILocal.R")
source("~/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Development/GeneralPurpose/R/lsos.R")
source("~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/R_scripts/prepare_ensembl_annotation.R")

# local functions ---------------------------------------------------------
lDir <- function(x, y, ...){
  paste(x, y, ..., sep = "/")
}

# get snakemake run configuration -----------------------------------------
runConfig <- jsonlite::fromJSON("~/Development/JCSMR-Tremethick-Lab/H2AZ_EMT/snakemake/configs/config.json")
runID <- names(runConfig$samples$`RNA-Seq`)[2]
refVersion <- "CanFam3.1"
annotationVersion <- runConfig$references[[refVersion]]$version
annotationVersion <- annotationVersion[1]

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

analysisDir <- lDir(pathPrefix, 
                    "Data/Tremethick/EMT/RNA-Seq/",
                    runID,
                    "/processed_data/CanFam3.1_ensembl84_ERCC/R_Analysis/")

if (dir.exists(analysisDir)){
  setwd(analysisDir)
} else {
  dir.create(analysisDir)
  setwd(analysisDir)
}

dataPath <- lDir(pathPrefix, 
                 "Data/Tremethick/EMT/RNA-Seq/",
                 runID,
                 "/processed_data/CanFam3.1_ensembl84_ERCC/HTSeq/count/")

devPath <- "~/Development"

files <- list.files(path = dataPath, full.names = T)
names(files) <- list.files(path = dataPath, full.names = F)


# preparing annotation data from Ensembl ------------------------------------
annotationDataPath <- analysisDir
prepare_ensembl_annotation(annotationDataPath, annotationVersion)

# load kallisto data with tximport and inspect via PCA -------------------------
base_dir <- paste(pathPrefix, 
                "Data/Tremethick/EMT/RNA-Seq",
                runID ,
                "processed_data",
                annotationVersion, 
                "kallisto", 
                sep = "/")
sample_id <- dir(base_dir)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
condition <- unlist(lapply(strsplit(sample_id, "_"), function(x) x[1]))
condition <- gsub("D6", "", condition)
files <- paste(kal_dirs, "abundance.h5", sep = "/")
names(files) <- sample_id
txi <- tximport::tximport(files, 
                          type = "kallisto",
                          geneIdCol = "ens_gene",
                          txIdCol = "target_id",
                          tx2gene = t2g)
save(txi, file = "kallisto_estimated_abundances.rda")

# exploratory analysis of the abundance data ------------------------------
# perform PCA for first inspection of data
pca1 <- ade4::dudi.pca(t(txi$abundance), scannf = F, nf = 6)
pdf("Exploratory_Analysis_sclass_plot.pdf")
ade4::s.class(pca1$li, fac = as.factor(condition))
dev.off()
pdf("Exploratory_Analysis_sarrow_plot.pdf")
ade4::s.arrow(pca1$li, clabel = 0.7)
dev.off()

pairs(log2(txi$abundance + 1))
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

condition <- as.factor(condition)
s2c <- data.frame(sample = sample_id, condition = condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$sample <- as.character(s2c$sample)
s2c.list <- list(MDCK = s2c)

sleuth_results_files <- "MDCK_kallisto_analysis_results_V2.rda"
# actual processing using sleuth------------------------------------------------
if (!file.exists(sleuth_results_files)){
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
      rt.gene <- sleuth::sleuth_results(so.gene, x, show_all = F)
      rt.gene <- data.table::data.table(rt.gene[order(rt.gene$qval),])
    })
    names(rt.gene.list) <- colnames(design)[grep("Intercept", colnames(design), invert = T)]
    # gene-level expression is summed from transcript level data (sum(TPM))
    kt_genes <- data.table::data.table(sleuth::kallisto_table(so.gene))
    kt_genes_wide <- tidyr::spread(kt_genes[, c("target_id", "sample", "tpm")], sample, tpm)
    kt_wide <- data.table::data.table(tidyr::spread(kt[, c("target_id", "sample", "tpm")], sample, tpm))
    return(list(sleuth_object = so,
                sleuth_object_genes = so.gene,
                sleuth_results = rt.list,
                kallisto_table = kt,
                kallisto_table_wide = kt_wide,
                sleuth_results.gene = rt.gene.list,
                kallisto_table_genes = kt_genes))
  })
  names(results) <- names(s2c.list)
  save(results, file = sleuth_results_files)
} else {
  load(sleuth_results_files)
}
names(results) <- names(s2c.list)

# re-formatting of list object --------------------------------------------
if (file.exists("resultsCompressed.rda")){
  load("resultsCompressed.rda")
} else {
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
}
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
df1 <- sleuth::sleuth_to_matrix(results[["MDCK"]]$sleuth_object, "obs_norm", "tpm")
df1 <- df1$data
tgfb <- grep("TGFb", colnames(df1))
shz <- grep("shZ", colnames(df1))
wt <- grep("TGFb|shZ", colnames(df1), invert = T)
df1 <- df1[, c(wt,shz,tgfb)]
dim(df1)
head(df1)
df2 <- log2(df1 + 1)
filter <- apply(df2, 1, function(y) length(y[y>2])>=0.1)
df2 <- as.matrix(df2[filter, ])
sd1 <- apply(df2, 1, sd)
length(sd1)
table(sd1 > 0.5)
pca1 <- ade4::dudi.pca(t(df2[sd1 > 0.5, ]), scannf = F, nf = 5, scale = T, center = T)
ade4::s.class(pca1$li, fac = as.factor(condition))
pdf("Heatmaps_MDCK_WT_TGFb_shZ.pdf")
heatmap.3(df2[sd1 > 0.5,], 
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
