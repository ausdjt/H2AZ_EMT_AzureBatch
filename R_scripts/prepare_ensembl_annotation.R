# function prepare_ensembl_annotation()

prepare_ensembl_annotation <- function(annotationDataPath = NULL, annotationVersion = NULL, ...){
  stopifnot(!is.null(annotationDataPath), !is.null(annotationVersion))
  adp = annotationDataPath
  av = annotationVersion
  ensGenes_file <- paste(adp, "ensGenes_", av, ".rda", sep = "")
  ensTranscripts_file <- paste(adp, "ensTranscripts_", av, ".rda", sep = "")
  t2g_file <- paste(adp, "t2g_", av, ".rda", sep = "")
  myLength_file <- paste(adp, "mylength_", av, ".rda", sep = "")
  myGC_file <- paste(adp, "myGC_", av, ".rda", sep = "")
  myBiotypes_file <- paste(adp, "myBiotypes_", av, ".rda", sep = "")
  myChroms_file <- paste(adp, "myChroms_", av, ".rda", sep = "")
  annotationFileList <- list(ensGenes_file, ensTranscripts_file, t2g_file, myLength_file, myGC_file, myBiotypes_file, myChroms_file) 
  
  if (!all(sapply(annotationFileList, file.exists))){
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
    save(ensGenes, file = ensGenes_file)
    
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
    save(ensTranscripts, file = ensTranscripts_file)
    
    # create t2g object
    t2g <- ensTranscripts[, c("ensembl_transcript_id", 
                              "ensembl_gene_id", 
                              "external_gene_name", 
                              "version", 
                              "transcript_version")]
    t2g$ensembl_transcript_id <- paste(t2g$ensembl_transcript_id, t2g$transcript_version, sep = ".")
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    save(t2g, file = t2g_file)
    
    mylength <- sapply(ensGenes$ensembl_gene_id, function(x){
      y <- ensTranscripts[which(ensTranscripts$ensembl_gene_id == x), ]
      y <- y[which.max(y$transcript_length), ]$transcript_length})
    save(mylength, file = myLength_file)
    mygc <- ensGenes$percentage_gc_content
    names(mygc) <- ensGenes$ensembl_gene_id
    save(mygc, file = myGC_file)
    mybiotypes <- ensGenes$gene_biotype
    names(mybiotypes) <- ensGenes$ensembl_gene_id
    save(mybiotypes, file = myBiotypes_file)
    mychroms <- data.frame(Chr = ensGenes$chromosome_name, GeneStart = ensGenes$start_position, GeneEnd = ensGenes$end_position)
    save(mychroms, file = myChroms_file)
  } else {
    load(ensGenes_file, .GlobalEnv)
    load(ensTranscripts_file, .GlobalEnv)
    load(myLength_file, .GlobalEnv)
    load(myGC_file, .GlobalEnv)
    load(myBiotypes_file, .GlobalEnv)
    load(myChroms_file, .GlobalEnv)
    load(t2g_file, .GlobalEnv)
  }
}
