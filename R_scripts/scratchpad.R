down.genes <- qPCRGenesTab[rownames(qDE.limma.tab[which(-log10(qDE.limma.tab$adj.P.Val) >= 1 & qDE.limma.tab$logFC < 0), ]), ]
up.genes <- qPCRGenesTab[rownames(qDE.limma.tab[which(-log10(qDE.limma.tab$adj.P.Val) >= 1 & qDE.limma.tab$logFC  > 0), ]), ]

down.genes.tss <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "transcription_start_site"), filters = "ensembl_gene_id", values = down.genes$ensembl_gene_id, dog)
up.genes.tss <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "transcription_start_site"), filters = "ensembl_gene_id", values = up.genes$ensembl_gene_id, dog)

gr.down.genes.tss <- GRanges(down.genes.tss$chromosome_name, IRanges(down.genes.tss$transcription_start_site, down.genes.tss$transcription_start_site + 1), strand = "*", down.genes.tss$hgnc_symbol, down.genes.tss$ensembl_gene_id)
gr.up.genes.tss <- GRanges(up.genes.tss$chromosome_name, IRanges(up.genes.tss$transcription_start_site, up.genes.tss$transcription_start_site + 1), strand = "*", up.genes.tss$hgnc_symbol, up.genes.tss$ensembl_gene_id)

export(gr.down.genes.tss, format = "BED", "~/Data/Annotations/CanFam3/EMT_down_genes_TSS.bed")
export(gr.up.genes.tss, format = "BED", "~/Data/Annotations/CanFam3/EMT_up_genes_TSS.bed")



TxDb.CanFam3.1 <- makeTxDbFromGFF(gzfile('~/Data/Annotations/CanFam3/Ensembl/Canis_familiaris.CanFam3.1.83.chr.gtf.gz'), dataSource = "Ensembl", organism = "Canis familiaris", circ_seqs=DEFAULT_CIRC_SEQS)
gr.genes <- genes(TxDb.CanFam3.1)

all.genes.tss <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "transcription_start_site"), 
                       filters = "ensembl_gene_id", 
                       values = gr.genes$gene_id,
                       dog)
gr.all.genes.tss <- GRanges(all.genes.tss$chromosome_name, IRanges(all.genes.tss$transcription_start_site, all.genes.tss$transcription_start_site + 1), strand = "*", all.genes.tss$hgnc_symbol, all.genes.tss$ensembl_gene_id)
gr.all.genes.tss <- reduce(gr.all.genes.tss)

export(gr.all.genes.tss, format = "BED", "~/Data/Annotations/CanFam3/Ensembl/all_genes_TSS.bed")
