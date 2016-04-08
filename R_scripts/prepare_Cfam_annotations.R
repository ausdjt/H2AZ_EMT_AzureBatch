# prepare Cfam annotations for deepTools plotting

# all genes
df <- data.frame(seqnames = seqnames(Cfam3.genes),
                 starts = start(Cfam3.genes) -1,
                 ends = end(Cfam3.genes),
                 names = Cfam3.genes$hgnc_symbol,
                 scores = ".",
                 strand = strand(Cfam3.genes))
write.table(df, file = "~/Data/Tremethick/EMT/GenomeWide/deepTools/regionFiles/Cfam3.genes.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# qPCR genes
df <- data.frame(seqnames = seqnames(gr.qPCRGenesPositions),
                 starts = start(gr.qPCRGenesPositions) -1,
                 ends = end(gr.qPCRGenesPositions),
                 names = gr.qPCRGenesPositions$hgnc_symbol,
                 scores = ".",
                 strand = strand(gr.qPCRGenesPositions))
write.table(df, file = "~/Data/Tremethick/EMT/GenomeWide/deepTools/regionFiles/Cfam3.qPCR.genes.bed", quote = F, sep = "\t", row.names = F, col.names = F)



