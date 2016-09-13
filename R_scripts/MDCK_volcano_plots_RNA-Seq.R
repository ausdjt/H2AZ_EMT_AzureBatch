# volcano plot - TGFb vs WT

datDE <- lapply(names(resultsCompressed), function(x){
                      data <- lapply(names(resultsCompressed[[x]]$sleuth_results.gene), function(y){
                        print(x)
                        print(y)
                        dat <- resultsCompressed[[x]]$sleuth_results.gene[[y]][which(resultsCompressed[[x]]$sleuth_results.gene[[y]]$target_id %in% cfam.qPCRGenesTab$ensembl_gene_id),]
                        dat <- merge(dat,
                                     ensGenes[,c("ensembl_gene_id", "external_gene_name")], 
                                     by.x = "target_id", 
                                     by.y = "ensembl_gene_id")
                        dat <- dat[order(dat$qval), ]
                        yAxisMax <- max(-log10(dat$qval))
                        xAxisMax <- max(dat$b)
                        pdf(paste("Volcano_plot_", x, y, ".pdf", sep = ""))
                        plot(dat$b,
                             -log10(dat$qval), 
                             axes = F, 
                             xlab = "", 
                             ylab = "", 
                             frame = F,
                             xlim = c(-10,10),
                             cex = 0.3,
                             pch = 16, main = paste("Volcano plot\n", x, " Condition: ", y, sep = ""))
                        points(dat[which(-log10(dat$qval) >= 1), "b"], 
                               -log10(dat[which(-log10(dat$qval) >= 1), "qval"]),
                               col = "red", 
                               pch = 16, 
                               cex = 1.2)
                        axis(2, 
                             pos = 0, 
                             lwd = 3, 
                             at = c(seq(0,yAxisMax,10)))
                        axis(1,
                             pos = 0, 
                             lwd = 3)
                        mtext("-log10(q-value)", 
                              side = 2)
                        mtext("beta value",
                              side = 1, 
                              line = 2)
                        abline(h = 1, col = "red", lty = 2, lwd = 2)
                        # only adding labels to genes with adjuste p-value <= 0.01
                        text(dat[which(-log10(dat$qval) >= 2), "b"], -log10(dat[which(-log10(dat$qval) >= 2), "qval"]), 
                             labels = dat[which(-log10(dat$qval) >= 2), "external_gene_name" ],
                             cex = 0.7,
                             pos = 4, offset = 0.3)
                        text(-4,1.2, "< 0.1 [adjusted p-value]", cex = 0.7)
                        dev.off()
                        return(dat)
                      })
    names(data) <- names(resultsCompressed[[x]]$sleuth_results.gene)
    return(data)
})
names(datDE) <- names(resultsCompressed)
