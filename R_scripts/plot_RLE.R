# plot_RLE.R
require("reshape2")

calculate_relative_log_expression <- function(mat = NULL, ...){
  stopifnot(!is.null(mat))
  med <- apply(mat, 1, median)
  rle <- log2(mat / med)
  return(rle)
}

# first for all transcripts
m1 <- calculate_relative_log_expression(df1)
m1 <- melt(m1)
p1 <- ggplot(m1, aes(x=Var2, y = value, fill = Var2))
p1 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# second for ERCCs only
erccs <- grep("ERCC", rownames(df1))
m2 <- calculate_relative_log_expression(df1[erccs,])
dim(m2)
m2 <- melt(m2)
p2 <- ggplot(m2, aes(x=Var2, y = value, fill = Var2))
p2 + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


