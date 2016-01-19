calculateCoverage <- function(step = 1, gr1 = NULL, cov1 = NULL, func = "mean"){
  if (is.null(range)) stop("GRanges object missing")
  if (is.null(cov1)) stop("Coverage object missing")
  gr.tiles <- unlist(tile(gr1, width = step))
  seqlevels(gr.tiles, force = T) <- seqlevels(gr.tiles)[order(seqlevels(gr.tiles))]
  cov1 <- cov1[which(names(cov1) %in% seqlevels(gr.tiles))]
  cov1 <- cov1[names(cov1)[order(names(cov1))]]
  bA <- binnedAverage(gr.tiles, cov1, func)
  return(bA)
}