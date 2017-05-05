#' Estimate saturation of genes based on rarefaction of reads
#' 
#' Estimate the saturation of gene detection based on rarefaction of the mapped read counts from each
#' library in a read counts object. This function takes the read counts for each library and
#' sequentially rarefies them at different levels to determine how thoroughly genes are being sampled.
#' Optional settings include the number of intermediate points to sample (default=6), the number of
#' times to sample at each depth (default=5), and the minimum number of counts for a gene to be
#' counted as "detected" (default=1).
#' @param counts a numeric matrix (or object that can be coerced to a matrix) containing read counts, with genes in rows and samples in columns.
#' @param max_reads the maximum number of reads to sample at. By default, this value is the maximum of total read counts across all libraries.
#' @param method character, either "division" or "sampling". Method "sampling" is slower but more realistic, and yields smoother curves. Method "division" is faster but more coarse and less realistic. See Details for more complete description
#' @param ndepths the number of depths to sample at. 0 is always included.
#' @param nreps the number of samples to take for each library at each depth. With well-sampled libraries, 1 should be sufficient. With poorly-sampled libraries, sampling variance may be substantial, requiring higher values.
#' @param min_counts the minimum number of counts for a gene to be counted as detected. Genes with sample counts >= this value are considered detected. Defaults to 1. Set to NULL to use min_cpm.
#' @param min_cpm the minimum counts per million for a gene to be counted as detected. Genes with sample counts >= this value are considered detected. Either this or min_count should be specified, but not both; including both yields an error. Defaults to NULL.
#' @param verbose logical, whether to output the status of the estimation.
#' @export
#' @details The \code{method} parameter determines the approach used to estimate the number of genes detected at different sequencing depths. Method "division" simply divides the counts for each gene by a series of scaling factors, then counts the genes whose adjusted counts exceed the detection threshold. Method "sampling" generates a number of sets (\code{nreps}) of simulated counts for each library at each sequencing depth, by probabilistically simulating counts using observed proportions. It then counts the number of genes that meet the detection threshold in each simulation, and takes the arithmetic mean of the values for each library at each depth.
#' @return A data frame containing \code{nrep * ndepths} rows, with one row for each sample at each depth. Columns include "sample" (the name of the sample identifier), "depth" (the depth value for that iteration), and "sat" (the number of genes detected at that depth for that sample).  For method "sampling", it includes an additional column with the variance of genes detected across all replicates of each sample at each depth.
#' @usage \code{
#' estimate_saturation(
#'   counts, max_reads=Inf,
#'   method="sampling",
#'   ndepths=6, nreps=5,
#'   min_counts=1, min_cpm=NULL,
#'   verbose=FALSE)}
estimate_saturation <-
  function(counts, max_reads=Inf,
           method="sampling",
           ndepths=6, nreps=5,
           min_counts=1, min_cpm=NULL,
           verbose=FALSE) {
    if (sum(!is.null(min_counts), !is.null(min_cpm)) != 1)
      stop("Either min_counts or min_cpm must be specified, but not")
    method <- match.arg(method, choices=c("division", "sampling"))
    if (!is.matrix(counts))
      counts <- as.matrix(counts)
    readsums <- colSums(counts)
    max_reads <- min(max(readsums), max_reads)
    depths <- round(seq(from=0, to=max_reads, length.out=ndepths+1))
    saturation <- data.frame(sample=as.vector(sapply(colnames(counts), FUN=rep, times=ndepths+1)),
                             depth=rep(depths, time=ncol(counts)))
    sat.estimates <- numeric()
    if (method=="sampling") sat.var.estimates <- numeric()
    for (i in 1:ncol(counts)) {
      if (verbose) cat("Working on library", i, "of", ncol(counts), "\n")
      
      # adjust to min_cpm if specified
      if (!is.null(min_cpm)) {
        min_counts.lib <- readsums[i] / 1000000
      } else {
        min_counts.lib <- min_counts
      }
      
      probs <- counts[,i] / readsums[i] # calculate gene probabilities for the library
      probs <- probs[probs > 0] # zero counts add nothing but computational time!
      ngenes <- length(probs)
      for (j in depths) {
        if (j > readsums[i]) {
          sat.estimates <- c(sat.estimates, NA)
          if (exists("sat.var.estimates"))
            sat.var.estimates <- c(sat.var.estimates, NA)
        } else if (method=="division") {
          sat.estimates <- c(sat.estimates,
                             sum((probs * j) >= min_counts.lib))
        } else if (method=="sampling") {
          est <- as.numeric(rep(NA, nreps))
          for (k in 1:nreps) {
            reads <- sample.int (n=ngenes, size=j, replace=TRUE, prob=probs)
            est[k] <- sum(table(reads) >= min_counts.lib)
          }
          sat.estimates <- c(sat.estimates, mean(est))
          sat.var.estimates <- c(sat.var.estimates, var(est))
        }
      }
    }
    saturation$sat <- sat.estimates
    if (exists("sat.var.estimates"))
      saturation$sat.var <- sat.var.estimates
    
    return(saturation)
  }
