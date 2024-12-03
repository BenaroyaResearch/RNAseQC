#' Estimate saturation of genes based on rarefaction of reads
#'
#' Estimate the saturation of gene detection based on rarefaction of the mapped read counts from each
#' library in a read counts object. This function takes the read counts for each library and
#' sequentially rarefies them at different levels to determine how thoroughly genes are being sampled.
#' Optional settings include the minimum number of counts for a gene to be counted as "detected"
#' (default=1), and if using the sampling method, the number of intermediate points to sample (default=6)
#' and the number of times to sample at each depth (default=5).
#' @param counts a numeric matrix (or object that can be coerced to a matrix) containing read counts, or an object from which counts can be extracted. Should have genes in rows and samples in columns.
#' @param max_reads the maximum number of reads to sample at. By default, this value is the maximum of total read counts across all libraries.
#' @param method character, either "division" or "sampling". Method "sampling" is slower but more realistic, and yields smoother curves. Method "division" is faster but more coarse and less realistic. See Details for more complete description
#' @param ndepths the number of depths to sample at. 0 is always included.
#' @param nreps the number of sampling iterations to take for each library at each depth. With well-sampled libraries, 1 should be sufficient. With poorly-sampled libraries, sampling variance may be substantial, requiring higher values.
#' @param min_counts the minimum number of counts for a gene to be counted as detected. Genes with sample counts >= this value are considered detected. Defaults to 1. Set to NULL to use min_cpm.
#' @param min_cpm the minimum counts per million for a gene to be counted as detected. Only relevant with \code{method = "sampling"}. Genes with sample cpm >= this value are considered detected. Either this or min_count should be specified, but not both; including both yields an error, as does specifying min_cpm with \code{method = "division"}. Defaults to NULL.
#' @param verbose logical, whether to output the status of the estimation.
#' @import countSubsetNorm
#' @importFrom bigtabulate bigtable
#' @export
#' @details The \code{method} parameter determines the approach used to estimate the number of genes detected at different sequencing depths. Method "division" simply divides the counts for each gene by a series of scaling factors, then counts the genes whose adjusted counts exceed the detection threshold. Method "sampling" generates a number of sets (\code{nreps}) of simulated counts for each library at each sequencing depth, by probabilistically simulating counts using observed proportions. It then counts the number of genes that meet the detection threshold in each simulation, and takes the arithmetic mean of the values for each library at each depth.
#' @return A data frame containing \code{ndepths x nSamples} rows, with one row for each sample at each depth. Columns include "sample" (the name of the sample identifier), "depth" (the depth value for that iteration), and "sat" (the number of genes detected at that depth for that sample).  For method "sampling", it includes an additional column with the variance of genes detected across all replicates of each sample at each depth.
estimate_saturation <-
  function(counts, max_reads = Inf,
           method = "sampling",
           ndepths = 6, nreps = 5,
           min_counts = 1, min_cpm = NULL,
           verbose = FALSE) {
    
    # extract counts and/or convert to matrix
    counts <-
      countSubsetNorm::extract_counts(counts, return_class = "matrix")
    
    # check inputs
    checkmate::assert_array(counts)
    checkmate::assert_numeric(counts)
    if (nrow(counts) < ncol(counts))
      warning("The input counts object has more columns than rows. ",
              "Are you sure the rows are genes and the columns are samples?")
    
    checkmate::assert_numeric(max_reads, lower = 0)
    checkmate::assert_character(method)
    method <- match.arg(method, choices = c("division", "sampling"))
    
    checkmate::assert_integerish(ndepths, lower = 1)
      
    checkmate::assert_numeric(min_counts, null.ok = TRUE)
    checkmate::assert_numeric(min_cpm, null.ok = TRUE)
    if (sum(!is.null(min_counts), !is.null(min_cpm)) != 1)
      stop("One of min_counts or min_cpm must be specified, but not both.")

    if (!is.null(min_cpm) && method == "division")
      stop("Thresholding by min_cpm is only relevant with method = 'sampling'.")
    
    checkmate::assert_logical(verbose)

    # calculate characteristics of libraries and calculations
    readsums <- colSums(counts)
    max_reads <- min(max(readsums), max_reads)
    depths <- round(seq(from = 0, to = max_reads, length.out = ndepths+1))
    
    # generate empty data frame to store results
    saturation <-
      data.frame(sample = as.vector(sapply(colnames(counts), FUN = rep, times = ndepths+1)),
                 depth = rep(depths, time = ncol(counts)))
    sat_estimates <- as.numeric(rep(NA, ncol(counts) * length(depths)))
    
    # create an iterator for rows of the data frame
    counter <- 0
    if (method == "sampling") {
      checkmate::assert_integerish(nreps, lower = 1)
      sat_var_estimates <- as.numeric(rep(NA, ncol(counts) * length(depths)))
    }
    
    # iterate over each column in counts, and each depth
    for (lib_current in seq_len(ncol(counts))) {
      if (verbose) cat("Working on library", lib_current, "of", ncol(counts), "\n")

      # calculate gene probabilities for the focal library
      probs <- counts[, lib_current, drop = TRUE] / readsums[lib_current]
      probs <- probs[probs > 0] # zero counts add nothing but computational time!
      ngenes <- length(probs)
      
      # calculate for each depth in the focal library
      for (depth_current in depths) {
        counter <- counter + 1
        if (depth_current == 0) {
          sat_estimates[counter] <- 0
          if (method == "sampling")
            sat_var_estimates[counter] <- 0
        } else if (depth_current > readsums[lib_current]) {
          sat_estimates[counter] <- NA
          if (method == "sampling")
            sat_var_estimates[counter] <- NA
        } else if (method == "division") {
          sat_estimates[counter] <- sum((probs * depth_current) >= min_counts)
        } else if (method == "sampling") {
          sat_estimate_lib_depth_current <- as.numeric(rep(NA, nreps))
          
          if (!is.null(min_cpm)) {
            min_counts_lib_current <- depth_current * min_cpm / 1E6
          } else min_counts_lib_current <- min_counts
          
          for (rep_current in seq_len(nreps)) {
            reads <- as.matrix(sample.int(n = ngenes, size = depth_current, replace = TRUE, prob = probs))
            sat_estimate_lib_depth_current[rep_current] <-
              sum(bigtabulate::bigtable(reads, ccol = 1) >= min_counts_lib_current)
          }
          sat_estimates[counter] <- mean(sat_estimate_lib_depth_current)
          sat_var_estimates[counter] <- var(sat_estimate_lib_depth_current)
        }
      }
    }
    saturation$sat <- sat_estimates
    if (method == "sampling")
      saturation$sat.var <- sat_var_estimates
    
    return(saturation)
  }
