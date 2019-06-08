#' Run PCA on gene expression count data
#'
#' This is a wrapper function for prcomp, designed for running PCA on RNAseq data. It transposes the input data, and adds two elements to the prcomp object, containing the following for each PC: the percent variance, and a character string for labeling plot axes.
#' @param counts matrix or data frame of gene expression counts, or an object from which counts can be extracted. Should have genes in rows and samples in columns
#' @param cpm logical, whether to transform the counts to counts-per-million (using \code{edgeR::cpm}) before running prcomp.
#' @param log2_transform logical, whether to transform the counts (using \code{log2(counts+1)}) before running prcomp.
#' @param ... (optional) additional arguments passed to prcomp, such as \code{center} and \code{scale}.
#' @export
#' @return a prcomp object, with additional list elements containing the percent variance and strings for labeling plot axes.
calc_PCAs <- function(counts, cpm=TRUE, log2_transform=FALSE,
                      ...) {
  counts <- extract_counts(counts)
  if (cpm) counts <- edgeR::cpm(counts)
  t.counts <- t(counts)
  if (log2_transform) t.counts <- log2(t.counts + 1)
  pcaAll <- prcomp(t.counts, ...)
  pcaAll$pvars <- (pcaAll$sdev^2)/sum(pcaAll$sdev^2)
  pcaAll$pvars.labs <- paste0("PC", 1:length(pcaAll$pvars), " (", round(pcaAll$pvars*100, 1), "%)")
  return(pcaAll)
}
