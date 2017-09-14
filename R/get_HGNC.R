#' Convert gene identifiers to HGNC gene symbols
#'
#' This is a simple wrapper for the more flexible function \code{convert_gene_names}. I created it to
#' facilitate compatability with old code. It simply calls \code{convert_gene_names} with the hard-coded
#' specification to return HGNC symbols, and some variable name translation.
#' It converts a vector of gene identifiers to HGNC gene symbols. The type of the input identifiers is
#' specified using \code{from}. By default, this function uses the \code{annotables} package. If annotables
#' is not installed, or if the user specifies not to use it, \code{biomaRt} is used instead. Genes can be
#' filtered by type; any genes not matching the specified type(s) will return NA.
#' @param genes character vector containing gene identifiers to convert
#' @param from the gene identifier class to convert from. Must match a corresponding variable in \code{annotables} or \code{biomaRt}.
#' @param type (optional) character vector containing types of genes to retain. If provided, genes with biotype not matching elements of this vector return NA. If NULL, all genes are retained.
#' @param use_annotables boolean, whether to use the annotables package. If annotables is not installed, the function defaults to using \code{biomaRt}.
#' @export
#' @return a vector of HGNC gene symbols, of length equivalent to the length of genes.
#' @usage \code{get_HGNC(genes, from="ensgene", type=NULL, use_annotables=TRUE)}
get_HGNC <- function(genes, from="ensgene", type=NULL, use_annotables=TRUE) {
  if (require(annotables) & use_annotables) {
    output_type <- "symbol"
  } else output_type <- "hgnc_symbol"
  genes_HGNC <- convert_gene_names(
    genes, input_type=from, output_type=output_type, biotype=type, use_annotables=use_annotables)
  return(genes_HGNC)
}
