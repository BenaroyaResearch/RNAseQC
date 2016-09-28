#' Convert gene identifiers to HGNC gene symbols
#'
#' Convert a vector of gene identifiers to HGNC gene symbols. The type of the input identifiers is
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
    cat('Package "annotables" detected. Using data from annotables instead of BioMart.\n')
    grch38.tmp <- if (is.null(type)) {grch38} else {grch38[grch38[["biotype"]] %in% type,]}
    genes_HGNC <- grch38.tmp[["symbol"]][match(genes, grch38.tmp[[from]])] # get the HGNC gene names
  } else {
    cat('Using BioMart to retrieve HGNC gene symbols. This could take a while.\n')
    if (from=="ensgene") from <- "ensembl_gene_id"
    ensembl <- biomaRt::useMart('ensembl', 'hsapiens_gene_ensembl')
    gene_ID_to_HGNC <-
      biomaRt::getBM(attributes=c(from, "hgnc_symbol","gene_biotype"),
                     filters=from, values=list(genes),
                     mart=ensembl)
    if (!is.null(type))
      gene_ID_to_HGNC <- gene_ID_to_HGNC[gene_ID_to_HGNC[,"gene_biotype"] %in% type,]
    genes_HGNC <- gene_ID_to_HGNC[match(genes, gene_ID_to_HGNC[,from]),"hgnc_symbol"] # get the HGNC gene names
  }
  return(genes_HGNC)
}