#' Convert gene identifiers
#'
#' Convert a vector of gene identifiers from one standard to another. The input type of the identifiers
#' is specified using \code{input_type}, and the output type is specified using \code{output_type}.
#' By default, this function uses the \code{annotables} package. If annotables
#' is not installed, or if the user specifies not to use it, \code{biomaRt} is used instead. Genes can be
#' filtered by biotype; any genes not matching the specified biotype(s) will return NA.
#' Currently, support is only provided for human gene identifiers.
#' @param genes character vector containing gene identifiers to convert
#' @param input_type the gene identifier type to convert from. Must match a corresponding identifier type in \code{annotables} or \code{biomaRt}.
#' @param output_type the gene identifier type to convert to Must match a corresponding identifier type in \code{annotables} or \code{biomaRt}.
#' @param biotype (optional) character vector containing biotypes of genes to retain. If provided, genes with biotype not matching elements of this vector return NA. If NULL, all genes are retained.
#' @param use_annotables boolean, whether to use the annotables package. If annotables is not installed, the function defaults to using \code{biomaRt}.
#' @export
#' @return a vector of gene symbols, of length equivalent to the length of genes.
#' @usage \code{convert_gene_names(genes, input_type, output_type, biotype=NULL, use_annotables=TRUE)}
convert_gene_names <- function(genes, input_type, output_type, biotype=NULL, use_annotables=TRUE) {
  if (require(annotables) & use_annotables) {
    cat('Package "annotables" detected. Using data from annotables instead of BioMart.\n')
    grch38.tmp <- if (is.null(biotype)) {grch38} else {grch38[grch38[["biotype"]] %in% biotype,]}
    if (!(output_type %in% colnames(grch38.tmp))) stop("Specified 'output_type' is not available.")
    genes_out <- grch38.tmp[[output_type]][match(genes, grch38.tmp[[input_type]])] # get the gene names
  } else {
    cat('Using BioMart to convert gene symbols. This could take a while.\n')
    if (input_type=="ensgene") input_type <- "ensembl_gene_id"
    ensembl <- biomaRt::useMart('ensembl', 'hsapiens_gene_ensembl')
    gene_ID_translate <-
      biomaRt::getBM(attributes=c(input_type, output_type, "gene_biotype"),
                     filters=input_type, values=list(genes),
                     mart=ensembl)
    if (!is.null(biotype))
      gene_ID_translate <- gene_ID_translate[gene_ID_translate[,"gene_biotype"] %in% biotype,]
    genes_out <- gene_ID_translate[match(genes, gene_ID_translate[,input_type]),outype_type] # get the gene names
  }
  return(genes_out)
}
