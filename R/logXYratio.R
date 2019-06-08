#' Calculate ratio of total counts of genes mapping to the X and Y chromosome
#'
#' This function calculates the number of reads that map to the X and Y
#' chromosomes, and returns the (natural) log of the ratio of total X reads
#' to total Y reads. It can be used to verify or infer sex using RNAseq data.
#' Large values are likely to be from female samples; small values from male samples.
#' @param counts a matrix or data frame containing the gene expression counts. Should have samples in columns and genes in rows. Row names must contain gene names, in the class standard matching \code{gene_ID}.
#' @param lib_cols a numeric vector, the indices of columns containing count data. Defaults to all columns in \code{counts}; can be adjusted if non-count columns are included in the \code{counts} object.
#' @param gene_ID the gene identifier class of the gene names. Must match a corresponding variable in \code{annotables} or \code{biomaRt}.
#' @param use_annotables boolean, whether to use the annotables package. If annotables is not installed, defaults to using \code{biomaRt}.
#' @export
#' @return a vector of ratios, with one element for each sample.
#' @details \code{counts} should be normalized or raw, but not log-transformed.
#'  It assumes that each column in the counts object corresponds to a library.
#'  If the counts object contains additional columns, the columns containing
#'  libraries must be indicated in \code{lib_cols}.
logXYratio <-
  function(counts, lib_cols=1:ncol(counts),
           gene_ID="symbol", use_annotables=TRUE) {
    if (!is.data.frame(counts)) {
      stop("Counts object must be a data frame.")
    }
    force(lib_cols)  # causes lib_cols to be evaluated, which is necessary for later use
    
    if (require(annotables) & use_annotables) {
      warning('Package "annotables" detected. Used data from annotables instead of BioMart.\n')
      if (gene_ID=="ensembl") gene_ID <- "ensgene"
      gene_ID_to_chrom_name <-
        data.frame(gene_ID = rownames(counts),
                   chromosome_name = grch38[["chr"]][match(rownames(counts), grch38[[gene_ID]])])
      
    } else if (require(biomaRt)) {
      warning('Used BioMart to retrieve chromosome identity for all gene IDs.\n')
      if (gene_ID=="ensembl") gene_ID <- "ensembl_gene_id"
      ensembl <- useMart('ensembl', 'hsapiens_gene_ensembl')
      gene_ID_to_chrom_name <-
        biomaRt::getBM(
          attributes=c(gene_ID, "chromosome_name"), filters=gene_ID,
          values=rownames(counts), mart=ensembl)
      
    } else {
      stop("Sorry, I'm missing some needed packages.\nPlease install either the annotables or biomaRt package.\n")
    }
    
    counts$chromosome_name <-
      gene_ID_to_chrom_name$chromosome_name[match(rownames(counts), gene_ID_to_chrom_name[,1])]
    x_counts <- colSums(counts[counts$chromosome_name=="X", lib_cols], na.rm=TRUE)
    y_counts <- colSums(counts[counts$chromosome_name=="Y", lib_cols], na.rm=TRUE)
    ratios <- log((x_counts+1) / (y_counts+1))  # calculate log-transformed ratios of X reads to Y reads; add 1 to each because Y counts can be 0, which yields infinite ratio
    names(ratios) <- colnames(counts)[lib_cols]  # name the vector of ratios with the library IDs
    
    return(ratios)
  }
