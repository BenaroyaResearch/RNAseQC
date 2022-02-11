#' Calculate ratio of total counts of genes mapping to the X and Y chromosome
#'
#' This function calculates the number of reads that map to the X and Y
#' chromosomes, and returns the (natural) log of the ratio of total X reads
#' to total Y reads. It can be used to verify or infer sex using RNAseq data.
#' Large values are likely to be from female samples; small values from male samples.
#' @param counts a matrix or data frame containing the gene expression counts. Should have samples in columns and genes in rows. Row names must contain gene names, in the class standard matching \code{gene_ID}.
#' @param lib_cols a numeric vector, the indices of columns containing count data. Defaults to all columns in \code{counts}; can be adjusted if non-count columns are included in the \code{counts} object.
#' @param gene_ID the gene identifier class of the gene names. Must match a corresponding variable in \code{annotables} or \code{biomaRt}.
#' @param species character, the species that the data derive from. Can be "human" or "mouse", or any species abbreviation used in the BioMart ensembl datasets (as long as the species has X and Y chromosomes). For a full list of possible species, use \code{biomaRt::listDatasets(biomaRt::useMart("ensembl"))$dataset)}. For backward compatibility, defaults to human.
#' @param use_annotables boolean, whether to use the annotables package. If annotables is not installed, defaults to using \code{biomaRt}.
#' @export
#' @return a vector of ratios, with one element for each sample.
#' @details \code{counts} should be normalized or raw, but not log-transformed.
#'  It assumes that each column in the counts object corresponds to a library.
#'  If the counts object contains additional columns, the columns containing
#'  libraries must be indicated in \code{lib_cols}.
logXYratio <-
  function(counts, lib_cols=1:ncol(counts),
           gene_ID="symbol", species="human",
           use_annotables=TRUE) {
    if (!(is.data.frame(counts) | is.matrix(counts) | is(counts, "dgCMatrix"))) {
      stop("Counts object must be a data frame, matrix, or sparse matrix of class 'dgCMatrix'.")
    }
    force(lib_cols)  # causes lib_cols to be evaluated, which is necessary for later use
    
    if (require(annotables) & use_annotables) {
      warning('Package "annotables" detected. Used data from annotables instead of BioMart.\n')
      if (gene_ID=="ensembl") gene_ID <- "ensgene"
      if (species=="human") {
        geneData <- grch38
      } else if (species=="mouse") {
        geneData <- grcm38
      } else stop('Package "annotables" can only be used with human or mouse genes. Please check species name or specify "use_annotables = FALSE".')
      gene_ID_to_chrom_name <-
        data.frame(gene_ID = rownames(counts),
                   chromosome_name = geneData[["chr"]][match(rownames(counts), geneData[[gene_ID]])])
      
    } else if (require(biomaRt)) {
      warning('Used BioMart to retrieve chromosome identity for all gene IDs.\n')
      
      # check parameters for biomaRt based on species
      if (gene_ID=="ensembl") gene_ID <- "ensembl_gene_id"
      if (species=="human") {
        ensemblDataset <- "hsapiens_gene_ensembl"
        if (gene_ID=="symbol") gene_ID <- "hgnc_smybol"
      } else if (species=="mouse") {
        ensemblDataset <- "mmusculus_gene_ensembl"
        if (gene_ID=="symbol") gene_ID <- "mgi_smybol"
      } else {
        ensemblDataset <- paste0(species, "_gene_ensembl")
        if (gene_ID == "symbol") {
          warning('In order to use BioMart with species other than human or mouse, gene IDs must be Ensembl IDs. Defaulting to gene_ID = "ensembl_gene_id".')
          gene_ID <- "ensembl_gene_id"
        }
      }
      
      # check that dataset is present in BioMart ensembl
      if (!ensemblDataset %in% listDatasets(useMart("ensembl"))$dataset) stop('Species "species" not found in BioMart ensembl datasets. Please check name against species abbreviations in listDatasets(useMart("ensembl"))$dataset.')
      
      ensembl <- useMart('ensembl', ensemblDataset)
      gene_ID_to_chrom_name <-
        biomaRt::getBM(
          attributes=c(gene_ID, "chromosome_name"), filters=gene_ID,
          values=rownames(counts), mart=ensembl)
      
    } else {
      stop("Sorry, I'm missing some needed packages.\nPlease install either the annotables or biomaRt package.\n")
    }
    
    chromosome_name <-
      gene_ID_to_chrom_name$chromosome_name[match(rownames(counts), gene_ID_to_chrom_name[,1])]
    x_counts <- Matrix::colSums(counts[which(chromosome_name=="X"), lib_cols], na.rm=TRUE)
    y_counts <- Matrix::colSums(counts[which(chromosome_name=="Y"), lib_cols], na.rm=TRUE)
    ratios <- log((x_counts+1) / (y_counts+1))  # calculate log-transformed ratios of X reads to Y reads; add 1 to each because Y counts can be 0, which yields infinite ratio
    names(ratios) <- colnames(counts)[lib_cols]  # name the vector of ratios with the library IDs
    
    return(ratios)
  }
