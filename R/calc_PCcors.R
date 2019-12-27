#' Calculate correlations of principal components with quality or annotation variables
#'
#' This is a helper function to facilitate calculating correlations of sample-specific variables with principal components (generally from RNAseq data). By default, it uses Spearman (rank) correlations for continuous variables and intraclass correlations (as implemented in \code{ICC::ICCbare}) for categorical variables. The correlation method can be changed for continuous variables, but not currently for categorical variables.
#' @param PCA_result result of a principal component analysis, generally of gene expression data. Typically the output of \code{prcomp} or \code{calc_PCAs}. Can also be a matrix with samples in rows and dimensions in columns.
#' @param annotation a data frame containing annotation data for the samples. May include clinical data, sample quality metrics, etc.
#' @param PCs numeric vector of principal component axes to include in correlation calculations. Defaults to 1:10, which will calculate for the first 10 PCs. Any PCs specified that are not found in the PCA object will not be compared.
#' @param id_col name or number of the column of \code{annotation} containing the library identifiers. These are matched to the rownames of \code{PCA_result}. Ignored if \code{PCA_result} does not have rownames.
#' @param var_cols numbers or names of columns to include in the correlation calculations. If not specified, all columns will be included, subject to other exclusion criteria.
#' @param ignore_unique_nonnumeric logical, whether to drop columns from annotation if they contain unique non-numeric values. Correlations for such variables are meaningless. Defaults to TRUE.
#' @param ignore_invariant logical, whether to drop columns from annotation if all non-NA values are identical. Correlations for such variables are meaningless. Defaults to TRUE.
#' @param date_as_numeric logical, whether to treat data of class "Date" and "POSIXt" as numeric. If set to FALSE, dates are treated as categorical variables.
#' @param min_libs number, the minimum number of libraries containing non-NA values for a variable. Variables in \code{annotation} with fewer non-NA values will be dropped. Defaults to 5.
#' @param cont_method character, the name of the correlation coefficient to use for continuous variables. Passed to \code{stats::cor}, and must be one of "pearson", "kendall", or "spearman", or abbreviations thereof. Defaults to "spearman".
#' @param cat_method character, the name of the correlation coefficient to use for categorical variables. Currently, the only acceptable option is "ICC", which uses the intraclass correlation coefficient as implemented in \code{ICC::ICCbare}.
#' @param ... (optional) additional arguments passed to \code{cor} or other functions.
#' @export
#' @return a matrix of correlation coefficients, wih the column and row names reflecting the PC axes and annotation variables for which correlations were calculated.
calc_PCcors <-
  function(PCA_result, annotation,
           PCs=1:10, id_col="libid",
           var_cols, ignore_unique_nonnumeric=TRUE, ignore_invariant=TRUE, date_as_numeric=TRUE,
           min_libs=5, cont_method="spearman", cat_method="ICC",
           ...) {
    if (inherits(PCA_result, "princomp")) {
      PCA_result <- PCA_result$scores
    } else if (inherits(PCA_result, "prcomp")) {
      PCA_result <- PCA_result$x
    } else if (!inherits(PCA_result, "matrix")) stop("Class of object PCA_result not recognized.")

    # drop objects from PCA matrix if not found in annotation object
    if (is.null(rownames(PCA_result))) {
      if (nrow(PCA_result) != nrow(annotation)) {
        stop("PCA_result must have rownames, or must have the same number of rows as annotation.")
      } else cat("No rownames in PCA_result. Assuming that PCA_result and annotation are in the same order.")
    } else {
      if (!(id_col) %in% colnames(annotation))
        stop(paste("Annotation object does not contain expected column of library identifiers:", id_col))
      PCA_result <- PCA_result[na.omit(match(annotation[,id_col, drop=TRUE], rownames(PCA_result))),]
      annotation <- annotation[na.omit(match(rownames(PCA_result), annotation[,id_col, drop=TRUE])),]
    }

    ## check that PCs specified are in PCA_result
    if (!all(PCs %in% 1:ncol(PCA_result))) {
      warning("Some specified PC axes were not found in PCA_result.")
      PCs <- intersect(PCs, 1:ncol(PCA_result))
    }

    ### drop columns from annotation object

    # keep columns specified by var_cols (or all if unspecified)
    if (exists("var_cols")) annotation <- annotation[,var_cols]

    # drop columns containing only unique non-numeric values
    if (ignore_unique_nonnumeric) {
      annotation <-
        annotation[
          , sapply(annotation,
                   function(x) {
                     is.numeric(x) |
                       length(unique(na.omit(x)))!=length(na.omit(x))})]
    }

    # drop columns containing all identical values
    if (ignore_invariant) {
      annotation <-
        annotation[
          , sapply(annotation,
                   function(x) {
                       length(unique(na.omit(x))) > 1})]
    }

    # convert Date columns to numeric (if specified)
    if (date_as_numeric) {
      for (
        j in
        colnames(annotation)[
          sapply(annotation, inherits, "Date") |
          sapply(annotation, inherits, "POSIXt")]) {
        annotation[,j] <- as.numeric(annotation[,j, drop=TRUE])
      }
    }

    # drop columns with too few non-NA values
    annotation <-
      annotation[
        , sapply(annotation,
                 function(x) {sum(!is.na(x)) >= min_libs})]

    # generate empty matrix to store correlation data
    if (is.null(colnames(PCA_result)))
      colnames(PCA_result) <- paste0("PC", 1:ncol(PCA_result))
    PCcors <-
      matrix(data=NA,
        nrow=length(PCs),
        ncol=ncol(annotation),
        dimnames=list(colnames(PCA_result)[PCs], colnames(annotation)))

    # calculate correlations of all variables with PCs, using specified methods for continuous
    for (i in colnames(PCcors)) {
      if (is.numeric(annotation[,i])) {
        PCcors[,i] <-
          cor(annotation[,i, drop=TRUE],
              PCA_result[,PCs],
              method=cont_method, use="pairwise")
      } else {
        if (cat_method=="ICC") {
          for (j in rownames(PCcors)) {
            data.tmp <-
              data.frame(
                annotation[,i], PCA_result[,j])
            colnames(data.tmp) <- c(i,j)

            capture.output(suppressWarnings(  # suppress output from ICCbare
              PCcors[j,i] <-
                ICC::ICCbare(data=data.tmp, x=i, y=j)))
          }
        } else stop ("Only method 'ICC' is currently supported for categorical variables.")
      }
    }

    PCcors
  }
