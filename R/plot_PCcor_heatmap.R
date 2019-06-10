#' Plot a heatmap of correlations between PCs and other sample variables
#'
#' This is a wrapper for \code{ComplexHeatmap::Heatmap}, designed to make plotting of PC correlations more
#' automated. It also provides a nicer default color palette than the \code{ComplexHeatmap::Heatmap} default.
#' 
#' @param PCcor_result a numeric matrix or data.frame, generally the result of \code{calc_PCcors}. Should have annotation variables as columns and PCs as rows.
#' @param filename a character string. If provided, the function outputs a pdf of the plot, named "{filename}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object (if outputting to a file)
#' @param my_heatmap_cols a vector of color names, typically the result of a call to \code{colorRampPalette}. Default is a blue (low) to white (middle) to red (high) palette.
#' @param center_colors_zero logical, whether to center the color scale on 0. Ignored if \code{my_heatmap_cols} is specified. Defaults to \code{TRUE}.
#' @param row_dendro,col_dendro logical, whether to include the row and/or columns dendrogram(s). Both default to \code{FALSE}.
#' @param remove_all_NA_cols logical, whether to remove variables where all the correlation values are NA. Defaults to \code{TRUE}.
#' @param ... (optional) additional arguments passed to \code{ComplexHeatmap::Heatmap}.
#' @export
plot_PCcor_heatmap <-
  function(PCcor_result, filename=NULL, plotdims=c(9,9),
           my_heatmap_cols, center_colors_zero=TRUE,
           key=TRUE, keysize=0.8, density.info="none", trace="none",
           row_dendro=FALSE, col_dendro=FALSE,
           remove_all_NA_cols=TRUE,
           ...) {
    PCcor_result <- as.matrix(PCcor_result)
    
    checkmate::assert(
      checkmate::check_numeric(PCcor_result)
    )
    
    if (remove_all_NA_cols)
      PCcor_result <- miscHelpers::remove_all_NA_rowcols(PCcor_result, rows=FALSE)
    
    if (missing(my_heatmap_cols)) {
      if (center_colors_zero) {
        val_max_abs <- max(abs(range(PCcor_result, na.rm=TRUE)))
        my_heatmap_cols <-
          circlize::colorRamp2(
            breaks=c(-val_max_abs, 0, val_max_abs),
            colors=c("blue", "white", "red"))
      } else {
        my_heatmap_cols <-
          circlize::colorRamp2(
            breaks=range(PCcor_result, na.rm=TRUE),
            colors=c("blue", "white", "red"))
      }
    }
    
    # open plotting device
    if (!is.null(filename)) {
      pdf(file=filename, w=plotdims[1], h=plotdims[2])
      on.exit(dev.off()) # close plotting device on exit
    }
    
    # generate heatmap
    heatmap_result <-
      ComplexHeatmap::Heatmap(
        matrix=PCcor_result,
        col=my_heatmap_cols,
        name="cor",
        cluster_rows=row_dendro,
        cluster_columns=col_dendro)
    
    print(heatmap_result)
  }
