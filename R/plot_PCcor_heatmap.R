#' Plot a heatmap of correlations between PCs and other sample variables
#'
#' This is a wrapper for \code{heatmap.2}, designed to make plotting of PC correlations more
#' automated. It also provides a nicer default color palette than the \code{heatmap.2} default.
#' 
#' @param PCcor_result a numeric matrix or data.frame, generally the result of \code{calc_PCcors}. Should have annotation variables as columns and PCs as rows.
#' @param filename a character string. If provided, the function outputs a pdf of the plot, named "{filename}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @param my_heatmap_cols a vector of color names, typically the result of a call to \code{colorRampPalette}. Default is a blue (low) to red (high) palette.
#' @param key,keysize,density.info,trace parameters passed to \code{heatmap.2}. Default values here are to override default behavior of \code{heatmap.2}.
#' @param row_dendro,col_dendro logical, whether to include the row and/or columns dendrogram(s). Both default to \code{FALSE}.
#' @param ... (optional) additional arguments passed to \code{heatmap.2}.
#' @export
#' @usage \code{
#' plot_PCcor_heatmap(
#'      PCcor_result, filename=NULL, plotdims=c(9,9),
#'      my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
#'      key=TRUE, keysize=0.8, density.info="none", trace="none",
#'      row_dendro=FALSE, col_dendro=FALSE,
#'      ...)}
plot_PCcor_heatmap <-
  function(PCcor_result, filename=NULL, plotdims=c(9,9),
           my_heatmap_cols=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           key=TRUE, keysize=0.8, density.info="none", trace="none",
           row_dendro=FALSE, col_dendro=FALSE,
           ...) {
    
    # open plotting device
    if (!is.null(filename)) {
      pdf(file=filename, w=plotdims[1], h=plotdims[2])
      on.exit(dev.off()) # close plotting device on exit
    } else quartz(plotdims[1],plotdims[2])
    
    # generate heatmap
    
    gplots::heatmap.2(
      PCcor_result, col=my_heatmap_cols,
      Rowv=row_dendro, Colv=col_dendro,
      key=key, keysize=keysize, margins=c(8,12),
      trace=trace, density.info=density.info,
      ...)
    
    
    
  }
