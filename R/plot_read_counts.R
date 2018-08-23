#' Plot read counts by library
#' 
#' Output plots of read counts by library, at coarse and fine scales. Libraries areordered by the total
#' number of reads. These plots can optionally be output to pdfs by specifing \code{file_prefix}. The first
#' plot shows all libraries, with a line at \code{threshold_line} million reads. The second plot zooms in
#' on the \code{n_lowcount} libraries with the lowest counts.
#' @param metrics matrix or data frame containing values of metrics. Should have metrics in columns and libraries in rows.
#' @param file_prefix a character string. If provided, the function outputs pdfs of the plots, named "{file_prefix}_{plot_name}.pdf". If not provided, the function prints to a plotting window.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdfs, or the size of the plotting windows.
#' @param threshold_line numeric, the values (in millions of reads) at which to plot a horizontal line.
#' @param n_lowcount numeric, the number of libraries to include in the plot of low-count libraries
#' @param id_col numeric or character, the number or name of the column containing the library identifiers. Used to plot identifiers of low-count libraries.
#' @export
#' @usage \code{
#' plot_read_counts(
#'      metrics,
#'      file_prefix=NULL, plotdims=c(9,6),
#'      threshold_line=5, n_lowcount=20,
#'      id_col="lib.id"
#'      )}
plot_read_counts <-
  function(metrics,
           file_prefix=NULL, plotdims=c(9,6),
           threshold_line=5, n_lowcount=20,
           id_col="lib.id") {
    metrics <- arrange(metrics, fastq_total_reads)
    
    if (!is.null(file_prefix)) {
      pdf(file=paste0(file_prefix, ".read_count_all_libs.pdf"), w=plotdims[1], h=plotdims[2])
      on.exit(while ("pdf" %in% names(dev.list())) dev.off()) # close plotting device on exit (mostly important for errors that could leave pdf output open)
    } else quartz(w=plotdims[1], h=plotdims[2])
    barplot(metrics[,"fastq_total_reads"]/10^6, main="Read count for all libraries",
            xlab="libraries", ylab = "total reads (in millions)")
    abline(h=threshold_line)
    
    if (!is.null(file_prefix)) {
      pdf(file=paste0(file_prefix, ".read_count_lowcount_libs.pdf"), w=plotdims[1], h=plotdims[2])
      on.exit(while ("pdf" %in% names(dev.list())) dev.off()) # close plotting device on exit (mostly important for errors that could leave pdf output open)
    } else quartz(w=plotdims[1], h=plotdims[2])
    barplot(metrics[1:n_lowcount,"fastq_total_reads"]/10^6, main="Read count for low-count libraries",
            names.arg = metrics[1:n_lowcount,id_col, drop=TRUE],
            ylab = "total reads (in millions)",las=2)
    abline(h=threshold_line)
  }
