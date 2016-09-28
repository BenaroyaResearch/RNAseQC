#' Plot saturation curves of genes based on reads object
#' 
#' Plot the saturation of gene detection based on rarefaction of the mapped read counts from a library.
#' This function takes the output of \code{estimate_saturation}.
#' @param saturation a data frame, typically the output of \code{estimate_saturation}. Should contain columns "sample" with sample identifiers, "depth" with the varied simulated sequencing depths, and "sat" with the number of genes detected at that depth.
#' @param plot_points logical, whether to include points in the plot.
#' @param plot_terminal_points, whether to plot the maximum-depth point for each library. This helps to visualize the actual distribution of observed counts by library.
#' @param plot_lines logical, whether to include (unsmoothed) lines in the plot. It is not recommended to include both \code{plot_lines} and \code{plot_smooth}, as this will make plots difficult to read.
#' @param plot_smooth logical, whether to include lowess curves in the plot. It is not recommended to include both \code{plot_lines} and \code{plot_smooth}, as this will make plots difficult to read.
#' @param color_points_by_var,color_lines_by_var,color_terminal_points_by_var,color_smooth_by_var variable to use for plotting points, lines, terminal points, and/or smooths by a discrete variable. Must be either the name of a column in \code{saturation}, or a named vector containing the values, with element names corresponding to values in \code{saturation$sample}.
#' @import ggplot2
#' @import dplyr
#' @export
#' @usage \code{
#' plot_saturation_curve(saturation,
#'                       plot_points=TRUE,
#'                       plot_lines=TRUE,
#'                       plot_terminal_points=TRUE,
#'                       plot_smooth=FALSE,
#'                       color_points_by_var=NULL, my_point_cols=NULL
#'                       color_lines_by_var=NULL, my_line_cols=NULL
#'                       color_terminal_points_by_var=NULL, my_terminal_point_cols=NULL
#'                       color_smooths_by_var=NULL, my_smooth_cols=NULL)}
plot_saturation_curve <- function(saturation,
                                  plot_points=TRUE,
                                  plot_lines=TRUE,
                                  plot_terminal_points=TRUE,
                                  plot_smooth=FALSE) {
  satplot <- ggplot(data=saturation, mapping=aes(x=depth, y=sat, group=sample, color=sample)) +
    labs(x="Reads", y="Genes above threshold") +
    guides(color=FALSE)
  
  if (plot_points) satplot <- satplot + geom_point()
  if (plot_lines) satplot <- satplot + geom_line()
  if (plot_smooth) satplot <- satplot + geom_smooth(se=FALSE)
  
  if (plot_terminal_points) {
    saturation.maxdepth <- saturation %>%
      filter(!is.na(sat)) %>%
      group_by(sample) %>%
      summarise(max.depth=max(depth))
    saturation.terminal <- saturation[sapply(1:nrow(saturation.maxdepth), function(x) {
      which((saturation$sample==saturation.maxdepth$sample[x]) &
              (saturation$depth == saturation.maxdepth$max.depth[x]))
      }),]
    satplot <- satplot + geom_point(data=saturation.terminal, mapping=aes(x=depth, y=sat),
                                    color="black")
  }
  print(satplot)
}
