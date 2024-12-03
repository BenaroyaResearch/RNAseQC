#' Plot saturation curves of genes based on reads object
#' 
#' Plot the saturation of gene detection based on rarefaction of the mapped read counts from a library.
#' This function takes the output of \code{estimate_saturation}.
#' @param saturation a data frame, typically the output of \code{estimate_saturation}. Should contain columns "sample" with sample identifiers, "depth" with the varied simulated sequencing depths, and "sat" with the number of genes detected at that depth.
#' @param design an optional data frame containing annotation data for use in coloring points, lines, and/or smooths. Must include a column matching column "sample" in \code{saturation}, and any columns specified by the "color_by" arguments.
#' @param design_id_col the number or name of the column containing identifiers that match the "sample" column in \code{saturation}. Defaults to "libid".
#' @param plot_points logical, whether to include points in the plot.
#' @param plot_terminal_points, whether to plot the maximum-depth point for each library. This helps to visualize the actual distribution of observed counts by library.
#' @param plot_lines logical, whether to include (unsmoothed) lines in the plot. It is not recommended to include both \code{plot_lines} and \code{plot_smooths}, as this will make plots difficult to read.
#' @param plot_smooths logical, whether to include lowess curves in the plot. It is not recommended to include both \code{plot_lines} and \code{plot_smooths}, as this will make plots difficult to read.
#' @param color_points_by_var,color_lines_by_var,color_terminal_points_by_var,color_smooths_by_var variable to use for plotting points, lines, terminal points, and/or smooths by a discrete variable. Must be either the name of a column in \code{saturation}, or a column in \code{design}.
#' @param log_transform_depth logical, whether to plot the read counts on a log10 scale. Defaults to FALSE.
#' @param log_transform_genes logical, whether to plot the number of genes detected on a log10 scale. Defaults to FALSE.
#' @param my_cols character vector of color names for use in plotting. If not specified, \code{ggthemes::colorblind_pal()} is used. If the number of colors provided is less than the number of colors needed, additional colors are interpolated using \code{colorRampPalette}.
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym
#' @export
plot_saturation_curve <-
  function(saturation,
           design, design_id_col = "libid",
           plot_points = FALSE, color_points_by_var = NULL,
           plot_lines = TRUE, color_lines_by_var = NULL,
           plot_terminal_points = TRUE, color_terminal_points_by_var = NULL,
           plot_smooths = FALSE, color_smooths_by_var = NULL,
           log_transform_depth = FALSE, log_transform_genes = FALSE,
           my_cols = NULL) {
    if (all(is.null(color_points_by_var), is.null(color_lines_by_var),
            is.null(color_terminal_points_by_var), is.null(color_smooths_by_var))) {
      satplot <- ggplot(data = saturation,
                        mapping = aes(x = depth, y = sat, group = sample, color = sample)) +
        labs(x = "Reads", y = "Genes above threshold")
    } else {
      saturation <-
        merge(saturation,
              design[
                ,c(design_id_col, color_points_by_var, color_lines_by_var,
                   color_terminal_points_by_var, color_smooths_by_var)],
              by.x = "sample", by.y = design_id_col, all.x = TRUE)
      satplot <- ggplot(data = saturation, mapping = aes(x = depth, y = sat, group = sample)) +
        labs(x = "Reads", y = "Genes above threshold")
    }
    
    n_col <- 8
    
    if (plot_points) {
      if (!is.null(color_points_by_var)) {
        satplot <- satplot + geom_point(mapping = aes(color = !!sym(color_points_by_var)))
        n_col <- max(n_col, length(unique(saturation[, color_points_by_var])))
      } else {
        satplot <- satplot + geom_point()
      }
    }
    
    if (plot_lines) {
      if (!is.null(color_lines_by_var)) {
        satplot <- satplot + geom_line(mapping = aes(color = !!sym(color_lines_by_var)))
        n_col <- max(n_col, length(unique(saturation[, color_lines_by_var])))
      } else {
        satplot <- satplot + geom_line()
      }
    }
    
    if (plot_smooths) {
      if (!is.null(color_smooths_by_var)) {
        satplot <- satplot + geom_smooth(mapping = aes(color = !!sym(color_smooths_by_var)))
        n_col <- max(n_col, length(unique(saturation[, color_smooths_by_var])))
      } else {
        satplot <- satplot + geom_smooth(se = FALSE)
      }
    }
    
    if (plot_terminal_points) {
      saturation.maxdepth <- saturation %>%
        filter(!is.na(sat)) %>%
        group_by(sample) %>%
        summarise(max.depth = max(depth))
      saturation.terminal <- saturation[sapply(1:nrow(saturation.maxdepth), function(x) {
        which((saturation$sample == saturation.maxdepth$sample[x]) &
                (saturation$depth == saturation.maxdepth$max.depth[x]))}),]
      
      if (!is.null(color_terminal_points_by_var)) {
        satplot <- satplot +
          geom_point(data = saturation.terminal,
                     mapping = aes(x = depth, y = sat, color = !!sym(color_terminal_points_by_var)))
        n_col <- max(n_col, length(unique(saturation.terminal[,color_terminal_points_by_var])))
      } else {
        satplot <- satplot +
          geom_point(data = saturation.terminal, mapping = aes(x = depth, y = sat),
                     color = "black")
      }
    }
    
    if (all(is.null(color_points_by_var), is.null(color_lines_by_var),
            is.null(color_terminal_points_by_var), is.null(color_smooths_by_var))) {
      satplot <- satplot + guides(color=FALSE)
    } else {
      if (is.null(my_cols)) {
        my_cols <-
          colorRampPalette(ggthemes::colorblind_pal()(8))(n_col)
      }
      satplot <- satplot + scale_color_manual(values=my_cols)
    }
    
    if (log_transform_depth) satplot <- satplot + scale_x_log10()
    if (log_transform_genes) satplot <- satplot + scale_y_log10()
    
    print(satplot)
  }
