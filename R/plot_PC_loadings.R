#' Plot PC loadings
#' 
#' Plot the loadings of each variable on each principal component. Each variable is plotted as a colored line.
#' In the future, I may expand this function to take prcomp and princomp objects, and to provide additional
#' handles for manipulating the plots (e.g. number of variables and/or number of PCs to plot), outputting
#' to pdf, etc.
#' @param rotation_matrix a PCA rotation or loading matrix. Should have genes in rows and PCs in columns.
#' @import ggplot2
#' @export
#' @usage \code{plot_PC_loadings(rotation_matrix)}
plot_PC_loadings <- function(rotation_matrix) {
  rotation_matrix <-
    reshape2::melt(rotation_matrix, varnames=c("gene", "PC"), value.name="loading")
  ggplot(data=rotation_matrix, mapping=aes(x=PC, y=loading, group=gene, color=gene)) +
    geom_line() + guides(color=FALSE)
}
