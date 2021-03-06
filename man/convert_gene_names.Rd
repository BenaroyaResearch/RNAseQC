% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_gene_names.R
\name{convert_gene_names}
\alias{convert_gene_names}
\title{Convert gene identifiers}
\usage{
convert_gene_names(
  genes,
  input_type,
  output_type,
  biotype = NULL,
  use_annotables = TRUE
)
}
\arguments{
\item{genes}{character vector containing gene identifiers to convert}

\item{input_type}{the gene identifier type to convert from. Must match a corresponding identifier type in \code{annotables} or \code{biomaRt}.}

\item{output_type}{the gene identifier type to convert to Must match a corresponding identifier type in \code{annotables} or \code{biomaRt}.}

\item{biotype}{(optional) character vector containing biotypes of genes to retain. If provided, genes with biotype not matching elements of this vector return NA. If NULL, all genes are retained.}

\item{use_annotables}{logical, whether to use the annotables package. If annotables is not installed, the function defaults to using \code{biomaRt}.}
}
\value{
a vector of gene symbols, of length equivalent to the length of genes.
}
\description{
Convert a vector of gene identifiers from one standard to another. The input type of the identifiers
is specified using \code{input_type}, and the output type is specified using \code{output_type}.
By default, this function uses the \code{annotables} package. If annotables
is not installed, or if the user specifies not to use it, \code{biomaRt} is used instead. Genes can be
filtered by biotype; any genes not matching the specified biotype(s) will return NA.
Currently, support is only provided for human gene identifiers.
}
