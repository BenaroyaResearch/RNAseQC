#' Convert gene identifiers from mouse to human or human to mouse, using orthologs
#'
#' Convert a vector of symbols from one species to another. Current options
#' exist only for mouse to human and  identifiers from one standard to another.
#' Mouse identifiers will be MGI symbols. Human identifiers will be HGNC
#' symbols. Multiple genes may be returned for any input gene, based on the
#' possible one-to-many mapping of gene orthologs. This function is adapted
#' from https://www.biostars.org/p/9567892/, and uses the MGI list of
#' orthologs at https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
#' @param genes character vector containing gene identifiers to convert
#' @param species_from character, the species of the input gene identifiers. Must be "mouse" or "human"; partial matches are accepted.
#' @param species_to character, the species of the output gene identifiers. Must be "mouse" or "human"; partial matches are accepted.
#' @param return_format character, the format of the output gene identifiers. Must be "orthologs_only" or "orig_and_ortholog"; partial matches are accepted.
#' @export
#' @return a vector of gene symbols, of length similar to the input genes, minus any genes that do not have orthologs in `species_to`.
convert_gene_symbols_species <-
  function(genes, species_from, species_to, return_format = "orthologs_only") {
    PATH_GENE_NAME_CONVERSION <-
      "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
    
    # check inputs
    checkmate::assert(
      checkmate::check_character(genes),
      checkmate::check_character(species_from),
      checkmate::check_character(species_to),
      checkmate::check_character(return_format)
    )
    
    species_from <- match.arg(species_from, c("mouse", "human"))
    species_to <- match.arg(species_to, c("mouse", "human"))
    if (species_from == species_to) {
      stop("species_from and species_to must be different")
    }
    
    return_format <-
      match.arg(return_format, c("orthologs_only", "orig_and_ortholog"))
    
    mouse_human_genes <-
      read.delim(PATH_GENE_NAME_CONVERSION, sep="\t")
    
    genes_out <-
      data.frame(
        gene_from = character(),
        gene_to = character())
    
    # replace "mouse" with "mouse, laboratory" to match the species name in the dataset
    species_from <- str_replace(species_from, "mouse", "mouse, laboratory")
    species_to <- str_replace(species_to, "mouse", "mouse, laboratory")
    
    for (gene_from in genes) {
      class_key <-
        mouse_human_genes %>%
        dplyr::filter(
          Symbol %in% gene_from,
          Common.Organism.Name %in% species_from) %>%
        dplyr::pull(DB.Class.Key)
      if (!identical(class_key, integer(0))) {
        genes_to <-
          mouse_human_genes %>%
          dplyr::filter(
            DB.Class.Key %in% class_key,
            Common.Organism.Name %in% species_to) %>%
          dplyr::pull(Symbol)
        for (gene_to in genes_to) {
          genes_out <-
            rbind(genes_out,
                  list("gene_from" = gene_from, "gene_to" = gene_to))
        }
      }
    }
    
    if (return_format == "orthologs_only")
      genes_out <- unique(genes_out[,2])
    
    return(genes_out)
  }
