#' singleton_to_data_frame
#'
#' Removes singletons from a dataset and prunes taxa from a phylogeny.
#'
#' @param data coordinates - in decimal degrees - with the following header: \strong{spp, lat,
#'  long (all lowercase)}
#' @param phylogeny a phylogenetic tree in Newick format
#' @return Both a new dataset as a data frame object and a modified tree without singletons.
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details It can be utilized even if the dataset does not contain any singletons, which may only be discovered later.
#'  It is also applicable in situations where a taxon present in the tree is absent from the dataset, but it is also applicable
#'    in situations where a taxon is absent from the tree and present in the dataset. Finally, it can be used without a phylogeny.
#'      The idea was based on the phytools package (Revell, 2024).
#' @references Revell, L.J., 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other
#'  things). \emph{PeerJ} \strong{12}: e16505.
#' @seealso
#' \code{\link[phytools]{drop.tip}} for pruning a taxon out of phylogenetic trees.
#'
#' \code{\link[PanBioGeo]{MST_node}} for estimating the minimum spanning tree (MST).
#'
#' @export
#' @examples
#' # load example data (coordinates and tree)
#' library(phytools)
#' data(lycipta)
#'
#' # The columns are in the following order: species (spp), longitude (long), latitude (lat)
#' lycipta.coords <- lycipta$coordinates
#' head(lycipta.coords, 20)
#'
#' # phylogenetic tree
#' lycipta.tree <- lycipta$phylogenetic_tree
#'
#' # removal of singletons with a tree
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)], phylogeny = lycipta.tree)
#'
#' # removal of singletons without a tree
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)])

singleton_to_data_frame <- function(data = NULL, phylogeny = NULL){
  # Input validation
  if(is.null(data)){
    stop("The 'data' argument cannot be NULL. Provide a data.frame with columns: species, latitude and longitude.")
  }
  spp <- as.character(data[, 1])
  lat <- data[, 2]
  long <- data[, 3]

  names(lat) <- spp
  names(long) <- spp
  locais <- cbind(lat, long)
  rownames(locais) <- spp

  spp_counts <- table(spp)
  singleton_taxa <- names(spp_counts[spp_counts == 1])

  build_data_df <- function(locais_obj) {
    data.frame(
      spp = rownames(locais_obj),
      long = locais_obj[, 2],
      lat = locais_obj[, 1],
      stringsAsFactors = FALSE
    )
  }

  # without a phylogeny:
  if(is.null(phylogeny)){
    if(length(singleton_taxa) == 0){
      data_df <- build_data_df(locais)
      message("The dataset is OK. No modification was necessary!")
      return(list(data_df = data_df))
    } else {
      locais_new <- locais[!(rownames(locais) %in% singleton_taxa), , drop = FALSE]
      data_df <- build_data_df(locais_new)
      message(paste("Removed", length(singleton_taxa), "species with only one occurrence."))
      return(list(data_df = data_df))
    }
  } else if(!is.null(phylogeny)){ # with a phylogeny
    # Tree validation
    if(class(phylogeny)[1] != "phylo"){
      stop("The 'phylogeny' argument must be an object of class 'phylo' from the 'ape' package.")
    }

    locais_new <- locais
    tree_new <- phylogeny
    n_singleton_removed <- 0L
    n_data_not_tree_removed <- 0L
    n_tree_not_data_removed <- 0L

    if (length(singleton_taxa) > 0) {
      locais_new <- locais_new[!(rownames(locais_new) %in% singleton_taxa), , drop = FALSE]
      drop_singletons <- intersect(singleton_taxa, tree_new$tip.label)
      if (length(drop_singletons) > 0) {
        tree_new <- ape::drop.tip(tree_new, drop_singletons)
      }
      n_singleton_removed <- length(singleton_taxa)
    }

    if (nrow(locais_new) == 0) {
      stop("All taxa were removed during singleton filtering. Please review the occurrence data.")
    }

    data_taxa <- unique(rownames(locais_new))
    tree_taxa <- tree_new$tip.label

    taxa_data_not_tree <- setdiff(data_taxa, tree_taxa)
    taxa_tree_not_data <- setdiff(tree_taxa, data_taxa)

    if (length(taxa_data_not_tree) > 0) {
      locais_new <- locais_new[!(rownames(locais_new) %in% taxa_data_not_tree), , drop = FALSE]
      n_data_not_tree_removed <- length(taxa_data_not_tree)
    }
    if (length(taxa_tree_not_data) > 0) {
      tree_new <- ape::drop.tip(tree_new, taxa_tree_not_data)
      n_tree_not_data_removed <- length(taxa_tree_not_data)
    }

    data_df <- build_data_df(locais_new)
    tree_modified <- n_tree_not_data_removed > 0 || length(intersect(singleton_taxa, phylogeny$tip.label)) > 0

    if (length(singleton_taxa) == 0 && !tree_modified) {
      message("The dataset is OK. No modification was necessary!")
      return(list(data_df = data_df, treeNonMod = tree_new))
    } else {
      message(paste0(
        "Singleton/tree harmonization complete. ",
        "Removed singleton taxa: ", n_singleton_removed, "; ",
        "removed data-only taxa: ", n_data_not_tree_removed, "; ",
        "removed tree-only taxa: ", n_tree_not_data_removed, "."
      ))
      return(list(data_df = data_df, treeMod = tree_new))
    }
  }
}
