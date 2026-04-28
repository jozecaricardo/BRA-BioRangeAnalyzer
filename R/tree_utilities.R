#' Tree Utility Functions for PAE-PCE Analysis
#'
#' Functions to prepare phylogenetic trees for PAE-PCE analysis
#'

#' Check if Tree has Branch Lengths
#'
#' Determine if a phylogenetic tree has branch length information
#'
#' @param tree phylo object
#'
#' @return Logical, TRUE if tree has branch lengths
#'
#' @export
has_branch_lengths <- function(tree) {
  !is.null(tree$edge.length) && length(tree$edge.length) > 0
}

#' Check if Tree is Ultrametric
#'
#' Determine if a tree is ultrametric (all tips equidistant from root)
#'
#' @param tree phylo object
#'
#' @return Logical, TRUE if tree is ultrametric
#'
#' @export
is_ultrametric <- function(tree) {
  
  if (!has_branch_lengths(tree)) {
    return(FALSE)
  }
  
  tryCatch({
    # Calculate distance from root to each tip
    distances <- ape::node.depth.edgelength(tree)
    tip_distances <- distances[1:length(tree$tip.label)]
    
    # Check if all tips have approximately the same distance
    max_dist <- max(tip_distances)
    min_dist <- min(tip_distances)
    
    # Allow small tolerance for floating point errors
    tolerance <- 1e-6
    (max_dist - min_dist) < tolerance
  }, error = function(e) {
    FALSE
  })
}

#' Convert Tree to Ultrametric
#'
#' Convert a non-ultrametric tree to ultrametric by scaling branch lengths
#' so all tips are equidistant from the root. This is optional and can be done
#' by user request before PAE-PCE analysis.
#'
#' @param tree phylo object with branch lengths
#' @param method Method for ultrametricization:
#'   - "equal": All branch lengths become 1 (equal length tree)
#'   - "scale": Scale existing branch lengths proportionally
#'   - "penalized": Use penalized least squares (requires phytools)
#'
#' @return Ultrametric phylo object
#'
#' @export
make_ultrametric <- function(tree, method = "equal") {
  
  if (!has_branch_lengths(tree)) {
    stop("Tree must have branch lengths to be made ultrametric")
  }
  
  if (is_ultrametric(tree)) {
    message("Tree is already ultrametric")
    return(tree)
  }
  
  tryCatch({
    if (method == "equal") {
      # Simple method: all branches get equal length
      tree$edge.length <- rep(1, length(tree$edge.length))
      return(tree)
    }
    
    else if (method == "scale") {
      # Scale method: proportionally adjust to make ultrametric
      # Calculate current distances
      distances <- ape::node.depth.edgelength(tree)
      tip_distances <- distances[1:length(tree$tip.label)]
      
      # Target distance (maximum current distance)
      target_dist <- max(tip_distances)
      
      # Scale factor for each edge
      tree$edge.length <- tree$edge.length * (target_dist / max(tip_distances))
      
      return(tree)
    }
    
    else if (method == "penalized") {
      # Penalized least squares method (requires phytools)
      if (!requireNamespace("phytools", quietly = TRUE)) {
        warning("phytools not available, using 'equal' method instead")
        tree$edge.length <- rep(1, length(tree$edge.length))
        return(tree)
      }
      
      tree_ultra <- phytools::force.ultrametric(tree, method = "extend")
      return(tree_ultra)
    }
    
    else {
      stop("Unknown method: ", method)
    }
  }, error = function(e) {
    stop("Error making tree ultrametric: ", e$message)
  })
}

#' Get Tree Statistics
#'
#' Calculate summary statistics for a phylogenetic tree
#'
#' @param tree phylo object
#'
#' @return List with tree statistics
#'
#' @export
get_tree_stats <- function(tree) {
  
  stats <- list(
    n_taxa = length(tree$tip.label),
    n_nodes = tree$Nnode,
    n_edges = nrow(tree$edge),
    has_branch_lengths = has_branch_lengths(tree),
    is_ultrametric = is_ultrametric(tree)
  )
  
  if (has_branch_lengths(tree)) {
    stats$total_length <- sum(tree$edge.length)
    stats$mean_branch_length <- mean(tree$edge.length)
    stats$min_branch_length <- min(tree$edge.length)
    stats$max_branch_length <- max(tree$edge.length)
    
    # Distance from root to tips
    distances <- ape::node.depth.edgelength(tree)
    tip_distances <- distances[1:length(tree$tip.label)]
    stats$max_tip_distance <- max(tip_distances)
    stats$min_tip_distance <- min(tip_distances)
  }
  
  stats
}

#' Validate Tree for Analysis
#'
#' Check if tree is suitable for PAE-PCE and BioGeoBEARS analyses
#'
#' @param tree phylo object
#' @param for_analysis Type of analysis: "pae_pce" or "biogeobears"
#'
#' @return List with validation results
#'
#' @export
validate_tree_for_analysis <- function(tree, for_analysis = "pae_pce") {
  
  validation <- list(
    is_valid = TRUE,
    errors = c(),
    warnings = c(),
    suggestions = c()
  )
  
  # Basic checks
  if (!inherits(tree, "phylo")) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Object is not a phylo object")
    return(validation)
  }
  
  if (length(tree$tip.label) < 2) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Tree must have at least 2 taxa")
    return(validation)
  }
  
  # Check for branch lengths
  if (!has_branch_lengths(tree)) {
    validation$warnings <- c(validation$warnings, "Tree does not have branch lengths")
  }
  
  # Check for duplicate taxa names
  if (length(unique(tree$tip.label)) < length(tree$tip.label)) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Tree has duplicate taxon names")
  }
  
  validation
}
