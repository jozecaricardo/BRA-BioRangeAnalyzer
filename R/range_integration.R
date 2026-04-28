#' Range Extrapolation Integration Functions
#'
#' Wrapper functions to integrate actual range extrapolation methods into the Shiny app
#'

#' Load Phylogenetic Tree
#'
#' Load a phylogenetic tree from file (Newick or Nexus format)
#'
#' @param tree_file Path to tree file
#'
#' @return phylo object
#'
#' @export
load_phylogenetic_tree <- function(tree_file) {
  
  tryCatch({
    # Try reading as Newick first
    tree <- tryCatch({
      ape::read.tree(tree_file)
    }, error = function(e) {
      # Try reading as Nexus
      ape::read.nexus(tree_file)
    })
    
    # Check if tree has branch lengths
    if (is.null(tree$edge.length)) {
      warning("Tree does not have branch lengths")
    }
    
    tree
  }, error = function(e) {
    stop("Error reading tree file: ", e$message)
  })
}

#' Validate Phylogenetic Tree
#'
#' Check if tree is valid
#'
#' @param tree phylo object
#'
#' @return List with validation results
#'
#' @export
validate_phylogenetic_tree <- function(tree) {
  
  list(
    is_valid = inherits(tree, "phylo"),
    n_taxa = length(tree$tip.label),
    has_branch_lengths = !is.null(tree$edge.length),
    taxa = tree$tip.label
  )
}

#' Prepare Occurrence Data
#'
#' Standardize column names for range functions
#'
#' @param occurrence_data Data frame with species, longitude, latitude
#'
#' @return Data frame with standardized columns
#'
#' @export
prepare_occurrence_data <- function(occurrence_data) {
  
  # Standardize column names
  names(occurrence_data) <- tolower(names(occurrence_data))
  
  # Rename to match function requirements
  if ("species" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "species"] <- "spp"
  }
  if ("longitude" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "longitude"] <- "long"
  }
  if ("latitude" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "latitude"] <- "lat"
  }
  
  # Return only required columns
  occurrence_data[, c("spp", "long", "lat")]
}

#' Extract Presence-Absence Matrix
#'
#' Extract pres_abs from range function results
#'
#' @param result Result from range extrapolation function
#'
#' @return Presence-absence matrix
#'
#' @export
extract_pres_abs <- function(result) {
  
  if (is.list(result) && "pres_abs" %in% names(result)) {
    return(result$pres_abs)
  }
  
  if (is.matrix(result)) {
    return(result)
  }
  
  stop("Cannot extract presence-absence matrix")
}

#' Extract Geometry from Range Results
#'
#' Extract spatial geometry (polygons) from range results
#'
#' @param result Result from range extrapolation function
#'
#' @return Spatial geometry object or NULL
#'
#' @export
extract_geometry <- function(result) {
  
  if (is.list(result) && "geometry" %in% names(result)) {
    return(result$geometry)
  }
  
  NULL
}

#' Create BioGeoBEARS Geography File
#'
#' Convert presence-absence matrix to BioGeoBEARS .data format
#'
#' @param pres_abs Presence-absence matrix
#' @param method Method name ("BUFF", "MPC", "MST")
#' @param output_dir Output directory
#'
#' @return Path to created file
#'
#' @export
create_biogeobears_file <- function(pres_abs, method = "BUFF", output_dir = ".") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Remove empty ranges (all zeros)
  pres_abs_clean <- pres_abs[rowSums(pres_abs) > 0, ]
  
  # Create output file
  output_file <- file.path(output_dir, paste0("pres_abs_", method, "_geog.data"))
  
  # Write header
  n_species <- ncol(pres_abs_clean)
  n_ranges <- nrow(pres_abs_clean)
  range_names <- paste(1:n_ranges, collapse = " ")
  
  # Convert matrix to binary strings
  binary_strings <- apply(pres_abs_clean, 1, function(row) paste(row, collapse = ""))
  
  # Write file
  writeLines(c(
    paste(n_species, n_ranges, range_names),
    paste(colnames(pres_abs_clean), binary_strings, sep = " ")
  ), output_file)
  
  output_file
}

#' Convert Geometry to SF for Leaflet
#'
#' Convert raster or terra geometry to sf polygons for Leaflet visualization
#'
#' @param geometry Spatial geometry object (raster, SpatRaster, or sf)
#'
#' @return sf object with polygons
#'
#' @export
convert_geometry_to_sf <- function(geometry) {
  
  if (is.null(geometry)) {
    return(NULL)
  }
  
  tryCatch({
    # If already sf, return as is
    if (inherits(geometry, "sf")) {
      return(geometry)
    }
    
    # If terra SpatRaster, convert to polygons
    if (inherits(geometry, "SpatRaster")) {
      # Convert raster to polygons
      poly <- terra::as.polygons(geometry)
      # Convert to sf
      sf_poly <- sf::st_as_sf(poly)
      return(sf_poly)
    }
    
    # If terra SpatVector, convert to sf
    if (inherits(geometry, "SpatVector")) {
      sf_poly <- sf::st_as_sf(geometry)
      return(sf_poly)
    }
    
    # If raster package raster, convert
    if (inherits(geometry, "RasterLayer") || inherits(geometry, "RasterStack")) {
      # Convert to polygons
      poly <- raster::rasterToPolygons(geometry)
      # Convert to sf
      sf_poly <- sf::st_as_sf(poly)
      return(sf_poly)
    }
    
    # If sp object, convert to sf
    if (inherits(geometry, "Spatial")) {
      sf_poly <- sf::st_as_sf(geometry)
      return(sf_poly)
    }
    
    warning("Unsupported geometry type")
    NULL
  }, error = function(e) {
    warning("Error converting geometry: ", e$message)
    NULL
  })
}

#' Add Extrapolation Polygons to Leaflet Map
#'
#' Add extrapolation polygons to an existing leaflet map
#'
#' @param map Leaflet map object
#' @param geometry Spatial geometry object
#' @param method Method name ("BUFF", "MPC", "MST")
#' @param color Color for polygons
#' @param opacity Opacity of polygons
#'
#' @return Updated leaflet map
#'
#' @export
add_extrapolation_polygons <- function(map, geometry, method = "BUFF", 
                                        color = "blue", opacity = 0.5) {
  
  if (is.null(geometry)) {
    return(map)
  }
  
  tryCatch({
    # Convert to sf if needed
    sf_geom <- convert_geometry_to_sf(geometry)
    
    if (is.null(sf_geom)) {
      return(map)
    }
    
    # Add polygons to map
    map <- map %>%
      leaflet::addPolygons(
        data = sf_geom,
        color = color,
        weight = 2,
        opacity = opacity,
        fillOpacity = opacity * 0.5,
        popup = paste(method, "extrapolation"),
        group = method
      )
    
    return(map)
  }, error = function(e) {
    warning("Error adding polygons: ", e$message)
    return(map)
  })
}
