#' Calculate Richness by Irregular Polygons from Extrapolations
#'
#' This function takes an extrapolation (as sf geometry) and calculates richness
#' by intersecting it with irregular polygons. Used for BioGeoBEARS analysis.
#'
#' @param extrapolation_sf sf object with extrapolation geometries (polygons or grid cells)
#' @param bins_shapefile sf object or file path with irregular polygons
#' @param bin_id_column Character name of the column containing polygon IDs
#' @param crs_input Integer CRS code (default: 4326 for WGS84)
#'
#' @return List containing:
#'   - bins_richness: Data frame with richness per bin
#'   - species_per_bin: Data frame with species presence per bin
#'
#' @export
calcRichness_extrapolation_irregular_bins <- function(extrapolation_sf,
                                                       bins_shapefile,
                                                       bin_id_column = "bin_id",
                                                       crs_input = 4326) {

  # Check for required packages
  required_packages <- c("sf", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' not found. Please install: install.packages('", pkg, "')"),
           call. = FALSE)
    }
  }

  library(sf)
  library(dplyr)

  cat("============================================================\n")
  cat("  RICHNESS BY IRREGULAR POLYGONS FROM EXTRAPOLATION\n")
  cat("============================================================\n\n")

  # ============================================================================
  # 1. VALIDATE AND LOAD INPUTS
  # ============================================================================

  cat("1) Validating inputs...\n")

  # Validate extrapolation_sf
  if (!inherits(extrapolation_sf, "sf")) {
    stop("extrapolation_sf must be an sf object")
  }

  if (nrow(extrapolation_sf) == 0) {
    stop("extrapolation_sf is empty")
  }

  # Extract species information from extrapolation_sf
  # Assuming columns are species names (excluding geometry)
  geom_col <- attr(extrapolation_sf, "sf_column")
  species_cols <- setdiff(names(extrapolation_sf), geom_col)
  
  if (length(species_cols) == 0) {
    stop("No species columns found in extrapolation_sf")
  }

  cat(paste0("   - Extrapolation geometries: ", nrow(extrapolation_sf), "\n"))
  cat(paste0("   - Species columns found: ", paste(species_cols, collapse = ", "), "\n\n"))

  # ============================================================================
  # 2. LOAD AND VALIDATE BINS SHAPEFILE
  # ============================================================================

  cat("2) Loading and validating bins shapefile...\n")

  # Universal shapefile converter
  if (is.character(bins_shapefile)) {
    bins_sf <- sf::st_read(bins_shapefile, quiet = TRUE)
    cat(paste0("   - Loaded from file: ", bins_shapefile, "\n"))
  } else if (inherits(bins_shapefile, "sf")) {
    bins_sf <- bins_shapefile
    cat("   - Using provided sf object\n")
  } else if (inherits(bins_shapefile, "SpatVector")) {
    bins_sf <- sf::st_as_sf(bins_shapefile)
    cat("   - Converted from terra::SpatVector to sf\n")
  } else if (inherits(bins_shapefile, "Spatial")) {
    bins_sf <- sf::st_as_sf(bins_shapefile)
    cat("   - Converted from sp::Spatial* to sf\n")
  } else {
    stop("bins_shapefile type not recognized")
  }

  # Check if bin_id_column exists
  if (!bin_id_column %in% names(bins_sf)) {
    stop(paste0("Column '", bin_id_column, "' not found in bins_shapefile"))
  }

  cat(paste0("   - Bin ID column: ", bin_id_column, "\n"))
  cat(paste0("   - Number of bins: ", nrow(bins_sf), "\n\n"))

  # ============================================================================
  # 3. TRANSFORM CRS TO MATCH BINS
  # ============================================================================

  cat("3) Transforming CRS...\n")

  extrapolation_transformed <- sf::st_transform(extrapolation_sf, sf::st_crs(bins_sf))
  cat(paste0("   - Target CRS: ", sf::st_crs(bins_sf)$input, "\n\n"))

  # ============================================================================
  # 4. INTERSECT EXTRAPOLATION WITH BINS
  # ============================================================================

  cat("4) Intersecting extrapolation with irregular polygons...\n")

  # For each bin, check which extrapolation cells intersect it
  bins_richness_list <- list()
  species_per_bin_list <- list()

  for (i in seq_len(nrow(bins_sf))) {
    bin_id <- bins_sf[[bin_id_column]][i]
    bin_geom <- bins_sf[i, ]

    # Find intersections with this bin
    intersects <- sf::st_intersects(extrapolation_transformed, bin_geom, sparse = FALSE)
    intersecting_cells <- which(intersects[, 1])

    if (length(intersecting_cells) > 0) {
      # Get the species data for intersecting cells
      intersecting_data <- extrapolation_transformed[intersecting_cells, species_cols, drop = FALSE]

      # Calculate richness (number of species with presence = 1)
      richness <- 0
      species_present <- character(0)

      for (sp in species_cols) {
        # Check if any cell has presence (1) for this species
        if (any(intersecting_data[[sp]] > 0, na.rm = TRUE)) {
          richness <- richness + 1
          species_present <- c(species_present, sp)
        }
      }

      bins_richness_list[[as.character(bin_id)]] <- data.frame(
        bin_id = bin_id,
        richness = richness,
        stringsAsFactors = FALSE
      )

      species_per_bin_list[[as.character(bin_id)]] <- data.frame(
        bin_id = bin_id,
        species = species_present,
        stringsAsFactors = FALSE
      )
    } else {
      # No intersections - richness = 0
      bins_richness_list[[as.character(bin_id)]] <- data.frame(
        bin_id = bin_id,
        richness = 0,
        stringsAsFactors = FALSE
      )

      species_per_bin_list[[as.character(bin_id)]] <- data.frame(
        bin_id = bin_id,
        species = character(0),
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine results
  bins_richness <- do.call(rbind, bins_richness_list)
  rownames(bins_richness) <- NULL

  species_per_bin <- do.call(rbind, species_per_bin_list)
  rownames(species_per_bin) <- NULL

  cat(paste0("   - Bins with species: ", sum(bins_richness$richness > 0), "\n"))
  cat(paste0("   - Mean richness: ", round(mean(bins_richness$richness), 2), "\n"))
  cat(paste0("   - Max richness: ", max(bins_richness$richness), "\n"))

  # ============================================================================
  # 5. CREATE PRESENCE-ABSENCE MATRIX
  # ============================================================================

  # Create a presence-absence matrix with bins as rows and species as columns
  all_species <- sort(unique(species_per_bin$species))
  all_bins <- sort(unique(bins_richness$bin_id))

  pres_abs_matrix <- matrix(0,
                            nrow = length(all_bins),
                            ncol = length(all_species),
                            dimnames = list(as.character(all_bins), all_species))

  # Fill in the matrix
  for (i in seq_len(nrow(species_per_bin))) {
    bin_id <- as.character(species_per_bin$bin_id[i])
    species <- species_per_bin$species[i]
    if (bin_id %in% rownames(pres_abs_matrix) && species %in% colnames(pres_abs_matrix)) {
      pres_abs_matrix[bin_id, species] <- 1
    }
  }

  cat("\n")
  cat("6) Creating presence-absence matrix...\n")
  cat(paste0("   - Matrix dimensions: ", nrow(pres_abs_matrix), " bins x ", ncol(pres_abs_matrix), " species\n"))
  cat(paste0("   - Total presences: ", sum(pres_abs_matrix), "\n\n"))

  # ============================================================================
  # 6. CREATE SF OBJECT WITH RICHNESS FOR MAPPING
  # ============================================================================

  cat("7) Creating sf object for mapping...\n")

  # Join richness data back to bins_sf
  bins_sf_with_richness <- bins_sf
  bins_sf_with_richness$n_species <- 0  # Initialize
  bins_sf_with_richness$richness <- 0   # Also add richness column
  bins_sf_with_richness$species_list <- ""  # For popup

  # Fill in richness values
  for (i in seq_len(nrow(bins_richness))) {
    bin_id <- bins_richness$bin_id[i]
    bin_idx <- which(bins_sf_with_richness[[bin_id_column]] == bin_id)
    if (length(bin_idx) > 0) {
      bins_sf_with_richness$n_species[bin_idx] <- bins_richness$richness[i]
      bins_sf_with_richness$richness[bin_idx] <- bins_richness$richness[i]
    }
  }

  # Add species list for popup
  for (bin_id in unique(species_per_bin$bin_id)) {
    species_in_bin <- species_per_bin$species[species_per_bin$bin_id == bin_id]
    species_str <- paste(species_in_bin, collapse = ", ")
    bin_idx <- which(bins_sf_with_richness[[bin_id_column]] == bin_id)
    if (length(bin_idx) > 0) {
      bins_sf_with_richness$species_list[bin_idx] <- species_str
    }
  }

  cat(paste0("   - SF object created with ", nrow(bins_sf_with_richness), " polygons\n\n"))

  # ============================================================================
  # 7. RETURN RESULTS
  # ============================================================================

  return(list(
    bins_richness = bins_richness,
    species_per_bin = species_per_bin,
    pres_abs = pres_abs_matrix,
    bins_richness_sf = bins_sf_with_richness
  ))
}
