#' Calculate Irregular Bins Richness from Extrapolation Geometries
#'
#' This function calculates species richness for irregular polygons based on
#' extrapolation geometries (Buffer, Convex Hull, MST) rather than just occurrence points.
#'
#' @param extrap_sf sf object with extrapolation geometries (e.g., grid cells from Buffer/Convex Hull/MST)
#' @param bins_shapefile SpatVector or sf object with irregular polygons
#' @param bin_id_column Character name of the column containing polygon IDs
#' @param crs_input Numeric EPSG code for input CRS (default: 4326)
#'
#' @return List containing:
#'   - bins_richness: Data frame with bin_id and richness (n_species)
#'   - species_per_bin: Data frame with bin_id and species names
#'   - pres_abs: Presence-absence matrix (bins x species)
#'
#' @keywords internal
#'
calcRange_irregular_bins_from_extrap <- function(
    extrap_sf,
    bins_shapefile,
    bin_id_column,
    crs_input = 4326) {

  # ============================================================================
  # 1. VALIDATE INPUTS
  # ============================================================================

  if (!inherits(extrap_sf, "sf")) {
    stop("extrap_sf must be an sf object")
  }

  if (nrow(extrap_sf) == 0) {
    stop("extrap_sf has no geometries")
  }

  # Convert bins_shapefile to sf if needed
  if (inherits(bins_shapefile, "SpatVector")) {
    bins_sf <- sf::st_as_sf(bins_shapefile)
  } else if (inherits(bins_shapefile, "sf")) {
    bins_sf <- bins_shapefile
  } else {
    stop("bins_shapefile must be SpatVector or sf object")
  }

  if (!(bin_id_column %in% names(bins_sf))) {
    stop(paste0("Column '", bin_id_column, "' not found in bins_shapefile"))
  }

  cat("\n")
  cat("Calculating irregular bins richness from extrapolation geometries...\n")
  cat("=" %+% rep("=", 60) %+% "\n\n")

  # ============================================================================
  # 2. ENSURE SAME CRS
  # ============================================================================

  cat("1) Checking coordinate reference systems...\n")

  # Transform extrap_sf to match bins_sf CRS
  if (!identical(sf::st_crs(extrap_sf), sf::st_crs(bins_sf))) {
    cat(paste0("   - Transforming extrapolation CRS from ", sf::st_crs(extrap_sf)$input, " to ", sf::st_crs(bins_sf)$input, "\n"))
    extrap_sf <- sf::st_transform(extrap_sf, sf::st_crs(bins_sf))
  }

  cat("   - CRS aligned\n\n")

  # ============================================================================
  # 3. EXTRACT SPECIES FROM EXTRAPOLATION GEOMETRIES
  # ============================================================================

  cat("2) Extracting species information from extrapolation geometries...\n")

  # The extrap_sf should have a 'species' or 'spp' column
  species_col <- NA
  if ("species" %in% names(extrap_sf)) {
    species_col <- "species"
  } else if ("spp" %in% names(extrap_sf)) {
    species_col <- "spp"
  } else if ("taxon" %in% names(extrap_sf)) {
    species_col <- "taxon"
  } else {
    # Try to find any column that looks like species
    possible_cols <- names(extrap_sf)[!names(extrap_sf) %in% c("geometry", "geom")]
    if (length(possible_cols) > 0) {
      species_col <- possible_cols[1]
      cat(paste0("   - Warning: Using column '", species_col, "' as species column\n"))
    }
  }

  if (is.na(species_col)) {
    stop("Could not find species column in extrap_sf. Expected 'species', 'spp', or 'taxon'")
  }

  all_species <- unique(extrap_sf[[species_col]])
  cat(paste0("   - Found ", length(all_species), " species in extrapolation\n\n"))

  # ============================================================================
  # 4. INTERSECT EXTRAPOLATION WITH BINS
  # ============================================================================

  cat("3) Intersecting extrapolation with irregular polygons...\n")

  bins_richness_list <- list()
  species_per_bin_list <- list()

  for (i in seq_len(nrow(bins_sf))) {
    bin_id <- bins_sf[[bin_id_column]][i]
    bin_geom <- bins_sf[i, ]

    # Find extrapolation cells that intersect this bin
    intersecting <- sf::st_intersects(extrap_sf, bin_geom, sparse = FALSE)[, 1]
    species_in_bin <- unique(extrap_sf[[species_col]][intersecting])
    species_in_bin <- species_in_bin[!is.na(species_in_bin)]

    if (length(species_in_bin) > 0) {
      for (sp in species_in_bin) {
        species_per_bin_list[[length(species_per_bin_list) + 1]] <- data.frame(
          bin_id = bin_id,
          species = sp,
          stringsAsFactors = FALSE
        )
      }

      bins_richness_list[[length(bins_richness_list) + 1]] <- data.frame(
        bin_id = bin_id,
        richness = length(species_in_bin),
        n_species = length(species_in_bin),
        stringsAsFactors = FALSE
      )
    } else {
      bins_richness_list[[length(bins_richness_list) + 1]] <- data.frame(
        bin_id = bin_id,
        richness = 0,
        n_species = 0,
        stringsAsFactors = FALSE
      )
    }
  }

  bins_richness <- do.call(rbind, bins_richness_list)
  rownames(bins_richness) <- NULL

  species_per_bin <- do.call(rbind, species_per_bin_list)
  rownames(species_per_bin) <- NULL

  cat(paste0("   - Bins with species: ", sum(bins_richness$richness > 0), "\n"))
  cat(paste0("   - Mean richness: ", round(mean(bins_richness$richness), 2), "\n"))
  cat(paste0("   - Max richness: ", max(bins_richness$richness), "\n\n"))

  # ============================================================================
  # 5. CREATE PRESENCE-ABSENCE MATRIX
  # ============================================================================

  cat("4) Creating presence-absence matrix...\n")

  all_species_sorted <- sort(unique(species_per_bin$species))
  all_bins_sorted <- sort(unique(bins_richness$bin_id))

  pres_abs_matrix <- matrix(0,
                            nrow = length(all_bins_sorted),
                            ncol = length(all_species_sorted),
                            dimnames = list(as.character(all_bins_sorted), all_species_sorted))

  for (i in seq_len(nrow(species_per_bin))) {
    bin_id <- as.character(species_per_bin$bin_id[i])
    species <- species_per_bin$species[i]
    if (bin_id %in% rownames(pres_abs_matrix) && species %in% colnames(pres_abs_matrix)) {
      pres_abs_matrix[bin_id, species] <- 1
    }
  }

  cat(paste0("   - Matrix dimensions: ", nrow(pres_abs_matrix), " bins x ", ncol(pres_abs_matrix), " species\n"))
  cat(paste0("   - Total presences: ", sum(pres_abs_matrix), "\n\n"))

  # ============================================================================
  # 6. RETURN RESULTS
  # ============================================================================

  return(list(
    bins_richness = bins_richness,
    species_per_bin = species_per_bin,
    pres_abs = pres_abs_matrix
  ))
}

# Helper operator for string concatenation
`%+%` <- function(x, y) paste0(x, y)
