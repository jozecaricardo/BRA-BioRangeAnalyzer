#' Calculate Irregular Bins Richness from Extrapolation Shapefiles
#'
#' This function loads extrapolation shapefiles and calculates species richness
#' for irregular polygons by intersecting the extrapolations with the polygons.
#'
#' @param extrap_method Character indicating extrapolation method ("buffer", "convex_hull", "mst")
#' @param output_dir Character path to the output directory containing extrapolation shapefiles
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
calcRange_irregular_bins_from_extrap_shapefiles <- function(
    extrap_method,
    output_dir,
    bins_shapefile,
    bin_id_column,
    crs_input = 4326) {

  # ============================================================================
  # 1. VALIDATE INPUTS
  # ============================================================================

  if (!dir.exists(output_dir)) {
    stop(paste0("Output directory not found: ", output_dir))
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
  cat("Calculating irregular bins richness from extrapolation shapefiles...\n")
  cat("=" %+% rep("=", 60) %+% "\n\n")

  # ============================================================================
  # 2. FIND AND LOAD EXTRAPOLATION SHAPEFILES
  # ============================================================================

  cat("1) Finding extrapolation shapefiles...\n")

  # Look for taxon-specific extrapolation shapefiles
  # Pattern: GRIDS_taxon_[SPECIES]_[METHOD]_q[RES].shp or
  #          [METHOD]_[SPECIES].shp or similar
  pattern <- paste0(".*", toupper(extrap_method), ".*\\.shp$")
  shapefiles <- list.files(output_dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)

  if (length(shapefiles) == 0) {
    stop(paste0("No extrapolation shapefiles found for method '", extrap_method, "' in ", output_dir))
  }

  cat(paste0("   - Found ", length(shapefiles), " extrapolation shapefile(s)\n"))

  # Load all shapefiles and combine
  all_extrap_sf <- NULL
  species_list <- character(0)

  for (shp_file in shapefiles) {
    tryCatch({
      extrap_sf <- sf::st_read(shp_file, quiet = TRUE)

      # Extract species name from filename
      filename <- basename(shp_file)
      # Try to extract species from patterns like "GRIDS_taxon_[SPECIES]_" or "[SPECIES]_"
      species_match <- regmatches(filename, regexpr("(?<=taxon_)[^_]+", filename, perl = TRUE))
      if (length(species_match) == 0) {
        species_match <- regmatches(filename, regexpr("^[^_]+", filename))
      }
      species_name <- if (length(species_match) > 0) species_match[1] else "unknown"

      # Add species column if not present
      if (!"species" %in% names(extrap_sf) && !"spp" %in% names(extrap_sf)) {
        extrap_sf$species <- species_name
      }

      species_list <- c(species_list, species_name)

      # Combine with other shapefiles
      if (is.null(all_extrap_sf)) {
        all_extrap_sf <- extrap_sf
      } else {
        # Ensure same columns
        common_cols <- intersect(names(all_extrap_sf), names(extrap_sf))
        if ("species" %in% common_cols || "spp" %in% common_cols) {
          all_extrap_sf <- rbind(all_extrap_sf[, common_cols], extrap_sf[, common_cols])
        }
      }
    }, error = function(e) {
      cat(paste0("   - Warning: Could not load ", shp_file, ": ", e$message, "\n"))
    })
  }

  if (is.null(all_extrap_sf) || nrow(all_extrap_sf) == 0) {
    stop("Could not load any valid extrapolation shapefiles")
  }

  cat(paste0("   - Loaded ", nrow(all_extrap_sf), " geometries from ", length(unique(species_list)), " species\n\n"))

  # ============================================================================
  # 3. ENSURE SAME CRS
  # ============================================================================

  cat("2) Checking coordinate reference systems...\n")

  if (!identical(sf::st_crs(all_extrap_sf), sf::st_crs(bins_sf))) {
    cat(paste0("   - Transforming extrapolation CRS to ", sf::st_crs(bins_sf)$input, "\n"))
    all_extrap_sf <- sf::st_transform(all_extrap_sf, sf::st_crs(bins_sf))
  }

  cat("   - CRS aligned\n\n")

  # ============================================================================
  # 4. INTERSECT EXTRAPOLATION WITH BINS
  # ============================================================================

  cat("3) Intersecting extrapolation with irregular polygons...\n")

  # Determine species column
  species_col <- NA
  if ("species" %in% names(all_extrap_sf)) {
    species_col <- "species"
  } else if ("spp" %in% names(all_extrap_sf)) {
    species_col <- "spp"
  }

  if (is.na(species_col)) {
    stop("Could not find species column in extrapolation shapefiles")
  }

  bins_richness_list <- list()
  species_per_bin_list <- list()

  for (i in seq_len(nrow(bins_sf))) {
    bin_id <- bins_sf[[bin_id_column]][i]
    bin_geom <- bins_sf[i, ]

    # Find extrapolation geometries that intersect this bin
    intersecting <- sf::st_intersects(all_extrap_sf, bin_geom, sparse = FALSE)[, 1]
    species_in_bin <- unique(all_extrap_sf[[species_col]][intersecting])
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
