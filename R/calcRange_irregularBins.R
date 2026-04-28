#' calcRange_irregularBins
#'
#' Extrapolates species distributions to irregular spatial units (biogeographic regions,
#' biomes, provinces, ecoregions, etc.) using spatial intersection. Creates a
#' presence-absence matrix and calculates species richness per bin.
#'
#' @param xy Data.frame or list containing species occurrence data with columns:
#'   - `spp`: species names
#'   - `Long` or `long`: longitude coordinates
#'   - `Lat` or `lat`: latitude coordinates
#' @param bins_shapefile Shapefile containing irregular bins. Accepts any of:
#'   - File path (character): "path/to/shapefile.shp"
#'   - sf object: from st_read()
#'   - terra::SpatVector: from vect()
#'   - sp::Spatial*: from readOGR() or shapefile()
#'   Must have a column identifying each bin (specified in `bin_id_column`).
#' @param bin_id_column Character. Name of the column in `bins_shapefile` that
#'   identifies each bin (e.g., "Dominio", "Biome", "Province"). Default: "bin_id"
#' @param resol Numeric vector of length 2. Resolution of grid cells in degrees
#'   (longitude, latitude) for generating grid shapefiles. Default: c(1, 1) (1 degrees x 1 degrees).
#'   Examples: c(0.5, 0.5) for 0.5 degrees cells, c(2, 2) for 2 degrees cells.
#' @param crs_input Integer or character. CRS of input coordinates. Default: 4326 (WGS84)
#' @param output_dir Character. Directory for output files. Default: "out_irregular_bins/"
#'
#' @return A list with four elements:
#'   - `bins_richness`: sf object with species richness per bin
#'   - `pres_abs`: presence-absence matrix (bins x species) with NUMERIC rownames
#'   - `species_per_bin`: data.frame with species count per bin
#'   - `bin_id_mapping`: data.frame mapping bin names to numeric IDs
#'
#' @details
#' This function performs the following steps:
#' 1. Converts occurrence points to sf object
#' 2. Validates and fixes invalid geometries in bins shapefile
#' 3. Transforms coordinates to match bins CRS
#' 4. Performs spatial join (st_join) to assign occurrences to bins
#' 5. Creates presence-absence matrix (bins x species)
#' 6. Calculates species richness per bin
#' 7. Saves outputs to specified directory
#'
#' **Input Validation:**
#' - Automatically fixes invalid geometries using `st_make_valid()`
#' - Handles different input formats (data.frame or list)
#' - Transforms CRS to match bins shapefile
#'
#' **Output Files:**
#' - `pres_abs_irregular_bins.txt`: presence-absence matrix (numeric rownames)
#' - `species_richness_bins.shp`: shapefile with richness per bin
#' - `species_per_bin.csv`: table with species count per bin
#' - `bin_id_mapping.csv`: mapping between bin names and numeric IDs
#'
#' @examples
#' \dontrun{
#' # Example with South American countries
#' library(sf)
#' library(rnaturalearthdata)
#'
#' # Load occurrence data
#' data <- data.frame(
#'   spp = c("sp1", "sp1", "sp2", "sp2", "sp3"),
#'   Long = c(-50, -51, -48, -49, -52),
#'   Lat = c(-10, -11, -12, -13, -14)
#' )
#'
#' # Get South America countries and Brazilian states as separate polygons
#' sa_countries <- ne_countries(continent = "South America", scale = "medium",
#'   returnclass = "sf")
#'
#' brazil_states <- ne_states(country = "Brazil", returnclass = "sf")
#'
#' # removing Brazil
#' sa_without_brazil <- sa_countries[sa_countries$iso_a3 != "BRA", ]
#'
#' # standardizing columns
#' countries <- sa_no_br[, c("name", "iso_a3")]
#' names(countries) <- c("region_name", "iso_code", "geometry")
#' countries$type <- "country"
#'
#' states <- brazil_states[, c("name", "iso_3166_2")]
#' names(states) <- c("region_name", "iso_code", "geometry")
#'
#' states$type <- "state"
#'
#' # validating geometries
#' countries_clean <- st_make_valid(countries)
#' states_clean <- st_make_valid(states)
#'
#' # Combining polygons
#' sa_bins <- rbind(countries_clean, states_clean)
#'
#' # Convert to terra SpatVector if needed
#' sa_bins_sv <- vect(sa_bins)
#'
#' # Calculate distribution using irregular bins
#' result <- calcRange_irregular_bins(
#'   xy = data,
#'   bins_shapefile = sa_bins,
#'   bin_id_column = "region_name",
#'   resol = c(1, 1)  # 1 degrees x 1 degrees grid cells
#' )
#'
#' # Access results
#' result$bins_richness        # sf object with richness
#' result$pres_abs             # presence-absence matrix (numeric rownames)
#' result$species_per_bin      # species count per bin
#' result$bin_id_mapping       # bin names <-> numeric IDs mapping
#'
#' # Plot richness map
#' library(ggplot2)
#' ggplot(result$bins_richness) +
#'   geom_sf(aes(fill = n_species)) +
#'   scale_fill_viridis_c() +
#'   theme_minimal()
#' }
#'
#' @export
calcRange_irregular_bins <- function(xy,
                                     bins_shapefile,
                                     bin_id_column = "bin_id",
                                     resol = c(1, 1),
                                     crs_input = 4326,
                                     output_dir = "out_irregular_bins/") {

  # Check for required packages
  required_packages <- c("sf", "dplyr", "lwgeom")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' not found. Please install: install.packages('", pkg, "')"),
           call. = FALSE)
    }
  }

  library(sf)
  library(dplyr)
  library(lwgeom)

  cat("============================================================\n")
  cat("  DISTRIBUTION EXTRAPOLATION USING IRREGULAR BINS\n")
  cat("============================================================\n\n")

  # ============================================================================
  # 1. PREPARE INPUT DATA
  # ============================================================================

  cat("1) Preparing occurrence data...\n")

  # Handle different input data types
  if (is.data.frame(xy)) {
    names(xy) <- tolower(names(xy))
    dat <- xy[, c("spp", "long", "lat")]
  } else if (is.list(xy) && "samples" %in% names(xy)) {
    dat <- xy$samples[, 1:3]
    names(dat) <- c("spp", "long", "lat")
  } else {
    stop("Input 'xy' must be a data.frame with columns 'spp', 'Long', 'Lat' or a list with 'samples'")
  }

  # Remove NA coordinates
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) == 0) {
    stop("No valid occurrence records found after removing NAs")
  }

  cat(paste0("   - Total occurrences: ", nrow(dat), "\n"))
  cat(paste0("   - Total species: ", length(unique(dat$spp)), "\n\n"))

  # ============================================================================
  # 2. CONVERT TO SF OBJECT
  # ============================================================================

  cat("2) Converting occurrences to sf object...\n")

  occurrences_sf <- st_as_sf(dat,
                             coords = c('long', 'lat'),
                             crs = crs_input)

  cat(paste0("   - CRS: ", st_crs(occurrences_sf)$input, "\n\n"))

  # ============================================================================
  # 3. LOAD AND VALIDATE BINS SHAPEFILE
  # ============================================================================

  cat("3) Loading and validating bins shapefile...\n")

  # Universal shapefile converter: accepts any type
  if (is.character(bins_shapefile)) {
    # File path provided
    bins_sf <- st_read(bins_shapefile, quiet = TRUE)
    cat(paste0("   - Loaded from file: ", bins_shapefile, "\n"))
  } else if (inherits(bins_shapefile, "sf")) {
    # Already sf object
    bins_sf <- bins_shapefile
    cat("   - Using provided sf object\n")
  } else if (inherits(bins_shapefile, "SpatVector")) {
    # terra::vect() object
    bins_sf <- st_as_sf(bins_shapefile)
    cat("   - Converted from terra::SpatVector to sf\n")
  } else if (inherits(bins_shapefile, "Spatial")) {
    # sp package object (SpatialPolygons*, etc.)
    bins_sf <- st_as_sf(bins_shapefile)
    cat("   - Converted from sp::Spatial* to sf\n")
  } else {
    stop(paste0("bins_shapefile type not recognized. ",
                "Accepted types: file path (character), sf, terra::SpatVector, or sp::Spatial*. ",
                "Provided type: ", class(bins_shapefile)[1]))
  }

  # Check if bin_id_column exists
  if (!bin_id_column %in% names(bins_sf)) {
    stop(paste0("Column '", bin_id_column, "' not found in bins_shapefile. ",
                "Available columns: ", paste(names(bins_sf), collapse = ", ")))
  }

  cat(paste0("   - Bin ID column: ", bin_id_column, "\n"))
  cat(paste0("   - Number of bins: ", nrow(bins_sf), "\n"))

  # Check for invalid geometries (safe mode)
  valid_probe <- st_is_valid(bins_sf, NA_on_exception = TRUE)
  invalid_geoms <- is.na(valid_probe) | !valid_probe

  if (any(invalid_geoms)) {
    cat(paste0("   - WARNING: Found ", sum(invalid_geoms), " invalid/problematic geometries\n"))
    cat("   - Fixing invalid geometries with robust sequence (st_make_valid -> lwgeom -> st_buffer(0))...\n")

    bins_sf <- tryCatch(st_make_valid(bins_sf), error = function(e) bins_sf)

    valid_after_make <- st_is_valid(bins_sf, NA_on_exception = TRUE)
    still_bad <- is.na(valid_after_make) | !valid_after_make

    if (any(still_bad) && requireNamespace("lwgeom", quietly = TRUE)) {
      repaired <- tryCatch(lwgeom::st_make_valid(bins_sf[still_bad, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(repaired)) {
        bins_sf[still_bad, ] <- repaired
      }
    }

    valid_after_lw <- st_is_valid(bins_sf, NA_on_exception = TRUE)
    still_bad <- is.na(valid_after_lw) | !valid_after_lw
    if (any(still_bad)) {
      buffered <- tryCatch(st_buffer(bins_sf[still_bad, , drop = FALSE], 0), error = function(e) NULL)
      if (!is.null(buffered)) {
        bins_sf[still_bad, ] <- buffered
      }
    }

    valid_final <- st_is_valid(bins_sf, NA_on_exception = TRUE)
    valid_final[is.na(valid_final)] <- FALSE
    if (!all(valid_final)) {
      warning("Some geometries could not be fixed. Filtering them out...")
      bins_sf <- bins_sf[valid_final, , drop = FALSE]
    }

    if (nrow(bins_sf) == 0) {
      stop("No valid geometries remained in bins shapefile after repair.")
    }
    cat("   - [OK] Geometry validation step finished\n")
  } else {
    cat("   - [OK] All geometries are valid\n")
  }

  cat("\n")

  # ============================================================================
  # 4. TRANSFORM CRS TO MATCH BINS
  # ============================================================================

  cat("4) Transforming coordinates to match bins CRS...\n")

  occurrences_transformed <- st_transform(occurrences_sf, st_crs(bins_sf))

  cat(paste0("   - Target CRS: ", st_crs(bins_sf)$input, "\n\n"))

  # ============================================================================
  # 5. SPATIAL JOIN (ASSIGN OCCURRENCES TO BINS)
  # ============================================================================

  cat("5) Performing spatial join (occurrences -> bins)...\n")

  intersections <- st_intersects(occurrences_transformed, bins_sf)
  n_assigned <- sum(lengths(intersections) > 0)
  n_outside <- nrow(occurrences_transformed) - n_assigned

  occurrences_with_bins <- st_join(occurrences_transformed, bins_sf)

  # Remove occurrences that didn't fall in any bin
  occurrences_with_bins <- occurrences_with_bins[!is.na(occurrences_with_bins[[bin_id_column]]), ]

  cat(paste0("   - Occurrences assigned to bins: ", n_assigned, "\n"))
  cat(paste0("   - Occurrences outside bins: ", n_outside, "\n\n"))

  if (nrow(occurrences_with_bins) == 0) {
    stop("No occurrences fell within any bin. Check CRS and spatial overlap.")
  }

  # ============================================================================
  # 6. CREATE PRESENCE-ABSENCE MATRIX
  # ============================================================================

  cat("6) Creating presence-absence matrix...\n")

  # Get unique species and bins
  all_species <- sort(unique(dat$spp))
  all_bins <- sort(unique(bins_sf[[bin_id_column]]))

  # Create mapping between bin names and numeric IDs
  bin_id_mapping <- data.frame(
    bin_name = all_bins,
    bin_number = 1:length(all_bins),
    stringsAsFactors = FALSE
  )

  # Initialize matrix with zeros (using AREA NAMES as rownames)
  pres_abs_matrix <- matrix(0,
                            nrow = length(all_bins),
                            ncol = length(all_species),
                            dimnames = list(all_bins, all_species))

  # Convert to data.frame for easier manipulation
  occ_df <- as.data.frame(occurrences_with_bins)
  occ_df <- occ_df[, c("spp", bin_id_column)]

  # Fill presence-absence matrix
  for (i in 1:nrow(occ_df)) {
    species <- as.character(occ_df$spp[i])
    bin_name <- as.character(occ_df[[bin_id_column]][i])

    if (bin_name %in% all_bins && species %in% all_species) {
      # Use bin_name directly as rowname
      pres_abs_matrix[bin_name, species] <- 1
    }
  }

  # Add ROOT row (for TNT compatibility)
  pres_abs_matrix <- rbind(pres_abs_matrix, ROOT = rep(0, ncol(pres_abs_matrix)))

  cat(paste0("   - Matrix dimensions: ", nrow(pres_abs_matrix) - 1, " bins x ",
             ncol(pres_abs_matrix), " species\n"))
  cat(paste0("   - Total presences: ", sum(pres_abs_matrix), "\n\n"))

  # ============================================================================
  # 7. CALCULATE SPECIES RICHNESS PER BIN
  # ============================================================================

  cat("7) Calculating species richness per bin...\n")

  # Count unique species per bin
  species_per_bin <- occurrences_with_bins %>%
    st_drop_geometry() %>%
    dplyr::group_by(!!rlang::sym(bin_id_column)) %>%
    dplyr::summarise(n_species = dplyr::n_distinct(spp),
                     n_occurrences = dplyr::n(),
                     species_list = paste(unique(spp), collapse = ", "),
                     .groups = "drop")

  cat(paste0("   - Bins with species: ", nrow(species_per_bin), "\n"))
  cat(paste0("   - Mean richness: ", round(mean(species_per_bin$n_species), 2), "\n"))
  cat(paste0("   - Max richness: ", max(species_per_bin$n_species), "\n\n"))

  # ============================================================================
  # 8. JOIN RICHNESS TO BINS SHAPEFILE
  # ============================================================================

  cat("8) Joining richness data to bins shapefile...\n")

  bins_with_richness <- bins_sf %>%
    left_join(species_per_bin, by = bin_id_column) %>%
    mutate(n_species = ifelse(is.na(n_species), 0, n_species),
           n_occurrences = ifelse(is.na(n_occurrences), 0, n_occurrences))

  cat("   - [OK] Richness data joined successfully\n\n")

  # ============================================================================
  # 9. SAVE OUTPUTS
  # ============================================================================

  cat("9) Saving outputs...\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("   - Created directory: ", output_dir, "\n"))
  }

  # Save presence-absence matrix
  write.table(pres_abs_matrix,
              file = file.path(output_dir, "pres_abs_irregular_bins.txt"),
              sep = "\t",
              quote = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "pres_abs_irregular_bins.txt"), "\n"))

  # Save species per bin table
  write.csv(species_per_bin,
            file = file.path(output_dir, "species_per_bin.csv"),
            row.names = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "species_per_bin.csv"), "\n"))

  # Save bin ID mapping (bin names <-> numeric IDs)
  write.csv(bin_id_mapping,
            file = file.path(output_dir, "bin_id_mapping.csv"),
            row.names = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "bin_id_mapping.csv"), "\n"))

  # Save shapefile with richness
  richness_shp <- file.path(output_dir, "species_richness_bins.shp")
  st_write(bins_with_richness,
           richness_shp,
           delete_layer = file.exists(richness_shp),
           quiet = TRUE)
  cat(paste0("   - Saved: ", richness_shp, "\n"))

  # Save individual species shapefiles
  cat("\n   - Saving individual species shapefiles...\n")

  for (species in all_species) {
    species_occ <- occurrences_with_bins[occurrences_with_bins$spp == species, ]

    if (nrow(species_occ) > 0) {
      occ_shp <- file.path(output_dir, paste0("pointshape_", species, ".shp"))
      st_write(species_occ,
               occ_shp,
               delete_layer = file.exists(occ_shp),
               quiet = TRUE)

      species_bins <- unique(species_occ[[bin_id_column]])
      species_bins_sf <- bins_sf[bins_sf[[bin_id_column]] %in% species_bins, ]

      bins_shp <- file.path(output_dir, paste0("BINS_", species, ".shp"))
      st_write(species_bins_sf,
               bins_shp,
               delete_layer = file.exists(bins_shp),
               quiet = TRUE)
    }
  }

  cat(paste0("   - Saved ", length(all_species), " species shapefiles\n"))
  cat("     * pointshape_[species].shp (occurrence points)\n")
  cat("     * BINS_[species].shp (bins where species occurs)\n\n")

  # Save grid shapefiles for each bin
  cat("   - Saving grid shapefiles for each bin...\n")

  # Check if raster package is available
  if (!requireNamespace("raster", quietly = TRUE)) {
    cat("     WARNING: Package 'raster' not found. Skipping grid shapefile generation.\n")
  } else {
    library(raster)

    # Create a grid based on bins extent with user-defined resolution
    bins_extent <- extent(as(bins_sf, "Spatial"))
    grid_raster <- raster(bins_extent, resolution = resol,
                          crs = CRS("+proj=longlat +datum=WGS84"))
    grid_polygons <- rasterToPolygons(grid_raster)
    grid_polygons_sf <- st_as_sf(grid_polygons)
    grid_polygons_sf$grid_id <- 1:nrow(grid_polygons_sf)

    # Spatial join: which grid cell belongs to which bin?
    grid_bin_join <- tryCatch({
      st_join(grid_polygons_sf, bins_sf[, c(bin_id_column, "geometry")],
              join = st_intersects, largest = TRUE)
    }, error = function(e) {
      cat(paste0("     WARNING: Grid-bin join failed with default geometry engine: ",
                 e$message, "\n"))
      cat("     - Retrying with validated geometries and s2 disabled...\n")

      tryCatch({
        old_s2 <- sf::sf_use_s2()
        on.exit(sf::sf_use_s2(old_s2), add = TRUE)
        sf::sf_use_s2(FALSE)

        bins_for_join <- bins_sf[, c(bin_id_column, "geometry")]
        bins_for_join <- tryCatch(st_make_valid(bins_for_join), error = function(e) bins_for_join)
        grid_for_join <- tryCatch(st_make_valid(grid_polygons_sf), error = function(e) grid_polygons_sf)

        st_join(grid_for_join, bins_for_join, join = st_intersects, largest = TRUE)
      }, error = function(e2) {
        warning(paste0("Could not assign grid cells to bins: ", e2$message,
                       ". Skipping grid shapefile generation."))
        NULL
      })
    })

    if (!is.null(grid_bin_join)) {
      # Save grid shapefile for each bin
      for (bin_name in all_bins) {
        bin_grids <- grid_bin_join[!is.na(grid_bin_join[[bin_id_column]]) &
                                     grid_bin_join[[bin_id_column]] == bin_name, ]

        if (nrow(bin_grids) > 0) {
          # Sanitize bin name for filename
          safe_bin_name <- gsub(" ", "_", bin_name)
          safe_bin_name <- gsub("[^A-Za-z0-9_]", "", safe_bin_name)

          grid_shp <- file.path(output_dir, paste0("GRIDS_", safe_bin_name, ".shp"))
          st_write(bin_grids,
                   grid_shp,
                   delete_layer = file.exists(grid_shp),
                   quiet = TRUE)
        }
      }

      cat(paste0("   - Saved grid shapefiles for ", length(all_bins), " bins\n"))
      cat("     * GRIDS_[bin_name].shp (grid cells within each bin)\n\n")
    } else {
      cat("   - WARNING: Grid shapefiles were skipped due to geometry issues\n\n")
    }
  }

  # ============================================================================
  # 10. SUMMARY
  # ============================================================================

  cat("============================================================\n")
  cat("  SUMMARY\n")
  cat("============================================================\n")
  cat(paste0("Total species: ", length(all_species), "\n"))
  cat(paste0("Total bins: ", length(all_bins), "\n"))
  cat(paste0("Bins with species: ", nrow(species_per_bin), "\n"))
  cat(paste0("Empty bins: ", length(all_bins) - nrow(species_per_bin), "\n"))
  cat(paste0("Mean richness per bin: ", round(mean(species_per_bin$n_species), 2), "\n"))
  cat(paste0("Max richness: ", max(species_per_bin$n_species),
             " (", species_per_bin[[bin_id_column]][which.max(species_per_bin$n_species)], ")\n"))
  cat("============================================================\n\n")

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  return(list(
    bins_richness = bins_with_richness,
    pres_abs = pres_abs_matrix,
    species_per_bin = species_per_bin,
    bin_id_mapping = bin_id_mapping
  ))
}
