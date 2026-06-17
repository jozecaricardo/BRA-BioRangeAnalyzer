###############################################################################################
###############################################################################################
#' BioGeoBEARS Range Calculation Methods
#' 
#' These functions calculate species ranges using different extrapolation methods
#' (MST, Convex Hull, Buffer) and format the output specifically for BioGeoBEARS analysis.
#' 
#' @import sf
#' @import dplyr
#' @import geosphere
#' @import fossil

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Create 2-letter area code
create_area_code <- function(area_name) {
  name_clean <- gsub(" ", "", area_name)
  first_letter <- toupper(substr(name_clean, 1, 1))
  
  if (nchar(name_clean) >= 2) {
    second_letter <- toupper(substr(name_clean, 2, 2))
  } else {
    second_letter <- "X"
  }
  
  paste0(first_letter, second_letter)
}

#' Resolve duplicate area codes
resolve_duplicate_codes <- function(area_names) {
  codes <- sapply(area_names, create_area_code, USE.NAMES = FALSE)
  
  df <- data.frame(
    original_name = area_names,
    code = codes,
    index = seq_along(area_names),
    stringsAsFactors = FALSE
  )
  
  duplicates <- duplicated(df$code) | duplicated(df$code, fromLast = TRUE)
  
  if (sum(duplicates) == 0) {
    return(df)
  }
  
  while (any(duplicated(df$code))) {
    for (code in unique(df$code[duplicated(df$code)])) {
      dup_indices <- which(df$code == code)
      
      if (length(dup_indices) > 1) {
        for (j in seq_along(dup_indices)) {
          idx <- dup_indices[j]
          name_clean <- gsub(" ", "", df$original_name[idx])
          
          # Try to use different letters from the name
          if (nchar(name_clean) >= 3) {
            for (letter_pos in 3:nchar(name_clean)) {
              new_code <- paste0(
                toupper(substr(name_clean, 1, 1)),
                toupper(substr(name_clean, letter_pos, letter_pos))
              )
              
              if (!(new_code %in% df$code[-idx])) {
                df$code[idx] <- new_code
                break
              }
            }
          }
          
          # If still duplicate or name too short, append index number
          if (df$code[idx] == code) {
            df$code[idx] <- paste0(substr(df$code[idx], 1, 1), j)
          }
        }
      }
    }
  }
  
  return(df)
}

#' Save matrix in BioGeoBEARS format
save_biogeobears_format <- function(pres_abs_matrix, output_dir, method_name, original_names = NULL, area_mapping = NULL) {
  # Remove empty species
  col_sums <- colSums(pres_abs_matrix)
  pres_abs_clean <- pres_abs_matrix[, col_sums > 0, drop = FALSE]
  
  # Extract info
  n_species <- ncol(pres_abs_clean)
  n_ranges <- nrow(pres_abs_clean)
  species_names <- colnames(pres_abs_clean)
  
  # Create sequential codes (A, B, C, ...) using the same logic as Step 5
  make_area_code <- function(idx) {
    if (!is.finite(idx) || idx < 1) return(NA_character_)
    idx <- as.integer(idx)
    out <- ""
    while (idx > 0) {
      r <- (idx - 1) %% 26
      out <- paste0(LETTERS[r + 1], out)
      idx <- (idx - 1) %/% 26
    }
    out
  }
  
  sequential_codes <- vapply(seq_len(n_ranges), make_area_code, character(1))
  
  # Create mapping dataframe: code => original_name
  if (!is.null(original_names)) {
    # Use provided original names
    mapping_df <- data.frame(
      code = sequential_codes,
      original_name = original_names,
      stringsAsFactors = FALSE
    )
  } else {
    # Fallback if not provided
    mapping_df <- data.frame(
      code = sequential_codes,
      original_name = rownames(pres_abs_clean),
      stringsAsFactors = FALSE
    )
  }
  
  # Create text
  texto <- character()
  header <- paste(n_species, n_ranges, 
                  paste0('(', paste(sequential_codes, collapse = ' '), ')'),
                  sep = ' ')
  texto <- c(texto, header)
  
  for (j in 1:n_species) {
    presence_string <- paste(pres_abs_clean[, j], collapse = '')
    sp_name <- gsub(" ", "_", species_names[j])
    texto <- c(texto, paste(sp_name, presence_string))
  }
  
  # Save file
  output_file <- file.path(output_dir, paste0("pres_abs_", method_name, "_geog.data"))
  writeLines(texto, con = output_file)
  
  # Save mapping file
  mapping_file <- file.path(output_dir, "area_code_mapping.csv")
  write.csv(mapping_df, mapping_file, row.names = FALSE)
  
  return(list(
    file = output_file,
    n_species = n_species,
    n_ranges = n_ranges,
    matrix = pres_abs_clean,
    mapping = mapping_df
  ))
}

# ============================================================================
# 1. MST METHOD
# ============================================================================

calcRange_mst_biogeobears <- function(occurrence_data, irregular_polygons, bin_id_column, crs_input = 4326, output_dir = "out_MST_biogeobears/") {
  
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  names(occurrence_data) <- tolower(names(occurrence_data))
  dat <- occurrence_data[, c("spp", "long", "lat")]
  dat <- dat[complete.cases(dat), ]
  
  occurrences_sf <- sf::st_as_sf(dat, coords = c("long", "lat"), crs = crs_input)
  
  # Shapefile should already be sf (converted in Step 1 for BioGeoBEARS)
  bins_sf <- irregular_polygons
  
  # Check and fix invalid geometries
  invalid_geoms <- !sf::st_is_valid(bins_sf)
  if (any(invalid_geoms)) {
    bins_sf <- sf::st_make_valid(bins_sf)
    
    # If still invalid after st_make_valid, filter them out
    invalid_after <- !sf::st_is_valid(bins_sf)
    if (any(invalid_after)) {
      bins_sf <- bins_sf[!invalid_after, ]
    } else {
    }
  }
  
  occurrences_transformed <- sf::st_transform(occurrences_sf, sf::st_crs(bins_sf))
  
  all_species <- sort(unique(dat$spp))
  all_bins <- sort(unique(bins_sf[[bin_id_column]]))
  
  pres_abs_matrix <- matrix(0, nrow = length(all_bins), ncol = length(all_species), dimnames = list(all_bins, all_species))
  
  for (sp in all_species) {
    sp_occurrences <- occurrences_transformed[occurrences_transformed$spp == sp, ]
    n_points <- nrow(sp_occurrences)
    
    if (n_points < 2) {
      next
    }
    
    # Calculate MST
    sp_coords <- sf::st_coordinates(sp_occurrences)
    distances <- geosphere::distm(sp_coords, fun = geosphere::distHaversine)
    distances_km <- distances / 1000
    
    mst_result <- fossil::dino.mst(distances_km)
    
    if (is.null(mst_result) || !is.matrix(mst_result)) next
    
    mst_edges <- which(mst_result > 0, arr.ind = TRUE)
    if (!is.matrix(mst_edges)) mst_edges <- matrix(mst_edges, nrow = 1)
    if (nrow(mst_edges) == 0) next
    
    mst_edges <- mst_edges[mst_edges[, 1] < mst_edges[, 2], ]
    if (!is.matrix(mst_edges)) mst_edges <- matrix(mst_edges, nrow = 1)
    if (nrow(mst_edges) == 0) next
    
    # Create MST geometry (multilinestring)
    mst_lines <- list()
    for (i in seq_len(nrow(mst_edges))) {
      idx1 <- mst_edges[i, 1]
      idx2 <- mst_edges[i, 2]
      mst_lines[[i]] <- sf::st_linestring(matrix(c(sp_coords[idx1, ], sp_coords[idx2, ]), nrow = 2, byrow = TRUE))
    }
    
    mst_geometry <- sf::st_multilinestring(mst_lines)
    mst_sf <- sf::st_sf(species = sp, geometry = sf::st_sfc(mst_geometry, crs = sf::st_crs(bins_sf)))
    
    # Use st_join to find which bins the MST intersects
    mst_joined <- sf::st_join(bins_sf[, bin_id_column], mst_sf, join = sf::st_intersects)
    bins_with_mst <- unique(mst_joined[[bin_id_column]][!is.na(mst_joined$species)])
    bins_with_mst <- bins_with_mst[!is.na(bins_with_mst)]
    
    # Also check which bins contain the occurrence points
    points_joined <- sf::st_join(bins_sf[, bin_id_column], sp_occurrences, join = sf::st_intersects)
    bins_with_points <- unique(points_joined[[bin_id_column]][!is.na(points_joined$spp)])
    bins_with_points <- bins_with_points[!is.na(bins_with_points)]
    
    # Combine both
    bins_with_presence <- sort(unique(c(bins_with_mst, bins_with_points)))
    
    # Fill matrix
    for (bin_name in bins_with_presence) {
      if (bin_name %in% rownames(pres_abs_matrix)) {
        pres_abs_matrix[bin_name, sp] <- 1
      }
    }
  }
  
  pres_abs_clean <- pres_abs_matrix[!is.na(rownames(pres_abs_matrix)), ]
  pres_abs_clean <- pres_abs_clean[rowSums(pres_abs_clean) > 0, ]
  
  # Store original names
  original_names <- rownames(pres_abs_clean)
  
  # Use sequential codes (A, B, C...) instead of 2-letter codes
  n_ranges <- nrow(pres_abs_clean)
  sequential_codes <- LETTERS[1:n_ranges]
  if (n_ranges > 26) {
    sequential_codes <- c(LETTERS, paste0(LETTERS, LETTERS[1:(n_ranges - 26)]))
  }
  rownames(pres_abs_clean) <- sequential_codes
  
  # Pass original names to save_biogeobears_format
  result <- save_biogeobears_format(pres_abs_clean, output_dir, "MST", original_names)
  return(result)
}

# ============================================================================
# 2. CONVEX HULL METHOD
# ============================================================================

calcRange_convexhull_biogeobears <- function(occurrence_data, irregular_polygons, bin_id_column, crs_input = 4326, output_dir = "out_MPC_biogeobears/") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  names(occurrence_data) <- tolower(names(occurrence_data))
  dat <- occurrence_data[, c("spp", "long", "lat")]
  dat <- dat[complete.cases(dat), ]
  
  occurrences_sf <- sf::st_as_sf(dat, coords = c("long", "lat"), crs = crs_input)
  
  # Shapefile should already be sf (converted in Step 1 for BioGeoBEARS)
  bins_sf <- irregular_polygons
  
  # Check and fix invalid geometries
  invalid_geoms <- !sf::st_is_valid(bins_sf)
  if (any(invalid_geoms)) {
    bins_sf <- sf::st_make_valid(bins_sf)
    
    # If still invalid after st_make_valid, filter them out
    invalid_after <- !sf::st_is_valid(bins_sf)
    if (any(invalid_after)) {
      bins_sf <- bins_sf[!invalid_after, ]
    } else {
    }
  }
  
  occurrences_transformed <- sf::st_transform(occurrences_sf, sf::st_crs(bins_sf))
  
  all_species <- sort(unique(dat$spp))
  all_bins <- sort(unique(bins_sf[[bin_id_column]]))
  
  pres_abs_matrix <- matrix(0, nrow = length(all_bins), ncol = length(all_species), dimnames = list(all_bins, all_species))
  
  for (sp in all_species) {
    sp_occurrences <- occurrences_transformed[occurrences_transformed$spp == sp, ]
    n_points <- nrow(sp_occurrences)
    
    if (n_points < 3) {
      next
    }
    
    # Calculate convex hull
    hull_geometry <- sf::st_convex_hull(sf::st_union(sp_occurrences))
    hull_sf <- sf::st_sf(species = sp, geometry = sf::st_sfc(hull_geometry, crs = sf::st_crs(bins_sf)))
    
    # Use st_join to find which bins the hull intersects
    hull_joined <- sf::st_join(bins_sf[, bin_id_column], hull_sf, join = sf::st_intersects)
    bins_with_hull <- unique(hull_joined[[bin_id_column]][!is.na(hull_joined$species)])
    bins_with_hull <- bins_with_hull[!is.na(bins_with_hull)]
    
    # Also check which bins contain the occurrence points
    points_joined <- sf::st_join(bins_sf[, bin_id_column], sp_occurrences, join = sf::st_intersects)
    bins_with_points <- unique(points_joined[[bin_id_column]][!is.na(points_joined$spp)])
    bins_with_points <- bins_with_points[!is.na(bins_with_points)]
    
    # Combine both
    bins_with_presence <- sort(unique(c(bins_with_hull, bins_with_points)))
    
    # Fill matrix
    for (bin_name in bins_with_presence) {
      if (bin_name %in% rownames(pres_abs_matrix)) {
        pres_abs_matrix[bin_name, sp] <- 1
      }
    }
  }
  
  pres_abs_clean <- pres_abs_matrix[!is.na(rownames(pres_abs_matrix)), ]
  pres_abs_clean <- pres_abs_clean[rowSums(pres_abs_clean) > 0, ]
  
  # Store original names before replacing with codes
  original_names <- rownames(pres_abs_clean)
  
  # Use sequential codes (A, B, C...) instead of 2-letter codes
  n_ranges <- nrow(pres_abs_clean)
  sequential_codes <- LETTERS[1:n_ranges]
  if (n_ranges > 26) {
    sequential_codes <- c(LETTERS, paste0(LETTERS, LETTERS[1:(n_ranges - 26)]))
  }
  rownames(pres_abs_clean) <- sequential_codes
  
  result <- save_biogeobears_format(pres_abs_clean, output_dir, "MPC", original_names)
  return(result)
}

# ============================================================================
# 3. BUFFER METHOD
# ============================================================================

calcRange_buffer_biogeobears <- function(occurrence_data, irregular_polygons, bin_id_column, buffer_width_km = 100, crs_input = 4326, output_dir = "out_BUFF_biogeobears/") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  names(occurrence_data) <- tolower(names(occurrence_data))
  dat <- occurrence_data[, c("spp", "long", "lat")]
  dat <- dat[complete.cases(dat), ]
  
  occurrences_sf <- sf::st_as_sf(dat, coords = c("long", "lat"), crs = crs_input)
  
  # Shapefile should already be sf (converted in Step 1 for BioGeoBEARS)
  bins_sf <- irregular_polygons
  
  # Check and fix invalid geometries
  invalid_geoms <- !sf::st_is_valid(bins_sf)
  if (any(invalid_geoms)) {
    bins_sf <- sf::st_make_valid(bins_sf)
    
    # If still invalid after st_make_valid, filter them out
    invalid_after <- !sf::st_is_valid(bins_sf)
    if (any(invalid_after)) {
      bins_sf <- bins_sf[!invalid_after, ]
    } else {
    }
  }
  
  occurrences_transformed <- sf::st_transform(occurrences_sf, sf::st_crs(bins_sf))
  
  all_species <- sort(unique(dat$spp))
  all_bins <- sort(unique(bins_sf[[bin_id_column]]))
  
  pres_abs_matrix <- matrix(0, nrow = length(all_bins), ncol = length(all_species), dimnames = list(all_bins, all_species))
  
  # Detect CRS and convert buffer width to appropriate units
  crs_obj <- sf::st_crs(bins_sf)
  if (is.na(crs_obj$proj4string) || grepl("longlat", crs_obj$proj4string)) {
    # Geographic CRS - project to Web Mercator for distance calculations
    bins_sf_proj <- sf::st_transform(bins_sf, 3857)
    occurrences_proj <- sf::st_transform(occurrences_transformed, 3857)
    buffer_dist_m <- buffer_width_km * 1000
  } else {
    # Already projected
    bins_sf_proj <- bins_sf
    occurrences_proj <- occurrences_transformed
    buffer_dist_m <- buffer_width_km * 1000
  }
  
  for (sp in all_species) {
    sp_occurrences <- occurrences_proj[occurrences_proj$spp == sp, ]
    n_points <- nrow(sp_occurrences)
    
    if (n_points < 1) {
      next
    }
    
    # Calculate buffer
    buffer_geometry <- sf::st_buffer(sf::st_union(sp_occurrences), dist = buffer_dist_m)
    buffer_sf <- sf::st_sf(species = sp, geometry = sf::st_sfc(buffer_geometry, crs = sf::st_crs(bins_sf_proj)))
    
    # Use st_join to find which bins the buffer intersects
    buffer_joined <- sf::st_join(bins_sf_proj[, bin_id_column], buffer_sf, join = sf::st_intersects)
    bins_with_buffer <- unique(buffer_joined[[bin_id_column]][!is.na(buffer_joined$species)])
    bins_with_buffer <- bins_with_buffer[!is.na(bins_with_buffer)]
    
    # Also check which bins contain the occurrence points
    points_joined <- sf::st_join(bins_sf_proj[, bin_id_column], sp_occurrences, join = sf::st_intersects)
    bins_with_points <- unique(points_joined[[bin_id_column]][!is.na(points_joined$spp)])
    bins_with_points <- bins_with_points[!is.na(bins_with_points)]
    
    # Combine both
    bins_with_presence <- sort(unique(c(bins_with_buffer, bins_with_points)))
    
    # Fill matrix
    for (bin_name in bins_with_presence) {
      if (bin_name %in% rownames(pres_abs_matrix)) {
        pres_abs_matrix[bin_name, sp] <- 1
      }
    }
  }
  
  pres_abs_clean <- pres_abs_matrix[!is.na(rownames(pres_abs_matrix)), ]
  pres_abs_clean <- pres_abs_clean[rowSums(pres_abs_clean) > 0, ]
  
  # Store original names
  original_names <- rownames(pres_abs_clean)
  
  # Use sequential codes (A, B, C...) instead of 2-letter codes
  n_ranges <- nrow(pres_abs_clean)
  sequential_codes <- LETTERS[1:n_ranges]
  if (n_ranges > 26) {
    sequential_codes <- c(LETTERS, paste0(LETTERS, LETTERS[1:(n_ranges - 26)]))
  }
  rownames(pres_abs_clean) <- sequential_codes
  
  # Pass original names to save_biogeobears_format
  result <- save_biogeobears_format(pres_abs_clean, output_dir, "BUFF", original_names)
  return(result)
}
