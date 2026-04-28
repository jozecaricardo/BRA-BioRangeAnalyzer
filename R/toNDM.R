#' toNDM
#'
#' This function converts species data into the XYD format required by NDM (VNDM)
#' software for identifying areas of endemism. The function accepts three types of input:
#' (1) occurrence data from CSV/data frames/matrices, (2) shapefiles from extrapolated
#' distributions (buffers, convex hulls, MST), or (3) a directory containing multiple
#' species shapefiles. For shapefiles, the function creates a grid with specified
#' resolution and uses grid cell centroids as occurrence points.
#'
#' @param input_data Either:
#'   \itemize{
#'     \item A character string with path to CSV file (occurrence data)
#'     \item A data frame with 3 columns: species, longitude, latitude
#'     \item A matrix/table with 3 columns: species, longitude, latitude
#'     \item A character string with path to shapefile directory (e.g., "out_buffers/")
#'     \item A character string with "SHAPEFILE" to trigger shapefile mode
#'   }
#' @param shapefile_dir Character string specifying directory containing species
#'   shapefiles (e.g., "out_buffers/", "out_MPC/", "out/"). Required when using
#'   shapefile mode. Each shapefile should be named with species identifier.
#' @param resolution Numeric vector of length 2 specifying grid resolution in degrees
#'   (e.g., c(5, 5) for 5x5 degree cells). Required for shapefile mode. This should
#'   match the resolution used in the extrapolation analysis.
#' @param shape_file Optional spatial object (sf, SpatVector, or Spatial*) representing
#'   the study region boundary. Used to crop the grid in shapefile mode.
#' @param output_file Character string specifying output XYD file name
#'   (default: "output_NDM.xyd").
#' @param separator Character string specifying field separator for CSV files (default: ",").
#' @param header Logical, whether CSV has column names (default: TRUE).
#' @param na_strings Character vector of strings to interpret as NA (default: c("?", "-", "NA", "")).
#' @param remove_na Logical, if TRUE removes rows with NA values (default: TRUE).
#' @param output_dir Character string specifying output directory (default: "out_NDM/").
#'
#' @return The function does not return a value. It creates an XYD-formatted file
#'   ready for NDM/VNDM analysis. Invisibly returns a list with conversion statistics.
#'
#' @details
#' The function operates in two modes:
#'
#' **Mode 1: Occurrence Data**
#' \itemize{
#'   \item Reads occurrence coordinates from CSV, data frame, or matrix
#'   \item Removes NA values and empty species
#'   \item Formats for NDM
#' }
#'
#' **Mode 2: Shapefile (Extrapolated Distributions)**
#' \itemize{
#'   \item Reads species shapefiles from specified directory
#'   \item Creates grid with specified resolution
#'   \item Identifies grid cells intersecting each species shapefile
#'   \item Calculates centroid coordinates for each grid cell
#'   \item Uses centroids as "occurrence" points for NDM
#'   \item Exports to XYD format
#' }
#'
#' @section Shapefile Mode Workflow:
#' When using extrapolated distributions (buffers, convex hulls, MST):
#' \enumerate{
#'   \item Run extrapolation analysis (e.g., \code{calcRange_buffers()})
#'   \item Shapefiles are saved in output directory (e.g., "out_buffers/")
#'   \item Use \code{toNDM()} with shapefile_dir and resolution parameters
#'   \item Function reads all shapefiles, creates grid, extracts centroids
#'   \item Generates XYD file with centroid coordinates
#'   \item Analyze in NDM using grid-based "occurrences"
#' }
#'
#' @section Input Formats:
#'
#' **Format 1: Occurrence Data (CSV/Data Frame)**
#' \preformatted{
#' species,longitude,latitude
#' Belostoma_amazonum,-60.5,-3.2
#' Belostoma_angustum,-58.7,-5.3
#' }
#'
#' **Format 2: Shapefiles in Directory**
#' \preformatted{
#' out_buffers/
#'    Belostoma_amazonum.shp
#'    Belostoma_angustum.shp
#'    Belostoma_anurum.shp
#' }
#'
#' @note
#' \itemize{
#'   \item For occurrence data: input must have 3 columns (species, lon, lat)
#'   \item For shapefiles: resolution parameter is required
#'   \item Shapefile names should contain species identifiers
#'   \item Grid resolution should match the one used in extrapolation analysis
#'   \item Centroids are calculated for grid cells intersecting species ranges
#'   \item Output directory is created automatically if needed
#' }
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # MODE 1: From occurrence data (original functionality)
#' # ============================================================================
#'
#' # Example 1: From CSV
#' toNDM(input_data = "occurrences.csv")
#'
#' # Example 2: From data frame
#' occurrences <- data.frame(
#'   species = c("Species_A", "Species_B"),
#'   longitude = c(-60.5, -58.7),
#'   latitude = c(-3.2, -5.3)
#' )
#' toNDM(input_data = occurrences)
#'
#' # ============================================================================
#' # MODE 2: From extrapolated distributions (NEW!)
#' # ============================================================================
#'
#' # Example 3: From buffer shapefiles
#' # Step 1: Run buffer analysis
#' calcRange_buffers(
#'   xy = lycipta_final,
#'   shape_file = neo,
#'   resol = 10,
#'   buffer.width = 500000
#' )
#' # This creates shapefiles in out_buffers/
#'
#' # Step 2: Convert shapefiles to NDM format
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(10, 10),
#'   shape_file = neo,
#'   output_file = "buffer_extrapolated_NDM.xyd"
#' )
#'
#' # Example 4: From convex hull shapefiles
#' # Step 1: Run convex hull analysis
#' calcRange_convexHull(
#'   xy = species_coords,
#'   shape_file = south_america,
#'   resol = 5
#' )
#' # This creates shapefiles in out_MPC/
#'
#' # Step 2: Convert to NDM
#' toNDM(
#'   shapefile_dir = "out_MPC/",
#'   resolution = c(5, 5),
#'   shape_file = south_america,
#'   output_file = "convexhull_extrapolated_NDM.xyd"
#' )
#'
#' # Example 5: From MST shapefiles
#' # Step 1: Run MST analysis
#' calcRange_mst(
#'   xy = belostomatidae_coords,
#'   shape_file = neotropics,
#'   resol = 10
#' )
#' # This creates shapefiles in out/
#'
#' # Step 2: Convert to NDM
#' toNDM(
#'   shapefile_dir = "out/",
#'   resolution = c(10, 10),
#'   shape_file = neotropics,
#'   output_file = "mst_extrapolated_NDM.xyd"
#' )
#'
#' # ============================================================================
#' # Example 6: Complete workflow - Buffer to NDM
#' # ============================================================================
#'
#' library(sf)
#'
#' # Step 1: Load occurrence data
#' occurrences <- read.csv("species_occurrences.csv")
#'
#' # Step 2: Load study area shapefile
#' study_area <- st_read("study_area.shp")
#'
#' # Step 3: Calculate buffer ranges
#' calcRange_buffers(
#'   xy = occurrences,
#'   shape_file = study_area,
#'   resol = 5,
#'   buffer.width = 300000  # 300 km
#' )
#'
#' # Step 4: Convert extrapolated ranges to NDM format
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = study_area,
#'   output_file = "buffer_300km_NDM.xyd"
#' )
#'
#' # Step 5: Analyze in NDM/VNDM
#' # - Open VNDM
#' # - Load out_NDM/buffer_300km_NDM.xyd
#' # - Set parameters
#' # - Run analysis
#'
#' # ============================================================================
#' # Example 7: Comparing raw occurrences vs extrapolated ranges
#' # ============================================================================
#'
#' # Analysis 1: Raw occurrences
#' toNDM(
#'   input_data = "occurrences.csv",
#'   output_file = "raw_occurrences_NDM.xyd",
#'   output_dir = "NDM_comparison/raw/"
#' )
#'
#' # Analysis 2: Buffer extrapolation (100 km)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 5, buffer.width = 100000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "buffer_100km_NDM.xyd",
#'   output_dir = "NDM_comparison/buffer_100/"
#' )
#'
#' # Analysis 3: Buffer extrapolation (500 km)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 5, buffer.width = 500000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "buffer_500km_NDM.xyd",
#'   output_dir = "NDM_comparison/buffer_500/"
#' )
#'
#' # Analysis 4: Convex hull extrapolation
#' calcRange_convexHull(xy = coords, shape_file = shape, resol = 5)
#' toNDM(
#'   shapefile_dir = "out_MPC/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "convexhull_NDM.xyd",
#'   output_dir = "NDM_comparison/hull/"
#' )
#'
#' # Compare all 4 analyses in NDM to assess method sensitivity
#'
#' # ============================================================================
#' # Example 8: Different grid resolutions
#' # ============================================================================
#'
#' # Coarse resolution (10 degrees)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 10, buffer.width = 300000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(10, 10),
#'   output_file = "buffer_10deg_NDM.xyd"
#' )
#'
#' # Fine resolution (2 degrees)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 2, buffer.width = 300000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(2, 2),
#'   output_file = "buffer_2deg_NDM.xyd"
#' )
#' }
#'
#' @seealso
#' \code{\link{calcRange_buffers}} for buffer-based range extrapolation
#' \code{\link{calcRange_convexHull}} for convex hull range extrapolation
#' \code{\link{calcRange_mst}} for MST-based range extrapolation
#'
#' @references
#' Szumik, C. A., & Goloboff, P. A. (2004). Areas of endemism: an improved
#' optimality criterion. Systematic Biology, 53(6), 968-977.
#'
#' Goloboff, P. A. (2004). NDM/VNDM, programs for identification of areas of
#' endemism. Program and documentation available at: www.lillo.org.ar/phylogeny/endemism
#'
#' @export
toNDM <- function(input_data = NULL,
                  shapefile_dir = NULL,
                  resolution = NULL,
                  shape_file = NULL,
                  output_file = "output_NDM.xyd",
                  separator = ",",
                  header = TRUE,
                  na_strings = c("?", "-", "NA", ""),
                  remove_na = TRUE,
                  output_dir = "out_NDM/") {

  ######################
  ## Input validation ##
  ######################

  # Check if at least one input is provided
  if (is.null(input_data) && is.null(shapefile_dir)) {
    stop("Error: Either 'input_data' or 'shapefile_dir' must be provided.")
  }

  ######################
  ## Detect mode #######
  ######################

  mode <- NULL

  if (!is.null(shapefile_dir)) {
    mode <- "shapefile"
    message("Mode: Shapefile (extrapolated distributions)")

    # Validate shapefile mode requirements
    if (is.null(resolution)) {
      stop("Error: 'resolution' parameter is required for shapefile mode.
           Example: resolution = c(5, 5)")
    }

    if (!dir.exists(shapefile_dir)) {
      stop(paste0("Error: Shapefile directory not found: ", shapefile_dir))
    }

  } else {
    mode <- "occurrences"
    message("Mode: Occurrence data")
  }

  ######################
  ## Process data ######
  ######################

  if (mode == "occurrences") {
    # MODE 1: Process occurrence data
    dados <- process_occurrences(input_data, separator, header, na_strings, remove_na)

  } else if (mode == "shapefile") {
    # MODE 2: Process shapefiles
    dados <- process_shapefiles(shapefile_dir, resolution, shape_file)
  }

  ######################
  ## Format for NDM ####
  ######################

  output <- format_for_ndm(dados)

  ######################
  ## Write output ######
  ######################

  write_ndm_file(output, output_file, output_dir)

  ######################
  ## Report results ####
  ######################

  report_conversion(dados, output_file, output_dir, mode)

  # Return summary invisibly
  invisible(list(
    output_file = file.path(output_dir, output_file),
    n_occurrences = nrow(dados),
    n_species = length(unique(dados[[1]])),
    mode = mode
  ))
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Process occurrence data
#' @keywords internal
process_occurrences <- function(input_data, separator, header, na_strings, remove_na) {

  if (is.null(input_data)) {
    stop("Error: 'input_data' is NULL.")
  }

  # Read data based on input type
  if (is.character(input_data)) {
    if (!file.exists(input_data)) {
      stop(paste0("Error: File not found: ", input_data))
    }
    message(paste0("Reading data from file: ", input_data))
    dados <- read.csv(
      file = input_data,
      sep = separator,
      header = header,
      na.strings = na_strings,
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(input_data)) {
    message("Using provided data frame")
    dados <- input_data
  } else if (is.matrix(input_data) || is.table(input_data)) {
    message("Converting matrix/table to data frame")
    dados <- as.data.frame(input_data, stringsAsFactors = FALSE)
  } else if(is.list(input_data)){ # Extract matrix from list if needed
    if("data_df" %in% names(input_data)){
      dados <- input_data$data_df
    } else {
      stop("Error: List does not contain 'data_df' component.")
    }
  } else {
    stop("Error: 'input_data' must be a file path, data frame, matrix, singleton_to_data_frame object, or table.")
  }

  # Validate 3 columns
  if (ncol(dados) != 3) {
    stop(paste0("Error: Data must have 3 columns (species, lon, lat). Found ", ncol(dados)))
  }

  # Convert to appropriate types
  dados[[1]] <- as.character(dados[[1]])
  dados[[2]] <- as.numeric(as.character(dados[[2]]))
  dados[[3]] <- as.numeric(as.character(dados[[3]]))

  # Remove NA
  if (remove_na) {
    dados <- na.omit(dados)
    rownames(dados) <- NULL
  }

  if (nrow(dados) == 0) {
    stop("Error: No valid data after removing NA values.")
  }

  # Sort by species
  dados <- dados[order(dados[[1]]), ]
  rownames(dados) <- NULL

  return(dados)
}


#' Process shapefiles from extrapolated distributions
#' @keywords internal
process_shapefiles <- function(shapefile_dir, resolution, shape_file) {

  # Load required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for shapefile mode. Install with: install.packages('sf')")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for shapefile mode. Install with: install.packages('terra')")
  }

  library(sf)
  library(terra)

  message(paste0("Reading shapefiles from: ", shapefile_dir))

  # Find all shapefiles in directory
  all_shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)

  # Filter to keep only range shapefiles (BUFF_, MPC_, MST_)
  # Exclude point shapefiles
  shapefiles <- all_shapefiles[grepl("(BUFF_|MPC_|MST_)", basename(all_shapefiles))]

  if (length(shapefiles) == 0) {
    stop(paste0("Error: No range shapefiles (BUFF_*, MPC_*, MST_*) found in directory: ", shapefile_dir,
                "\nFound ", length(all_shapefiles), " total shapefiles, but none with BUFF_, MPC_, or MST_ prefix."))
  }

  message(paste0("Found ", length(shapefiles), " shapefiles"))

  # Load study area if provided
  if (!is.null(shape_file)) {
    if (inherits(shape_file, "sf")) {
      study_area <- shape_file
    } else if (inherits(shape_file, "SpatVector")) {
      study_area <- st_as_sf(shape_file)
    } else {
      study_area <- st_read(shape_file, quiet = TRUE)
    }
    bbox <- st_bbox(study_area)
  } else {
    # Use bounding box from all shapefiles
    all_bbox <- lapply(shapefiles, function(f) st_bbox(st_read(f, quiet = TRUE)))
    bbox <- c(
      xmin = min(sapply(all_bbox, function(b) b["xmin"])),
      ymin = min(sapply(all_bbox, function(b) b["ymin"])),
      xmax = max(sapply(all_bbox, function(b) b["xmax"])),
      ymax = max(sapply(all_bbox, function(b) b["ymax"]))
    )
  }

  # Create grid
  message(paste0("Creating grid with resolution: ", resolution[1], " x ", resolution[2], " degrees"))

  grid_polygon <- st_make_grid(
    st_as_sfc(bbox),
    cellsize = resolution,
    what = "polygons"
  )
  grid_sf <- st_sf(geometry = grid_polygon)
  grid_sf$grid_id <- 1:nrow(grid_sf)

  # Calculate centroids
  centroids <- st_centroid(grid_sf)
  centroid_coords <- st_coordinates(centroids)
  grid_sf$centroid_lon <- centroid_coords[, 1]
  grid_sf$centroid_lat <- centroid_coords[, 2]

  # Process each shapefile
  occurrences_list <- list()

  for (shp_file in shapefiles) {
    # Extract species name from filename (remove BUFF_, MPC_, or MST_ prefix)
    species_name <- tools::file_path_sans_ext(basename(shp_file))
    species_name <- gsub("^(BUFF_|MPC_|MST_)", "", species_name)  # Remove prefix
    species_name <- gsub("_", " ", species_name)  # Replace underscores with spaces

    message(paste0("Processing: ", species_name))

    # Read shapefile
    species_range <- st_read(shp_file, quiet = TRUE)

    # Find intersecting grid cells
    intersects <- st_intersects(grid_sf, species_range, sparse = FALSE)
    intersecting_cells <- grid_sf[apply(intersects, 1, any), ]

    if (nrow(intersecting_cells) > 0) {
      # Create occurrence records from centroids
      for (i in 1:nrow(intersecting_cells)) {
        occurrences_list[[length(occurrences_list) + 1]] <- data.frame(
          species = species_name,
          longitude = intersecting_cells$centroid_lon[i],
          latitude = intersecting_cells$centroid_lat[i],
          stringsAsFactors = FALSE
        )
      }
    } else {
      warning(paste0("No grid cells found for species: ", species_name))
    }
  }

  # Combine all occurrences
  if (length(occurrences_list) == 0) {
    stop("Error: No valid occurrences generated from shapefiles.")
  }

  dados <- do.call(rbind, occurrences_list)
  rownames(dados) <- NULL

  # Sort by species
  dados <- dados[order(dados$species), ]

  message(paste0("Generated ", nrow(dados), " grid cell centroids from ",
                 length(unique(dados$species)), " species"))

  return(dados)
}


#' Format data for NDM
#' @keywords internal
format_for_ndm <- function(dados) {

  output <- list()
  n <- 0
  i <- 1
  t <- 1
  name <- ""
  size <- nrow(dados)
  final_species <- length(unique(dados[[1]]))

  # Header
  output$c1[t] <- 'longlat'
  output$c2[t] <- ''
  t <- t + 1

  output$c1[t] <- paste('spp', final_species)
  output$c2[t] <- ''
  t <- t + 1

  output$c1[t] <- 'xydata'
  output$c2[t] <- ''
  t <- t + 1

  # Data
  while(i <= size) {
    if (as.character(dados[[1]][i]) != name) {
      output$c1[t] <- paste0("sp ", n)
      output$c2[t] <- paste0("[", as.character(dados[[1]][i]), "]")
      name <- as.character(dados[[1]][i])
      n <- n + 1
      t <- t + 1
    }

    output$c1[t] <- paste(as.character(dados[[2]][i]), ',', sep = '')
    output$c2[t] <- as.character(dados[[3]][i])
    t <- t + 1
    i <- i + 1
  }

  output <- na.omit(output)
  rownames(output) <- NULL

  return(output)
}


#' Write NDM file
#' @keywords internal
write_ndm_file <- function(output, output_file, output_dir) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste0("Created output directory: ", output_dir))
  }

  output_path <- file.path(output_dir, output_file)

  write.table(
    output,
    file = output_path,
    sep = "",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


#' Report conversion results
#' @keywords internal
report_conversion <- function(dados, output_file, output_dir, mode) {

  output_path <- file.path(output_dir, output_file)

  message("\n[OK] NDM/VNDM XYD file created successfully")
  message(paste0("  Output file: ", output_path))
  message(paste0("  Mode: ", mode))

  if (mode == "shapefile") {
    message(paste0("  Grid cell centroids: ", nrow(dados)))
  } else {
    message(paste0("  Occurrences: ", nrow(dados)))
  }

  message(paste0("  Species: ", length(unique(dados[[1]]))))
  message("\nNext steps:")
  message("  1. Open NDM/VNDM software")
  message(paste0("  2. Load file: ", output_path))
  message("  3. Set parameters and run analysis")
}
