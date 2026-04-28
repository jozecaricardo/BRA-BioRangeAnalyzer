#' Wrapper for Range Extrapolation Functions
#'
#' Calls the appropriate range extrapolation function based on method selection
#' with proper error handling and fallback to mock data
#'
#' @param occurrence_data Data frame with columns: spp, long, lat
#' @param method Character string: "BUFF", "MPC", or "MST"
#' @param grid_resolution Numeric, grid resolution in degrees
#' @param buffer_width Numeric, buffer width in meters (for BUFF method)
#'
#' @return List with $pres_abs (matrix) and $geometry (sf object or NULL)
#'
#' @export
run_extrapolation <- function(occurrence_data, method = "BUFF", 
                              grid_resolution = 1, buffer_width = 100000) {
  
  tryCatch({
    # Ensure we have the required columns
    if (!all(c("spp", "long", "lat") %in% names(occurrence_data))) {
      stop("Data must have columns: spp, long, lat")
    }
    
    # Call the appropriate function based on method
    if (method == "BUFF") {
      result <- calcRange_buffers(
        xy = occurrence_data,
        buffer.width = buffer_width,
        shape_file = NULL,
        resol = grid_resolution,
        mean_dist = FALSE
      )
    } else if (method == "MPC") {
      result <- calcRange_convexHull(
        xy = occurrence_data,
        shape_file = NULL,
        resol = grid_resolution
      )
    } else if (method == "MST") {
      result <- calcRange_irregular_bins(
        xy = occurrence_data,
        bins_shapefile = NULL,
        resol = c(grid_resolution, grid_resolution)
      )
    } else {
      stop("Unknown method: ", method)
    }
    
    # Ensure result has the expected structure
    if (!is.list(result)) {
      stop("Function did not return a list")
    }
    
    # Extract pres_abs
    pres_abs <- if ("pres_abs" %in% names(result)) {
      result$pres_abs
    } else if ("pres.abs" %in% names(result)) {
      result$pres.abs
    } else if (is.matrix(result)) {
      result
    } else {
      stop("Could not extract presence-absence matrix from result")
    }
    
    # Extract geometry if available
    geometry <- if ("geometry" %in% names(result)) {
      result$geometry
    } else {
      NULL
    }
    
    # Try to convert geometry to sf if it's not NULL
    if (!is.null(geometry)) {
      tryCatch({
        geometry <- convert_geometry_to_sf(geometry)
      }, error = function(e) {
        warning("Could not convert geometry to sf: ", e$message)
        geometry <<- NULL
      })
    }
    
    # Return standardized result
    list(
      pres_abs = pres_abs,
      geometry = geometry,
      method = method,
      success = TRUE
    )
    
  }, error = function(e) {
    # Return error result
    warning("Extrapolation failed: ", e$message)
    list(
      pres_abs = NULL,
      geometry = NULL,
      method = method,
      success = FALSE,
      error = e$message
    )
  })
}

#' Create Mock Extrapolation Result
#'
#' Creates a mock result for testing when real functions fail
#'
#' @param occurrence_data Data frame with species occurrences
#' @param method Character string: "BUFF", "MPC", or "MST"
#'
#' @return List with $pres_abs (matrix) and $geometry (sf object with mock polygons)
#'
#' @export
create_mock_extrapolation <- function(occurrence_data, method = "BUFF") {
  
  tryCatch({
    # Create a simple presence-absence matrix
    species <- unique(occurrence_data$spp)
    n_species <- length(species)
    
    # Create a grid based on the extent of the data
    lon_range <- range(occurrence_data$long, na.rm = TRUE)
    lat_range <- range(occurrence_data$lat, na.rm = TRUE)
    
    # Create a simple grid
    lon_seq <- seq(lon_range[1] - 1, lon_range[2] + 1, by = 1)
    lat_seq <- seq(lat_range[1] - 1, lat_range[2] + 1, by = 1)
    
    n_cells <- length(lon_seq) * length(lat_seq)
    
    # Create random presence-absence matrix
    pres_abs <- matrix(
      sample(0:1, n_cells * n_species, replace = TRUE, prob = c(0.7, 0.3)),
      nrow = n_cells,
      ncol = n_species,
      dimnames = list(
        paste0("Cell_", 1:n_cells),
        species
      )
    )
    
    # Create mock geometry as sf polygons
    geometry <- NULL
    tryCatch({
      # Try to create sf polygons
      if (requireNamespace("sf", quietly = TRUE)) {
        polygons <- list()
        
        for (i in 1:(length(lon_seq) - 1)) {
          for (j in 1:(length(lat_seq) - 1)) {
            # Create a simple square polygon
            coords <- matrix(c(
              lon_seq[i], lat_seq[j],
              lon_seq[i+1], lat_seq[j],
              lon_seq[i+1], lat_seq[j+1],
              lon_seq[i], lat_seq[j+1],
              lon_seq[i], lat_seq[j]
            ), ncol = 2, byrow = TRUE)
            
            polygons[[length(polygons) + 1]] <- sf::st_polygon(list(coords))
          }
        }
        
        # Create sf object
        geometry <- sf::st_sfc(polygons, crs = 4326)
        geometry <- sf::st_sf(geometry = geometry)
      }
    }, error = function(e) {
      warning("Could not create sf geometry: ", e$message)
    })
    
    list(
      pres_abs = pres_abs,
      geometry = geometry,
      method = method,
      is_mock = TRUE
    )
    
  }, error = function(e) {
    warning("Could not create mock extrapolation: ", e$message)
    list(
      pres_abs = NULL,
      geometry = NULL,
      method = method,
      is_mock = TRUE,
      error = e$message
    )
  })
}
