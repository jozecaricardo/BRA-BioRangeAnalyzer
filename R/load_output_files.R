#' Load Shapefiles from Output Directories
#'
#' Searches for and loads all shapefiles from output directories created by
#' range extrapolation functions (out_buffers, out_MCP, out_MST, etc.)
#'
#' @return A list containing:
#'   - shapefiles: named list of loaded shapefiles (terra objects)
#'   - rasters: named list of loaded rasters
#'   - directories: list of directories searched
#'
#' @export
load_output_shapefiles <- function() {
  
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out_mst", "out", "out_domains", "out_irregular_bins")
  
  shapefiles <- list()
  rasters <- list()
  found_dirs <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      next
    }
    
    found_dirs <- c(found_dirs, dir)
    
    # Load shapefiles
    shp_files <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)
    
    for (shp_file in shp_files) {
      tryCatch({
        shp_name <- tools::file_path_sans_ext(basename(shp_file))
        shapefiles[[paste0(dir, "/", shp_name)]] <- terra::vect(shp_file)
      }, error = function(e) {
        warning("Could not load shapefile: ", shp_file, " - ", e$message)
      })
    }
    
    # Load rasters
    tif_files <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
    
    for (tif_file in tif_files) {
      tryCatch({
        tif_name <- tools::file_path_sans_ext(basename(tif_file))
        rasters[[paste0(dir, "/", tif_name)]] <- terra::rast(tif_file)
      }, error = function(e) {
        warning("Could not load raster: ", tif_file, " - ", e$message)
      })
    }
  }
  
  return(list(
    shapefiles = shapefiles,
    rasters = rasters,
    directories = found_dirs
  ))
}

#' Convert Shapefile to Leaflet-compatible format
#'
#' Converts terra SpatVector to sf for use with leaflet
#'
#' @param shapefile A terra SpatVector object
#'
#' @return An sf object
#'
#' @export
shapefile_to_sf <- function(shapefile) {
  tryCatch({
    sf::st_as_sf(shapefile)
  }, error = function(e) {
    NULL
  })
}

#' Convert Raster to Leaflet-compatible format
#'
#' Converts terra raster to RasterLayer for use with leaflet
#'
#' @param raster A terra raster object
#'
#' @return A RasterLayer object
#'
#' @export
raster_to_raster_layer <- function(raster) {
  tryCatch({
    raster::raster(raster)
  }, error = function(e) {
    NULL
  })
}

#' Get list of available output directories
#'
#' @return Character vector of directories that exist
#'
#' @export
get_available_output_dirs <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out_mst", "out", "out_domains", "out_irregular_bins")
  available <- output_dirs[sapply(output_dirs, dir.exists)]
  return(available)
}

#' Get list of shapefile names from output directories
#'
#' @return Character vector of shapefile names
#'
#' @export
get_shapefile_names <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out_mst", "out", "out_domains", "out_irregular_bins")
  
  shp_names <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) next
    
    shp_files <- list.files(dir, pattern = "\\.shp$")
    shp_names <- c(shp_names, tools::file_path_sans_ext(shp_files))
  }
  
  return(shp_names)
}

#' Get list of raster names from output directories
#'
#' @return Character vector of raster names
#'
#' @export
get_raster_names <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out_mst", "out", "out_domains", "out_irregular_bins")
  
  tif_names <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) next
    
    tif_files <- list.files(dir, pattern = "\\.tif$")
    tif_names <- c(tif_names, tools::file_path_sans_ext(tif_files))
  }
  
  return(tif_names)
}
