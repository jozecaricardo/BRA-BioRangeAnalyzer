#' extrapolation_plot_from_shapeFiles
#'
#' Searches for, loads, and plots species distribution shapefiles from specified directories.
#' Can work with either a list of SpatVector objects or directly from directories containing
#' shapefiles (e.g., minimum convex polygons or buffers). Automatically assigns colors and
#' generates legends for multiple species.
#'
#' @param species_names Character vector with one or more species names to plot. If NULL,
#'   lists all available species in the directory without plotting
#' @param spatvector_list Optional list of SpatVector objects to search through. If provided,
#'   the function searches this list instead of loading from directories
#' @param directory Character string or vector specifying directory(ies) to search (default: "out_MCP").
#'   Can be "out_MCP", "out_buffers", or any custom directory path. If a vector is provided,
#'   searches all specified directories
#' @param file_type Character string specifying the file type pattern to match:
#'   "MCP" for minimum convex polygons, "BUFF" for buffer files, 'mst' for MSTs, or a custom pattern
#'   (default: "MCP")
#' @param add Logical, if TRUE adds to existing plot, if FALSE creates new plot (default: TRUE)
#' @param colors Named vector of colors for each species, or NULL for automatic colors using viridis palette
#' @param alpha Numeric value between 0 and 1 for transparency (default: 0.5)
#' @param legend Logical, if TRUE shows legend with species names (default: TRUE)
#' @param legend_position Character string specifying legend position: "topright", "topleft",
#'   "bottomright", "bottomleft", "top", "bottom", "left", "right", or "center" (default: "topright")
#' @param legend_inset Numeric value indicating the percentage of the space between the margin of the legend and
#'   the border of the plot
#' @param list_only Logical, if TRUE only lists available species without plotting (default: FALSE)
#' @param ... Additional graphical parameters passed to \code{\link[terra]{plot}}
#' @return If \code{list_only = TRUE}, returns a data frame with available species.
#'   Otherwise, invisibly returns a named list of loaded SpatVector objects
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function provides a flexible interface for working with species distribution
#' shapefiles. It can operate in three modes:
#'
#' \strong{Mode 1: List species} - When \code{species_names = NULL} or \code{list_only = TRUE},
#' lists all available species in the specified directory(ies).
#'
#' \strong{Mode 2: Search in list} - When \code{spatvector_list} is provided, searches for
#' species by name in the provided list of SpatVector objects (useful when you've already
#' loaded multiple shapefiles into a list).
#'
#' \strong{Mode 3: Load from directory} - When \code{spatvector_list = NULL}, searches for
#' and loads shapefiles directly from the default directory "out_MCP". Alternative directories are
#' 'out_buffers' for extrapolation done with defined area buffers and 'out' for MSTs.
#'
#' The function automatically handles multiple directories, allowing comparison of different
#' range estimation methods (e.g., MCP vs. buffers). Colors are automatically assigned using
#' the viridis color palette unless custom colors are provided. A legend is generated showing
#' all plotted species.
#'
#' This function is particularly useful for visualizing species distributions from PAE-PCE
#' analyses or other biogeographic studies, where multiple species need to be plotted together
#' with consistent styling. In this case, the file_type argument is 'mst'.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025.
#'   Panbiotracks: software for track analysis. \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'   Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1-17.
#' @seealso
#' \code{\link[terra]{vect}} for loading spatial vector data,
#' \code{\link[terra]{plot}} for plotting spatial data,
#' \code{\link[viridis]{viridis}} for color palettes.
#' @export
#' @examples
#' # Example 1: List available species in a directory (default: "out_MCP")
#' available <- extrapolation_plot_from_shapeFiles(list_only = TRUE)
#' print(available)
#'
#' # Example 2: Plot a single species from directory (default: "out_MCP")
#' # shape file of South America
#' lycipta.asul <- lycipta$polygon
#' plot(lycipta.asul, axes = TRUE, las = 1)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", legend_inset = 0.05) #  Legenda 5% para dentro
#'
#' # Example 3: Plot multiple species with automatic colors (default: "out_MCP")
#' plot(lycipta.asul, axes = TRUE, las = 1)
#' extrapolation_plot_from_shapeFiles(c("E_L_triangulator", "E_L_aceratos", "E_L_imitator")) # legenda cortando o mapa
#'
#' # Example 4: Plot with custom colors
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles(c("E_L_imitator", "E_L_triangulator"),
#'              colors = c("E_L_imitator" = "red", "E_L_triangulator" = "blue"))
#'
#' # Example 5: Plot from buffers directory
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", directory = "out_buffers", file_type = "BUFF")
#'
#' # Example 6: Plot from multiple directories (MCP + buffers)
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", directory = c("out_MCP", "out_buffers"))
#'
#' # Example 7: Search in a pre-loaded list
#' MCP_shapes_n <- list()
#' for(k in 1:length(species_names)){
#'   MCP_shapes_n[[k]] <- terra::vect(paste0('out_MCP/MCP_', species_names[k], '.shp'))
#' }
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", spatvector_list = MCP_shapes_n)
#'
#' # Example 8: Integration with PAE-PCE results
#' pae1 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10), gridView = TRUE, labelGrid = TRUE,
#'    nonHomoplasticSpeciesList = FALSE, N = 10)
#' pae1 # data frame showing the homoplastic and non-homoplastic species (spp) and their grid numbers (grid_n),
#' # from the anterior iterations (PAE-PCE)
#'
#' especies_sinap <- unique(pae1$nonHomoplastic_species[[1]]$spp)
#' plot(lycipta.asul, axes = TRUE, las = 1, main = "PAE-PCE: Individual tracks of synapomorphic species")
#' extrapolation_plot_from_shapeFiles(especies_sinap, alpha = 0.6, directory = 'out', file_type = 'mst')
#'

extrapolation_plot_from_shapeFiles <- function(species_names = NULL,
                                               spatvector_list = NULL,
                                               directory = "out_MCP",
                                               file_type = "MCP",
                                               add = TRUE,
                                               colors = NULL,
                                               alpha = 0.5,
                                               legend = TRUE,
                                               legend_position = "topright",
                                               legend_inset = 0.02,
                                               list_only = FALSE,
                                               ...) {

  # ============================================================================
  # MODE 1: LIST SPECIES ONLY
  # ============================================================================

  if(is.null(species_names) || list_only) {
    all_species_df <- data.frame()

    for(dir in directory) {
      if(!dir.exists(dir)) {
        warning("Directory '", dir, "' does not exist. Skipping.")
        next
      }

      all_shp <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)

      if(length(all_shp) == 0) {
        message("No shapefiles found in directory '", dir, "'")
        next
      }

      pattern <- switch(file_type,
                        "MCP" = "MCP_",
                        "buffer" = "BUFF_",
                        file_type)

      filtered_shp <- grep(pattern, all_shp, value = TRUE)

      if(length(filtered_shp) == 0) {
        filtered_shp <- all_shp
      }

      basenames <- basename(filtered_shp)
      species_names_found <- gsub(paste0("^", pattern), "", basenames)
      species_names_found <- gsub("\\.shp$", "", species_names_found)

      df <- data.frame(
        Directory = dir,
        Species = species_names_found,
        Filename = basenames,
        FullPath = filtered_shp,
        stringsAsFactors = FALSE
      )

      all_species_df <- rbind(all_species_df, df)
    }

    if(nrow(all_species_df) > 0) {
      all_species_df$Index <- seq_len(nrow(all_species_df))
      all_species_df <- all_species_df[, c("Index", "Directory", "Species", "Filename", "FullPath")]
    }

    return(all_species_df)
  }

  # ============================================================================
  # MODE 2: SEARCH IN PROVIDED LIST
  # ============================================================================

  if(!is.null(spatvector_list)) {
    all_sources <- sapply(spatvector_list, terra::sources)
    loaded_vectors <- list()
    found_species <- character()

    for(sp in species_names) {
      # Use exact match: species name must be followed by .shp or end of string
      pattern_exact <- paste0(sp, "(\\.shp)?$")
      index <- grep(pattern_exact, all_sources)

      if(length(index) == 0) {
        warning("No SpatVector found with pattern: ", sp)
        next
      }

      if(length(index) > 1) {
        warning("Multiple matches found for '", sp, "'. Using first one.")
      }

      loaded_vectors[[sp]] <- spatvector_list[[index[1]]]
      found_species <- c(found_species, sp)
      cat("Found:", basename(all_sources[index[1]]), "\n")
    }

    if(length(loaded_vectors) == 0) {
      stop("No species could be found in the provided list")
    }

  } else {

    # ============================================================================
    # MODE 3: LOAD FROM DIRECTORY
    # ============================================================================

    loaded_vectors <- list()
    found_species <- character()

    for(dir in directory) {
      if(!dir.exists(dir)) {
        warning("Directory '", dir, "' does not exist. Skipping.")
        next
      }

      all_shp <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)

      if(length(all_shp) == 0) {
        warning("No shapefiles found in directory '", dir, "'")
        next
      }

      pattern <- switch(file_type,
                        "MCP" = "MCP_",
                        "buffer" = "BUFF_",
                        file_type)

      filtered_shp <- grep(pattern, all_shp, value = TRUE)

      if(length(filtered_shp) == 0) {
        warning("No files matching pattern '", pattern, "' found in '", dir, "'")
        filtered_shp <- all_shp
      }

      for(sp in species_names) {
        # Skip if already found
        if(sp %in% found_species) next

        # Use exact match: species name must be followed by .shp or end of string
        pattern_exact <- paste0(sp, "(\\.shp)?$")
        matching_files <- grep(pattern_exact, filtered_shp, value = TRUE)

        if(length(matching_files) == 0) {
          next
        }

        if(length(matching_files) > 1) {
          warning("Multiple files found for '", sp, "' in '", dir, "'. Using first one.")
        }

        tryCatch({
          v <- terra::vect(matching_files[1])
          loaded_vectors[[sp]] <- v
          found_species <- c(found_species, sp)
          cat("Loaded:", basename(matching_files[1]), "\n")
        }, error = function(e) {
          warning("Failed to load shapefile for '", sp, "': ", e$message)
        })
      }
    }

    if(length(loaded_vectors) == 0) {
      stop("No shapefiles could be loaded for the specified species")
    }

    # Warn about species not found
    not_found <- setdiff(species_names, found_species)
    if(length(not_found) > 0) {
      warning("Species not found: ", paste(not_found, collapse = ", "))
    }
  }

  # ============================================================================
  # PLOT ALL LOADED SPECIES
  # ============================================================================

  # Generate colors if not provided
  if(is.null(colors)) {
    colors <- setNames(viridis::viridis(length(found_species), option = 'D'),
                       found_species)
  } else if(length(colors) != length(found_species)) {
    colors <- setNames(viridis::viridis(length(found_species), option = 'D'),
                       found_species)
  } else if(is.null(names(colors))) {
    colors <- setNames(colors, found_species)
  }

  # Plot each species
  for(sp in found_species) {
    v <- loaded_vectors[[sp]]
    terra::plot(v, add = add, col = colors[sp], alpha = alpha,
                legend = FALSE, ...)
    add <- TRUE  # After first plot, always add
  }

  # Add legend
  if(legend && length(found_species) > 0) {
    par(xpd = TRUE)
    legend(legend_position,
           legend = found_species,
           fill = colors[found_species],
           title = "Species",
           bg = "white",
           cex = 0.7,
           inset = legend_inset,
           box.lwd = 1)
    par(xpd = FALSE)
  }

  invisible(loaded_vectors)
}
