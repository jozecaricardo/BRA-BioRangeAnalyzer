#' Visualization Functions for BioGeoBEARS Shiny
#'
#' Functions for creating interactive maps, phylogenetic trees, and model comparison plots
#'

#' Plot Phylogenetic Tree
#'
#' Create a phylogenetic tree visualization with branch lengths
#'
#' @param tree phylo object
#' @param title Title for the plot
#'
#' @return NULL (creates plot as side effect)
#'
#' @export
plot_phylogenetic_tree <- function(tree, title = "Phylogenetic Tree") {
  
  plot(tree, main = title, cex = 0.8)
  
  # Add axis information
  axisPhylo()
  
  invisible(NULL)
}

#' Plot Species Distribution Map
#'
#' Create an interactive map showing species occurrences with optional extrapolation polygons
#'
#' @param occurrence_data Data frame with species, longitude, latitude columns
#' @param geometry Optional spatial geometry object (polygons from extrapolation)
#' @param method Optional method name ("BUFF", "MPC", "MST") for labeling
#'
#' @return leaflet map object
#'
#' @export
plot_distribution_map <- function(occurrence_data, geometry = NULL, method = NULL) {
  
  # Create base map
  m <- leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::setView(lng = mean(occurrence_data$long, na.rm = TRUE),
                     lat = mean(occurrence_data$lat, na.rm = TRUE),
                     zoom = 4)
  
  # Add occurrence points
  species_unique <- unique(occurrence_data$spp)
  n_species <- length(species_unique)
  
  # Choose appropriate color palette based on number of species
  if (n_species <= 3) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set1")
  } else if (n_species <= 9) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set1")
  } else if (n_species <= 12) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set3")
  } else {
    # For more than 12 species, use a continuous palette
    colors <- grDevices::hcl.colors(n_species, "Spectral")
  }
  
  for (i in seq_along(species_unique)) {
    sp <- species_unique[i]
    sp_data <- occurrence_data[occurrence_data$spp == sp, ]
    
    m <- m %>%
      leaflet::addCircleMarkers(
        data = sp_data,
        lng = ~long,
        lat = ~lat,
        radius = 5,
        color = colors[i],
        fillOpacity = 0.7,
        popup = ~paste(spp, "<br>Lon:", round(long, 2), "<br>Lat:", round(lat, 2)),
        group = sp
      )
  }
  
  # Add extrapolation polygons if provided
  if (!is.null(geometry)) {
    m <- biogeoshiny::add_extrapolation_polygons(
      map = m,
      geometry = geometry,
      method = ifelse(is.null(method), "Extrapolation", method),
      color = "blue",
      opacity = 0.5
    )
  }
  
  # Prepare layer groups for control
  layer_groups <- c(species_unique)
  if (!is.null(geometry)) {
    layer_groups <- c(layer_groups, ifelse(is.null(method), "Extrapolation", method))
  }
  
  # Add layer control
  m <- m %>%
    leaflet::addLayersControl(
      overlayGroups = layer_groups,
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )
  
  m
}

#' Plot Extrapolation Polygons on Map
#'
#' Add extrapolation polygons (buffers, convex hulls, MST) to a leaflet map
#'
#' @param map leaflet map object
#' @param geometry Spatial geometry object (SpatVector or sf)
#' @param species_names Character vector of species names
#'
#' @return Updated leaflet map object
#'
#' @export
plot_extrapolation_polygons <- function(map, geometry, species_names) {
  
  if (is.null(geometry)) {
    return(map)
  }
  
  tryCatch({
    # Handle different geometry types
    if (inherits(geometry, "SpatVector")) {
      # Convert terra SpatVector to sf
      geometry <- sf::st_as_sf(geometry)
    }
    
    if (inherits(geometry, "sf") || inherits(geometry, "data.frame")) {
      # Add polygons to map
      colors <- RColorBrewer::brewer.pal(min(length(species_names), 12), "Set1")
      
      for (i in seq_along(species_names)) {
        if (i <= nrow(geometry)) {
          map <- map %>%
            leaflet::addPolygons(
              data = geometry[i, ],
              color = colors[i],
              fillColor = colors[i],
              fillOpacity = 0.3,
              weight = 2,
              popup = paste("Species:", species_names[i]),
              group = paste("Polygon:", species_names[i])
            )
        }
      }
    }
    
    map
  }, error = function(e) {
    warning("Could not add polygons to map: ", e$message)
    map
  })
}

#' Plot Model Comparison by AICc
#'
#' Create an interactive bar plot comparing models by AICc
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_model_comparison <- function(comparison_table) {
  
  # Sort by AICc
  comparison_table <- comparison_table[order(comparison_table$AICc), ]
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~reorder(Model, AICc),
    y = ~AICc,
    type = "bar",
    marker = list(color = ~DeltaAICc),
    text = ~paste("AICc:", round(AICc, 2), "<br>ΔAICc:", round(DeltaAICc, 2)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = "Model Comparison by AICc",
      xaxis = list(title = "Model"),
      yaxis = list(title = "AICc"),
      hovermode = "closest"
    )
  
  p
}

#' Plot Model Likelihood
#'
#' Create a plot showing log-likelihood values for different models
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_model_likelihood <- function(comparison_table) {
  
  # Sort by LnL
  comparison_table <- comparison_table[order(comparison_table$LnL, decreasing = TRUE), ]
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~reorder(Model, LnL),
    y = ~LnL,
    type = "bar",
    marker = list(color = "steelblue"),
    text = ~paste("LnL:", round(LnL, 2)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = "Model Log-Likelihood Comparison",
      xaxis = list(title = "Model"),
      yaxis = list(title = "Log-Likelihood"),
      hovermode = "closest"
    )
  
  p
}

#' Plot Parameter Comparison
#'
#' Create a plot comparing parameter values across models
#'
#' @param comparison_table Data frame with model comparison results
#' @param parameter Parameter to plot ("d", "e", or "j")
#'
#' @return plotly object
#'
#' @export
plot_parameter_comparison <- function(comparison_table, parameter = "d") {
  
  if (!(parameter %in% c("d", "e", "j"))) {
    stop("Parameter must be 'd', 'e', or 'j'")
  }
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~Model,
    y = as.formula(paste("~", parameter)),
    type = "bar",
    marker = list(color = "coral"),
    text = ~paste(parameter, "=", get(parameter)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = paste("Parameter", parameter, "Comparison"),
      xaxis = list(title = "Model"),
      yaxis = list(title = paste("Parameter", parameter)),
      hovermode = "closest"
    )
  
  p
}

#' Plot Parameter Heatmap
#'
#' Create a heatmap showing parameter values across models
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_parameter_heatmap <- function(comparison_table) {
  
  # Extract parameter values
  params_matrix <- as.matrix(comparison_table[, c("d", "e")])
  rownames(params_matrix) <- comparison_table$Model
  
  # Create heatmap
  p <- plotly::plot_ly(
    z = params_matrix,
    x = colnames(params_matrix),
    y = rownames(params_matrix),
    type = "heatmap",
    colorscale = "Viridis"
  ) %>%
    plotly::layout(
      title = "Parameter Heatmap",
      xaxis = list(title = "Parameter"),
      yaxis = list(title = "Model")
    )
  
  p
}
