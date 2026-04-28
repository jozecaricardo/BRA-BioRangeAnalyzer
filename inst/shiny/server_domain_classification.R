# ===== DOMAIN CLASSIFICATION TAB HANDLERS =====

# Reactive storage for domain data
domain_store <- reactiveValues(
  shapefile = NULL,
  domain_column = NULL,
  classification = NULL
)

# Load domain shapefile
observeEvent(input$load_domain_shapefile, {
  tryCatch({
    if (is.null(input$domain_shapefile)) {
      output$domain_shapefile_status <- renderPrint({
        cat("Error: Please select a shapefile\n")
      })
    } else {
      # Read the shapefile using terra
      shapefile <- terra::vect(input$domain_shapefile$datapath)
      
      # Set CRS if needed
      if (is.na(terra::crs(shapefile))) {
        terra::crs(shapefile) <- "+proj=longlat +datum=WGS84"
      }
      
      domain_store$shapefile <- shapefile
      domain_store$domain_column <- input$domain_column
      
      output$domain_shapefile_status <- renderPrint({
        cat("✓ Shapefile loaded successfully!\n")
        cat("Number of features:", nrow(shapefile), "\n")
        cat("CRS:", terra::crs(shapefile), "\n")
        cat("Columns:", paste(names(shapefile), collapse = ", "), "\n")
      })
    }
  }, error = function(e) {
    output$domain_shapefile_status <- renderPrint({
      cat("Error loading shapefile:\n")
      cat(e$message, "\n")
    })
  })
})

# Classify occurrences by domain
observeEvent(input$classify_by_domain, {
  if (is.null(data_store$occurrence)) {
    output$classification_status <- renderPrint({
      cat("Error: Please load occurrence data first\n")
    })
  } else if (is.null(domain_store$shapefile)) {
    output$classification_status <- renderPrint({
      cat("Error: Please load a domain shapefile first\n")
    })
  } else {
    tryCatch({
      # Prepare occurrence data
      occ_data <- data_store$occurrence
      names(occ_data) <- c("species", "long", "lat")
      
      # Convert shapefile to sf for SpGeoCod
      shapefile_sf <- sf::st_as_sf(domain_store$shapefile)
      
      # Convert occurrence data to sf
      occ_sf <- sf::st_as_sf(occ_data, coords = c("long", "lat"), crs = 4326)
      
      # Perform spatial join
      occ_classified <- sf::st_join(occ_sf, shapefile_sf)
      
      # Store classification
      domain_store$classification <- occ_classified
      
      # Count species per domain
      domain_col <- domain_store$domain_column
      
      if (domain_col %in% names(occ_classified)) {
        species_per_domain <- occ_classified %>%
          sf::st_drop_geometry() %>%
          dplyr::group_by(!!rlang::sym(domain_col)) %>%
          dplyr::summarise(n_species = n_distinct(species), .groups = "drop")
        
        output$classification_status <- renderPrint({
          cat("✓ Classification completed!\n")
          cat("Total occurrences:", nrow(occ_classified), "\n")
          cat("Unique species:", n_distinct(occ_classified$species), "\n")
          cat("\nSpecies per domain:\n")
          print(species_per_domain)
        })
      } else {
        output$classification_status <- renderPrint({
          cat("Error: Domain column '", domain_col, "' not found in shapefile\n")
          cat("Available columns:", paste(names(shapefile_sf), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$classification_status <- renderPrint({
        cat("Error during classification:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# Domain map visualization
output$domain_map <- leaflet::renderLeaflet({
  if (is.null(domain_store$shapefile)) {
    leaflet::leaflet() %>%
      leaflet::addTiles() %>%
      leaflet::setView(lng = -60, lat = -15, zoom = 4)
  } else {
    # Convert to sf for leaflet
    shapefile_sf <- sf::st_as_sf(domain_store$shapefile)
    
    # Create color palette
    n_domains <- nrow(shapefile_sf)
    colors <- if (n_domains <= 9) {
      RColorBrewer::brewer.pal(n_domains, "Set1")
    } else {
      grDevices::hcl.colors(n_domains, "Spectral")
    }
    
    # Create map
    map <- leaflet::leaflet(shapefile_sf) %>%
      leaflet::addTiles() %>%
      leaflet::fitBounds(
        sf::st_bbox(shapefile_sf)["xmin"],
        sf::st_bbox(shapefile_sf)["ymin"],
        sf::st_bbox(shapefile_sf)["xmax"],
        sf::st_bbox(shapefile_sf)["ymax"]
      )
    
    # Add polygons
    for (i in 1:nrow(shapefile_sf)) {
      map <- map %>%
        leaflet::addPolygons(
          data = shapefile_sf[i, ],
          color = colors[i],
          fillOpacity = 0.5,
          popup = paste(names(shapefile_sf), shapefile_sf[i, ], sep = ": ")
        )
    }
    
    map
  }
})

# Species richness by domain plot
output$richness_by_domain_plot <- renderPlot({
  if (is.null(domain_store$classification)) {
    plot(1, type = "n", main = "Species Richness by Domain", xlab = "", ylab = "")
    text(1, 1, "Classify occurrences first", cex = 1.5, col = "red")
  } else {
    tryCatch({
      domain_col <- domain_store$domain_column
      
      # Count species per domain
      richness_data <- domain_store$classification %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(!!rlang::sym(domain_col)) %>%
        dplyr::summarise(n_species = n_distinct(species), .groups = "drop") %>%
        dplyr::arrange(desc(n_species))
      
      # Create bar plot
      barplot(richness_data$n_species,
             names.arg = richness_data[[domain_col]],
             main = "Species Richness by Domain",
             xlab = domain_col,
             ylab = "Number of Species",
             col = "steelblue",
             las = 2)
    }, error = function(e) {
      plot(1, type = "n", main = "Error", xlab = "", ylab = "")
      text(1, 1, paste("Error:", e$message), cex = 1, col = "red")
    })
  }
})

# Export to NEXUS
observeEvent(input$export_domain_nexus, {
  if (is.null(domain_store$classification)) {
    output$domain_nexus_status <- renderPrint({
      cat("Error: Please classify occurrences first\n")
    })
  } else {
    tryCatch({
      domain_col <- domain_store$domain_column
      
      # Create presence-absence matrix by domain
      occ_data <- domain_store$classification %>%
        sf::st_drop_geometry()
      
      # Create matrix: domains as rows, species as columns
      domains <- unique(occ_data[[domain_col]])
      species <- unique(occ_data$species)
      
      pres_abs <- matrix(0, nrow = length(domains), ncol = length(species),
                        dimnames = list(domains, species))
      
      for (i in 1:nrow(occ_data)) {
        domain <- occ_data[i, domain_col][[1]]
        sp <- occ_data[i, "species"][[1]]
        pres_abs[domain, sp] <- 1
      }
      
      # Create output directory
      dir.create("out_domains", showWarnings = FALSE)
      
      # Write NEXUS file
      nexus_file <- file.path("out_domains", "domains_classification.nex")
      
      # Simple NEXUS format
      write("#NEXUS", nexus_file)
      write("begin data;", nexus_file, append = TRUE)
      write(paste("dimensions ntax=", nrow(pres_abs), " nchar=", ncol(pres_abs), ";", sep = ""),
           nexus_file, append = TRUE)
      write("format datatype=standard symbols=\"01\";", nexus_file, append = TRUE)
      write("matrix", nexus_file, append = TRUE)
      
      for (i in 1:nrow(pres_abs)) {
        line <- paste(rownames(pres_abs)[i], paste(pres_abs[i, ], collapse = ""))
        write(line, nexus_file, append = TRUE)
      }
      
      write(";", nexus_file, append = TRUE)
      write("end;", nexus_file, append = TRUE)
      
      output$domain_nexus_status <- renderPrint({
        cat("✓ NEXUS file exported!\n")
        cat("File:", nexus_file, "\n")
        cat("Dimensions:", nrow(pres_abs), "domains x", ncol(pres_abs), "species\n")
      })
    }, error = function(e) {
      output$domain_nexus_status <- renderPrint({
        cat("Error exporting NEXUS:\n")
        cat(e$message, "\n")
      })
    })
  }
})

# Export to BioGeoBEARS
observeEvent(input$export_domain_biogeobears, {
  if (is.null(domain_store$classification)) {
    output$domain_biogeobears_status <- renderPrint({
      cat("Error: Please classify occurrences first\n")
    })
  } else {
    tryCatch({
      domain_col <- domain_store$domain_column
      
      # Create presence-absence matrix by domain
      occ_data <- domain_store$classification %>%
        sf::st_drop_geometry()
      
      # Create matrix: domains as rows, species as columns
      domains <- unique(occ_data[[domain_col]])
      species <- unique(occ_data$species)
      
      pres_abs <- matrix(0, nrow = length(domains), ncol = length(species),
                        dimnames = list(domains, species))
      
      for (i in 1:nrow(occ_data)) {
        domain <- occ_data[i, domain_col][[1]]
        sp <- occ_data[i, "species"][[1]]
        pres_abs[domain, sp] <- 1
      }
      
      # Create output directory
      dir.create("out_domains", showWarnings = FALSE)
      
      # Write BioGeoBEARS format
      biogeobears_file <- file.path("out_domains", "domains_BioGeoBEARS.txt")
      
      # Write as tab-separated file
      write.table(pres_abs, biogeobears_file, sep = "\t", quote = FALSE)
      
      output$domain_biogeobears_status <- renderPrint({
        cat("✓ BioGeoBEARS file exported!\n")
        cat("File:", biogeobears_file, "\n")
        cat("Dimensions:", nrow(pres_abs), "domains x", ncol(pres_abs), "species\n")
      })
    }, error = function(e) {
      output$domain_biogeobears_status <- renderPrint({
        cat("Error exporting BioGeoBEARS:\n")
        cat(e$message, "\n")
      })
    })
  }
})

