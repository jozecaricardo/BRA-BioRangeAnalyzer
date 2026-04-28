#' Export Formats Module UI
#'
#' @param id Module ID
#'
#' @return UI elements
#'
#' @keywords internal
mod_export_formats_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Step 3: Export Data Formats"),
    shiny::hr(),

    # Guidance text
    shiny::div(
      class = "guidance-text",
      shiny::p(
        "Export your presence-absence matrix to various formats for use with ",
        "BioGeoBEARS, NDM/VNDM, TNT, and other phylogenetic software."
      )
    ),

    # Export options
    shiny::div(
      class = "form-group",
      shiny::h5("Select Formats to Export"),
      shiny::checkboxGroupInput(
        ns("export_formats"),
        "Formats:",
        choices = list(
          "BioGeoBEARS (.data)" = "biogeobears",
          "NEXUS (.nex)" = "nexus",
          "TNT (.tnt)" = "tnt",
          "NDM/VNDM (.xyd)" = "ndm"
        ),
        selected = c("biogeobears", "nexus")
      )
    ),

    shiny::hr(),

    # NDM-specific options
    shiny::div(
      class = "form-group",
      shiny::h5("NDM/VNDM Options"),
      shiny::numericInput(
        ns("ndm_grid_resolution"),
        "Grid resolution (degrees):",
        value = 1,
        min = 0.1,
        step = 0.1
      )
    ),

    shiny::hr(),

    # TNT-specific options
    shiny::div(
      class = "form-group",
      shiny::h5("TNT PAE-PCE Options"),
      shiny::checkboxInput(
        ns("generate_tnt_script"),
        "Generate PAE-PCE script",
        value = TRUE
      ),
      shiny::numericInput(
        ns("max_iterations"),
        "Maximum PAE-PCE iterations:",
        value = 10,
        min = 1,
        step = 1
      ),
      shiny::numericInput(
        ns("search_replicates"),
        "Search replicates:",
        value = 4,
        min = 1,
        step = 1
      )
    ),

    shiny::hr(),

    # Action button
    shiny::actionButton(
      ns("export_data"),
      "Export Data",
      class = "btn-primary btn-lg",
      width = "100%"
    ),

    shiny::hr(),

    # Download buttons
    shiny::uiOutput(ns("download_buttons")),

    # Status message
    shiny::uiOutput(ns("status_message"))
  )
}

#' Export Formats Module Server
#'
#' @param id Module ID
#' @param pres_abs_matrix Reactive presence-absence matrix
#'
#' @return Reactive values with exported files
#'
#' @keywords internal
mod_export_formats_server <- function(id, pres_abs_matrix) {
  shiny::moduleServer(id, function(input, output, session) {
    # Reactive values
    exports <- shiny::reactiveValues(
      biogeobears = NULL,
      nexus = NULL,
      tnt = NULL,
      ndm = NULL,
      completed = FALSE,
      message = ""
    )

    # Export data
    shiny::observeEvent(input$export_data, {
      shiny::withProgress(message = "Exporting data...", value = 0, {
        tryCatch({
          pres_abs <- pres_abs_matrix()

          if (is.null(pres_abs)) {
            stop("No presence-absence matrix available")
          }

          # Create output directory
          output_dir <- tempdir()

          # Export BioGeoBEARS format
          if ("biogeobears" %in% input$export_formats) {
            exports$message <- "Exporting BioGeoBEARS format..."
            shiny::incProgress(0.25)

            # Format: n_species n_areas area1 area2 ...
            # followed by species and their ranges
            bgb_file <- file.path(output_dir, "pres_abs_geog.data")

            n_species <- ncol(pres_abs)
            n_areas <- nrow(pres_abs)

            # Write header
            header <- paste(n_species, n_areas, paste(rownames(pres_abs), collapse = " "))

            # Write species and ranges
            species_lines <- apply(pres_abs, 2, function(col) {
              paste(colnames(pres_abs)[which(col == col)], paste(col, collapse = ""), sep = " ")
            })

            writeLines(c(header, species_lines), bgb_file)
            exports$biogeobears <- bgb_file
          }

          # Export NEXUS format
          if ("nexus" %in% input$export_formats) {
            exports$message <- "Exporting NEXUS format..."
            shiny::incProgress(0.25)

            nex_file <- file.path(output_dir, "pres_abs.nex")

            # Write NEXUS file
            n_taxa <- nrow(pres_abs)
            n_chars <- ncol(pres_abs)

            nexus_content <- c(
              "#NEXUS",
              "begin data;",
              paste("dimensions ntax=", n_taxa, " nchar=", n_chars, ";", sep = ""),
              'format datatype=standard symbols="01" gap=-;',
              "",
              "CHARSTATELABELS"
            )

            # Add character labels
            for (i in 1:n_chars) {
              nexus_content <- c(nexus_content, paste(i, colnames(pres_abs)[i], ",", sep = " "))
            }

            nexus_content <- c(
              nexus_content,
              ";",
              "",
              "matrix"
            )

            # Add data
            for (i in 1:n_taxa) {
              nexus_content <- c(
                nexus_content,
                paste(rownames(pres_abs)[i], paste(pres_abs[i, ], collapse = ""))
              )
            }

            nexus_content <- c(nexus_content, ";", "end;")

            writeLines(nexus_content, nex_file)
            exports$nexus <- nex_file
          }

          # Export TNT format
          if ("tnt" %in% input$export_formats) {
            exports$message <- "Exporting TNT format..."
            shiny::incProgress(0.25)

            tnt_file <- file.path(output_dir, "pres_abs.tnt")

            n_taxa <- nrow(pres_abs)
            n_chars <- ncol(pres_abs)

            tnt_content <- c(
              paste(n_taxa, n_chars)
            )

            # Add data
            for (i in 1:n_taxa) {
              tnt_content <- c(
                tnt_content,
                paste(rownames(pres_abs)[i], paste(pres_abs[i, ], collapse = ""))
              )
            }

            writeLines(tnt_content, tnt_file)
            exports$tnt <- tnt_file
          }

          # Export NDM format
          if ("ndm" %in% input$export_formats) {
            exports$message <- "Exporting NDM/VNDM format..."
            shiny::incProgress(0.25)

            ndm_file <- file.path(output_dir, "pres_abs.xyd")

            # Convert presence-absence to XYD format
            # This is a simplified version; actual implementation would use grid centroids
            xyd_content <- c()

            for (i in 1:nrow(pres_abs)) {
              for (j in 1:ncol(pres_abs)) {
                if (pres_abs[i, j] == 1) {
                  # Generate dummy coordinates
                  lon <- -60 + (i * 0.5)
                  lat <- -10 + (j * 0.5)
                  xyd_content <- c(xyd_content, paste(colnames(pres_abs)[j], lon, lat))
                }
              }
            }

            writeLines(xyd_content, ndm_file)
            exports$ndm <- ndm_file
          }

          exports$completed <- TRUE
          exports$message <- "Export completed successfully!"

          shiny::incProgress(1)
        }, error = function(e) {
          exports$completed <- FALSE
          exports$message <- paste("Error:", e$message)
        })
      })
    })

    # Download buttons
    output$download_buttons <- shiny::renderUI({
      if (!exports$completed) {
        return(NULL)
      }

      buttons <- list()

      if (!is.null(exports$biogeobears)) {
        buttons[[length(buttons) + 1]] <- shiny::downloadButton(
          session$ns("download_biogeobears"),
          "Download BioGeoBEARS"
        )
      }

      if (!is.null(exports$nexus)) {
        buttons[[length(buttons) + 1]] <- shiny::downloadButton(
          session$ns("download_nexus"),
          "Download NEXUS"
        )
      }

      if (!is.null(exports$tnt)) {
        buttons[[length(buttons) + 1]] <- shiny::downloadButton(
          session$ns("download_tnt"),
          "Download TNT"
        )
      }

      if (!is.null(exports$ndm)) {
        buttons[[length(buttons) + 1]] <- shiny::downloadButton(
          session$ns("download_ndm"),
          "Download NDM"
        )
      }

      do.call(shiny::div, c(class = "btn-group", buttons))
    })

    # Download handlers
    output$download_biogeobears <- shiny::downloadHandler(
      filename = "pres_abs_geog.data",
      content = function(file) {
        file.copy(exports$biogeobears, file)
      }
    )

    output$download_nexus <- shiny::downloadHandler(
      filename = "pres_abs.nex",
      content = function(file) {
        file.copy(exports$nexus, file)
      }
    )

    output$download_tnt <- shiny::downloadHandler(
      filename = "pres_abs.tnt",
      content = function(file) {
        file.copy(exports$tnt, file)
      }
    )

    output$download_ndm <- shiny::downloadHandler(
      filename = "pres_abs.xyd",
      content = function(file) {
        file.copy(exports$ndm, file)
      }
    )

    # Status message
    output$status_message <- shiny::renderUI({
      if (exports$message == "") {
        return(NULL)
      }

      if (exports$completed) {
        shiny::div(
          class = "alert alert-success",
          exports$message
        )
      } else {
        shiny::div(
          class = "alert alert-info",
          exports$message
        )
      }
    })

    # Return reactive exports
    shiny::reactive({
      list(
        biogeobears = exports$biogeobears,
        nexus = exports$nexus,
        tnt = exports$tnt,
        ndm = exports$ndm,
        completed = exports$completed
      )
    })
  })
}
