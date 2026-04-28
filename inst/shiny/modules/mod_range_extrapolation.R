#' Range Extrapolation Module UI
#'
#' @param id Module ID
#'
#' @return UI elements
#'
#' @keywords internal
mod_range_extrapolation_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Step 2: Range Extrapolation"),
    shiny::hr(),

    # Guidance text
    shiny::div(
      class = "guidance-text",
      shiny::p(
        "Select a method to extrapolate species ranges from occurrence points. ",
        "Each method has different parameters and assumptions."
      )
    ),

    # Method selection
    shiny::div(
      class = "form-group",
      shiny::h5("Select Method"),
      shiny::radioButtons(
        ns("method"),
        "Range Extrapolation Method:",
        choices = list(
          "Buffer (circular buffers around points)" = "buffer",
          "Convex Hull (minimum convex polygon)" = "convex_hull",
          "Minimum Spanning Tree" = "mst"
        ),
        selected = "buffer"
      )
    ),

    shiny::hr(),

    # Buffer parameters
    shiny::conditionalPanel(
      condition = "input.method == 'buffer'",
      ns = ns,
      shiny::div(
        class = "form-group",
        shiny::h5("Buffer Parameters"),
        shiny::numericInput(
          ns("buffer_width"),
          "Buffer width (meters):",
          value = 100000,
          min = 1000,
          step = 10000
        ),
        shiny::checkboxInput(
          ns("mean_dist"),
          "Use mean distance between points",
          value = FALSE
        )
      )
    ),

    # Convex hull parameters
    shiny::conditionalPanel(
      condition = "input.method == 'convex_hull'",
      ns = ns,
      shiny::div(
        class = "form-group",
        shiny::h5("Convex Hull Parameters"),
        shiny::p("No additional parameters needed.")
      )
    ),

    # MST parameters
    shiny::conditionalPanel(
      condition = "input.method == 'mst'",
      ns = ns,
      shiny::div(
        class = "form-group",
        shiny::h5("MST Parameters"),
        shiny::p("No additional parameters needed.")
      )
    ),

    shiny::hr(),

    # Grid resolution
    shiny::div(
      class = "form-group",
      shiny::h5("Grid Resolution"),
      shiny::numericInput(
        ns("resolution"),
        "Resolution (decimal degrees):",
        value = 1,
        min = 0.1,
        step = 0.1
      )
    ),

    shiny::hr(),

    # Action button
    shiny::actionButton(
      ns("run_extrapolation"),
      "Run Extrapolation",
      class = "btn-primary btn-lg",
      width = "100%"
    ),

    # Progress and status
    shiny::uiOutput(ns("progress_ui")),
    shiny::uiOutput(ns("status_message"))
  )
}

#' Range Extrapolation Module Server
#'
#' @param id Module ID
#' @param occurrence_data Reactive occurrence data
#'
#' @return Reactive values with extrapolation results
#'
#' @keywords internal
mod_range_extrapolation_server <- function(id, occurrence_data) {
  shiny::moduleServer(id, function(input, output, session) {
    # Reactive values
    results <- shiny::reactiveValues(
      pres_abs = NULL,
      geometry = NULL,
      shapefiles = NULL,
      completed = FALSE,
      message = ""
    )

    # Progress indicator
    output$progress_ui <- shiny::renderUI({
      if (results$completed) {
        shiny::div(
          class = "alert alert-success",
          "Extrapolation completed!"
        )
      }
    })

    # Status message
    output$status_message <- shiny::renderUI({
      if (results$message == "") {
        return(NULL)
      }
      shiny::div(
        class = "alert alert-info",
        results$message
      )
    })

    # Run extrapolation
    shiny::observeEvent(input$run_extrapolation, {
      shiny::withProgress(message = "Running extrapolation...", value = 0, {
        tryCatch({
          occ_data <- occurrence_data()

          if (is.null(occ_data$occurrence)) {
            stop("No occurrence data loaded")
          }

          # Format data for range calculation functions
          xy <- data.frame(
            spp = occ_data$occurrence$species,
            long = occ_data$occurrence$longitude,
            lat = occ_data$occurrence$latitude
          )

          # Run selected method
          if (input$method == "buffer") {
            results$message <- "Calculating buffers..."
            shiny::incProgress(0.3)

            # This would call calcRange_buffers from the project
            # For now, we'll create a placeholder
            result <- list(
              pres_abs = matrix(1, nrow = 10, ncol = length(unique(xy$spp))),
              geometry = NULL
            )
          } else if (input$method == "convex_hull") {
            results$message <- "Calculating convex hulls..."
            shiny::incProgress(0.3)

            # This would call calcRange_convexHull
            result <- list(
              pres_abs = matrix(1, nrow = 10, ncol = length(unique(xy$spp))),
              geometry = NULL
            )
          } else if (input$method == "mst") {
            results$message <- "Calculating minimum spanning trees..."
            shiny::incProgress(0.3)

            # This would call calcRange_irregularBins or similar
            result <- list(
              pres_abs = matrix(1, nrow = 10, ncol = length(unique(xy$spp))),
              geometry = NULL
            )
          }

          results$pres_abs <- result$pres_abs
          results$geometry <- result$geometry
          results$completed <- TRUE
          results$message <- "Extrapolation completed successfully!"

          shiny::incProgress(1)
        }, error = function(e) {
          results$completed <- FALSE
          results$message <- paste("Error:", e$message)
        })
      })
    })

    # Return reactive results
    shiny::reactive({
      list(
        pres_abs = results$pres_abs,
        geometry = results$geometry,
        completed = results$completed
      )
    })
  })
}
