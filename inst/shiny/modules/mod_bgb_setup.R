#' BioGeoBEARS Setup Module UI
#'
#' @param id Module ID
#'
#' @return UI elements
#'
#' @keywords internal
mod_bgb_setup_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Step 4: BioGeoBEARS Configuration"),
    shiny::hr(),

    # Guidance text
    shiny::div(
      class = "guidance-text",
      shiny::p(
        "Configure the BioGeoBEARS models and parameters for ancestral range estimation. ",
        "You can run multiple models for comparison."
      )
    ),

    # Model selection
    shiny::div(
      class = "form-group",
      shiny::h5("Select Models"),
      shiny::checkboxGroupInput(
        ns("models"),
        "Models to run:",
        choices = list(
          "DEC (Dispersal-Extinction-Cladogenesis)" = "DEC",
          "DEC+J (DEC with jump dispersal)" = "DEC+J",
          "DIVALIKE (Dispersal-Vicariance)" = "DIVALIKE",
          "DIVALIKE+J" = "DIVALIKE+J",
          "BAYAREALIKE" = "BAYAREALIKE",
          "BAYAREALIKE+J" = "BAYAREALIKE+J"
        ),
        selected = c("DEC", "DEC+J")
      )
    ),

    shiny::hr(),

    # Parameter bounds
    shiny::div(
      class = "form-group",
      shiny::h5("Parameter Bounds"),

      shiny::div(
        class = "row",
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("d_min"),
            "d (dispersal) minimum:",
            value = 0.0001,
            min = 0,
            step = 0.0001
          )
        ),
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("d_max"),
            "d maximum:",
            value = 100,
            min = 0,
            step = 1
          )
        )
      ),

      shiny::div(
        class = "row",
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("e_min"),
            "e (extinction) minimum:",
            value = 0.0001,
            min = 0,
            step = 0.0001
          )
        ),
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("e_max"),
            "e maximum:",
            value = 100,
            min = 0,
            step = 1
          )
        )
      ),

      shiny::div(
        class = "row",
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("j_min"),
            "j (jump dispersal) minimum:",
            value = 0.0001,
            min = 0,
            step = 0.0001
          )
        ),
        shiny::div(
          class = "col-md-6",
          shiny::numericInput(
            ns("j_max"),
            "j maximum:",
            value = 1,
            min = 0,
            step = 0.01
          )
        )
      )
    ),

    shiny::hr(),

    # Advanced options
    shiny::div(
      class = "form-group",
      shiny::h5("Advanced Options"),

      shiny::checkboxInput(
        ns("time_stratification"),
        "Use time stratification",
        value = FALSE
      ),

      shiny::conditionalPanel(
        condition = "input.time_stratification",
        ns = ns,
        shiny::fileInput(
          ns("time_strat_file"),
          "Time stratification file",
          accept = c(".txt", ".csv")
        )
      ),

      shiny::checkboxInput(
        ns("distance_matrix"),
        "Use distance matrix (x parameter)",
        value = FALSE
      ),

      shiny::conditionalPanel(
        condition = "input.distance_matrix",
        ns = ns,
        shiny::fileInput(
          ns("distance_file"),
          "Distance matrix file",
          accept = c(".txt", ".csv")
        )
      ),

      shiny::checkboxInput(
        ns("dispersal_multiplier"),
        "Use dispersal multiplier matrix",
        value = FALSE
      ),

      shiny::conditionalPanel(
        condition = "input.dispersal_multiplier",
        ns = ns,
        shiny::fileInput(
          ns("dispersal_file"),
          "Dispersal multiplier file",
          accept = c(".txt", ".csv")
        )
      )
    ),

    shiny::hr(),

    # Optimization options
    shiny::div(
      class = "form-group",
      shiny::h5("Optimization Settings"),

      shiny::selectInput(
        ns("optimizer"),
        "Optimizer:",
        choices = list(
          "GenSA (default)" = "GenSA",
          "optim" = "optim"
        ),
        selected = "GenSA"
      ),

      shiny::numericInput(
        ns("max_iterations_optim"),
        "Maximum iterations:",
        value = 1000,
        min = 100,
        step = 100
      )
    ),

    shiny::hr(),

    # Action button
    shiny::actionButton(
      ns("setup_complete"),
      "Proceed to Analysis",
      class = "btn-primary btn-lg",
      width = "100%"
    ),

    # Status message
    shiny::uiOutput(ns("status_message"))
  )
}

#' BioGeoBEARS Setup Module Server
#'
#' @param id Module ID
#' @param pres_abs_matrix Reactive presence-absence matrix
#' @param tree Reactive phylogenetic tree
#'
#' @return Reactive values with configuration
#'
#' @keywords internal
mod_bgb_setup_server <- function(id, pres_abs_matrix, tree) {
  shiny::moduleServer(id, function(input, output, session) {
    # Reactive values
    config <- shiny::reactiveValues(
      models = NULL,
      parameters = NULL,
      advanced_options = NULL,
      optimizer = NULL,
      completed = FALSE,
      message = ""
    )

    # Setup complete button
    shiny::observeEvent(input$setup_complete, {
      tryCatch({
        if (length(input$models) == 0) {
          stop("Please select at least one model")
        }

        # Store configuration
        config$models <- input$models

        config$parameters <- list(
          d = c(min = input$d_min, max = input$d_max),
          e = c(min = input$e_min, max = input$e_max),
          j = c(min = input$j_min, max = input$j_max)
        )

        config$advanced_options <- list(
          time_stratification = input$time_stratification,
          distance_matrix = input$distance_matrix,
          dispersal_multiplier = input$dispersal_multiplier
        )

        config$optimizer <- list(
          method = input$optimizer,
          max_iterations = input$max_iterations_optim
        )

        config$completed <- TRUE
        config$message <- "Configuration saved successfully!"
      }, error = function(e) {
        config$completed <- FALSE
        config$message <- paste("Error:", e$message)
      })
    })

    # Status message
    output$status_message <- shiny::renderUI({
      if (config$message == "") {
        return(NULL)
      }

      if (config$completed) {
        shiny::div(
          class = "alert alert-success",
          config$message
        )
      } else {
        shiny::div(
          class = "alert alert-danger",
          config$message
        )
      }
    })

    # Return reactive configuration
    shiny::reactive({
      list(
        models = config$models,
        parameters = config$parameters,
        advanced_options = config$advanced_options,
        optimizer = config$optimizer,
        completed = config$completed
      )
    })
  })
}
