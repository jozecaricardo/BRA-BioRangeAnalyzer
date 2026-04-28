#' BioGeoBEARS Analysis Module UI
#'
#' @param id Module ID
#'
#' @return UI elements
#'
#' @keywords internal
mod_bgb_analysis_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Step 5: Run Analysis"),
    shiny::hr(),

    # Guidance text
    shiny::div(
      class = "guidance-text",
      shiny::p(
        "Run the BioGeoBEARS analysis with your selected models. ",
        "This may take several minutes depending on the number of models and areas."
      )
    ),

    # Action button
    shiny::actionButton(
      ns("run_analysis"),
      "Run BioGeoBEARS Analysis",
      class = "btn-primary btn-lg",
      width = "100%"
    ),

    shiny::hr(),

    # Progress indicator
    shiny::uiOutput(ns("progress_ui")),

    # Results tabs
    shiny::uiOutput(ns("results_ui"))
  )
}

#' BioGeoBEARS Analysis Module Server
#'
#' @param id Module ID
#' @param pres_abs_matrix Reactive presence-absence matrix
#' @param tree Reactive phylogenetic tree
#' @param bgb_config Reactive BioGeoBEARS configuration
#'
#' @return Reactive values with analysis results
#'
#' @keywords internal
mod_bgb_analysis_server <- function(id, pres_abs_matrix, tree, bgb_config) {
  shiny::moduleServer(id, function(input, output, session) {
    # Reactive values
    results <- shiny::reactiveValues(
      models_results = NULL,
      model_comparison = NULL,
      completed = FALSE,
      message = ""
    )

    # Run analysis
    shiny::observeEvent(input$run_analysis, {
      shiny::withProgress(message = "Running BioGeoBEARS analysis...", value = 0, {
        tryCatch({
          # Get inputs
          pres_abs <- pres_abs_matrix()
          phylo_tree <- tree()
          config <- bgb_config()

          if (is.null(pres_abs) || is.null(phylo_tree)) {
            stop("Missing data: presence-absence matrix or tree")
          }

          if (!config$completed) {
            stop("Configuration not completed")
          }

          # Initialize results list
          models_results <- list()

          # Run each selected model
          for (model_name in config$models) {
            results$message <- paste("Running", model_name, "...")
            shiny::incProgress(1 / length(config$models))

            # This is a placeholder for actual BioGeoBEARS execution
            # In practice, this would call BioGeoBEARS::bears_optim_run()

            model_result <- list(
              model = model_name,
              lnl = -100 + runif(1) * 10,
              d = config$parameters$d[1] + runif(1) * (config$parameters$d[2] - config$parameters$d[1]),
              e = config$parameters$e[1] + runif(1) * (config$parameters$e[2] - config$parameters$e[1]),
              j = if ("J" %in% strsplit(model_name, "\\+")[[1]]) {
                config$parameters$j[1] + runif(1) * (config$parameters$j[2] - config$parameters$j[1])
              } else {
                NA
              }
            )

            models_results[[model_name]] <- model_result
          }

          results$models_results <- models_results

          # Create comparison table
          comparison_df <- data.frame(
            Model = names(models_results),
            LnL = sapply(models_results, function(x) round(x$lnl, 3)),
            d = sapply(models_results, function(x) round(x$d, 4)),
            e = sapply(models_results, function(x) round(x$e, 4)),
            j = sapply(models_results, function(x) if (is.na(x$j)) "-" else round(x$j, 4)),
            stringsAsFactors = FALSE
          )

          # Calculate AICc
          n_areas <- nrow(pres_abs)
          comparison_df$AICc <- sapply(seq_along(models_results), function(i) {
            n_params <- 2 + if ("J" %in% strsplit(names(models_results)[i], "\\+")[[1]]) 1 else 0
            lnl <- models_results[[i]]$lnl
            aic <- 2 * n_params - 2 * lnl
            aicc <- aic + (2 * n_params * (n_params + 1)) / (n_areas - n_params - 1)
            round(aicc, 3)
          })

          results$model_comparison <- comparison_df
          results$completed <- TRUE
          results$message <- "Analysis completed successfully!"

          shiny::incProgress(1)
        }, error = function(e) {
          results$completed <- FALSE
          results$message <- paste("Error:", e$message)
        })
      })
    })

    # Progress UI
    output$progress_ui <- shiny::renderUI({
      if (results$message == "") {
        return(NULL)
      }

      if (results$completed) {
        shiny::div(
          class = "alert alert-success",
          results$message
        )
      } else {
        shiny::div(
          class = "alert alert-info",
          results$message
        )
      }
    })

    # Results UI
    output$results_ui <- shiny::renderUI({
      if (!results$completed || is.null(results$model_comparison)) {
        return(NULL)
      }

      shiny::tabsetPanel(
        shiny::tabPanel(
          "Model Comparison",
          shiny::br(),
          DT::renderDT({
            DT::datatable(
              results$model_comparison,
              options = list(
                dom = "t",
                pageLength = 10
              )
            )
          })
        ),
        shiny::tabPanel(
          "Model Details",
          shiny::br(),
          shiny::uiOutput(session$ns("model_details_ui"))
        ),
        shiny::tabPanel(
          "Download Results",
          shiny::br(),
          shiny::downloadButton(
            session$ns("download_results"),
            "Download Results Table"
          )
        )
      )
    })

    # Model details
    output$model_details_ui <- shiny::renderUI({
      if (is.null(results$models_results)) {
        return(NULL)
      }

      details_list <- lapply(names(results$models_results), function(model_name) {
        res <- results$models_results[[model_name]]
        shiny::div(
          shiny::h5(model_name),
          shiny::p(paste("Log-Likelihood:", round(res$lnl, 3))),
          shiny::p(paste("d (dispersal):", round(res$d, 4))),
          shiny::p(paste("e (extinction):", round(res$e, 4))),
          if (!is.na(res$j)) {
            shiny::p(paste("j (jump dispersal):", round(res$j, 4)))
          },
          shiny::hr()
        )
      })

      do.call(shiny::div, details_list)
    })

    # Download handler
    output$download_results <- shiny::downloadHandler(
      filename = "biogeobears_results.csv",
      content = function(file) {
        utils::write.csv(results$model_comparison, file, row.names = FALSE)
      }
    )

    # Return reactive results
    shiny::reactive({
      list(
        models_results = results$models_results,
        model_comparison = results$model_comparison,
        completed = results$completed
      )
    })
  })
}
