#' Data Input Module UI
#'
#' @param id Module ID
#'
#' @return UI elements
#'
#' @keywords internal
mod_data_input_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Step 1: Upload Data"),
    shiny::hr(),

    # Guidance text
    shiny::div(
      class = "guidance-text",
      shiny::p(
        "Upload your occurrence data (CSV) and phylogenetic tree (Newick or Nexus format). ",
        "The occurrence data should have columns: species, longitude, latitude."
      )
    ),

    # Occurrence data upload
    shiny::div(
      class = "form-group",
      shiny::h5("Occurrence Data"),
      shiny::fileInput(
        ns("occurrence_file"),
        "Choose CSV/TXT file",
        accept = c(".csv", ".txt")
      ),
      shiny::uiOutput(ns("occurrence_preview"))
    ),

    shiny::hr(),

    # Tree file upload
    shiny::div(
      class = "form-group",
      shiny::h5("Phylogenetic Tree"),
      shiny::fileInput(
        ns("tree_file"),
        "Choose Newick/Nexus file",
        accept = c(".tre", ".nex", ".nexus", ".txt")
      ),
      shiny::uiOutput(ns("tree_preview"))
    ),

    shiny::hr(),

    # Optional: shapefile upload
    shiny::div(
      class = "form-group",
      shiny::h5("Optional: Study Area Boundary (Shapefile)"),
      shiny::p(
        "Upload a shapefile to define the study area boundary. ",
        "Include all files (.shp, .shx, .dbf, .prj)."
      ),
      shiny::fileInput(
        ns("shapefile"),
        "Choose shapefile files",
        multiple = TRUE,
        accept = c(".shp", ".shx", ".dbf", ".prj")
      )
    ),

    shiny::hr(),

    # Action button
    shiny::actionButton(
      ns("load_data"),
      "Load Data",
      class = "btn-primary btn-lg",
      width = "100%"
    ),

    # Status message
    shiny::uiOutput(ns("status_message"))
  )
}

#' Data Input Module Server
#'
#' @param id Module ID
#'
#' @return Reactive values with loaded data
#'
#' @keywords internal
mod_data_input_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    read_occurrence_file <- function(file_input) {
      file_ext <- tolower(tools::file_ext(file_input$name))

      if (file_ext == "txt") {
        tryCatch(
          utils::read.delim(file_input$datapath, stringsAsFactors = FALSE),
          error = function(e) {
            utils::read.table(
              file_input$datapath,
              header = TRUE,
              sep = ",",
              stringsAsFactors = FALSE
            )
          }
        )
      } else {
        utils::read.csv(file_input$datapath, stringsAsFactors = FALSE)
      }
    }

    # Reactive values
    data <- shiny::reactiveValues(
      occurrence = NULL,
      tree = NULL,
      shapefile = NULL,
      loaded = FALSE,
      message = ""
    )

    # Preview occurrence data
    output$occurrence_preview <- shiny::renderUI({
      if (is.null(input$occurrence_file)) {
        return(NULL)
      }

      tryCatch({
        occ_data <- read_occurrence_file(input$occurrence_file)
        preview_data <- utils::head(occ_data)

        shiny::div(
          shiny::p(
            "Preview (head):",
            class = "text-muted"
          ),
          DT::datatable(preview_data, options = list(dom = "t"), rownames = FALSE)
        )
      }, error = function(e) {
        shiny::div(
          class = "alert alert-danger",
          paste("Error reading file:", e$message)
        )
      })
    })

    # Preview tree
    output$tree_preview <- shiny::renderUI({
      if (is.null(input$tree_file)) {
        return(NULL)
      }

      tryCatch({
        tree <- ape::read.tree(input$tree_file$datapath)
        shiny::div(
          shiny::p(
            paste(
              "Tree loaded:",
              length(tree$tip.label),
              "taxa"
            ),
            class = "text-success"
          )
        )
      }, error = function(e) {
        shiny::div(
          class = "alert alert-danger",
          paste("Error reading tree:", e$message)
        )
      })
    })

    # Load data button
    shiny::observeEvent(input$load_data, {
      tryCatch({
        # Load occurrence data
        if (!is.null(input$occurrence_file)) {
          occ_data <- read_occurrence_file(input$occurrence_file)

          # Validate
          if (!all(c("species", "longitude", "latitude") %in% names(occ_data))) {
            stop("Occurrence data must have columns: species, longitude, latitude")
          }

          data$occurrence <- occ_data
        }

        # Load tree
        if (!is.null(input$tree_file)) {
          data$tree <- ape::read.tree(input$tree_file$datapath)
        }

        data$loaded <- TRUE
        data$message <- "Data loaded successfully!"
      }, error = function(e) {
        data$loaded <- FALSE
        data$message <- paste("Error:", e$message)
      })
    })

    # Status message
    output$status_message <- shiny::renderUI({
      if (data$message == "") {
        return(NULL)
      }

      if (data$loaded) {
        shiny::div(
          class = "alert alert-success",
          data$message
        )
      } else {
        shiny::div(
          class = "alert alert-danger",
          data$message
        )
      }
    })

    # Return reactive data
    shiny::reactive({
      list(
        occurrence = data$occurrence,
        tree = data$tree,
        loaded = data$loaded
      )
    })
  })
}
