#' Run the BioGeoBEARS Shiny Application
#'
#' Launches the BioGeoBEARS Shiny application for ancestral range estimation
#' and biogeographic analysis.
#'
#' @param ... Additional arguments passed to [shiny::shinyApp()]
#'
#' @return Invisibly returns the Shiny application object
#'
#' @examples
#' \dontrun{
#' run_biogeoshiny()
#' }
#'
#' @export
run_biogeoshiny <- function(...) {
  app_dir <- system.file("shiny", package = "biogeoshiny")
  if (app_dir == "") {
    stop("Could not find shiny app directory. Try reinstalling biogeoshiny.",
         call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal", ...)
}
