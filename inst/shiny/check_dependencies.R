# Check and install missing dependencies for BioGeoBEARS Shiny App

check_and_install_dependencies <- function() {
  # List of required packages
  required_packages <- c(
    # Shiny and UI
    "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinythemes",
    
    # Data manipulation
    "dplyr", "tidyr",
    
    # Visualization
    "ggplot2", "plotly", "leaflet", "RColorBrewer", "viridis",
    
    # Data tables
    "DT",
    
    # Phylogenetics
    "ape", "phangorn", "TreeSearch",
    
    # Spatial analysis
    "sp", "sf", "raster", "terra", "geosphere",
    
    # Statistics
    "stats"
  )
  
  # Check which packages are missing
  installed <- rownames(installed.packages())
  missing <- required_packages[!(required_packages %in% installed)]
  
  if (length(missing) == 0) {
    return(list(
      all_installed = TRUE,
      missing_packages = NULL,
      message = "✓ All dependencies are installed!"
    ))
  }
  
  return(list(
    all_installed = FALSE,
    missing_packages = missing,
    message = paste("Missing packages:", paste(missing, collapse = ", "))
  ))
}

install_missing_dependencies <- function(packages) {
  if (length(packages) == 0) {
    return(TRUE)
  }
  
  tryCatch({
    cat("Installing missing packages...\n")
    install.packages(packages, dependencies = TRUE, quiet = TRUE)
    cat("✓ Installation complete!\n")
    return(TRUE)
  }, error = function(e) {
    cat("✗ Error installing packages:", e$message, "\n")
    return(FALSE)
  })
}
