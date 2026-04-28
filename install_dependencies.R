# Install all dependencies for BioGeoBEARS Shiny App
# Run this script before using the app for the first time

cat("Installing BioGeoBEARS Shiny App dependencies...\n\n")

# List of required packages
packages <- c(
  # Shiny and UI
  "shiny",
  "shinydashboard",
  "shinyjs",
  "shinyWidgets",
  "shinythemes",
  
  # Data manipulation
  "dplyr",
  "tidyr",
  "data.table",
  
  # Visualization
  "ggplot2",
  "plotly",
  "leaflet",
  "RColorBrewer",
  "viridis",
  
  # Data tables
  "DT",
  
  # Phylogenetics
  "ape",
  "phytools",
  
  # Spatial analysis
  "sp",
  "sf",
  "raster",
  "terra",
  "rgeos",
  "rgdal",
  "geosphere",
  "mapdata",
  
  # Statistics
  "vegan",
  "spdep",
  
  # Utilities
  "devtools"
)

# Install missing packages
missing_packages <- packages[!(packages %in% rownames(installed.packages()))]

if (length(missing_packages) > 0) {
  cat("Installing", length(missing_packages), "missing packages...\n")
  install.packages(missing_packages, dependencies = TRUE, quiet = TRUE)
  cat("✓ Installation complete!\n\n")
} else {
  cat("✓ All packages are already installed!\n\n")
}

# Verify installation
cat("Verifying installation...\n")
all_installed <- all(packages %in% rownames(installed.packages()))

if (all_installed) {
  cat("✓ All dependencies are installed and ready to use!\n")
  cat("You can now run: devtools::install_local('biogeoshiny')\n")
  cat("Then: biogeoshiny::run_biogeoshiny()\n")
} else {
  missing <- packages[!(packages %in% rownames(installed.packages()))]
  cat("✗ Some packages failed to install:\n")
  cat(paste("  -", missing, collapse = "\n"), "\n")
  cat("Please install them manually or check your internet connection.\n")
}
