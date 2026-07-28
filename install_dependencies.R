# ==============================================================================
# Install all dependencies for BioRangeAnalyzer
# ==============================================================================
#
# Run this script BEFORE using the app for the first time, or AFTER updating R.
#
# Some packages are only available from Bioconductor or GitHub (not CRAN).
# This script handles all sources automatically.
#
# Usage:
#   source("install_dependencies.R")
#
# ==============================================================================

cat("=================================================================\n")
cat("  BioRangeAnalyzer - Installing Dependencies\n")
cat("=================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: CRAN packages
# ------------------------------------------------------------------------------

cat("STEP 1: Installing CRAN packages...\n")

cran_packages <- c(
  # Shiny and UI
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinythemes",
  
  # Interactive graphics
  "ggiraph",
  
  # Data manipulation
  "dplyr", "tidyr", "data.table",
  
  # Visualization
  "ggplot2", "plotly", "leaflet", "RColorBrewer", "viridis",
  

  # Data tables
  "DT",
  
  # Phylogenetics (CRAN)
  "ape", "phangorn", "phytools", "TreeSearch",
  
  # Spatial analysis
  "sp", "sf", "raster", "terra", "geosphere",
  
  # Optimization
  "GenSA", "FD",
  
  # Parallel computing
  "snow", "parallel",
  
  # Utilities
  "devtools", "remotes", "markdown", "knitr"
)

missing_cran <- cran_packages[!(cran_packages %in% rownames(installed.packages()))]

if (length(missing_cran) > 0) {
  cat("  Installing", length(missing_cran), "missing CRAN packages...\n")
  install.packages(missing_cran, dependencies = TRUE)
  cat("  Done.\n\n")
} else {
  cat("  All CRAN packages already installed.\n\n")
}

# ------------------------------------------------------------------------------
# STEP 2: Bioconductor packages (ggtree, rexpokit, cladoRcpp)
# ------------------------------------------------------------------------------

cat("STEP 2: Installing Bioconductor packages...\n")
cat("  (ggtree, rexpokit, cladoRcpp)\n")
cat("  NOTE: These are NOT on CRAN. They require BiocManager.\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  Installing BiocManager...\n")
  install.packages("BiocManager")
}

bioc_packages <- c("ggtree", "rexpokit", "cladoRcpp")
missing_bioc <- bioc_packages[!(bioc_packages %in% rownames(installed.packages()))]

if (length(missing_bioc) > 0) {
  cat("  Installing", length(missing_bioc), "Bioconductor packages...\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  cat("  Done.\n\n")
} else {
  cat("  All Bioconductor packages already installed.\n\n")
}

# ------------------------------------------------------------------------------
# STEP 3: GitHub packages (BioGeoBEARS)
# ------------------------------------------------------------------------------

cat("STEP 3: Installing GitHub packages...\n")
cat("  (BioGeoBEARS)\n")
cat("  NOTE: BioGeoBEARS is NOT on CRAN or Bioconductor.\n")
cat("  It must be installed from GitHub: nmatzke/BioGeoBEARS\n\n")

if (!requireNamespace("BioGeoBEARS", quietly = TRUE)) {
  cat("  Installing BioGeoBEARS from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("nmatzke/BioGeoBEARS", upgrade = "never")
  cat("  Done.\n\n")
} else {
  cat("  BioGeoBEARS already installed.\n\n")
}

# ------------------------------------------------------------------------------
# STEP 4: Verify installation
# ------------------------------------------------------------------------------

cat("=================================================================\n")
cat("  VERIFICATION\n")
cat("=================================================================\n\n")

all_packages <- c(cran_packages, bioc_packages, "BioGeoBEARS")
installed <- all_packages %in% rownames(installed.packages())

if (all(installed)) {
  cat("  ALL DEPENDENCIES INSTALLED SUCCESSFULLY!\n\n")
  cat("  Next steps:\n")
  cat("    1. devtools::document()           # Generate documentation\n")
  cat("    2. devtools::install_local('.')    # Install the package\n")
  cat("    3. BioRangeAnalyzer::run_biogeoshiny()  # Run the app\n")
} else {
  missing <- all_packages[!installed]
  cat("  WARNING: Some packages failed to install:\n")
  for (pkg in missing) {
    cat(paste0("    - ", pkg, "\n"))
  }
  cat("\n  Troubleshooting:\n")
  cat("    - Check your internet connection\n")
  cat("    - For Bioconductor packages: BiocManager::install('package_name')\n")
  cat("    - For BioGeoBEARS: remotes::install_github('nmatzke/BioGeoBEARS')\n")
  cat("    - For ggtree: BiocManager::install('ggtree')\n")
  cat("    - On macOS, you may need Xcode: xcode-select --install\n")
  cat("    - On Linux, you may need: sudo apt install libgdal-dev libproj-dev\n")
}

cat("\n=================================================================\n")
cat("  IMPORTANT: After updating R to a new major version,\n")
cat("  you MUST re-run this script to reinstall all packages.\n")
cat("=================================================================\n")
