# Global configuration and shared resources for BRA (BioRangeAnalyzer Shiny)

library(raster)  # Load raster before terra to avoid conflicts
library(terra)   # Load terra after raster to avoid method conflicts
# Check dependencies first
source("check_dependencies.R", local = TRUE)

# Check if all dependencies are installed
dep_check <- check_and_install_dependencies()

if (!dep_check$all_installed) {
  cat("\n=== MISSING DEPENDENCIES ===")
  cat("\nThe following packages are missing:")
  cat(paste("  -", dep_check$missing_packages, collapse = "\n"), "\n\n")
  cat("Installing missing packages...\n")
  
  if (!install_missing_dependencies(dep_check$missing_packages)) {
    cat("\n✗ Failed to install some packages.")
    cat("\nPlease install them manually using:")
    cat("\ninstall.packages(c(", paste(paste0('"', dep_check$missing_packages, '"'), collapse = ", "), "))\n")
  }
}

# Source all custom functions
tryCatch({
  source("all_functions.R", local = FALSE)
  source("biogeobears_range_methods.R", local = FALSE)
}, error = function(e) {
  warning("Could not source custom functions: ", e$message)
})

bind_optional_fun <- function(path, fun_name, fallback) {
  stopifnot(is.function(fallback))
  env <- new.env(parent = globalenv())

  if (!file.exists(path)) {
    return(fallback)
  }

  sys.source(path, envir = env, chdir = TRUE)
  fun <- get0(fun_name, envir = env, inherits = FALSE)
  if (is.function(fun)) fun else fallback
}

tryCatch({
  pae_path_candidates <- c(
    file.path("R", "pae_pce.R"),
    file.path("..", "R", "pae_pce.R"),
    file.path("..", "..", "R", "pae_pce.R"),
    file.path(getwd(), "R", "pae_pce.R"),
    file.path(getwd(), "biogeoshiny", "R", "pae_pce.R"),
    file.path(getwd(), "biogeoshiny_improved_0.1.0", "biogeoshiny", "R", "pae_pce.R")
  )

  resolved_pae_path <- NULL
  for (cand in unique(pae_path_candidates)) {
    if (file.exists(cand)) {
      resolved_pae_path <- cand
      break
    }
  }

  if (!is.null(resolved_pae_path)) {
    fallback_fun <- get0("pae_pce", mode = "function", inherits = TRUE)
    if (!is.function(fallback_fun)) {
      fallback_fun <- function(...) {
        stop("pae_pce() is unavailable in runtime.")
      }
    }

    assign(
      "pae_pce",
      bind_optional_fun(
        path = resolved_pae_path,
        fun_name = "pae_pce",
        fallback = fallback_fun
      ),
      envir = .GlobalEnv
    )
    message("PAE-PCE function bound from: ", normalizePath(resolved_pae_path, winslash = "/", mustWork = FALSE))
  }
}, error = function(e) {
  warning("Could not pre-bind pae_pce(): ", e$message)
})

# Load required packages with error handling
required_packages <- c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinythemes",
  "leaflet", "DT", "ape", "phangorn", "TreeSearch", "dplyr", "tidyr", "ggplot2", "plotly", 
  "RColorBrewer", "stats", "sp", "raster", "terra", "sf", "viridis", "geosphere", "fossil"
)

for (pkg in required_packages) {
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
  }, error = function(e) {
    stop(paste("CRITICAL: Could not load package:", pkg, "-", e$message))
  })
}

# Optional packages (loaded if available)
optional_packages <- c("BioGeoBEARS", "phytools", "vegan", "spdep", "rgeos", "rgdal")

for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      library(pkg, character.only = TRUE, quietly = TRUE)
    }, error = function(e) {
      # Silent fail for optional packages
    })
  }
}

if (!requireNamespace("BioGeoBEARS", quietly = TRUE)) {
  warning(
    paste(
      "Optional package 'BioGeoBEARS' is not installed.",
      "BioGeoBEARS analysis/ancestral-range plotting tabs will not run until installed.",
      "Install using:",
      "install.packages('devtools', repos='https://cloud.r-project.org')",
      "devtools::install_github('nmatzke/BioGeoBEARS', INSTALL_opts='--byte-compile', upgrade='never')",
      sep = "\n"
    )
  )
}

# Set options
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB max file size
options(stringsAsFactors = FALSE)

# Define color palette
app_colors <- list(
  primary = "#3c8dbc",
  secondary = "#367fa9",
  success = "#00a65a",
  danger = "#dd4b39",
  warning = "#f39c12",
  info = "#00c0ef"
)

# Define model descriptions
model_descriptions <- list(
  DEC = "Dispersal-Extinction-Cladogenesis (2 parameters: d, e)",
  "DEC+J" = "DEC with jump dispersal (3 parameters: d, e, j)",
  DIVALIKE = "Dispersal-Vicariance model (2 parameters: d, e)",
  "DIVALIKE+J" = "DIVALIKE with jump dispersal (3 parameters: d, e, j)",
  BAYAREALIKE = "BayArea model (2 parameters: d, e)",
  "BAYAREALIKE+J" = "BAYAREALIKE with jump dispersal (3 parameters: d, e, j)"
)

# Define parameter descriptions
parameter_descriptions <- list(
  d = "Dispersal rate: rate of range expansion/colonization",
  e = "Extinction rate: rate of range contraction/local extinction",
  j = "Jump dispersal: probability of founder-event speciation",
  x = "Distance exponent: parameter for distance-dependent dispersal"
)

# Function to format numbers
format_number <- function(x, digits = 3) {
  if (is.na(x)) return("NA")
  format(round(x, digits), nsmall = digits)
}

# Function to validate CSV
validate_csv <- function(file_path) {
  tryCatch({
    data <- utils::read.csv(file_path)
    if (nrow(data) == 0) {
      return(list(valid = FALSE, message = "CSV file is empty"))
    }
    list(valid = TRUE, data = data)
  }, error = function(e) {
    list(valid = FALSE, message = paste("Error reading CSV:", e$message))
  })
}

# Function to validate tree file
validate_tree <- function(file_path) {
  tryCatch({
    # Try reading as Newick first
    tree <- tryCatch({
      ape::read.tree(file_path)
    }, error = function(e) {
      # Try reading as Nexus
      ape::read.nexus(file_path)
    })
    list(valid = TRUE, tree = tree, n_taxa = length(tree$tip.label))
  }, error = function(e) {
    list(valid = FALSE, message = paste("Error reading tree:", e$message))
  })
}

# Function to create presence-absence matrix
create_pres_abs_matrix <- function(occurrence_data, method = "buffer", ...) {
  species <- unique(occurrence_data$spp)
  areas <- 1:10
  
  matrix(
    sample(0:1, length(species) * length(areas), replace = TRUE),
    nrow = length(areas),
    ncol = length(species),
    dimnames = list(paste0("Area_", areas), species)
  )
}

# Function to calculate model statistics
calculate_model_stats <- function(lnl, n_params, n_areas) {
  aic <- 2 * n_params - 2 * lnl
  aicc <- aic + (2 * n_params * (n_params + 1)) / (n_areas - n_params - 1)
  
  list(
    aic = aic,
    aicc = aicc,
    lnl = lnl,
    n_params = n_params
  )
}

# Function to perform LRT
perform_lrt <- function(lnl_simple, lnl_complex, df) {
  chisq <- 2 * (lnl_complex - lnl_simple)
  pval <- stats::pchisq(chisq, df = df, lower.tail = FALSE)
  
  list(
    chisq = chisq,
    df = df,
    pval = pval,
    significant = pval < 0.05
  )
}


# Detect and set the project directory
# Try to find the RStudio project directory
project_dir <- NULL

# Method 1: Check if RStudio is active and get project directory
if (exists(".rs.getProjectDirectory")) {
  project_dir <- tryCatch({
    .rs.getProjectDirectory()
  }, error = function(e) NULL)
}

# Method 2: Look for .Rproj file in parent directories
if (is.null(project_dir)) {
  current_dir <- getwd()
  for (i in 1:10) {  # Check up to 10 levels up
    rproj_files <- list.files(current_dir, pattern = "\\.Rproj$", full.names = FALSE)
    if (length(rproj_files) > 0) {
      project_dir <- current_dir
      break
    }
    parent_dir <- dirname(current_dir)
    if (parent_dir == current_dir) break  # Reached root
    current_dir <- parent_dir
  }
}

# Method 3: Use current working directory if no project found
if (is.null(project_dir)) {
  project_dir <- getwd()
}

# Set working directory to project directory
setwd(project_dir)

# Create output directories for extrapolation results
output_dirs <- c(
  "out_buffers",
  "out_MCP",
  "out_MST",
  "out_grid",
  "out_irregular_bins",
  "out_TNT"
)

for (dir_name in output_dirs) {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
  }
}

cat("\n=== BRA (BioRangeAnalyzer Shiny) Started ===")
cat("\n✓ Output directories created in:", getwd())
cat("\n  Project directory detected:", project_dir)
cat("\n  Working directory:", getwd())
cat("\n  Files will be saved in subdirectories of:", getwd(), "\n\n")

# Function to load shapefiles from output directories
load_output_shapefiles <- function(extrap_method = NULL, since_time = NULL, taxa_filter = NULL) {
  shapefiles <- list()
  rasters <- list()
  shapefile_info <- list()
  if (is.null(taxa_filter)) {
    taxa_filter <- character(0)
  }
  taxa_filter <- unique(as.character(taxa_filter))
  taxa_filter <- taxa_filter[nzchar(taxa_filter)]

  file_is_current <- function(path, since_time_value) {
    if (is.null(since_time_value)) {
      return(TRUE)
    }
    mtime <- tryCatch(file.info(path)$mtime, error = function(e) NA)
    if (is.na(mtime)) {
      return(FALSE)
    }

    # Some filesystems record mtime with coarse precision (e.g., 1 second).
    # Use a small tolerance to avoid dropping outputs created in the same second.
    tol_secs <- 2
    mtime >= (since_time_value - tol_secs)
  }

  keep_by_taxa <- function(base_name, inferred_species, taxa_values) {
    if (length(taxa_values) == 0) {
      return(TRUE)
    }
    has_species_tag <- grepl("^(BUFF_|MCP_|mst_|MST_|presence_BUFF_|presence_MCP_|presence_mst_|presence_MST_)", base_name)
    if (!has_species_tag || is.na(inferred_species) || !nzchar(inferred_species)) {
      return(TRUE)
    }
    inferred_species %in% taxa_values
  }

  infer_species_id <- function(base_name) {
    patterns <- c("^BUFF_", "^MCP_", "^mst_", "^MST_", "^BINS_", "^presence_BUFF_", "^presence_MCP_", "^presence_mst_", "^presence_MST_")
    for (pat in patterns) {
      if (grepl(pat, base_name)) {
        return(sub(pat, "", base_name))
      }
    }
    NA_character_
  }
  
  # Define output directories and their properties
  output_dirs_info <- list(
    list(dir = "out_buffers", name = "Buffer", color = "#FF6B6B"),
    list(dir = "out_MCP", name = "Convex Hull", color = "#4ECDC4"),
    list(dir = "out_MST", name = "MST", color = "#45B7D1"),
    list(dir = "out_grid", name = "Regular Grid", color = "#7D8EA3"),
    list(dir = "out_irregular_bins", name = "Irregular Bins", color = "#96CEB4")
  )

  if (!is.null(extrap_method) && nzchar(extrap_method)) {
    method_dir_map <- c(
      buffer = "out_buffers",
      convex_hull = "out_MCP",
      mst = "out_MST",
      occurrence_only = "out_irregular_bins",
      irregular_bins = "out_irregular_bins"
    )
    if (extrap_method %in% names(method_dir_map)) {
      selected_dir <- unname(method_dir_map[extrap_method])
      selected_dirs <- unique(c(selected_dir, "out_grid"))
      output_dirs_info <- Filter(function(x) x$dir %in% selected_dirs, output_dirs_info)
    } else {
      output_dirs_info <- list()
    }
  }
  
  for (dir_info in output_dirs_info) {
    dir_path <- file.path(getwd(), dir_info$dir)
    
    if (dir.exists(dir_path)) {
      # Look for .shp files
      shp_files <- list.files(dir_path, pattern = "\\.shp$", full.names = TRUE)
      
      for (shp_file in shp_files) {
        tryCatch({
          # Skip pointshape files - they should not be displayed on the map
          base_name <- tools::file_path_sans_ext(basename(shp_file))
          if (grepl("^pointshape", base_name)) {
            next  # Skip this file
          }

          if (!file_is_current(shp_file, since_time)) {
            next
          }

          # Read shapefile using sf
          shp <- sf::st_read(shp_file, quiet = TRUE)
          species_id <- infer_species_id(base_name)
          if (!keep_by_taxa(base_name, species_id, taxa_filter)) {
            next
          }
          
          # Create a unique name
          unique_name <- paste0(dir_info$name, "_", base_name)
          
          # Store shapefile with metadata
          shapefiles[[unique_name]] <- shp
          shapefile_info[[unique_name]] <- list(
            name = unique_name,
            method = dir_info$name,
            species_id = species_id,
            color = dir_info$color,
            path = shp_file
          )
        }, error = function(e) {
          warning("Could not read shapefile ", shp_file, ": ", e$message)
        })
      }
      
      # Look for raster files (.tif, .asc, .grd)
      raster_files <- list.files(dir_path, pattern = "\\.(tif|asc|grd)$", full.names = TRUE)
      
      for (raster_file in raster_files) {
        tryCatch({
          base_name <- tools::file_path_sans_ext(basename(raster_file))

          if (!file_is_current(raster_file, since_time)) {
            next
          }

          species_id <- infer_species_id(base_name)
          if (!keep_by_taxa(base_name, species_id, taxa_filter)) {
            next
          }

          # Read raster using terra
          rst <- terra::rast(raster_file)
          
          # Create a unique name
          unique_name <- paste0(dir_info$name, "_", base_name)
          
          # Store raster
          rasters[[unique_name]] <- rst

          shapefile_info[[unique_name]] <- list(
            name = unique_name,
            method = dir_info$name,
            species_id = species_id,
            color = dir_info$color,
            path = raster_file,
            is_raster = TRUE
          )
        }, error = function(e) {
          warning("Could not read raster ", raster_file, ": ", e$message)
        })
      }
    }
  }
  
  list(
    shapefiles = shapefiles,
    rasters = rasters,
    shapefile_info = shapefile_info
  )
}
