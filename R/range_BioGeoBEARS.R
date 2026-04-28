#' range_BioGeoBEARS
#'
#' This function converts a presence-absence matrix into a BioGeoBEARS-formatted
#' geography file (.data format) for historical biogeography analysis. The function
#' automatically removes empty ranges (areas with no taxa) and formats the data
#' according to BioGeoBEARS requirements. It works with range data from different
#' methods including buffers, convex hulls (MPC), and minimum spanning trees (MST).
#'
#' @param pres_abs Either a matrix object, a list containing \code{$pres_abs}
#'   component, or a character string specifying the path to a text file containing
#'   the presence-absence matrix. This is typically the output from range calculation
#'   functions such as \code{calcRange_buffers()}, \code{calcRange_convexHull()},
#'   or \code{calcRange_mst()}.
#' @param meto Character string specifying the method used to calculate ranges.
#'   Options are:
#'   \itemize{
#'     \item \code{"BUFF"} - Buffer method (from \code{calcRange_buffers()})
#'     \item \code{"MPC"} - Minimum Polygon Convex/Convex Hull method (from \code{calcRange_convexHull()})
#'     \item \code{"MST"} - Minimum Spanning Tree method (from \code{calcRange_mst()})
#'   }
#'   This parameter determines the output directory and file name.
#'
#' @return The function does not return a value. It creates a BioGeoBEARS-formatted
#'   geography file in the appropriate output directory:
#'   \itemize{
#'     \item \code{out_buffers/pres_abs_BUFF_geog.data} for buffer method
#'     \item \code{out_MPC/pres_abs_MPC_geog.data} for convex hull method
#'     \item \code{out/pres_abs_MST_geog.data} for MST method
#'   }
#'   Empty ranges (rows with all zeros) are automatically removed from the output.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Reads the presence-absence matrix from file or R object
#'   \item Removes ROOT taxon (last row) if present
#'   \item **Filters out empty ranges** (rows where no species occur)
#'   \item Replaces dots and spaces with underscores in species names
#'   \item Formats data according to BioGeoBEARS specifications
#'   \item Creates output directory if it doesn't exist
#'   \item Exports formatted file ready for BioGeoBEARS analysis
#' }
#'
#' @section BioGeoBEARS File Format:
#' The output file follows the BioGeoBEARS geography file format:
#' \preformatted{
#' n_species n_ranges (range1 range2 range3 ...)
#' Species_A 101010
#' Species_B 110001
#' Species_C 001100
#' }
#'
#' Where:
#' \itemize{
#'   \item First line: number of species, number of ranges, and range identifiers
#'   \item Following lines: species name followed by binary presence-absence string
#'   \item 1 = species present in that range, 0 = species absent
#' }
#'
#' @section Input Matrix Structure:
#' The presence-absence matrix should be structured as:
#' \itemize{
#'   \item **Rows**: Geographic ranges (grid cells, areas, or operational geographic units)
#'   \item **Columns**: Species (taxa)
#'   \item **Values**: Binary (0 = absence, 1 = presence)
#'   \item **Row names**: Range identifiers (e.g., grid cell numbers)
#'   \item **Column names**: Species names
#'   \item **Last row**: ROOT (automatically removed by the function)
#' }
#'
#' @section Compatible Range Methods:
#' This function works with presence-absence matrices generated from:
#' \itemize{
#'   \item **Buffer method** (\code{calcRange_buffers()}): Creates buffer zones around
#'     occurrence points with specified width
#'   \item **Convex hull method** (\code{calcRange_convexHull()}): Creates minimum
#'     convex polygons around occurrence points
#'   \item **Minimum Spanning Tree** (\code{calcRange_mst()}): Connects occurrence
#'     points with minimum total distance
#' }
#'
#' @section BioGeoBEARS Analysis Workflow:
#' After generating the geography file:
#' \enumerate{
#'   \item Load the geography file in R
#'   \item Prepare phylogenetic tree in Newick format
#'   \item Set up BioGeoBEARS models (DEC, DIVALIKE, BAYAREALIKE, etc.)
#'   \item Run biogeographical analysis
#'   \item Visualize ancestral range reconstructions
#' }
#'
#' @note
#' \itemize{
#'   \item The function automatically creates output directories if they don't exist
#'   \item Empty ranges (with no species occurrences) are automatically removed
#'   \item The ROOT taxon (last row) is automatically excluded from the output
#'   \item Species names with dots or spaces are converted to underscores
#'   \item The function will overwrite existing files with the same name
#'   \item Ensure the \code{meto} parameter matches the method used to generate the ranges
#' }
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # Example 1: Buffer method
#' # ============================================================================
#'
#' # Load required functions
#' source('calcRange_buffers.R')
#' source('range_BioGeoBears.R')
#'
#' # Calculate ranges using buffer method
#' rang_buffer <- calcRange_buffers(
#'   xy = lycipta_final,
#'   shape_file = neo,
#'   resol = 10,
#'   buffer.width = 500000  # 500 km buffer
#' )
#'
#' # Generate BioGeoBEARS file
#' range_BioGeoBears(pres_abs = rang_buffer, meto = 'BUFF')
#' # Output: out_buffers/pres_abs_BUFF_geog.data
#'
#' # ============================================================================
#' # Example 2: Convex hull method (MPC)
#' # ============================================================================
#'
#' # Calculate ranges using convex hull method
#' rang_hull <- calcRange_convexHull(
#'   xy = species_coords,
#'   shape_file = south_america,
#'   resol = 5
#' )
#'
#' # Generate BioGeoBEARS file
#' range_BioGeoBears(pres_abs = rang_hull, meto = 'MPC')
#' # Output: out_MPC/pres_abs_MPC_geog.data
#'
#' # ============================================================================
#' # Example 3: Minimum Spanning Tree method
#' # ============================================================================
#'
#' # Calculate ranges using MST method
#' rang_mst <- calcRange_mst(
#'   xy = belostomatidae_coords,
#'   shape_file = neotropics,
#'   resol = 10
#' )
#'
#' # Generate BioGeoBEARS file
#' range_BioGeoBears(pres_abs = rang_mst, meto = 'MST')
#' # Output: out/pres_abs_MST_geog.data
#'
#' # ============================================================================
#' # Example 4: From matrix file
#' # ============================================================================
#'
#' # If you have a presence-absence matrix saved as text file
#' range_BioGeoBears(
#'   pres_abs = "path/to/presence_absence_matrix.txt",
#'   meto = 'BUFF'
#' )
#'
#' # ============================================================================
#' # Example 5: From matrix object in R
#' # ============================================================================
#'
#' # Create example presence-absence matrix
#' pa_matrix <- matrix(
#'   c(1, 0, 1, 0, 1,
#'     1, 1, 0, 0, 0,
#'     0, 1, 1, 1, 0,
#'     0, 0, 0, 1, 1),
#'   nrow = 4, byrow = TRUE,
#'   dimnames = list(
#'     c("Range_1", "Range_2", "Range_3", "Range_4"),
#'     c("Species_A", "Species_B", "Species_C", "Species_D", "Species_E")
#'   )
#' )
#'
#' # Add ROOT row (will be automatically removed)
#' pa_matrix <- rbind(pa_matrix, ROOT = rep(0, ncol(pa_matrix)))
#'
#' # Generate BioGeoBEARS file
#' range_BioGeoBears(pres_abs = pa_matrix, meto = 'BUFF')
#'
#' # ============================================================================
#' # Example 6: Complete workflow with BioGeoBEARS analysis
#' # ============================================================================
#'
#' library(BioGeoBEARS)
#'
#' # Step 1: Calculate ranges with buffer method
#' rang2 <- calcRange_buffers(
#'   xy = lycipta_final,
#'   shape_file = neo,
#'   resol = 10,
#'   buffer.width = 500000
#' )
#'
#' # Step 2: Generate BioGeoBEARS geography file
#' range_BioGeoBears(pres_abs = rang2, meto = 'BUFF')
#'
#' # Step 3: Load geography file in BioGeoBEARS
#' geogfn <- "out_buffers/pres_abs_BUFF_geog.data"
#' tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
#'
#' # Step 4: Load phylogenetic tree
#' trfn <- "path/to/phylogeny.newick"
#' tr <- read.tree(trfn)
#'
#' # Step 5: Set up and run BioGeoBEARS analysis
#' BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
#' BioGeoBEARS_run_object$trfn <- trfn
#' BioGeoBEARS_run_object$geogfn <- geogfn
#'
#' # Run DEC model
#' resDEC <- bears_optim_run(BioGeoBEARS_run_object)
#'
#' # Step 6: Visualize results
#' plot_BioGeoBEARS_results(
#'   results_object = resDEC,
#'   analysis_titletxt = "DEC Analysis",
#'   addl_params = list("j"),
#'   plotwhat = "text",
#'   label.offset = 0.45
#' )
#'
#' # ============================================================================
#' # Example 7: Different buffer widths comparison
#' # ============================================================================
#'
#' # Small buffer (100 km)
#' rang_small <- calcRange_buffers(
#'   xy = species_coords,
#'   shape_file = study_area,
#'   resol = 5,
#'   buffer.width = 100000
#' )
#' range_BioGeoBears(pres_abs = rang_small, meto = 'BUFF')
#'
#' # Large buffer (500 km)
#' rang_large <- calcRange_buffers(
#'   xy = species_coords,
#'   shape_file = study_area,
#'   resol = 5,
#'   buffer.width = 500000
#' )
#' range_BioGeoBears(pres_abs = rang_large, meto = 'BUFF')
#'
#' # Compare results in BioGeoBEARS to assess sensitivity
#'
#' # ============================================================================
#' # Example 8: Handling list output
#' # ============================================================================
#'
#' # If calcRange function returns a list with multiple components
#' range_list <- calcRange_buffers(
#'   xy = coords,
#'   shape_file = shape,
#'   resol = 10,
#'   buffer.width = 300000
#' )
#'
#' # Function automatically extracts $pres_abs component
#' range_BioGeoBears(pres_abs = range_list, meto = 'BUFF')
#'
#' # Or extract manually
#' range_BioGeoBears(pres_abs = range_list$pres_abs, meto = 'BUFF')
#' }
#'
#' @seealso
#' \code{\link{calcRange_buffers}} for calculating ranges using buffer method
#' \code{\link{calcRange_convexHull}} for calculating ranges using convex hull method
#' \code{\link{calcRange_mst}} for calculating ranges using minimum spanning tree method
#' \code{\link[BioGeoBEARS]{define_BioGeoBEARS_run}} for setting up BioGeoBEARS analysis
#'
#' @references
#' Matzke, N. J. (2013). Probabilistic historical biogeography: new models for
#' founder-event speciation, imperfect detection, and fossils allow improved
#' accuracy and model-testing. Frontiers of Biogeography, 5(4), 242-248.
#'
#' Matzke, N. J. (2014). Model selection in historical biogeography reveals that
#' founder-event speciation is a crucial process in island clades. Systematic
#' Biology, 63(6), 951-970.
#'
#' Ree, R. H., & Smith, S. A. (2008). Maximum likelihood inference of geographic
#' range evolution by dispersal, local extinction, and cladogenesis. Systematic
#' Biology, 57(1), 4-14.
#'
#' @export
range_BioGeoBears <- function(pres_abs, meto = NULL){

  ######################
  ## Input validation ##
  ######################

  # Check if pres_abs is provided
  if (missing(pres_abs)) {
    stop("Error: 'pres_abs' argument is required.
         Provide output from calcRange_buffers(), calcRange_convexHull(), or calcRange_mst().")
  }

  # Check if NULL
  if (is.null(pres_abs)) {
    stop("Error: 'pres_abs' is NULL.
         Provide a valid presence-absence matrix or range calculation output.")
  }

  # Check if meto is provided
  if (is.null(meto)) {
    warning("Warning: 'meto' parameter not specified. Using default output directory 'out/'.")
    meto <- "MST"  # Default to MST
  }

  # Validate meto parameter
  valid_methods <- c("BUFF", "MPC", "MST")
  if (!(meto %in% valid_methods)) {
    stop(paste0("Error: 'meto' must be one of: ", paste(valid_methods, collapse = ", ")))
  }

  ######################
  ## Read input data ###
  ######################

  # Extract matrix from list if needed
  if(is.list(pres_abs)){
    if("pres_abs" %in% names(pres_abs)){
      pres_abs <- pres_abs$pres_abs
    } else {
      stop("Error: List does not contain 'pres_abs' component.")
    }
  }

  # Read matrix from file or object
  if(class(pres_abs)[1] == 'matrix'){
    mat <- pres_abs
  } else if(class(pres_abs)[1] == 'character'){
    if (!file.exists(pres_abs)) {
      stop(paste0("Error: File not found: ", pres_abs))
    }
    mat <- read.table(file = pres_abs, sep = '', dec = '.', row.names = 1)
  } else {
    stop("Error: 'pres_abs' must be a matrix, list with $pres_abs, or file path.")
  }

  # Store original dimensions
  original_n_ranges <- nrow(mat)
  original_n_species <- ncol(mat)

  ######################
  ## Remove ROOT #######
  ######################

  # Remove ROOT (last row) if present
  mat <- mat[-dim(mat)[1], ]

  ######################
  ## Filter empty ranges
  ######################

  # Calculate row sums (total occurrences per range)
  row_sums <- rowSums(mat)

  # Identify non-empty rows (at least one species present)
  non_empty_rows <- row_sums > 0

  # Store removed ranges for reporting
  removed_ranges <- rownames(mat)[!non_empty_rows]

  # Filter matrix to keep only non-empty ranges
  mat <- mat[non_empty_rows, , drop = FALSE]

  ######################
  ## Extract info ######
  ######################

  # Extract range identifiers (grid cell numbers or areas)
  taxa <- row.names(mat)

  # Get matrix dimensions (after filtering)
  n.cara <- ncol(mat)  # number of species
  n.taxa <- nrow(mat)  # number of ranges (after removing empty ones)

  # Replace dots with underscores in species names
  colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))

  ######################
  ## Generate output ###
  ######################

  nomes <- NULL
  conta <- 0

  # Capture R output to text connection
  zz <- textConnection('texto', "w")
  sink(zz)

  # Header line: n_species n_ranges (range_list)
  cat(paste(n.cara, n.taxa, paste0('(', paste(taxa, collapse = ' '), ')'), sep = ' '), sep = '\n')

  # Format species names (replace spaces with underscores)
  species <- colnames(mat)
  species <- gsub(pattern = " ", replacement = "_", x = species)

  # Write presence-absence data for each species
  for(j in 1:n.cara){
    cat(paste(mat[, j], collapse = ''), fill = TRUE,
        labels = paste0(species[j]))
  }

  sink()
  close(zz)

  ######################
  ## Write output file #
  ######################

  # Determine output directory and file name based on method
  if(meto == 'MPC'){
    if (!file.exists('out_MPC/'))
      dir.create('out_MPC/')
    exte <- 'out_MPC/pres_abs_MPC_geog.data'
  } else if(meto == 'BUFF'){
    if (!file.exists('out_buffers/'))
      dir.create('out_buffers/')
    exte <- 'out_buffers/pres_abs_BUFF_geog.data'
  } else if(meto == 'MST'){
    if (!file.exists('out/'))
      dir.create('out/')
    exte <- 'out/pres_abs_MST_geog.data'
  }

  # Write file
  tempor <- file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)

  ######################
  ## Report results ####
  ######################

  message("[OK] BioGeoBEARS geography file created successfully")
  message(paste0("  Output file: ", exte))
  message(paste0("  Method: ", meto))
  message(paste0("  Number of species: ", n.cara))
  message(paste0("  Number of ranges (after filtering): ", n.taxa))
  message(paste0("  Original ranges (before filtering): ", original_n_ranges - 1, " (ROOT excluded)"))

  if(length(removed_ranges) > 0){
    message(paste0("  Empty ranges removed: ", length(removed_ranges)))
    if(length(removed_ranges) <= 20){
      message(paste0("    Removed: ", paste(removed_ranges, collapse = ", ")))
    } else {
      message(paste0("    Removed: ", paste(removed_ranges[1:20], collapse = ", "), " ... (",
                     length(removed_ranges) - 20, " more)"))
    }
  } else {
    message("  No empty ranges found (all ranges have at least one species)")
  }

  message("\nNext steps:")
  message("  1. Load geography file in BioGeoBEARS:")
  message(paste0("     tipranges <- getranges_from_LagrangePHYLIP('", exte, "')"))
  message("  2. Prepare phylogenetic tree in Newick format")
  message("  3. Set up and run BioGeoBEARS models (DEC, DIVALIKE, BAYAREALIKE)")
}
