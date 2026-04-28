#' range_nexus
#'
#' This function converts a presence-absence matrix into NEXUS format for
#' phylogenetic and biogeographic analysis. The function works with matrices
#' generated from different range extrapolation methods including buffers,
#' convex hulls (MPC), and minimum spanning trees (MST). Empty ranges (areas
#' with no species) are automatically removed.
#'
#' @param pres_abs Either:
#'   \itemize{
#'     \item A matrix object with presence-absence data
#'     \item A character string specifying path to a text file with the matrix
#'     \item A list containing \code{$pres_abs} component (output from calcRange functions)
#'   }
#'   The matrix should have ranges/areas as rows and species as columns, with
#'   binary values (0 = absence, 1 = presence).
#' @param meto Character string specifying the method used. Options are:
#'   \itemize{
#'     \item \code{"BUFF"} - Buffer method (from \code{calcRange_buffers()})
#'     \item \code{"MPC"} - Minimum Polygon Convex/Convex Hull (from \code{calcRange_convexHull()})
#'     \item \code{"MST"} - Minimum Spanning Tree (from \code{calcRange_mst()})
#'   }
#'   This parameter determines the output directory and file name.
#'
#' @return The function does not return a value. It creates a NEXUS-formatted
#'   file in the appropriate output directory:
#'   \itemize{
#'     \item \code{out_buffers/pres_abs_BUFF.nex} for buffer method
#'     \item \code{out_MPC/pres_abs_MPC.nex} for convex hull method
#'     \item \code{out/pres_abs_MST.nex} for MST method
#'   }
#'   Empty ranges (rows with all zeros) are automatically removed from the output.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Reads presence-absence matrix from file, matrix object, or list
#'   \item Removes ROOT taxon (last row) if present
#'   \item **Filters out empty ranges** (rows where no species occur)
#'   \item Replaces dots and spaces with underscores in species names
#'   \item Formats data according to NEXUS specifications
#'   \item Creates character state labels for each species
#'   \item Exports formatted NEXUS file ready for phylogenetic analysis
#' }
#'
#' @section NEXUS File Format:
#' The output file follows the standard NEXUS format:
#' \preformatted{
#' #NEXUS
#' begin data;
#' dimensions ntax=10 nchar=25;
#' format datatype=standard symbols="01" gap=-;
#'
#' CHARSTATELABELS
#' 1 Species_A,
#' 2 Species_B,
#' ...
#' 25 Species_Z;
#'
#' matrix
#' Species_A 1010101010
#' Species_B 1100110011
#' ...
#' ;
#' end;
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
#' @section Compatible Methods:
#' This function works with presence-absence matrices generated from:
#' \itemize{
#'   \item **Buffer method** (\code{calcRange_buffers()}): Creates buffer zones
#'     around occurrence points
#'   \item **Convex hull method** (\code{calcRange_convexHull()}): Creates minimum
#'     convex polygons around occurrence points
#'   \item **Minimum Spanning Tree** (\code{calcRange_mst()}): Connects occurrence
#'     points with minimum total distance
#' }
#'
#' @note
#' \itemize{
#'   \item The function automatically creates output directories if they don't exist
#'   \item Empty ranges (with no species occurrences) are automatically removed
#'   \item The ROOT taxon (last row) is automatically excluded from the output
#'   \item Species names with dots or spaces are converted to underscores
#'   \item The function will overwrite existing files with the same name
#'   \item Ensure the \code{meto} parameter matches the method used to generate ranges
#'   \item Fixed loop index error in CHARSTATELABELS section
#'   \item Fixed matrix reading from list (moved before matrix operations)
#' }
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # Example 1: From buffer analysis
#' # ============================================================================
#'
#' # Load functions
#' source('calcRange_buffers.R')
#' source('range_nexus.R')
#'
#' # Calculate ranges using buffer method
#' rang_buffer <- calcRange_buffers(
#'   xy = lycipta_final,
#'   shape_file = neo,
#'   resol = 10,
#'   buffer.width = 500000
#' )
#'
#' # Generate NEXUS file
#' range_nexus(pres_abs = rang_buffer, meto = 'BUFF')
#' # Output: out_buffers/pres_abs_BUFF.nex
#'
#' # ============================================================================
#' # Example 2: From convex hull analysis
#' # ============================================================================
#'
#' # Calculate ranges using convex hull
#' rang_hull <- calcRange_convexHull(
#'   xy = species_coords,
#'   shape_file = south_america,
#'   resol = 5
#' )
#'
#' # Generate NEXUS file
#' range_nexus(pres_abs = rang_hull, meto = 'MPC')
#' # Output: out_MPC/pres_abs_MPC.nex
#'
#' # ============================================================================
#' # Example 3: From MST analysis
#' # ============================================================================
#'
#' # Calculate ranges using MST
#' rang_mst <- calcRange_mst(
#'   xy = belostomatidae_coords,
#'   shape_file = neotropics,
#'   resol = 10
#' )
#'
#' # Generate NEXUS file
#' range_nexus(pres_abs = rang_mst, meto = 'MST')
#' # Output: out/pres_abs_MST.nex
#'
#' # ============================================================================
#' # Example 4: From matrix file
#' # ============================================================================
#'
#' # If you have a presence-absence matrix saved as text file
#' range_nexus(
#'   pres_abs = "path/to/presence_absence_matrix.txt",
#'   meto = 'BUFF'
#' )
#'
#' # ============================================================================
#' # Example 5: From matrix object
#' # ============================================================================
#'
#' # Create example matrix
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
#' # Generate NEXUS file
#' range_nexus(pres_abs = pa_matrix, meto = 'BUFF')
#'
#' # ============================================================================
#' # Example 6: Complete workflow for phylogenetic analysis
#' # ============================================================================
#'
#' # Step 1: Calculate ranges with buffer method
#' rang_buffer <- calcRange_buffers(
#'   xy = species_coords,
#'   shape_file = study_area,
#'   resol = 5,
#'   buffer.width = 300000
#' )
#'
#' # Step 2: Generate NEXUS file
#' range_nexus(pres_abs = rang_buffer, meto = 'BUFF')
#'
#' # Step 3: Load in phylogenetic software
#' # - Open PAUP*, MrBayes, BEAST, or other software
#' # - Load out_buffers/pres_abs_BUFF.nex
#' # - Run phylogenetic analysis
#' # - Infer biogeographic patterns from tree topology
#'
#' # ============================================================================
#' # Example 7: Comparing different methods
#' # ============================================================================
#'
#' # Method 1: Buffer
#' rang_buff <- calcRange_buffers(xy = coords, shape_file = shape,
#'                                 resol = 5, buffer.width = 200000)
#' range_nexus(pres_abs = rang_buff, meto = 'BUFF')
#'
#' # Method 2: Convex hull
#' rang_hull <- calcRange_convexHull(xy = coords, shape_file = shape, resol = 5)
#' range_nexus(pres_abs = rang_hull, meto = 'MPC')
#'
#' # Method 3: MST
#' rang_mst <- calcRange_mst(xy = coords, shape_file = shape, resol = 5)
#' range_nexus(pres_abs = rang_mst, meto = 'MST')
#'
#' # Compare resulting NEXUS files in phylogenetic analysis
#' }
#'
#' @seealso
#' \code{\link{calcRange_buffers}} for buffer-based range calculation
#' \code{\link{calcRange_convexHull}} for convex hull range calculation
#' \code{\link{calcRange_mst}} for MST-based range calculation
#' \code{\link{range_BioGeoBears}} for BioGeoBEARS format export
#' \code{\link{tnt_matrix}} for TNT format export
#'
#' @references
#' Maddison, D. R., Swofford, D. L., & Maddison, W. P. (1997). NEXUS: an
#' extensible file format for systematic information. Systematic Biology,
#' 46(4), 590-621.
#'
#' @export
range_nexus <- function(pres_abs, meto = NULL){

  ######################
  ## Input validation ##
  ######################

  # Check if pres_abs is provided
  if (missing(pres_abs)) {
    stop("Error: 'pres_abs' argument is required.
         Provide output from calcRange functions or a presence-absence matrix.")
  }

  # Check if NULL
  if (is.null(pres_abs)) {
    stop("Error: 'pres_abs' is NULL.
         Provide a valid presence-absence matrix or range calculation output.")
  }

  # Check if meto is provided
  if (is.null(meto)) {
    warning("Warning: 'meto' parameter not specified. Using default 'MST'.")
    meto <- "MST"
  }

  # Validate meto parameter
  valid_methods <- c("BUFF", "MPC", "MST")
  if (!(meto %in% valid_methods)) {
    stop(paste0("Error: 'meto' must be one of: ", paste(valid_methods, collapse = ", ")))
  }

  ######################
  ## Read input data ###
  ######################

  # Extract matrix from list BEFORE any matrix operations
  if (is.list(pres_abs)) {
    if ("pres_abs" %in% names(pres_abs)) {
      mat <- pres_abs$pres_abs
      message("Extracted presence-absence matrix from list")
    } else {
      stop("Error: List does not contain 'pres_abs' component.")
    }
  } else if (class(pres_abs)[1] == 'matrix') {
    mat <- pres_abs
  } else if (class(pres_abs)[1] == 'character') {
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
  mat <- mat[-nrow(mat), ]

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

  # Extract range identifiers
  taxa <- row.names(mat)

  # Get matrix dimensions (after filtering)
  n.cara <- ncol(mat)  # number of characters (species)
  n.taxa <- nrow(mat)  # number of taxa (ranges)

  # Replace dots with underscores in species names
  colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))

  ######################
  ## Generate NEXUS ####
  ######################

  nomes <- NULL
  conta <- 0

  # Capture R output to text connection
  zz <- textConnection('texto', "w")
  sink(zz)

  # NEXUS header
  cat('#NEXUS', sep = '\n')
  cat('begin data;', sep = '\n')
  cat(paste0('dimensions ntax=', n.taxa, ' nchar=', n.cara, ';'), sep = '\n')
  cat('format datatype=standard symbols="01" gap=-;', sep = '\n')
  cat('\n')

  # Character state labels (species names)
  cat('CHARSTATELABELS', sep = '\n')

  species <- colnames(mat)
  species <- gsub(pattern = " ", replacement = "_", x = species)

  for(k in 1:(n.cara - 1)){
    cat(sprintf("%s %s,", k, species[k]), sep = '\n')
  }
  # Last species without comma
  cat(sprintf("%s %s;", n.cara, species[n.cara]), sep = '\n')
  cat('\n')

  # Matrix section
  cat('matrix', sep = '\n')

  for(j in 1:n.cara){
    cat(paste(mat[, j], collapse = ''), fill = TRUE,
        labels = paste0(species[j]))
  }

  cat(';', sep = '\n')
  cat('end;')

  sink()
  close(zz)

  ######################
  ## Write output file #
  ######################

  # Determine output directory and file name based on method
  if(meto == 'MPC'){
    if (!file.exists('out_MPC/'))
      dir.create('out_MPC/')
    exte <- 'out_MPC/pres_abs_MPC.nex'
  } else if(meto == 'BUFF'){
    if (!file.exists('out_buffers/'))
      dir.create('out_buffers/')
    exte <- 'out_buffers/pres_abs_BUFF.nex'
  } else if(meto == 'MST'){
    if (!file.exists('out/'))
      dir.create('out/')
    exte <- 'out/pres_abs_MST.nex'
  }

  # Write file
  tempor <- file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)

  ######################
  ## Report results ####
  ######################

  message("[OK] NEXUS file created successfully")
  message(paste0("  Output file: ", exte))
  message(paste0("  Method: ", meto))
  message(paste0("  Number of species (characters): ", n.cara))
  message(paste0("  Number of ranges (taxa, after filtering): ", n.taxa))
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
  message("  1. Load NEXUS file in phylogenetic software (PAUP*, MrBayes, BEAST)")
  message(paste0("  2. File location: ", exte))
  message("  3. Run phylogenetic analysis")
  message("  4. Interpret biogeographic patterns from tree topology")

  # Return summary invisibly
  invisible(list(
    output_file = exte,
    n_species = n.cara,
    n_ranges = n.taxa,
    removed_ranges = removed_ranges
  ))
}
