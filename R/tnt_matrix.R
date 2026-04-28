#' tnt_matrix
#'
#' This function converts a presence-absence matrix into TNT (Tree analysis using
#' New Technology) format for parsimony analysis in Goloboff's TNT software. The
#' function works with matrices generated from different PAE-PCE methods including
#' minimum spanning trees (MST), buffers, and convex hulls.
#'
#' @param pres_abs Either a matrix object or a character string specifying the path
#'   to a text file containing the presence-absence matrix. The matrix should have
#'   taxa (grid cells or operational geographic units) as rows and species as columns,
#'   with binary values (0 = absence, 1 = presence). This can be the output from
#'   \code{pae_pce()} function using any method (raster, MST, buffer, or convex hull).
#'
#' @return The function does not return a value. It creates a TNT-formatted file
#'   named \code{pres_abs.tnt} in the \code{out/} directory. This file is ready
#'   to be loaded and analyzed in TNT software.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Reads the presence-absence matrix (from file or R object)
#'   \item Converts the matrix to TNT format with proper syntax
#'   \item Adds character names (species names) to the matrix
#'   \item Generates TNT commands for parsimony analysis
#'   \item Exports the formatted matrix to \code{out/pres_abs.tnt}
#' }
#'
#' The TNT file structure includes:
#' \itemize{
#'   \item \code{nstates num 32;} - Sets maximum number of character states
#'   \item \code{xread} - Command to read data matrix
#'   \item Matrix dimensions (number of characters and taxa)
#'   \item Binary presence-absence data for each taxon
#'   \item Character names (species names) with indices
#'   \item \code{proc/;} - Command to process the data
#' }
#'
#' @section Input Matrix Format:
#' The presence-absence matrix should be structured as:
#' \itemize{
#'   \item **Rows**: Taxa (grid cells, operational geographic units, or areas)
#'   \item **Columns**: Species (characters in parsimony analysis)
#'   \item **Values**: Binary (0 = absence, 1 = presence)
#'   \item **Row names**: Taxon identifiers (e.g., grid cell numbers)
#'   \item **Column names**: Species names
#' }
#'
#' @section Compatible Methods:
#' This function works with presence-absence matrices generated from:
#' \itemize{
#'   \item **Minimum Spanning Trees (MST)**: Network-based species ranges
#'   \item **Buffer method**: Buffer zones around occurrence points
#'   \item **Convex hulls**: Polygon-based species ranges
#' }
#'
#' @section TNT Analysis Workflow:
#' After generating the TNT file:
#' \enumerate{
#'   \item Open TNT software
#'   \item Load the file: \code{proc out/pres_abs.tnt;}
#'   \item Run parsimony analysis: \code{ienum;} or \code{mult;}
#'   \item View consensus tree: \code{nelsen;}
#'   \item Identify generalized tracks from tree topology
#'   \item Map results back to geographic space using \code{paeGridNumbers()}
#' }
#'
#' @note
#' \itemize{
#'   \item The output directory \code{out/} must exist before running the function
#'   \item The function will overwrite existing \code{pres_abs.tnt} file
#'   \item TNT software is required to analyze the generated file
#'   \item Row names and column names in the matrix are preserved in the TNT file
#'   \item The function uses \code{nstates num 32;} which allows up to 32 character states
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: From minimum spanning trees (MST)
#'
#' # todos os taxons:
#' mst_all_taxa <- terminal_node(coordin = resul1a$data_df, shape_file = neo,
#'  sobrepo = FALSE, caption = TRUE, resol = c(10, 10), seephylog = FALSE)
#'
#' # Generate TNT matrix
#' tnt_matrix(mst_all_taxa) # Output: out/pres_abs.tnt
#'
#' # Example 2: From convex hull analysis
#' rang1 <- calcRange_convexHull(x = lycipta_final,
#' shape_file = neo, resol = 5) # quadriculas de 5 graus
#'
#' tnt_matrix(rang1) # Output: out/pres_abs.tnt
#'
#' # Example 3: From buffer analysis
#' rang2 <- calcRange_buffers(xy = data, buffer.width = 10000, shape_file = neo,
#'  resol = 5)
#'
#' tnt_matrix(rang2) # Output: out/pres_abs.tnt
#'
#' # Example 4: From a matrix file
#' tnt_matrix("path/to/presence_absence_matrix.txt")
#'
#' # Example 5: From a matrix object in R
#' # Create example matrix
#' pa_matrix <- matrix(
#'   c(1, 0, 1, 0,
#'     1, 1, 0, 0,
#'     0, 1, 1, 1),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(
#'     c("Grid_1", "Grid_2", "Grid_3"),
#'     c("Species_A", "Species_B", "Species_C", "Species_D")
#'   )
#' )
#'
#' tnt_matrix(pa_matrix)
#'
#' # Complete workflow example:
#'
#' # Step 1: Extrapolate the distribution and generate TNT matrix
#'
#' # Step 2: Analyze in TNT software
#' # - Open TNT
#' # - Load file: proc out/pres_abs.tnt;
#' # - Run analysis: ienum; or mult;
#' # - Get consensus: nelsen;
#' # - Identify generalized tracks (clades in the tree)
#'
#' # Step 3: Visualize results back in R
#' # Assuming TNT identified grids 48 and 47 as a generalized track
#' paeGridNumbers(
#'   pae_pce_result = pae_result,
#'   shape_file = south_america,
#'   resol = c(10, 10),
#'   gridCell = c(48, 47),
#'   transp = 0.8
#' )
#' }
#'
#' @section TNT File Format:
#' The generated TNT file has the following structure:
#' \preformatted{
#' nstates num 32;
#' xread
#' 4 3
#' 1010 Grid_1
#' 1100 Grid_2
#' 0111 Grid_3
#' ;
#'
#' cnames
#' {0 Species_A;}
#' {1 Species_B;}
#' {2 Species_C;}
#' {3 Species_D;}
#' ;
#'
#' proc/;
#' }
#'
#' @seealso
#' \code{\link{pae_pce}} for running PAE-PCE analysis with different methods
#' \code{\link{paeGridNumbers}} for visualizing TNT results on maps
#'
#' @references
#' Goloboff, P. A., Farris, J. S., & Nixon, K. C. (2008). TNT, a free program
#' for phylogenetic analysis. Cladistics, 24(5), 774-786.
#'
#' Morrone, J. J. (1994). On the identification of areas of endemism.
#' Systematic Biology, 43(3), 438-441.
#'
#' Craw, R. C., Grehan, J. R., & Heads, M. J. (1999). Panbiogeography:
#' Tracking the history of life. Oxford University Press.
#'
#' Escalante, T. (2009). Un ensayo sobre regionalizacion biogeografica.
#' Revista Mexicana de Biodiversidad, 80(3), 551-560.
#'
#' @export

tnt_matrix <- function(pres_abs){
  # pres_abs = .txt file or an object which is the output from the node_terminal function

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

  ######################
  ## Temporarily remove ROOT #######
  ######################

  # Remove ROOT (last row) if present
  matTemp <- mat[-dim(mat)[1], ]

  ######################
  ## Filter empty ranges
  ######################

  # Calculate row sums (total occurrences per range)
  row_sums <- rowSums(matTemp)

  # Identify non-empty rows (at least one species present)
  non_empty_rows <- row_sums > 0

  # Filter matrix to keep only non-empty ranges
  matTemp1 <- mat[non_empty_rows, , drop = FALSE]
  matTemp2 <- rbind(matTemp1, rep(0, ncol(matTemp1)))
  rownames(matTemp2) <- c(rownames(matTemp1), 'ROOT')

  # ROOT returns...
  mat <- matTemp2

  # extracting the raster cell numbers
  taxa <- row.names(mat)

  # producing the *.txt files
  n.cara <- ncol(mat)
  n.taxa <- nrow(mat)

  nomes <- NULL

  conta <- 0

  # capture R output:
  zz <- textConnection('texto', "w")
  sink(zz)
  cat('nstates num 32;',sep = '\n')
  cat('xread',sep = '\n')
  #cat('/* Matrix for a PAE-PCE analysis */', sep = '\n')
  cat(n.cara, n.taxa, fill = T)

  # placing the character matrix pres_abs...
  for(j in 1:n.taxa){
    nomes[j] <- paste(c(taxa[j]),collapse = ', ')
    cat(paste(mat[j,], collapse = ''), fill = T,
        labels = paste0(nomes[j]))
  }

  cat(';',sep = '\n')

  # taxon names:
  cat(sep = '\n')
  cat('cnames', sep = '\n')
  for(i in 1:n.cara){
    cat(paste0(c('{',conta, ' ', colnames(mat)[i], ';'), collapse = ''), sep = '\n')
    conta <- conta + 1
  }

  cat(';',sep = '\n')
  cat(sep = '\n')
  # cat('cc+ 0.6;',sep = '\n')
  cat('proc/;',sep = '\n')
  sink()
  close(zz)

  # opening, writing and closing the file...
  if (!file.exists('out/'))
    dir.create('out/') # output directory for results

  exte <- 'out/pres_abs.tnt'
  tempor = file(exte, "wt")
  writeLines(texto, con = tempor)
  close(tempor)
}
