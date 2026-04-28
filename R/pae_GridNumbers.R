#' pae_GridNumbers
#'
#' This function plots the results of a PAE-PCE analysis, displaying a shapefile
#' with grid cells overlaid. Specified grid cells (generalized tracks or endemic areas)
#' are highlighted and numbered. A georeferenced raster file is automatically generated
#' in the output directory. The function also displays a legend with the species that
#' characterize the highlighted areas (synapomorphic species).
#'
#' @param pae_pce_result The complete output list from the \code{pae_pce()} function,
#'   containing at minimum \code{$nonHomoplastic_species} and \code{$homoplastic_species} or
#'   a \emph{list} with the taxa and the number of grids from a TNT analysis with the following
#'   column titles: \emph{spp} and \emph{grid_n}
#' @param shape_file A spatial object (sf, SpatVector, or Spatial*) representing
#'   the study region used in the PAE-PCE analysis
#' @param resol Numeric vector of length 2 specifying the grid resolution
#'   in degrees (e.g., \code{c(10, 10)} for 10x10 degree cells). Must match the
#'   resolution used in the PAE-PCE analysis.
#' @param gridCell Numeric vector of grid cell numbers to be highlighted and numbered.
#'   These typically represent generalized tracks or endemic areas identified from
#'   TNT analysis results (Round-1, Round-2, etc.)
#' @param transp Numeric, transparency level for the highlighted grid cells
#'   (default: 0.8). Values range from 0 (fully transparent) to 1 (fully opaque).
#' @param xmin Numeric, minimum x-axis limit for plot extent (default: NULL, auto-calculated)
#' @param xmax Numeric, maximum x-axis limit for plot extent (default: NULL, auto-calculated)
#' @param ymin Numeric, minimum y-axis limit for plot extent (default: NULL, auto-calculated)
#' @param ymax Numeric, maximum y-axis limit for plot extent (default: NULL, auto-calculated)
#'
#' @return The function creates a plot and invisibly returns NULL. As a side effect,
#'   it generates a georeferenced raster file named \code{resultPAE_PCE_X_Y.tif} in
#'   the \code{out/} directory, where X_Y represents the highlighted grid cell numbers
#'   separated by underscores.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Creates a grid overlay on the shapefile with the specified resolution
#'   \item Extracts non-homoplastic (synapomorphic) species from PAE-PCE results
#'   \item Filters species that occur in the specified grid cells
#'   \item Highlights the specified grid cells with semi-transparent shading
#'   \item Labels the highlighted cells with their grid numbers
#'   \item Displays a legend showing which species characterize those areas
#'   \item Exports a georeferenced raster file (.tif) with the highlighted cells
#' }
#'
#' The grid numbering system follows a sequential pattern starting from the top-left
#' corner, proceeding left-to-right and top-to-bottom across the study region.
#'
#' The PAE-PCE result object must contain:
#' \itemize{
#'   \item \code{$nonHomoplastic_species}: List of data frames with columns 'spp' and 'grid_n'
#'   \item \code{$homoplastic_species}: List of data frames with columns 'spp' and 'grid_n'
#' }
#'
#' @note
#' \itemize{
#'   \item The output directory \code{out/} must exist before running the function
#'   \item The resolution (\code{resol}) must match the resolution used in the PAE-PCE analysis
#'   \item Grid cell numbers are determined by the raster cell indexing system
#'   \item The coordinate reference system (CRS) is automatically inherited from the shapefile
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Example 1: Visualize generalized tracks from PAE-PCE Round-1 (TNT results)
#' # Highlighting grid cells 48 and 47
#' pae_GridNumbers(
#'   pae_pce_result = my_pae_result,
#'   shape_file = asul,
#'   resol = c(10, 10),
#'   gridCell = c(48, 47),
#'   transp = 0.8
#' )
#'
#' # Example 2: Visualize generalized tracks from PAE-PCE Round-2 (TNT results)
#' # Highlighting grid cells 54, 48, and 47
#' pae_GridNumbers(
#'   pae_pce_result = my_pae_result,
#'   shape_file = asul,
#'   resol = c(10, 10),
#'   gridCell = c(54, 48, 47),
#'   transp = 0.8
#' )
#'
#' # Example 3: Custom plot extent with specific coordinates
#' pae_GridNumbers(
#'   pae_pce_result = my_pae_result,
#'   shape_file = asul,
#'   resol = c(10, 10),
#'   gridCell = c(48, 47),
#'   transp = 0.8,
#'   xmin = -80,
#'   xmax = -40,
#'   ymin = -30,
#'   ymax = 10
#' )
#'
#' # Example 4: Higher transparency for subtle highlighting
#' pae_GridNumbers(
#'   pae_pce_result = my_pae_result,
#'   shape_file = asul,
#'   resol = c(10, 10),
#'   gridCell = c(48, 47),
#'   transp = 0.5
#' )
#'
#' # Complete workflow example:
#' # Step 1: Run PAE-PCE analysis
#' result <- pae_pce(
#'   shape = asul,
#'   coords = lycipta.coords,
#'   resolution = c(10, 10)
#' )
#'
#' # Step 2: Generate presence/absence matrix for TNT
#' tnt_matrix(result)
#'
#' # Step 3: Run TNT analysis (external software)
#' # Identify generalized tracks from TNT output
#'
#' # Step 4: Visualize results in R
#' cells <- unique(as.numeric(pae1$nonHomoplastic_species[[1]]$grid_n)) # run 1
#' pae_GridNumbers(
#'   pae_pce_result = result,
#'   shape_file = asul,
#'   resol = c(10, 10),
#'   gridCell = cells,
#'   transp = 0.8
#' )
#'
#' # Example 5: when pae_pce-result is not a list and comes from TNT software:
#' library(rnaturalearth)
#' library(sf)
#' library(terra)
#' # Obter Americas (Norte + Sul)
#' americas <- ne_countries(continent = c("North America", "South America"),
#'                          returnclass = "sf")
#'
#' americas_vect <- vect(americas)
#' crs(americas_vect) <- "+proj=longlat +datum=WGS84"
#'
#' resulPae <- list()
#' # For example: species A, B and C in grids 156 and 157 # results from a TNT parsimony analysis
#' resulPae[[1]] <- data.frame(spp = c('A', 'A', 'B', 'B'), grid_n = rep(c(156, 157), 2))
#'
#' paeGridNumbers(shape_file = americas_vect,
#'  resol = c(10, 10), pae_pce_result = resulPae,
#'  gridCell = c(156, 157), transp = 0.8)
#' }
#'
#' @seealso
#' \code{\link{pae_pce}} for running PAE-PCE analysis
#' \code{\link{tnt_matrix}} for generating presence/absence matrices for TNT
#'
#' @references
#' Morrone, J. J. (1994). On the identification of areas of endemism.
#' Systematic Biology, 43(3), 438-441.
#'
#' Escalante, T. (2009). Un ensayo sobre regionalizacion biogeografica.
#' Revista Mexicana de Biodiversidad, 80(3), 551-560.
#'
#' @export

pae_GridNumbers <- function(pae_pce_result, shape_file, resol, gridCell, transp = 0.8, xmin = NULL,
                            xmax = NULL, ymin = NULL, ymax = NULL){
  library(raster)

  library(raster)
  ######################
  ##### shape file #####
  # shapeFile <- asul
  resolut <- resol

  # Verifica se e lista
  if (!is.list(pae_pce_result)) {
    stop("Error: 'pae_pce_result' must be a list...")
  }


  # loop:
  tabela_completa <- list()
  # Verifica componente nonHomoplastic_species
  if (!("nonHomoplastic_species" %in% names(pae_pce_result))){
    warning("'pae_pce_result' does not contain 'nonHomoplastic_species'...")
    for(li in 1:length(pae_pce_result)){
      resulTemp <- pae_pce_result[[li]]
      resulTemp <- resulTemp %>%
        filter(grid_n %in% gridCell)
      tabela_completa[[li]] <- resulTemp
    }
  } else {
    for(li in 1:length(pae_pce_result$nonHomoplastic_species)){
      resulTemp <- pae_pce_result$nonHomoplastic_species[[li]]
      resulTemp <- resulTemp %>%
        filter(grid_n %in% gridCell)
      tabela_completa[[li]] <- resulTemp
    }
  }
  tabela_total <- bind_rows(tabela_completa)

  speciesN <- unique(tabela_total$spp)

  cols1 <- setNames(n = viridis(length(speciesN)),
                    unique(speciesN))

  xmin = xmin; xmax = xmax; ymin = ymin; ymax = ymax

  grid <- raster(extent(as(shape_file, 'Spatial')),
                 resolution = resolut)
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  # proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  if(is.null(crs(shape_file))){
    crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  } else {
    crs(gridPolygon) <- crs(shape_file)
  }

  # clipping the intersected cells:
  cropped_map <- raster::intersect(gridPolygon, as(shape_file, 'Spatial'))
  # plot(cropped_map, xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = T)
  mask.raster <- raster(extent(as(shape_file, 'Spatial')), resolution = resolut,
                        crs = crs(shape_file))
  r <- rasterize(as(shape_file, 'Spatial'), mask.raster, fun = 'first')
  proj4string(r) <- crs(shape_file)
  r <- merge(r, mask.raster)
  ncellR <- ncell(r)
  r[r > 0] <- NA
  for(i in 1:length(gridCell)){
    for(j in 1:ncellR){
      if(j == gridCell[i]){
        values(r)[j] <- 1
      }
    }
  }

  # Set plot limits if not provided
  if (is.null(xmin)) xmin <- extent(cropped_map)[1]
  if (is.null(xmax)) xmax <- extent(cropped_map)[2]
  if (is.null(ymin)) ymin <- extent(cropped_map)[3]
  if (is.null(ymax)) ymax <- extent(cropped_map)[4]

  map.r <- raster::as.data.frame(raster::rasterToPoints(r))
  pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, mask = T)
  # plot(pontosRaster)
  gridNumber <- which(pontosRaster@data@values == 1)
  # map.r$gridNumber <- which(pontosRaster@data@values == 1)

  writeRaster(pontosRaster,
              paste0('out/resultPAE_PCE_', paste(gridNumber,
                                                 collapse = '_'), '.tif'), format = "GTiff",
              overwrite = TRUE)

  colores <- hcl.colors(1)
  brks <- seq(0, 1, by=1)

  plot(cropped_map, axes = TRUE, main = 'PAE-PCE final result with numbered grids',
       sub = paste0(c('resolution: ', resolut[1], ' x ', resolut[2]), collapse = ''), cex.main = 0.9,
       cex.sub = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))

  nspeciesN <- length(speciesN)
  legend(x = "bottomleft", legend = speciesN, pch = 16, col = names(cols1[1:length(speciesN)]),
         title = 'Species',  title.col = 'red', pt.cex = 0.6, cex = 0.6)

  plot(pontosRaster, add = T, axes = F, legend = F, col = colores, alpha = 0.40)

  map.r$gridNumber <- which(pontosRaster@data@values == 1)

  text(subset(map.r[,c(1, 2)], map.r$layer == 1), labels = gridNumber, cex = (0.8 * resolut[1]) / resolut[1],
       col = adjustcolor(col = 'black', alpha.f = 0.9), font = 2)
}
