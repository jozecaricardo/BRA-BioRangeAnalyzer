#' speciesGeneTrackPlot
#'
#' Produces a map - to be chosen by the user (in shapefile format) - and plots estimated generalized tracks
#' with PAE-PCE, in the form of grid cells with the same resolution selected in the analysis.
#'
#' @param PAE_PCE object produced by the function \code{\link[PanBioGeo]{pae_pce}} with the list of species (non-homoplastic
#' synapomorphies) of each clade found
#' @param resol a vector with resolution of quadrats
#' @param spvect a polygon of some place already with datum to be used
#' @param runPaePce a vector indicating which PAE-PCE run is to be indicated on the map
#' @return Map with hatched and properly numbered grid cells indicating the estimated generalized tracks and the
#' occurrence points of species (synapomorphies) associated with their MSTs located in the /out directory
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function plots the generalized tracks produced by the pae_pce function - one for each run (when applicable) -
#' for a given set of species, whose occurrence points and MSTs are already available in the /out directory.
#' The produced generalized track is also found in the /out directory with the name "GeneralizedTrack_#.tif", where #
#' represents the run number of the PAE-PCE analysis.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025. Panbiotracks: software for track analysis.
#'  \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'  Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1.
#' @seealso
#' \code{\link[phytools]{phytools}} for phylogenetic trees,
#' \code{\link[terra]{terra}} for spatial data handling,
#' \code{\link[PanBioGeo]{MST_node}} for the production of presence-absence matrices.
#' \code{\link[PanBioGeo]{pae_pce}} for the production of generalized tracks.
#' @export
#' @examples
#' # First, determine the minimum spanning trees of all taxa and produce a presence-absence matrix
#' # ---------------------------- with all taxa ------------------------------
#' # ----------- with a tree and resolution of 10 X 10 degrees ---------------
#' mst_all_taxa <- MST_node(coordin = resul$data_df, shape_file = lycipta.asul, sobrepo = F, caption = T, resol = c(10, 10), tree = resul$treeMod)
#'
#' # --------------------------------- generalized track ---------------------
#' # Now, produce the generalized track and a list with grid numbers and the characters (i.e., taxaset)
#' # Please try to run at least more than one round, with N > 1...
#' # With a new plot using the polygon of South America (asul)
#' pae1 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10), gridView = TRUE, labelGrid = TRUE, nonHomoplasticSpeciesList = FALSE, N = 10)
#' pae1 # data frame showing the homoplastic and non-homoplastic species (spp) and their grid numbers (grid_n), from the anterior iterations (PAE-PCE)
#'
#' speciesPlot(PAE_PCE = pae1, resol = c(10, 10), spvect = neo, runPaePce = 1)

speciesGeneTrackPlot <- function(PAE_PCE, resol, spvect, runPaePce){
  require(viridis)
  require(dplyr)
  require(terra)
  require(sf)
  require(sp)
  require(raster)

  resol <- resol
  vec_sp <- as(spvect, 'Spatial')
  datum <- "+proj=longlat +datum=WGS84"

  # Create base raster (EXACTLY as pae_pce)
  mask.raster <- raster::raster(raster::extent(vec_sp),
                                resolution = resol,
                                crs = raster::crs(datum))
  r <- raster::rasterize(vec_sp, mask.raster)
  raster::crs(r) <- raster::crs(datum)
  r <- raster::merge(r, mask.raster)

  # Create grid polygons (EXACTLY as pae_pce)
  grid_base <- raster::raster(raster::extent(vec_sp),
                              resolution = resol,
                              crs = raster::crs(datum))
  grid_base <- raster::rasterize(vec_sp, grid_base, field = 1)
  grid_base[is.na(grid_base)] <- 0
  grid_base[grid_base == 0] <- NA
  cropped_map <- raster::rasterToPolygons(grid_base, dissolve = FALSE)
  raster::crs(cropped_map) <- datum

  # Load the SAVED raster file (pontosRaster)
  geneTrack_file <- list.files('out/',
                               pattern = '\\.tif$',
                               full.names = TRUE)
  pontosRaster_loaded <- raster::raster(geneTrack_file[runPaePce])

  # Create grid numbering (EXACTLY as pae_pce line 707)
  map.r <- as.data.frame(raster::rasterToPoints(pontosRaster_loaded))
  map.r$gridNumber <- raster::cellFromXY(r, map.r[, c("x", "y")])

  # Extract synapomorphic grids
  synapomorphic_grids <- unique(PAE_PCE$nonHomoplastic_species[[runPaePce]]$grid_n)

  # PLOT
  plot(vec_sp, axes = TRUE, las = 1,
       main = "PAE-PCE: Generalized Track with Synapomorphic Grids")

  # Plot loaded raster (this is what pae_pce created)
  plot(pontosRaster_loaded, add = TRUE, col = gray.colors(100, start = 0.3, end = 0.7),
       alpha = 0.6, legend = FALSE)

  # Add grid
  plot(cropped_map, add = TRUE, border = "black", lwd = 0.8)

  # Add labels
  map.r_filtered <- map.r[map.r$gridNumber %in% synapomorphic_grids, ]
  text(map.r_filtered[, c("x", "y")],
       labels = map.r_filtered$gridNumber,
       cex = 0.8, col = 'red', font = 2)

  # pontos e MSTs:
  especies <- unique(PAE_PCE$nonHomoplastic_species[[runPaePce]]$spp)
  coresSpp <- setNames(object = viridis(n = length(especies),
                                        option = 'B'), nm = especies)
  for(sp in especies){
    pontos_MST <- vect(paste0('out/pointshape_', sp, '.shp'), crs = "+proj=longlat +datum=WGS84")
    MSTs <- vect(paste0('out/mst_', sp, '.shp'), crs = "+proj=longlat +datum=WGS84")
    plot(pontos_MST, add = T, col = coresSpp[sp])
    plot(MSTs, add = T, col = coresSpp[sp])
  }

  # Add legend
  par(xpd = TRUE)
  legend("bottomleft",
         legend = especies,
         pch = 19,
         col = unname(coresSpp[especies]),
         title = 'Non-homoplastic species',
         title.col = 'red',
         pt.cex = 0.8,
         cex = 0.7,
         bg = "white",
         box.lwd = 1)
  par(xpd = FALSE)
}
