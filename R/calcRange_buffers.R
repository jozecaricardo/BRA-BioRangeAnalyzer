#' calcRange_buffers
#'
#' Produces buffers of known area for one or more taxa on a map, generating presence-absence
#' matrices based on quadrats of user-defined resolution constructed within each buffer polygon.
#' Also produces a raster with species richness values based on polygon overlap.
#'
#' @param xy A data frame with occurrence points for one or more taxa, containing columns
#'   'spp', 'long', and 'lat' in decimal degrees; or an object produced by spatial
#'   analysis functions with a 'samples' component
#' @param buffer.width Numeric value specifying the buffer width in meters. If NULL and
#'   \code{mean_dist = TRUE}, uses the mean distance between points (default: NULL)
#' @param shape_file A polygon shapefile of a region with defined coordinate reference system (CRS)
#'   to be used as the study area boundary
#' @param resol Numeric value or vector specifying the resolution of quadrats in decimal degrees
#' @param mean_dist Logical, if TRUE calculates buffer width as the mean distance between points
#'   for each species; if FALSE uses the value provided in \code{buffer.width} (default: FALSE)
#' @return A list with two elements:
#'   \describe{
#'     \item{geometry}{A SpatRaster object with species richness values per cell}
#'     \item{pres_abs}{A presence-absence matrix based on quadrats of user-defined resolution,
#'       with species as columns and grid cells as rows}
#'   }
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function implements range estimation using circular buffers around occurrence
#' points. It provides an alternative to minimum convex polygons (MCP) by allowing explicit
#' control over the buffer area, which can be based on ecological or biological criteria.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters species with fewer than 3 occurrence records
#'   \item Removes species where all coordinates are identical (preventing degenerate geometries)
#'   \item Creates circular buffers around occurrence points for each species
#'   \item Converts buffers to a grid-based presence-absence matrix
#'   \item Calculates species richness by overlaying all buffer polygons
#' }
#'
#' Buffer width can be specified in two ways:
#' \itemize{
#'   \item \strong{Fixed width}: Set \code{buffer.width} to a specific value in meters
#'   \item \strong{Mean distance}: Set \code{mean_dist = TRUE} to use the mean distance
#'     between occurrence points as buffer width (species-specific)
#' }
#'
#' The function creates an 'out_buffers' directory containing:
#' \itemize{
#'   \item \strong{BUFF_*.shp}: Buffer polygon shapefiles for each species
#'   \item \strong{pointshape_*.shp}: Point occurrence shapefiles for each species
#'   \item \strong{presence_BUFF_*.tif}: Raster files with species presence
#'   \item \strong{pres_abs_BUFF_quad.txt}: Presence-absence matrix in text format
#' }
#'
#' \strong{Note:} The function assumes WGS84 coordinate reference system (EPSG:4326) and
#' issues a warning if coordinates are not in this format.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025.
#'   Panbiotracks: software for track analysis. \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'   Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1-17.
#'
#'   Gaston, K.J. & Fuller, R.A., 2009. The sizes of species' geographic ranges.
#'   \emph{Journal of Applied Ecology} \strong{46}: 1-9.
#' @seealso
#' \code{\link[terra]{buffer}}, \code{\link[raster]{rasterize}}, \code{\link[geosphere]{distm}}, \code{\link{extrapolation_plot_from_shapeFiles}}
#' @export
#' @examples
#' # Example 1: Basic usage with fixed buffer width (10 km)
#' data <- data.frame(
#'   species = rep(c("E_L_triangulator", "E_L_aceratos"), each = 5),
#'   long = c(-50, -51, -52, -53, -54, -45, -46, -47, -48, -49),
#'   lat = c(-10, -11, -12, -13, -14, -15, -16, -17, -18, -19)
#' )
#' result <- calcRange_buffers(data, buffer.width = 10000, shape_file = neo, resol = 15)
#'
#' # View results
#' plot(result$geometry, main = "Species Richness (Buffer Method)")
#' print(result$pres_abs)
#' #######################################################################################
#' # Example 2: Using mean distance between points as buffer width
#' result2 <- calcRange_buffers(data, shape_file = neo, resol = 15, mean_dist = TRUE)
#' #######################################################################################
#' # Example 3: Integration with extrapolation_plot_from_shapeFiles
#' result3 <- calcRange_buffers(data, buffer.width = 10000, shape_file = neo, resol = 15)
#' plot(neo, axes = TRUE, las = 1, main = "Species Buffers")
#' extrapolation_plot_from_shapeFiles(c("E_L_triangulator", "E_L_aceratos"),
#'              directory = "out_buffers", file_type = "BUFF")
#' #######################################################################################
#' # Doing a PAE-PCE analysis using buffers ...
#' rang_buff <- calcRange_buffers(xy = lycipta_final, shape_file = neo, resol = 10, buffer.width = 500000) # with 500 km buffers and 10 by 10 degree grid cells
#' pae_buff <- pae_pce(preabsMat = rang_buff$pres_abs, shapeFile = neo, resol = c(10, 10),
#'  gridView = TRUE, labelGrid = TRUE, nonHomoplasticSpeciesList = TRUE, N = 10, sobrepo = FALSE)
#'
#' # Species E. L. machadus and E. L. triangulator defined an area formed by grid cells 48 and 55
#' extrapolation_plot_from_shapeFiles(c("E_L_machadus", "E_L_triangulator"), directory = "out_buffers", file_type = "BUFF")
#'
#' # To view the occurrence points on the same map:
#' nomes <- unique(pae_buff$nonHomoplastic_species[[1]]$spp)
#'
#' # plotting the points
#' points_shapes_n <- list()
#' for(k in 1:length(nomes)){
#'   points_shapes_n[[k]] <- vect(paste0('out_buffers/pointshape_', nomes[k], '.shp'))
#' }
#' points_shapes_n
#'
#' # plotando os shapefiles no mapa
#' for(i in 1:length(points_shapes_n)){
#'   plot(points_shapes_n[[i]], add = T, col = cores)
#' }
#'
#' # Or...
#' extrapolation_plot_from_shapeFiles(c("E_L_machadus", "E_L_triangulator"), directory = "out_buffers",
#' file_type = "pointshape")
#'
#'
#'


calcRange_buffers <- function(xy, buffer.width = NULL, shape_file, resol, mean_dist = FALSE) {
  # projection
  wgs84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  warning("Assuming lat/long wgs84 coordinates")

  # check for geosphere package
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package 'geosphere' not found. Please install the package.", call. = FALSE)
  }

  make_valid_study_area <- function(shape_file_obj, context_label = "study area shapefile") {
    shape_sf <- if (inherits(shape_file_obj, "sf")) {
      shape_file_obj
    } else {
      sf::st_as_sf(shape_file_obj)
    }

    old_s2 <- sf::sf_use_s2()
    on.exit(sf::sf_use_s2(old_s2), add = TRUE)
    sf::sf_use_s2(FALSE)

    shape_sf <- tryCatch(
      sf::st_make_valid(shape_sf),
      error = function(e) {
        if (requireNamespace("lwgeom", quietly = TRUE)) {
          lwgeom::st_make_valid(shape_sf)
        } else {
          stop(
            paste0(
              "Unable to fix invalid ", context_label,
              ". Install package 'lwgeom' or repair the shapefile externally."
            )
          )
        }
      }
    )

    shape_sf <- suppressWarnings(sf::st_collection_extract(shape_sf, "POLYGON", warn = FALSE))
    if (nrow(shape_sf) == 0) {
      stop(paste0("No polygon geometries remain after validation of ", context_label, "."))
    }

    shape_sf <- suppressWarnings(sf::st_cast(shape_sf, "MULTIPOLYGON", warn = FALSE))
    as(shape_sf, "Spatial")
  }

  # fix different input data types data.frame
  if (is.data.frame(xy)) {
    names(xy) <- tolower(names(xy))
    dat <- xy[, c("spp", "long", "lat")]
  }
  ## spgeoOUt
  if (is.list(xy) && "samples" %in% names(xy)) {
    dat <- xy$samples[, 1:3]
    names(dat) <- c("spp", "long", "lat")
  }

  # check for species with less than 3 records
  filt <- table(dat$spp)
  # filt
  sortout <- names(filt[filt <= 2])
  # sortout
  filt <- filt[filt > 2]

  dat.filt <- droplevels(subset(dat, dat$spp %in% as.character(names(filt))))

  # check for species where all lat or long ar identical, or almost identical,
  # to prevent line polygons longitude
  test <- split(dat.filt, f = dat.filt$spp)
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallongitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$spp %in% sortout))

  # latitude
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallatitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$spp %in% sortout))
  #
  # # test for almost perfect fit
  # test2 <- sapply(test, function(k) {
  #   round(abs(cor(k[, "decimallongitude"], k[, "decimallatitude"])), 6)
  # })
  # sortout2 <- names(test2[test2 == 1])
  # sortout <- c(sortout, sortout2)
  # dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))
  #
  # sortout <- sortout[!is.na(sortout)]
  #
  # if (length(sortout) > 0) {
  #   warning("found species with < 3 occurrences:", paste("\n", sortout))
  # }

  print('1) preparing the coordinates...')
  coordin <- matrix(as.matrix(dat.filt[, c(2, 3)]), nrow(dat.filt),
                    2, dimnames = list(dat.filt[,1], colnames(dat.filt)[c(2, 3)]))

  cols1 <- setNames(viridis(n = length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))

  # preparing species' buffers...
  #######
  ### vector for accumulating specimens ###
  species <- unique(rownames(coordin))
  qde <- 0
  for(i in species){qde[i] <- length(which(rownames(coordin) == i))}
  qde <- cumsum(qde)
  tabela <- matrix(NA, nrow = dim(coordin), nc = 4)
  #######

  print('2) calculating buffers...')
  tempo <- coordin
  shape_spatial <- make_valid_study_area(shape_file)

  # tempo <- as.data.frame(na.exclude(tempo))
  # colnames(tempo) <- c('species', colnames(coordin[, c(1, 2)]))

  Long <- tempo[,1]
  Lat <- tempo[, 2]
  names(Long) <- rownames(tempo)

  tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(tempo), 2, dimnames = list(rownames(tempo),
                                                                                            colnames(tempo)[c(1,2)]))
  #######
  # preparando o raster com grid do shapefile escolhido
  grid <- raster(extent(shape_spatial), resolution = resol,
                 crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  #########
  # producing a raster of the shapefile
  mask.raster <- raster(extent(shape_spatial), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(shape_spatial, mask.raster)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  # mask.raster[is.na(mask.raster)] <- 0
  r <- merge(r, mask.raster)
  ncellras <- ncell(r)

  # create a folder for the output files of this course (if it doesn't already exist)
  if (!file.exists("out_buffers/")) dir.create("out_buffers/")

  ##############
  conta <- 0
  coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))),
                   dimnames = list(seq(1:ncellras), unique(rownames(tempo)))
                   [c(1:2)])
  lista_r <- list()
  current_run_tifs <- character(0)

  for(j in unique(rownames(tempo))){
    conta <- conta + 1
    tempoo <- subset(tempo, rownames(tempo) == j)

    #### shapefile ###
    tempo.d <- as.data.frame(tempoo)
    tempo_pts <- data.frame(
      long = as.numeric(tempo.d[, 1]),
      lat = as.numeric(tempo.d[, 2])
    )
    resul1_shape <- terra::vect(tempo_pts, geom = c("long", "lat"), crs = "EPSG:4326")
    writeVector(resul1_shape, paste0(c('out_buffers/pointshape_', j, '.shp'), collapse = ''), overwrite = TRUE)
    # resul1_shape <- vect(resul1_shape, crs = "+proj=longlat")

    ##### BUFFERS #####
    if(mean_dist == TRUE){
      dmat <- terra::distance(resul1_shape)
      dvals <- as.numeric(dmat)
      dvals <- dvals[is.finite(dvals) & dvals > 0]
      mean_dista <- if (length(dvals) > 0) mean(dvals) else NA_real_
      if (!is.finite(mean_dista) || mean_dista <= 0) {
        mean_dista <- if (!is.null(buffer.width) && is.finite(buffer.width) && buffer.width > 0) buffer.width else resol
      }
      obj_buff <- aggregate(buffer(resul1_shape, width = mean_dista))
    } else {
      obj_buff <- aggregate(buffer(resul1_shape, width = buffer.width))
    }
    # write.shapefile(hull, paste0(c('out_buffers/MCP_', j), collapse = ''))
    writeVector(obj_buff, paste0(c('out_buffers/BUFF_', j, '.shp'), collapse = ''), overwrite = T)

    print(paste0(c(conta + 2, ') calculating buffer... Done'), collapse = ''))

    # preparing shapes
    pasta <- paste0('BUFF_', j)
    pontos_linha <- shapefile(paste0(c('out_buffers/BUFF_', j, '.shp'), collapse = ''),
                              warnPRJ = FALSE)
    proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
    # plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
    # plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)

    idx <- which(names(qde) == j)
    tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- tempo[which(rownames(tempo) == j), 1]
    tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- tempo[which(rownames(tempo) == j), 2]
    tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))

    # back-transforming lines in points:
    # suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

    cover_r <- raster::rasterize(pontos_linha, r, getCover = TRUE)
    pres_vals <- raster::getValues(cover_r)
    valid_cells <- which(!is.na(r[]))
    hit_cells <- which(!is.na(pres_vals) & pres_vals > 0)

    coor.l[valid_cells, conta] <- 0
    if (length(hit_cells) > 0) {
      coor.l[hit_cells, conta] <- 1
    }

    lista_r[[conta]] <- r
    raster::values(lista_r[[conta]])[valid_cells] <- 0
    if (length(hit_cells) > 0) {
      raster::values(lista_r[[conta]])[hit_cells] <- 1
    }
    proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    tif_path <- paste0(c("out_buffers/presence_BUFF_q", resol, j, ".tif"), collapse = '')
    writeRaster(lista_r[[conta]], tif_path, overwrite = TRUE)
    current_run_tifs <- c(current_run_tifs, tif_path)

    # teste <- tabelao[j,j] # positions of taxa in the node
    teste <- conta
    # text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])

    # presence-absence matrix was filled by polygon-grid coverage intersection
  }

  coor.l <- na.exclude(coor.l)
  coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
  rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
  write.table(x = coor.ll, file = 'out_buffers/pres_abs_BUFF_quad.txt', sep = '\t')

  ##################
  print(paste0(conta + 3, ') Loading the raster files and adding the presences...'))
  current_run_tifs <- unique(current_run_tifs[file.exists(current_run_tifs)])
  if (length(current_run_tifs) == 0) {
    stop("No buffer raster outputs were generated in the current run.")
  }
  layers_MPC <- terra::rast(current_run_tifs)
  # Sum all layers to obtain species richness per cell
  MPC_soma <- sum(layers_MPC, na.rm = TRUE)
  return(list(geometry = MPC_soma, pres_abs = coor.ll))
}
