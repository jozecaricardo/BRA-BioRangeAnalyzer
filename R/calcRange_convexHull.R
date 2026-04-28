#' calcRange_convexHull
#'
#' Extrapolates distribution by generating minimum convex polygons from taxon occurrence points,
#' and produces a presence/absence matrix by cell.
#'
#' @param xy object produced by the function \code{\link[PanBioGeo]{singleton_to_data_frame}} or
#' a table with occurrence points of one or more taxa with columns 'spp', 'long' and 'lat' in decimal degrees
#' @param shape_file a polygon of some place already with datum to be used
#' @param resol a vector with resolution of quadrats
#' @return Presence/absence matrix based on quadrats - with user-chosen resolution - built
#' within each minimum convex polygon and a raster with richness values, based on polygon overlap.
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function promotes the generalization of distributions, from minimum convex polygons converted
#' into quadrats with user-chosen resolution. An 'out_MCP' directory is built and .shp files ('MCP_' for
#' convex hulls and 'pointshape_' for points), .tif (with richness) and .txt (the binary matrix)
#' with the polygons and occurrence points of taxa are placed in it.
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
#' \code{\link[PanBioGeo]{extrapolation_plot_from_shapeFiles}} for plotting species distributions.
#' @export
#' @examples
#' library(devtools)
#' library(tidyverse)
#' library(lwgeom)
#' library(tmap)
#' library(viridis)
#' library(readr)
#' library(sf)
#' library(phytools)
#' library(terra)
#'
#' # load example data (coordinates and tree)
#' data(lycipta)
#'
#' # The columns are in the following order: species, longitude, latitude
#' lycipta.coords <- lycipta$coordinates
#' head(lycipta.coords, 20)
#'
#' # shape file of South America
#' lycipta.asul <- lycipta$polygon
#'
#' # removal of singletons: the columns should be in the following order: species, latitude,
#' # longitude (without a tree)
#' resul <- singleton_to_data_frame_without_tree(spp = lycipta.coords$spp,
#' lat = lycipta.coords$Lat, long = lycipta.coords$Long)
#'
#' # important to change column names to long and lat
#' colnames(resul$data_df) <- c('species', 'long', 'lat')
#' rang1 <- calcRange_convexHull(xy = resul$data_df, shape_file = lycipta.asul,
#'                                resolut = c(1, 1)) # 1 degree quadrats
#' rang1$pres_abs # presence and absence matrix
#' rang1$geometry # raster with values
#'
#' # plot raster
#' cores <- viridis(n = 15, alpha = 0.5, option = 'B')
#' plot(lycipta.asul)
#' plot(rang1$geometry, col = cores, add = TRUE)
#'
#' # adding points
#' nomes <- unique(resul$data_df$species)
#' points_shapes_n <- list()
#' pasta <- 'pointshape_'
#' for(k in 1:length(nomes)){
#'   if(resul[k] < 3){ # ensuring we keep only taxa with three or more records
#'       next
#'   }
#'
#'   points_shapes_n[[k]] <- vect(paste0('out_MCP/pointshape_', names(resul)[k], '.shp'))
#' }
#' points_shapes_n <- points_shapes_n[-which(sapply(points_shapes_n, is.null))]
#'
#' # plotting shapefiles on the map
#' for(i in 1:length(points_shapes_n)){
#'   plot(points_shapes_n[[i]], add = TRUE, col = cores)
#' }
#'
#' # adding all minimum convex polygon shapefiles
#' MCP_shapes_n <- list()
#' pasta <- 'MCP_'
#' for(k in 1:length(nomes)){
#'   if(resul[k] < 3){ # ensuring we keep only taxa with three or more records
#'       next
#'   }
#'   MCP_shapes_n[[k]] <- vect(paste0('out_MCP/MCP_', names(resul)[k], '.shp'))
#' }
#' MCP_shapes_n <- MCP_shapes_n[-which(sapply(MCP_shapes_n, is.null))]
#'
#' # plotting shapefiles on the map
#' for(i in 1:length(MCP_shapes_n)){
#'   plot(MCP_shapes_n[[i]], add = TRUE, col = cores)
#' }
#' #######################################################################################
#' # -------- adding user-chosen polygon (e.g., E. L. imitator) ----------
#' plot(lycipta.asul)
#' # Search by part of the name
#' v_imitator_MCP <- MCP_shapes_n[[grep("E_L_imitator", sapply(MCP_shapes_n, sources))[1]]] # MCP
#' v_imitator_points <- points_shapes_n[[grep("E_L_imitator", sapply(points_shapes_n, sources))[1]]] # points
#'
#' plot(neo)
#' plot(v_imitator_MCP, add = TRUE, col = 'grey')
#' plot(v_imitator_points, add = TRUE, col = 'green')
#' ########################################################################################
#' # Integration with extrapolation_plot_from_shapeFiles
#' rang1 <- calcRange_convexHull(xy = resul$data_df, shape_file = lycipta.asul,
#'                                resol = c(1, 1)) # quadrats
#' plot(neo, axes = TRUE, las = 1, main = "Species MCPs")
#' extrapolation_plot_from_shapeFiles(c("E_L_triangulator", "E_L_aceratos"),
#'              directory = "out_MCP", file_type = "MCP")
#' ########################################################################################
#' # Doing a PAE-PCE analysis using convex hulls...
#' rang1 <- calcRange_convexHull(xy = resul$data_df, shape_file = lycipta.asul,
#'                                resol = c(5, 5)) # quadrats
#' pae_convexHulls <- pae_pce(preabsMat = rang1$pres_abs, shapeFile = lycipta.asul,
#'                            resol = c(5, 5), gridView = TRUE, labelGrid = TRUE,
#'                            nonHomoplasticSpeciesList = TRUE, N = 10, sobrepo = FALSE)
#'

calcRange_convexHull <- function(xy, shape_file, resol) {
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
    dat <- xy[, c("species", "long", "lat")]
  }
  ## spgeoOUt
  if (is.list(xy) && "samples" %in% names(xy)) {
    dat <- xy$samples[, 1:3]
    names(dat) <- c("species", "long", "lat")
  }

  # check for species with less than 3 records
  filt <- table(dat$species)
  # filt
  sortout <- names(filt[filt <= 2])
  # sortout
  filt <- filt[filt > 2]

  dat.filt <- droplevels(subset(dat, dat$species %in% as.character(names(filt))))

  # check for species where all lat or long are identical, or almost identical,
  # to prevent line polygons longitude
  test <- split(dat.filt, f = dat.filt$species)
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallongitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))

  # latitude
  test2 <- sapply(test, function(k) {
    length(unique(k$decimallatitude))
  })
  sortout2 <- names(test2[test2 == 1])
  sortout <- c(sortout, sortout2)
  dat.filt <- droplevels(subset(dat.filt, !dat.filt$species %in% sortout))

  print('1) preparing the coordinates...')
  coordin <- matrix(as.matrix(dat.filt[, c(2, 3)]), nrow(dat.filt),
                    2, dimnames = list(dat.filt[,1], colnames(dat.filt)[c(2, 3)]))

  cols1 <- setNames(viridis(n = length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))

  # preparing species' MCPs...
  #######
  ### accumulation vector of specimens ###
  species <- unique(rownames(coordin))
  qde <- 0
  for(i in species){qde[i] <- length(which(rownames(coordin) == i))}
  qde <- cumsum(qde)
  tabela <- matrix(NA, nrow = dim(coordin), ncol = 4)
  #######

  print('2) calculating MCP...')
  tempo <- coordin
  shape_spatial <- make_valid_study_area(shape_file)
  Long <- tempo[,1]
  Lat <- tempo[, 2]
  names(Long) <- rownames(tempo)

  tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(tempo), 2,
                  dimnames = list(rownames(tempo), colnames(tempo)[c(1,2)]))
  #######
  # producing the exact study-area raster used as the Convex Hull grid basis
  mask.raster <- raster(extent(shape_spatial), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  proj4string(mask.raster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  r <- raster::rasterize(shape_spatial, mask.raster, getCover = TRUE)
  r[is.na(r) | r <= 0] <- NA
  r[!is.na(r)] <- 1
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  r <- raster::merge(r, mask.raster)
  ncellras <- ncell(r)
  grid_index_r <- raster::setValues(raster::raster(r), seq_len(ncellras))
  gridPolygon <- raster::rasterToPolygons(grid_index_r)
  names(gridPolygon) <- "grid_id"
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  valid_cells <- which(!is.na(r[]))
  valid_grid <- gridPolygon[gridPolygon$grid_id %in% valid_cells, ]
  valid_grid_sf <- sf::st_as_sf(valid_grid)

  # create a folder for the output files (if it doesn't already exist)
  if (!file.exists("out_MCP/")) dir.create("out_MCP/")

  ##############
  conta <- 0
  coor.l <- matrix(NA, nrow = ncellras, ncol = length(unique(rownames(tempo))),
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
    writeVector(resul1_shape, paste0(c('out_MCP/pointshape_', j, '.shp'), collapse = ''), overwrite = TRUE)
    # resul1_shape <- vect(resul1_shape, crs = "+proj=longlat")

    ##### MCP #####
    hull <- convHull(resul1_shape)
    # write.shapefile(hull, paste0(c('out_MCP/MCP_', j), collapse = ''))
    writeVector(hull, paste0(c('out_MCP/MCP_', j, '.shp'), collapse = ''), overwrite = TRUE)

    print(paste0(c(conta + 2, ') calculating MCP... Done'), collapse = ''))

    # preparing shapes
    pasta <- paste0('MCP_', j)
    pontos_linha <- shapefile(paste0(c('out_MCP/MCP_', j, '.shp'), collapse = ''),
                              warnPRJ = FALSE)
    proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
    # plot(pontos_linha, col = cols1[j], lwd = 3, add = TRUE)
    # plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = TRUE)

    idx <- which(names(qde) == j)
    tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- tempo[which(rownames(tempo) == j), 1]
    tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- tempo[which(rownames(tempo) == j), 2]
    tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))

    # back-transforming lines in points:
    # suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

    hull_sf <- sf::st_as_sf(pontos_linha)
    old_s2_loop <- sf::sf_use_s2()
    on.exit(sf::sf_use_s2(old_s2_loop), add = TRUE)
    sf::sf_use_s2(FALSE)
    hull_sf <- suppressWarnings(sf::st_make_valid(hull_sf))

    hit_cells <- tryCatch({
      hit_idx <- lengths(sf::st_intersects(valid_grid_sf, hull_sf)) > 0
      valid_grid_sf$grid_id[hit_idx]
    }, error = function(e) integer(0))
    if (length(hit_cells) == 0) {
      cover_r <- tryCatch(
        raster::rasterize(pontos_linha, r, getCover = TRUE),
        error = function(e) NULL
      )
      if (!is.null(cover_r)) {
        pres_vals <- raster::getValues(cover_r)
        hit_cells <- which(!is.na(pres_vals) & pres_vals > 0)
      }
    }

    point_cells <- raster::cellFromXY(r, as.matrix(tempo_pts[, c("long", "lat")]))
    point_cells <- unique(point_cells[!is.na(point_cells)])
    hit_cells <- sort(unique(c(hit_cells, point_cells)))
    hit_cells <- hit_cells[hit_cells %in% valid_cells]

    coor.l[valid_cells, conta] <- 0
    if (length(hit_cells) > 0) {
      coor.l[as.character(hit_cells), conta] <- 1
    }

    lista_r[[conta]] <- r
    raster::values(lista_r[[conta]])[valid_cells] <- 0
    if (length(hit_cells) > 0) {
      raster::values(lista_r[[conta]])[hit_cells] <- 1
    }
    proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    tif_path <- paste0(c("out_MCP/presence_MCP_q", resol, j, ".tif"), collapse = '')
    writeRaster(lista_r[[conta]], tif_path, overwrite = TRUE)
    current_run_tifs <- c(current_run_tifs, tif_path)

    # teste <- tabelao[j,j] # taxon positions of the node
    teste <- conta
    # text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])

    # presence-absence matrix was filled by polygon-grid coverage intersection
  }

  coor.l <- na.exclude(coor.l)
  coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
  rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
  write.table(x = coor.ll, file = 'out_MCP/pres_abs_MCP_quad.txt', sep = '\t')

  occupied_ids <- suppressWarnings(as.integer(rownames(coor.l)[rowSums(coor.l, na.rm = TRUE) > 0]))
  occupied_ids <- unique(occupied_ids[!is.na(occupied_ids)])
  occupied_grid <- if (length(occupied_ids) > 0) {
    valid_grid[valid_grid$grid_id %in% occupied_ids, ]
  } else {
    valid_grid[0, ]
  }

  resol_tag <- gsub("[^0-9A-Za-z]+", "_", format(resol, scientific = FALSE, trim = TRUE))
  if (!nzchar(resol_tag)) resol_tag <- "grid"
  terra::writeVector(terra::vect(valid_grid), paste0('out_MCP/GRIDS_allcells_convex_hull_q', resol_tag, '.shp'), overwrite = TRUE)
  if (nrow(occupied_grid) > 0) {
    terra::writeVector(terra::vect(occupied_grid), paste0('out_MCP/GRIDS_presence_convex_hull_q', resol_tag, '.shp'), overwrite = TRUE)
  }

  ##################
  print(paste0(conta + 3, ') Loading the raster files and adding the presences...'))
  current_run_tifs <- unique(current_run_tifs[file.exists(current_run_tifs)])
  if (length(current_run_tifs) == 0) {
    stop("No convex-hull raster outputs were generated in the current run.")
  }
  layers_MCP <- terra::rast(current_run_tifs)
  # Sum all layers to get species richness per cell
  MCP_soma <- sum(layers_MCP, na.rm = TRUE)
  return(list(
    geometry = MCP_soma,
    pres_abs = coor.ll,
    grid_all_sf = sf::st_as_sf(valid_grid),
    grid_presence_sf = sf::st_as_sf(occupied_grid)
  ))
}
