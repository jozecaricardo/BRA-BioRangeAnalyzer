#' MST_node
#'
#' Estimates the minimum spanning tree (MST) for a given taxon, taxa of a given node or nodes,
#' both internal and terminal. It also estimates the MST of a given species regardless of its presence
#' on a tree.
#'
#' @param coordin coordinates - in decimal degrees - with the following header: spp, lat,
#'  long (all lowercase)
#' @param tree a phylogenetic tree in Newick format
#' @param shape_file a polygon of some place already with datum
#' @param seeres logical value that enables seeing the resolution of quadrats on the map
#' @param resol a vector with resolution of quadrats
#' @param sobrepo logical value that enables plotting on a map a new minimum spanning trees (MSTs) have already
#' been produced
#' @param mintreeall logical value that enables plotting a single minimum spanning tree (MST) that represents
#' all species
#' @param caption logical value that enables the labeling of a plot based on the specified species
#' @param nodes a vector with the node(s) of a tree that allows for plotting a single minimum spanning tree (MST)
#' representing all species related to that (those) node(s) (optional)
#' @param taxon a vector containing the name(s) of terminal taxon (taxa) for creating a minimum spanning tree (MST)
#' for each of them
#' @param transp a value (from 0 to 1) representing the degree of transparency of a polygon
#' @param seephylog a logical option to visualize the phylogeny with nodes on the map (optional)
#' @param ... optional plot arguments
#' @return Matrix with a number of rows equal to the number of cells of the polygon, determined by
#' the resol parameter, and with a number of columns equal to the number of taxa.
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function produces individual tracks, which are the basic unit of a panbiogeographical
#'  analysis (Page, 1987; Castillo-Garcia \emph{et al.}, 2025). It defines a minimum spanning tree,
#'  which is a graph with \emph{n} localities connected through \emph{n}-1 edges, producing a shortest
#'  total length graph.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025. Panbiotracks: software for track analysis.
#'  \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'  Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1.
#' @seealso
#' \code{\link[phytools]{phytools}} for phylogenetic trees,
#' \code{\link[rnaturalearthdata]{rnaturalearthdata}} for global geographic data,
#' \code{\link[terra]{terra}} for spatial data handling.
#' @export
#' @examples
#' # load example data (coordinates and tree)
#' library(phytools)
#' library(terra)
#' library(sf)
#' data(lycipta)
#'
#' # The columns are in the following order: species, longitude, latitude
#' lycipta.coords <- lycipta$coordinates
#' head(lycipta.coords, 20)
#'
#' # phylogenetic tree
#' lycipta.tree <- lycipta$phylogenetic_tree
#'
#' # removal of singletons: the columns should be in the following order: species, latitude,
#' # longitude (with a tree)
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)], phylogeny = lycipta.tree)
#'
#' #' # removal of singletons: the columns should be in the following order: species, latitude,
#' # longitude (without a tree)
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)])
#'###############################################################################################
#' # Download South America countries (low or high resolutions)
#' sa_countries <- rnaturalearth::ne_countries(continent = "South America", scale = "medium",
#'   returnclass = "sf")
#'
#' sa_countries <- rnaturalearth::ne_countries(continent = "South America", scale = "large",
#'   returnclass = "sf")
#'
#' brazil_states <- rnaturalearth::ne_states(country = "Brazil", returnclass = "sf")
#'
#' # removing Brazil
#' sa_no_br <- sa_countries[sa_countries$iso_a3 != "BRA", ]
#'
#' # standardizing columns (keeping geometry column intact)
#' countries <- sa_no_br[, c("name", "iso_a3")]
#' names(countries)[1:2] <- c("region_name", "iso_code")
#' countries$type <- "country"
#'
#' states <- brazil_states[, c("name", "iso_3166_2")]
#' names(states)[1:2] <- c("region_name", "iso_code")
#' states$type <- "state"
#'
#' st_is_valid(countries)
#' st_is_valid(states) # there is a problem here!
#'
#' # Validating gemometries
#' countries_clean <- st_make_valid(countries)
#' states_clean <- st_make_valid(states)
#'
#' # Combining
#' sa_bins <- rbind(countries_clean, states_clean)
#'
#' # Convert to terra SpatVector if needed
#' sa_bins_sv <- terra::vect(sa_bins)
#'
#'
#' # ---------------------------- with all taxa ------------------------------
#' # ----------- with a tree and resolution of 10 X 10 degrees ----------
#' mst_all_taxa <- MST_node(coordin = resul$data_df, shape_file = sa_bins_sv, sobrepo = FALSE, caption = TRUE,
#'   resol = c(10, 10), tree = resul$treeMod)
#'
#' # -------- without a tree and resolution of 10 X 10 degrees ----------
#' mst_all_taxa <- MST_node(coordin = resul$data_df, shape_file = sa_bins_sv, sobrepo = FALSE, caption = TRUE,
#'   resol = c(10, 10))
#' #################################################################################################
#' # ---------------- with terminals: E. l. triangulator ----------------
#' resul_triangulator <- MST_node(taxon = 'E_L_triangulator', coordin = resul$data_df, shape_file = sa_bins_sv,
#'   sobrepo = FALSE, caption = TRUE, resol = c(10, 10), tree = resul$treeMod)
#'
#' # --------------- with terminals: E. l. triangulator + E. l. aceratos
#' resul_twoTaxa <- MST_node(taxon = c('E_L_triangulator', 'E_L_aceratos'), coordin = resul$data_df,
#'   shape_file = sa_bins_sv, sobrepo = FALSE, caption = TRUE, resol = c(10, 10), tree = resul$treeMod)
#' #################################################################################################
#' # ----------------------------- with internal nodes --------------------------------
#' # To determine the minimum spanning tree of members in an internal node, we first need to know the number of internal
#' # nodes in the tree.
#'  plotTree(resul$treeMod)
#'  labelnodes(1:(Ntip(resul$treeMod) + resul$treeMod$Nnode), 1:(Ntip(resul$treeMod) + resul$treeMod$Nnode),
#'   interactive = F, cex = 0.6, circle.exp = 0.4)
#'
#' # Now that we know the numbers of the internal nodes, we will ask R to calculate the individual tracks of all
#' # members of any monophyletic group.
#' # internal nodes: node number 23
#' resul_node_23 <- MST_node(coordin = resul$data_df, tree = resul$treeMod, shape_file = sa_bins_sv,
#'  caption = TRUE, resol = c(10, 10), nodes = 23, sobrepo = FALSE)
#' ##################################################################################################
#' # ------------------- internal nodes + terminal nodes (E_L_imitator) ---------------
#' resul_node23_imitator <- MST_node(coordin = resulcTree$data_df, tree = resulcTree$treeMod, shape_file = asul,
#' caption = TRUE, resol = c(10, 10), nodes = 23, taxon = 'E_L_imitator')
#' ##################################################################################################
#' # the map has been plotted! Now, let's ask R to calculate the individual track of a taxon at a node
#' resul_node20_sobrepo <- MST_node(coordin = resul$data_df, tree = resul$treeMod, shape_file = asul,
#' caption = TRUE, resol = c(10, 10), sobrepo = T, nodes = 20)
#'


MST_node <- function(coordin, tree = NULL, shape_file, resol, seeres = FALSE, sobrepo = FALSE, mintreeall = FALSE, caption = TRUE,
                     nodes = NULL, taxon = NULL, pol = FALSE, transp = 0.5, seephylog = FALSE, xmin = NULL, xmax = NULL,
                     ymin = NULL, ymax = NULL){

  print('1) loading functions...')
  library(geiger); library(devtools)
  library(phytools)
  library(letsR); library(fossil)  # Package maptools was removed from the CRAN repository.
  library(viridis); library(mapdata); library(vegan); library(spdep)
  library(dismo); library(spatstat); library(terra)

  validate_mst_matrix <- function(mst_obj, labels, context_label = "selected subset") {
    if (is.null(dim(mst_obj)) || length(dim(mst_obj)) != 2 || any(dim(mst_obj) < 2)) {
      stop(paste0(
        "MST could not be generated for ", context_label,
        ": fewer than two distinct occurrence points are available after filtering."
      ))
    }
    if (length(labels) != nrow(mst_obj)) {
      stop(paste0(
        "MST failed for ", context_label,
        ": point labels are inconsistent with the MST matrix dimensions."
      ))
    }
    rownames(mst_obj) <- labels
    colnames(mst_obj) <- labels
    mst_obj
  }

  make_valid_polygon_sf <- function(x, context_label = "geometry") {
    old_s2 <- sf::sf_use_s2()
    on.exit(sf::sf_use_s2(old_s2), add = TRUE)
    sf::sf_use_s2(FALSE)

    x <- tryCatch(
      sf::st_make_valid(x),
      error = function(e) {
        if (requireNamespace("lwgeom", quietly = TRUE)) {
          lwgeom::st_make_valid(x)
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

    x <- suppressWarnings(sf::st_collection_extract(x, "POLYGON", warn = FALSE))
    if (nrow(x) == 0) {
      stop(paste0("No polygon geometries remain after validation of ", context_label, "."))
    }

    x <- suppressWarnings(sf::st_cast(x, "MULTIPOLYGON", warn = FALSE))
    x
  }

  # Store original shapefile at the beginning (before any plot calls)
  shape_file_original <- shape_file

  if(is.null(taxon) == FALSE){
    if(length(taxon) < 2){
      warning('It is not possible with only one species to produce generalized tracks by the PAE-PCE analysis!')
    }
  }

  if(sobrepo == FALSE){
    par(mfrow = c(1, 1))
    # Use original shapefile for plotting (shows all polygons)
    plot(shape_file_original, type = 'n', xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  }

  print('2) preparing the coordinates...')
  coordin <- matrix(as.matrix(coordin[, c(2, 3)]), nrow(coordin),
                    2, dimnames = list(coordin[,1], colnames(coordin)[c(2, 3)]))

  cols1 <- setNames(viridis(n = length(unique(rownames(coordin)))),
                    unique(rownames(coordin)))

  if(seephylog == TRUE && sobrepo == TRUE){
    stop('This is not possible! Please try to use the each argument separately.')
  }


  if(seeres == TRUE){
    layout(matrix(1, nr = 1, nc = 1, byrow = F))
  }
  if(seephylog == TRUE && seeres == TRUE){
    layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = F))
  }

  if(!is.null(tree)){
    print('3) preparing the tree...')
    # phylogenetic trees and nodes:
    if(is.null(tree$edge.length) == TRUE){
      print('...creating the branch lengths of a tree equal to one...')
      tree <- compute.brlen(tree, 1)
    }

    print('...checking names from dataset and names in the species tree.')
    if(name.check(tree, coordin)[1] != 'OK'){
      chk <- name.check(tree, coordin)
      tree <- drop.tip(tree, chk$tree_not_data)
    }
  }

  # Validate and fix geometries to prevent TopologyException
  print('Validating shapefile geometries...')

  if(inherits(shape_file, 'SpatVector')){
    # Convert to sf for validation
    shape_sf <- sf::st_as_sf(shape_file)
    shape_sf <- make_valid_polygon_sf(shape_sf, context_label = "input shapefile")
    shape_file <- terra::vect(shape_sf)
    print(paste('  Geometries validated. Features kept for raster operations:', nrow(shape_sf)))
  }

  # plotting the tree and coordinates:
  if(is.na(proj4string(as(shape_file, 'Spatial'))) == TRUE){
    print('Attention: your shapefile has not been associated to any datum and this routine is going
            to associate it to WGS84 datum!')
    # proj4string(shape_file) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    # projection(as(shape_file, 'Spatial')) <- CRS("+proj=longlat +datum=WGS84")
    crs(shape_file) <- "+proj=longlat +datum=WGS84"
  } else if(proj4string(as(shape_file, 'Spatial')) != "+proj=longlat +datum=WGS84"){
    # proj4string(shape_file) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    # projection(as(shape_file, 'Spatial')) <- CRS("+proj=longlat +datum=WGS84")
    crs(shape_file) <- "+proj=longlat +datum=WGS84"
  }

  if(seephylog == TRUE){
    obj2 <- phylo.to.map(tree = tree, coords = coordin[,2:1],
                         database = as(shape_file, 'Spatial'), plot = F, rotate = F)
    plot(obj2, colors = cols1, ftype = 'i', fsize = 0.6, cex.points = c(0.7, 1.2), pts = F,
         direction = "rightwards")
    labelnodes(1:(Ntip(tree) + tree$Nnode), 1:(Ntip(tree) + tree$Nnode),
               interactive = F, cex = .6, circle.exp = 0.4)
    mtext('Phylogeny on the map', side = 3, line = -5)
  }

  if(!is.null(tree)){
    print('4) dividing polygon into grid and producing a raster...')
  } else {
    print('3) dividing polygon into grid and producing a raster...')
  }

  if(!is.null(tree)){
    arv <- reorder(tree, "postorder") # reordering the levels
    e1 <- arv$edge[, 1] # internal nodes
    e2 <- arv$edge[, 2] # terminal nodes and root
    EL <- arv$edge.length # branch lengths

    tabelao <- mrca(phy = tree, full = F)
  }


  grid <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                 crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84"


  # clipping the intersected cells:
  # Convert to Spatial and validate geometries to prevent TopologyException
  shape_spatial <- as(shape_file, 'Spatial')

  # Validate Spatial object using rgeos (if available) or sf
  if(requireNamespace("rgeos", quietly = TRUE)){
    print('  Cleaning Spatial geometries with rgeos...')
    shape_spatial <- rgeos::gBuffer(shape_spatial, byid = TRUE, width = 0)
  } else {
    print('  Cleaning Spatial geometries with sf...')
    shape_sf_temp <- sf::st_as_sf(shape_spatial)
    shape_sf_temp <- make_valid_polygon_sf(shape_sf_temp, context_label = "raster mask shapefile")
    shape_sf_temp <- sf::st_buffer(shape_sf_temp, dist = 0)
    shape_spatial <- as(shape_sf_temp, 'Spatial')
  }

  # Build a single grid mask using polygon coverage (> 0) to keep consistency
  # across matrix and map (also keeps tiny polygon slivers represented).
  mask.raster <- raster(extent(as(shape_file, 'Spatial')), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  cover.r <- rasterize(shape_spatial, mask.raster, getCover = TRUE)
  cover.r <- merge(cover.r, mask.raster)
  cover.r[is.na(cover.r) | cover.r <= 0] <- NA
  cover.r[!is.na(cover.r)] <- 1
  r <- cover.r
  cropped_map <- rasterToPolygons(r, dissolve = FALSE)
  if (seeres == TRUE){
    plot(cropped_map, xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = T)
    proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    ncellR <- ncell(r)
    map.r <- raster::as.data.frame(raster::rasterToPoints(r))
    pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r, field = 1)
    proj4string(pontosRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84(pontosRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
    # plot(pontosRaster, add = T)
    map.r$gridNumber <- which(pontosRaster@data@values == 1)
    text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8,
         col = adjustcolor(col = 'red', alpha = 1), font = 2)
  }

  # r already built above from polygon coverage mask (>= tiny intersections)
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  ncellras <- ncell(r)
  grid_index_r <- raster::setValues(raster::raster(r), seq_len(ncellras))
  gridPolygon_idx <- raster::rasterToPolygons(grid_index_r)
  names(gridPolygon_idx) <- "grid_id"
  crs(gridPolygon_idx) <- "+proj=longlat +datum=WGS84"
  valid_cells <- which(!is.na(r[]))
  valid_grid <- gridPolygon_idx[gridPolygon_idx$grid_id %in% valid_cells, ]
  valid_grid_sf <- sf::st_as_sf(valid_grid)

  mst_hit_cells <- function(track_spatial, point_coords = NULL, context_label = "MST geometry") {
    hit_cells <- integer(0)

    old_s2_loop <- sf::sf_use_s2()
    on.exit(sf::sf_use_s2(old_s2_loop), add = TRUE)
    sf::sf_use_s2(FALSE)

    track_sf <- tryCatch(
      sf::st_as_sf(track_spatial),
      error = function(e) NULL
    )

    if (!is.null(track_sf)) {
      track_sf <- tryCatch(
        sf::st_make_valid(track_sf),
        error = function(e) track_sf
      )
      hit_cells <- tryCatch({
        hit_idx <- lengths(sf::st_intersects(valid_grid_sf, track_sf)) > 0
        valid_grid_sf$grid_id[hit_idx]
      }, error = function(e) integer(0))
    }

    if (length(hit_cells) == 0) {
      cover_r <- tryCatch(
        raster::rasterize(track_spatial, r, getCover = TRUE),
        error = function(e) NULL
      )
      if (!is.null(cover_r)) {
        pres_vals <- raster::getValues(cover_r)
        hit_cells <- which(!is.na(pres_vals) & pres_vals > 0)
      }
    }

    if (!is.null(point_coords)) {
      point_xy <- as.matrix(point_coords)
      if (!is.null(dim(point_xy)) && ncol(point_xy) >= 2) {
        point_xy <- point_xy[, 1:2, drop = FALSE]
        point_cells <- raster::cellFromXY(r, point_xy)
        point_cells <- unique(point_cells[!is.na(point_cells)])
        hit_cells <- sort(unique(c(hit_cells, point_cells)))
      }
    }

    hit_cells <- hit_cells[hit_cells %in% valid_cells]
    hit_cells
  }

  build_mst_presence_raster <- function(hit_cells) {
    out_r <- r
    raster::values(out_r)[valid_cells] <- 0
    if (length(hit_cells) > 0) {
      raster::values(out_r)[hit_cells] <- 1
    }
    proj4string(out_r) <- CRS("+proj=longlat +datum=WGS84")
    out_r
  }

  if (!file.exists("out_MST/")) dir.create("out_MST/") # temporary folder

  if(sobrepo == FALSE){
    if(seephylog == TRUE){
      if(is.null(nodes)){
        plot(shape_file_original, axes = TRUE, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
        # abline(h = 0, col = 'red', lty = 2) # equator
        mtext('Minimum spanning trees of the terminal nodes', side = 3, line = 2)
      } else if(!is.null(nodes)){
        if(length(nodes) == 1){
          plot(shape_file_original, axes = TRUE, cex.main = 0.7, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
          mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
          # abline(h = 0, col = 'red', lty = 2) # equator
          mtext(paste0(c('Mapping a minimum spanning tree of the internal node:',
                         nodes), collapse = ' '), side = 3, line = 2)
        } else if(length(nodes) > 1){
          plot(shape_file_original, axes = TRUE, cex.main = 0.7,
               xlim = c(xmin, xmax), ylim = c(ymin, ymax))
          mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
          # abline(h = 0, col = 'red', lty = 2) # equator
          mtext(paste0(c('Mapping a minimum spanning tree of the internal nodes:',
                         nodes), collapse = ' '), side = 3, line = 2)
        }
      }
    } else {
      if(is.null(nodes) && !is.null(taxon)){
        plot(shape_file_original, axes = TRUE, main = 'Minimum spanning tree(s) of the terminal node(s)',
             cex.main = 0.9, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
        # abline(h = 0, col = 'red', lty = 2) # equator
      } else if(!is.null(nodes) || !is.null(taxon)){
        plot(shape_file_original, axes = TRUE, main = paste0(c('Mapping a minimum spanning tree of the node(s):',
                                                               nodes), collapse = ' '), cex.main = 0.9, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
        # abline(h = 0, col = 'red', lty = 2) # equator
      } else if(is.null(nodes) && is.null(taxon)){
        plot(shape_file_original, axes = TRUE, main = 'Mapping a minimum spanning tree of the terminal nodes!',
             cex.main = 0.9, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
        mtext(paste0('resolution: ', resol[1], ' x ', resol[2]), side = 1, line = 3, cex = 0.7)
        # abline(h = 0, col = 'red', lty = 2) # equator
      }
    }
  }


  # preparing species' msts...
  #######
  ### vector for accumulating specimens ###
  species <- unique(rownames(coordin))
  qde <- 0
  for(i in species){qde[i] <- length(which(rownames(coordin) == i))}
  qde <- cumsum(qde)
  tabela <- matrix(NA, nrow = dim(coordin), nc = 4)
  #######

  if(!is.null(tree)){
    print('5) calculating mst...')
  } else {
    print('4) calculating mst...')
  }
  if(is.null(nodes)){

    if(!is.null(taxon)){

      if(mintreeall == TRUE){

        for(j in taxon){
          # sppp <- which(unique(rownames(coords)) == taxon[j])

          ##### If we want to make MST of all specimens TOGETHER from the chosen species...
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
        }
        frame <- as.data.frame(na.exclude(tabela))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))

        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]

        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2, dimnames = list(frame[,1],
                                                                                                  colnames(frame)[c(2,3)]))

        #### shapefile ###
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        # dir.create('out_MST/')
        write.shapefile(tempo_shape, 'out_MST/pointsshape_mintreeall')
        # resul1_shape <- rgdal::readOGR(dsn = 'out_MST/pointsshape_mintreeall.shp', verbose = FALSE)
        # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
        resul1_shape <- vect('out_MST/pointsshape_mintreeall.shp',
                             crs = "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        terra::project(resul1_shape, "+proj=longlat +datum=WGS84")


        ##### MST based on geographic distance #####
        # rownames(resul1_shape@coords) <- rownames(tempo.d)
        # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
        mst2 <- dino.mst(dista)
        mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = "mintreeall selected taxa")
        lats <- cbind(resul1_shape$long,
                      resul1_shape$lat)
        rownames(lats) <- rownames(tempo.d)
        colnames(lats) <- c('longitude', 'latitude')
        mst_shape <- msn2Shape(msn = mst2, lats = lats)
        # mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
        write.shapefile(mst_shape, 'out_MST/mst_mintreeall')

        if(!is.null(tree)){
          print('5) calculating mst... Done')
        } else {
          print('4) calculating mst... Done')
        }

        # plotting shapes
        pontos_linha <- shapefile('out_MST/mst_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        crs(pontos_linha) <- "+proj=longlat +datum=WGS84"

        plot(pontos_linha, col = cols1[1], lwd = 3, add = T)
        plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[1:length(taxon)], add = T)

        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor(cols1[taxon[1]], transp))
          conv <- st_as_sf(conv)
          conv <- as_Spatial(conv)
          conv <- vect(conv, crs = "+proj=longlat +datum=WGS84")
          writeVector(conv, 'out_MST/mcp_mintreeall', overwrite = T)
        }

        # legend
        if(caption == TRUE){
          x <- taxon
          if(sobrepo == TRUE){
            legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x, pch = 19, col = cols1[x], bty = 'n',
                   pt.cex = 0.6, cex = 0.6, title = 'Adding species...', title.col = 'red')
          } else {
            legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = 'Members of the minimum spanning tree',
                   title.col = 'red', pt.cex = 0.6, cex = 0.6)
          }
        }

        hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo, context_label = "mintreeall selected taxa")

        # presence-absence matrix:
        ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(taxon), dimnames = list(seq(1:ncellras),
                                                                                unique(rownames(tempo))[c(1:length(taxon))]))  # tabela com todas as celulas
        linhasRaster <- build_mst_presence_raster(hit_cells)

        writeRaster(linhasRaster, "out_MST/presence_mintreeall.tif",
                    overwrite = TRUE)

        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
             alpha = transp)

        coor.l[valid_cells, ] <- 0
        if (length(hit_cells) > 0) {
          coor.l[hit_cells, ] <- 1
        }

        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

        return(coor.ll)


      } else if (mintreeall == FALSE){

        conta <- 0
        coor.l <- matrix(NA, nr = ncellras, nc = length(taxon), dimnames = list(seq(1:ncellras),
                                                                                taxon)[c(1:2)])  # tabela com todas as celulas
        lista_r <- list()

        for(j in taxon){
          conta <- conta + 1
          tempo <- subset(coordin[, 1:2], rownames(coordin) == j)

          if (nrow(tempo) < 2 || nrow(unique(tempo[, 1:2, drop = FALSE])) < 2) {
            warning(paste0("Skipping taxon '", j, "': fewer than two distinct occurrence points for MST."))
            coor.l[, conta] <- 0
            next
          }

          #### shapefile ###
          tempo.d <- as.data.frame(tempo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          # dir.create('out_MST/')
          write.shapefile(tempo_shape, paste0(c('out_MST/pointshape_', j), collapse = ''))
          resul1_shape <- vect(paste0(c('out_MST/pointshape_', j, '.shp'),
                                      collapse = ''), crs = "+proj=longlat +datum=WGS84")
          # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
          # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84

          ##### MST based on geographic distance #####
          # rownames(resul1_shape@coords) <- rownames(tempo.d)
          # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
          mst2 <- dino.mst(dista)
          mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = paste0("taxon '", j, "'"))
          lats <- cbind(resul1_shape$long,
                        resul1_shape$lat)
          rownames(lats) <- rownames(tempo.d)
          colnames(lats) <- c('longitude', 'latitude')
          mst_shape <- msn2Shape(msn = mst2, lats = lats)
          # mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
          write.shapefile(mst_shape, paste0(c('out_MST/mst_', j), collapse = ''))

          if(!is.null(tree)){
            print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
          } else {
            print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
          }

          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out_MST/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
          plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)

          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempo)
            conv <- st_as_sf(conv)
            conv <- as_Spatial(conv)
            conv <- vect(conv, crs = '+proj=longlat +datum=WGS84')
            writeVector(conv, paste0(c('out_MST/mcp_', j), collapse = ''), overwrite = T)
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
          }

          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))


          hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo.d[, 1:2, drop = FALSE], context_label = paste0("taxon '", j, "'"))

          lista_r[[conta]] <- build_mst_presence_raster(hit_cells)
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE,
               col = cols1[j], alpha = transp)
          # colocando labels:
          teste <- tabelao[j,j] # posicoes dos taxons do no
          text(lista_r[[conta]], labels = rep(teste, dim(tempo)[1]), cex = 0.8, pos = 2, col = cols1[j])

          coor.l[valid_cells, conta] <- 0
          if (length(hit_cells) > 0) {
            coor.l[hit_cells, conta] <- 1
          }
        }

        if(caption == TRUE){
          x <- taxon
          if(sobrepo == TRUE){
            legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x, pch = 19, col = cols1[x],
                   pt.cex = 0.6, cex = 0.6, title = 'Adding species...', title.col = 'red')
          } else {
            # Create vector with track numbers for each taxon
            track_numbers <- sapply(taxon, function(j) tabelao[j,j])
            # Create legends with species name and track number in parentheses
            legend_labels <- paste0(taxon, " (", track_numbers, ")")

            legend(x = "bottomleft", legend = legend_labels, pch = 19, col = cols1[taxon],
                   title = 'Members of the minimum spanning tree(s)',
                   title.col = 'red', pt.cex = 0.6, cex = 0.6)
          }
        }

        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

        return(coor.ll)

      }
    } else if(is.null(taxon)){ # null taxa

      tempo <- coordin

      # tempo <- as.data.frame(na.exclude(tempo))
      # colnames(tempo) <- c('species', colnames(coordin[, c(1, 2)]))

      Long <- tempo[,1]
      Lat <- tempo[, 2]
      names(Long) <- rownames(tempo)

      tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(tempo), 2, dimnames = list(rownames(tempo),
                                                                                                colnames(tempo)[c(1,2)]))


      if(mintreeall == TRUE){

        #### shapefile ###
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        write.shapefile(tempo_shape, 'out_MST/pointsshape_mintreeall')
        resul1_shape <- vect('out_MST/pointsshape_mintreeall.shp',
                             crs = "+proj=longlat +datum=WGS84")
        # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
        terra::project(resul1_shape, "+proj=longlat +datum=WGS84") # datum WGS84

        ##### MST based on geographic distance #####
        # rownames(resul1_shape@coords) <- rownames(tempo.d)
        # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
        mst2 <- dino.mst(dista)
        mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = "mintreeall all taxa")
        lats <- cbind(resul1_shape$long,
                      resul1_shape$lat)
        rownames(lats) <- rownames(tempo.d)
        colnames(lats) <- c('longitude', 'latitude')
        mst_shape <- msn2Shape(msn = mst2, lats = lats)
        write.shapefile(mst_shape, 'out_MST/mst_mintreeall')

        if(!is.null(tree)){
          print('5) calculating mst... Done')
        } else {
          print('4) calculating mst... Done')
        }

        # plotting shapes
        pontos_linha <- shapefile('out_MST/mst_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        plot(pontos_linha, col = cols1[unique(rownames(tempo))[1]], lwd = 3, add = T)
        plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[unique(rownames(tempo))[1]], add = T)

        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor(cols1[unique(rownames(tempo))[1]], transp))
          conv <- st_as_sf(conv)
          conv <- as_Spatial(conv)
          conv <- vect(conv, crs = '+proj=longlat +datum=WGS84')
          writeVector(conv, 'out_MST/mcp_mintreeall', overwrite = T)
        }

        # legend
        if(caption == TRUE){
          x <- 'All species'
          legend(x = "bottomleft", legend = x, pch = 19, col = cols1[unique(rownames(tempo))[1]],
                 bty = 'o', pt.cex = 0.6, cex = 0.6)
        }

        hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo, context_label = "mintreeall all taxa")

        # presence-absence matrix:
        # ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                                                                                                  unique(rownames(tempo))))
        linhasRaster <- build_mst_presence_raster(hit_cells)
        writeRaster(linhasRaster, "out_MST/presence_mintreeall.tif",
                    overwrite = TRUE)
        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1], alpha = transp)

        coor.l[valid_cells, ] <- 0
        if (length(hit_cells) > 0) {
          coor.l[hit_cells, ] <- 1
        }


        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

        return(coor.ll)


      } else if (mintreeall == FALSE){

        conta <- 0
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))),
                         dimnames = list(seq(1:ncellras), unique(rownames(tempo)))
                         [c(1:2)])
        lista_r <- list()
        track_numbers <- setNames(rep(NA, length(unique(rownames(tempo)))), unique(rownames(tempo)))
        for(j in unique(rownames(tempo))){
          conta <- conta + 1
          tempoo <- subset(tempo[, 1:2], rownames(tempo) == j)

          if (nrow(tempoo) < 2 || nrow(unique(tempoo[, 1:2, drop = FALSE])) < 2) {
            warning(paste0("Skipping taxon '", j, "': fewer than two distinct occurrence points for MST."))
            coor.l[, conta] <- 0
            next
          }

          #### shapefile ###
          tempo.d <- as.data.frame(tempoo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          write.shapefile(tempo_shape, paste0(c('out_MST/pointshape_', j), collapse = ''))
          resul1_shape <- vect(paste0(c('out_MST/pointshape_', j, '.shp'),
                                      collapse = ''), crs = "+proj=longlat +datum=WGS84")
          # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
          # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          terra::project(resul1_shape, "+proj=longlat +datum=WGS84")
          # resul1_shape <- vect(resul1_shape, crs = "+proj=longlat")

          ##### MST based on geographic distance #####
          # rownames(resul1_shape@coords) <- rownames(tempo.d)
          # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
          mst2 <- dino.mst(dista)
          mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = paste0("taxon '", j, "'"))
          lats <- cbind(resul1_shape$long,
                        resul1_shape$lat)
          rownames(lats) <- rownames(tempo.d)
          colnames(lats) <- c('longitude', 'latitude')
          mst_shape <- msn2Shape(msn = mst2, lats = lats)
          write.shapefile(mst_shape, paste0(c('out_MST/mst_', j), collapse = ''))

          if(!is.null(tree)){
            print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
          } else {
            print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
          }

          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out_MST/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, add = T)
          plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)

          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempoo)
            write.shapefile(conv, paste0(c('out_MST/mcp_', j), collapse = ''))
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
          }

          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- tempo[which(rownames(tempo) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- tempo[which(rownames(tempo) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))


          hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo.d[, 1:2, drop = FALSE], context_label = paste0("taxon '", j, "'"))

          lista_r[[conta]] <- build_mst_presence_raster(hit_cells)
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)

          # teste <- tabelao[j,j] # posicoes dos taxons do no
          teste <- conta
          track_numbers[j] <- teste  # where 'teste' is the calculated track number
          text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])

          coor.l[valid_cells, conta] <- 0
          if (length(hit_cells) > 0) {
            coor.l[hit_cells, conta] <- 1
          }
        }

        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo))
          # Create legends with track numbers
          legend_labels <- paste0(x, " (", track_numbers[x], ")")
          legend(x = "bottomleft", legend = legend_labels, pch = 19, col = cols1[x],
                 bty = 'o', pt.cex = 0.6, cex = 0.6)
        }

        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

        return(coor.ll)


      }
    }
    ########### Choosing a node or nodes... ############
  } else if(!is.null(nodes)){

    anc <- nodes
    no_num <- 0

    ## looping the nodes:
    for(k in anc){
      no_num <- no_num + 1

      lis <- arv$tip.label[arv$edge[which(arv$edge[,1] == k), 2]]

      # criando a condicao de checagem dos taxons terminais:
      if (any(is.na(lis)) == FALSE){
        for(j in lis){
          # sppp <- which(unique(rownames(coords)) == taxon[j])

          ##### If we want to make MST of all specimens TOGETHER from the chosen species...
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          tabela[(qde[idx - 1] + 1) : qde[idx], 4] <- rep(no_num, (qde[idx] - qde[idx - 1]))
        }

      } else if(any(is.na(lis)) == TRUE){

        # tabelao <- mrca(phy = tree, full = F)
        teste <- which(tabelao %in% anc[no_num]) # posicoes dos taxons do no
        tabelao.v <- as.vector(tabelao)
        names(tabelao.v) <- rep(colnames(tabelao), dim(tabelao)[1])
        clado <- unique(names(tabelao.v)[teste])

        for(j in clado){
          # sppp <- which(unique(rownames(coords)) == taxon[j])

          ##### If we want to make MST of all specimens TOGETHER from the chosen species...
          idx <- which(names(qde) == j)
          tabela[(qde[idx - 1] + 1) : qde[idx], 2] <- coordin[which(rownames(coordin) == j), 1]
          tabela[(qde[idx - 1] + 1) : qde[idx], 3] <- coordin[which(rownames(coordin) == j), 2]
          tabela[(qde[idx - 1] + 1) : qde[idx], 1] <- rep(j, (qde[idx] - qde[idx - 1]))
          tabela[(qde[idx - 1] + 1) : qde[idx], 4] <- rep(no_num, (qde[idx] - qde[idx - 1]))
        }
      }
    }

    # creating a table just for those terminal species
    if(!is.null(taxon)){
      # preparing species' msts...
      #######
      ### vector for accumulating specimens ###
      qde_sp <- 0
      for(b in species){qde_sp[b] <- length(which(rownames(coordin) == b))}
      qde_sp <- cumsum(qde_sp)
      tab_sp <- matrix(NA, nrow = dim(coordin), nc = 4)
      no_num <- 0
      #######

      for(z in taxon){
        no_num <- no_num + 1
        # sppp <- which(unique(rownames(coordin)) == taxon[j])

        ##### Se quisermos fazer MST de todos os especimes JUNTOS das especies escolhidas...
        idx <- which(names(qde_sp) == z)
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 2] <- coordin[which(rownames(coordin) == z), 1]
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 3] <- coordin[which(rownames(coordin) == z), 2]
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 1] <- rep(z, (qde_sp[idx] - qde_sp[idx - 1]))
        tab_sp[(qde_sp[idx - 1] + 1) : qde_sp[idx], 4] <- rep(no_num, (qde_sp[idx] - qde_sp[idx - 1]))
      }
      tabela <- tabela
    } else if(is.null(taxon)){
      tabela <- tabela
    }

    if(mintreeall == TRUE){

      # internal nodes:
      print('adding internal nodes')
      frame.n <- as.data.frame(na.exclude(tabela))
      colnames(frame.n) <- c('species', colnames(coordin[, c(1, 2)]))

      Long_n <- frame.n[,2]
      Lat_n <- frame.n[, 3]
      names(Long_n) <- frame.n[,1]

      tempo.temp2 <- matrix(cbind(as.numeric(Long_n), as.numeric(Lat_n)), nrow(frame.n), 2,
                            dimnames = list(frame.n[,1], colnames(frame.n)[c(2,3)]))
      #######

      if(!is.null(taxon)){

        # terminal nodes:
        print('adding terminal nodes...')
        frame <- as.data.frame(na.exclude(tab_sp))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))

        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]

        tempo.temp <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2,
                             dimnames = list(frame[,1], colnames(frame)[c(2,3)]))
        #######

        # putting in there alltogether:
        total <- rbind(frame.n, frame)

        Long <- total[,2]
        Lat <- total[, 3]
        names(Long) <- total[,1]

        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(total), 2,
                        dimnames = list(total[,1], colnames(total)[c(2,3)]))

        #### shapefiles ###
        ## terminal + nodes
        tempo.d <- as.data.frame(tempo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        # dir.create('temp/')
        write.shapefile(tempo_shape, 'out_MST/ancterminal_points_mintreeall')
        resul1_shape <- vect('out_MST/ancterminal_points_mintreeall.shp',
                             crs = "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        terra::project(resul1_shape, "+proj=longlat +datum=WGS84")

        #only terminals:
        tempo_temp.d <- as.data.frame(tempo.temp)
        tempo_shape.temp <- lats2Shape(lats = tempo_temp.d)
        write.shapefile(tempo_shape.temp, 'out_MST/points_mintreeall_onlyterminal')
        resul1_shape.temp <- vect('out_MST/points_mintreeall_onlyterminal.shp',
                                  crs = "+proj=longlat +datum=WGS84")
        terra::project(resul1_shape.temp, "+proj=longlat +datum=WGS84")
        # proj4string(resul1_shape.temp) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        # projection(resul1_shape.temp) <- CRS("+proj=longlat +datum=WGS84")
      }

      #only internal nodes:
      tempo_temp.d2 <- as.data.frame(tempo.temp2)
      tempo_shape.temp <- lats2Shape(lats = tempo_temp.d2)
      write.shapefile(tempo_shape.temp, 'out_MST/points_mintreeall_onlyinternal')
      resul2_shape.temp <- vect('out_MST/points_mintreeall_onlyinternal.shp',
                                crs = "+proj=longlat +datum=WGS84")
      terra::project(resul2_shape.temp, "+proj=longlat +datum=WGS84")
      # proj4string(resul2_shape.temp) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
      # projection(resul2_shape.temp) <- CRS("+proj=longlat +datum=WGS84")

      if(!is.null(taxon)){
        ##### total MST based on geographic distance #####
        # rownames(resul1_shape@coords) <- rownames(tempo.d)
        # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
        mst2 <- dino.mst(dista)
        mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = "internal+terminal mintreeall")
        lats <- cbind(resul1_shape$long,
                      resul1_shape$lat)
        rownames(lats) <- rownames(tempo.d)
        colnames(lats) <- c('longitude', 'latitude')
        mst_shape <- msn2Shape(msn = mst2, lats = lats)
        # mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
        write.shapefile(mst_shape, 'out_MST/mst_ancterminal_mintreeall')

        if(!is.null(tree)){
          print('5) calculating total mst... Done')
        } else {
          print('4) calculating total mst... Done')
        }

        # plotting shapes
        pontos_linha <- shapefile('out_MST/mst_ancterminal_mintreeall.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        crs(pontos_linha) <- "+proj=longlat +datum=WGS84"
        # terra::project(pontos_linha, "+proj=longlat +datum=WGS84")

        plot(pontos_linha, col = 'red', lwd = 3, lty = 2, add = T)
        plot(resul2_shape.temp, cex = 1.5, pch = 'o', add = T)
        # colocando labels:
        text(resul2_shape.temp, labels = rep(nodes, dim(tempo)[1]), cex = 0.8, pos = 1,
             col = cols1[lis])

        # plotting only terminal:
        plot(resul1_shape.temp, cex = 1.5, pch = 21, bg = cols1[taxon], add = T)
        for(nomes in taxon){
          teste <- tabelao[nomes,nomes] # posicoes dos taxons do no
          text(tempo_temp.d, labels = rep(teste, dim(tempo_temp.d)[1]), cex = 0.8, pos = 2, col = cols1[nomes])
        }


        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo)
          plot(conv, add = T, col = adjustcolor('gray', transp))
          conv <- st_as_sf(conv)
          conv <- as_Spatial(conv)
          conv <- vect(conv, crs = '+proj=longlat +datum=WGS84')
          writeVector(conv, 'out_MST/mcp_mintreeall_ancterminal', overwrite = T)
        }

        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo))
          xx <- unique(rownames(tempo.temp))
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 'o', col = cols1[lis], title = paste0(c('Adding descendents from node ',
                                                                   nodes, '...'), collapse = ''), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(as(shape_file, 'Spatial'))[1], y = extent(as(shape_file, 'Spatial'))[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...',
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                     pch = 'o', col = cols1[lis], title = paste0(c('Descendents from node ', nodes),
                                                                 collapse = ''), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(as(shape_file, 'Spatial'))[1], y = extent(as(shape_file, 'Spatial'))[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...',
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 'o', col = cols1[lis], title = paste0(c('Adding descendents from nodes',
                                                                   nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(as(shape_file, 'Spatial'))[1], y = extent(as(shape_file, 'Spatial'))[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...',
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                     pch = 'o', col = cols1[lis], title = paste0(c('Descendents from nodes', nodes),
                                                                 collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = extent(as(shape_file, 'Spatial'))[1], y = extent(as(shape_file, 'Spatial'))[4], legend = xx,
                     pch = 'o', col = cols1[xx], title = 'Adding terminal(s)...',
                     title.col = 'blue', pt.cex = 0.6, cex = 0.6)
            }
          }
        }

        hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo, context_label = "internal+terminal mintreeall")

        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                                                                                                  unique(rownames(tempo))))  # tabela com todas as celulas
        linhasRaster <- build_mst_presence_raster(hit_cells)
        writeRaster(linhasRaster, "out_MST/presence_mst_mintreeall_ancterminal.tif",
                    overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)

        coor.l[valid_cells, ] <- 0
        if (length(hit_cells) > 0) {
          coor.l[hit_cells, ] <- 1
        }

      } else {
        ##### MST of the internal nodes based on geographic distance #####
        # rownames(resul2_shape.temp@coords) <- rownames(tempo_temp.d2)
        # colnames(resul2_shape.temp@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = as(resul2_shape.temp, 'Spatial'))
        mst2 <- dino.mst(dista)
        mst2 <- validate_mst_matrix(mst2, rownames(tempo_temp.d2), context_label = "internal nodes mintreeall")
        lats <- cbind(resul2_shape.temp$long,
                      resul2_shape.temp$lat)
        rownames(lats) <- rownames(tempo_temp.d2)
        colnames(lats) <- c('longitude', 'latitude')
        mst_shape <- msn2Shape(msn = mst2, lats = lats)
        write.shapefile(mst_shape, 'out_MST/mst_mintreeall_onlyinternal')

        if(!is.null(tree)){
          print('5) calculating mst from internal node(s)... Done')
        } else {
          print('4) calculating mst from internal node(s)... Done')
        }

        # plotting shapes
        pontos_linha <- shapefile('out_MST/mst_mintreeall_onlyinternal.shp', warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        crs(pontos_linha) <- "+proj=longlat +datum=WGS84"
        # terra::project(pontos_linha, "+proj=longlat +datum=WGS84")

        plot(pontos_linha, col = 'red', lwd = 2, lty = 2, add = T)
        plot(resul2_shape.temp, cex = 1.5, pch = 'o', add = T)
        # colocando labels:
        text(resul2_shape.temp, labels = rep(nodes, dim(tempo_temp.d2)[1]), cex = 0.8, pos = 1,
             col = cols1[lis])


        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempo_temp.d2)
          plot(conv, add = T, col = adjustcolor('gray', transp))
          conv <- st_as_sf(conv)
          conv <- as_Spatial(conv)
          conv <- vect(conv, crs = '+proj=longlat +datum=WGS84')
          writeVector(conv, 'out_MST/mcp_mintreeall_ancterminal', overwrite = T)
        }


        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.temp2))
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 'o', title = paste0(c('Adding descendents from node',
                                                 nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                     pch = 'o', title = paste0(c('Descendents from node', nodes),
                                               collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 'o', title = paste0(c('Adding descendents from nodes',
                                                 nodes, '...'), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x,
                     pch = 'o', title = paste0(c('Descendents from nodes', nodes),
                                               collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)

            }
          }
        }

        hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo_temp.d2[, 1:2, drop = FALSE], context_label = "internal nodes mintreeall")

        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo.temp2))), dimnames = list(seq(1:ncellras),
                                                                                                        unique(rownames(tempo.temp2)))[c(1:2)])  # tabela com todas as celulas
        linhasRaster <- build_mst_presence_raster(hit_cells)
        writeRaster(linhasRaster, "out_MST/presence_mst_mintreeall_ancterminal.tif",
                    overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)

        coor.l[valid_cells, ] <- 0
        if (length(hit_cells) > 0) {
          coor.l[hit_cells, ] <- 1
        }
      }

      coor.l <- na.exclude(coor.l)
      coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
      rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
      write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

      return(coor.ll)


    } else if(mintreeall == FALSE){
      # internal nodes:
      print('adding internal nodes...')
      frame.n <- as.data.frame(na.exclude(tabela))
      # frame.n <- as.data.frame(tabela)
      colnames(frame.n) <- c('species', colnames(coordin[, c(1, 2)]), 'nodes')

      Long_n <- frame.n[,2]
      Lat_n <- frame.n[, 3]
      nodd <- frame.n[, 4]
      names(Long_n) <- frame.n[,1]
      #######

      # terminal nodes:
      if(!is.null(taxon)){
        print('adding terminal nodes...')
        frame <- as.data.frame(na.exclude(tab_sp))
        colnames(frame) <- c('species', colnames(coordin[, c(1, 2)]))

        Long <- frame[,2]
        Lat <- frame[, 3]
        names(Long) <- frame[,1]
      }
      #######

      # internal nodes:
      tempo.n <- matrix(cbind(as.numeric(Long_n), as.numeric(Lat_n), as.numeric(nodd)), nrow(frame.n), 3,
                        dimnames = list(frame.n[,1], colnames(frame.n)[c(2:4)]))

      if(!is.null(taxon)){
        #terminal nodes:
        tempo <- matrix(cbind(as.numeric(Long), as.numeric(Lat)), nrow(frame), 2,
                        dimnames = list(frame[,1], colnames(frame)[c(2,3)]))
      }

      if(!is.null(taxon)){
        #let's put in there alltogether!
        coor.l <- matrix(NA, nr = ncellras, nc = sum(length(unique(rownames(tempo.n))), length(unique(rownames(tempo)))),
                         dimnames = list(seq(1:ncellras), c(unique(rownames(tempo.n)), unique(rownames(tempo))))[c(1:2)])  # tabela com todas as celulas
      } else {
        #let's put in there alltogether!
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo.n))),
                         dimnames = list(seq(1:ncellras), unique(rownames(tempo.n)))[c(1:2)])  # tabela com todas as celulas

      }

      # first, to internal nodes...
      lista_r <- list()
      conta <- 0

      for(j in unique(rownames(tempo.n))){
        conta <- conta + 1
        tempoo <- subset(tempo.n[, 1:2], rownames(tempo.n) == j)

        #### shapefile ###
        tempo.d <- as.data.frame(tempoo)
        tempo_shape <- lats2Shape(lats = tempo.d)
        # dir.create('out_MST/')
        write.shapefile(tempo_shape, paste0(c('out_MST/pointshape_', j), collapse = ''))
        resul1_shape <- vect(paste0(c('out_MST/pointshape_', j, '.shp'),
                                    collapse = ''),
                             crs = "+proj=longlat +datum=WGS84")
        # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
        # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        terra::project(resul1_shape, "+proj=longlat +datum=WGS84")

        ##### MST based on geographic distance #####
        # rownames(resul1_shape@coords) <- rownames(tempo.d)
        # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
        dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
        mst2 <- dino.mst(dista)
        mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = paste0("node-descendent taxon '", j, "'"))
        lats <- cbind(resul1_shape$long, resul1_shape$lat)
        rownames(lats) <- rownames(tempo.d)
        colnames(lats) <- c('longitude', 'latitude')
        mst_shape <- msn2Shape(msn = mst2, lats = lats)
        write.shapefile(mst_shape, paste0(c('out_MST/mst_', j), collapse = ''))

        if(!is.null(tree)){
          print(paste0(c(conta + 5, ') calculating mst... Done'), collapse = ''))
          contass <- conta + 5
        } else {
          print(paste0(c(conta + 4, ') calculating mst... Done'), collapse = ''))
          contass <- conta + 4
        }

        # plotting shapes
        pontos_linha <- shapefile(paste0(c('out_MST/mst_', j, '.shp'), collapse = ''),
                                  warnPRJ = FALSE)
        proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
        crs(pontos_linha) <- "+proj=longlat +datum=WGS84"
        plot(pontos_linha, col = 'red', lwd = 3, lty = 2, add = T)
        plot(resul1_shape, cex = 1.1, pch = 21, bg = cols1[j], add = T)

        # labels:
        text(resul1_shape, labels = rep(nodes[unique(tempo.n[which(rownames(tempo.n) == j), 3])],
                                        dim(tempoo)[1]), cex = 0.6, pos = 1, col = 'black')

        # minimum convex polygon
        if(pol == TRUE){
          conv <- convexhull.xy(tempoo)
          conv <- st_as_sf(conv)
          conv <- as_Spatial(conv)
          conv <- vect(conv, crs = "+proj=longlat +datum=WGS84")
          plot(conv, add = T, col = adjustcolor(cols1[j], transp))
          writeVector(conv, paste0(c('out_MST/mcp_', j), collapse = ''), overwrite = T)
        }

        hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo.d[, 1:2, drop = FALSE], context_label = paste0("node-descendent taxon '", j, "'"))

        lista_r[[conta]] <- build_mst_presence_raster(hit_cells)
        writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                             collapse = ''), overwrite = TRUE)

        plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE,
             col = cols1[unique(tempo.n[which(rownames(tempo.n) == j), 3])], alpha = transp)

        coor.l[valid_cells, conta] <- 0
        if (length(hit_cells) > 0) {
          coor.l[hit_cells, conta] <- 1
        }
      }


      if(!is.null(taxon)){
        #finally, the terminal nodes:
        contas <- 0
        for(j in taxon){
          contas <- contas + 1
          conta <- conta + 1
          tempoo <- subset(tempo[, 1:2], rownames(tempo) == j)

          #### shapefile ###
          tempo.d <- as.data.frame(tempoo)
          tempo_shape <- lats2Shape(lats = tempo.d)
          # dir.create('temp/')
          write.shapefile(tempo_shape, paste0(c('out_MST/pointshape_', j), collapse = ''))
          resul1_shape <- vect(paste0(c('out_MST/pointshape_', j, '.shp'),
                                      collapse = ''), crs = "+proj=longlat +datum=WGS84")
          # resul1_shape <- readShapeSpatial('tempshape1_out_mst.shp') # shapefile
          # proj4string(resul1_shape) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          terra::project(resul1_shape, "+proj=longlat +datum=WGS84")

          ##### MST based on geographic distance #####
          # rownames(resul1_shape@coords) <- rownames(tempo.d)
          # colnames(resul1_shape@coords) <- c('longitude', 'latitude')
          dista <- earth.dist(lats = as(resul1_shape, 'Spatial'))
          mst2 <- dino.mst(dista)
          mst2 <- validate_mst_matrix(mst2, rownames(tempo.d), context_label = paste0("terminal taxon '", j, "'"))
          lats <- cbind(resul1_shape$long, resul1_shape$lat)
          rownames(lats) <- rownames(tempo.d)
          colnames(lats) <- c('longitude', 'latitude')
          mst_shape <- msn2Shape(msn = mst2, lats = lats)
          # mst_shape <- msn2Shape(msn = mst2, lats = resul1_shape, dist = NULL)
          write.shapefile(mst_shape, paste0(c('out_MST/mst_', j), collapse = ''))

          print(paste0(c(contas + contass, ') calculating mst... Done'), collapse = ''))

          # plotting shapes
          pontos_linha <- shapefile(paste0(c('out_MST/mst_', j, '.shp'), collapse = ''),
                                    warnPRJ = FALSE)
          proj4string(pontos_linha) <- CRS("+proj=longlat +datum=WGS84") # wgs84 datum
          plot(pontos_linha, col = cols1[j], lwd = 3, lty = 2, add = T)
          # plotting only terminal:
          plot(resul1_shape, cex = 1.5, pch = 21, bg = cols1[j], add = T)

          teste <- tabelao[j, j] # posicoes dos taxons do no
          text(tempoo, labels = rep(teste, dim(tempoo)[1]), cex = 0.8, pos = 2, col = cols1[j])

          # minimum convex polygon
          if(pol == TRUE){
            conv <- convexhull.xy(tempoo)
            plot(conv, add = T, col = adjustcolor(cols1[j], transp))
            conv <- st_as_sf(conv)
            conv <- as_Spatial(conv)
            conv <- vect(conv, crs = "+proj=longlat +datum=WGS84")
            writeVector(conv, paste0(c('out_MST/mcp_', j), collapse = ''), overwrite = T)
          }

          hit_cells <- mst_hit_cells(pontos_linha, point_coords = tempo.d[, 1:2, drop = FALSE], context_label = paste0("terminal taxon '", j, "'"))

          lista_r[[conta]] <- build_mst_presence_raster(hit_cells)
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)

          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)

          coor.l[valid_cells, conta] <- 0
          if (length(hit_cells) > 0) {
            coor.l[hit_cells, conta] <- 1
          }
        }
      }

      if(!is.null(taxon)){
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.n))
          xx <- unique(rownames(tempo))
          xxx <- unique(tempo.n[,3])
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the node',
                                                                nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                     pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the node',
                                                                                              nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                     pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                     title.col = 'blue', pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the nodes',
                                                                nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                     pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                     title.col = 'blue', pt.cex = 1.2, cex = 0.6, bty = 'n')
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the nodes',
                                                                                              nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'bottomright', legend = xx,
                     pch = 19, col = cols1[xx], title = 'Adding terminal node(s)...',
                     title.col = 'blue', pt.cex = 1.2, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)
            }
          }
        }
      } else {
        # legend
        if(caption == TRUE){
          x <- unique(rownames(tempo.n))
          xxx <- unique(tempo.n[,3])
          if(length(anc) == 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the node',
                                                                nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the node',
                                                                                              nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            }
          } else if(length(anc) > 1){
            if(sobrepo == TRUE){
              legend(x = extent(as(shape_file, 'Spatial'))[2] - 30, y = extent(as(shape_file, 'Spatial'))[4], legend = x,
                     pch = 19, col = cols1[x], title = paste0(c('Adding descendents from the nodes',
                                                                nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            } else if(sobrepo == FALSE){
              legend(x = "bottomleft", legend = x, pch = 19, col = cols1[x], title = paste0(c('Descendents from the nodes',
                                                                                              nodes), collapse = ' '), title.col = 'red', pt.cex = 0.6, cex = 0.6)
              legend(x = 'topright', legend = nodes, pch = 15, col = cols1[xxx], title = 'Adding the node(s)...', title.col = 'blue',
                     pt.cex = 1.2, cex = 0.6)

            }
          }
        }

      }

      coor.l <- na.exclude(coor.l)
      coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
      rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
      write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

      return(coor.ll)

    }
  }
}
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
  # producing a raster of the shapefile using polygon coverage so border
  # cells intersecting the study area remain valid for buffer occupancy
  mask.raster <- raster(extent(shape_spatial), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  cover.r <- rasterize(shape_spatial, mask.raster, getCover = TRUE)
  cover.r <- merge(cover.r, mask.raster)
  cover.r[is.na(cover.r) | cover.r <= 0] <- NA
  cover.r[!is.na(cover.r)] <- 1
  r <- cover.r
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
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
  # preparing the raster with grid of chosen shapefile
  grid <- raster(extent(shape_spatial), resolution = resol,
                 crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  # suppressWarnings(proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84")) # datum WGS84
  #proj4string(gridPolygon) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84"
  #########
  # producing a raster of the shapefile using polygon coverage so the valid grid
  # matches the study area boundary instead of the full bounding box
  mask.raster <- raster(extent(shape_spatial), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  cover.r <- rasterize(shape_spatial, mask.raster, getCover = TRUE)
  cover.r <- merge(cover.r, mask.raster)
  cover.r[is.na(cover.r) | cover.r <= 0] <- NA
  cover.r[!is.na(cover.r)] <- 1
  r <- cover.r
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  ncellras <- ncell(r)
  grid_index_r <- raster::setValues(raster::raster(r), seq_len(ncellras))
  raster::values(grid_index_r)[is.na(raster::values(r))] <- NA
  grid_all_sf <- sf::st_as_sf(rasterToPolygons(grid_index_r, dissolve = FALSE))
  names(grid_all_sf)[names(grid_all_sf) != attr(grid_all_sf, "sf_column")] <- "grid_id"
  grid_all_sf$grid_id <- as.integer(grid_all_sf$grid_id)
  grid_all_sf <- grid_all_sf[, c("grid_id", attr(grid_all_sf, "sf_column"))]

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

  ##################
  print(paste0(conta + 3, ') Loading the raster files and adding the presences...'))
  current_run_tifs <- unique(current_run_tifs[file.exists(current_run_tifs)])
  if (length(current_run_tifs) == 0) {
    stop("No convex-hull raster outputs were generated in the current run.")
  }
  layers_MCP <- terra::rast(current_run_tifs)
  # Sum all layers to get species richness per cell
  MCP_soma <- sum(layers_MCP, na.rm = TRUE)
  return(list(geometry = MCP_soma, pres_abs = coor.ll, grid_all_sf = grid_all_sf))
}
#' calcRange_irregularBins
#'
#' Extrapolates species distributions to irregular spatial units (biogeographic regions,
#' biomes, provinces, ecoregions, etc.) using spatial intersection. Creates a
#' presence-absence matrix and calculates species richness per bin.
#'
#' @param xy Data.frame or list containing species occurrence data with columns:
#'   - `spp`: species names
#'   - `Long` or `long`: longitude coordinates
#'   - `Lat` or `lat`: latitude coordinates
#' @param bins_shapefile Shapefile containing irregular bins. Accepts any of:
#'   - File path (character): "path/to/shapefile.shp"
#'   - sf object: from st_read()
#'   - terra::SpatVector: from vect()
#'   - sp::Spatial*: from readOGR() or shapefile()
#'   Must have a column identifying each bin (specified in `bin_id_column`).
#' @param bin_id_column Character. Name of the column in `bins_shapefile` that
#'   identifies each bin (e.g., "Dominio", "Biome", "Province"). Default: "bin_id"
#' @param resol Numeric vector of length 2. Resolution of grid cells in degrees
#'   (longitude, latitude) for generating grid shapefiles. Default: c(1, 1) (1 degrees x 1 degrees).
#'   Examples: c(0.5, 0.5) for 0.5 degrees cells, c(2, 2) for 2 degrees cells.
#' @param crs_input Integer or character. CRS of input coordinates. Default: 4326 (WGS84)
#' @param output_dir Character. Directory for output files. Default: "out_irregular_bins/"
#'
#' @return A list with four elements:
#'   - `bins_richness`: sf object with species richness per bin
#'   - `pres_abs`: presence-absence matrix (bins x species) with NUMERIC rownames
#'   - `species_per_bin`: data.frame with species count per bin
#'   - `bin_id_mapping`: data.frame mapping bin names to numeric IDs
#'
#' @details
#' This function performs the following steps:
#' 1. Converts occurrence points to sf object
#' 2. Validates and fixes invalid geometries in bins shapefile
#' 3. Transforms coordinates to match bins CRS
#' 4. Performs spatial join (st_join) to assign occurrences to bins
#' 5. Creates presence-absence matrix (bins x species)
#' 6. Calculates species richness per bin
#' 7. Saves outputs to specified directory
#'
#' **Input Validation:**
#' - Automatically fixes invalid geometries using `st_make_valid()`
#' - Handles different input formats (data.frame or list)
#' - Transforms CRS to match bins shapefile
#'
#' **Output Files:**
#' - `pres_abs_irregular_bins.txt`: presence-absence matrix (numeric rownames)
#' - `species_richness_bins.shp`: shapefile with richness per bin
#' - `species_per_bin.csv`: table with species count per bin
#' - `bin_id_mapping.csv`: mapping between bin names and numeric IDs
#'
#' @examples
#' \dontrun{
#' # Example with South American countries
#' library(sf)
#' library(rnaturalearthdata)
#'
#' # Load occurrence data
#' data <- data.frame(
#'   spp = c("sp1", "sp1", "sp2", "sp2", "sp3"),
#'   Long = c(-50, -51, -48, -49, -52),
#'   Lat = c(-10, -11, -12, -13, -14)
#' )
#'
#' # Get South America countries and Brazilian states as separate polygons
#' sa_countries <- ne_countries(continent = "South America", scale = "medium",
#'   returnclass = "sf")
#'
#' brazil_states <- ne_states(country = "Brazil", returnclass = "sf")
#'
#' # removing Brazil
#' sa_without_brazil <- sa_countries[sa_countries$iso_a3 != "BRA", ]
#'
#' # standardizing columns
#' countries <- sa_no_br[, c("name", "iso_a3")]
#' names(countries) <- c("region_name", "iso_code", "geometry")
#' countries$type <- "country"
#'
#' states <- brazil_states[, c("name", "iso_3166_2")]
#' names(states) <- c("region_name", "iso_code", "geometry")
#'
#' states$type <- "state"
#'
#' # validating geometries
#' countries_clean <- st_make_valid(countries)
#' states_clean <- st_make_valid(states)
#'
#' # Combining polygons
#' sa_bins <- rbind(countries_clean, states_clean)
#'
#' # Convert to terra SpatVector if needed
#' sa_bins_sv <- vect(sa_bins)
#'
#' # Calculate distribution using irregular bins
#' result <- calcRange_irregular_bins(
#'   xy = data,
#'   bins_shapefile = sa_bins,
#'   bin_id_column = "region_name",
#'   resol = c(1, 1)  # 1 degrees x 1 degrees grid cells
#' )
#'
#' # Access results
#' result$bins_richness        # sf object with richness
#' result$pres_abs             # presence-absence matrix (numeric rownames)
#' result$species_per_bin      # species count per bin
#' result$bin_id_mapping       # bin names <-> numeric IDs mapping
#'
#' # Plot richness map
#' library(ggplot2)
#' ggplot(result$bins_richness) +
#'   geom_sf(aes(fill = n_species)) +
#'   scale_fill_viridis_c() +
#'   theme_minimal()
#' }
#'
#' @export
calcRange_irregular_bins <- function(xy,
                                     bins_shapefile,
                                     bin_id_column = "bin_id",
                                     resol = c(1, 1),
                                     crs_input = 4326,
                                     output_dir = "out_irregular_bins/") {

  # Check for required packages
  required_packages <- c("sf", "dplyr", "lwgeom")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' not found. Please install: install.packages('", pkg, "')"),
           call. = FALSE)
    }
  }

  library(sf)
  library(dplyr)
  library(lwgeom)

  cat("============================================================\n")
  cat("  DISTRIBUTION EXTRAPOLATION USING IRREGULAR BINS\n")
  cat("============================================================\n\n")

  # ============================================================================
  # 1. PREPARE INPUT DATA
  # ============================================================================

  cat("1) Preparing occurrence data...\n")

  # Handle different input data types
  if (is.data.frame(xy)) {
    names(xy) <- tolower(names(xy))
    dat <- xy[, c("spp", "long", "lat")]
  } else if (is.list(xy) && "samples" %in% names(xy)) {
    dat <- xy$samples[, 1:3]
    names(dat) <- c("spp", "long", "lat")
  } else {
    stop("Input 'xy' must be a data.frame with columns 'spp', 'Long', 'Lat' or a list with 'samples'")
  }

  # Remove NA coordinates
  dat <- dat[complete.cases(dat), ]

  if (nrow(dat) == 0) {
    stop("No valid occurrence records found after removing NAs")
  }

  cat(paste0("   - Total occurrences: ", nrow(dat), "\n"))
  cat(paste0("   - Total species: ", length(unique(dat$spp)), "\n\n"))

  # ============================================================================
  # 2. CONVERT TO SF OBJECT
  # ============================================================================

  cat("2) Converting occurrences to sf object...\n")

  occurrences_sf <- st_as_sf(dat,
                             coords = c('long', 'lat'),
                             crs = crs_input)

  cat(paste0("   - CRS: ", st_crs(occurrences_sf)$input, "\n\n"))

  # ============================================================================
  # 3. LOAD AND VALIDATE BINS SHAPEFILE
  # ============================================================================

  cat("3) Loading and validating bins shapefile...\n")

  # Universal shapefile converter: accepts any type
  if (is.character(bins_shapefile)) {
    # File path provided
    bins_sf <- st_read(bins_shapefile, quiet = TRUE)
    cat(paste0("   - Loaded from file: ", bins_shapefile, "\n"))
  } else if (inherits(bins_shapefile, "sf")) {
    # Already sf object
    bins_sf <- bins_shapefile
    cat("   - Using provided sf object\n")
  } else if (inherits(bins_shapefile, "SpatVector")) {
    # terra::vect() object
    bins_sf <- st_as_sf(bins_shapefile)
    cat("   - Converted from terra::SpatVector to sf\n")
  } else if (inherits(bins_shapefile, "Spatial")) {
    # sp package object (SpatialPolygons*, etc.)
    bins_sf <- st_as_sf(bins_shapefile)
    cat("   - Converted from sp::Spatial* to sf\n")
  } else {
    stop(paste0("bins_shapefile type not recognized. ",
                "Accepted types: file path (character), sf, terra::SpatVector, or sp::Spatial*. ",
                "Provided type: ", class(bins_shapefile)[1]))
  }

  # Check if bin_id_column exists
  if (!bin_id_column %in% names(bins_sf)) {
    stop(paste0("Column '", bin_id_column, "' not found in bins_shapefile. ",
                "Available columns: ", paste(names(bins_sf), collapse = ", ")))
  }

  cat(paste0("   - Bin ID column: ", bin_id_column, "\n"))
  cat(paste0("   - Number of bins: ", nrow(bins_sf), "\n"))

  # Check for invalid geometries
  invalid_geoms <- !st_is_valid(bins_sf)

  if (any(invalid_geoms)) {
    cat(paste0("   - WARNING: Found ", sum(invalid_geoms), " invalid geometries\n"))
    cat("   - Fixing invalid geometries with st_make_valid()...\n")
    bins_sf <- st_make_valid(bins_sf)

    # Verify fix
    if (all(st_is_valid(bins_sf))) {
      cat("   - [OK] All geometries are now valid\n")
    } else {
      warning("Some geometries could not be fixed. Filtering them out...")
      bins_sf <- bins_sf[st_is_valid(bins_sf), ]
    }
  } else {
    cat("   - [OK] All geometries are valid\n")
  }

  cat("\n")

  # ============================================================================
  # 4. TRANSFORM CRS TO MATCH BINS
  # ============================================================================

  cat("4) Transforming coordinates to match bins CRS...\n")

  occurrences_transformed <- st_transform(occurrences_sf, st_crs(bins_sf))

  cat(paste0("   - Target CRS: ", st_crs(bins_sf)$input, "\n\n"))

  # ============================================================================
  # 5. SPATIAL JOIN (ASSIGN OCCURRENCES TO BINS)
  # ============================================================================

  cat("5) Performing spatial join (occurrences -> bins)...\n")

  intersections <- st_intersects(occurrences_transformed, bins_sf)
  n_assigned <- sum(lengths(intersections) > 0)
  n_outside <- nrow(occurrences_transformed) - n_assigned

  occurrences_with_bins <- st_join(occurrences_transformed, bins_sf)

  # Remove occurrences that didn't fall in any bin
  occurrences_with_bins <- occurrences_with_bins[!is.na(occurrences_with_bins[[bin_id_column]]), ]

  cat(paste0("   - Occurrences assigned to bins: ", n_assigned, "\n"))
  cat(paste0("   - Occurrences outside bins: ", n_outside, "\n\n"))

  if (nrow(occurrences_with_bins) == 0) {
    stop("No occurrences fell within any bin. Check CRS and spatial overlap.")
  }

  # ============================================================================
  # 6. CREATE PRESENCE-ABSENCE MATRIX
  # ============================================================================

  cat("6) Creating presence-absence matrix...\n")

  # Get unique species and bins
  all_species <- sort(unique(dat$spp))
  all_bins <- sort(unique(bins_sf[[bin_id_column]]))

  # Create mapping between bin names and numeric IDs
  bin_id_mapping <- data.frame(
    bin_name = all_bins,
    bin_number = 1:length(all_bins),
    stringsAsFactors = FALSE
  )

  # Initialize matrix with zeros (using AREA NAMES as rownames)
  pres_abs_matrix <- matrix(0,
                            nrow = length(all_bins),
                            ncol = length(all_species),
                            dimnames = list(all_bins, all_species))

  # Convert to data.frame for easier manipulation
  occ_df <- as.data.frame(occurrences_with_bins)
  occ_df <- occ_df[, c("spp", bin_id_column)]

  # Fill presence-absence matrix
  for (i in 1:nrow(occ_df)) {
    species <- as.character(occ_df$spp[i])
    bin_name <- as.character(occ_df[[bin_id_column]][i])

    if (bin_name %in% all_bins && species %in% all_species) {
      # Use bin_name directly as rowname
      pres_abs_matrix[bin_name, species] <- 1
    }
  }

  # Add ROOT row (for TNT compatibility)
  pres_abs_matrix <- rbind(pres_abs_matrix, ROOT = rep(0, ncol(pres_abs_matrix)))

  cat(paste0("   - Matrix dimensions: ", nrow(pres_abs_matrix) - 1, " bins x ",
             ncol(pres_abs_matrix), " species\n"))
  cat(paste0("   - Total presences: ", sum(pres_abs_matrix), "\n\n"))

  # ============================================================================
  # 7. CALCULATE SPECIES RICHNESS PER BIN
  # ============================================================================

  cat("7) Calculating species richness per bin...\n")

  # Count unique species per bin
  species_per_bin <- occurrences_with_bins %>%
    st_drop_geometry() %>%
    dplyr::group_by(!!rlang::sym(bin_id_column)) %>%
    dplyr::summarise(n_species = dplyr::n_distinct(spp),
                     n_occurrences = dplyr::n(),
                     species_list = paste(unique(spp), collapse = ", "),
                     .groups = "drop")

  cat(paste0("   - Bins with species: ", nrow(species_per_bin), "\n"))
  cat(paste0("   - Mean richness: ", round(mean(species_per_bin$n_species), 2), "\n"))
  cat(paste0("   - Max richness: ", max(species_per_bin$n_species), "\n\n"))

  # ============================================================================
  # 8. JOIN RICHNESS TO BINS SHAPEFILE
  # ============================================================================

  cat("8) Joining richness data to bins shapefile...\n")

  bins_with_richness <- bins_sf %>%
    left_join(species_per_bin, by = bin_id_column) %>%
    mutate(n_species = ifelse(is.na(n_species), 0, n_species),
           n_occurrences = ifelse(is.na(n_occurrences), 0, n_occurrences))

  cat("   - [OK] Richness data joined successfully\n\n")

  # ============================================================================
  # 9. SAVE OUTPUTS
  # ============================================================================

  cat("9) Saving outputs...\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(paste0("   - Created directory: ", output_dir, "\n"))
  }

  # Save presence-absence matrix
  write.table(pres_abs_matrix,
              file = file.path(output_dir, "pres_abs_irregular_bins.txt"),
              sep = "\t",
              quote = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "pres_abs_irregular_bins.txt"), "\n"))

  # Save species per bin table
  write.csv(species_per_bin,
            file = file.path(output_dir, "species_per_bin.csv"),
            row.names = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "species_per_bin.csv"), "\n"))

  # Save bin ID mapping (bin names <-> numeric IDs)
  write.csv(bin_id_mapping,
            file = file.path(output_dir, "bin_id_mapping.csv"),
            row.names = FALSE)
  cat(paste0("   - Saved: ", file.path(output_dir, "bin_id_mapping.csv"), "\n"))

  # Save shapefile with richness
  richness_shp <- file.path(output_dir, "species_richness_bins.shp")
  st_write(bins_with_richness,
           richness_shp,
           delete_layer = file.exists(richness_shp),
           quiet = TRUE)
  cat(paste0("   - Saved: ", richness_shp, "\n"))

  # Save individual species shapefiles
  cat("\n   - Saving individual species shapefiles...\n")

  for (species in all_species) {
    species_occ <- occurrences_with_bins[occurrences_with_bins$spp == species, ]

    if (nrow(species_occ) > 0) {
      occ_shp <- file.path(output_dir, paste0("pointshape_", species, ".shp"))
      st_write(species_occ,
               occ_shp,
               delete_layer = file.exists(occ_shp),
               quiet = TRUE)

      species_bins <- unique(species_occ[[bin_id_column]])
      species_bins_sf <- bins_sf[bins_sf[[bin_id_column]] %in% species_bins, ]

      bins_shp <- file.path(output_dir, paste0("BINS_", species, ".shp"))
      st_write(species_bins_sf,
               bins_shp,
               delete_layer = file.exists(bins_shp),
               quiet = TRUE)
    }
  }

  cat(paste0("   - Saved ", length(all_species), " species shapefiles\n"))
  cat("     * pointshape_[species].shp (occurrence points)\n")
  cat("     * BINS_[species].shp (bins where species occurs)\n\n")

  # Save grid shapefiles for each bin
  cat("   - Saving grid shapefiles for each bin...\n")

  # Check if raster package is available
  if (!requireNamespace("raster", quietly = TRUE)) {
    cat("     WARNING: Package 'raster' not found. Skipping grid shapefile generation.\n")
  } else {
    library(raster)

    # Create a grid based on bins extent with user-defined resolution
    bins_extent <- extent(as(bins_sf, "Spatial"))
    grid_raster <- raster(bins_extent, resolution = resol,
                          crs = CRS("+proj=longlat +datum=WGS84"))
    grid_polygons <- rasterToPolygons(grid_raster)
    grid_polygons_sf <- st_as_sf(grid_polygons)
    grid_polygons_sf$grid_id <- 1:nrow(grid_polygons_sf)

    # Spatial join: which grid cell belongs to which bin?
    grid_bin_join <- tryCatch({
      st_join(grid_polygons_sf, bins_sf[, c(bin_id_column, "geometry")],
              join = st_intersects, largest = TRUE)
    }, error = function(e) {
      cat(paste0("     WARNING: Grid-bin join failed with default geometry engine: ",
                 e$message, "\n"))
      cat("     - Retrying with validated geometries and s2 disabled...\n")

      tryCatch({
        old_s2 <- sf::sf_use_s2()
        on.exit(sf::sf_use_s2(old_s2), add = TRUE)
        sf::sf_use_s2(FALSE)

        bins_for_join <- bins_sf[, c(bin_id_column, "geometry")]
        bins_for_join <- tryCatch(st_make_valid(bins_for_join), error = function(e) bins_for_join)
        grid_for_join <- tryCatch(st_make_valid(grid_polygons_sf), error = function(e) grid_polygons_sf)

        st_join(grid_for_join, bins_for_join, join = st_intersects, largest = TRUE)
      }, error = function(e2) {
        warning(paste0("Could not assign grid cells to bins: ", e2$message,
                       ". Skipping grid shapefile generation."))
        NULL
      })
    })

    if (!is.null(grid_bin_join)) {
      # Save grid shapefile for each bin
      for (bin_name in all_bins) {
        bin_grids <- grid_bin_join[!is.na(grid_bin_join[[bin_id_column]]) &
                                     grid_bin_join[[bin_id_column]] == bin_name, ]

        if (nrow(bin_grids) > 0) {
          # Sanitize bin name for filename
          safe_bin_name <- gsub(" ", "_", bin_name)
          safe_bin_name <- gsub("[^A-Za-z0-9_]", "", safe_bin_name)

          grid_shp <- file.path(output_dir, paste0("GRIDS_", safe_bin_name, ".shp"))
          st_write(bin_grids,
                   grid_shp,
                   delete_layer = file.exists(grid_shp),
                   quiet = TRUE)
        }
      }

      cat(paste0("   - Saved grid shapefiles for ", length(all_bins), " bins\n"))
      cat("     * GRIDS_[bin_name].shp (grid cells within each bin)\n\n")
    } else {
      cat("   - WARNING: Grid shapefiles were skipped due to geometry issues\n\n")
    }
  }

  # ============================================================================
  # 10. SUMMARY
  # ============================================================================

  cat("============================================================\n")
  cat("  SUMMARY\n")
  cat("============================================================\n")
  cat(paste0("Total species: ", length(all_species), "\n"))
  cat(paste0("Total bins: ", length(all_bins), "\n"))
  cat(paste0("Bins with species: ", nrow(species_per_bin), "\n"))
  cat(paste0("Empty bins: ", length(all_bins) - nrow(species_per_bin), "\n"))
  cat(paste0("Mean richness per bin: ", round(mean(species_per_bin$n_species), 2), "\n"))
  cat(paste0("Max richness: ", max(species_per_bin$n_species),
             " (", species_per_bin[[bin_id_column]][which.max(species_per_bin$n_species)], ")\n"))
  cat("============================================================\n\n")

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  return(list(
    bins_richness = bins_with_richness,
    pres_abs = pres_abs_matrix,
    species_per_bin = species_per_bin,
    bin_id_mapping = bin_id_mapping
  ))
}
#' extrapolation_plot_from_shapeFiles
#'
#' Searches for, loads, and plots species distribution shapefiles from specified directories.
#' Can work with either a list of SpatVector objects or directly from directories containing
#' shapefiles (e.g., minimum convex polygons or buffers). Automatically assigns colors and
#' generates legends for multiple species.
#'
#' @param species_names Character vector with one or more species names to plot. If NULL,
#'   lists all available species in the directory without plotting
#' @param spatvector_list Optional list of SpatVector objects to search through. If provided,
#'   the function searches this list instead of loading from directories
#' @param directory Character string or vector specifying directory(ies) to search (default: "out_MCP").
#'   Can be "out_MCP", "out_buffers", or any custom directory path. If a vector is provided,
#'   searches all specified directories
#' @param file_type Character string specifying the file type pattern to match:
#'   "MCP" for minimum convex polygons, "BUFF" for buffer files, 'mst' for MSTs, or a custom pattern
#'   (default: "MCP")
#' @param add Logical, if TRUE adds to existing plot, if FALSE creates new plot (default: TRUE)
#' @param colors Named vector of colors for each species, or NULL for automatic colors using viridis palette
#' @param alpha Numeric value between 0 and 1 for transparency (default: 0.5)
#' @param legend Logical, if TRUE shows legend with species names (default: TRUE)
#' @param legend_position Character string specifying legend position: "topright", "topleft",
#'   "bottomright", "bottomleft", "top", "bottom", "left", "right", or "center" (default: "topright")
#' @param legend_inset Numeric value indicating the percentage of the space between the margin of the legend and
#'   the border of the plot
#' @param list_only Logical, if TRUE only lists available species without plotting (default: FALSE)
#' @param ... Additional graphical parameters passed to \code{\link[terra]{plot}}
#' @return If \code{list_only = TRUE}, returns a data frame with available species.
#'   Otherwise, invisibly returns a named list of loaded SpatVector objects
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function provides a flexible interface for working with species distribution
#' shapefiles. It can operate in three modes:
#'
#' \strong{Mode 1: List species} - When \code{species_names = NULL} or \code{list_only = TRUE},
#' lists all available species in the specified directory(ies).
#'
#' \strong{Mode 2: Search in list} - When \code{spatvector_list} is provided, searches for
#' species by name in the provided list of SpatVector objects (useful when you've already
#' loaded multiple shapefiles into a list).
#'
#' \strong{Mode 3: Load from directory} - When \code{spatvector_list = NULL}, searches for
#' and loads shapefiles directly from the default directory "out_MCP". Alternative directories are
#' 'out_buffers' for extrapolation done with defined area buffers and 'out' for MSTs.
#'
#' The function automatically handles multiple directories, allowing comparison of different
#' range estimation methods (e.g., MCP vs. buffers). Colors are automatically assigned using
#' the viridis color palette unless custom colors are provided. A legend is generated showing
#' all plotted species.
#'
#' This function is particularly useful for visualizing species distributions from PAE-PCE
#' analyses or other biogeographic studies, where multiple species need to be plotted together
#' with consistent styling. In this case, the file_type argument is 'mst'.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025.
#'   Panbiotracks: software for track analysis. \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'   Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1-17.
#' @seealso
#' \code{\link[terra]{vect}} for loading spatial vector data,
#' \code{\link[terra]{plot}} for plotting spatial data,
#' \code{\link[viridis]{viridis}} for color palettes.
#' @export
#' @examples
#' # Example 1: List available species in a directory (default: "out_MCP")
#' available <- extrapolation_plot_from_shapeFiles(list_only = TRUE)
#' print(available)
#'
#' # Example 2: Plot a single species from directory (default: "out_MCP")
#' # shape file of South America
#' lycipta.asul <- lycipta$polygon
#' plot(lycipta.asul, axes = TRUE, las = 1)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", legend_inset = 0.05) #  Legenda 5% para dentro
#'
#' # Example 3: Plot multiple species with automatic colors (default: "out_MCP")
#' plot(lycipta.asul, axes = TRUE, las = 1)
#' extrapolation_plot_from_shapeFiles(c("E_L_triangulator", "E_L_aceratos", "E_L_imitator")) # legenda cortando o mapa
#'
#' # Example 4: Plot with custom colors
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles(c("E_L_imitator", "E_L_triangulator"),
#'              colors = c("E_L_imitator" = "red", "E_L_triangulator" = "blue"))
#'
#' # Example 5: Plot from buffers directory
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", directory = "out_buffers", file_type = "BUFF")
#'
#' # Example 6: Plot from multiple directories (MCP + buffers)
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", directory = c("out_MCP", "out_buffers"))
#'
#' # Example 7: Search in a pre-loaded list
#' MCP_shapes_n <- list()
#' for(k in 1:length(species_names)){
#'   MCP_shapes_n[[k]] <- terra::vect(paste0('out_MCP/MCP_', species_names[k], '.shp'))
#' }
#' plot(lycipta.asul)
#' extrapolation_plot_from_shapeFiles("E_L_triangulator", spatvector_list = MCP_shapes_n)
#'
#' # Example 8: Integration with PAE-PCE results
#' pae1 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10), gridView = TRUE, labelGrid = TRUE,
#'    nonHomoplasticSpeciesList = FALSE, N = 10)
#' pae1 # data frame showing the homoplastic and non-homoplastic species (spp) and their grid numbers (grid_n),
#' # from the anterior iterations (PAE-PCE)
#'
#' especies_sinap <- unique(pae1$nonHomoplastic_species[[1]]$spp)
#' plot(lycipta.asul, axes = TRUE, las = 1, main = "PAE-PCE: Individual tracks of synapomorphic species")
#' extrapolation_plot_from_shapeFiles(especies_sinap, alpha = 0.6, directory = 'out', file_type = 'mst')
#'

extrapolation_plot_from_shapeFiles <- function(species_names = NULL,
                                               spatvector_list = NULL,
                                               directory = "out_MCP",
                                               file_type = "MCP",
                                               add = TRUE,
                                               colors = NULL,
                                               alpha = 0.5,
                                               legend = TRUE,
                                               legend_position = "topright",
                                               legend_inset = 0.02,
                                               list_only = FALSE,
                                               ...) {

  # ============================================================================
  # MODE 1: LIST SPECIES ONLY
  # ============================================================================

  if(is.null(species_names) || list_only) {
    all_species_df <- data.frame()

    for(dir in directory) {
      if(!dir.exists(dir)) {
        warning("Directory '", dir, "' does not exist. Skipping.")
        next
      }

      all_shp <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)

      if(length(all_shp) == 0) {
        message("No shapefiles found in directory '", dir, "'")
        next
      }

      pattern <- switch(file_type,
                        "MCP" = "MCP_",
                        "buffer" = "BUFF_",
                        file_type)

      filtered_shp <- grep(pattern, all_shp, value = TRUE)

      if(length(filtered_shp) == 0) {
        filtered_shp <- all_shp
      }

      basenames <- basename(filtered_shp)
      species_names_found <- gsub(paste0("^", pattern), "", basenames)
      species_names_found <- gsub("\\.shp$", "", species_names_found)

      df <- data.frame(
        Directory = dir,
        Species = species_names_found,
        Filename = basenames,
        FullPath = filtered_shp,
        stringsAsFactors = FALSE
      )

      all_species_df <- rbind(all_species_df, df)
    }

    if(nrow(all_species_df) > 0) {
      all_species_df$Index <- seq_len(nrow(all_species_df))
      all_species_df <- all_species_df[, c("Index", "Directory", "Species", "Filename", "FullPath")]
    }

    return(all_species_df)
  }

  # ============================================================================
  # MODE 2: SEARCH IN PROVIDED LIST
  # ============================================================================

  if(!is.null(spatvector_list)) {
    all_sources <- sapply(spatvector_list, terra::sources)
    loaded_vectors <- list()
    found_species <- character()

    for(sp in species_names) {
      # Use exact match: species name must be followed by .shp or end of string
      pattern_exact <- paste0(sp, "(\\.shp)?$")
      index <- grep(pattern_exact, all_sources)

      if(length(index) == 0) {
        warning("No SpatVector found with pattern: ", sp)
        next
      }

      if(length(index) > 1) {
        warning("Multiple matches found for '", sp, "'. Using first one.")
      }

      loaded_vectors[[sp]] <- spatvector_list[[index[1]]]
      found_species <- c(found_species, sp)
      cat("Found:", basename(all_sources[index[1]]), "\n")
    }

    if(length(loaded_vectors) == 0) {
      stop("No species could be found in the provided list")
    }

  } else {

    # ============================================================================
    # MODE 3: LOAD FROM DIRECTORY
    # ============================================================================

    loaded_vectors <- list()
    found_species <- character()

    for(dir in directory) {
      if(!dir.exists(dir)) {
        warning("Directory '", dir, "' does not exist. Skipping.")
        next
      }

      all_shp <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)

      if(length(all_shp) == 0) {
        warning("No shapefiles found in directory '", dir, "'")
        next
      }

      pattern <- switch(file_type,
                        "MCP" = "MCP_",
                        "buffer" = "BUFF_",
                        file_type)

      filtered_shp <- grep(pattern, all_shp, value = TRUE)

      if(length(filtered_shp) == 0) {
        warning("No files matching pattern '", pattern, "' found in '", dir, "'")
        filtered_shp <- all_shp
      }

      for(sp in species_names) {
        # Skip if already found
        if(sp %in% found_species) next

        # Use exact match: species name must be followed by .shp or end of string
        pattern_exact <- paste0(sp, "(\\.shp)?$")
        matching_files <- grep(pattern_exact, filtered_shp, value = TRUE)

        if(length(matching_files) == 0) {
          next
        }

        if(length(matching_files) > 1) {
          warning("Multiple files found for '", sp, "' in '", dir, "'. Using first one.")
        }

        tryCatch({
          v <- terra::vect(matching_files[1])
          loaded_vectors[[sp]] <- v
          found_species <- c(found_species, sp)
          cat("Loaded:", basename(matching_files[1]), "\n")
        }, error = function(e) {
          warning("Failed to load shapefile for '", sp, "': ", e$message)
        })
      }
    }

    if(length(loaded_vectors) == 0) {
      stop("No shapefiles could be loaded for the specified species")
    }

    # Warn about species not found
    not_found <- setdiff(species_names, found_species)
    if(length(not_found) > 0) {
      warning("Species not found: ", paste(not_found, collapse = ", "))
    }
  }

  # ============================================================================
  # PLOT ALL LOADED SPECIES
  # ============================================================================

  # Generate colors if not provided
  if(is.null(colors)) {
    colors <- setNames(viridis::viridis(length(found_species), option = 'D'),
                       found_species)
  } else if(length(colors) != length(found_species)) {
    colors <- setNames(viridis::viridis(length(found_species), option = 'D'),
                       found_species)
  } else if(is.null(names(colors))) {
    colors <- setNames(colors, found_species)
  }

  # Plot each species
  for(sp in found_species) {
    v <- loaded_vectors[[sp]]
    terra::plot(v, add = add, col = colors[sp], alpha = alpha,
                legend = FALSE, ...)
    add <- TRUE  # After first plot, always add
  }

  # Add legend
  if(legend && length(found_species) > 0) {
    par(xpd = TRUE)
    legend(legend_position,
           legend = found_species,
           fill = colors[found_species],
           title = "Species",
           bg = "white",
           cex = 0.7,
           inset = legend_inset,
           box.lwd = 1)
    par(xpd = FALSE)
  }

  invisible(loaded_vectors)
}
#' Wrapper for Range Extrapolation Functions
#'
#' Calls the appropriate range extrapolation function based on method selection
#' with proper error handling and fallback to mock data
#'
#' @param occurrence_data Data frame with columns: spp, long, lat
#' @param method Character string: "BUFF", "MPC", or "MST"
#' @param grid_resolution Numeric, grid resolution in degrees
#' @param buffer_width Numeric, buffer width in meters (for BUFF method)
#'
#' @return List with $pres_abs (matrix) and $geometry (sf object or NULL)
#'
#' @export
run_extrapolation <- function(occurrence_data, method = "BUFF", 
                              grid_resolution = 1, buffer_width = 100000) {
  
  tryCatch({
    # Ensure we have the required columns
    if (!all(c("spp", "long", "lat") %in% names(occurrence_data))) {
      stop("Data must have columns: spp, long, lat")
    }
    
    # Call the appropriate function based on method
    if (method == "BUFF") {
      result <- calcRange_buffers(
        xy = occurrence_data,
        buffer.width = buffer_width,
        shape_file = NULL,
        resol = grid_resolution,
        mean_dist = FALSE
      )
    } else if (method == "MPC") {
      result <- calcRange_convexHull(
        xy = occurrence_data,
        shape_file = NULL,
        resol = grid_resolution
      )
    } else if (method == "MST") {
      result <- calcRange_irregular_bins(
        xy = occurrence_data,
        bins_shapefile = NULL,
        resol = c(grid_resolution, grid_resolution)
      )
    } else {
      stop("Unknown method: ", method)
    }
    
    # Ensure result has the expected structure
    if (!is.list(result)) {
      stop("Function did not return a list")
    }
    
    # Extract pres_abs
    pres_abs <- if ("pres_abs" %in% names(result)) {
      result$pres_abs
    } else if ("pres.abs" %in% names(result)) {
      result$pres.abs
    } else if (is.matrix(result)) {
      result
    } else {
      stop("Could not extract presence-absence matrix from result")
    }
    
    # Extract geometry if available
    geometry <- if ("geometry" %in% names(result)) {
      result$geometry
    } else {
      NULL
    }
    
    # Try to convert geometry to sf if it's not NULL
    if (!is.null(geometry)) {
      tryCatch({
        geometry <- convert_geometry_to_sf(geometry)
      }, error = function(e) {
        warning("Could not convert geometry to sf: ", e$message)
        geometry <<- NULL
      })
    }
    
    # Return standardized result
    list(
      pres_abs = pres_abs,
      geometry = geometry,
      method = method,
      success = TRUE
    )
    
  }, error = function(e) {
    # Return error result
    warning("Extrapolation failed: ", e$message)
    list(
      pres_abs = NULL,
      geometry = NULL,
      method = method,
      success = FALSE,
      error = e$message
    )
  })
}

#' Create Mock Extrapolation Result
#'
#' Creates a mock result for testing when real functions fail
#'
#' @param occurrence_data Data frame with species occurrences
#' @param method Character string: "BUFF", "MPC", or "MST"
#'
#' @return List with $pres_abs (matrix) and $geometry (sf object with mock polygons)
#'
#' @export
create_mock_extrapolation <- function(occurrence_data, method = "BUFF") {
  
  tryCatch({
    # Create a simple presence-absence matrix
    species <- unique(occurrence_data$spp)
    n_species <- length(species)
    
    # Create a grid based on the extent of the data
    lon_range <- range(occurrence_data$long, na.rm = TRUE)
    lat_range <- range(occurrence_data$lat, na.rm = TRUE)
    
    # Create a simple grid
    lon_seq <- seq(lon_range[1] - 1, lon_range[2] + 1, by = 1)
    lat_seq <- seq(lat_range[1] - 1, lat_range[2] + 1, by = 1)
    
    n_cells <- length(lon_seq) * length(lat_seq)
    
    # Create random presence-absence matrix
    pres_abs <- matrix(
      sample(0:1, n_cells * n_species, replace = TRUE, prob = c(0.7, 0.3)),
      nrow = n_cells,
      ncol = n_species,
      dimnames = list(
        paste0("Cell_", 1:n_cells),
        species
      )
    )
    
    # Create mock geometry as sf polygons
    geometry <- NULL
    tryCatch({
      # Try to create sf polygons
      if (requireNamespace("sf", quietly = TRUE)) {
        polygons <- list()
        
        for (i in 1:(length(lon_seq) - 1)) {
          for (j in 1:(length(lat_seq) - 1)) {
            # Create a simple square polygon
            coords <- matrix(c(
              lon_seq[i], lat_seq[j],
              lon_seq[i+1], lat_seq[j],
              lon_seq[i+1], lat_seq[j+1],
              lon_seq[i], lat_seq[j+1],
              lon_seq[i], lat_seq[j]
            ), ncol = 2, byrow = TRUE)
            
            polygons[[length(polygons) + 1]] <- sf::st_polygon(list(coords))
          }
        }
        
        # Create sf object
        geometry <- sf::st_sfc(polygons, crs = 4326)
        geometry <- sf::st_sf(geometry = geometry)
      }
    }, error = function(e) {
      warning("Could not create sf geometry: ", e$message)
    })
    
    list(
      pres_abs = pres_abs,
      geometry = geometry,
      method = method,
      is_mock = TRUE
    )
    
  }, error = function(e) {
    warning("Could not create mock extrapolation: ", e$message)
    list(
      pres_abs = NULL,
      geometry = NULL,
      method = method,
      is_mock = TRUE,
      error = e$message
    )
  })
}

# FunÃ§Ãµes auxiliares para o BioGeoBEARS App

# FunÃ§Ã£o para formatar resultados
format_results <- function(results) {
  # ImplementaÃ§Ã£o futura
  return(results)
}

# FunÃ§Ã£o para calcular teste de razÃ£o de verossimilhanÃ§a
calc_LRT <- function(res1, res2) {
  # res1 deve ser o modelo mais simples (aninhado)
  # res2 deve ser o modelo mais complexo
  
  LnL1 <- res1$total_loglikelihood
  LnL2 <- res2$total_loglikelihood
  
  # DiferenÃ§a em nÃºmero de parÃ¢metros
  df <- length(res2$inputs$BioGeoBEARS_model_object@params_table[res2$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"]) - 
        length(res1$inputs$BioGeoBEARS_model_object@params_table[res1$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])
  
  # EstatÃ­stica de teste
  chisq <- 2 * (LnL2 - LnL1)
  
  # Valor p
  pval <- pchisq(chisq, df = df, lower.tail = FALSE)
  
  # Retornar resultados
  result <- list(
    model1 = res1$inputs$description,
    model2 = res2$inputs$description,
    LnL1 = LnL1,
    LnL2 = LnL2,
    df = df,
    chisq = chisq,
    pval = pval,
    significant = pval < 0.05
  )
  
  return(result)
}

# Outras funÃ§Ãµes auxiliares seriam adicionadas aqui

#' generate_tnt_pae_pce
#'
#' Creates a TNT script file for automated PAE-PCE (Parsimony Analysis of Endemicity
#' with Parsimony Constraint) analysis with multiple rounds of character removal.
#'
#' @param matrix_file Character. Path to the TNT matrix file (default: "out/pres_abs.tnt")
#' @param output_file Character. Name of the output TNT script file (default: "PAE_PCE_tntAUTO.run")
#' @param max_iterations Integer. Maximum number of PAE-PCE rounds to perform (default: 10)
#' @param search_replicates Integer. Number of search replicates to perform (default: 4)
#' @param search_method Character. Search method to use: "traditional" or "new_technology" (default: "new_technology")
#' @param output_dir Character. Directory for output files (default: "out_TNT/")
#'
#' @details
#' This function generates a TNT script that performs iterative PAE-PCE analysis:
#'
#' **Search Methods:**
#' - `"traditional"`: Uses traditional search with TBR branch swapping
#'   - Command: `mult 1000 = hold 100 tbr; bb=tbr;`
#' - `"new_technology"`: Uses new technology search with sectorial searches
#'   - Command: `xmult= hits 7 rep 5 rss level 10;`
#'
#' **Iterative Process:**
#' 1. Performs parsimony analysis
#' 2. Calculates CI (Consistency Index) and RI (Retention Index) for all characters
#' 3. Identifies non-homoplastic synapomorphies (CI=1, RI=1)
#' 4. Deactivates these characters cumulatively
#' 5. Repeats until no more synapomorphies are found or max_iterations is reached
#'
#' **Output Files (per round N):**
#' - `Round_OutputN.txt`: Complete log of the round
#' - `Round_CI_RI_TableN.txt`: CI/RI table for all characters
#' - `Round_SynapomorphiesN.txt`: List of non-homoplastic synapomorphies
#' - `Round_MPTN.tre`: Most parsimonious trees
#' - `Round_ConsensusN.tmp`: Consensus tree (TNT format)
#' - `Round_TreeN.tre`: Consensus tree (Newick format)
#' - `Round_MapN.svg`: Tree visualization with mapped synapomorphies
#'
#' **Summary Files:**
#' - `PAE_PCE_Summary.txt`: Summary of all rounds
#' - `resul_CI_RI.log`: CI/RI statistics for all rounds
#'
#' @return Invisibly returns the path to the generated script file
#'
#' @examples
#' \dontrun{
#' # Generate script with default settings (4 replicates, new technology)
#' generate_tnt_pae_pce()
#'
#' # Generate script with 20 replicates and traditional search
#' generate_tnt_pae_pce(
#'   search_replicates = 20,
#'   search_method = "traditional"
#' )
#'
#' # Generate script with custom settings
#' generate_tnt_pae_pce(
#'   matrix_file = "data/my_matrix.tnt",
#'   output_file = "my_pae_pce.run",
#'   max_iterations = 15,
#'   search_replicates = 10,
#'   search_method = "new_technology"
#' )
#' }
#'
#' @export
generate_tnt_pae_pce <- function(matrix_file = "out/pres_abs.tnt",
                                 output_file = "PAE_PCE_tntAUTO.run",
                                 max_iterations = 10,
                                 search_replicates = 4,
                                 search_method = "new_technology",
                                 output_dir = "out_TNT/") {

  ######################
  ## Input validation ##
  ######################

  if (!file.exists(matrix_file)) {
    stop(paste0("Error: Matrix file not found: ", matrix_file))
  }

  valid_methods <- c("traditional", "new_technology")
  if (!(search_method %in% valid_methods)) {
    stop(paste0("Error: search_method must be one of: ", paste(valid_methods, collapse = ", ")))
  }

  if (!is.numeric(max_iterations) || max_iterations < 1) {
    stop("Error: max_iterations must be a positive integer")
  }

  if (!is.numeric(search_replicates) || search_replicates < 1) {
    stop("Error: search_replicates must be a positive integer")
  }

  ######################
  ## Generate script ###
  ######################

  message("Generating TNT PAE-PCE script based on template...")
  message(paste0("  Matrix: ", matrix_file))
  message(paste0("  Output: ", output_file))
  message(paste0("  Max iterations: ", max_iterations))
  message(paste0("  Search replicates: ", search_replicates))
  message(paste0("  Search method: ", search_method))

  # Open connection to write script
  con <- file(output_file, "w", encoding = "UTF-8")

  # Write header (lines 1-56 from template)
  cat("macro= ;\n", file = con)
  cat("clb;\n", file = con)
  cat("mxram 1000;\n", file = con)
  cat("macro * 15 1000;\n", file = con)
  cat("macfloat 3;\n", file = con)
  cat("piwe-;\n", file = con)
  cat("\n", file = con)
  cat("report-;\n", file = con)
  cat("\n", file = con)
  cat(paste0("proc ", matrix_file, " ;\n"), file = con)
  cat("\n", file = con)
  cat("hold 99999;\n", file = con)
  cat("\n", file = con)
  cat(paste0("log ", output_dir, "PAE_PCE_Summary.txt ;\n"), file = con)
  cat("quote\n", file = con)
  cat("================================================================================\n", file = con)
  cat("           partially AUTOMATED PAE-PCE ANALYSIS - SUMMARY\n", file = con)
  cat("================================================================================\n", file = con)
  cat(paste0("Matrix: ", matrix_file, "\n"), file = con)
  cat(paste0("Max iterations: ", max_iterations, " (by author)\n"), file = con)
  cat("================================================================================\n", file = con)
  cat(";\n", file = con)
  cat("log/ ;\n", file = con)
  cat("\n", file = con)
  cat("\n", file = con)
  cat("\n", file = con)
  cat(paste0("log+ ", output_dir, "resul_CI_RI.log;\n"), file = con)
  cat("quote RUN, CI, RI, tree_length;\n", file = con)
  cat("log/;\n", file = con)
  cat("\n", file = con)
  cat("var: removed_count;\n", file = con)
  cat("\n", file = con)
  cat("/* ========================================================================== */\n", file = con)
  cat("/* MAIN ITERATION LOOP                                                        */\n", file = con)
  cat("/* ========================================================================== */\n", file = con)
  cat("\n", file = con)
  cat("/* Define variables for automation */\n", file = con)
  cat("var: total_removed conta contagem;\n", file = con)
  cat("var: theminchar themaxchar ci_val[(nchar + 1 )] ri_val[(nchar + 1)] rc_val[(nchar + 1)] char_i num_chars;\n", file = con)
  cat("var: THEMINtree THEMAXtree; \n", file = con)
  cat("var: THIStree CItree RItree RCtree nome list_idx;\n", file = con)
  cat("\n", file = con)
  cat("/* [OK] NEW: Global array to store ALL deactivated characters */\n", file = con)
  cat("var: all_deactivated[(nchar + 1)] total_deactivated;\n", file = con)
  cat("\n", file = con)
  cat("set num_chars nchar;\n", file = con)
  cat("set total_removed 0;\n", file = con)
  cat("set total_deactivated 0;\n", file = con)
  cat("\n", file = con)
  cat("outgroup ROOT;\n", file = con)
  cat("\n", file = con)
  cat("set theminchar 0;\n", file = con)
  cat("set themaxchar 1;\n", file = con)
  cat("\n", file = con)
  cat("cnames=;\n", file = con)
  cat("\n", file = con)

  # Generate rounds (replicate the structure search_replicates times)
  for (round_num in 1:search_replicates) {

    # Round header
    cat("/* --------------ROUNDS------------------- */\n", file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("                ROUND ", round_num, " - Starting iteration\n"), file = con)
    cat("============================================================================\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat(paste0("log ", output_dir, "Round_Output", round_num, ".txt;\n"), file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("                ROUND ", round_num, " - PAE-PCE Iteration\n"), file = con)
    cat("============================================================================\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)

    # Clear memory
    cat("/* Clear memory */\n", file = con)
    cat("keep 0 ;\n", file = con)
    cat("ttags - ;\n", file = con)
    cat("rseed 0 ; \n", file = con)
    cat("\n", file = con)

    # Parsimony search - DIFFERENT BASED ON METHOD
    if (search_method == "traditional") {
      cat("/* Traditional search */  \n", file = con)
      cat("quote Running parsimony analysis... ;\n", file = con)
      cat("mult 1000 = hold 100 tbr;\n", file = con)
      cat("bb=tbr;\n", file = con)
      cat("reroot [;\n", file = con)
    } else {
      cat("/* New tech search */  \n", file = con)
      cat("quote Running parsimony analysis... ;\n", file = con)
      cat("xmult= hits 7 rep 5 rss level 10;\n", file = con)
      cat("reroot [;\n", file = con)
    }
    cat("\n", file = con)

    # Save most parsimonious trees
    cat("/* Save most parsimonious trees */\n", file = con)
    cat(paste0("tsave *", output_dir, "Round_MPT", round_num, ".tre ;\n"), file = con)
    cat("save ;\n", file = con)
    cat("tsave / ;\n", file = con)
    cat("quote Most parsimonious trees saved ;\n", file = con)
    cat("\n", file = con)

    # Calculate strict consensus
    cat("/* Calculate strict consensus */\n", file = con)
    cat("if ((ntrees+1) > 1)\n", file = con)
    cat("    quote Calculating strict consensus... ;\n", file = con)
    cat("    nelsen*;      /*  calculate strict consensus   */\n", file = con)
    cat("    tchoose {strict};\n", file = con)
    cat(paste0("    tsave ", output_dir, "Round_Consensus", round_num, ".tmp ;\n"), file = con)
    cat("    save;\n", file = con)
    cat("    tsave/;      /*  close file with saved consensus   */\n", file = con)
    cat("    tt-;\n", file = con)
    cat("    ttags = ;\n", file = con)
    cat("    apo ] ;\n", file = con)
    cat(paste0("    ttag & ", output_dir, "Round_Map", round_num, ".svg ;\n"), file = con)
    cat("    taxname=;\n", file = con)
    cat(paste0("    export - ", output_dir, "Round_Tree", round_num, ".tre ;\n"), file = con)
    cat("    tt-;    \n", file = con)
    cat("    quote Consensus tree saved ;\n", file = con)
    cat("    keep 0;\n", file = con)
    cat("else\n", file = con)
    cat("    quote Last tree... ;\n", file = con)
    cat(paste0("    tsave ", output_dir, "Round_Consensus", round_num, ".tmp ;\n"), file = con)
    cat("    save 0;\n", file = con)
    cat("    tsave/;      /*  close file with saved consensus   */\n", file = con)
    cat("    tt-;\n", file = con)
    cat("    ttags = ;\n", file = con)
    cat("    apo ] ;\n", file = con)
    cat(paste0("    ttag & ", output_dir, "Round_Map", round_num, ".svg ;\n"), file = con)
    cat("    taxname=;\n", file = con)
    cat(paste0("    export - ", output_dir, "Round_Tree", round_num, ".tre ;\n"), file = con)
    cat("    tt-;    \n", file = con)
    cat("    quote Last tree saved ;\n", file = con)
    cat("    keep 0;\n", file = con)
    cat("end;    \n", file = con)
    cat("\n", file = con)

    # Calculate and save complete CI/RI statistics
    cat("/* ====================================================================== */\n", file = con)
    cat("/* CALCULATE AND SAVE COMPLETE CI/RI STATISTICS                          */\n", file = con)
    cat("/* ====================================================================== */\n", file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("----------------------------------------------------------------------------\n", file = con)
    cat("Tree Statistics:\n", file = con)
    cat("----------------------------------------------------------------------------\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat("/* get the tree */\n", file = con)
    cat(paste0("short ", output_dir, "Round_Consensus", round_num, ".tmp;\n"), file = con)
    cat("\n", file = con)
    cat("set THEMINtree minsteps;\n", file = con)
    cat("set THEMAXtree maxsteps;\n", file = con)
    cat("set THIStree length[0];\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "Round_CI_RI_Table", round_num, ".txt ;\n"), file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("        ROUND ", round_num, " - COMPLETE CI AND RI TABLE FOR ALL CHARACTERS\n"), file = con)
    cat("============================================================================\n", file = con)
    cat("\n", file = con)
    cat("Character | CI    | RI    | RC    | Status\n", file = con)
    cat("----------|-------|-------|-------|----------------------------------\n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)

    # Loop through all characters and print CI/RI
    cat("/* Loop through all characters and print CI/RI */\n", file = con)
    cat("loop 0 ('num_chars' - 1) ;    \n", file = con)
    cat("    set char_i length[0 #1];\n", file = con)
    cat("\n", file = con)
    cat("    if ('char_i' > 0)\n", file = con)
    cat("        set themaxchar maxsteps[#1] ;\n", file = con)
    cat("        set theminchar minsteps[#1] ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate CI */\n", file = con)
    cat("        if ('char_i' == 0)\n", file = con)
    cat("            set ci_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ci_val[#1] 'theminchar'/'char_i' ;\n", file = con)
    cat("        end;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RI */\n", file = con)
    cat("        if (('themaxchar'-'theminchar') == 0)\n", file = con)
    cat("            set ri_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ri_val[#1] ('themaxchar'-'char_i')/('themaxchar'-'theminchar') ;\n", file = con)
    cat("        end;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RC */\n", file = con)
    cat("        set rc_val[#1]  'ri_val[#1]'*'ci_val[#1]';\n", file = con)
    cat("\n", file = con)
    cat("        /* Print character info */\n", file = con)
    cat("        if ('ci_val[#1]' == 1)\n", file = con)
    cat("            if ('ri_val[#1]' == 1)\n", file = con)
    cat("                quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | NON-HOMOPLASTIC SYNAPOMORPHY ;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | Homoplastic ;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            quote #1 | 'ci_val[#1]' | 'ri_val[#1]' | 'rc_val[#1]' | Homoplastic ;\n", file = con)
    cat("        end;\n", file = con)
    cat("    else\n", file = con)
    cat("        /* Character is INACTIVE - skip */\n", file = con)
    cat("        quote #1 | INACTIVE | INACTIVE | INACTIVE | Character deactivated in previous round;\n", file = con)
    cat("    end;\n", file = con)
    cat("stop;\n", file = con)
    cat("\n", file = con)

    # Tree statistics
    cat("/* tree statistics (usando valores capturados anteriormente) */\n", file = con)
    cat("set CItree 'THEMINtree'/'THIStree';\n", file = con)
    cat("set RItree ('THEMAXtree'-'THIStree')/('THEMAXtree'-'THEMINtree');\n", file = con)
    cat("set RCtree 'RItree'*'CItree';\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "resul_CI_RI.log;\n"), file = con)
    cat("\n", file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("RUN_", round_num, ", 'CItree', 'RItree', 'RCtree', 'THIStree';\n"), file = con)
    cat("\n", file = con)
    cat("log/ ;\n", file = con)
    cat("\n", file = con)

    # Identify and list non-homoplastic synapomorphies
    cat("/* ====================================================================== */\n", file = con)
    cat("/* IDENTIFY AND LIST NON-HOMOPLASTIC SYNAPOMORPHIES                      */\n", file = con)
    cat("/* ====================================================================== */\n", file = con)
    cat("\n", file = con)
    cat(paste0("log+ ", output_dir, "Round_Synapomorphies", round_num, ".txt ;\n"), file = con)
    cat("quote\n", file = con)
    cat("============================================================================\n", file = con)
    cat(paste0("ROUND ", round_num, " - NON-HOMOPLASTIC SYNAPOMORPHIES (CI=1, RI=1)\n"), file = con)
    cat("============================================================================\n", file = con)
    cat("\n", file = con)
    cat("These characters (taxa/species) will be removed for the next round:\n", file = con)
    cat("\n", file = con)
    cat("Character(Taxon) | Run\n", file = con)
    cat("------------------------------------------------------------    \n", file = con)
    cat(";\n", file = con)
    cat("\n", file = con)
    cat("set removed_count 0;\n", file = con)
    cat("set conta 0;\n", file = con)
    cat("set contagem 0;\n", file = con)
    cat("\n", file = con)

    # Loop to identify and list synapomorphies
    cat("/* Loop to identify and list synapomorphies */\n", file = con)
    cat("loop 0 ('num_chars' - 1) ;\n", file = con)
    cat("    set char_i length[0 #1];\n", file = con)
    cat("\n", file = con)
    cat("    if ('char_i' > 0)\n", file = con)
    cat("        set themaxchar maxsteps[#1] ;\n", file = con)
    cat("        set theminchar minsteps[#1] ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate CI */\n", file = con)
    cat("        if ('char_i' == 0)\n", file = con)
    cat("            set ci_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ci_val[#1] 'theminchar'/'char_i' ;\n", file = con)
    cat("        end ;\n", file = con)
    cat("\n", file = con)
    cat("        /* Calculate RI */\n", file = con)
    cat("        if ( ('themaxchar'-'theminchar') == 0 )\n", file = con)
    cat("            set ri_val[#1] 1 ;\n", file = con)
    cat("        else\n", file = con)
    cat("            set ri_val[#1] ('themaxchar'-'char_i')/('themaxchar'-'theminchar') ;\n", file = con)
    cat("        end ;\n", file = con)
    cat("\n", file = con)
    cat("        set nome #1;\n", file = con)
    cat("\n", file = con)
    cat("        /* If CI=1 AND RI=1, it's a non-homoplastic synapomorphy */\n", file = con)
    cat("        if ('ci_val[#1]' == 1)\n", file = con)
    cat("            if ('ri_val[#1]' == 1)\n", file = con)
    cat(paste0("                quote Character 'nome' | Run ", round_num, ";\n"), file = con)
    cat("                set conta ++;\n", file = con)
    cat("                set removed_count ++ ;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote It is a homoplastic synapomorphy;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            quote It is a homoplastic synapomorphy;\n", file = con)
    cat("        end;\n", file = con)
    cat("    else\n", file = con)
    cat("        /* Character is INACTIVE - skip */\n", file = con)
    cat("        quote Character #1 is INACTIVE (deactivated in previous round);\n", file = con)
    cat("    end;            \n", file = con)
    cat("stop;\n", file = con)
    cat("\n", file = con)

    # If synapomorphies found, deactivate them
    cat("if ('removed_count' > 0)\n", file = con)
    cat("    var: caraExtract['conta'];\n", file = con)
    cat("    loop 0 ('num_chars' - 1);\n", file = con)
    cat("        set char_i length[0 #1];\n", file = con)
    cat("        if ('char_i' > 0)\n", file = con)
    cat("            if ('ci_val[#1]' == 1)\n", file = con)
    cat("                if ('ri_val[#1]' == 1)\n", file = con)
    cat("                    set caraExtract['contagem'] #1;\n", file = con)
    cat("                    set all_deactivated['total_deactivated'] #1;\n", file = con)
    cat("                    set total_deactivated ++;\n", file = con)
    cat("                    set contagem ++;\n", file = con)
    cat("                    set total_removed ++;\n", file = con)
    cat("                else\n", file = con)
    cat("                    quote Not counted;\n", file = con)
    cat("                end;\n", file = con)
    cat("            else\n", file = con)
    cat("                quote Not counted;\n", file = con)
    cat("            end;\n", file = con)
    cat("        else\n", file = con)
    cat("            /* Character is INACTIVE - skip */\n", file = con)
    cat("            quote Character #1 is INACTIVE - not counted;\n", file = con)
    cat("        end;        \n", file = con)
    cat("    stop;\n", file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat("    Total synapomorphies found: 'removed_count'\n", file = con)
    cat("    Characters removed:;\n", file = con)
    cat("    loop 0 ('conta' - 1) ;\n", file = con)
    cat("        quote character 'caraExtract[#1]';\n", file = con)
    cat("    stop;\n", file = con)
    cat("    \n", file = con)
    cat("    quote\n", file = con)
    cat("    ============================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("\n", file = con)
    cat("\n", file = con)
    cat("    set list_idx 0;\n", file = con)
    cat("    loop 0 ('total_deactivated' - 1) ;\n", file = con)
    cat("        ccode] 'all_deactivated['list_idx']';\n", file = con)
    cat("        quote Deactivating character 'all_deactivated['list_idx']';\n", file = con)
    cat("        set list_idx ++;\n", file = con)
    cat("    stop;\n", file = con)
    cat("\n", file = con)
    cat("    /* [OK] PROTECTION 2: Update num_chars after deactivation */\n", file = con)
    cat("    set num_chars nchar;\n", file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat(paste0("    Round ", round_num, " Summary:\n"), file = con)
    cat("    Characters removed this round: 'all_deactivated[0-('total_deactivated' - 1) &45]'\n", file = con)
    cat("    Total characters removed so far: 'total_removed'\n", file = con)

    # Add additional info for rounds after the first
    if (round_num > 1) {
      cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
      cat("    Active characters remaining: 'num_chars'\n", file = con)
    }

    cat("    ----------------------------------------------------------------------------\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)
    cat("\n", file = con)
    cat(paste0("    log+ ", output_dir, "PAE_PCE_Summary.txt;\n"), file = con)
    cat("\n", file = con)
    cat("    quote\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("                        REACHED ITERATION NUMBER ", round_num, "\n"), file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("    Maximum iterations (", max_iterations, ")\n"), file = con)
    cat("    Total characters removed: 'total_removed'\n", file = con)
    cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
    cat("\n", file = con)

    # Different message for last round vs. intermediate rounds
    if (round_num == search_replicates) {
      cat("    Note: Analysis has stopped due to iteration limit.\n", file = con)
    } else {
      cat("    Note: Analysis MIGHT BE stopped due to iteration limit.\n", file = con)
    }

    cat("    Consider increasing max_iterations if more rounds are needed.\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)

    # If no synapomorphies found, stop
    cat("else\n", file = con)
    cat("    log/;\n", file = con)
    cat("    quote\n", file = con)
    cat("    ============================================================================\n", file = con)
    cat(paste0("    STOPPING: No more non-homoplastic synapomorphies found in Round ", round_num, "\n"), file = con)
    cat("    ============================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    \n", file = con)
    cat(paste0("    log+ ", output_dir, "PAE_PCE_Summary.txt;\n"), file = con)
    cat("    quote\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("                        ANALYSIS COMPLETED SUCCESSFULLY\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat(paste0("    Stopped at iteration: ", round_num, "\n"), file = con)
    cat("    Reason: No more non-homoplastic synapomorphies (CI=1, RI=1)\n", file = con)
    cat("    Total characters removed: 'total_removed'\n", file = con)
    cat("    Total cumulative deactivated: 'total_deactivated'\n", file = con)
    cat("    ================================================================================\n", file = con)
    cat("    ;\n", file = con)
    cat("    log/;\n", file = con)
    cat("    silent-;\n", file = con)
    cat("    report=;\n", file = con)
    cat("    proc/;\n", file = con)
    cat("end;\n", file = con)

    # Add separator between rounds (except after last round)
    if (round_num < search_replicates) {
      cat("\n", file = con)
      cat(paste0("/* ---------------- RUN ", round_num + 1, " ----------------- */\n"), file = con)
      cat("set theminchar 0;\n", file = con)
      cat("set themaxchar 1;\n", file = con)
      cat("\n", file = con)
      cat("cnames=;\n", file = con)
      cat("\n", file = con)
    }
  }

  # Final closing
  cat("silent-;\n", file = con)
  cat("report=;\n", file = con)
  cat("proc/;", file = con)

  # Close connection
  close(con)

  message(paste0("[OK] TNT script generated successfully: ", output_file))
  message(paste0("   Total rounds: ", search_replicates))
  message(paste0("   Script will stop automatically when no more synapomorphies are found"))

  invisible(output_file)
}
#' Load Shapefiles from Output Directories
#'
#' Searches for and loads all shapefiles from output directories created by
#' range extrapolation functions (out_buffers, out_MCP, out_MST, etc.)
#'
#' @return A list containing:
#'   - shapefiles: named list of loaded shapefiles (terra objects)
#'   - rasters: named list of loaded rasters
#'   - directories: list of directories searched
#'
#' @export
load_output_shapefiles <- function() {
  
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out", "out_domains", "out_irregular_bins")
  
  shapefiles <- list()
  rasters <- list()
  found_dirs <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) {
      next
    }
    
    found_dirs <- c(found_dirs, dir)
    
    # Load shapefiles
    shp_files <- list.files(dir, pattern = "\\.shp$", full.names = TRUE)
    
    for (shp_file in shp_files) {
      tryCatch({
        shp_name <- tools::file_path_sans_ext(basename(shp_file))
        shapefiles[[paste0(dir, "/", shp_name)]] <- terra::vect(shp_file)
      }, error = function(e) {
        warning("Could not load shapefile: ", shp_file, " - ", e$message)
      })
    }
    
    # Load rasters
    tif_files <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
    
    for (tif_file in tif_files) {
      tryCatch({
        tif_name <- tools::file_path_sans_ext(basename(tif_file))
        rasters[[paste0(dir, "/", tif_name)]] <- terra::rast(tif_file)
      }, error = function(e) {
        warning("Could not load raster: ", tif_file, " - ", e$message)
      })
    }
  }
  
  return(list(
    shapefiles = shapefiles,
    rasters = rasters,
    directories = found_dirs
  ))
}

#' Convert Shapefile to Leaflet-compatible format
#'
#' Converts terra SpatVector to sf for use with leaflet
#'
#' @param shapefile A terra SpatVector object
#'
#' @return An sf object
#'
#' @export
shapefile_to_sf <- function(shapefile) {
  tryCatch({
    sf::st_as_sf(shapefile)
  }, error = function(e) {
    NULL
  })
}

#' Convert Raster to Leaflet-compatible format
#'
#' Converts terra raster to RasterLayer for use with leaflet
#'
#' @param raster A terra raster object
#'
#' @return A RasterLayer object
#'
#' @export
raster_to_raster_layer <- function(raster) {
  tryCatch({
    raster::raster(raster)
  }, error = function(e) {
    NULL
  })
}

#' Get list of available output directories
#'
#' @return Character vector of directories that exist
#'
#' @export
get_available_output_dirs <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out", "out_domains", "out_irregular_bins")
  available <- output_dirs[sapply(output_dirs, dir.exists)]
  return(available)
}

#' Get list of shapefile names from output directories
#'
#' @return Character vector of shapefile names
#'
#' @export
get_shapefile_names <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out", "out_domains", "out_irregular_bins")
  
  shp_names <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) next
    
    shp_files <- list.files(dir, pattern = "\\.shp$")
    shp_names <- c(shp_names, tools::file_path_sans_ext(shp_files))
  }
  
  return(shp_names)
}

#' Get list of raster names from output directories
#'
#' @return Character vector of raster names
#'
#' @export
get_raster_names <- function() {
  output_dirs <- c("out_buffers", "out_MCP", "out_MST", "out", "out_domains", "out_irregular_bins")
  
  tif_names <- c()
  
  for (dir in output_dirs) {
    if (!dir.exists(dir)) next
    
    tif_files <- list.files(dir, pattern = "\\.tif$")
    tif_names <- c(tif_names, tools::file_path_sans_ext(tif_files))
  }
  
  return(tif_names)
}
#' Calculate AIC and AICc
#'
#' Calculate Akaike Information Criterion and corrected AIC
#'
#' @param lnl Log-likelihood value
#' @param k Number of parameters
#' @param n Sample size (number of areas)
#'
#' @return List with aic and aicc values
#'
#' @export
calc_aic_aicc <- function(lnl, k, n) {
  aic <- 2 * k - 2 * lnl
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  
  list(
    aic = aic,
    aicc = aicc
  )
}

#' Calculate Delta AIC and Delta AICc
#'
#' Calculate differences in AIC and AICc relative to the best model
#'
#' @param aic Vector of AIC values
#' @param aicc Vector of AICc values
#'
#' @return List with delta_aic and delta_aicc vectors
#'
#' @export
calc_delta_aic <- function(aic, aicc) {
  delta_aic <- aic - min(aic, na.rm = TRUE)
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  
  list(
    delta_aic = delta_aic,
    delta_aicc = delta_aicc
  )
}

#' Perform Likelihood Ratio Test
#'
#' Test if a complex model is significantly better than a simple model
#'
#' @param lnl_simple Log-likelihood of simple model
#' @param lnl_complex Log-likelihood of complex model
#' @param df Degrees of freedom (difference in number of parameters)
#'
#' @return List with test statistic, p-value, and significance
#'
#' @export
perform_lrt <- function(lnl_simple, lnl_complex, df) {
  # LRT statistic: 2 * (lnl_complex - lnl_simple)
  # Follows chi-square distribution with df degrees of freedom
  
  if (lnl_complex < lnl_simple) {
    warning("Complex model has lower likelihood than simple model")
  }
  
  test_stat <- 2 * (lnl_complex - lnl_simple)
  pval <- stats::pchisq(test_stat, df = df, lower.tail = FALSE)
  
  list(
    test_statistic = test_stat,
    df = df,
    p_value = pval,
    significant = pval < 0.05
  )
}

#' Create Model Comparison Table
#'
#' Create a comprehensive comparison table for multiple models
#'
#' @param models Character vector of model names
#' @param lnl Numeric vector of log-likelihood values
#' @param params List of parameter values for each model
#' @param n_areas Number of areas (for AICc calculation)
#'
#' @return Data frame with model comparison
#'
#' @export
create_model_comparison <- function(models, lnl, params, n_areas) {
  
  # Calculate number of parameters for each model
  n_params <- sapply(params, function(x) {
    n <- 2  # d and e
    if (!is.na(x$j) && x$j > 0) n <- n + 1  # Add j if present
    n
  })
  
  # Calculate AIC and AICc
  aic <- sapply(seq_along(models), function(i) {
    calc_aic_aicc(lnl[i], n_params[i], n_areas)$aic
  })
  
  aicc <- sapply(seq_along(models), function(i) {
    calc_aic_aicc(lnl[i], n_params[i], n_areas)$aicc
  })
  
  # Calculate delta AIC and delta AICc
  delta_aic <- aic - min(aic, na.rm = TRUE)
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  
  # Create comparison table
  comparison <- data.frame(
    Model = models,
    LnL = round(lnl, 3),
    K = n_params,
    AIC = round(aic, 2),
    DeltaAIC = round(delta_aic, 2),
    AICc = round(aicc, 2),
    DeltaAICc = round(delta_aicc, 2),
    stringsAsFactors = FALSE
  )
  
  # Add parameter values
  comparison$d <- round(sapply(params, function(x) x$d), 4)
  comparison$e <- round(sapply(params, function(x) x$e), 4)
  comparison$j <- sapply(params, function(x) {
    if (is.na(x$j)) "-" else round(x$j, 4)
  })
  
  # Sort by AICc
  comparison <- comparison[order(comparison$AICc), ]
  rownames(comparison) <- NULL
  
  comparison
}

#' Calculate Akaike Weights
#'
#' Calculate Akaike weights for model comparison
#'
#' @param aicc Vector of AICc values
#'
#' @return Vector of Akaike weights
#'
#' @export
calc_akaike_weights <- function(aicc) {
  delta_aicc <- aicc - min(aicc, na.rm = TRUE)
  L <- exp(-0.5 * delta_aicc)
  w <- L / sum(L, na.rm = TRUE)
  w
}

#' Generate Mock BioGeoBEARS Results
#'
#' Generate realistic mock results for demonstration
#'
#' @param models Character vector of model names
#' @param n_areas Number of areas
#' @param d_bounds Numeric vector with min and max for d parameter
#' @param e_bounds Numeric vector with min and max for e parameter
#' @param j_bounds Numeric vector with min and max for j parameter
#'
#' @return List with lnl and params
#'
#' @export
generate_mock_results <- function(models, n_areas, d_bounds, e_bounds, j_bounds) {
  
  lnl <- numeric(length(models))
  params <- list()
  
  for (i in seq_along(models)) {
    # Generate realistic log-likelihood values
    # Better models have higher (less negative) lnl
    base_lnl <- -50 - (i - 1) * 5
    lnl[i] <- base_lnl + rnorm(1, 0, 2)
    
    # Generate parameter values
    d <- runif(1, d_bounds[1], d_bounds[2])
    e <- runif(1, e_bounds[1], e_bounds[2])
    j <- NA
    
    if (grepl("\\+J", models[i])) {
      j <- runif(1, j_bounds[1], j_bounds[2])
    }
    
    params[[i]] <- list(d = d, e = e, j = j)
  }
  
  list(
    lnl = lnl,
    params = params
  )
}
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
#' pae_pce
#'
#' Produces generalized tracks with PAE-PCE, from the presence-absence matrix produced by the MST_node function.
#'
#' @param preabsMat the presence-absence matrix from MST_node function
#' @param shapeFile a polygon of some place already with datum to be used
#' @param resol a vector with resolution of quadrats
#' @param sobrepo logical value that enables plotting on a map with minimum spanning trees (MSTs) already
#' produced
#' @param nonHomoplasticSpeciesList logical value that enables the production of a table (in .txt) showing
#' the presence of non-homoplastic species and its referred number of grid in a /out directory automatically produced
#' @param gridView logical value that enables viewing a grid with a given resolution provided by
#' the resol argument
#' @param labelGrid logical value that enables viewing the numbers of a grid with a given resolution provided by
#' the resol argument
#' @param N number of iterations in order to get the most parsimonious tree
#' @param ... optional plot arguments
#' @return Matrix with homoplastic and non-homoplastic species and their grid numbers associated, as well as a plot
#' with the polygon indicated by the shapeFile argument and the hatched grids indicating the generalized track
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details This function produces generalized tracks, which can be used as primary biogeographical homology
#' hypothesis in Biogeography (Page, 1987; Castillo-Garcia \emph{et al.}, 2025). It is defined as the superimposition
#' of at least two minimum spanning trees, which are graphs with \emph{n} localities connected through \emph{n}-1 edges, producing shortest
#' total length graphs.
#' @references Castillo-Garcia, C.F., Morrone, J.J., Salgado-Ugarte, I.H. & D. Espinosa, 2025. Panbiotracks: software for track analysis.
#'  \emph{Revista Mexicana de Biodiversidad} \strong{96}: e965429.
#'
#'  Page, R.D.M., 1987. Graphs and generalized tracks: quantifying Croizat's panbiogeography.
#'   \emph{Systematic Zoology} \strong{36}: 1.
#' @seealso
#' \code{\link[phytools]{phytools}} for phylogenetic trees,
#' \code{\link[terra]{terra}} for spatial data handling,
#' \code{\link[PanBioGeo]{MST_node}} for the production of presence-absence matrices.
#' @export
#' @examples
#' # First, determine the minimum spanning trees of all taxa and produce a presence-absence matrix
#' # ---------------------------- with all taxa ------------------------------
#' # ----------- with a tree and resolution of 10 X 10 degrees ----------
#' mst_all_taxa <- MST_node(coordin = resul$data_df, shape_file = lycipta.asul, sobrepo = F, caption = T,
#'   resol = c(10, 10), tree = resul$treeMod)
#'
#' # Now, produce the generalized track and a list with grid numbers and the characters (i.e., taxaset)
#' # Please try to run at least more than one round, with N > 1...
#' # With a new plot using the polygon of South America (asul)
#' pae1 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10), gridView = TRUE, labelGrid = TRUE,
#'    nonHomoplasticSpeciesList = FALSE, N = 10)
#' pae1 # data frame showing the homoplastic and non-homoplastic species (spp) and their grid numbers (grid_n),
#' # from the anterior iterations (PAE-PCE)
#'
#' # Plotting onto a map already plotted with minimum spanning trees
#' pae2 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10),
#'                gridView = TRUE, labelGrid = TRUE, nonHomoplasticSpeciesList = TRUE, sobrepo = TRUE,
#'                N = 10)
#' pae2
#'
#' # Doing a zoom
#' pae3 <- pae_pce(preabsMat = mst_all_taxa, shapeFile = asul, resol = c(10, 10),
#'                gridView = TRUE, labelGrid = TRUE, nonHomoplasticSpeciesList = TRUE, sobrepo = TRUE,
#'                N = 10, xmin = min(lycipta.coords$Long) - 5,
#'                xmax = max(lycipta.coords$Long), ymin = min(lycipta.coords$Lat) - 7,
#'                ymax = max(lycipta.coords$Lat))
#' pae3
#'

pae_pce <- function(preabsMat, shapeFile, resol, N = NULL,
                    gridView = FALSE, labelGrid = FALSE, nonHomoplasticSpeciesList = TRUE, sobrepo = FALSE,
                    xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL){

  library(phangorn)
  library(raster)
  library(viridis)
  library(TreeSearch) # TBR search

  ## HELPER FUNCTION FOR ELEGANT TERMINATIONS ##
  finalizar_analise_elegante <- function(reason, details = NULL, iteration = 1, plot_grid = TRUE) {
    # Function to gracefully finalize the analysis with informative reports

    print("")
    print("")
    print(paste(rep("=", 60), collapse=""))
    print("FINALIZED PAE-PCE ANALYSIS")
    print(paste(rep("=", 60), collapse=""))

    print(paste0("Stopping iteration: ", iteration))
    print(paste0("Reason for termination: ", reason))

    if(!is.null(details)) {
      print("Details:")
      for(detail in details) {
        print(paste0("  * ", detail))
      }
    }

    print("")
    print("Scientific interpretation:")
    if(grepl("synapomorphy", reason, ignore.case = TRUE)) {
      print("* This is a normal completion of the PAE-PCE analysis")
      print("* Indicates that there are no more informative characters for analysis")
      print("* All non-homoplastic synapomorphies were identified")
    } else if(grepl("overlap", reason, ignore.case = TRUE)) {
      print("* The analyzed species do not have sufficiently overlapping distributions")
      print("* Consider adjusting the resolution of the grids or the overlap criteria")
      print("* It may indicate that the species have naturally distinct distributions")
    } else {
      print("* This is a normal termination based on the analysis criteria")
      print("* The data do not meet the minimum requirements for generalized traits")
    }

    print("")
    print("Recommendations:")
    print("* Check the quality of the input data")
    print("* Consider adjusting resolution parameters")
    print("* Assess whether the data are suitable for PAE-PCE analysis")

    print(paste(rep("=", 60), collapse=""))

    # Plot grid if requested
    if(plot_grid && exists("cropped_map")) {
      plot(cropped_map, add = TRUE, border = "gray", lwd = 0.5)
    }

    # Return structured information instead of throwing error
    return(list(
      status = "finished_without_results",
      reason = reason,
      details = details,
      iteration_parada = iteration,
      timestamp = Sys.time(),
      message = "PAE-PCE analysis completed - please refer to the report for details"
    ))
  }

  ## HELPER FUNCTION TO CHECK GRID OVERLAP - CRITERION 1.0 ##
  verificar_sobreposicao_grids <- function(matTemp2, min_sobreposicao = 1.0) {
    # Identifies groups of species with PERFECT overlap (1.0 = identical distributions)

    print("=== CHECKING GROUPS WITH PERFECT OVERLAP ===")

    # Get list of grids for each species
    especies_grids <- list()
    for(especie in colnames(matTemp2)) {
      grids_especie <- rownames(matTemp2)[matTemp2[, especie] == 1]
      especies_grids[[especie]] <- grids_especie
      print(paste0(especie, ": grids ", paste(grids_especie, collapse = ", ")))
    }

    especies <- names(especies_grids)
    n_especies <- length(especies)

    if(n_especies < 2) {
      print("Only 1 species - automatically accepted")
      return(list(especies_validas = especies, sobreposicao_ok = TRUE))
    }

    # Calculate overlap matrix among all species
    sobreposicoes <- matrix(0, nrow = n_especies, ncol = n_especies,
                            dimnames = list(especies, especies))

    for(i in 1:n_especies) {
      sobreposicoes[i, i] <- 1  # Self-overlap = 1
    }

    for(i in 1:(n_especies-1)) {
      for(j in (i+1):n_especies) {
        esp1 <- especies[i]
        esp2 <- especies[j]
        grids1 <- especies_grids[[esp1]]
        grids2 <- especies_grids[[esp2]]

        # Calculate overlap (Jaccard index)
        intersecao <- length(intersect(grids1, grids2))
        uniao <- length(union(grids1, grids2))
        sobreposicao <- if(uniao > 0) intersecao / uniao else 0

        sobreposicoes[i, j] <- sobreposicao
        sobreposicoes[j, i] <- sobreposicao

        print(paste0("Overlap ", esp1, " vs ", esp2, ": ", round(sobreposicao, 3)))

        # STRICT CRITERION: Only overlap = 1.0 is accepted
        if(sobreposicao == 1.0) {
          print(paste0("  [OK] PERFECT OVERLAP (= 1.0) - Identical distributions"))
        } else {
          print(paste0("  [X] Imperfect overlap (", round(sobreposicao, 3), " != 1.0) - Different distributions"))
        }
        print(paste0("  - Common grids: ", paste(intersect(grids1, grids2), collapse = ", ")))
        if(sobreposicao < 1.0) {
          print(paste0("  - Unique grids of ", esp1, ": ", paste(setdiff(grids1, grids2), collapse = ", ")))
          print(paste0("  - Unique grids of ", esp2, ": ", paste(setdiff(grids2, grids1), collapse = ", ")))
        }
      }
    }

    # LOGIC: Find groups of species with PERFECT overlap (1.0)
    print("")
    print("=== IDENTIFYING GROUPS WITH IDENTICAL DISTRIBUTIONS ===")

    # Create connection graph (species connected ONLY if overlap = 1.0)
    conexoes <- sobreposicoes == 1.0

    # Find connected components (species groups)
    visitado <- rep(FALSE, n_especies)
    grupos <- list()
    grupo_id <- 1

    for(i in 1:n_especies) {
      if(!visitado[i]) {
        # Depth-first search to find all connected species
        grupo_atual <- c()
        pilha <- c(i)

        while(length(pilha) > 0) {
          atual <- pilha[length(pilha)]
          pilha <- pilha[-length(pilha)]

          if(!visitado[atual]) {
            visitado[atual] <- TRUE
            grupo_atual <- c(grupo_atual, atual)

            # Add unvisited neighbors to stack
            vizinhos <- which(conexoes[atual, ] & !visitado)
            pilha <- c(pilha, vizinhos)
          }
        }

        if(length(grupo_atual) >= 2) {
          grupos[[grupo_id]] <- especies[grupo_atual]
          print(paste0("Group ", grupo_id, " (", length(grupo_atual), " species with identical distributions): ",
                       paste(especies[grupo_atual], collapse = ", ")))
          grupo_id <- grupo_id + 1
        } else {
          print(paste0("Isolated species (unique distribution): ", especies[grupo_atual]))
        }
      }
    }

    # Determine valid species (those belonging to groups >=2 species with overlap = 1.0)
    especies_validas <- c()
    for(grupo in grupos) {
      especies_validas <- c(especies_validas, grupo)
    }

    if(length(especies_validas) >= 2) {
      print(paste0("", "[OK] FOUND ", length(grupos), " GROUP(S) WITH IDENTICAL DISTRIBUTIONS"))
      print(paste0("Accepted species (identical distributions): ", paste(especies_validas, collapse = ", ")))
      return(list(especies_validas = especies_validas, sobreposicao_ok = TRUE,
                  sobreposicoes = sobreposicoes, grupos = grupos))
    } else {
      print("")
      print("[X] NO VALID GROUP FOUND (no group with >=2 species of identical distributions)")
      return(list(especies_validas = c(), sobreposicao_ok = FALSE,
                  sobreposicoes = sobreposicoes, grupos = list()))
    }
  }

  # ============================================================================
  # DETECT GRID NAME TYPE AND CREATE MAPPING IF NEEDED
  # ============================================================================

  # Check if rownames are numeric or strings
  grid_names <- rownames(preabsMat)
  grid_names_no_root <- grid_names[grid_names != "ROOT"]

  # Detect if grid names are numeric
  is_numeric_grid <- all(suppressWarnings(!is.na(as.numeric(grid_names_no_root))))

  if(is_numeric_grid) {
    cat("Grid names detected as NUMERIC\n")
    use_grid_mapping <- FALSE
    grid_mapping <- NULL
  } else {
    cat("Grid names detected as STRINGS (area names)\n")
    cat("Creating internal mapping: area names <-> numeric IDs\n")
    use_grid_mapping <- TRUE

    # Create mapping: area_name -> numeric_id
    grid_mapping <- data.frame(
      area_name = grid_names_no_root,
      grid_id = 1:length(grid_names_no_root),
      stringsAsFactors = FALSE
    )

    cat(paste0("Mapped ", nrow(grid_mapping), " areas to numeric IDs\n"))
    cat("First 5 mappings:\n")
    print(head(grid_mapping, 5))
  }

  # ============================================================================
  # Convert the character matrix to a phyDat object
  # in order to infer a tree. By assigning the matrix
  # type as 'user' we can specify the components of the
  # matrix (in this case, binary 1s and 0s).
  # ============================================================================

  tempMatrix <- as.phyDat(preabsMat, type = 'USER', levels = c(0, 1))

  colu <- ncol(preabsMat)

  # Creating a restriction:
  ciVec <- rep(1.0, colu)
  riVec <- rep(1.0, colu)


  lista <- list()
  listaR <- list()

  print('Please, remember that the resulting rasters will be kept in your setted directory!')

  if(sobrepo == TRUE){
    print('Please do not forgive that grids are in different colors!')
    if(is.null(N)){
      # ELEGANT STOP: Mandatory parameter not provided
      return(finalizar_analise_elegante(
        reason = "Parameter N not specified",
        details = c("Parameter N (number of iterations) is mandatory when sobrepo=TRUE",
                    "Specify N=10 or another suitable value"),
        iteration = 0,
        plot_grid = FALSE
      ))
    } else {
      print('The N argument is equal to 1 as default!')
    }
  }

  ######################
  ##### shape file #####

  # Store original shapefile for irregular bins
  shapeFile_sf <- if(inherits(shapeFile, "sf")) shapeFile else st_as_sf(shapeFile)

  # If using irregular bins (string grid names), find matching column in shapefile
  bin_column_name <- NULL
  grid_to_area_mapping <- NULL

  if(use_grid_mapping) {
    cat("\n=== IRREGULAR BINS MODE ===")
    cat("\nSearching for shapefile column matching matrix rownames...\n")

    # Try to find column in shapefile that matches grid names
    for(col in names(shapeFile_sf)) {
      if(col == "geometry") next

      shapefile_values <- as.character(shapeFile_sf[[col]])
      matrix_rownames <- grid_mapping$area_name

      # Check if there's overlap
      overlap <- sum(matrix_rownames %in% shapefile_values)

      if(overlap > 0) {
        cat(paste0("  Found matching column: '", col, "' (", overlap, "/", length(matrix_rownames), " matches)\n"))
        bin_column_name <- col
        break
      }
    }

    if(is.null(bin_column_name)) {
      cat("  WARNING: Could not find shapefile column matching matrix rownames.\n")
      cat("  The shapefile will be used only for masking, not for irregular bins association.\n")
      cat("  Proceeding with regular grid mode...\n")
      cat("===========================\n\n")
      # Disable irregular bins mode
      use_grid_mapping <- FALSE
      grid_mapping <- NULL
    } else {
      cat(paste0("Using shapefile column: '", bin_column_name, "'\n"))
      cat("Regular grids will be created and associated with irregular bins\n")
      cat("===========================\n\n")
    }
  }

  shapeFile <- as(shapeFile, 'Spatial')

  # Create regular grid
  grid <- raster(extent(shapeFile), resolution = resol, crs = CRS("+proj=longlat +datum=WGS84"))
  grid <- raster::extend(grid, c(1, 1))
  gridPolygon <- rasterToPolygons(grid)
  crs(gridPolygon) <- "+proj=longlat +datum=WGS84" # datum WGS84

  # If using irregular bins, create mapping between grid cells and areas
  if(use_grid_mapping) {
    cat("Creating spatial association between grid cells and irregular bins...\n")

    # Convert gridPolygon to sf
    gridPolygon_sf <- st_as_sf(gridPolygon)
    gridPolygon_sf$grid_id <- 1:nrow(gridPolygon_sf)

    # Spatial join: which grid cell intersects which polygon?
    grid_area_join <- st_join(gridPolygon_sf, shapeFile_sf[, c(bin_column_name, "geometry")],
                              join = st_intersects, largest = TRUE)

    # Create mapping: grid_id -> area_name
    grid_to_area_mapping <- data.frame(
      grid_id = grid_area_join$grid_id,
      area_name = as.character(grid_area_join[[bin_column_name]]),
      stringsAsFactors = FALSE
    )

    # Remove NAs (grids outside all polygons)
    grid_to_area_mapping <- grid_to_area_mapping[!is.na(grid_to_area_mapping$area_name), ]

    cat(paste0("  Associated ", nrow(grid_to_area_mapping), " grid cells to irregular bins\n"))
    cat(paste0("  Unique areas: ", length(unique(grid_to_area_mapping$area_name)), "\n"))

    # Show summary
    area_counts <- table(grid_to_area_mapping$area_name)
    cat("\n  Grid cells per area:\n")
    for(area in names(area_counts)) {
      cat(paste0("    ", area, ": ", area_counts[area], " cells\n"))
    }
    cat("\n")
  }


  # Build one canonical grid mask from polygon coverage (>0)
  # to keep matrix and map using identical grid membership.
  mask.raster <- raster(extent(shapeFile), resolution = resol,
                        crs = CRS("+proj=longlat +datum=WGS84"))
  r <- rasterize(shapeFile, mask.raster, getCover = TRUE)
  r[is.na(r) | r <= 0] <- NA
  r[!is.na(r)] <- 1
  proj4string(r) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
  r <- merge(r, mask.raster)

  cropped_map <- rasterToPolygons(r, dissolve = FALSE)
  crs(cropped_map) <- "+proj=longlat +datum=WGS84"


  if(sobrepo == FALSE){
    colores <- viridis(n = length(r[!is.na(r)]), option = 'A')
  } else if(sobrepo == TRUE){
    colores <- viridis(n = N, option = 'A')
  }

  # Plotting the shapefile
  if(is.null(xmin)){
    plot(shapeFile, axes = TRUE, las = 1)
  } else if(!is.null(xmin)){
    plot(shapeFile, axes = TRUE, las = 1, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  }

  if (!file.exists('out/'))
    dir.create('out/') # output directory for results

  #---------------------------------- #### PAE-PCE ### -------------------------------------#
  contagem <- 1
  conta <- 0

  # MAIN CORRECTION: Variable to control when to stop the loop
  deve_continuar <- TRUE

  # CORRECTION: More robust while condition
  while(deve_continuar && length(which((ciVec == 1) == TRUE)) > 1 && dim(preabsMat)[2] > 1){

    print(paste0("=== STARTING ITERATION ", contagem, " ==="))
    print(paste0("Remaining species in matrix: ", dim(preabsMat)[2]))
    print(paste0("Grids in matrix: ", dim(preabsMat)[1]))

    # We can get a random starting tree for a parsimony search using the
    # random.addition function:

    if(dim(as.data.frame(tempMatrix))[1] < 2){
      # ELEGANT STOP: Insufficient data
      return(finalizar_analise_elegante(
        reason = "Insufficient data for phylogenetic analysis",
        details = c("Less than 2 terminals in the data matrix",
                    "Parsimony analysis requires at least 2 terminals"),
        iteration = contagem
      ))
    }

    ### TBR + RATCHET searches
    dm <- dist.hamming(tempMatrix)
    ra.tre <- NJ(dm)

    # Stores the score (length) of the best tree found so far
    best_score <- parsimony(ra.tre, tempMatrix)
    cat(paste("Score inicial:", best_score, "\n"))

    # Use the pratchet function (parsimony ratchet) to find the most
    # parsimonious tree. Specifying k means the algorithm will search
    # through k possible trees to find the most parsimonious solution.
    treeIt1 <- list()

    get_single_tree <- function(tree_obj, data_mat) {
      if (inherits(tree_obj, "multiPhylo")) {
        if (length(tree_obj) == 0) return(NULL)
        scores_tmp <- sapply(tree_obj, function(tt) parsimony(tt, data_mat))
        return(tree_obj[[which.min(scores_tmp)[1]]])
      }
      tree_obj
    }

    for(rou in 1:N){
      teste2 <- pratchet(tempMatrix,
                         start = ra.tre,
                         maxit = 1000,
                         minit = 100,
                         k = 20,
                         trace = 1,
                         method = 'fitch',
                         perturbation = 'ratchet',
                         rearrangements = 'TBR')

      teste2 <- optim.parsimony(tree = teste2,
                                data = tempMatrix,
                                method = 'fitch',
                                rearrangements = 'TBR')

      teste2 <- get_single_tree(teste2, tempMatrix)
      if (is.null(teste2)) {
        next
      }

      # Evaluate the score of the new tree
      current_score <- attr(teste2, "pscore")
      if(is.null(current_score) || length(current_score) == 0 || !is.finite(current_score)) {
        current_score <- parsimony(teste2, tempMatrix)
      }
      current_score <- as.numeric(current_score)[1]

      cat(paste("Best score found in this iteration:", current_score, "\n"))

      # Use current tree as next start to improve search trajectory
      ra.tre <- teste2

      # If we find a better tree (lower score), reset list of best trees
      if(current_score < best_score){
        cat(paste("Novo melhor score encontrado! Atualizando de", best_score, "para", current_score, "\n"))
        best_score <- current_score
        treeIt1 <- list(teste2)
      } else if(current_score == best_score){
        cat("Score equal to the best already found. Adding trees to the list.\n")
        treeIt1[[length(treeIt1) + 1]] <- teste2
      }
    }
    cat('\nGenerating majority consensus tree...\n')

    treeIt1 <- treeIt1[!vapply(treeIt1, is.null, logical(1))]
    if (length(treeIt1) == 0) {
      treeIt1 <- list(ra.tre)
    }

    class(treeIt1) <- c("multiPhylo", "list")
    treeIt1 <- tryCatch(unique(treeIt1), error = function(e) treeIt1)

    if (length(treeIt1) > 1) {
      consensus_tree <- consensus(treeIt1, p = 0.5)
      write.tree(phy = consensus_tree, file = 'out/consensus_pae_pce.tre',
                 append = TRUE, tree.names = paste0('run_', contagem))
      cat(paste("Consensus generated from", length(treeIt1), "trees\n"))
    } else {
      consensus_tree <- treeIt1[[1]]
      write.tree(phy = consensus_tree, file = 'out/consensus_pae_pce.tre',
                 append = TRUE, tree.names = paste0('run_', contagem))
      cat("Only one tree found. Returning this tree.\n")
    }

    if(is.null(treeIt1)){
      # ELEGANT STOP: Tree construction failure
      return(finalizar_analise_elegante(
        reason = "Failure in phylogenetic topology construction",
        details = c("Search algorithm could not generate valid trees",
                    "May indicate problems in the input data"),
        iteration = contagem
      ))
    }

    treeT <- consensus_tree # 'majority rule': change p to 0.5

    # Root the tree by the designated outgroup
    # (write the species name as it appears in the tree,
    # but add an underscore to fill the space).
    strictTree <- root(treeT, outgroup = "ROOT", resolve.root = TRUE)

    # CI of each character (i.e., taxon):
    ciVec <- CI(strictTree, data = tempMatrix, sitewise = TRUE)

    # RI of each character (i.e., taxon):
    riVec <- RI(strictTree, data = tempMatrix, sitewise = TRUE)
    riVec[which(is.na(riVec))] <- 0

    print(paste0("CI calculated for ", length(ciVec), " characters"))
    print(paste0("RI calculated for ", length(riVec), " characters"))
    print(paste0("Characters with CI=1: ", length(which(ciVec == 1))))
    print(paste0("Characters with RI=1: ", length(which(riVec == 1))))

    if(all(is.na(riVec) == TRUE)){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }

    if(any(riVec == 1.0) == FALSE){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }


    ## First load packages
    require(phytools)
    tree_temp <- strictTree
    if(is.null(tree_temp$edge.length) == TRUE){
      print('...creating the branch lengths of a tree equal to one...')
      tree_temp <- compute.brlen(tree_temp, 1)
    }

    # Discrete mapping
    inter1 <- which(ciVec == 1)
    inter2 <- which(riVec == 1)
    interT1 <- intersect(inter1, inter2)

    print(paste0("Characters with CI=1 and RI=1 (non-homoplastic synapomorphies): ", length(interT1)))

    if(length(interT1) == 0 && contagem == 1){
      # ELEGANT STOP: No non-homoplastic synapomorphy in first iteration
      return(finalizar_analise_elegante(
        reason = "No non-homoplastic synapomorphy found",
        details = c("First iteration did not identify characters with CI=1 and RI=1",
                    "Data may not be suitable for PAE-PCE analysis",
                    "Consider checking the quality of the input data"),
        iteration = contagem
      ))
    }

    if(length(interT1) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }

    # Things that had IC and IR equal to 1...
    lista_T <- list()
    matTemp2 <- preabsMat[, interT1] # temporary matrix with non-homoplastic synapomorphies
    matTemp2 <- as.matrix(matTemp2)

    nomesCOLsim <- colnames(matTemp2) # species which are synapomorphies
    print(paste0("Synapomorphic species identified: ", paste(nomesCOLsim, collapse = ", ")))

    if(length(nomesCOLsim) < 2){
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break # remaining body of loop will not be executed!
    }

    resulT <- dim(matTemp2)[2]

    # NEW VERIFICATION: Strict overlap of distributions
    print(paste0('Iteration ', contagem, ': Identifying groups of species with identical distributions...'))

    verificacao <- verificar_sobreposicao_grids(matTemp2, min_sobreposicao)

    # CRITICAL CORRECTION: Stop iteration if no perfect overlap found
    if(!verificacao$sobreposicao_ok || length(verificacao$especies_validas) < 2) {
      # No valid overlap group found
      if(contagem == 1) {
        # First iteration without perfect overlap - INVALIDATE
        return(finalizar_analise_elegante(
          reason = "No valid group of species with overlapping distributions",
          details = c("First iteration found no species with identical distributions",
                      paste0("Species analyzed: ", paste(colnames(matTemp2), collapse = ", ")),
                      "Consider adjusting grid resolution or overlap criteria"),
          iteration = contagem
        ))
      } else {
        # Subsequent iterations (2+) - ACCEPT previous iterations as VALID and STOP
        print("")
        print(paste(rep("=", 60), collapse=""))
        print(paste0("ITERATION ", contagem, ": NO PERFECT OVERLAP FOUND"))
        print(paste(rep("=", 60), collapse=""))
        print("No group with identical distributions (overlap = 1.0) was found.")
        print("Previous iterations are VALID. Stopping analysis here.")
        print(paste0("Species analyzed: ", paste(colnames(matTemp2), collapse = ", ")))
        print(paste(rep("=", 60), collapse=""))

        # Incrementing contagem to signal that at least one iteration was successful:
        # contagem <- contagem + 1
        deve_continuar <- FALSE
        break  # STOP, but previous iterations remain VALID
      }
    }

    # Filter only species with identical distributions
    especies_validas <- verificacao$especies_validas
    matTemp2 <- matTemp2[, especies_validas, drop = FALSE]

    print(paste0('Accepted species (identical distributions): ', paste(especies_validas, collapse = ", ")))

    # Continue with normal processing

    # Common grids (synapomorphies)
    roundRes1 <- rownames(matTemp2[rowSums(matTemp2) <= resulT, , drop = F])
    roundRes2 <- rownames(matTemp2[rowSums(matTemp2) > 1, , drop = F])
    roundRes <- intersect(roundRes1, roundRes2)

    print(paste0("Common grids identified: ", length(roundRes)))
    if(length(roundRes) > 0) {
      print(paste0("Grids: ", paste(roundRes, collapse = ", ")))
    }

    if(is.null(roundRes) && contagem == 1){
      # ELEGANT STOP: No common grid in first iteration
      return(finalizar_analise_elegante(
        reason = "No common grid identified as synapomorphy",
        details = c("First iteration found no shared grids among species",
                    "Species may have very dispersed distributions",
                    "Consider adjusting grid resolution"),
        iteration = contagem
      ))
    } else if(is.null(roundRes) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }

    ################################################################
    ###########################
    ## To see the mapping by locality of the species that form
    # the generalized track(s):
    # After the first round and exclusion of characters with IC and IR equal to 1...
    lista_TT <- list()

    # Common species which are synapomorphies
    nLi <- dim(matTemp2)[1]
    roundResCol1 <- colnames(matTemp2[, colSums(matTemp2) <= nLi , drop = F])
    roundResCol2 <- colnames(matTemp2[, colSums(matTemp2) > 1 , drop = F])
    roundResCol <- intersect(roundResCol1, roundResCol2)

    print(paste0("Common synapomorphic species identified: ", length(roundResCol)))
    if(length(roundResCol) > 0) {
      print(paste0("Species: ", paste(roundResCol, collapse = ", ")))
    }

    if(is.null(roundResCol) && contagem == 1){
      # ELEGANT STOP: No common synapomorphic species in first iteration
      return(finalizar_analise_elegante(
        reason = "No common synapomorphic species identified",
        details = c("First iteration found no shared species among grids",
                    "May indicate that species have very specific distributions",
                    "Consider checking species distribution data"),
        iteration = contagem
      ))
    } else if(is.null(roundResCol) && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }

    if(length(roundResCol) < 2 && contagem == 1){
      # ELEGANT STOP: Less than 2 synapomorphic species in first iteration
      return(finalizar_analise_elegante(
        reason = "Insufficient synapomorphic species for analysis",
        details = c("First iteration found less than 2 synapomorphic species",
                    paste0("Species found: ", length(roundResCol)),
                    "PAE-PCE analysis requires at least 2 synapomorphic species"),
        iteration = contagem
      ))
    } else if(length(roundResCol) < 2 && contagem > 1){
      print(paste0(c('The iteration number ', contagem, ' has finished without results.'), collapse = ''))
      print('This analysis had to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      break
    }

    quantum <- 1
    for(j in 1:length(roundResCol)){
      k <- roundResCol[quantum]
      discre <- matTemp2[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species

      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, 'is outside of map! Please check it!'))
        next
      }

      discre <- as.factor(discre)  ### species presence vector

      lista_TT[[k]] <- names(discre)[which(discre == k)]

      quantum <- quantum + 1
    }

    subtrair <- NULL
    for(i in colnames(matTemp2)){
      subtrair[i] <- which(colnames(preabsMat) == i)
    }

    print(paste0("Removing ", length(subtrair), " species (characters) from matrix for next iteration"))
    print(paste0("Species (characters) to remove: ", paste(names(subtrair), collapse = ", ")))

    # Changing the matrix for the next iteration:
    tempMatrix <- subset(tempMatrix, select=-subtrair, site.pattern = FALSE)

    print(paste0("New tempMatrix: ", dim(as.data.frame(tempMatrix))[1], " characters and ",
                 dim(as.data.frame(tempMatrix))[2], " terminals"))

    # CRUCIAL VERIFICATION: If there are not enough species, stop
    if(dim(as.data.frame(tempMatrix))[1] < 3) {
      print("Less than 3 species remaining in matrix - stopping analysis")
      deve_continuar <- FALSE
      break
    }

    speciesNamesT <- names(lista_TT)

    qde <- 0
    idx <- 1
    for(i in speciesNamesT){
      idx <- idx + 1
      qde[idx] <- length(unlist(lista_TT[speciesNamesT[idx - 1]]))
    }

    qdeCum <- cumsum(qde)

    tabelaT <- matrix(NA, nrow = qdeCum[length(qde)], nc = 2)

    indice <- 0
    for(j in speciesNamesT){
      indice <- indice + 1
      tabelaT[(qdeCum[indice] + 1): qdeCum[indice + 1], 1] <- as.character(rep(j, qde[indice + 1]))
      if(is.numeric(unlist(lista_TT[indice])) == TRUE){
        tabelaT[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.numeric(unlist(lista_TT[[indice]]))
      } else {
        tabelaT[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.character(unlist(lista_TT[[indice]]))
      }
    }

    if(dim(tabelaT)[1] == 0){
      # ELEGANT STOP: Empty results table
      return(finalizar_analise_elegante(
        reason = "No results generated in the iteration",
        details = c("Species and grids table resulted empty",
                    "May indicate natural endpoint of PAE-PCE analysis"),
        iteration = contagem
      ))
    }
    # Data frame with the results...
    # synapomorphies
    frameTempT <- data.frame(spp = tabelaT[, 1], grid_n = tabelaT[, 2], row.names = 1:dim(tabelaT)[1])

    # Matrix with species names and grid numbers:
    resulPaeRasterT <- matrix(as.matrix(frameTempT[,2]), nrow(frameTempT),
                              1, dimnames = list(frameTempT[,1], colnames(frameTempT)[2]))


    ## The number of species supporting the grid number
    speciesNumberT <- data.frame(table(resulPaeRasterT))

    n_occurT <- data.frame(table(unlist(lista_TT))) # frequency of each grid

    similarNames <- list()
    # Get unique grid values from speciesNumberT
    unique_grids <- as.character(speciesNumberT$resulPaeRasterT)

    for(volta in 1:length(unique_grids)){
      grid_value <- unique_grids[volta]
      # Find all species that have this grid value
      matching_rows <- which(as.character(resulPaeRasterT[,1]) == grid_value)

      if(length(matching_rows) > 0) {
        similarNames[[volta]] <- resulPaeRasterT[matching_rows, ]
        # Get species names for this grid
        species_names <- subset(frameTempT, as.character(grid_n) == grid_value)[,1]
        names(similarNames[[volta]]) <- species_names
      }
    }

    ## Grid numbers that are generalized tracks in this iteration
    if(is.numeric(unlist(similarNames)) == TRUE){
      gridIt <- unique(as.numeric(unlist(similarNames)))
    } else {
      gridIt <- unique(as.character(unlist(similarNames)))
    }


    if(length(unique(rownames(resulPaeRasterT))) == 0){
      # ELEGANT STOP: No unique species in result
      return(finalizar_analise_elegante(
        reason = "No unique species identified in result",
        details = c("Iteration result contains no distinct species",
                    "May indicate end of analysis or data problems"),
        iteration = contagem
      ))
    }

    if(length(unique(frameTempT$spp)) > 1){

      if(use_grid_mapping) {
        # ===== IRREGULAR BINS MODE: Plot grids associated with endemic areas =====
        cat("\nPlotting grids associated with irregular bins...\n")

        # Get area names with frequency > 1 (generalized tracks)
        endemic_areas <- c()
        for(j in 1:dim(speciesNumberT)[1]){
          if(speciesNumberT$Freq[j] > 1){
            area_name <- as.character(speciesNumberT[j,1])
            endemic_areas <- c(endemic_areas, area_name)
          }
        }

        cat(paste0("Endemic areas identified: ", paste(endemic_areas, collapse = ", "), "\n"))

        # Find all grid cells that belong to endemic areas
        endemic_grid_ids <- grid_to_area_mapping$grid_id[
          grid_to_area_mapping$area_name %in% endemic_areas
        ]

        cat(paste0("Grid cells to plot: ", length(endemic_grid_ids), "\n"))

        if(length(endemic_grid_ids) > 0) {
          # Get the polygons of endemic grids
          endemic_grids_sf <- gridPolygon_sf[gridPolygon_sf$grid_id %in% endemic_grid_ids, ]

          # Convert to Spatial for plotting
          endemic_grids_sp <- as(endemic_grids_sf, "Spatial")

          # Plot ALL grid polygons directly with hatching
          if(sobrepo == TRUE){
            plot(endemic_grids_sp, add = TRUE,
                 col = adjustcolor(colores[contagem], alpha.f = 0.3),
                 border = colores[contagem], lwd = 1.5,
                 density = 15, angle = 45)  # Hatching pattern
          } else {
            plot(endemic_grids_sp, add = TRUE,
                 col = adjustcolor(colores[1], alpha.f = 0.3),
                 border = colores[1], lwd = 1.5,
                 density = 15, angle = 45)  # Hatching pattern
          }

          cat(paste0("Plotted ", nrow(endemic_grids_sf), " grid cells (all cells of endemic areas)\n"))
        } else {
          cat("No endemic areas to plot\n")
        }

        # Save raster output
        if(length(endemic_grid_ids) > 0) {
          # Rasterize endemic grids for output file
          r[r > 0] <- NA
          r_endemic <- rasterize(endemic_grids_sp, r, field = 1)
          r_masked <- mask(r_endemic, shapeFile)

          map.r <- raster::as.data.frame(raster::rasterToPoints(r_masked))
          if(nrow(map.r) > 0) {
            pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r_masked, field = 1)
            proj4string(pontosRaster) <- CRS("+proj=longlat +datum=WGS84")
            writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
                        overwrite = TRUE)
          }
        }

      } else {
        # ===== REGULAR GRID MODE: Use raster as before =====
        # Set the cells associated with the shapefile to the specified value
        r[r > 0] <- NA

        for(j in 1:dim(speciesNumberT)[1]){
          if(speciesNumberT$Freq[j] > 1){
            grid_value <- as.character(speciesNumberT[j,1])
            gTrack <- as.numeric(grid_value)
            values(r)[gTrack] <- 1
          }
        }

        # CORRECTION: Cut plotted cells at shapefile boundaries
        r_masked <- mask(r, shapeFile)

        if(sobrepo == TRUE){
          plot(r_masked, axes = FALSE, legend = FALSE, add = TRUE, col = colores[contagem],
               alpha = 0.60)
        } else if(sobrepo == FALSE){
          plot(r_masked, axes = FALSE, legend = FALSE, add = TRUE, col = colores[1],
               alpha = 0.50)
        }

        # Producing a raster:
        # Convert the raster to points for plotting the number of a grid
        map.r <- raster::as.data.frame(raster::rasterToPoints(r_masked))
        pontosRaster <- rasterize(cbind(map.r$x, map.r$y), r_masked, field = 1) # raster with presences
        proj4string(pontosRaster) <- CRS("+proj=longlat +datum=WGS84")
        writeRaster(pontosRaster, paste0(c('out/generalizedTrack_', contagem, '.tif'), collapse = ''),
                    overwrite = TRUE)
      }

      # Grid view and labeling
      if(!use_grid_mapping) {
        ## Note: gridPolygon@data has all polygons established with a specific resolution

        if(gridView == TRUE){
          plot(cropped_map, add = TRUE, border = "gray", lwd = 0.5)

          map.r$gridNumber <- which(pontosRaster@data@values == 1)

          if(labelGrid == TRUE){
            text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8, col = 'black', font = 2)
          }
        }
      } else {
        # For irregular bins, show grid boundaries and optionally label
        if(gridView == TRUE){
          plot(cropped_map, add = TRUE, border = "gray", lwd = 0.5)

          if(labelGrid == TRUE){
            # Label grid cells with abbreviated area names
            if(exists("endemic_grid_ids") && length(endemic_grid_ids) > 0) {
              # Get centroids of endemic grids
              endemic_grids_sf <- gridPolygon_sf[gridPolygon_sf$grid_id %in% endemic_grid_ids, ]
              centroids <- st_centroid(st_geometry(endemic_grids_sf))
              coords <- st_coordinates(centroids)

              # Get area names for these grids
              area_labels_full <- grid_to_area_mapping$area_name[
                match(endemic_grid_ids, grid_to_area_mapping$grid_id)
              ]

              # Create abbreviations (2-3 letters)
              # Store unique area names and their abbreviations
              if(!exists("area_abbreviations")) {
                unique_areas <- unique(area_labels_full)
                area_abbreviations <<- data.frame(
                  full_name = unique_areas,
                  abbreviation = sapply(unique_areas, function(name) {
                    # Split by spaces and take first letters
                    words <- strsplit(name, " ")[[1]]
                    if(length(words) >= 2) {
                      # Use first letter of each word (max 3)
                      abbr <- paste0(substr(words[1:min(3, length(words))], 1, 1), collapse = "")
                    } else {
                      # Use first 3 letters of single word
                      abbr <- substr(name, 1, 3)
                    }
                    return(toupper(abbr))
                  }),
                  stringsAsFactors = FALSE
                )

                # Print legend to console
                cat("\n=== AREA ABBREVIATIONS LEGEND ===\n")
                for(i in 1:nrow(area_abbreviations)) {
                  cat(paste0(area_abbreviations$abbreviation[i], " = ",
                             area_abbreviations$full_name[i], "\n"))
                }
                cat("=================================\n\n")
              }

              # Get abbreviations for labels
              area_labels_abbr <- area_abbreviations$abbreviation[
                match(area_labels_full, area_abbreviations$full_name)
              ]

              text(coords[,1], coords[,2],
                   labels = area_labels_abbr,
                   cex = 0.7, col = 'black', font = 2)
            }
          }
        }
      }

      if(nonHomoplasticSpeciesList == TRUE){
        print('Please, check the nonHomoplasticSpeciesList.txt file in the out/ directory...')
        if (!file.exists('out/'))
          dir.create('out/')
        logfile <- "out/nonHomoplasticSpeciesList.txt"
        cat(c("grid_identifier", "species", "\n"), file = logfile, sep="\t")

        # Iterate over unique grid values
        for(i in unique_grids){
          grid_label <- if(use_grid_mapping) {
            # Use area name for string grids
            i
          } else {
            # Use numeric ID
            paste0('grid_', i)
          }

          # Find species with this grid value
          species_in_grid <- rownames(subset(resulPaeRasterT, as.character(resulPaeRasterT[,'grid_n']) == i))

          if(length(species_in_grid) > 0) {
            cat(c(grid_label, species_in_grid, '\n'), sep = '\t', file = logfile, append = TRUE)
          }
        }
      }
    }

    # ##################################
    # ## Changing matTemp:
    ######################################################
    matTemp <- preabsMat[, -subtrair] # temporary matrix with homoplastic synapomorphies
    matTemp <- as.matrix(matTemp)
    #######################################################


    nomesCOL <- colnames(matTemp) # species which are not synapomorphies

    print(paste0("Remaining species (non-synapomorphic): ", length(nomesCOL)))
    if(length(nomesCOL) > 0) {
      print(paste0("Species: ", paste(nomesCOL, collapse = ", ")))
    }

    quantum <- 1
    for(j in 1:length(nomesCOL)){
      k <- nomesCOL[quantum]
      discre <- matTemp[, as.character(k)]
      for(a in 1:length(discre)){
        ifelse(discre[a] == 1, discre[a] <- k, discre[a] <- 'absent')
      } # vector with the presence and absence of one species

      if(all(discre == 'absent')){
        quantum <- quantum + 1
        print(paste0(k, ' is outside of map! This species is absent from the given set of grids!'))
        next
      }

      discre <- as.factor(discre)  ### species presence vector

      lista_T[[k]] <- names(discre)[which(discre == k)]

      quantum <- quantum + 1
    }

    speciesNames <- names(lista_T)

    qde <- 0
    idx <- 1
    for(i in speciesNames){
      idx <- idx + 1
      qde[idx] <- length(unlist(lista_T[speciesNames[idx - 1]]))
    }

    qdeCum <- cumsum(qde)

    tabela <- matrix(NA, nrow = qdeCum[length(qde)], nc = 2)

    indice <- 0
    for(j in speciesNames){
      indice <- indice + 1
      tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 1] <- as.character(rep(j, qde[indice + 1]))
      if(is.numeric(unlist(lista_T[[indice]])) == TRUE){
        tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.numeric(unlist(lista_T[[indice]]))
      } else {
        tabela[(qdeCum[indice] + 1): qdeCum[indice + 1], 2] <- as.character(unlist(lista_T[[indice]]))
      }
    }

    if(dim(tabela)[1] == 0){
      print('No more iterations are needed.')
      deve_continuar <- FALSE
      break
    }

    # Data frame with the results... Homoplastic species
    frameTemp <- data.frame(spp = tabela[, 1], grid_n = tabela[, 2], row.names = 1:dim(tabela)[1])

    # Matrix with species names and grid numbers:
    resulPaeRaster <- matrix(as.matrix(frameTemp[,2]), nrow(frameTemp),
                             1, dimnames = list(frameTemp[,1], colnames(frameTemp)[2]))


    ## The number of species supporting the grid number
    speciesNumber <- data.frame(table(resulPaeRaster))

    n_occur <- data.frame(table(unlist(lista_T))) # frequency of each grid


    syn_grids <- unlist(lista_T)[unlist(lista_T) %in% n_occur$Var1[n_occur$Freq > 1]]
    if(length(syn_grids) < 2){
      print('This iteration has to be stopped because there is no more synapomorphies!')
      deve_continuar <- FALSE
      next
    }

    # ESSENTIAL: Save species lists for each iteration
    lista[[contagem]] <- frameTemp  # species which are not synapomorphies (homoplastic_species)
    listaR[[contagem]] <- frameTempT # species which are synapomorphies (nonHomoplastic_species)


    preabsMat <- matTemp

    print(paste0("New preabsMat matrix: ", dim(preabsMat)[2], " species, ", dim(preabsMat)[1], " grids"))

    # FINAL VERIFICATION: If there are not enough species to continue
    if(dim(preabsMat)[2] < 3) {
      print("Less than 3 species remaining - stopping analysis")
      deve_continuar <- FALSE
      break
    }

    print(paste0(c('The iteration number ', contagem, ' has finished.'), collapse = ''))

    contagem <- contagem + 1
    conta[contagem] <- contagem

    print(paste0("=== END OF ITERATION ", contagem - 1, " ==="))
    print("")
  } # close while looping

  print("=== WHILE LOOP FINALIZED ===")
  print(paste0("Reason for exit: deve_continuar = ", deve_continuar))
  print(paste0("Characters with CI=1: ", length(which((ciVec == 1) == TRUE))))
  print(paste0("Columns in preabsMat: ", dim(preabsMat)[2]))

  # ELEGANT TERMINATIONS FOR END OF LOOP
  if (exists("matTemp") == FALSE && length(conta) == 1){
    # ELEGANT STOP: Temporary matrix does not exist
    return(finalizar_analise_elegante(
      reason = "Temporary matrix was not created",
      details = c("Variable matTemp does not exist after first iteration",
                  "May indicate failure in analysis initialization"),
      iteration = length(conta)
    ))
  }

  if(length(colnames(matTemp)) == 0 && length(conta) == 1){
    # ELEGANT STOP: Empty matrix in first iteration
    return(finalizar_analise_elegante(
      reason = "Data matrix empty after first iteration",
      details = c("No columns remaining in matrix after first iteration",
                  "All characters were removed or there is no valid data"),
      iteration = length(conta)
    ))
  } else if(length(colnames(matTemp)) == 0 && length(conta) > 1){
    xis <- seq(1:length(listaR))
  }

  if(is.null(nomesCOL) || dim(as.data.frame(tempMatrix))[1] == 0){
    xis <- seq(1:length(listaR))
  }


  # Only label grids if there are results
  if(labelGrid == TRUE && length(listaR) > 0 && exists("map.r")){
    if(sobrepo == TRUE){
      print('Bold numbers refer to PAE analysis using generalized tracks!')
    }
    text(map.r[,c(1, 2)], labels = map.r$gridNumber, cex = 0.8, col = 'black', font = 2)
  }

  conta <- NULL
  # Only show legend if there are results
  if(length(listaR) > 0) {
    if(sobrepo == FALSE){
      xis <- seq(1:length(listaR))
      legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
             title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
    } else if(sobrepo == TRUE){
      xis <- seq(1:length(listaR))
      legend(x = 'topright', legend = xis, pch = 15, col = colores[xis],
             title = 'Adding generalized track of the iterations...', title.col = 'red', pt.cex = 1.5, cex = 0.8)
    }
  }

  if(contagem == 1 && (length(listaR) == 0 || is.null(listaR[[1]]))){
    # ELEGANT TERMINATION: No results
    return(finalizar_analise_elegante(
      reason = "Analysis completed without results",
      details = c("No iteration produced valid results",
                  "Data may not be suitable for PAE-PCE analysis"),
      iteration = contagem,
      plot_grid = FALSE
    ))
  } else if(contagem == 2 && length(listaR) != 0){
    # NORMAL TERMINATION: Successful analysis
    print("Only one iteration!")
    print("")
    print(paste(rep("=", 60), collapse=""))
    print("        PAE-PCE ANALYSIS COMPLETED SUCCESSFULLY")
    print(paste(rep("=", 60), collapse=""))
    print(paste0("Total iterations: ", contagem))
    print(paste(rep("=", 60), collapse=""))

    return(list(homoplastic_species = lista,
                nonHomoplastic_species = listaR))
  } else if(contagem > 2){
    # NORMAL TERMINATION: Successful analysis
    print("")
    print("")
    print(paste(rep("=", 60), collapse=""))
    print("        PAE-PCE ANALYSIS COMPLETED SUCCESSFULLY")
    print(paste(rep("=", 60), collapse=""))
    print(paste0("Total iterations: ", contagem - 1))
    print(paste(rep("=", 60), collapse=""))

    return(list(homoplastic_species = lista,
                nonHomoplastic_species = listaR))
  }
} # close if... else... condition
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
#' Range Extrapolation Integration Functions
#'
#' Wrapper functions to integrate actual range extrapolation methods into the Shiny app
#'

#' Load Phylogenetic Tree
#'
#' Load a phylogenetic tree from file (Newick or Nexus format)
#'
#' @param tree_file Path to tree file
#'
#' @return phylo object
#'
#' @export
load_phylogenetic_tree <- function(tree_file) {
  
  tryCatch({
    # Try reading as Newick first
    tree <- tryCatch({
      ape::read.tree(tree_file)
    }, error = function(e) {
      # Try reading as Nexus
      ape::read.nexus(tree_file)
    })
    
    # Check if tree has branch lengths
    if (is.null(tree$edge.length)) {
      warning("Tree does not have branch lengths")
    }
    
    tree
  }, error = function(e) {
    stop("Error reading tree file: ", e$message)
  })
}

#' Validate Phylogenetic Tree
#'
#' Check if tree is valid
#'
#' @param tree phylo object
#'
#' @return List with validation results
#'
#' @export
validate_phylogenetic_tree <- function(tree) {
  
  list(
    is_valid = inherits(tree, "phylo"),
    n_taxa = length(tree$tip.label),
    has_branch_lengths = !is.null(tree$edge.length),
    taxa = tree$tip.label
  )
}

#' Prepare Occurrence Data
#'
#' Standardize column names for range functions
#'
#' @param occurrence_data Data frame with species, longitude, latitude
#'
#' @return Data frame with standardized columns
#'
#' @export
prepare_occurrence_data <- function(occurrence_data) {
  
  # Standardize column names
  names(occurrence_data) <- tolower(names(occurrence_data))
  
  # Rename to match function requirements
  if ("species" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "species"] <- "spp"
  }
  if ("longitude" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "longitude"] <- "long"
  }
  if ("latitude" %in% names(occurrence_data)) {
    names(occurrence_data)[names(occurrence_data) == "latitude"] <- "lat"
  }
  
  # Return only required columns
  occurrence_data[, c("spp", "long", "lat")]
}

#' Extract Presence-Absence Matrix
#'
#' Extract pres_abs from range function results
#'
#' @param result Result from range extrapolation function
#'
#' @return Presence-absence matrix
#'
#' @export
extract_pres_abs <- function(result) {
  
  if (is.list(result) && "pres_abs" %in% names(result)) {
    return(result$pres_abs)
  }
  
  if (is.matrix(result)) {
    return(result)
  }
  
  stop("Cannot extract presence-absence matrix")
}

#' Extract Geometry from Range Results
#'
#' Extract spatial geometry (polygons) from range results
#'
#' @param result Result from range extrapolation function
#'
#' @return Spatial geometry object or NULL
#'
#' @export
extract_geometry <- function(result) {
  
  if (is.list(result) && "geometry" %in% names(result)) {
    return(result$geometry)
  }
  
  NULL
}

#' Create BioGeoBEARS Geography File
#'
#' Convert presence-absence matrix to BioGeoBEARS .data format
#'
#' @param pres_abs Presence-absence matrix
#' @param method Method name ("BUFF", "MPC", "MST")
#' @param output_dir Output directory
#'
#' @return Path to created file
#'
#' @export
create_biogeobears_file <- function(pres_abs, method = "BUFF", output_dir = ".") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Remove empty ranges (all zeros)
  pres_abs_clean <- pres_abs[rowSums(pres_abs) > 0, ]
  
  # Create output file
  output_file <- file.path(output_dir, paste0("pres_abs_", method, "_geog.data"))
  
  # Write header
  n_species <- ncol(pres_abs_clean)
  n_ranges <- nrow(pres_abs_clean)
  range_names <- paste(1:n_ranges, collapse = " ")
  
  # Convert matrix to binary strings
  binary_strings <- apply(pres_abs_clean, 1, function(row) paste(row, collapse = ""))
  
  # Write file
  writeLines(c(
    paste(n_species, n_ranges, range_names),
    paste(colnames(pres_abs_clean), binary_strings, sep = " ")
  ), output_file)
  
  output_file
}

#' Convert Geometry to SF for Leaflet
#'
#' Convert raster or terra geometry to sf polygons for Leaflet visualization
#'
#' @param geometry Spatial geometry object (raster, SpatRaster, or sf)
#'
#' @return sf object with polygons
#'
#' @export
convert_geometry_to_sf <- function(geometry) {
  
  if (is.null(geometry)) {
    return(NULL)
  }
  
  tryCatch({
    # If already sf, return as is
    if (inherits(geometry, "sf")) {
      return(geometry)
    }
    
    # If terra SpatRaster, convert to polygons
    if (inherits(geometry, "SpatRaster")) {
      # Convert raster to polygons
      poly <- terra::as.polygons(geometry)
      # Convert to sf
      sf_poly <- sf::st_as_sf(poly)
      return(sf_poly)
    }
    
    # If terra SpatVector, convert to sf
    if (inherits(geometry, "SpatVector")) {
      sf_poly <- sf::st_as_sf(geometry)
      return(sf_poly)
    }
    
    # If raster package raster, convert
    if (inherits(geometry, "RasterLayer") || inherits(geometry, "RasterStack")) {
      # Convert to polygons
      poly <- raster::rasterToPolygons(geometry)
      # Convert to sf
      sf_poly <- sf::st_as_sf(poly)
      return(sf_poly)
    }
    
    # If sp object, convert to sf
    if (inherits(geometry, "Spatial")) {
      sf_poly <- sf::st_as_sf(geometry)
      return(sf_poly)
    }
    
    warning("Unsupported geometry type")
    NULL
  }, error = function(e) {
    warning("Error converting geometry: ", e$message)
    NULL
  })
}

#' Add Extrapolation Polygons to Leaflet Map
#'
#' Add extrapolation polygons to an existing leaflet map
#'
#' @param map Leaflet map object
#' @param geometry Spatial geometry object
#' @param method Method name ("BUFF", "MPC", "MST")
#' @param color Color for polygons
#' @param opacity Opacity of polygons
#'
#' @return Updated leaflet map
#'
#' @export
add_extrapolation_polygons <- function(map, geometry, method = "BUFF", 
                                        color = "blue", opacity = 0.5) {
  
  if (is.null(geometry)) {
    return(map)
  }
  
  tryCatch({
    # Convert to sf if needed
    sf_geom <- convert_geometry_to_sf(geometry)
    
    if (is.null(sf_geom)) {
      return(map)
    }
    
    # Add polygons to map
    map <- map %>%
      leaflet::addPolygons(
        data = sf_geom,
        color = color,
        weight = 2,
        opacity = opacity,
        fillOpacity = opacity * 0.5,
        popup = paste(method, "extrapolation"),
        group = method
      )
    
    return(map)
  }, error = function(e) {
    warning("Error adding polygons: ", e$message)
    return(map)
  })
}
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
#' Run the BioGeoBEARS Shiny Application
#'
#' Launches the BioGeoBEARS Shiny application for ancestral range estimation
#' and biogeographic analysis.
#'
#' @param ... Additional arguments passed to [shiny::shinyApp()]
#'
#' @return Invisibly returns the Shiny application object
#'
#' @examples
#' \dontrun{
#' run_biogeoshiny()
#' }
#'
#' @export
run_biogeoshiny <- function(...) {
  app_dir <- system.file("shiny", package = "biogeoshiny")
  if (app_dir == "") {
    stop("Could not find shiny app directory. Try reinstalling biogeoshiny.",
         call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal", ...)
}
#' singleton_to_data_frame
#'
#' Removes singletons from a dataset and prunes taxa from a phylogeny.
#'
#' @param data coordinates - in decimal degrees - with the following header: \strong{spp, lat,
#'  long (all lowercase)}
#' @param phylogeny a phylogenetic tree in Newick format
#' @return Both a new dataset as a data frame object and a modified tree without singletons.
#' @author Jose Ricardo Inacio Ribeiro \email{joseribeiro@@unipampa.edu.br}
#'
#' Augusto Ferrari \email{ferrariaugusto@@gmail.com}
#' @details It can be utilized even if the dataset does not contain any singletons, which may only be discovered later.
#'  It is also applicable in situations where a taxon present in the tree is absent from the dataset, but it is also applicable
#'    in situations where a taxon is absent from the tree and present in the dataset. Finally, it can be used without a phylogeny.
#'      The idea was based on the phytools package (Revell, 2024).
#' @references Revell, L.J., 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other
#'  things). \emph{PeerJ} \strong{12}: e16505.
#' @seealso
#' \code{\link[phytools]{drop.tip}} for pruning a taxon out of phylogenetic trees.
#'
#' \code{\link[PanBioGeo]{MST_node}} for estimating the minimum spanning tree (MST).
#'
#' @export
#' @examples
#' # load example data (coordinates and tree)
#' library(phytools)
#' data(lycipta)
#'
#' # The columns are in the following order: species (spp), longitude (long), latitude (lat)
#' lycipta.coords <- lycipta$coordinates
#' head(lycipta.coords, 20)
#'
#' # phylogenetic tree
#' lycipta.tree <- lycipta$phylogenetic_tree
#'
#' # removal of singletons with a tree
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)], phylogeny = lycipta.tree)
#'
#' # removal of singletons without a tree
#' resul <- singleton_to_data_frame(data = lycipta.coords[, c(1, 3, 2)])

singleton_to_data_frame <- function(data = NULL, phylogeny = NULL){
  # Input validation
  if(is.null(data)){
    stop("The 'data' argument cannot be NULL. Provide a data.frame with columns: species, latitude and longitude.")
  }
  spp <- as.character(data[, 1])
  lat <- data[, 2]
  long <- data[, 3]

  names(lat) <- spp
  names(long) <- spp
  locais <- cbind(lat, long)
  rownames(locais) <- spp

  spp_counts <- table(spp)
  singleton_taxa <- names(spp_counts[spp_counts == 1])

  build_data_df <- function(locais_obj) {
    data.frame(
      spp = rownames(locais_obj),
      long = locais_obj[, 2],
      lat = locais_obj[, 1],
      stringsAsFactors = FALSE
    )
  }

  # without a phylogeny:
  if(is.null(phylogeny)){
    if(length(singleton_taxa) == 0){
      data_df <- build_data_df(locais)
      message("The dataset is OK. No modification was necessary!")
      return(list(data_df = data_df))
    } else {
      locais_new <- locais[!(rownames(locais) %in% singleton_taxa), , drop = FALSE]
      data_df <- build_data_df(locais_new)
      message(paste("Removed", length(singleton_taxa), "species with only one occurrence."))
      return(list(data_df = data_df))
    }
  } else if(!is.null(phylogeny)){ # with a phylogeny
    # Tree validation
    if(class(phylogeny)[1] != "phylo"){
      stop("The 'phylogeny' argument must be an object of class 'phylo' from the 'ape' package.")
    }

    locais_new <- locais
    tree_new <- phylogeny
    n_singleton_removed <- 0L
    n_data_not_tree_removed <- 0L
    n_tree_not_data_removed <- 0L

    if (length(singleton_taxa) > 0) {
      locais_new <- locais_new[!(rownames(locais_new) %in% singleton_taxa), , drop = FALSE]
      drop_singletons <- intersect(singleton_taxa, tree_new$tip.label)
      if (length(drop_singletons) > 0) {
        tree_new <- ape::drop.tip(tree_new, drop_singletons)
      }
      n_singleton_removed <- length(singleton_taxa)
    }

    if (nrow(locais_new) == 0) {
      stop("All taxa were removed during singleton filtering. Please review the occurrence data.")
    }

    data_taxa <- unique(rownames(locais_new))
    tree_taxa <- tree_new$tip.label

    taxa_data_not_tree <- setdiff(data_taxa, tree_taxa)
    taxa_tree_not_data <- setdiff(tree_taxa, data_taxa)

    if (length(taxa_data_not_tree) > 0) {
      locais_new <- locais_new[!(rownames(locais_new) %in% taxa_data_not_tree), , drop = FALSE]
      n_data_not_tree_removed <- length(taxa_data_not_tree)
    }
    if (length(taxa_tree_not_data) > 0) {
      tree_new <- ape::drop.tip(tree_new, taxa_tree_not_data)
      n_tree_not_data_removed <- length(taxa_tree_not_data)
    }

    data_df <- build_data_df(locais_new)
    tree_modified <- n_tree_not_data_removed > 0 || length(intersect(singleton_taxa, phylogeny$tip.label)) > 0

    if (length(singleton_taxa) == 0 && !tree_modified) {
      message("The dataset is OK. No modification was necessary!")
      return(list(data_df = data_df, treeNonMod = tree_new))
    } else {
      message(paste0(
        "Singleton/tree harmonization complete. ",
        "Removed singleton taxa: ", n_singleton_removed, "; ",
        "removed data-only taxa: ", n_data_not_tree_removed, "; ",
        "removed tree-only taxa: ", n_tree_not_data_removed, "."
      ))
      return(list(data_df = data_df, treeMod = tree_new))
    }
  }
}
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
#' toNDM
#'
#' This function converts species data into the XYD format required by NDM (VNDM)
#' software for identifying areas of endemism. The function accepts three types of input:
#' (1) occurrence data from CSV/data frames/matrices, (2) shapefiles from extrapolated
#' distributions (buffers, convex hulls, MST), or (3) a directory containing multiple
#' species shapefiles. For shapefiles, the function creates a grid with specified
#' resolution and uses grid cell centroids as occurrence points.
#'
#' @param input_data Either:
#'   \itemize{
#'     \item A character string with path to CSV file (occurrence data)
#'     \item A data frame with 3 columns: species, longitude, latitude
#'     \item A matrix/table with 3 columns: species, longitude, latitude
#'     \item A character string with path to shapefile directory (e.g., "out_buffers/")
#'     \item A character string with "SHAPEFILE" to trigger shapefile mode
#'   }
#' @param shapefile_dir Character string specifying directory containing species
#'   shapefiles (e.g., "out_buffers/", "out_MPC/", "out/"). Required when using
#'   shapefile mode. Each shapefile should be named with species identifier.
#' @param resolution Numeric vector of length 2 specifying grid resolution in degrees
#'   (e.g., c(5, 5) for 5x5 degree cells). Required for shapefile mode. This should
#'   match the resolution used in the extrapolation analysis.
#' @param shape_file Optional spatial object (sf, SpatVector, or Spatial*) representing
#'   the study region boundary. Used to crop the grid in shapefile mode.
#' @param output_file Character string specifying output XYD file name
#'   (default: "output_NDM.xyd").
#' @param separator Character string specifying field separator for CSV files (default: ",").
#' @param header Logical, whether CSV has column names (default: TRUE).
#' @param na_strings Character vector of strings to interpret as NA (default: c("?", "-", "NA", "")).
#' @param remove_na Logical, if TRUE removes rows with NA values (default: TRUE).
#' @param output_dir Character string specifying output directory (default: "out_NDM/").
#'
#' @return The function does not return a value. It creates an XYD-formatted file
#'   ready for NDM/VNDM analysis. Invisibly returns a list with conversion statistics.
#'
#' @details
#' The function operates in two modes:
#'
#' **Mode 1: Occurrence Data**
#' \itemize{
#'   \item Reads occurrence coordinates from CSV, data frame, or matrix
#'   \item Removes NA values and empty species
#'   \item Formats for NDM
#' }
#'
#' **Mode 2: Shapefile (Extrapolated Distributions)**
#' \itemize{
#'   \item Reads species shapefiles from specified directory
#'   \item Creates grid with specified resolution
#'   \item Identifies grid cells intersecting each species shapefile
#'   \item Calculates centroid coordinates for each grid cell
#'   \item Uses centroids as "occurrence" points for NDM
#'   \item Exports to XYD format
#' }
#'
#' @section Shapefile Mode Workflow:
#' When using extrapolated distributions (buffers, convex hulls, MST):
#' \enumerate{
#'   \item Run extrapolation analysis (e.g., \code{calcRange_buffers()})
#'   \item Shapefiles are saved in output directory (e.g., "out_buffers/")
#'   \item Use \code{toNDM()} with shapefile_dir and resolution parameters
#'   \item Function reads all shapefiles, creates grid, extracts centroids
#'   \item Generates XYD file with centroid coordinates
#'   \item Analyze in NDM using grid-based "occurrences"
#' }
#'
#' @section Input Formats:
#'
#' **Format 1: Occurrence Data (CSV/Data Frame)**
#' \preformatted{
#' species,longitude,latitude
#' Belostoma_amazonum,-60.5,-3.2
#' Belostoma_angustum,-58.7,-5.3
#' }
#'
#' **Format 2: Shapefiles in Directory**
#' \preformatted{
#' out_buffers/
#'    Belostoma_amazonum.shp
#'    Belostoma_angustum.shp
#'    Belostoma_anurum.shp
#' }
#'
#' @note
#' \itemize{
#'   \item For occurrence data: input must have 3 columns (species, lon, lat)
#'   \item For shapefiles: resolution parameter is required
#'   \item Shapefile names should contain species identifiers
#'   \item Grid resolution should match the one used in extrapolation analysis
#'   \item Centroids are calculated for grid cells intersecting species ranges
#'   \item Output directory is created automatically if needed
#' }
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # MODE 1: From occurrence data (original functionality)
#' # ============================================================================
#'
#' # Example 1: From CSV
#' toNDM(input_data = "occurrences.csv")
#'
#' # Example 2: From data frame
#' occurrences <- data.frame(
#'   species = c("Species_A", "Species_B"),
#'   longitude = c(-60.5, -58.7),
#'   latitude = c(-3.2, -5.3)
#' )
#' toNDM(input_data = occurrences)
#'
#' # ============================================================================
#' # MODE 2: From extrapolated distributions (NEW!)
#' # ============================================================================
#'
#' # Example 3: From buffer shapefiles
#' # Step 1: Run buffer analysis
#' calcRange_buffers(
#'   xy = lycipta_final,
#'   shape_file = neo,
#'   resol = 10,
#'   buffer.width = 500000
#' )
#' # This creates shapefiles in out_buffers/
#'
#' # Step 2: Convert shapefiles to NDM format
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(10, 10),
#'   shape_file = neo,
#'   output_file = "buffer_extrapolated_NDM.xyd"
#' )
#'
#' # Example 4: From convex hull shapefiles
#' # Step 1: Run convex hull analysis
#' calcRange_convexHull(
#'   xy = species_coords,
#'   shape_file = south_america,
#'   resol = 5
#' )
#' # This creates shapefiles in out_MPC/
#'
#' # Step 2: Convert to NDM
#' toNDM(
#'   shapefile_dir = "out_MPC/",
#'   resolution = c(5, 5),
#'   shape_file = south_america,
#'   output_file = "convexhull_extrapolated_NDM.xyd"
#' )
#'
#' # Example 5: From MST shapefiles
#' # Step 1: Run MST analysis
#' calcRange_mst(
#'   xy = belostomatidae_coords,
#'   shape_file = neotropics,
#'   resol = 10
#' )
#' # This creates shapefiles in out/
#'
#' # Step 2: Convert to NDM
#' toNDM(
#'   shapefile_dir = "out/",
#'   resolution = c(10, 10),
#'   shape_file = neotropics,
#'   output_file = "mst_extrapolated_NDM.xyd"
#' )
#'
#' # ============================================================================
#' # Example 6: Complete workflow - Buffer to NDM
#' # ============================================================================
#'
#' library(sf)
#'
#' # Step 1: Load occurrence data
#' occurrences <- read.csv("species_occurrences.csv")
#'
#' # Step 2: Load study area shapefile
#' study_area <- st_read("study_area.shp")
#'
#' # Step 3: Calculate buffer ranges
#' calcRange_buffers(
#'   xy = occurrences,
#'   shape_file = study_area,
#'   resol = 5,
#'   buffer.width = 300000  # 300 km
#' )
#'
#' # Step 4: Convert extrapolated ranges to NDM format
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = study_area,
#'   output_file = "buffer_300km_NDM.xyd"
#' )
#'
#' # Step 5: Analyze in NDM/VNDM
#' # - Open VNDM
#' # - Load out_NDM/buffer_300km_NDM.xyd
#' # - Set parameters
#' # - Run analysis
#'
#' # ============================================================================
#' # Example 7: Comparing raw occurrences vs extrapolated ranges
#' # ============================================================================
#'
#' # Analysis 1: Raw occurrences
#' toNDM(
#'   input_data = "occurrences.csv",
#'   output_file = "raw_occurrences_NDM.xyd",
#'   output_dir = "NDM_comparison/raw/"
#' )
#'
#' # Analysis 2: Buffer extrapolation (100 km)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 5, buffer.width = 100000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "buffer_100km_NDM.xyd",
#'   output_dir = "NDM_comparison/buffer_100/"
#' )
#'
#' # Analysis 3: Buffer extrapolation (500 km)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 5, buffer.width = 500000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "buffer_500km_NDM.xyd",
#'   output_dir = "NDM_comparison/buffer_500/"
#' )
#'
#' # Analysis 4: Convex hull extrapolation
#' calcRange_convexHull(xy = coords, shape_file = shape, resol = 5)
#' toNDM(
#'   shapefile_dir = "out_MPC/",
#'   resolution = c(5, 5),
#'   shape_file = shape,
#'   output_file = "convexhull_NDM.xyd",
#'   output_dir = "NDM_comparison/hull/"
#' )
#'
#' # Compare all 4 analyses in NDM to assess method sensitivity
#'
#' # ============================================================================
#' # Example 8: Different grid resolutions
#' # ============================================================================
#'
#' # Coarse resolution (10 degrees)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 10, buffer.width = 300000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(10, 10),
#'   output_file = "buffer_10deg_NDM.xyd"
#' )
#'
#' # Fine resolution (2 degrees)
#' calcRange_buffers(xy = coords, shape_file = shape,
#'                   resol = 2, buffer.width = 300000)
#' toNDM(
#'   shapefile_dir = "out_buffers/",
#'   resolution = c(2, 2),
#'   output_file = "buffer_2deg_NDM.xyd"
#' )
#' }
#'
#' @seealso
#' \code{\link{calcRange_buffers}} for buffer-based range extrapolation
#' \code{\link{calcRange_convexHull}} for convex hull range extrapolation
#' \code{\link{calcRange_mst}} for MST-based range extrapolation
#'
#' @references
#' Szumik, C. A., & Goloboff, P. A. (2004). Areas of endemism: an improved
#' optimality criterion. Systematic Biology, 53(6), 968-977.
#'
#' Goloboff, P. A. (2004). NDM/VNDM, programs for identification of areas of
#' endemism. Program and documentation available at: www.lillo.org.ar/phylogeny/endemism
#'
#' @export
toNDM <- function(input_data = NULL,
                  shapefile_dir = NULL,
                  resolution = NULL,
                  shape_file = NULL,
                  output_file = "output_NDM.xyd",
                  separator = ",",
                  header = TRUE,
                  na_strings = c("?", "-", "NA", ""),
                  remove_na = TRUE,
                  output_dir = "out_NDM/") {

  ######################
  ## Input validation ##
  ######################

  # Check if at least one input is provided
  if (is.null(input_data) && is.null(shapefile_dir)) {
    stop("Error: Either 'input_data' or 'shapefile_dir' must be provided.")
  }

  ######################
  ## Detect mode #######
  ######################

  mode <- NULL

  if (!is.null(shapefile_dir)) {
    mode <- "shapefile"
    message("Mode: Shapefile (extrapolated distributions)")

    # Validate shapefile mode requirements
    if (is.null(resolution)) {
      stop("Error: 'resolution' parameter is required for shapefile mode.
           Example: resolution = c(5, 5)")
    }

    if (!dir.exists(shapefile_dir)) {
      stop(paste0("Error: Shapefile directory not found: ", shapefile_dir))
    }

  } else {
    mode <- "occurrences"
    message("Mode: Occurrence data")
  }

  ######################
  ## Process data ######
  ######################

  if (mode == "occurrences") {
    # MODE 1: Process occurrence data
    dados <- process_occurrences(input_data, separator, header, na_strings, remove_na)

  } else if (mode == "shapefile") {
    # MODE 2: Process shapefiles
    dados <- process_shapefiles(shapefile_dir, resolution, shape_file)
  }

  ######################
  ## Format for NDM ####
  ######################

  output <- format_for_ndm(dados)

  ######################
  ## Write output ######
  ######################

  write_ndm_file(output, output_file, output_dir)

  ######################
  ## Report results ####
  ######################

  report_conversion(dados, output_file, output_dir, mode)

  # Return summary invisibly
  invisible(list(
    output_file = file.path(output_dir, output_file),
    n_occurrences = nrow(dados),
    n_species = length(unique(dados[[1]])),
    mode = mode
  ))
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Process occurrence data
#' @keywords internal
process_occurrences <- function(input_data, separator, header, na_strings, remove_na) {

  if (is.null(input_data)) {
    stop("Error: 'input_data' is NULL.")
  }

  # Read data based on input type
  if (is.character(input_data)) {
    if (!file.exists(input_data)) {
      stop(paste0("Error: File not found: ", input_data))
    }
    message(paste0("Reading data from file: ", input_data))
    dados <- read.csv(
      file = input_data,
      sep = separator,
      header = header,
      na.strings = na_strings,
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(input_data)) {
    message("Using provided data frame")
    dados <- input_data
  } else if (is.matrix(input_data) || is.table(input_data)) {
    message("Converting matrix/table to data frame")
    dados <- as.data.frame(input_data, stringsAsFactors = FALSE)
  } else if(is.list(input_data)){ # Extract matrix from list if needed
    if("data_df" %in% names(input_data)){
      dados <- input_data$data_df
    } else {
      stop("Error: List does not contain 'data_df' component.")
    }
  } else {
    stop("Error: 'input_data' must be a file path, data frame, matrix, singleton_to_data_frame object, or table.")
  }

  # Validate 3 columns
  if (ncol(dados) != 3) {
    stop(paste0("Error: Data must have 3 columns (species, lon, lat). Found ", ncol(dados)))
  }

  # Convert to appropriate types
  dados[[1]] <- as.character(dados[[1]])
  dados[[2]] <- as.numeric(as.character(dados[[2]]))
  dados[[3]] <- as.numeric(as.character(dados[[3]]))

  # Remove NA
  if (remove_na) {
    dados <- na.omit(dados)
    rownames(dados) <- NULL
  }

  if (nrow(dados) == 0) {
    stop("Error: No valid data after removing NA values.")
  }

  # Sort by species
  dados <- dados[order(dados[[1]]), ]
  rownames(dados) <- NULL

  return(dados)
}


#' Process shapefiles from extrapolated distributions
#' @keywords internal
process_shapefiles <- function(shapefile_dir, resolution, shape_file) {

  # Load required packages
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for shapefile mode. Install with: install.packages('sf')")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for shapefile mode. Install with: install.packages('terra')")
  }

  library(sf)
  library(terra)

  message(paste0("Reading shapefiles from: ", shapefile_dir))

  # Find all shapefiles in directory
  all_shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)

  # Filter to keep only range shapefiles (BUFF_, MPC_, MST_)
  # Exclude point shapefiles
  shapefiles <- all_shapefiles[grepl("(BUFF_|MPC_|MST_)", basename(all_shapefiles))]

  if (length(shapefiles) == 0) {
    stop(paste0("Error: No range shapefiles (BUFF_*, MPC_*, MST_*) found in directory: ", shapefile_dir,
                "\nFound ", length(all_shapefiles), " total shapefiles, but none with BUFF_, MPC_, or MST_ prefix."))
  }

  message(paste0("Found ", length(shapefiles), " shapefiles"))

  # Load study area if provided
  if (!is.null(shape_file)) {
    if (inherits(shape_file, "sf")) {
      study_area <- shape_file
    } else if (inherits(shape_file, "SpatVector")) {
      study_area <- st_as_sf(shape_file)
    } else {
      study_area <- st_read(shape_file, quiet = TRUE)
    }
    bbox <- st_bbox(study_area)
  } else {
    # Use bounding box from all shapefiles
    all_bbox <- lapply(shapefiles, function(f) st_bbox(st_read(f, quiet = TRUE)))
    bbox <- c(
      xmin = min(sapply(all_bbox, function(b) b["xmin"])),
      ymin = min(sapply(all_bbox, function(b) b["ymin"])),
      xmax = max(sapply(all_bbox, function(b) b["xmax"])),
      ymax = max(sapply(all_bbox, function(b) b["ymax"]))
    )
  }

  # Create grid
  message(paste0("Creating grid with resolution: ", resolution[1], " x ", resolution[2], " degrees"))

  grid_polygon <- st_make_grid(
    st_as_sfc(bbox),
    cellsize = resolution,
    what = "polygons"
  )
  grid_sf <- st_sf(geometry = grid_polygon)
  grid_sf$grid_id <- 1:nrow(grid_sf)

  # Calculate centroids
  centroids <- st_centroid(grid_sf)
  centroid_coords <- st_coordinates(centroids)
  grid_sf$centroid_lon <- centroid_coords[, 1]
  grid_sf$centroid_lat <- centroid_coords[, 2]

  # Process each shapefile
  occurrences_list <- list()

  for (shp_file in shapefiles) {
    # Extract species name from filename (remove BUFF_, MPC_, or MST_ prefix)
    species_name <- tools::file_path_sans_ext(basename(shp_file))
    species_name <- gsub("^(BUFF_|MPC_|MST_)", "", species_name)  # Remove prefix
    species_name <- gsub("_", " ", species_name)  # Replace underscores with spaces

    message(paste0("Processing: ", species_name))

    # Read shapefile
    species_range <- st_read(shp_file, quiet = TRUE)

    # Find intersecting grid cells
    intersects <- st_intersects(grid_sf, species_range, sparse = FALSE)
    intersecting_cells <- grid_sf[apply(intersects, 1, any), ]

    if (nrow(intersecting_cells) > 0) {
      # Create occurrence records from centroids
      for (i in 1:nrow(intersecting_cells)) {
        occurrences_list[[length(occurrences_list) + 1]] <- data.frame(
          species = species_name,
          longitude = intersecting_cells$centroid_lon[i],
          latitude = intersecting_cells$centroid_lat[i],
          stringsAsFactors = FALSE
        )
      }
    } else {
      warning(paste0("No grid cells found for species: ", species_name))
    }
  }

  # Combine all occurrences
  if (length(occurrences_list) == 0) {
    stop("Error: No valid occurrences generated from shapefiles.")
  }

  dados <- do.call(rbind, occurrences_list)
  rownames(dados) <- NULL

  # Sort by species
  dados <- dados[order(dados$species), ]

  message(paste0("Generated ", nrow(dados), " grid cell centroids from ",
                 length(unique(dados$species)), " species"))

  return(dados)
}


#' Format data for NDM
#' @keywords internal
format_for_ndm <- function(dados) {

  output <- list()
  n <- 0
  i <- 1
  t <- 1
  name <- ""
  size <- nrow(dados)
  final_species <- length(unique(dados[[1]]))

  # Header
  output$c1[t] <- 'longlat'
  output$c2[t] <- ''
  t <- t + 1

  output$c1[t] <- paste('spp', final_species)
  output$c2[t] <- ''
  t <- t + 1

  output$c1[t] <- 'xydata'
  output$c2[t] <- ''
  t <- t + 1

  # Data
  while(i <= size) {
    if (as.character(dados[[1]][i]) != name) {
      output$c1[t] <- paste0("sp ", n)
      output$c2[t] <- paste0("[", as.character(dados[[1]][i]), "]")
      name <- as.character(dados[[1]][i])
      n <- n + 1
      t <- t + 1
    }

    output$c1[t] <- paste(as.character(dados[[2]][i]), ',', sep = '')
    output$c2[t] <- as.character(dados[[3]][i])
    t <- t + 1
    i <- i + 1
  }

  output <- na.omit(output)
  rownames(output) <- NULL

  return(output)
}


#' Write NDM file
#' @keywords internal
write_ndm_file <- function(output, output_file, output_dir) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste0("Created output directory: ", output_dir))
  }

  output_path <- file.path(output_dir, output_file)

  write.table(
    output,
    file = output_path,
    sep = "",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


#' Report conversion results
#' @keywords internal
report_conversion <- function(dados, output_file, output_dir, mode) {

  output_path <- file.path(output_dir, output_file)

  message("\n[OK] NDM/VNDM XYD file created successfully")
  message(paste0("  Output file: ", output_path))
  message(paste0("  Mode: ", mode))

  if (mode == "shapefile") {
    message(paste0("  Grid cell centroids: ", nrow(dados)))
  } else {
    message(paste0("  Occurrences: ", nrow(dados)))
  }

  message(paste0("  Species: ", length(unique(dados[[1]]))))
  message("\nNext steps:")
  message("  1. Open NDM/VNDM software")
  message(paste0("  2. Load file: ", output_path))
  message("  3. Set parameters and run analysis")
}
#' Tree Utility Functions for PAE-PCE Analysis
#'
#' Functions to prepare phylogenetic trees for PAE-PCE analysis
#'

#' Check if Tree has Branch Lengths
#'
#' Determine if a phylogenetic tree has branch length information
#'
#' @param tree phylo object
#'
#' @return Logical, TRUE if tree has branch lengths
#'
#' @export
has_branch_lengths <- function(tree) {
  !is.null(tree$edge.length) && length(tree$edge.length) > 0
}

#' Check if Tree is Ultrametric
#'
#' Determine if a tree is ultrametric (all tips equidistant from root)
#'
#' @param tree phylo object
#'
#' @return Logical, TRUE if tree is ultrametric
#'
#' @export
is_ultrametric <- function(tree) {
  
  if (!has_branch_lengths(tree)) {
    return(FALSE)
  }
  
  tryCatch({
    # Calculate distance from root to each tip
    distances <- ape::node.depth.edgelength(tree)
    tip_distances <- distances[1:length(tree$tip.label)]
    
    # Check if all tips have approximately the same distance
    max_dist <- max(tip_distances)
    min_dist <- min(tip_distances)
    
    # Allow small tolerance for floating point errors
    tolerance <- 1e-6
    (max_dist - min_dist) < tolerance
  }, error = function(e) {
    FALSE
  })
}

#' Convert Tree to Ultrametric
#'
#' Convert a non-ultrametric tree to ultrametric by scaling branch lengths
#' so all tips are equidistant from the root. This is optional and can be done
#' by user request before PAE-PCE analysis.
#'
#' @param tree phylo object with branch lengths
#' @param method Method for ultrametricization:
#'   - "equal": All branch lengths become 1 (equal length tree)
#'   - "scale": Scale existing branch lengths proportionally
#'   - "penalized": Use penalized least squares (requires phytools)
#'
#' @return Ultrametric phylo object
#'
#' @export
make_ultrametric <- function(tree, method = "equal") {
  
  if (!has_branch_lengths(tree)) {
    stop("Tree must have branch lengths to be made ultrametric")
  }
  
  if (is_ultrametric(tree)) {
    message("Tree is already ultrametric")
    return(tree)
  }
  
  tryCatch({
    if (method == "equal") {
      # Simple method: all branches get equal length
      tree$edge.length <- rep(1, length(tree$edge.length))
      return(tree)
    }
    
    else if (method == "scale") {
      # Scale method: proportionally adjust to make ultrametric
      # Calculate current distances
      distances <- ape::node.depth.edgelength(tree)
      tip_distances <- distances[1:length(tree$tip.label)]
      
      # Target distance (maximum current distance)
      target_dist <- max(tip_distances)
      
      # Scale factor for each edge
      tree$edge.length <- tree$edge.length * (target_dist / max(tip_distances))
      
      return(tree)
    }
    
    else if (method == "penalized") {
      # Penalized least squares method (requires phytools)
      if (!requireNamespace("phytools", quietly = TRUE)) {
        warning("phytools not available, using 'equal' method instead")
        tree$edge.length <- rep(1, length(tree$edge.length))
        return(tree)
      }
      
      tree_ultra <- phytools::force.ultrametric(tree, method = "extend")
      return(tree_ultra)
    }
    
    else {
      stop("Unknown method: ", method)
    }
  }, error = function(e) {
    stop("Error making tree ultrametric: ", e$message)
  })
}

#' Get Tree Statistics
#'
#' Calculate summary statistics for a phylogenetic tree
#'
#' @param tree phylo object
#'
#' @return List with tree statistics
#'
#' @export
get_tree_stats <- function(tree) {
  
  stats <- list(
    n_taxa = length(tree$tip.label),
    n_nodes = tree$Nnode,
    n_edges = nrow(tree$edge),
    has_branch_lengths = has_branch_lengths(tree),
    is_ultrametric = is_ultrametric(tree)
  )
  
  if (has_branch_lengths(tree)) {
    stats$total_length <- sum(tree$edge.length)
    stats$mean_branch_length <- mean(tree$edge.length)
    stats$min_branch_length <- min(tree$edge.length)
    stats$max_branch_length <- max(tree$edge.length)
    
    # Distance from root to tips
    distances <- ape::node.depth.edgelength(tree)
    tip_distances <- distances[1:length(tree$tip.label)]
    stats$max_tip_distance <- max(tip_distances)
    stats$min_tip_distance <- min(tip_distances)
  }
  
  stats
}

#' Validate Tree for Analysis
#'
#' Check if tree is suitable for PAE-PCE and BioGeoBEARS analyses
#'
#' @param tree phylo object
#' @param for_analysis Type of analysis: "pae_pce" or "biogeobears"
#'
#' @return List with validation results
#'
#' @export
validate_tree_for_analysis <- function(tree, for_analysis = "pae_pce") {
  
  validation <- list(
    is_valid = TRUE,
    errors = c(),
    warnings = c(),
    suggestions = c()
  )
  
  # Basic checks
  if (!inherits(tree, "phylo")) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Object is not a phylo object")
    return(validation)
  }
  
  if (length(tree$tip.label) < 2) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Tree must have at least 2 taxa")
    return(validation)
  }
  
  # Check for branch lengths
  if (!has_branch_lengths(tree)) {
    validation$warnings <- c(validation$warnings, "Tree does not have branch lengths")
  }
  
  # Check for duplicate taxa names
  if (length(unique(tree$tip.label)) < length(tree$tip.label)) {
    validation$is_valid <- FALSE
    validation$errors <- c(validation$errors, "Tree has duplicate taxon names")
  }
  
  validation
}
#' Utility Functions for BioGeoBEARS Shiny App
#'
#' Collection of helper functions for data processing and validation

#' Calculate Likelihood Ratio Test
#'
#' Performs likelihood ratio test between two nested models
#'
#' @param res_simple Results from simpler (nested) model
#' @param res_complex Results from more complex model
#'
#' @return List with LRT statistics
#'
#' @keywords internal
calc_lrt <- function(res_simple, res_complex) {
  lnl_simple <- res_simple$total_loglikelihood
  lnl_complex <- res_complex$total_loglikelihood

  # Calculate degrees of freedom
  df_simple <- length(res_simple$inputs$BioGeoBEARS_model_object@params_table[
    res_simple$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])
  df_complex <- length(res_complex$inputs$BioGeoBEARS_model_object@params_table[
    res_complex$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])

  df <- df_complex - df_simple

  # Test statistic
  chisq <- 2 * (lnl_complex - lnl_simple)

  # P-value
  pval <- stats::pchisq(chisq, df = df, lower.tail = FALSE)

  list(
    model_simple = res_simple$inputs$description,
    model_complex = res_complex$inputs$description,
    lnl_simple = lnl_simple,
    lnl_complex = lnl_complex,
    df = df,
    chisq = chisq,
    pval = pval,
    significant = pval < 0.05
  )
}

#' Calculate AICc
#'
#' Calculate Akaike Information Criterion corrected for small samples
#'
#' @param lnl Log-likelihood value
#' @param k Number of parameters
#' @param n Sample size
#'
#' @return AICc value
#'
#' @keywords internal
calc_aicc <- function(lnl, k, n) {
  aic <- 2 * k - 2 * lnl
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  aicc
}

#' Validate Occurrence Data
#'
#' Check if occurrence data has required columns and valid format
#'
#' @param data Data frame with occurrence data
#'
#' @return Logical, TRUE if valid
#'
#' @keywords internal
validate_occurrence_data <- function(data) {
  required_cols <- c("species", "longitude", "latitude")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: species, longitude, latitude")
  }

  if (nrow(data) == 0) {
    stop("Data is empty")
  }

  if (any(is.na(data$longitude)) || any(is.na(data$latitude))) {
    stop("Coordinates contain NA values")
  }

  TRUE
}

#' Validate Tree File
#'
#' Check if tree file is valid Newick format
#'
#' @param tree_file Path to tree file
#'
#' @return Logical, TRUE if valid
#'
#' @keywords internal
validate_tree_file <- function(tree_file) {
  tryCatch({
    tree <- ape::read.tree(tree_file)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Format Model Results
#'
#' Format BioGeoBEARS results for display
#'
#' @param results List of BioGeoBEARS results
#'
#' @return Data frame with formatted results
#'
#' @keywords internal
format_results <- function(results) {
  if (is.null(results) || length(results) == 0) {
    return(data.frame())
  }

  # Extract key statistics for each model
  result_list <- lapply(results, function(res) {
    data.frame(
      Model = res$inputs$description,
      LnL = round(res$total_loglikelihood, 3),
      d = round(res$inputs$BioGeoBEARS_model_object@params_table["d", "est"], 4),
      e = round(res$inputs$BioGeoBEARS_model_object@params_table["e", "est"], 4),
      j = if ("j" %in% rownames(res$inputs$BioGeoBEARS_model_object@params_table)) {
        round(res$inputs$BioGeoBEARS_model_object@params_table["j", "est"], 4)
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind, result_list)
  rownames(results_df) <- NULL
  results_df
}
#' Visualization Functions for BioGeoBEARS Shiny
#'
#' Functions for creating interactive maps, phylogenetic trees, and model comparison plots
#'

#' Plot Phylogenetic Tree
#'
#' Create a phylogenetic tree visualization with branch lengths
#'
#' @param tree phylo object
#' @param title Title for the plot
#'
#' @return NULL (creates plot as side effect)
#'
#' @export
plot_phylogenetic_tree <- function(tree, title = "Phylogenetic Tree") {
  
  plot(tree, main = title, cex = 0.8)
  
  # Add axis information
  axisPhylo()
  
  invisible(NULL)
}

#' Plot Species Distribution Map
#'
#' Create an interactive map showing species occurrences with optional extrapolation polygons
#'
#' @param occurrence_data Data frame with species, longitude, latitude columns
#' @param geometry Optional spatial geometry object (polygons from extrapolation)
#' @param method Optional method name ("BUFF", "MPC", "MST") for labeling
#'
#' @return leaflet map object
#'
#' @export
plot_distribution_map <- function(occurrence_data, geometry = NULL, method = NULL) {
  
  # Create base map
  m <- leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::setView(lng = mean(occurrence_data$long, na.rm = TRUE),
                     lat = mean(occurrence_data$lat, na.rm = TRUE),
                     zoom = 4)
  
  # Add occurrence points
  species_unique <- unique(occurrence_data$spp)
  n_species <- length(species_unique)
  
  # Choose appropriate color palette based on number of species
  if (n_species <= 3) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set1")
  } else if (n_species <= 9) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set1")
  } else if (n_species <= 12) {
    colors <- RColorBrewer::brewer.pal(n_species, "Set3")
  } else {
    # For more than 12 species, use a continuous palette
    colors <- grDevices::hcl.colors(n_species, "Spectral")
  }
  
  for (i in seq_along(species_unique)) {
    sp <- species_unique[i]
    sp_data <- occurrence_data[occurrence_data$spp == sp, ]
    
    m <- m %>%
      leaflet::addCircleMarkers(
        data = sp_data,
        lng = ~long,
        lat = ~lat,
        radius = 5,
        color = colors[i],
        fillOpacity = 0.7,
        popup = ~paste(spp, "<br>Lon:", round(long, 2), "<br>Lat:", round(lat, 2)),
        group = sp
      )
  }
  
  # Add extrapolation polygons if provided
  if (!is.null(geometry)) {
    m <- biogeoshiny::add_extrapolation_polygons(
      map = m,
      geometry = geometry,
      method = ifelse(is.null(method), "Extrapolation", method),
      color = "blue",
      opacity = 0.5
    )
  }
  
  # Prepare layer groups for control
  layer_groups <- c(species_unique)
  if (!is.null(geometry)) {
    layer_groups <- c(layer_groups, ifelse(is.null(method), "Extrapolation", method))
  }
  
  # Add layer control
  m <- m %>%
    leaflet::addLayersControl(
      overlayGroups = layer_groups,
      options = leaflet::layersControlOptions(collapsed = FALSE)
    )
  
  m
}

#' Plot Extrapolation Polygons on Map
#'
#' Add extrapolation polygons (buffers, convex hulls, MST) to a leaflet map
#'
#' @param map leaflet map object
#' @param geometry Spatial geometry object (SpatVector or sf)
#' @param species_names Character vector of species names
#'
#' @return Updated leaflet map object
#'
#' @export
plot_extrapolation_polygons <- function(map, geometry, species_names) {
  
  if (is.null(geometry)) {
    return(map)
  }
  
  tryCatch({
    # Handle different geometry types
    if (inherits(geometry, "SpatVector")) {
      # Convert terra SpatVector to sf
      geometry <- sf::st_as_sf(geometry)
    }
    
    if (inherits(geometry, "sf") || inherits(geometry, "data.frame")) {
      # Add polygons to map
      colors <- RColorBrewer::brewer.pal(min(length(species_names), 12), "Set1")
      
      for (i in seq_along(species_names)) {
        if (i <= nrow(geometry)) {
          map <- map %>%
            leaflet::addPolygons(
              data = geometry[i, ],
              color = colors[i],
              fillColor = colors[i],
              fillOpacity = 0.3,
              weight = 2,
              popup = paste("Species:", species_names[i]),
              group = paste("Polygon:", species_names[i])
            )
        }
      }
    }
    
    map
  }, error = function(e) {
    warning("Could not add polygons to map: ", e$message)
    map
  })
}

#' Plot Model Comparison by AICc
#'
#' Create an interactive bar plot comparing models by AICc
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_model_comparison <- function(comparison_table) {
  
  # Sort by AICc
  comparison_table <- comparison_table[order(comparison_table$AICc), ]
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~reorder(Model, AICc),
    y = ~AICc,
    type = "bar",
    marker = list(color = ~DeltaAICc),
    text = ~paste("AICc:", round(AICc, 2), "<br>Î”AICc:", round(DeltaAICc, 2)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = "Model Comparison by AICc",
      xaxis = list(title = "Model"),
      yaxis = list(title = "AICc"),
      hovermode = "closest"
    )
  
  p
}

#' Plot Model Likelihood
#'
#' Create a plot showing log-likelihood values for different models
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_model_likelihood <- function(comparison_table) {
  
  # Sort by LnL
  comparison_table <- comparison_table[order(comparison_table$LnL, decreasing = TRUE), ]
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~reorder(Model, LnL),
    y = ~LnL,
    type = "bar",
    marker = list(color = "steelblue"),
    text = ~paste("LnL:", round(LnL, 2)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = "Model Log-Likelihood Comparison",
      xaxis = list(title = "Model"),
      yaxis = list(title = "Log-Likelihood"),
      hovermode = "closest"
    )
  
  p
}

#' Plot Parameter Comparison
#'
#' Create a plot comparing parameter values across models
#'
#' @param comparison_table Data frame with model comparison results
#' @param parameter Parameter to plot ("d", "e", or "j")
#'
#' @return plotly object
#'
#' @export
plot_parameter_comparison <- function(comparison_table, parameter = "d") {
  
  if (!(parameter %in% c("d", "e", "j"))) {
    stop("Parameter must be 'd', 'e', or 'j'")
  }
  
  # Create bar plot
  p <- plotly::plot_ly(
    data = comparison_table,
    x = ~Model,
    y = as.formula(paste("~", parameter)),
    type = "bar",
    marker = list(color = "coral"),
    text = ~paste(parameter, "=", get(parameter)),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      title = paste("Parameter", parameter, "Comparison"),
      xaxis = list(title = "Model"),
      yaxis = list(title = paste("Parameter", parameter)),
      hovermode = "closest"
    )
  
  p
}

#' Plot Parameter Heatmap
#'
#' Create a heatmap showing parameter values across models
#'
#' @param comparison_table Data frame with model comparison results
#'
#' @return plotly object
#'
#' @export
plot_parameter_heatmap <- function(comparison_table) {
  
  # Extract parameter values
  params_matrix <- as.matrix(comparison_table[, c("d", "e")])
  rownames(params_matrix) <- comparison_table$Model
  
  # Create heatmap
  p <- plotly::plot_ly(
    z = params_matrix,
    x = colnames(params_matrix),
    y = rownames(params_matrix),
    type = "heatmap",
    colorscale = "Viridis"
  ) %>%
    plotly::layout(
      title = "Parameter Heatmap",
      xaxis = list(title = "Parameter"),
      yaxis = list(title = "Model")
    )
  
  p
}

.resolve_pae_pce_impl <- local({
  cached_fun <- NULL

  function(force_reload = FALSE) {
    if (!force_reload && is.function(cached_fun)) {
      return(cached_fun)
    }

    candidate_paths <- unique(c(
      file.path(getwd(), "R", "pae_pce.R"),
      file.path(getwd(), "biogeoshiny", "R", "pae_pce.R"),
      file.path(getwd(), "biogeoshiny_improved_0.1.0", "biogeoshiny", "R", "pae_pce.R"),
      file.path(dirname(getwd()), "R", "pae_pce.R"),
      file.path(dirname(getwd()), "biogeoshiny", "R", "pae_pce.R"),
      file.path(dirname(getwd()), "biogeoshiny_improved_0.1.0", "biogeoshiny", "R", "pae_pce.R"),
      file.path("..", "R", "pae_pce.R")
    ))

    for (candidate in candidate_paths) {
      if (!file.exists(candidate)) {
        next
      }

      temp_env <- new.env(parent = globalenv())
      source(candidate, local = temp_env)
      candidate_fun <- get0("pae_pce", envir = temp_env, mode = "function", inherits = FALSE)

      if (is.function(candidate_fun)) {
        cached_fun <<- candidate_fun
        return(cached_fun)
      }
    }

    if (requireNamespace("biogeoshiny", quietly = TRUE)) {
      candidate_fun <- getExportedValue("biogeoshiny", "pae_pce")
      if (is.function(candidate_fun)) {
        cached_fun <<- candidate_fun
        return(cached_fun)
      }
    }

    stop("Could not resolve a pae_pce() implementation.")
  }
})

pae_pce <- function(...) {
  .resolve_pae_pce_impl()(...)
}

