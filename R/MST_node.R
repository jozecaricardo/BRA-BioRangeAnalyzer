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

        # back-transforming lines in points:
        suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

        # presence-absence matrix:
        ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(taxon), dimnames = list(seq(1:ncellras),
                                                                                unique(rownames(tempo))[c(1:length(taxon))]))  # tabela com todas as celulas
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
        proj4string(linhasRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        linhasRaster <- mask(linhasRaster, as(shape_file, 'Spatial'))

        writeRaster(linhasRaster, "out_MST/presence_mintreeall.tif",
                    overwrite = TRUE)

        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
             alpha = transp)

        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
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


          # back-transforming lines in points:
          suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))


          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
          proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          lista_r[[conta]] <- mask(lista_r[[conta]], as(shape_file, 'Spatial'))
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE,
               col = cols1[j], alpha = transp)
          # colocando labels:
          teste <- tabelao[j,j] # posicoes dos taxons do no
          text(lista_r[[conta]], labels = rep(teste, dim(tempo)[1]), cex = 0.8, pos = 2, col = cols1[j])

          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
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

        # back-transforming lines in points:
        suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

        # presence-absence matrix:
        # ncellras <- ncell(r)
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                                                                                                  unique(rownames(tempo))))
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
        proj4string(linhasRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        linhasRaster <- mask(linhasRaster, as(shape_file, 'Spatial'))
        writeRaster(linhasRaster, "out_MST/presence_mintreeall.tif",
                    overwrite = TRUE)
        plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1], alpha = transp)

        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
        }

        coor.l <- na.exclude(coor.l)
        coor.ll <- rbind(coor.l, rep(0, ncol(coor.l)))
        rownames(coor.ll) <- c(rownames(coor.l), 'ROOT')
        write.table(x = coor.ll, file = 'out_MST/pres_abs.txt', sep = '\t')

        return(coor.ll)


      } else if (mintreeall == FALSE){

        conta <- 0
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))),
                         dimnames = list(seq(1:ncellras), unique(rownames(tempo))))
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


          # back-transforming lines in points:
          suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))


          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
          proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          lista_r[[conta]] <- mask(lista_r[[conta]], as(shape_file, 'Spatial'))
          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)

          # teste <- tabelao[j,j] # posicoes dos taxons do no
          teste <- conta
          track_numbers[j] <- teste  # where 'teste' is the calculated track number
          text(tempoo, labels = rep(teste, dim(tempo)[1]), cex = 0.5, pos = 2, col = cols1[j])

          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
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

        # back-transforming lines in points:
        suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo))), dimnames = list(seq(1:ncellras),
                                                                                                  unique(rownames(tempo))))  # tabela com todas as celulas
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
        proj4string(linhasRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        linhasRaster <- mask(linhasRaster, as(shape_file, 'Spatial'))
        writeRaster(linhasRaster, "out_MST/presence_mst_mintreeall_ancterminal.tif",
                    overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)

        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
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

        # back-transforming lines in points:
        suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

        # presence-absence matrix:
        coor.l <- matrix(NA, nr = ncellras, nc = length(unique(rownames(tempo.temp2))), dimnames = list(seq(1:ncellras),
                                                                                                        unique(rownames(tempo.temp2)))[c(1:2)])  # tabela com todas as celulas
        linhasRaster <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
        proj4string(linhasRaster) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        linhasRaster <- mask(linhasRaster, as(shape_file, 'Spatial'))
        writeRaster(linhasRaster, "out_MST/presence_mst_mintreeall_ancterminal.tif",
                    overwrite = TRUE)
        # plot(linhasRaster, axes = FALSE, legend = FALSE, add = TRUE, col = cols1[1],
        #     alpha = transp)

        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == FALSE){
            coor.l[i, ] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(linhasRaster[i]) == TRUE){
            coor.l[i, ] <- 0
          }
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

        # back-transforming lines in points:
        suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

        lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
        proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
        lista_r[[conta]] <- mask(lista_r[[conta]], as(shape_file, 'Spatial'))
        writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                             collapse = ''), overwrite = TRUE)

        plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE,
             col = cols1[unique(tempo.n[which(rownames(tempo.n) == j), 3])], alpha = transp)

        # presence-absence matrix:
        # preparing the matrix
        for(i in 1:ncellras){
          if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
            coor.l[i, conta] <- 1
          } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
            coor.l[i, conta] <- 0
          }
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

          # back-transforming lines in points:
          suppressWarnings(pontos_linha2 <- spsample(pontos_linha, n = 100, type = 'regular'))

          lista_r[[conta]] <- rasterize(pontos_linha2, r, field = 1) # raster com as presencas
          proj4string(lista_r[[conta]]) <- CRS("+proj=longlat +datum=WGS84") # datum WGS84
          lista_r[[conta]] <- mask(lista_r[[conta]], as(shape_file, 'Spatial'))
          writeRaster(lista_r[[conta]], paste0(c("out_MST/presence_mst_", j, ".tif"),
                                               collapse = ''), overwrite = TRUE)

          plot(lista_r[[conta]], axes = FALSE, legend = FALSE, add = TRUE, col = cols1[j],
               alpha = transp)

          # presence-absence matrix:
          # preparing the matrix
          for(i in 1:ncellras){
            if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == FALSE){
              coor.l[i, conta] <- 1
            } else if(is.na(r[i]) == FALSE && is.na(lista_r[[conta]][i]) == TRUE){
              coor.l[i, conta] <- 0
            }
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

