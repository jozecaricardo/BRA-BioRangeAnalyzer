function(input, output, session) {
  
  # Reactive storage for data
  data_store <- reactiveValues(
    occurrence = NULL,
    tree = NULL,
    tree_ultrametric = NULL,
    pres_abs = NULL,
    pres_abs_regular = NULL,
    pres_abs_irregular = NULL,
    geometry = NULL,
    extrap_method = NULL,
    loaded_shapefiles = list(),
    loaded_rasters = list(),
    shapefile_info = list(),
    species_colors = list(),
    distance_matrix = NULL,
    dispersal_matrix = NULL,
    extinction_matrix = NULL,
    dispersal_multipliers = NULL,
    env_distance_matrix = NULL,
    areas_allowed_matrix = NULL,
    areas_adjacency_matrix = NULL,
    time_periods = NULL,
    mst_pres_abs = NULL,
    mst_context = NULL,
    irregular_bins_richness = NULL,
    irregular_bins_species_table = NULL,
    irregular_bins_id_column = NULL,
    irregular_bins_bin_id_mapping = NULL,
    irregular_bins_method = NULL,
    matrix_is_irregular_aggregated = FALSE,
    visual_layers_method = NULL,
    problematic_taxa = NULL,
    singleton_taxa = NULL,
    doubleton_taxa = NULL,
    duplicate_rows = NULL,
    duplicate_taxa = NULL,
    secondary_study_areas = list(),
    regular_grid_all_sf = NULL,
    regular_grid_presence_sf = NULL,
    regular_grid_presence_points_sf = NULL,
    regular_grid_presence_extrap_sf = NULL,
    regular_grid_presence_points_by_taxon = list(),
    regular_grid_presence_extrap_by_taxon = list(),
    regular_grid_presence_extrap_label = NULL,
    regular_grid_presence_by_taxon = list(),
    regular_grid_res = NULL,
    extrap_log = character(0),
    last_extrap_started_at = NULL,
    last_extrap_taxa = character(0),
    analysis_occurrence = NULL,
    analysis_tree = NULL
  )

  `%||%` <- function(x, y) if (is.null(x)) y else x

  reset_extrap_log <- function(header = NULL) {
    if (is.null(header)) {
      data_store$extrap_log <- character(0)
    } else {
      data_store$extrap_log <- as.character(header)
    }
  }

  append_extrap_log <- function(...) {
    parts <- unlist(list(...), use.names = FALSE)
    parts <- as.character(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) {
      return(invisible(NULL))
    }
    data_store$extrap_log <- c(data_store$extrap_log, parts)
    invisible(NULL)
  }

  capture_analysis_output <- function(expr) {
    collected <- character(0)
    result <- withCallingHandlers(
      {
        tmp <- NULL
        text_out <- capture.output(tmp <- force(expr), type = "output")
        if (length(text_out) > 0) {
          collected <- c(collected, text_out)
        }
        tmp
      },
      message = function(m) {
        collected <<- c(collected, paste0("[message] ", conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        collected <<- c(collected, paste0("[warning] ", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )

    list(result = result, log = collected)
  }

  extract_clade_subset <- function(occ_data, tree, node_id = NULL) {
    if (is.null(tree)) {
      stop("Clade filter requires a loaded phylogenetic tree.")
    }
    if (is.null(node_id) || !nzchar(as.character(node_id))) {
      stop("Please select an internal node for clade filtering.")
    }

    node_num <- suppressWarnings(as.integer(node_id))
    if (is.na(node_num)) {
      stop("Invalid clade node ID.")
    }

    n_tips <- ape::Ntip(tree)
    n_nodes <- tree$Nnode
    valid_internal <- (n_tips + 1):(n_tips + n_nodes)
    if (!(node_num %in% valid_internal)) {
      stop("Selected node is not a valid internal node in the loaded tree.")
    }

    clade_tree <- ape::extract.clade(tree, node = node_num)
    clade_tips <- clade_tree$tip.label
    occ_subset <- occ_data[occ_data$spp %in% clade_tips, , drop = FALSE]
    if (nrow(occ_subset) == 0) {
      stop("No occurrence records match the selected clade tips.")
    }

    list(
      occurrence = occ_subset,
      tree = clade_tree,
      tips = clade_tips,
      node = node_num
    )
  }

  show_tree_node_labels <- reactiveVal(FALSE)
  
  # Sync checkbox with reactive value
  observeEvent(input$show_node_labels, {
    show_tree_node_labels(isTRUE(input$show_node_labels))
  })

  safe_valid_sf <- function(sf_obj, context_label = "geometry") {
    n_before <- nrow(sf_obj)
    sf_obj <- tryCatch(sf::st_make_valid(sf_obj), error = function(e) {
      if (requireNamespace("lwgeom", quietly = TRUE)) {
        tryCatch(lwgeom::st_make_valid(sf_obj), error = function(e2) sf_obj)
      } else {
        sf_obj
      }
    })
    valid_idx <- tryCatch(sf::st_is_valid(sf_obj, NA_on_exception = TRUE), error = function(e) rep(FALSE, nrow(sf_obj)))
    invalid_idx <- is.na(valid_idx) | !valid_idx

    if (any(invalid_idx)) {
      geom_types <- as.character(sf::st_geometry_type(sf_obj, by_geometry = TRUE))
      bufferable <- geom_types %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION", "GEOMETRY")
      repair_idx <- which(invalid_idx & bufferable)
      if (length(repair_idx) > 0) {
        repaired <- tryCatch(sf::st_buffer(sf_obj[repair_idx, , drop = FALSE], dist = 0), error = function(e) NULL)
        if (!is.null(repaired)) {
          sf_obj[repair_idx, ] <- repaired
        }
      }
      valid_idx <- tryCatch(sf::st_is_valid(sf_obj, NA_on_exception = TRUE), error = function(e) rep(FALSE, nrow(sf_obj)))
      invalid_idx <- is.na(valid_idx) | !valid_idx
    }

    if (length(valid_idx) == nrow(sf_obj) && any(valid_idx %in% FALSE | is.na(valid_idx))) {
      sf_obj <- sf_obj[valid_idx %in% TRUE, , drop = FALSE]
    }

    if (nrow(sf_obj) > 0) {
      empty_idx <- tryCatch(sf::st_is_empty(sf_obj), error = function(e) rep(FALSE, nrow(sf_obj)))
      if (length(empty_idx) == nrow(sf_obj) && any(empty_idx)) {
        sf_obj <- sf_obj[!empty_idx, , drop = FALSE]
      }
    }
    n_after <- nrow(sf_obj)
    if (n_after < n_before) {
      message(sprintf("%s cleanup removed %d invalid/empty feature(s).", context_label, n_before - n_after))
    }
    sf_obj
  }

  normalize_to_wgs84_sf <- function(obj) {
    sf_obj <- sf::st_as_sf(obj)
    sf_obj <- safe_valid_sf(sf_obj, context_label = "normalize_to_wgs84_sf")
    if (nrow(sf_obj) == 0) {
      return(sf_obj)
    }
    crs_obj <- sf::st_crs(sf_obj)
    if (is.na(crs_obj)) {
      sf::st_crs(sf_obj) <- 4326
      return(sf_obj)
    }
    if (!isTRUE(sf::st_is_longlat(sf_obj)) || !isTRUE(crs_obj$epsg == 4326)) {
      sf_obj <- sf::st_transform(sf_obj, 4326)
    }
    sf_obj
  }

  build_bbox_grid_shape <- function(shape_file, pad_degrees = 0) {
    shp_sf <- normalize_to_wgs84_sf(shape_file)
    bb <- sf::st_bbox(shp_sf)
    if (pad_degrees > 0) {
      bb["xmin"] <- bb["xmin"] - pad_degrees
      bb["xmax"] <- bb["xmax"] + pad_degrees
      bb["ymin"] <- bb["ymin"] - pad_degrees
      bb["ymax"] <- bb["ymax"] + pad_degrees
    }
    bb_poly <- sf::st_as_sfc(bb)
    bb_sf <- sf::st_sf(id = "bbox_grid", geometry = bb_poly, crs = 4326)
    terra::vect(bb_sf)
  }

  safe_grid_resolution_tag <- function(grid_res) {
    tag <- gsub("[^0-9A-Za-z]+", "_", format(grid_res, scientific = FALSE, trim = TRUE))
    if (!nzchar(tag)) tag <- "grid"
    tag
  }

  build_regular_grid_sf <- function(shape_file, grid_res) {
    if (is.null(shape_file) || !is.finite(grid_res) || grid_res <= 0) {
      stop("Grid resolution must be > 0 to build regular grid layers.")
    }

    study_sf <- safe_valid_sf(normalize_to_wgs84_sf(shape_file))
    if (nrow(study_sf) == 0) {
      stop("Study area has no valid geometry for regular grid generation.")
    }

    bb <- sf::st_bbox(study_sf)
    grid_r <- raster::raster(
      xmn = bb["xmin"],
      xmx = bb["xmax"],
      ymn = bb["ymin"],
      ymx = bb["ymax"],
      resolution = c(grid_res, grid_res),
      crs = sp::CRS("+proj=longlat +datum=WGS84")
    )

    grid_sf <- sf::st_as_sf(raster::rasterToPolygons(grid_r))
    grid_sf$grid_id <- seq_len(nrow(grid_sf))
    grid_sf <- grid_sf[, c("grid_id", "geometry")]
    normalize_to_wgs84_sf(grid_sf)
  }

  align_regular_matrix_to_grid <- function(pres_abs, grid_sf) {
    if (is.null(pres_abs) || !is.matrix(pres_abs)) {
      stop("Regular-grid matrix is missing or invalid.")
    }
    if (is.null(grid_sf) || nrow(grid_sf) == 0) {
      stop("Regular grid polygons are missing.")
    }

    mat <- pres_abs
    if ("ROOT" %in% rownames(mat)) {
      mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
    }

    full_ids <- as.integer(grid_sf$grid_id)
    if (nrow(mat) == 0) {
      full_mat <- matrix(0, nrow = length(full_ids), ncol = ncol(pres_abs), dimnames = list(as.character(full_ids), colnames(pres_abs)))
      return(rbind(full_mat, ROOT = rep(0, ncol(full_mat))))
    }

    row_ids <- suppressWarnings(as.integer(rownames(mat)))
    if (all(is.na(row_ids)) && nrow(mat) == length(full_ids)) {
      row_ids <- full_ids
    }

    full_mat <- matrix(0, nrow = length(full_ids), ncol = ncol(mat), dimnames = list(as.character(full_ids), colnames(mat)))
    valid <- !is.na(row_ids) & as.character(row_ids) %in% rownames(full_mat)
    if (any(valid)) {
      full_mat[as.character(row_ids[valid]), ] <- mat[valid, , drop = FALSE]
    }

    full_mat[is.na(full_mat)] <- 0
    full_mat <- ifelse(full_mat > 0, 1, 0)
    full_mat <- as.matrix(full_mat)
    colnames(full_mat) <- colnames(mat)
    rownames(full_mat) <- as.character(full_ids)

    rbind(full_mat, ROOT = rep(0, ncol(full_mat)))
  }

  write_regular_grid_layers <- function(method, grid_sf, regular_matrix, grid_res) {
    if (is.null(grid_sf) || nrow(grid_sf) == 0 || is.null(regular_matrix) || !is.matrix(regular_matrix)) {
      return(invisible(NULL))
    }

    method_dir <- c(
      occurrence_only = "out_grid",
      buffer = "out_grid",
      convex_hull = "out_grid",
      mst = "out_grid"
    )
    out_subdir <- unname(method_dir[method])
    if (is.na(out_subdir) || !nzchar(out_subdir)) {
      return(invisible(NULL))
    }

    out_dir <- file.path(getwd(), out_subdir)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    delete_shapefile_set <- function(shp_path) {
      base <- tools::file_path_sans_ext(shp_path)
      exts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix")
      files <- paste0(base, exts)
      existing <- files[file.exists(files)]
      if (length(existing) > 0) {
        invisible(file.remove(existing))
      }
    }

    safe_write_shp <- function(obj, shp_path, label) {
      delete_shapefile_set(shp_path)
      tryCatch({
        sf::st_write(obj, shp_path, driver = "ESRI Shapefile", quiet = TRUE)
        append_extrap_log(paste0("Wrote regular grid layer (", label, "): ", basename(shp_path)))
        TRUE
      }, error = function(e) {
        append_extrap_log(paste0("[error] Could not write regular grid layer (", label, "): ", e$message))
        FALSE
      })
    }

    res_tag <- safe_grid_resolution_tag(grid_res)
    all_grid_path <- file.path(out_dir, paste0("GRIDS_allcells_", method, "_q", res_tag, ".shp"))
    safe_write_shp(grid_sf, all_grid_path, "all cells")

    mat <- regular_matrix
    if ("ROOT" %in% rownames(mat)) {
      mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
    }
    if (nrow(mat) == 0) {
      return(invisible(NULL))
    }

    occupied_ids <- suppressWarnings(as.integer(rownames(mat)))
    occupied_ids <- occupied_ids[rowSums(mat, na.rm = TRUE) > 0]
    occupied_ids <- unique(occupied_ids[!is.na(occupied_ids)])
    if (length(occupied_ids) == 0) {
      return(invisible(NULL))
    }

    occupied_grid <- grid_sf[grid_sf$grid_id %in% occupied_ids, c("grid_id", "geometry")]
    if (nrow(occupied_grid) == 0) {
      return(invisible(NULL))
    }
    occupied_path <- file.path(out_dir, paste0("GRIDS_presence_", method, "_q", res_tag, ".shp"))
    safe_write_shp(occupied_grid, occupied_path, "presence")

    invisible(NULL)
  }

  layer_species_name <- function(layer_info, layer_key) {
    if (!is.null(layer_info) && !is.null(layer_info$species_id) && !is.na(layer_info$species_id)) {
      return(as.character(layer_info$species_id))
    }
    basename(layer_key)
  }

  is_helper_points_layer <- function(base_name) {
    base_name <- as.character(base_name)
    grepl("^(MST_)?pointshape_", base_name, ignore.case = TRUE) ||
      grepl("^(MST_)?pointsshape_", base_name, ignore.case = TRUE) ||
      grepl("^(MST_)?points_.*onlyinternal$", base_name, ignore.case = TRUE)
  }

  is_mst_internal_helper_layer <- function(base_name, species_name = "") {
    base_name <- as.character(base_name)
    species_name <- as.character(species_name)
    grepl("mintreeall_onlyinternal", base_name, ignore.case = TRUE) ||
      grepl("mintreeall_onlyinternal", species_name, ignore.case = TRUE)
  }

  is_mst_onlyterminal_points_layer <- function(base_name, species_name = "") {
    base_name <- as.character(base_name)
    species_name <- as.character(species_name)
    grepl("points_mintreeall_onlyterminal", base_name, ignore.case = TRUE) ||
      grepl("points_mintreeall_onlyterminal", species_name, ignore.case = TRUE)
  }

  mst_grid_from_presence_raster <- function(loaded_rasters) {
    if (is.null(loaded_rasters) || length(loaded_rasters) == 0) {
      return(NULL)
    }

    rk <- names(loaded_rasters)
    target_idx <- grep("presence_mst_mintreeall_ancterminal", rk, ignore.case = TRUE)
    if (length(target_idx) == 0) {
      return(NULL)
    }

    raster_file <- loaded_rasters[[target_idx[1]]]
    if (is.null(raster_file) || !file.exists(raster_file)) {
      return(NULL)
    }

    rr <- terra::rast(raster_file)
    grid_vect <- terra::as.polygons(rr, na.rm = TRUE, dissolve = FALSE)
    if (is.null(grid_vect) || nrow(grid_vect) == 0) {
      return(NULL)
    }

    grid_sf <- sf::st_as_sf(grid_vect)
    normalize_to_wgs84_sf(grid_sf)
  }

  ensure_mst_node_intersection_outputs <- function(shape_file, grid_res) {
    if (is.null(shape_file) || !is.finite(grid_res) || grid_res <= 0) {
      return(FALSE)
    }

    out_dir <- file.path(getwd(), "out_MST")
    if (!dir.exists(out_dir)) {
      return(FALSE)
    }

    preferred_mst <- file.path(out_dir, "mst_ancterminal_mintreeall.shp")
    shp_files <- if (file.exists(preferred_mst)) {
      preferred_mst
    } else {
      list.files(out_dir, pattern = "^mst_.*\\.shp$", full.names = TRUE)
    }
    if (length(shp_files) == 0) {
      return(FALSE)
    }

    base_names <- tools::file_path_sans_ext(basename(shp_files))
    keep <- !grepl("mintreeall_onlyinternal", base_names, ignore.case = TRUE)
    shp_files <- shp_files[keep]
    base_names <- base_names[keep]
    if (length(shp_files) == 0) {
      return(FALSE)
    }

    ranked <- order(!grepl("ancterminal", base_names, ignore.case = TRUE), -as.numeric(file.info(shp_files)$mtime))
    mst_file <- shp_files[ranked[1]]

    mst_sf <- suppressWarnings(sf::st_read(mst_file, quiet = TRUE))
    mst_sf <- normalize_to_wgs84_sf(mst_sf)
    if (nrow(mst_sf) == 0) {
      return(FALSE)
    }

    gtypes <- as.character(sf::st_geometry_type(mst_sf, by_geometry = TRUE))
    line_idx <- !is.na(gtypes) & gtypes %in% c("LINESTRING", "MULTILINESTRING")
    coll_idx <- !is.na(gtypes) & gtypes %in% c("GEOMETRYCOLLECTION", "GEOMETRY")

    line_parts <- list()
    if (any(line_idx)) {
      line_parts[[length(line_parts) + 1]] <- mst_sf[line_idx, , drop = FALSE]
    }
    if (any(coll_idx)) {
      coll_lines <- suppressWarnings(sf::st_collection_extract(mst_sf[coll_idx, , drop = FALSE], "LINESTRING"))
      if (!is.null(coll_lines) && nrow(coll_lines) > 0) {
        line_parts[[length(line_parts) + 1]] <- coll_lines
      }
    }
    if (length(line_parts) == 0) {
      return(FALSE)
    }

    mst_lines <- do.call(rbind, line_parts)
    mst_lines <- safe_valid_sf(mst_lines, context_label = "mst_lines")
    if (nrow(mst_lines) == 0) {
      return(FALSE)
    }

    study_sf <- normalize_to_wgs84_sf(shape_file)
    bb <- sf::st_bbox(study_sf)
    xmin <- floor(as.numeric(bb["xmin"]) / grid_res) * grid_res
    xmax <- ceiling(as.numeric(bb["xmax"]) / grid_res) * grid_res
    ymin <- floor(as.numeric(bb["ymin"]) / grid_res) * grid_res
    ymax <- ceiling(as.numeric(bb["ymax"]) / grid_res) * grid_res

    grid_geom <- sf::st_make_grid(
      sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = sf::st_crs(study_sf))),
      cellsize = c(grid_res, grid_res),
      square = TRUE
    )
    if (length(grid_geom) == 0) {
      return(FALSE)
    }

    grid_sf <- sf::st_sf(cell_id = seq_along(grid_geom), geometry = grid_geom, crs = sf::st_crs(study_sf))
    line_union <- sf::st_union(sf::st_geometry(mst_lines))
    hit_idx <- lengths(sf::st_intersects(grid_sf, line_union)) > 0
    hit_grid <- grid_sf[hit_idx, , drop = FALSE]

    mst_target <- file.path(out_dir, "mst_ancterminal_mintreeall.shp")
    suppressWarnings(sf::st_write(mst_lines, mst_target, delete_layer = TRUE, quiet = TRUE))

    if (nrow(hit_grid) > 0) {
      grid_target <- file.path(out_dir, "GRIDS_mintreeall_ancterminal.shp")
      suppressWarnings(sf::st_write(hit_grid, grid_target, delete_layer = TRUE, quiet = TRUE))

      r_template <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, resolution = grid_res, crs = "EPSG:4326")
      hit_vect <- terra::vect(hit_grid)
      r_presence <- terra::rasterize(hit_vect, r_template, field = 1, background = NA)
      terra::writeRaster(r_presence, file.path(out_dir, "presence_mst_mintreeall_ancterminal.tif"), overwrite = TRUE)
    }

    TRUE
  }

  filter_visual_outputs <- function(output_files) {
    if (!is.null(output_files$shapefiles) && length(output_files$shapefiles) > 0) {
      keep_shp <- vapply(names(output_files$shapefiles), function(nm) {
        base_name <- tools::file_path_sans_ext(basename(nm))
        !is_helper_points_layer(base_name) && !is_mst_internal_helper_layer(base_name)
      }, logical(1))
      output_files$shapefiles <- output_files$shapefiles[keep_shp]
    }

    if (!is.null(output_files$shapefile_info) && length(output_files$shapefile_info) > 0) {
      kept_keys <- unique(c(names(output_files$shapefiles), names(output_files$rasters)))
      output_files$shapefile_info <- output_files$shapefile_info[intersect(names(output_files$shapefile_info), kept_keys)]
    }

    output_files
  }

  ensure_species_palette <- function(extra_species = character(0)) {
    species_from_occ <- if (!is.null(data_store$occurrence) && nrow(data_store$occurrence) > 0) {
      unique(as.character(data_store$occurrence$spp))
    } else {
      character(0)
    }

    species_from_layers <- character(0)
    if (!is.null(data_store$loaded_shapefiles) && length(data_store$loaded_shapefiles) > 0) {
      species_from_layers <- unique(vapply(seq_along(data_store$loaded_shapefiles), function(i) {
        layer_key <- names(data_store$loaded_shapefiles)[i]
        layer_info <- data_store$shapefile_info[[layer_key]]
        layer_species_name(layer_info, layer_key)
      }, character(1)))
    }

    species_list <- unique(c(species_from_occ, species_from_layers, as.character(extra_species)))
    species_list <- species_list[!is.na(species_list) & nzchar(species_list)]
    if (length(species_list) == 0) {
      return(invisible(NULL))
    }

    current_colors <- data_store$species_colors
    if (is.null(current_colors) || length(current_colors) == 0) {
      current_colors <- setNames(rep(NA_character_, length(species_list)), species_list)
    }

    missing_names <- setdiff(species_list, names(current_colors))
    if (length(missing_names) > 0) {
      current_colors <- c(current_colors, setNames(rep(NA_character_, length(missing_names)), missing_names))
    }

    palette_slice <- current_colors[species_list]
    na_idx <- is.na(palette_slice)
    if (any(na_idx)) {
      palette_slice[na_idx] <- grDevices::hcl.colors(sum(na_idx), "Dark 3")
      current_colors[species_list] <- palette_slice
    }

    data_store$species_colors <- current_colors
    invisible(NULL)
  }

  get_species_color_map <- function(species_values) {
    species_values <- unique(as.character(species_values))
    species_values <- species_values[!is.na(species_values) & nzchar(species_values)]
    if (length(species_values) == 0) {
      return(setNames(character(0), character(0)))
    }

    base_colors <- data_store$species_colors[species_values]
    missing_idx <- is.na(base_colors) | !nzchar(base_colors)
    if (any(missing_idx)) {
      fallback <- grDevices::hcl.colors(sum(missing_idx), "Dark 3")
      base_colors[missing_idx] <- fallback
    }
    setNames(as.character(base_colors), species_values)
  }

  clear_visual_output_files <- function(method = NULL) {
    method_dir_map <- c(
      buffer = "out_buffers",
      convex_hull = "out_MCP",
      mst = "out_MST",
      occurrence_only = "out_irregular_bins",
      irregular_bins = "out_irregular_bins"
    )

    target_dirs <- if (is.null(method)) {
      unique(c(unname(method_dir_map), "out_grid"))
    } else {
      dir_name <- unname(method_dir_map[method])
      if (is.na(dir_name) || !nzchar(dir_name)) {
        character(0)
      } else {
        unique(c(dir_name, "out_grid"))
      }
    }

    removed <- character(0)
    for (dir_name in target_dirs) {
      dir_path <- file.path(getwd(), dir_name)
      if (!dir.exists(dir_path)) next
      files <- list.files(
        dir_path,
        pattern = "\\.(shp|shx|dbf|prj|cpg|qix|tif|asc|grd)$",
        full.names = TRUE,
        ignore.case = TRUE
      )
      if (length(files) == 0) next
      keep <- basename(files) %in% c(".gitkeep")
      files <- files[!keep]
      if (length(files) == 0) next
      ok <- file.remove(files)
      removed <- c(removed, files[ok])
    }

    removed
  }

  startup_removed_files <- clear_visual_output_files(method = NULL)
  if (length(startup_removed_files) > 0) {
    message("Startup cleanup removed ", length(startup_removed_files), " old visual output files.")
  }

  make_area_code <- function(idx) {
    if (!is.finite(idx) || idx < 1) return(NA_character_)
    idx <- as.integer(idx)
    out <- ""
    while (idx > 0) {
      r <- (idx - 1) %% 26
      out <- paste0(LETTERS[r + 1], out)
      idx <- (idx - 1) %/% 26
    }
    out
  }

  abbreviate_area_labels <- function(labels, label_lookup = NULL) {
    labels <- as.character(labels)
    keep_root <- labels == "ROOT"
    original_non_root <- unique(labels[!keep_root])
    original_non_root <- original_non_root[!is.na(original_non_root) & nzchar(original_non_root)]

    if (length(original_non_root) == 0) {
      return(list(labels = labels, mapping = data.frame(original = character(0), code = character(0), stringsAsFactors = FALSE)))
    }

    codes <- vapply(seq_along(original_non_root), make_area_code, character(1))
    code_map <- setNames(codes, original_non_root)
    new_labels <- labels
    idx_non_root <- which(!keep_root)
    mapped <- unname(code_map[new_labels[idx_non_root]])
    mapped[is.na(mapped)] <- new_labels[idx_non_root][is.na(mapped)]
    new_labels[idx_non_root] <- mapped

    mapping_df <- data.frame(original = original_non_root, code = codes, stringsAsFactors = FALSE)
    if (!is.null(label_lookup) && length(label_lookup) > 0) {
      friendly <- unname(label_lookup[mapping_df$original])
      mapping_df$label <- ifelse(is.na(friendly) | !nzchar(friendly), mapping_df$original, friendly)
      mapping_df <- mapping_df[, c("code", "original", "label")]
    } else {
      mapping_df <- mapping_df[, c("code", "original")]
    }

    list(
      labels = new_labels,
      mapping = mapping_df
    )
  }

  get_area_label_lookup <- function() {
    bm <- data_store$irregular_bins_bin_id_mapping
    if (is.null(bm) || !is.data.frame(bm) || nrow(bm) == 0) {
      return(setNames(character(0), character(0)))
    }

    cn <- names(bm)
    id_col <- cn[grep("number|id", cn, ignore.case = TRUE)][1]
    name_col <- cn[grep("name|label|bin", cn, ignore.case = TRUE)][1]
    if (is.na(id_col) || is.na(name_col) || !nzchar(id_col) || !nzchar(name_col)) {
      return(setNames(character(0), character(0)))
    }

    ids <- as.character(bm[[id_col]])
    labs <- as.character(bm[[name_col]])
    ok <- !is.na(ids) & nzchar(ids) & !is.na(labs) & nzchar(labs)
    setNames(labs[ok], ids[ok])
  }

  get_pae_pce_function <- function() {
    resolve_from_file <- function(path) {
      if (!file.exists(path)) return(NULL)
      temp_env <- new.env(parent = globalenv())
      tryCatch(sys.source(path, envir = temp_env, chdir = TRUE), error = function(e) NULL)
      get0("pae_pce", envir = temp_env, mode = "function", inherits = FALSE)
    }

    candidate_paths <- unique(c(
      file.path("..", "..", "R", "pae_pce.R"),
      file.path("..", "R", "pae_pce.R"),
      file.path(getwd(), "R", "pae_pce.R"),
      file.path(getwd(), "biogeoshiny", "R", "pae_pce.R"),
      file.path(getwd(), "biogeoshiny_improved_0.1.0", "biogeoshiny", "R", "pae_pce.R")
    ))

    for (candidate in candidate_paths) {
      candidate_fun <- resolve_from_file(candidate)
      if (is.function(candidate_fun)) {
        return(candidate_fun)
      }
    }

    pae_fun <- get0("pae_pce", mode = "function", inherits = TRUE)
    if (is.function(pae_fun)) {
      return(pae_fun)
    }

    stop("PAE-PCE function could not be resolved in the current runtime.")
  }
  
  # ===== DATA INPUT TAB =====
  
  # Function to standardize and validate occurrence data
  standardize_occurrence_data <- function(data) {
    # Convert all column names to lowercase
    names(data) <- tolower(names(data))
    
    # Check for required columns (case-insensitive)
    required_cols <- c("spp", "long", "lat")
    available_cols <- tolower(names(data))
    
    if (!all(required_cols %in% available_cols)) {
      missing <- setdiff(required_cols, available_cols)
      stop(paste("Missing columns:", paste(missing, collapse=", ")))
    }
    
    # Select and reorder columns
    data <- data[, required_cols]
    
    # Ensure correct data types
    data$spp <- as.character(data$spp)
    data$long <- as.numeric(data$long)
    data$lat <- as.numeric(data$lat)
    
    # Remove rows with NA values
    data <- na.omit(data)
    
    if (nrow(data) == 0) {
      stop("No valid data rows after removing NAs")
    }
    
    return(data)
  }

  read_occurrence_table <- function(file_info) {
    if (is.null(file_info)) {
      stop("No file selected")
    }

    ext <- tolower(tools::file_ext(file_info$name %||% ""))
    path <- file_info$datapath

    read_attempts <- list()
    if (ext == "csv") {
      read_attempts <- list(
        function() utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        function() utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
        function() utils::read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
      )
    } else {
      read_attempts <- list(
        function() utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
        function() utils::read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE),
        function() utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
      )
    }

    last_err <- NULL
    for (reader in read_attempts) {
      out <- tryCatch(reader(), error = function(e) {
        last_err <<- e
        NULL
      })
      if (!is.null(out) && ncol(out) >= 3) {
        return(out)
      }
    }

    if (!is.null(last_err)) {
      stop(last_err$message)
    }
    stop("Could not read occurrence table. Please provide a valid CSV/TXT with header.")
  }
  
  observeEvent(input$load_occurrence, {
    tryCatch({
      if (is.null(input$occurrence_file)) {
        output$data_status <- renderPrint({
          cat("Please select a CSV/TXT file\n")
        })
      } else {
        data <- read_occurrence_table(input$occurrence_file)
        
        # Standardize and validate the data
        data <- standardize_occurrence_data(data)
        data_store$occurrence <- data
        data_store$analysis_occurrence <- NULL
        data_store$analysis_tree <- NULL
        data_store$problematic_taxa <- NULL
        data_store$singleton_taxa <- NULL
        data_store$doubleton_taxa <- NULL
        data_store$duplicate_rows <- NULL
        data_store$duplicate_taxa <- NULL
        data_store$singleton_result <- NULL
        data_store$pres_abs <- NULL
        data_store$pres_abs_regular <- NULL
        data_store$pres_abs_irregular <- NULL
        data_store$mst_pres_abs <- NULL
        data_store$mst_context <- NULL
                data_store$loaded_shapefiles <- list()
                data_store$loaded_rasters <- list()
                data_store$shapefile_info <- list()
                data_store$regular_grid_all_sf <- NULL
                data_store$regular_grid_presence_sf <- NULL
                data_store$regular_grid_presence_points_sf <- NULL
                data_store$regular_grid_presence_extrap_sf <- NULL
                data_store$regular_grid_presence_points_by_taxon <- list()
                data_store$regular_grid_presence_extrap_by_taxon <- list()
                data_store$regular_grid_presence_extrap_label <- NULL
                data_store$regular_grid_presence_by_taxon <- list()
                data_store$regular_grid_res <- NULL
                data_store$last_extrap_started_at <- NULL
                data_store$last_extrap_taxa <- character(0)
                data_store$extrap_log <- character(0)

        if (!is.null(data_store$tree)) {
          data_store$tree <- NULL
          show_tree_node_labels(FALSE)
          updateCheckboxInput(session, "use_clade_filter", value = FALSE)
          updateSelectizeInput(session, "clade_node_id", choices = character(0), selected = character(0))
          output$tree_load_status <- renderPrint({
            cat("Previous tree was cleared to avoid mismatch with the newly loaded occurrence dataset.\n")
            cat("Please reload the phylogenetic tree if you need tree-based analysis.\n")
          })
        }
        
        output$data_status <- renderPrint({
          cat("✓ Occurrence data loaded successfully!\n")
          cat("Rows:", nrow(data_store$occurrence), "\n")
          cat("Species:", length(unique(data_store$occurrence$spp)), "\n")
          cat("Tree state: reset (reload tree if needed)\n")
        })
      }
    }, error = function(e) {
      output$data_status <- renderPrint({
        cat("✗ Error loading data:\n", e$message, "\n")
      })
    })
  })
  
  observeEvent(input$load_tree, {
    tryCatch({
      if (is.null(input$tree_file)) {
        output$tree_load_status <- renderPrint({
          cat("Please select a tree file\n")
        })
      } else {
        # Try to read the tree file with different formats
        tree <- NULL
        file_path <- input$tree_file$datapath
        file_ext <- tolower(tools::file_ext(input$tree_file$name))
        
        # Try different reading methods
        if (file_ext %in% c("tre", "nwk", "newick", "txt")) {
          # Try as Newick format
          tryCatch({
            tree <- ape::read.tree(file_path)
          }, error = function(e) {
            NULL
          })
        } else if (file_ext %in% c("nex", "nexus")) {
          # Try as Nexus format
          tryCatch({
            tree <- ape::read.nexus(file_path)
          }, error = function(e) {
            NULL
          })
        }
        
        # If still NULL, try read.tree as fallback
        if (is.null(tree)) {
          tree <- ape::read.tree(file_path)
        }
        
        if (inherits(tree, "multiPhylo")) {
          if (length(tree) == 0) {
            stop("The uploaded tree file has no trees.")
          }
          tree <- tree[[1]]
        }

        data_store$tree <- tree
        shiny::updateCheckboxInput(session, "show_node_labels", value = FALSE)

        output$tree_load_status <- renderPrint({
          cat("✓ Tree loaded successfully!\n")
          cat("Number of tips:", length(tree$tip.label), "\n")
          cat("Has branch lengths:", !is.null(tree$edge.length), "\n")
        })
      }
    }, error = function(e) {
      output$tree_load_status <- renderPrint({
        cat("✗ Error loading tree:\n", e$message, "\n")
      })
    })
  })
  
  # ===== TREE VALIDATION TAB =====
  
  observeEvent(input$validate_tree, {
    if (is.null(data_store$tree)) {
      output$tree_validation_output <- renderPrint({
        cat("Error: Please load a tree first\n")
      })
    } else {
      tree <- data_store$tree
      shiny::updateCheckboxInput(session, "show_node_labels", value = TRUE)
      output$tree_validation_output <- renderPrint({
        n_tips <- ape::Ntip(tree)
        n_nodes <- tree$Nnode %||% 0
        cat("✓ Tree inspection:\n")
        cat("Number of tips:", length(tree$tip.label), "\n")
        cat("Number of internal nodes:", n_nodes, "\n")
        if (n_nodes > 0) {
          cat("Internal node IDs:", (n_tips + 1), "to", (n_tips + n_nodes), "\n")
        }
        cat("Has branch lengths:", !is.null(tree$edge.length), "\n")
        if (!is.null(tree$edge.length)) {
          cat("Branch length range:", min(tree$edge.length), "-", max(tree$edge.length), "\n")
        } else {
          cat("Observation: branch lengths are missing (this does not block visualization).\n")
        }
      })
    }
  })
  
  output$tree_plot <- ggiraph::renderGirafe({
    if (is.null(data_store$tree)) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No tree loaded", size = 6) +
        ggplot2::theme_void()
      return(ggiraph::girafe(ggobj = p, options = list(
        ggiraph::opts_hover(css = "opacity:0.7;"),
        ggiraph::opts_sizing(rescale = TRUE, width = 1)
      )))
    }
    
    tree <- data_store$tree
    n_tips <- ape::Ntip(tree)
    n_nodes <- tree$Nnode %||% 0
    show_labels <- isTRUE(input$show_node_labels)
    
    node_size <- if (n_tips > 200) 2 else if (n_tips > 100) 2.5 else 3
    text_size <- if (n_tips > 200) 2.5 else if (n_tips > 100) 3 else 3.5
    
    p <- ggtree::ggtree(tree, layout = "rectangular", size = 0.5, color = "#1f1f1f") +
      ggplot2::ggtitle("Phylogenetic Tree (Interactive)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    if (show_labels && n_nodes > 0) {
      internal_nodes <- (n_tips + 1):(n_tips + n_nodes)
      
      p <- p +
        ggiraph::geom_point_interactive(
          ggplot2::aes(
            x = x, y = y,
            tooltip = paste("Node ID:", node, "\nClick to see terminal taxa"),
            data_id = node
          ),
          data = function(d) d[d$node %in% internal_nodes, ],
          size = node_size,
          color = "#a61c1c",
          fill = "#a61c1c",
          shape = 21,
          stroke = 1.5,
          hover_nearest = TRUE
        ) +
        ggiraph::geom_text_interactive(
          ggplot2::aes(
            x = x, y = y,
            label = node,
            tooltip = paste("Node ID:", node, "\nClick to see terminal taxa"),
            data_id = node
          ),
          data = function(d) d[d$node %in% internal_nodes, ],
          size = text_size * 0.7,
          color = "white",
          fontface = "bold",
          hjust = 0.5,
          vjust = 0.5
        )
    }
    
    ggiraph::girafe(
      ggobj = p,
      options = list(
        ggiraph::opts_hover(css = "opacity:0.8; stroke-width:2px;"),
        ggiraph::opts_selection(type = "single", css = "stroke-width:3px;"),
        ggiraph::opts_zoom(min = 0.5, max = 5),
        ggiraph::opts_toolbar(position = "topright"),
        ggiraph::opts_sizing(rescale = TRUE, width = 1)
      )
    )
  })
  
  # Handle node selection to show terminal taxa
  observeEvent(input$tree_plot_selected, {
    if (is.null(data_store$tree) || is.null(input$tree_plot_selected) || length(input$tree_plot_selected) == 0) return()
    
    tree <- data_store$tree
    clicked_node <- as.numeric(input$tree_plot_selected[1])
    
    if (is.na(clicked_node)) return()
    
    tryCatch({
      # Get descendants of the clicked node
      descendants <- ape::extract.clade(tree, clicked_node)
      terminal_taxa <- descendants$tip.label
      
      # Create output for modal
      output$node_taxa_list <- renderPrint({
        cat("Node ID:", clicked_node, "\n")
        cat("Number of terminal taxa:", length(terminal_taxa), "\n\n")
        cat("Terminal taxa:\n")
        cat(paste("-", terminal_taxa, collapse = "\n"), "\n")
      })
      
      # Show modal
      shiny::showModal(shiny::modalDialog(
        title = paste("Terminal Taxa in Node", clicked_node),
        div(
          verbatimTextOutput("node_taxa_list"),
          style = "max-height: 500px; overflow-y: auto; font-size: 12px;"
        ),
        easyClose = TRUE,
        footer = shiny::modalButton("Close")
      ))
    }, error = function(e) {
      shiny::showModal(shiny::modalDialog(
        title = "Error",
        paste("Could not extract clade:", e$message),
        easyClose = TRUE,
        footer = shiny::modalButton("Close")
      ))
    })
  })
  
  # Handle clade node selection to show terminal taxa
  observeEvent(input$clade_tree_plot_selected, {
    if (is.null(data_store$tree) || is.null(input$clade_node_id) || input$clade_node_id == "") return()
    if (is.null(input$clade_tree_plot_selected) || length(input$clade_tree_plot_selected) == 0) return()
    
    tree <- data_store$tree
    clade_root <- as.numeric(input$clade_node_id)
    clicked_node_in_clade <- as.numeric(input$clade_tree_plot_selected[1])
    
    if (is.na(clicked_node_in_clade)) return()
    
    tryCatch({
      # Extract the clade first
      clade_tree <- ape::extract.clade(tree, clade_root)
      clade_n_tips <- ape::Ntip(clade_tree)
      
      # Map the clicked node back to the original tree numbering
      # The clade_tree has its own numbering, but we need to show the original tree node numbers
      if (clicked_node_in_clade <= clade_n_tips) {
        # It's a tip - show just that tip
        terminal_taxa <- clade_tree$tip.label[clicked_node_in_clade]
        node_label <- paste("Tip:", terminal_taxa)
        original_node_id <- "(terminal taxon)"
      } else {
        # It's an internal node in the clade - get its descendants
        descendants <- ape::extract.clade(clade_tree, clicked_node_in_clade)
        terminal_taxa <- descendants$tip.label
        # For internal nodes, we keep the clade's node ID since it's what the user sees
        original_node_id <- clicked_node_in_clade
        node_label <- paste("Node ID:", original_node_id)
      }
      
      # Create output for modal
      output$clade_node_taxa_list <- renderPrint({
        cat(node_label, "\n")
        cat("Number of terminal taxa:", length(terminal_taxa), "\n\n")
        cat("Terminal taxa:\n")
        cat(paste("-", terminal_taxa, collapse = "\n"), "\n")
      })
      
      # Show modal
      shiny::showModal(shiny::modalDialog(
        title = paste("Terminal Taxa in Node", original_node_id),
        div(
          verbatimTextOutput("clade_node_taxa_list"),
          style = "max-height: 500px; overflow-y: auto; font-size: 12px;"
        ),
        easyClose = TRUE,
        footer = shiny::modalButton("Close")
      ))
    }, error = function(e) {
      shiny::showModal(shiny::modalDialog(
        title = "Error",
        paste("Could not extract clade:", e$message),
        easyClose = TRUE,
        footer = shiny::modalButton("Close")
      ))
    })
  })
  
  # Render clade tree when a clade node is selected (with reactive dependency on checkbox)
  output$clade_tree_plot <- ggiraph::renderGirafe({
    # Trigger reactivity when checkbox changes
    input$show_clade_node_labels
    
    if (is.null(data_store$tree) || is.null(input$clade_node_id) || input$clade_node_id == "") {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Select a node to preview", size = 5) +
        ggplot2::theme_void()
      return(ggiraph::girafe(ggobj = p, options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1)
      )))
    }
    
    tree <- data_store$tree
    clade_node <- as.numeric(input$clade_node_id)
    
    tryCatch({
      # Extract the clade
      clade_tree <- ape::extract.clade(tree, clade_node)
      
      # Ensure the tree has tip labels
      if (is.null(clade_tree$tip.label)) {
        clade_tree$tip.label <- paste0("Tip_", 1:ape::Ntip(clade_tree))
      }
      
      # Build ggtree plot for the clade with tip labels
      clade_n_tips <- ape::Ntip(clade_tree)
      clade_n_nodes <- clade_tree$Nnode %||% 0
      show_clade_labels <- isTRUE(input$show_clade_node_labels)
      
      clade_node_size <- if (clade_n_tips > 200) 2 else if (clade_n_tips > 100) 2.5 else 3
      clade_text_size <- if (clade_n_tips > 200) 2.5 else if (clade_n_tips > 100) 3 else 3.5
      
      p <- ggtree::ggtree(clade_tree, layout = "rectangular", size = 0.5, color = "#1f1f1f") +
        ggtree::geom_tiplab(size = 1.8, hjust = -0.05, color = "#333333") +
        ggplot2::ggtitle(paste("Clade rooted at Node", clade_node)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(10, 80, 10, 10)
        )
      
      # Add internal node labels if checkbox is checked
      if (show_clade_labels && clade_n_nodes > 0) {
        clade_internal_nodes <- (clade_n_tips + 1):(clade_n_tips + clade_n_nodes)
        
        # Create a mapping from clade node numbers to original tree node numbers
        # This is done by matching tips and working backwards
        clade_tip_names <- clade_tree$tip.label
        original_tip_indices <- match(clade_tip_names, tree$tip.label)
        
        p <- p +
          ggiraph::geom_point_interactive(
            ggplot2::aes(
              x = x, y = y,
              tooltip = paste("Node ID:", node, "\nClick to see terminal taxa"),
              data_id = node
            ),
            data = function(d) d[d$node %in% clade_internal_nodes, ],
            size = clade_node_size,
            color = "#a61c1c",
            fill = "#a61c1c",
            shape = 21,
            stroke = 1.5,
            hover_nearest = TRUE
          ) +
          ggiraph::geom_text_interactive(
            ggplot2::aes(
              x = x, y = y,
              label = node,
              tooltip = paste("Node ID:", node, "\nClick to see terminal taxa"),
              data_id = node
            ),
            data = function(d) d[d$node %in% clade_internal_nodes, ],
            size = clade_text_size * 0.7,
            color = "white",
            fontface = "bold",
            hjust = 0.5,
            vjust = 0.5
          )
      }
      
      # Extend x-axis to accommodate tip labels
      tryCatch({
        x_range <- ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
        p <- p + ggplot2::xlim(x_range[1], x_range[2] * 1.5)
      }, error = function(e) {
        # If xlim fails, just continue without it
        NULL
      })
      
      # Create interactive girafe
      ggiraph::girafe(
        ggobj = p,
        options = list(
          ggiraph::opts_hover(css = "opacity:0.8;"),
          ggiraph::opts_selection(type = "single", css = "stroke-width:3px;"),
          ggiraph::opts_zoom(min = 0.5, max = 5),
          ggiraph::opts_toolbar(position = "topright"),
          ggiraph::opts_sizing(rescale = TRUE, width = 1)
        )
      )
    }, error = function(e) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Error extracting clade", size = 5, color = "red") +
        ggplot2::theme_void()
      ggiraph::girafe(ggobj = p, options = list(
        ggiraph::opts_sizing(rescale = TRUE, width = 1)
      ))
    })
  })
  
  # ===== DATA PREPROCESSING TAB =====

  output$problematic_recommendation <- renderPrint({
    method <- input$extrap_method %||% "buffer"
    if (identical(method, "buffer") || identical(method, "occurrence_only")) {
      cat("Method guidance:\n")
      cat("- Buffer/Occurrence-only can use singletons and doubletons.\n")
      cat("- If you want to keep all taxa, run detection but skip removal.\n")
    } else if (identical(method, "mst")) {
      cat("Method guidance:\n")
      cat("- MST requires >= 2 points per taxon.\n")
      cat("- Recommended strategy: remove only singletons.\n")
    } else if (identical(method, "convex_hull")) {
      cat("Method guidance:\n")
      cat("- Convex hull needs >= 3 points per taxon.\n")
      cat("- Recommended strategy: remove singletons + doubletons.\n")
    }
  })
  
  observeEvent(input$detect_problems, {
    if (is.null(data_store$occurrence)) {
      output$problem_detection_output <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
    } else {
      tryCatch({
        occ_data <- data_store$occurrence
        species_counts <- table(occ_data$spp)
        singletons <- names(species_counts[species_counts == 1])
        doubletons <- names(species_counts[species_counts == 2])
        dup_idx <- duplicated(occ_data[, c("spp", "long", "lat"), drop = FALSE])
        duplicate_taxa <- sort(unique(occ_data$spp[dup_idx]))
        duplicate_n <- sum(dup_idx)
        
        output$problem_detection_output <- renderPrint({
          cat("=== Problematic Taxa Detection ===\n\n")
          if (length(singletons) > 0) {
            cat("Singletons (1 point):\n")
            for (sp in singletons) cat("  -", sp, "\n")
            cat("\n")
          }
          if (length(doubletons) > 0) {
            cat("Doubletons (2 points):\n")
            for (sp in doubletons) cat("  -", sp, "\n")
            cat("\n")
          }
          if (duplicate_n > 0) {
            cat("Duplicate occurrence records (same spp + long + lat):", duplicate_n, "\n")
            cat("Affected taxa:\n")
            for (sp in duplicate_taxa) cat("  -", sp, "\n")
            cat("\n")
          }

          if (!is.null(data_store$tree)) {
            taxa_data <- unique(as.character(occ_data$spp))
            taxa_tree <- data_store$tree$tip.label
            taxa_data_only <- setdiff(taxa_data, taxa_tree)
            taxa_tree_only <- setdiff(taxa_tree, taxa_data)

            if (length(taxa_data_only) > 0) {
              cat("Taxa in data but not in tree:", length(taxa_data_only), "\n")
            }
            if (length(taxa_tree_only) > 0) {
              cat("Taxa in tree but not in data:", length(taxa_tree_only), "\n")
            }
            if (length(taxa_data_only) > 0 || length(taxa_tree_only) > 0) {
              cat("These mismatches will be harmonized when removing taxa.\n\n")
            }
          }

          problematic_union <- unique(c(singletons, doubletons, duplicate_taxa))
          if (length(singletons) == 0 && length(doubletons) == 0 && duplicate_n == 0) {
            cat("✓ No problematic taxa detected!\n")
          } else {
            cat("Total problematic taxa:", length(problematic_union), "\n")
          }
        })
        data_store$singleton_taxa <- singletons
        data_store$doubleton_taxa <- doubletons
        data_store$problematic_taxa <- unique(c(singletons, doubletons))
        data_store$duplicate_rows <- dup_idx
        data_store$duplicate_taxa <- duplicate_taxa
      }, error = function(e) {
        output$problem_detection_output <- renderPrint({
          cat("Error:\n", e$message, "\n")
        })
      })
    }
  })
  
  observeEvent(input$remove_problems, {
    if (is.null(data_store$occurrence)) {
      output$removal_output <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
    } else {
      tryCatch({
        occ_data <- data_store$occurrence
        tree <- data_store$tree
        removal_mode <- input$problematic_removal_mode %||% "singletons_doubletons"
        singletons <- data_store$singleton_taxa %||% character(0)
        doubletons <- data_store$doubleton_taxa %||% character(0)

        problematic <- switch(
          removal_mode,
          singletons_only = singletons,
          doubletons_only = doubletons,
          singletons_doubletons = unique(c(singletons, doubletons)),
          unique(c(singletons, doubletons))
        )

        original_count <- nrow(occ_data)
        occ_data_clean <- if (length(problematic) > 0) {
          occ_data[!(occ_data$spp %in% problematic), , drop = FALSE]
        } else {
          occ_data
        }
        removed_problematic_rows <- original_count - nrow(occ_data_clean)

        tree_tips_removed_problematic <- character(0)
        tree_tips_removed_tree_only <- character(0)
        taxa_data_only <- character(0)

        if (!is.null(tree) && nrow(occ_data_clean) > 0) {
          tree_clean <- tree

          tree_tips_removed_problematic <- intersect(problematic, tree_clean$tip.label)
          if (length(tree_tips_removed_problematic) > 0) {
            if (length(tree_tips_removed_problematic) >= length(tree_clean$tip.label)) {
              tree_clean <- NULL
            } else {
              tree_clean <- ape::drop.tip(tree_clean, tree_tips_removed_problematic)
            }
          }

          if (!is.null(tree_clean)) {
            taxa_data <- unique(as.character(occ_data_clean$spp))
            taxa_tree <- tree_clean$tip.label

            taxa_data_only <- setdiff(taxa_data, taxa_tree)
            if (length(taxa_data_only) > 0) {
              occ_data_clean <- occ_data_clean[!(occ_data_clean$spp %in% taxa_data_only), , drop = FALSE]
            }

            taxa_tree_only <- setdiff(tree_clean$tip.label, unique(as.character(occ_data_clean$spp)))
            tree_tips_removed_tree_only <- taxa_tree_only
            if (length(taxa_tree_only) > 0) {
              if (length(taxa_tree_only) >= length(tree_clean$tip.label)) {
                tree_clean <- NULL
              } else {
                tree_clean <- ape::drop.tip(tree_clean, taxa_tree_only)
              }
            }
          }

          if (nrow(occ_data_clean) == 0) {
            stop("All occurrences were removed after harmonizing data and tree. Choose a less restrictive option.")
          }
        } else {
          if (!is.null(tree)) {
            tree_tips_removed_problematic <- intersect(problematic, tree$tip.label)
            if (length(tree_tips_removed_problematic) > 0 && length(tree_tips_removed_problematic) < length(tree$tip.label)) {
              tree_clean <- ape::drop.tip(tree, tree_tips_removed_problematic)
            } else if (length(tree_tips_removed_problematic) >= length(tree$tip.label)) {
              tree_clean <- NULL
            } else {
              tree_clean <- tree
            }
          } else {
            tree_clean <- NULL
          }

          if (nrow(occ_data_clean) == 0) {
            stop("All occurrences were removed by the selected strategy. Choose a less restrictive option.")
          }
        }

        data_store$occurrence <- occ_data_clean
        data_store$tree <- tree_clean
        data_store$singleton_result <- NULL
        data_store$problematic_taxa <- NULL
        data_store$singleton_taxa <- NULL
        data_store$doubleton_taxa <- NULL
        data_store$duplicate_rows <- NULL
        data_store$duplicate_taxa <- NULL
        data_store$pres_abs <- NULL
        data_store$pres_abs_regular <- NULL
        data_store$pres_abs_irregular <- NULL
        data_store$mst_pres_abs <- NULL
        data_store$mst_context <- NULL

        output$removal_output <- renderPrint({
          cat("=== Removal Complete ===\n\n")
          cat("Selection:", removal_mode, "\n")
          cat("Removed taxa by selection:", if (length(problematic) == 0) "none" else paste(problematic, collapse = ", "), "\n")
          cat("Rows removed by selection:", removed_problematic_rows, "\n")
          cat("Data-only taxa removed during harmonization:", length(taxa_data_only), "\n")
          cat("Tree-only taxa removed during harmonization:", length(tree_tips_removed_tree_only), "\n")
          cat("Remaining rows:", nrow(occ_data_clean), "\n")
          if (!is.null(tree)) {
            cat("Tree tips removed by selection:", length(tree_tips_removed_problematic), "\n")
            cat("Tree tips remaining:", if (is.null(tree_clean)) 0 else length(tree_clean$tip.label), "\n")
          }
        })

        output$cleaned_data_summary <- renderPrint({
          cat("Cleaned Data Summary:\n")
          cat("Total occurrences:", nrow(occ_data_clean), "\n")
          cat("Total species:", length(unique(occ_data_clean$spp)), "\n")
          cat("Species list:", paste(unique(occ_data_clean$spp), collapse = ", "), "\n")
        })

        output$updated_tree_plot <- renderPlot({
          tree_display <- data_store$tree
          if (is.null(tree_display)) {
            plot(1, type = "n", main = "No tree available")
          } else {
            plot(tree_display, main = "Updated Phylogenetic Tree")
          }
        })
      }, error = function(e) {
        output$removal_output <- renderPrint({
          cat("Error:\n", e$message, "\n")
        })
      })
    }
  })

  observeEvent(input$remove_duplicates, {
    if (is.null(data_store$occurrence)) {
      output$duplicates_output <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
    } else {
      tryCatch({
        occ_data <- data_store$occurrence
        dup_idx <- duplicated(occ_data[, c("spp", "long", "lat"), drop = FALSE])
        duplicate_n <- sum(dup_idx)

        if (duplicate_n == 0) {
          output$duplicates_output <- renderPrint({
            cat("=== Duplicate Removal ===\n\n")
            cat("✓ No duplicate occurrence records were found.\n")
          })
          return()
        }

        duplicate_taxa <- sort(unique(occ_data$spp[dup_idx]))
        occ_clean <- occ_data[!dup_idx, , drop = FALSE]

        data_store$occurrence <- occ_clean
        data_store$problematic_taxa <- NULL
        data_store$singleton_taxa <- NULL
        data_store$doubleton_taxa <- NULL
        data_store$duplicate_rows <- NULL
        data_store$duplicate_taxa <- NULL
        data_store$singleton_result <- NULL

        # downstream results become stale after changing occurrence matrix
        data_store$pres_abs <- NULL
        data_store$pres_abs_regular <- NULL
        data_store$pres_abs_irregular <- NULL
        data_store$mst_pres_abs <- NULL
        data_store$mst_context <- NULL

        output$duplicates_output <- renderPrint({
          cat("=== Duplicate Removal Complete ===\n\n")
          cat("Removed duplicate occurrence records:", duplicate_n, "\n")
          cat("Affected taxa:", paste(duplicate_taxa, collapse = ", "), "\n")
          cat("Remaining occurrences:", nrow(occ_clean), "\n")
        })

        output$cleaned_data_summary <- renderPrint({
          cat("Cleaned Data Summary:\n")
          cat("Total occurrences:", nrow(occ_clean), "\n")
          cat("Total species:", length(unique(occ_clean$spp)), "\n")
          cat("Species list:", paste(unique(occ_clean$spp), collapse = ", "), "\n")
        })
      }, error = function(e) {
        output$duplicates_output <- renderPrint({
          cat("Error:\n", e$message, "\n")
        })
      })
    }
  })

  observeEvent(input$harmonize_tree_data, {
    if (is.null(data_store$occurrence)) {
      output$harmonization_output <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
    } else if (is.null(data_store$tree)) {
      output$harmonization_output <- renderPrint({
        cat("No tree loaded. Nothing to harmonize on tree side.\n")
      })
    } else {
      tryCatch({
        occ_data <- data_store$occurrence
        tree <- data_store$tree

        taxa_data <- unique(as.character(occ_data$spp))
        taxa_tree <- tree$tip.label

        taxa_data_only <- setdiff(taxa_data, taxa_tree)
        taxa_tree_only <- setdiff(taxa_tree, taxa_data)

        occ_data_clean <- if (length(taxa_data_only) > 0) {
          occ_data[!(occ_data$spp %in% taxa_data_only), , drop = FALSE]
        } else {
          occ_data
        }

        if (nrow(occ_data_clean) == 0) {
          stop("No overlapping taxa between occurrence data and tree. Harmonization would remove all occurrence rows.")
        }

        tree_clean <- tree
        if (length(taxa_tree_only) > 0) {
          if (length(taxa_tree_only) >= length(tree$tip.label)) {
            tree_clean <- NULL
          } else {
            tree_clean <- ape::drop.tip(tree, taxa_tree_only)
          }
        }

        data_store$occurrence <- occ_data_clean
        data_store$tree <- tree_clean
        data_store$singleton_result <- NULL
        data_store$problematic_taxa <- NULL
        data_store$singleton_taxa <- NULL
        data_store$doubleton_taxa <- NULL
        data_store$duplicate_rows <- NULL
        data_store$duplicate_taxa <- NULL
        data_store$pres_abs <- NULL
        data_store$pres_abs_regular <- NULL
        data_store$pres_abs_irregular <- NULL
        data_store$mst_pres_abs <- NULL
        data_store$mst_context <- NULL

        output$harmonization_output <- renderPrint({
          cat("=== Tree/Data Harmonization Complete ===\n\n")
          cat("Data-only taxa removed:", length(taxa_data_only), "\n")
          cat("Tree-only taxa removed:", length(taxa_tree_only), "\n")
          cat("Remaining occurrence rows:", nrow(occ_data_clean), "\n")
          cat("Remaining occurrence taxa:", length(unique(occ_data_clean$spp)), "\n")
          cat("Remaining tree tips:", if (is.null(tree_clean)) 0 else length(tree_clean$tip.label), "\n")
        })

        output$cleaned_data_summary <- renderPrint({
          cat("Cleaned Data Summary:\n")
          cat("Total occurrences:", nrow(occ_data_clean), "\n")
          cat("Total species:", length(unique(occ_data_clean$spp)), "\n")
          cat("Species list:", paste(unique(occ_data_clean$spp), collapse = ", "), "\n")
        })

        output$updated_tree_plot <- renderPlot({
          tree_display <- data_store$tree
          if (is.null(tree_display)) {
            plot(1, type = "n", main = "No tree available")
          } else {
            plot(tree_display, main = "Updated Phylogenetic Tree")
          }
        })
      }, error = function(e) {
        output$harmonization_output <- renderPrint({
          cat("Error:\n", e$message, "\n")
        })
      })
    }
  })
  
  output$cleaned_data_summary <- renderPrint({
    cat("Load and process data to see summary\n")
  })

  output$duplicates_output <- renderPrint({
    cat("Detect problematic taxa first, then click 'Remove Duplicates' if needed.\n")
  })

  output$harmonization_output <- renderPrint({
    cat("Click 'Harmonize Tree <-> Data' to reconcile taxa between occurrence data and tree.\n")
  })
  
  output$updated_tree_plot <- renderPlot({
    plot(1, type = "n", main = "No tree loaded")
  })

  output$download_pruned_tree <- downloadHandler(
    filename = function() {
      tr <- data_store$tree
      stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (is.null(tr)) {
        return(paste0("pruned_tree_unavailable_", stamp, ".txt"))
      }
      paste0("pruned_tree_", stamp, ".newick")
    },
    content = function(file) {
      tr <- data_store$tree
      if (is.null(tr)) {
        writeLines(
          c(
            "No pruned tree is currently available.",
            "Load a tree and run a pruning step in Data Preprocessing before downloading."
          ),
          con = file
        )
        return(invisible(NULL))
      }

      if (inherits(tr, "multiPhylo")) {
        if (length(tr) == 0) {
          writeLines("Tree object is empty.", con = file)
          return(invisible(NULL))
        }
        tr <- tr[[1]]
      }

      ape::write.tree(tr, file = file)
    }
  )
  

  # ===== RANGE EXTRAPOLATION TAB =====
  
  # Function to load shapefile from multiple files
  load_shapefile_from_files <- function(file_df) {
    tryCatch({
      # Find the .shp file (case-insensitive)
      shp_idx <- grep("\\.shp$", file_df$name, ignore.case = TRUE)
      if (length(shp_idx) == 0) {
        stop("No .shp file found. Please upload the .shp file.")
      }
      
      # Get base name from the .shp file (remove extension)
      shp_basename <- file_df$name[shp_idx[1]]
      shp_name <- sub("\\.[^.]*$", "", shp_basename)  # Remove any extension
      
      # Create a unique temporary directory for this shapefile
      temp_dir <- file.path(tempdir(), paste0("shp_", format(Sys.time(), "%s")))
      dir.create(temp_dir, showWarnings = FALSE)
      
      # Copy all related files to temp directory with same base name
      # Map of required extensions
      required_exts <- c(".shp", ".shx", ".dbf")
      optional_exts <- c(".prj", ".cpg")
      all_exts <- c(required_exts, optional_exts)
      
      files_copied <- list()
      for (i in seq_len(nrow(file_df))) {
        # Get extension (case-insensitive)
        file_name <- file_df$name[i]
        file_ext <- tolower(sub("^.*\\.", ".", file_name))
        
        # Check if this is a shapefile component
        if (file_ext %in% all_exts) {
          # Copy with lowercase extension
          new_name <- paste0(shp_name, file_ext)
          src_path <- file_df$datapath[i]
          dst_path <- file.path(temp_dir, new_name)
          
          if (file.copy(src_path, dst_path, overwrite = TRUE)) {
            files_copied[[file_ext]] <- TRUE
          }
        }
      }
      
      # Check if all required files were copied
      missing <- setdiff(required_exts, names(files_copied))
      if (length(missing) > 0) {
        stop(paste("Missing required files:", paste(missing, collapse = ", ")))
      }
      
      # Load the shapefile
      shp_path <- file.path(temp_dir, paste0(shp_name, ".shp"))
      
      if (!file.exists(shp_path)) {
        stop(paste("Shapefile not found at:", shp_path))
      }
      
      # Try loading with terra first
      if (requireNamespace("terra", quietly = TRUE)) {
        shape <- terra::vect(shp_path)
        return(shape)
      } else if (requireNamespace("sf", quietly = TRUE)) {
        shape <- sf::st_read(shp_path, quiet = TRUE)
        return(shape)
      } else if (requireNamespace("raster", quietly = TRUE)) {
        shape <- raster::shapefile(shp_path)
        return(shape)
      } else {
        stop("No spatial package available (terra, sf, or raster)")
      }
    }, error = function(e) {
      stop(paste("Error loading shapefile:", e$message))
    })
  }

  available_id_columns_from_upload <- function(file_df) {
    shp <- load_shapefile_from_files(file_df)
    shp_sf <- normalize_to_wgs84_sf(shp)
    cols <- setdiff(names(shp_sf), "geometry")
    if (length(cols) == 0) return(character(0))
    cols[sapply(cols, function(cc) sum(!is.na(shp_sf[[cc]])) > 0)]
  }

  observeEvent(input$irregular_bins_shapefile, {
    if (is.null(input$irregular_bins_shapefile) || nrow(input$irregular_bins_shapefile) == 0) {
      updateSelectInput(session, "irregular_bin_id_column", choices = c("(Auto-detect)" = ""), selected = "")
      return()
    }
    cols <- tryCatch(available_id_columns_from_upload(input$irregular_bins_shapefile), error = function(e) character(0))
    choices <- c("(Auto-detect)" = "", stats::setNames(cols, cols))
    updateSelectInput(session, "irregular_bin_id_column", choices = choices, selected = "")
  })

  observeEvent(input$irregular_richness_shapefile, {
    if (is.null(input$irregular_richness_shapefile) || nrow(input$irregular_richness_shapefile) == 0) {
      updateSelectInput(session, "irregular_richness_id_column", choices = c("(Auto-detect)" = ""), selected = "")
      return()
    }
    cols <- tryCatch(available_id_columns_from_upload(input$irregular_richness_shapefile), error = function(e) character(0))
    choices <- c("(Auto-detect)" = "", stats::setNames(cols, cols))
    updateSelectInput(session, "irregular_richness_id_column", choices = choices, selected = "")
  })

  observeEvent(input$points_irregular_bins_shapefile, {
    if (is.null(input$points_irregular_bins_shapefile) || nrow(input$points_irregular_bins_shapefile) == 0) {
      updateSelectInput(session, "points_irregular_bin_id_column", choices = c("(Auto-detect)" = ""), selected = "")
      return()
    }
    cols <- tryCatch(available_id_columns_from_upload(input$points_irregular_bins_shapefile), error = function(e) character(0))
    choices <- c("(Auto-detect)" = "", stats::setNames(cols, cols))
    updateSelectInput(session, "points_irregular_bin_id_column", choices = choices, selected = "")
  })

  compute_irregular_richness_from_layers <- function(loaded_shapefiles,
                                                     shapefile_info,
                                                     bins_shape,
                                                     bin_id_column = NULL) {
    if (is.null(loaded_shapefiles) || length(loaded_shapefiles) == 0) {
      stop("No extrapolated layers available for irregular-polygon richness.")
    }

    bins_sf <- safe_valid_sf(normalize_to_wgs84_sf(bins_shape))
    if (nrow(bins_sf) == 0) {
      stop("No valid geometries found in irregular polygons shapefile.")
    }

    detect_bin_col <- function(sf_obj, candidate) {
      if (!is.null(candidate) && nzchar(candidate) && candidate %in% names(sf_obj)) {
        return(candidate)
      }
      cols <- setdiff(names(sf_obj), "geometry")
      if (length(cols) == 0) {
        stop("No attribute columns found in irregular polygons shapefile.")
      }
      nn <- vapply(cols, function(cc) sum(!is.na(sf_obj[[cc]])), numeric(1))
      uniq <- vapply(cols, function(cc) length(unique(as.character(sf_obj[[cc]]))), numeric(1))
      score <- nn + uniq
      cols[which.max(score)]
    }

    bin_col <- detect_bin_col(bins_sf, bin_id_column)
    bins_sf[[bin_col]] <- as.character(bins_sf[[bin_col]])

    species_layers <- list()

    for (i in seq_along(loaded_shapefiles)) {
      layer_key <- names(loaded_shapefiles)[i]
      layer_info <- if (!is.null(shapefile_info)) shapefile_info[[layer_key]] else NULL
      species_name <- layer_species_name(layer_info, layer_key)
      if (!nzchar(species_name)) {
        next
      }

      shp_sf <- tryCatch(normalize_to_wgs84_sf(loaded_shapefiles[[i]]), error = function(e) NULL)
      if (is.null(shp_sf) || nrow(shp_sf) == 0) {
        next
      }

      shp_sf <- safe_valid_sf(shp_sf)
      if (nrow(shp_sf) == 0) {
        next
      }

      shp_sf <- sf::st_transform(shp_sf, sf::st_crs(bins_sf))
      shp_sf$species_id <- species_name
      species_layers[[length(species_layers) + 1]] <- shp_sf[, c("species_id", "geometry")]
    }

    if (length(species_layers) == 0) {
      stop("No valid extrapolated geometries were available for irregular-polygon intersection.")
    }

    species_sf <- do.call(rbind, species_layers)

    # Spatial intersection assignment analogous to st_join(pot_sf, morrone_sf)
    bins_joined <- sf::st_join(
      bins_sf[, c(bin_col, "geometry")],
      species_sf[, c("species_id", "geometry")],
      join = sf::st_intersects,
      left = TRUE
    )

    species_per_bin <- bins_joined %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(.data[[bin_col]]) %>%
      dplyr::summarise(
        n_species = dplyr::n_distinct(species_id, na.rm = TRUE),
        species_list = paste(sort(unique(stats::na.omit(species_id))), collapse = ", "),
        .groups = "drop"
      )

    names(species_per_bin)[1] <- "bin_id"
    species_per_bin$species_list[species_per_bin$species_list == ""] <- NA_character_

    bins_with_richness <- bins_sf %>%
      dplyr::left_join(species_per_bin, by = setNames("bin_id", bin_col)) %>%
      dplyr::mutate(
        n_species = ifelse(is.na(n_species), 0L, as.integer(n_species)),
        species_list = ifelse(is.na(species_list), "", species_list)
      )

    bin_id_mapping <- data.frame(
      bin_name = unique(as.character(bins_sf[[bin_col]])),
      stringsAsFactors = FALSE
    )
    bin_id_mapping$bin_number <- seq_len(nrow(bin_id_mapping))

    list(
      bins_richness = bins_with_richness,
      species_per_bin = species_per_bin,
      bin_id_mapping = bin_id_mapping,
      bin_id_column = bin_col
    )
  }

  build_points_presence_absence <- function(occ_data, shape_file, grid_res) {
    if (is.null(occ_data) || nrow(occ_data) == 0) {
      stop("No occurrence points available to build point-based matrix.")
    }

    if (!all(c("spp", "long", "lat") %in% names(occ_data))) {
      stop("Occurrence data must contain columns: spp, long, lat.")
    }

    study_sf <- safe_valid_sf(normalize_to_wgs84_sf(shape_file))
    if (nrow(study_sf) == 0) {
      stop("Study area has no valid geometry for point-based matrix generation.")
    }

    occ_df <- data.frame(
      spp = as.character(occ_data$spp),
      long = as.numeric(occ_data$long),
      lat = as.numeric(occ_data$lat),
      stringsAsFactors = FALSE
    )
    occ_df <- occ_df[stats::complete.cases(occ_df), , drop = FALSE]
    if (nrow(occ_df) == 0) {
      stop("No valid occurrence points remained after filtering NAs.")
    }

    occ_sf <- sf::st_as_sf(occ_df, coords = c("long", "lat"), crs = 4326)
    occ_sf <- sf::st_transform(occ_sf, sf::st_crs(study_sf))

    within_idx <- sf::st_intersects(occ_sf, sf::st_geometry(study_sf), sparse = TRUE)
    keep_occ <- lengths(within_idx) > 0
    occ_sf <- occ_sf[keep_occ, , drop = FALSE]
    if (nrow(occ_sf) == 0) {
      stop("No occurrence points fall inside the study area shapefile.")
    }

    grid_sf <- build_regular_grid_sf(shape_file = shape_file, grid_res = grid_res)

    grid_hits <- sf::st_intersects(occ_sf, sf::st_geometry(grid_sf), sparse = TRUE)
    point_grid <- vapply(grid_hits, function(ids) if (length(ids) > 0) ids[1] else NA_integer_, integer(1))

    point_df <- data.frame(
      spp = as.character(occ_sf$spp),
      grid_id = point_grid,
      stringsAsFactors = FALSE
    )
    point_df <- point_df[!is.na(point_df$grid_id), , drop = FALSE]

    if (nrow(point_df) == 0) {
      stop("Could not assign occurrence points to grid cells.")
    }

    all_species <- sort(unique(point_df$spp))
    all_bins <- sort(unique(grid_sf$grid_id))

    pres_abs <- matrix(
      0,
      nrow = length(all_bins),
      ncol = length(all_species),
      dimnames = list(as.character(all_bins), all_species)
    )

    for (i in seq_len(nrow(point_df))) {
      pres_abs[as.character(point_df$grid_id[i]), point_df$spp[i]] <- 1
    }

    occupied_grid <- grid_sf[grid_sf$grid_id %in% unique(point_df$grid_id), c("grid_id", "geometry")]

    pres_abs <- rbind(pres_abs, ROOT = rep(0, ncol(pres_abs)))
    list(
      pres_abs = pres_abs,
      geometry = occupied_grid
    )
  }

  build_occurrence_irregular_matrix <- function(occ_data, bins_shape, bin_col, mode = c("direct", "grid_overlay"), grid_res = NULL) {
    mode <- match.arg(mode)

    bins_sf <- safe_valid_sf(normalize_to_wgs84_sf(bins_shape))
    if (!(bin_col %in% names(bins_sf))) {
      stop("Invalid irregular polygon ID column for occurrence workflow.")
    }

    occ_df <- data.frame(
      spp = as.character(occ_data$spp),
      long = as.numeric(occ_data$long),
      lat = as.numeric(occ_data$lat),
      stringsAsFactors = FALSE
    )
    occ_df <- occ_df[stats::complete.cases(occ_df), , drop = FALSE]
    if (nrow(occ_df) == 0) stop("No valid occurrence points remained after filtering NAs.")

    occ_sf <- sf::st_as_sf(occ_df, coords = c("long", "lat"), crs = 4326)
    occ_sf <- sf::st_transform(occ_sf, sf::st_crs(bins_sf))

    all_bins <- sort(unique(as.character(bins_sf[[bin_col]])))
    all_species <- sort(unique(occ_df$spp))

    if (mode == "direct") {
      occ_bins <- sf::st_join(occ_sf, bins_sf[, c(bin_col, "geometry")], join = sf::st_intersects, left = FALSE)
      if (nrow(occ_bins) == 0) stop("No occurrence points intersect the uploaded irregular polygons.")

      pairs <- unique(data.frame(
        bin_id = as.character(occ_bins[[bin_col]]),
        spp = as.character(occ_bins$spp),
        stringsAsFactors = FALSE
      ))

      bin_summary <- occ_bins %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(.data[[bin_col]]) %>%
        dplyr::summarise(
          n_species = dplyr::n_distinct(spp),
          n_occurrences = dplyr::n(),
          species_list = paste(sort(unique(spp)), collapse = ", "),
          .groups = "drop"
        )
    } else {
      if (is.null(grid_res) || !is.finite(grid_res) || grid_res <= 0) {
        stop("Grid resolution must be > 0 for irregular polygons + regular grid workflow.")
      }

      bb <- sf::st_bbox(bins_sf)
      grid_r <- raster::raster(
        xmn = bb["xmin"],
        xmx = bb["xmax"],
        ymn = bb["ymin"],
        ymx = bb["ymax"],
        resolution = c(grid_res, grid_res),
        crs = sp::CRS("+proj=longlat +datum=WGS84")
      )

      grid_sf <- sf::st_as_sf(raster::rasterToPolygons(grid_r))
      grid_sf$grid_id <- seq_len(nrow(grid_sf))

      grid_keep <- lengths(sf::st_intersects(grid_sf, sf::st_geometry(bins_sf))) > 0
      grid_sf <- grid_sf[grid_keep, , drop = FALSE]
      if (nrow(grid_sf) == 0) stop("No regular grid cells intersect the uploaded irregular polygons.")

      point_grid_hits <- sf::st_intersects(occ_sf, sf::st_geometry(grid_sf), sparse = TRUE)
      point_grid <- vapply(point_grid_hits, function(ids) {
        if (length(ids) > 0) grid_sf$grid_id[ids[1]] else NA_integer_
      }, integer(1))

      point_grid_df <- data.frame(
        spp = as.character(occ_sf$spp),
        grid_id = point_grid,
        stringsAsFactors = FALSE
      )
      point_grid_df <- point_grid_df[!is.na(point_grid_df$grid_id), , drop = FALSE]
      if (nrow(point_grid_df) == 0) stop("No occurrence points were assigned to regular grid cells inside irregular polygons.")

      grid_points <- suppressWarnings(sf::st_point_on_surface(grid_sf))
      grid_bins <- sf::st_join(
        grid_points,
        bins_sf[, c(bin_col, "geometry")],
        join = sf::st_within,
        left = TRUE
      )

      missing_bin <- is.na(grid_bins[[bin_col]])
      if (any(missing_bin)) {
        fallback <- sf::st_join(
          grid_points[missing_bin, , drop = FALSE],
          bins_sf[, c(bin_col, "geometry")],
          join = sf::st_intersects,
          left = TRUE,
          largest = TRUE
        )
        grid_bins[[bin_col]][missing_bin] <- fallback[[bin_col]]
      }

      grid_bin_map <- data.frame(
        grid_id = grid_sf$grid_id,
        bin_id = as.character(grid_bins[[bin_col]]),
        stringsAsFactors = FALSE
      )
      grid_bin_map <- grid_bin_map[!is.na(grid_bin_map$bin_id), , drop = FALSE]

      point_grid_df <- dplyr::inner_join(point_grid_df, grid_bin_map, by = "grid_id")
      if (nrow(point_grid_df) == 0) stop("No grid cells with occurrences could be linked to irregular polygons.")

      pairs <- unique(data.frame(
        bin_id = point_grid_df$bin_id,
        spp = point_grid_df$spp,
        stringsAsFactors = FALSE
      ))

      bin_summary <- point_grid_df %>%
        dplyr::group_by(bin_id) %>%
        dplyr::summarise(
          n_species = dplyr::n_distinct(spp),
          n_occurrences = dplyr::n_distinct(grid_id),
          species_list = paste(sort(unique(spp)), collapse = ", "),
          .groups = "drop"
        )
      names(bin_summary)[1] <- bin_col
    }

    pres_abs <- matrix(0, nrow = length(all_bins), ncol = length(all_species), dimnames = list(all_bins, all_species))
    if (nrow(pairs) > 0) {
      for (i in seq_len(nrow(pairs))) {
        b <- pairs$bin_id[i]
        s <- pairs$spp[i]
        if (!is.na(b) && !is.na(s) && b %in% rownames(pres_abs) && s %in% colnames(pres_abs)) {
          pres_abs[b, s] <- 1
        }
      }
    }
    pres_abs <- rbind(pres_abs, ROOT = rep(0, ncol(pres_abs)))

    bin_summary <- as.data.frame(bin_summary)
    names(bin_summary)[1] <- bin_col
    bins_with_richness <- bins_sf %>%
      dplyr::left_join(bin_summary, by = bin_col) %>%
      dplyr::mutate(
        n_species = ifelse(is.na(n_species), 0L, as.integer(n_species)),
        n_occurrences = ifelse(is.na(n_occurrences), 0L, as.integer(n_occurrences)),
        species_list = ifelse(is.na(species_list), "", species_list)
      )

    species_per_bin <- bins_with_richness %>%
      sf::st_drop_geometry() %>%
      dplyr::select(dplyr::all_of(bin_col), n_species, n_occurrences, species_list)

    bin_id_mapping <- data.frame(bin_name = all_bins, stringsAsFactors = FALSE)
    bin_id_mapping$bin_number <- seq_len(nrow(bin_id_mapping))

    list(
      pres_abs = pres_abs,
      bins_richness = bins_with_richness,
      species_per_bin = species_per_bin,
      bin_id_mapping = bin_id_mapping,
      bin_id_column = bin_col,
      geometry = bins_with_richness
    )
  }

  aggregate_regular_matrix_to_irregular_bins <- function(pres_abs, study_shape, bins_shape, bin_col, grid_res) {
    if (is.null(pres_abs) || !is.matrix(pres_abs)) {
      stop("Regular-grid matrix is missing or invalid for irregular aggregation.")
    }
    if (is.null(grid_res) || !is.finite(grid_res) || grid_res <= 0) {
      stop("Grid resolution must be > 0 to aggregate matrix by irregular polygons.")
    }

    bins_sf <- safe_valid_sf(normalize_to_wgs84_sf(bins_shape))
    if (!(bin_col %in% names(bins_sf))) {
      stop("Invalid irregular polygon ID column for matrix aggregation.")
    }

    study_sf <- safe_valid_sf(normalize_to_wgs84_sf(study_shape))
    if (nrow(study_sf) == 0) {
      stop("Study-area shapefile has no valid geometry for irregular aggregation.")
    }

    mat <- pres_abs
    if ("ROOT" %in% rownames(mat)) {
      mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
    }
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      stop("No matrix content available after removing ROOT for irregular aggregation.")
    }

    cell_ids <- suppressWarnings(as.integer(rownames(mat)))
    if (any(is.na(cell_ids))) {
      stop("Current matrix rows are not regular-grid cell IDs, so irregular aggregation is not possible for this run.")
    }

    study_sp <- as(study_sf, "Spatial")
    mask_raster <- raster::raster(
      raster::extent(study_sp),
      resolution = c(grid_res, grid_res),
      crs = sp::CRS("+proj=longlat +datum=WGS84")
    )
    area_raster <- raster::rasterize(study_sp, mask_raster)
    area_raster <- raster::merge(area_raster, mask_raster)

    max_cell <- raster::ncell(area_raster)
    in_bounds <- cell_ids >= 1 & cell_ids <= max_cell
    if (!all(in_bounds)) {
      mat <- mat[in_bounds, , drop = FALSE]
      cell_ids <- cell_ids[in_bounds]
    }
    if (nrow(mat) == 0) {
      stop("No matrix cells could be aligned with the generated regular grid.")
    }

    grid_sp <- raster::rasterToPolygons(area_raster)
    grid_sf <- sf::st_as_sf(grid_sp)
    grid_centers <- sf::st_coordinates(sf::st_centroid(grid_sf$geometry))
    grid_sf$grid_id <- raster::cellFromXY(area_raster, grid_centers)
    grid_sf <- grid_sf[!is.na(grid_sf$grid_id) & grid_sf$grid_id %in% unique(cell_ids), c("grid_id", "geometry")]

    if (nrow(grid_sf) == 0) {
      stop("No regular-grid polygons could be recovered for matrix aggregation.")
    }

    grid_points <- suppressWarnings(sf::st_centroid(grid_sf))
    grid_bins <- sf::st_join(
      grid_points,
      bins_sf[, c(bin_col, "geometry")],
      join = sf::st_within,
      left = TRUE
    )

    missing_bin <- is.na(grid_bins[[bin_col]])
    if (any(missing_bin)) {
      fallback <- sf::st_join(
        grid_points[missing_bin, , drop = FALSE],
        bins_sf[, c(bin_col, "geometry")],
        join = sf::st_intersects,
        left = TRUE,
        largest = TRUE
      )
      grid_bins[[bin_col]][missing_bin] <- fallback[[bin_col]]
    }

    grid_bin_map <- unique(data.frame(
      grid_id = as.integer(grid_bins$grid_id),
      bin_id = as.character(grid_bins[[bin_col]]),
      stringsAsFactors = FALSE
    ))
    grid_bin_map <- grid_bin_map[!is.na(grid_bin_map$grid_id) & !is.na(grid_bin_map$bin_id), , drop = FALSE]

    if (nrow(grid_bin_map) == 0) {
      stop("No matrix grid cells intersect the uploaded irregular polygons.")
    }

    mat_df <- as.data.frame(mat, stringsAsFactors = FALSE)
    mat_df$grid_id <- cell_ids
    mat_df <- dplyr::left_join(mat_df, grid_bin_map, by = "grid_id")
    mat_df <- mat_df[!is.na(mat_df$bin_id), , drop = FALSE]

    if (nrow(mat_df) == 0) {
      stop("No grid cells from the matrix could be assigned to irregular polygons.")
    }

    all_bins <- sort(unique(as.character(bins_sf[[bin_col]])))
    all_species <- colnames(mat)

    pres_abs_irregular <- matrix(
      0,
      nrow = length(all_bins),
      ncol = length(all_species),
      dimnames = list(all_bins, all_species)
    )

    by_bin <- split(mat_df, mat_df$bin_id)
    for (bin_name in names(by_bin)) {
      block <- by_bin[[bin_name]][, all_species, drop = FALSE]
      pres_abs_irregular[bin_name, ] <- as.integer(colSums(as.matrix(block), na.rm = TRUE) > 0)
    }

    pres_abs_irregular <- rbind(pres_abs_irregular, ROOT = rep(0, ncol(pres_abs_irregular)))

    agg_no_root <- pres_abs_irregular[rownames(pres_abs_irregular) != "ROOT", , drop = FALSE]
    species_list <- vapply(seq_len(nrow(agg_no_root)), function(i) {
      spp <- colnames(agg_no_root)[agg_no_root[i, ] > 0]
      if (length(spp) == 0) "" else paste(spp, collapse = ", ")
    }, character(1))

    n_occ_per_bin <- vapply(all_bins, function(x) {
      length(unique(mat_df$grid_id[mat_df$bin_id == x]))
    }, integer(1))

    species_per_bin <- data.frame(
      n_species = as.integer(rowSums(agg_no_root, na.rm = TRUE)),
      n_occurrences = as.integer(n_occ_per_bin),
      species_list = species_list,
      stringsAsFactors = FALSE
    )
    species_per_bin[[bin_col]] <- rownames(agg_no_root)
    species_per_bin <- species_per_bin[, c(bin_col, "n_species", "n_occurrences", "species_list")]

    bins_with_richness <- bins_sf %>%
      dplyr::left_join(species_per_bin, by = bin_col) %>%
      dplyr::mutate(
        n_species = ifelse(is.na(n_species), 0L, as.integer(n_species)),
        n_occurrences = ifelse(is.na(n_occurrences), 0L, as.integer(n_occurrences)),
        species_list = ifelse(is.na(species_list), "", species_list)
      )

    bin_id_mapping <- data.frame(bin_name = all_bins, stringsAsFactors = FALSE)
    bin_id_mapping$bin_number <- seq_len(nrow(bin_id_mapping))

    list(
      pres_abs = pres_abs_irregular,
      bins_richness = bins_with_richness,
      species_per_bin = species_per_bin,
      bin_id_mapping = bin_id_mapping,
      bin_id_column = bin_col,
      geometry = bins_with_richness
    )
  }

  # Update taxa choices when data is loaded
  observe({
    taxa_from_occ <- if (!is.null(data_store$occurrence)) unique(data_store$occurrence$spp) else character(0)

    taxa_from_layers <- character(0)
    if (!is.null(data_store$shapefile_info) && length(data_store$shapefile_info) > 0) {
      taxa_from_layers <- unique(vapply(data_store$shapefile_info, function(info) {
        sp <- info$species_id %||% NA_character_
        ifelse(is.na(sp), "", as.character(sp))
      }, character(1)))
      taxa_from_layers <- taxa_from_layers[nzchar(taxa_from_layers)]
    }

    taxa_list <- sort(unique(c(taxa_from_occ, taxa_from_layers)))

    updateSelectizeInput(session, "mst_taxa_select", choices = taxa_from_occ)
  })

  observe({
    if (!is.null(data_store$tree)) {
      n_tips <- ape::Ntip(data_store$tree)
      n_nodes <- data_store$tree$Nnode
      if (!is.null(n_nodes) && n_nodes > 0) {
        internal_nodes <- as.character((n_tips + 1):(n_tips + n_nodes))
        updateSelectizeInput(session, "mst_node_ids", choices = internal_nodes, selected = character(0))
        updateSelectizeInput(session, "clade_node_id", choices = internal_nodes, selected = character(0))
      } else {
        updateSelectizeInput(session, "mst_node_ids", choices = character(0), selected = character(0))
        updateSelectizeInput(session, "clade_node_id", choices = character(0), selected = character(0))
      }
    } else {
      updateSelectizeInput(session, "mst_node_ids", choices = character(0), selected = character(0))
      updateSelectizeInput(session, "clade_node_id", choices = character(0), selected = character(0))
    }
  })

  output$clade_filter_status <- renderPrint({
    if (!isTRUE(input$use_clade_filter)) {
      cat("Clade filter is disabled.\n")
      return()
    }
    if (is.null(data_store$tree)) {
      cat("Load a phylogenetic tree to enable clade filtering.\n")
      return()
    }
    if (is.null(input$clade_node_id) || !nzchar(input$clade_node_id)) {
      cat("Select an internal node to filter taxa for extrapolation.\n")
      return()
    }

    node_num <- suppressWarnings(as.integer(input$clade_node_id))
    if (is.na(node_num)) {
      cat("Invalid node ID.\n")
      return()
    }

    n_tips <- ape::Ntip(data_store$tree)
    n_nodes <- data_store$tree$Nnode
    valid_internal <- (n_tips + 1):(n_tips + n_nodes)
    if (!(node_num %in% valid_internal)) {
      cat("Selected node is not a valid internal node.\n")
      return()
    }

    clade_tree <- ape::extract.clade(data_store$tree, node = node_num)
    clade_tips <- clade_tree$tip.label

    cat("Clade filter enabled\n")
    cat("Node:", node_num, "\n")
    cat("Descendant taxa:", length(clade_tips), "\n")
  })
  
  # Render tree with node numbers for MST node selection
  output$mst_tree_plot <- renderPlot({
    if (is.null(data_store$tree)) {
      plot(1, type = "n", main = "No tree loaded", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "Please load a phylogenetic tree in Step 1", cex = 1.2)
    } else {
      tree <- data_store$tree
      
      # Plot tree with nice styling
      par(mar = c(2.5, 1.5, 5, 1))
      phytools::plotTree(tree, fsize = 1.05, lwd = 1.8, color = "#1f1f1f", 
                         offset = 0.7, main = "Phylogenetic Tree with Node Numbers")
      
      # Add node labels
      n_tips <- ape::Ntip(tree)
      n_nodes <- tree$Nnode
      all_nodes <- 1:(n_tips + n_nodes)
      
      # Only label internal nodes to avoid cluttering
      internal_nodes <- (n_tips + 1):(n_tips + n_nodes)
      
      phytools::labelnodes(text = internal_nodes, node = internal_nodes,
                           interactive = FALSE, cex = 1.15,
                           circle.exp = 1.7, col = "white", bg = "#a61c1c")
      
      legend("bottomleft", legend = "Internal Nodes", pch = 21, 
             pt.bg = "#d9534f", col = "white", pt.cex = 2.2, bty = "n")
    }
  })

  # Handle shapefile upload
  observeEvent(input$study_area_shapefile, {
    tryCatch({
      if (is.null(input$study_area_shapefile) || nrow(input$study_area_shapefile) == 0) {
        output$shapefile_status <- renderPrint({
          cat("No shapefile loaded\n")
        })
      } else {
        # Load shapefile from multiple files
        shape <- load_shapefile_from_files(input$study_area_shapefile)
        data_store$study_area_shapefile <- shape
        
        output$shapefile_status <- renderPrint({
          cat("✓ Shapefile loaded successfully!\n")
          cat("Files uploaded:", nrow(input$study_area_shapefile), "\n")
          for (i in seq_len(nrow(input$study_area_shapefile))) {
            cat("  ✓", input$study_area_shapefile$name[i], "\n")
          }
        })
      }
    }, error = function(e) {
      output$shapefile_status <- renderPrint({
        cat("✗ Error loading shapefile:\n", e$message, "\n\n")
        cat("Troubleshooting:\n")
        cat("1. Make sure you selected ALL files together (use Ctrl+Click)\n")
        cat("2. Required files: .shp, .shx, .dbf\n")
        cat("3. All files must have the SAME base name (e.g., America_Sul.*)\n")
        cat("4. Try uploading again with all files selected at once\n")
      })
    })
  })

  observeEvent(input$plot_points_only, {
    if (is.null(data_store$occurrence) || nrow(data_store$occurrence) == 0) {
      reset_extrap_log("Points-only map not started: occurrence data is missing.")
      output$extrap_status <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
      return(invisible(NULL))
    }

    tryCatch({
      run_started_at <- Sys.time()
      reset_extrap_log(c(
        "=== Points-Only Map Log ===",
        paste0("Started: ", format(run_started_at, "%Y-%m-%d %H:%M:%S"))
      ))

      occ_data <- data_store$occurrence
      tree <- data_store$tree

      if (isTRUE(input$use_clade_filter)) {
        clade_subset <- extract_clade_subset(
          occ_data = occ_data,
          tree = tree,
          node_id = input$clade_node_id
        )
        occ_data <- clade_subset$occurrence
        tree <- clade_subset$tree
        append_extrap_log(paste0("Clade filter active (node ", clade_subset$node, ")."))
        append_extrap_log(paste0("Rows after clade filter: ", nrow(occ_data)))
      }

      data_store$analysis_occurrence <- occ_data
      data_store$analysis_tree <- tree
      data_store$loaded_shapefiles <- list()
      data_store$loaded_rasters <- list()
      data_store$shapefile_info <- list()
      data_store$regular_grid_all_sf <- NULL
      data_store$regular_grid_presence_sf <- NULL
      data_store$regular_grid_presence_points_sf <- NULL
      data_store$regular_grid_presence_extrap_sf <- NULL
      data_store$regular_grid_presence_points_by_taxon <- list()
      data_store$regular_grid_presence_extrap_by_taxon <- list()
      data_store$regular_grid_presence_extrap_label <- NULL
      data_store$regular_grid_presence_by_taxon <- list()
      data_store$regular_grid_res <- NULL
      data_store$secondary_study_areas <- list()
      data_store$visual_layers_method <- "occurrence_only"
      data_store$last_extrap_started_at <- run_started_at
      data_store$last_extrap_taxa <- unique(as.character(occ_data$spp))
      distribution_visible_groups(NULL)

      ensure_species_palette(extra_species = unique(as.character(occ_data$spp)))

      append_extrap_log("Loaded occurrence points directly into visualization map.")
      append_extrap_log("No extrapolation was executed.")
      append_extrap_log(paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      output$extrap_status <- renderPrint({
        cat("✓ Points loaded for visualization (no extrapolation)\n")
        cat("Rows:", nrow(occ_data), "\n")
        cat("Species:", length(unique(occ_data$spp)), "\n")
        cat("Map mode: occurrence points only\n")
      })
    }, error = function(e) {
      append_extrap_log(paste0("[error] ", e$message))
      output$extrap_status <- renderPrint({
        cat("Error loading points-only map:\n", e$message, "\n")
      })
    })
  })
  
  observeEvent(input$run_extrap, {
    if (is.null(data_store$occurrence)) {
      reset_extrap_log("No run started: occurrence data is missing.")
      output$extrap_status <- renderPrint({
        cat("Error: Please load occurrence data first\n")
      })
    } else {
      tryCatch({
        run_started_at <- Sys.time()
        reset_extrap_log(c(
          "=== Range Extrapolation Log ===",
          paste0("Started: ", format(run_started_at, "%Y-%m-%d %H:%M:%S"))
        ))

        # Use already preprocessed data from the preprocessing tab
        occ_data <- data_store$occurrence
        tree <- data_store$tree
        append_extrap_log(paste0("Loaded occurrence rows: ", nrow(occ_data)))
        append_extrap_log(paste0("Loaded occurrence taxa: ", length(unique(occ_data$spp))))
        append_extrap_log(paste0("Tree loaded: ", if (is.null(tree)) "no" else paste0("yes (", length(tree$tip.label), " tips)")))
        
        # Data should already be cleaned from preprocessing tab
        if (is.null(occ_data) || nrow(occ_data) == 0) {
          append_extrap_log("No occurrence data available after preprocessing.")
          output$extrap_status <- renderPrint({
            cat("Error: No occurrence data available. Please run Data Preprocessing first.
")
          })
          return()
        }

        if (isTRUE(input$use_clade_filter)) {
          clade_subset <- extract_clade_subset(
            occ_data = occ_data,
            tree = tree,
            node_id = input$clade_node_id
          )
          occ_data <- clade_subset$occurrence
          tree <- clade_subset$tree
          append_extrap_log(paste0("Clade filter active (node ", clade_subset$node, ")."))
          append_extrap_log(paste0("Clade descendant taxa: ", length(clade_subset$tips)))
          append_extrap_log(paste0("Rows after clade filter: ", nrow(occ_data)))
        }

        method <- input$extrap_method
        removed_before_run <- clear_visual_output_files(method = method)
        if (length(removed_before_run) > 0) {
          append_extrap_log(paste0("Cleared ", length(removed_before_run), " previous visual output file(s) for method '", method, "'."))
        }

        if (!all(c("spp", "long", "lat") %in% names(occ_data))) {
          stop("Occurrence data must include columns spp, long and lat.")
        }

        if (identical(method, "mst")) {
          occ_data_ordered <- data.frame(
            spp = occ_data$spp,
            lat = occ_data$lat,
            long = occ_data$long,
            stringsAsFactors = FALSE
          )

          singleton_result <- singleton_to_data_frame(data = occ_data_ordered, phylogeny = tree)
          occ_data <- singleton_result$data_df
          tree <- singleton_result$treeMod %||% singleton_result$treeNonMod %||% tree
          data_store$singleton_result <- singleton_result

          append_extrap_log("Applied MST singleton/tree harmonization.")
          append_extrap_log(paste0("Rows after MST harmonization: ", nrow(occ_data)))
          append_extrap_log(paste0("Taxa after MST harmonization: ", length(unique(occ_data$spp))))
          append_extrap_log(paste0("Tree after MST harmonization: ", if (is.null(tree)) "none" else paste0(length(tree$tip.label), " tips")))
        } else if (identical(method, "convex_hull")) {
          spp_counts <- table(occ_data$spp)
          lt3_taxa <- names(spp_counts[spp_counts < 3])
          if (length(lt3_taxa) > 0) {
            occ_data <- occ_data[!(occ_data$spp %in% lt3_taxa), , drop = FALSE]
            append_extrap_log(paste0("Convex hull fallback removed ", length(lt3_taxa), " taxa with < 3 points."))
            if (nrow(occ_data) == 0) {
              stop("No taxa with >= 3 points are available for convex hull extrapolation.")
            }
          }
          data_store$singleton_result <- NULL
        } else {
          data_store$singleton_result <- NULL
        }
        points_mode <- if (identical(method, "occurrence_only")) (input$points_occurrence_mode %||% "regular_grid") else "regular_grid"
        points_use_irregular <- identical(method, "occurrence_only") && points_mode %in% c("irregular_direct", "irregular_with_grid")
        points_use_grid_overlay <- identical(method, "occurrence_only") && identical(points_mode, "irregular_with_grid")
        grid_res <- input$grid_resolution
        data_store$irregular_bins_richness <- NULL
        data_store$irregular_bins_species_table <- NULL
        data_store$irregular_bins_id_column <- NULL
        data_store$irregular_bins_bin_id_mapping <- NULL
        data_store$irregular_bins_method <- NULL
        data_store$matrix_is_irregular_aggregated <- FALSE
        data_store$secondary_study_areas <- list()
        data_store$regular_grid_all_sf <- NULL
        data_store$regular_grid_presence_sf <- NULL
        data_store$regular_grid_presence_points_sf <- NULL
        data_store$regular_grid_presence_extrap_sf <- NULL
        data_store$regular_grid_presence_points_by_taxon <- list()
        data_store$regular_grid_presence_extrap_by_taxon <- list()
        data_store$regular_grid_presence_by_taxon <- list()
        data_store$regular_grid_res <- NULL
        
        # Check if shapefile is loaded
        if (is.null(data_store$study_area_shapefile)) {
          append_extrap_log("Study-area shapefile not loaded.")
          output$extrap_status <- renderPrint({
            cat("Error: Please upload a shapefile first!\n")
          })
          return()
        }
        
        shape_file <- data_store$study_area_shapefile
        
        # Ensure occurrence data has the correct column names for the functions
        # calcRange_convexHull expects 'species', 'long', 'lat'
        if ("spp" %in% names(occ_data) && !("species" %in% names(occ_data))) {
          occ_data$species <- occ_data$spp
        }

        if (method == "mst" && !is.null(tree)) {
          taxa_data <- unique(as.character(occ_data$spp))
          taxa_tree_only <- setdiff(tree$tip.label, taxa_data)
          taxa_data_only <- setdiff(taxa_data, tree$tip.label)

          if (length(taxa_tree_only) > 0) {
            tree <- ape::drop.tip(tree, taxa_tree_only)
            append_extrap_log(paste0("Pruned ", length(taxa_tree_only), " tree-only taxa before MST."))
          }
          if (length(taxa_data_only) > 0) {
            occ_data <- occ_data[!(occ_data$spp %in% taxa_data_only), , drop = FALSE]
            append_extrap_log(paste0("Removed ", length(taxa_data_only), " data-only taxa before MST."))
          }
        }

        append_extrap_log(paste0("Method selected: ", method))
        mst_level_selected <- if (identical(method, "mst")) (input$mst_level %||% "single_all") else NULL
        
        # Call the appropriate extrapolation function and capture console logs
        method_run <- capture_analysis_output({
          if (method == "buffer") {
            buffer_width <- input$buffer_width * 1000
            calcRange_buffers(
              xy = occ_data,
              buffer.width = buffer_width,
              shape_file = shape_file,
              resol = grid_res,
              mean_dist = input$use_mean_dist
            )
          } else if (method == "convex_hull") {
            calcRange_convexHull(
              xy = occ_data,
              shape_file = shape_file,
              resol = grid_res
            )
          } else if (method == "mst") {
            # Get MST specific parameters
            mst_level <- mst_level_selected

            # Prepare parameters for MST_node
            nodes_param <- NULL
            taxon_param <- NULL
            mintreeall_param <- FALSE

            if (mst_level == "single_all") {
              mintreeall_param <- TRUE
            } else if (mst_level == "node") {
              if (is.null(tree)) {
                stop("A phylogenetic tree is required for node-based MST analysis. Please load a tree in Step 1.")
              }
              if (is.null(input$mst_node_ids) || length(input$mst_node_ids) == 0) {
                stop("Please select at least one internal node for node-based MST analysis.")
              }
              nodes_param <- as.integer(input$mst_node_ids)
              mintreeall_param <- TRUE
            } else if (mst_level == "node_taxa") {
              if (is.null(tree)) {
                stop("A phylogenetic tree is required for node + taxa MST analysis. Please load a tree in Step 1.")
              }
              if (is.null(input$mst_node_ids) || length(input$mst_node_ids) == 0) {
                stop("Please select at least one internal node for node + taxa MST analysis.")
              }
              nodes_param <- as.integer(input$mst_node_ids)
              if (!is.null(input$mst_taxa_select) && length(input$mst_taxa_select) > 0) {
                taxon_param <- input$mst_taxa_select
              }
              mintreeall_param <- TRUE
            } else if (mst_level == "taxa") {
              if (is.null(input$mst_taxa_select) || length(input$mst_taxa_select) == 0) {
                stop("Please select at least one taxon for MST analysis.")
              }
              taxon_param <- input$mst_taxa_select
              mintreeall_param <- TRUE
            }

            # Use singleton_result data for MST
            mst_data <- occ_data
            mst_tree <- tree

            if ((mst_level %in% c("node", "node_taxa")) && !is.null(mst_tree) && !is.null(nodes_param) && length(nodes_param) > 0) {
              descendant_tips <- unique(unlist(lapply(nodes_param, function(nn) {
                tryCatch(ape::extract.clade(mst_tree, node = as.integer(nn))$tip.label, error = function(e) character(0))
              })))
              descendant_tips <- descendant_tips[nzchar(descendant_tips)]
              if (length(descendant_tips) > 0) {
                taxon_param <- unique(c(as.character(taxon_param), descendant_tips))
                append_extrap_log(paste0("MST node mode uses descendant tips as taxa set: ", length(taxon_param), " taxa."))
              }
            }

            mst_shape_for_grid <- shape_file
            if (is.finite(grid_res) && grid_res > 0) {
              mst_shape_for_grid <- tryCatch(
                build_bbox_grid_shape(shape_file = shape_file, pad_degrees = 0),
                error = function(e) shape_file
              )
            }
            append_extrap_log("MST grid extent uses study-area bbox for full border cells.")

            # Call MST_node function
            MST_node(
              coordin = mst_data,
              shape_file = mst_shape_for_grid,
              resol = c(grid_res, grid_res),
              tree = mst_tree,
              nodes = nodes_param,
              taxon = taxon_param,
              mintreeall = mintreeall_param,
              caption = TRUE,
              sobrepo = FALSE,
              seeres = FALSE
            )
          } else if (method == "occurrence_only") {
            if (points_use_irregular) {
              if (is.null(input$points_irregular_bins_shapefile) || nrow(input$points_irregular_bins_shapefile) == 0) {
                stop("Please upload the irregular polygons shapefile to build the point-based matrix.")
              }

              bins_shape <- load_shapefile_from_files(input$points_irregular_bins_shapefile)
              bins_sf <- sf::st_as_sf(bins_shape)
              bin_col <- input$points_irregular_bin_id_column

              if (!is.null(bin_col) && nzchar(bin_col) && !(bin_col %in% names(bins_sf))) {
                stop(paste0("Invalid irregular polygon ID column: '", bin_col,
                            "'. Available columns: ", paste(setdiff(names(bins_sf), "geometry"), collapse = ", ")))
              }

              if (is.null(bin_col) || !nzchar(bin_col)) {
                candidate_cols <- setdiff(names(bins_sf), "geometry")
                valid_cols <- candidate_cols[sapply(candidate_cols, function(x) sum(!is.na(bins_sf[[x]])) > 0)]
                if (length(valid_cols) == 0) {
                  stop("Could not detect a valid ID column in the irregular polygons shapefile.")
                }
                nn <- sapply(valid_cols, function(x) sum(!is.na(bins_sf[[x]])))
                bin_col <- valid_cols[which.max(nn)]
              }

              out <- if (points_use_grid_overlay) {
                build_occurrence_irregular_matrix(
                  occ_data = occ_data,
                  bins_shape = bins_shape,
                  bin_col = bin_col,
                  mode = "grid_overlay",
                  grid_res = grid_res
                )
              } else {
                build_occurrence_irregular_matrix(
                  occ_data = occ_data,
                  bins_shape = bins_shape,
                  bin_col = bin_col,
                  mode = "direct"
                )
              }

              secondary_label <- if (points_use_grid_overlay) {
                "Irregular Polygons (Points + Regular Grid Overlay)"
              } else {
                "Irregular Polygons (Points Direct)"
              }
              data_store$secondary_study_areas <- c(data_store$secondary_study_areas, setNames(list(bins_shape), secondary_label))
              out
            } else {
              build_points_presence_absence(
                occ_data = occ_data,
                shape_file = shape_file,
                grid_res = grid_res
              )
            }
          } else if (method == "irregular_bins") {
            if (is.null(input$irregular_bins_shapefile) || nrow(input$irregular_bins_shapefile) == 0) {
              stop("Please upload a shapefile for irregular bins")
            }

            bins_shape <- load_shapefile_from_files(input$irregular_bins_shapefile)
            bins_sf <- sf::st_as_sf(bins_shape)
            bin_col <- input$irregular_bin_id_column

            if (!is.null(bin_col) && nzchar(bin_col) && !(bin_col %in% names(bins_sf))) {
              stop(paste0("Invalid irregular bins ID column: '", bin_col,
                          "'. Available columns: ", paste(setdiff(names(bins_sf), "geometry"), collapse = ", ")))
            }

            if (is.null(bin_col) || !nzchar(bin_col)) {
              candidate_cols <- setdiff(names(bins_sf), "geometry")
              valid_cols <- candidate_cols[sapply(candidate_cols, function(x) sum(!is.na(bins_sf[[x]])) > 0)]
              if (length(valid_cols) == 0) {
                stop("Could not detect a valid bin ID column in irregular bins shapefile.")
              }
              # choose the column with the largest number of non-NA entries
              nn <- sapply(valid_cols, function(x) sum(!is.na(bins_sf[[x]])))
              bin_col <- valid_cols[which.max(nn)]
            }

            calcRange_irregular_bins(
              xy = occ_data,
              bins_shapefile = bins_shape,
              bin_id_column = bin_col,
              resol = c(grid_res, grid_res)
            )
          }
        })

        result <- method_run$result
        if (length(method_run$log) > 0) {
          append_extrap_log(method_run$log)
        }

        if (identical(method, "mst") && (mst_level_selected %in% c("node", "node_taxa"))) {
          rebuilt <- tryCatch(
            ensure_mst_node_intersection_outputs(shape_file = shape_file, grid_res = grid_res),
            error = function(e) {
              append_extrap_log(paste0("[warning] Could not rebuild MST node intersection outputs: ", e$message))
              FALSE
            }
          )
          if (isTRUE(rebuilt)) {
            append_extrap_log("Rebuilt MST node outputs from MST-grid intersections (line + grid + presence raster).")
          }

          mst_outputs_now <- list.files(file.path(getwd(), "out_MST"), pattern = "(mst_.*\\.shp$|GRIDS_.*\\.shp$|presence_.*\\.tif$)", full.names = FALSE)
          append_extrap_log(paste0("MST output files available: ", length(mst_outputs_now)))
          if (length(mst_outputs_now) > 0) {
            append_extrap_log(paste0("MST files: ", paste(sort(mst_outputs_now), collapse = ", ")))
          }
        }

        if (is.matrix(result)) {
          data_store$pres_abs <- result
          data_store$geometry <- NULL
        } else {
          data_store$pres_abs <- result$pres_abs
          data_store$geometry <- result$geometry %||% result$bins_richness

          if (!is.null(result$bins_richness) && !is.null(result$species_per_bin)) {
            data_store$irregular_bins_richness <- result$bins_richness
            data_store$irregular_bins_species_table <- result$species_per_bin
            data_store$irregular_bins_bin_id_mapping <- result$bin_id_mapping %||% NULL

            result_id_col <- names(result$species_per_bin)[1]
            if (is.null(result_id_col) || !nzchar(result_id_col)) {
              result_id_col <- NA_character_
            }
            data_store$irregular_bins_id_column <- result_id_col
            data_store$irregular_bins_method <- if (method == "occurrence_only" && points_use_irregular) {
              if (points_use_grid_overlay) "occurrence_only_irregular_grid_overlay" else "occurrence_only_irregular_direct"
            } else {
              method
            }
          }
        }

        data_store$pres_abs_regular <- NULL
        data_store$pres_abs_irregular <- NULL
        if (identical(method, "occurrence_only") && points_use_irregular) {
          data_store$pres_abs_irregular <- data_store$pres_abs
        } else if (identical(method, "irregular_bins")) {
          data_store$pres_abs_irregular <- data_store$pres_abs
        } else {
          data_store$pres_abs_regular <- data_store$pres_abs
        }

        if (identical(method, "occurrence_only") && isTRUE(points_use_grid_overlay) &&
            (is.null(data_store$pres_abs_regular) || !is.matrix(data_store$pres_abs_regular))) {
          regular_overlay <- tryCatch(
            build_points_presence_absence(
              occ_data = occ_data,
              shape_file = shape_file,
              grid_res = grid_res
            ),
            error = function(e) NULL
          )
          if (!is.null(regular_overlay) && is.list(regular_overlay) && !is.null(regular_overlay$pres_abs)) {
            data_store$pres_abs_regular <- regular_overlay$pres_abs
            append_extrap_log("Built auxiliary regular-grid matrix for occurrence points + irregular overlay workflow.")
          }
        }

        regular_grid_needed <- !identical(method, "irregular_bins") && !(
          identical(method, "occurrence_only") && identical(points_mode, "irregular_direct")
        )
        if (isTRUE(regular_grid_needed) && is.finite(grid_res) && grid_res > 0) {
          grid_sf <- tryCatch({
            if (identical(method, "convex_hull") && is.list(result) && !is.null(result$grid_all_sf)) {
              normalize_to_wgs84_sf(result$grid_all_sf)
            } else {
              build_regular_grid_sf(shape_file = shape_file, grid_res = grid_res)
            }
          }, error = function(e) NULL)
          extrap_methods_with_dual_grid <- c("convex_hull", "buffer", "mst")
          extrap_method_label <- switch(method,
            convex_hull = "Convex Hull",
            buffer = "Buffer",
            mst = "MST",
            "Extrapolation"
          )
          if (!is.null(grid_sf) && nrow(grid_sf) > 0 && !is.null(data_store$pres_abs_regular)) {
            geom_col <- attr(grid_sf, "sf_column")
            if (!is.null(geom_col) && nzchar(geom_col) && geom_col %in% names(grid_sf) && !identical(geom_col, "geometry")) {
              names(grid_sf)[names(grid_sf) == geom_col] <- "geometry"
              sf::st_geometry(grid_sf) <- "geometry"
            }
            if (!("grid_id" %in% names(grid_sf)) && nrow(grid_sf) > 0) {
              grid_sf$grid_id <- seq_len(nrow(grid_sf))
            }
            grid_sf <- grid_sf[, c("grid_id", "geometry")]
            data_store$pres_abs_regular <- align_regular_matrix_to_grid(data_store$pres_abs_regular, grid_sf)
            mat_no_root <- data_store$pres_abs_regular
            if ("ROOT" %in% rownames(mat_no_root)) {
              mat_no_root <- mat_no_root[rownames(mat_no_root) != "ROOT", , drop = FALSE]
            }
            occ_ids <- suppressWarnings(as.integer(rownames(mat_no_root)))
            occ_ids <- unique(occ_ids[!is.na(occ_ids) & rowSums(mat_no_root, na.rm = TRUE) > 0])
            data_store$regular_grid_all_sf <- grid_sf
            matrix_presence_sf <- if (length(occ_ids) > 0) {
              grid_sf[grid_sf$grid_id %in% occ_ids, c("grid_id", "geometry")]
            } else {
              grid_sf[0, c("grid_id", "geometry")]
            }
            data_store$regular_grid_presence_sf <- matrix_presence_sf
            data_store$regular_grid_presence_points_sf <- NULL
            data_store$regular_grid_presence_extrap_sf <- NULL
            data_store$regular_grid_presence_points_by_taxon <- list()
            data_store$regular_grid_presence_extrap_by_taxon <- list()
            presence_by_taxon <- list()
            taxa_cols <- colnames(mat_no_root)
            for (tx in taxa_cols) {
              tx_hits <- suppressWarnings(as.integer(rownames(mat_no_root)[mat_no_root[, tx] > 0]))
              tx_hits <- unique(tx_hits[!is.na(tx_hits)])
              if (length(tx_hits) == 0) next
              presence_by_taxon[[tx]] <- grid_sf[grid_sf$grid_id %in% tx_hits, c("grid_id", "geometry")]
            }
            if (method %in% extrap_methods_with_dual_grid) {
              data_store$regular_grid_presence_extrap_sf <- matrix_presence_sf
              data_store$regular_grid_presence_extrap_by_taxon <- presence_by_taxon
              data_store$regular_grid_presence_extrap_label <- extrap_method_label
              points_overlay <- tryCatch(
                build_points_presence_absence(
                  occ_data = occ_data,
                  shape_file = shape_file,
                  grid_res = grid_res
                ),
                error = function(e) NULL
              )
              if (!is.null(points_overlay) && is.list(points_overlay) && !is.null(points_overlay$pres_abs)) {
                points_mat <- tryCatch(align_regular_matrix_to_grid(points_overlay$pres_abs, grid_sf), error = function(e) NULL)
                if (!is.null(points_mat)) {
                  points_no_root <- points_mat
                  if ("ROOT" %in% rownames(points_no_root)) {
                    points_no_root <- points_no_root[rownames(points_no_root) != "ROOT", , drop = FALSE]
                  }
                  point_ids <- suppressWarnings(as.integer(rownames(points_no_root)))
                  point_ids <- unique(point_ids[!is.na(point_ids) & rowSums(points_no_root, na.rm = TRUE) > 0])
                  data_store$regular_grid_presence_points_sf <- if (length(point_ids) > 0) {
                    grid_sf[grid_sf$grid_id %in% point_ids, c("grid_id", "geometry")]
                  } else {
                    grid_sf[0, c("grid_id", "geometry")]
                  }
                  points_by_taxon <- list()
                  point_taxa_cols <- colnames(points_no_root)
                  for (tx in point_taxa_cols) {
                    tx_hits <- suppressWarnings(as.integer(rownames(points_no_root)[points_no_root[, tx] > 0]))
                    tx_hits <- unique(tx_hits[!is.na(tx_hits)])
                    if (length(tx_hits) == 0) next
                    points_by_taxon[[tx]] <- grid_sf[grid_sf$grid_id %in% tx_hits, c("grid_id", "geometry")]
                  }
                  data_store$regular_grid_presence_points_by_taxon <- points_by_taxon
                }
              }
            } else {
              data_store$regular_grid_presence_extrap_label <- NULL
            }
            data_store$regular_grid_presence_by_taxon <- presence_by_taxon
            data_store$regular_grid_res <- grid_res
            if (!identical(method, "occurrence_only") || !points_use_irregular) {
              data_store$pres_abs <- data_store$pres_abs_regular
            }
            write_regular_grid_layers(
              method = method,
              grid_sf = grid_sf,
              regular_matrix = data_store$pres_abs_regular,
              grid_res = grid_res
            )
          }
        }

        data_store$mst_pres_abs <- if (method == "mst") data_store$pres_abs else NULL
        data_store$mst_context <- if (method == "mst") list(
          preabsMat = data_store$pres_abs,
          shapeFile = shape_file,
          resol = c(grid_res, grid_res),
          taxa = colnames(data_store$pres_abs),
          timestamp = Sys.time()
        ) else NULL
        data_store$extrap_method <- method
        data_store$analysis_occurrence <- occ_data
        data_store$analysis_tree <- tree
        data_store$extrap_shapefile <- shape_file
        data_store$study_area <- shape_file
        
        if (isTRUE(input$enable_irregular_richness) && method %in% c("buffer", "convex_hull", "mst")) {
          if (is.null(input$irregular_richness_shapefile) || nrow(input$irregular_richness_shapefile) == 0) {
            stop("To compute diversity by irregular polygons, upload the SECOND subdivision shapefile (.shp, .shx, .dbf).")
          }
        }

        # Auto-load output files for visualization
        load_method <- if (method == "occurrence_only" && points_use_irregular) {
          "irregular_bins"
        } else {
          method
        }
        taxa_filter_for_layers <- if (identical(method, "mst")) character(0) else unique(as.character(occ_data$spp))
        output_files <- load_output_shapefiles(
          extrap_method = load_method,
          since_time = NULL,
          taxa_filter = taxa_filter_for_layers
        )
        output_files <- filter_visual_outputs(output_files)
        data_store$loaded_shapefiles <- output_files$shapefiles
        data_store$loaded_rasters <- output_files$rasters
        data_store$shapefile_info <- output_files$shapefile_info
        data_store$visual_layers_method <- if (method == "occurrence_only") "occurrence_only" else method
        data_store$last_extrap_started_at <- run_started_at
        data_store$last_extrap_taxa <- if (identical(method, "mst")) character(0) else unique(as.character(occ_data$spp))
        distribution_visible_groups(NULL)

        loaded_keys <- names(data_store$loaded_shapefiles)
        grid_loaded <- loaded_keys[grepl("GRIDS_", loaded_keys, ignore.case = TRUE)]
        append_extrap_log(paste0("Loaded shapefile layers: ", length(loaded_keys)))
        if (length(grid_loaded) > 0) {
          append_extrap_log(paste0("Loaded grid layers: ", paste(grid_loaded, collapse = ", ")))
        } else {
          append_extrap_log("Loaded grid layers: none found in output loader")
        }

        ensure_species_palette(extra_species = unique(occ_data$spp))

        if (isTRUE(input$enable_irregular_richness) && method %in% c("buffer", "convex_hull", "mst")) {
          bins_shape <- load_shapefile_from_files(input$irregular_richness_shapefile)
          richness_result <- compute_irregular_richness_from_layers(
            loaded_shapefiles = data_store$loaded_shapefiles,
            shapefile_info = data_store$shapefile_info,
            bins_shape = bins_shape,
            bin_id_column = input$irregular_richness_id_column
          )

          data_store$irregular_bins_richness <- richness_result$bins_richness
          data_store$irregular_bins_species_table <- richness_result$species_per_bin
          data_store$irregular_bins_id_column <- richness_result$bin_id_column
          data_store$irregular_bins_method <- method

          matrix_result <- aggregate_regular_matrix_to_irregular_bins(
            pres_abs = data_store$pres_abs_regular %||% data_store$pres_abs,
            study_shape = shape_file,
            bins_shape = bins_shape,
            bin_col = richness_result$bin_id_column,
            grid_res = grid_res
          )

          data_store$pres_abs <- matrix_result$pres_abs
          data_store$pres_abs_irregular <- matrix_result$pres_abs
          data_store$geometry <- matrix_result$geometry
          data_store$irregular_bins_richness <- matrix_result$bins_richness
          data_store$irregular_bins_species_table <- matrix_result$species_per_bin
          data_store$irregular_bins_bin_id_mapping <- matrix_result$bin_id_mapping
          data_store$irregular_bins_id_column <- matrix_result$bin_id_column
          data_store$irregular_bins_method <- paste0(method, "_matrix_aggregated_to_irregular")
          data_store$matrix_is_irregular_aggregated <- TRUE

          append_extrap_log(
            paste0("Presence-absence matrix aggregated from regular grid to irregular polygons (", data_store$irregular_bins_id_column, ").")
          )

          data_store$secondary_study_areas <- c(
            data_store$secondary_study_areas,
            list("Irregular Polygons (Diversity)" = bins_shape)
          )
        }
        
        output$extrap_status <- renderPrint({
          cat("✓ Extrapolation completed!\n")
          cat("Method:", method, "\n")
          if (!identical(method, "occurrence_only") || !identical(points_mode, "irregular_direct")) {
            cat("Grid resolution:", grid_res, "degrees\n")
          }
          cat("Presence-absence matrix:", nrow(data_store$pres_abs), "x", ncol(data_store$pres_abs), "\n")
          if (!is.null(data_store$pres_abs_regular)) {
            cat("Regular-grid matrix available:", nrow(data_store$pres_abs_regular), "x", ncol(data_store$pres_abs_regular), "\n")
          }
          if (!is.null(data_store$pres_abs_irregular)) {
            cat("Irregular-polygon matrix available:", nrow(data_store$pres_abs_irregular), "x", ncol(data_store$pres_abs_irregular), "\n")
          }
          if (method == "occurrence_only") {
            if (identical(points_mode, "irregular_direct")) {
              cat("Matrix basis: occurrence points counted directly in user irregular polygons (no regular grid)\n")
            } else if (identical(points_mode, "irregular_with_grid")) {
              cat("Matrix basis: occurrence points assigned to regular grid cells, then aggregated to user irregular polygons\n")
            } else {
              cat("Matrix basis: occurrence points assigned to regular grid cells\n")
            }
          } else if (isTRUE(data_store$matrix_is_irregular_aggregated)) {
            cat("Matrix basis: regular-grid presence aggregated to user irregular polygons\n")
          }
          if (!is.null(data_store$irregular_bins_richness)) {
            cat("\n✓ Diversity by irregular polygons computed\n")
            cat("  Source method:", data_store$irregular_bins_method, "\n")
            cat("  Subdivision ID column:", data_store$irregular_bins_id_column, "\n")
            cat("  Total subdivisions:", nrow(data_store$irregular_bins_richness), "\n")
          }
          if (isTRUE(data_store$matrix_is_irregular_aggregated)) {
            cat("\n✓ Matrix areas were reduced to irregular subdivisions for downstream exports/BioGeoBEARS\n")
          }
          cat("\nOutput files saved in:\n")
          cat("  - /out_buffers/ (if using buffers)\n")
          cat("  - /out_MCP/ (if using convex hulls)\n")
          cat("  - /out_MST/ (if using MST)\n")
          cat("  - /out_grid/ (regular-grid outputs for grid-based workflows)\n")
          cat("  - /out_irregular_bins/ (if using irregular bins)\n")
        })
        append_extrap_log(paste0("Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      }, error = function(e) {
        append_extrap_log(paste0("[error] ", e$message))
        output$extrap_status <- renderPrint({
          cat("Error in extrapolation:\n", e$message, "\n")
        })
      })
    }
  })

  output$extrap_console_log <- renderText({
    if (is.null(data_store$extrap_log) || length(data_store$extrap_log) == 0) {
      return("Run an extrapolation to see the execution log here.")
    }
    paste(data_store$extrap_log, collapse = "\n")
  })

  output$download_extrap_log <- downloadHandler(
    filename = function() {
      paste0("extrapolation_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      log_lines <- data_store$extrap_log
      if (is.null(log_lines) || length(log_lines) == 0) {
        log_lines <- c(
          "=== Range Extrapolation Log ===",
          "No extrapolation run has been executed yet."
        )
      }
      writeLines(log_lines, con = file, useBytes = TRUE)
    }
  )
  
  # ===== VISUALIZATIONS TAB =====

  observeEvent(input$clear_visual_outputs, {
    removed <- clear_visual_output_files(method = NULL)
    data_store$loaded_shapefiles <- list()
    data_store$loaded_rasters <- list()
    data_store$shapefile_info <- list()
    data_store$regular_grid_all_sf <- NULL
    data_store$regular_grid_presence_sf <- NULL
    data_store$regular_grid_presence_points_sf <- NULL
    data_store$regular_grid_presence_extrap_sf <- NULL
    data_store$regular_grid_presence_points_by_taxon <- list()
    data_store$regular_grid_presence_extrap_by_taxon <- list()
    data_store$regular_grid_presence_by_taxon <- list()
    data_store$regular_grid_res <- NULL
    distribution_visible_groups(NULL)

    output$polygons_status <- renderPrint({
      cat("✓ Visual output layers cleared\n")
      cat("Files removed:", length(removed), "\n")
      cat("Run Step 3 again to regenerate fresh layers.\n")
    })
  })
  
  observeEvent(input$load_polygons, {
    tryCatch({
      overlay_all_methods <- isTRUE(input$viz_overlay_all_methods)
      reload_method <- if (overlay_all_methods) NULL else (data_store$visual_layers_method %||% data_store$extrap_method)
      reload_since_time <- NULL
      reload_taxa <- if (overlay_all_methods) character(0) else data_store$last_extrap_taxa
      output_files <- load_output_shapefiles(
        extrap_method = reload_method,
        since_time = reload_since_time,
        taxa_filter = reload_taxa
      )
      output_files <- filter_visual_outputs(output_files)
      data_store$loaded_shapefiles <- output_files$shapefiles
      data_store$loaded_rasters <- output_files$rasters
      data_store$shapefile_info <- output_files$shapefile_info

      species_from_shapes <- character(0)
      if (!is.null(output_files$shapefiles) && length(output_files$shapefiles) > 0) {
        species_from_shapes <- unique(vapply(seq_along(output_files$shapefiles), function(i) {
          layer_key <- names(output_files$shapefiles)[i]
          layer_info <- output_files$shapefile_info[[layer_key]]
          layer_species_name(layer_info, layer_key)
        }, character(1)))
      }

      species_from_occ <- if (!is.null(data_store$occurrence) && nrow(data_store$occurrence) > 0) {
        unique(as.character(data_store$occurrence$spp))
      } else {
        character(0)
      }

      species_list <- unique(c(species_from_shapes, species_from_occ))
      ensure_species_palette(extra_species = species_list)
      
      output$polygons_status <- renderPrint({
        cat("✓ Output files loaded!\n")
        cat("Shapefiles:", length(output_files$shapefiles), "\n")
        cat("Rasters:", length(output_files$rasters), "\n")
      })
    }, error = function(e) {
      output$polygons_status <- renderPrint({
        cat("Error loading files:\n", e$message, "\n")
      })
    })
  })

  observeEvent(input$viz_overlay_all_methods, {
    tryCatch({
      overlay_all_methods <- isTRUE(input$viz_overlay_all_methods)
      reload_method <- if (overlay_all_methods) NULL else (data_store$visual_layers_method %||% data_store$extrap_method)
      reload_since_time <- NULL
      reload_taxa <- if (overlay_all_methods) character(0) else data_store$last_extrap_taxa
      output_files <- load_output_shapefiles(
        extrap_method = reload_method,
        since_time = reload_since_time,
        taxa_filter = reload_taxa
      )
      output_files <- filter_visual_outputs(output_files)
      data_store$loaded_shapefiles <- output_files$shapefiles
      data_store$loaded_rasters <- output_files$rasters
      data_store$shapefile_info <- output_files$shapefile_info
    }, error = function(e) {
      NULL
    })
  }, ignoreInit = TRUE)

  output$distribution_map_panel <- renderUI({
    leaflet::leafletOutput("distribution_map", height = "600px")
  })

  output$distribution_map_static <- renderPlot({
    req(identical(data_store$extrap_method, "mst"))

    is_mst_node_mode <- ((input$mst_level %||% "") %in% c("node", "node_taxa"))
    poly_opacity <- polygon_opacity_debounced()

    map_occurrence <- data_store$analysis_occurrence %||% data_store$occurrence

    legend_species <- character(0)
    legend_colors <- character(0)

    preferred_taxa <- character(0)
    if (!is.null(data_store$mst_context) && !is.null(data_store$mst_context$taxa)) {
      preferred_taxa <- c(preferred_taxa, as.character(data_store$mst_context$taxa))
    }
    if (!is.null(data_store$pres_abs) && ncol(data_store$pres_abs) > 0) {
      preferred_taxa <- c(preferred_taxa, colnames(data_store$pres_abs))
    }
    if (!is.null(map_occurrence) && nrow(map_occurrence) > 0) {
      preferred_taxa <- c(preferred_taxa, as.character(unique(map_occurrence$spp)))
    }
    preferred_taxa <- unique(preferred_taxa[nzchar(preferred_taxa)])

    if (length(preferred_taxa) > 0) {
      base_colors <- unname(data_store$species_colors[preferred_taxa])
      if (length(base_colors) != length(preferred_taxa)) {
        base_colors <- rep(NA_character_, length(preferred_taxa))
      }
      missing_base <- is.na(base_colors)
      if (any(missing_base)) {
        base_colors[missing_base] <- grDevices::hcl.colors(sum(missing_base), "Dark 3")
      }
      legend_species <- preferred_taxa
      legend_colors <- base_colors
    }

    if (!is.null(data_store$study_area_shapefile)) {
      plot(data_store$study_area_shapefile, col = "gray95", border = "gray35", axes = TRUE)
    } else if (!is.null(map_occurrence) && nrow(map_occurrence) > 0) {
      plot(
        range(map_occurrence$long, na.rm = TRUE),
        range(map_occurrence$lat, na.rm = TRUE),
        type = "n", xlab = "Longitude", ylab = "Latitude", axes = TRUE
      )
    } else {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    }

    if (isTRUE(input$viz_show_grid) && !is.null(data_store$loaded_rasters) && length(data_store$loaded_rasters) > 0) {
      for (i in seq_along(data_store$loaded_rasters)) {
        tryCatch({
          raster_key <- names(data_store$loaded_rasters)[i]
          rr <- raster::raster(data_store$loaded_rasters[[i]])
          raster::plot(rr, add = TRUE, legend = FALSE, col = grDevices::adjustcolor("#FDBB84", alpha.f = 0.45))
        }, error = function(e) NULL)
      }
    }

    if (!is.null(data_store$loaded_shapefiles) && length(data_store$loaded_shapefiles) > 0) {
      for (i in seq_along(data_store$loaded_shapefiles)) {
        tryCatch({
          layer_key <- names(data_store$loaded_shapefiles)[i]
          layer_info <- data_store$shapefile_info[[layer_key]]
          species_name <- layer_species_name(layer_info, layer_key)
          source_name <- if (!is.null(layer_info) && !is.null(layer_info$path)) {
            basename(layer_info$path)
          } else {
            basename(layer_key)
          }
          base_name <- tools::file_path_sans_ext(source_name)
          if (is_helper_points_layer(base_name)) {
            next
          }
          if (is_mst_internal_helper_layer(base_name, species_name)) {
            next
          }
          if (is_mst_onlyterminal_points_layer(base_name, species_name) && !identical(input$mst_level %||% "", "node_taxa")) {
            next
          }
          is_grid_layer <- grepl("^q[0-9]+$", base_name, ignore.case = TRUE) || grepl("^GRIDS_", base_name, ignore.case = TRUE)
          species_color <- unname(data_store$species_colors[species_name])
          if (is.na(species_color)) species_color <- "#FF6B6B"

          shp_sf <- normalize_to_wgs84_sf(data_store$loaded_shapefiles[[i]])
          empty_idx <- tryCatch(sf::st_is_empty(shp_sf), error = function(e) rep(FALSE, nrow(shp_sf)))
          keep_idx <- is.na(empty_idx) | !empty_idx
          if (!all(keep_idx)) shp_sf <- shp_sf[keep_idx, , drop = FALSE]
          if (nrow(shp_sf) == 0) next

          geom_types <- as.character(sf::st_geometry_type(shp_sf, by_geometry = TRUE))
          poly_idx <- !is.na(geom_types) & geom_types %in% c("POLYGON", "MULTIPOLYGON")
          line_idx <- !is.na(geom_types) & geom_types %in% c("LINESTRING", "MULTILINESTRING")
          geom_collection_idx <- !is.na(geom_types) & geom_types %in% c("GEOMETRYCOLLECTION", "GEOMETRY")

          if (any(poly_idx)) {
            plot(sf::st_geometry(shp_sf[poly_idx, , drop = FALSE]), add = TRUE, border = species_color, col = grDevices::adjustcolor(species_color, alpha.f = poly_opacity), lwd = 1)
          }
          if (any(line_idx)) {
            plot(sf::st_geometry(shp_sf[line_idx, , drop = FALSE]), add = TRUE, col = species_color, lwd = 3)
            if (!(species_name %in% legend_species)) {
              legend_species <- c(legend_species, species_name)
              legend_colors <- c(legend_colors, species_color)
            }
          }

          if (any(geom_collection_idx)) {
            extracted_line <- suppressWarnings(
              sf::st_collection_extract(shp_sf[geom_collection_idx, , drop = FALSE], "LINESTRING")
            )
            if (!is.null(extracted_line) && nrow(extracted_line) > 0) {
              plot(sf::st_geometry(extracted_line), add = TRUE, col = species_color, lwd = 3)
              if (!(species_name %in% legend_species)) {
                legend_species <- c(legend_species, species_name)
                legend_colors <- c(legend_colors, species_color)
              }
            }
          }
        }, error = function(e) NULL)
      }
    }

    if (!is.null(map_occurrence) && nrow(map_occurrence) > 0) {
      species_colors <- get_species_color_map(unique(map_occurrence$spp))
      for (sp in unique(map_occurrence$spp)) {
        sp_data <- map_occurrence[map_occurrence$spp == sp, , drop = FALSE]
        sp_color <- unname(species_colors[sp])
        if (is.na(sp_color)) sp_color <- "#2C7FB8"
        points(sp_data$long, sp_data$lat, pch = 21, cex = 0.8, bg = sp_color, col = "black")
      }
    }

    if (length(legend_species) == 0) {
      fallback_taxa <- character(0)

      if (length(preferred_taxa) > 0) {
        fallback_taxa <- preferred_taxa
      }

      if (length(fallback_taxa) == 0 && !is.null(data_store$loaded_shapefiles) && length(data_store$loaded_shapefiles) > 0) {
        fallback_taxa <- unique(vapply(seq_along(data_store$loaded_shapefiles), function(i) {
          layer_key <- names(data_store$loaded_shapefiles)[i]
          layer_info <- data_store$shapefile_info[[layer_key]]
          layer_species_name(layer_info, layer_key)
        }, character(1)))
        fallback_taxa <- fallback_taxa[!grepl("mintreeall|ancterminal", fallback_taxa, ignore.case = TRUE)]
      }

      fallback_taxa <- fallback_taxa[nzchar(fallback_taxa)]
      if (length(fallback_taxa) > 0) {
        palette_colors <- unname(data_store$species_colors[fallback_taxa])
        if (length(palette_colors) != length(fallback_taxa)) {
          palette_colors <- rep(NA_character_, length(fallback_taxa))
        }
        missing_idx <- is.na(palette_colors)
        if (any(missing_idx)) {
          gen_cols <- grDevices::hcl.colors(sum(missing_idx), "Dark 3")
          palette_colors[missing_idx] <- gen_cols
        }
        legend_species <- fallback_taxa
        legend_colors <- palette_colors
      }
    }

    if (length(legend_species) > 0) {
      keep <- !duplicated(legend_species)
      legend(
        x = "bottomleft",
        legend = legend_species[keep],
        col = legend_colors[keep],
        lwd = 3,
        bty = "n",
        cex = 1.05,
        title = "MST tracks"
      )
    }
  })
  
  distribution_visible_groups <- reactiveVal(NULL)
  polygon_opacity_debounced <- shiny::debounce(reactive(input$polygon_opacity %||% 0.5), millis = 250)

  observeEvent(input$distribution_map_groups, {
    distribution_visible_groups(input$distribution_map_groups)
  }, ignoreInit = TRUE)

  output$distribution_map <- renderLeaflet({
    poly_opacity <- polygon_opacity_debounced()
    method_label_map <- c(
      occurrence_only = "Occurrence Only",
      buffer = "Buffer",
      convex_hull = "Convex Hull",
      mst = "MST",
      irregular_bins = "Irregular Bins"
    )
    active_render_method <- data_store$visual_layers_method %||% data_store$extrap_method
    active_method_label <- unname(method_label_map[active_render_method])
    overlay_all_methods <- isTRUE(input$viz_overlay_all_methods)
    is_mst_node_mode <- identical(active_render_method, "mst") && ((input$mst_level %||% "") %in% c("node", "node_taxa"))

    # Create base map
    m <- leaflet::leaflet() %>%
      leaflet::addProviderTiles("OpenStreetMap.Mapnik") %>%
      leaflet::addMapPane("presence-rasters", zIndex = 350) %>%
      leaflet::addMapPane("range-polygons", zIndex = 420) %>%
      leaflet::addMapPane("mst-traces", zIndex = 450) %>%
      leaflet::addMapPane("occurrence-points", zIndex = 460)

    overlay_groups <- c("Study Area")

    selected_species <- character(0)
    map_occurrence <- data_store$analysis_occurrence %||% data_store$occurrence
    
    # Add study area shapefile if available
    if (!is.null(data_store$study_area_shapefile)) {
      tryCatch({
        shapefile_sf <- normalize_to_wgs84_sf(data_store$study_area_shapefile)
        m <- m %>%
          leaflet::addPolygons(
            data = shapefile_sf,
            color = "black",
            weight = 2,
            fillColor = "lightblue",
            fillOpacity = 0.3,
            popup = "Study Area",
            options = leaflet::pathOptions(pane = "range-polygons"),
            group = "Study Area"
          )
      }, error = function(e) {
        warning("Could not add study area shapefile: ", e$message)
      })
    }

    # Add optional secondary study areas (irregular polygons uploaded by user)
    if (!is.null(data_store$secondary_study_areas) && length(data_store$secondary_study_areas) > 0) {
      for (secondary_name in names(data_store$secondary_study_areas)) {
        tryCatch({
          secondary_sf <- normalize_to_wgs84_sf(data_store$secondary_study_areas[[secondary_name]])
          if (nrow(secondary_sf) == 0) next

          group_name <- secondary_name
          m <- m %>%
            leaflet::addPolygons(
              data = secondary_sf,
              color = "#1F78B4",
              weight = 2,
              dashArray = "6,4",
              fill = FALSE,
              popup = secondary_name,
              options = leaflet::pathOptions(pane = "range-polygons"),
              group = group_name
            )

          overlay_groups <- c(overlay_groups, group_name)
        }, error = function(e) {
          warning("Could not add secondary study area shapefile: ", e$message)
        })
      }
    }

    grid_res_label <- data_store$regular_grid_res
    grid_res_txt <- if (!is.null(grid_res_label) && is.finite(grid_res_label)) {
      format(grid_res_label, scientific = FALSE, trim = TRUE)
    } else {
      "custom"
    }
    grid_group_all <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - All cells")
    grid_group_presence <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - Occupied")
    extrap_grid_label <- data_store$regular_grid_presence_extrap_label %||% active_method_label %||% "Extrapolation"
    grid_group_presence_points <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - All taxa - Occurrence points")
    grid_group_presence_extrap <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - All taxa - ", extrap_grid_label)

    if (!is.null(data_store$regular_grid_all_sf) && nrow(data_store$regular_grid_all_sf) > 0) {
      m <- m %>% leaflet::addPolygons(
        data = data_store$regular_grid_all_sf,
        color = "#4d4d4d",
        weight = 1,
        fill = FALSE,
        popup = ~paste("Grid cell:", grid_id),
        options = leaflet::pathOptions(pane = "range-polygons"),
        group = grid_group_all
      )
      overlay_groups <- c(overlay_groups, grid_group_all)
    }

    if (!is.null(data_store$regular_grid_presence_extrap_sf) && nrow(data_store$regular_grid_presence_extrap_sf) > 0) {
      m <- m %>% leaflet::addPolygons(
        data = data_store$regular_grid_presence_extrap_sf,
        color = "#1f9e89",
        weight = 1,
        fillColor = "#1f9e89",
        fillOpacity = 0.15,
        popup = ~paste0(extrap_grid_label, " occupied grid cell: ", grid_id),
        options = leaflet::pathOptions(pane = "range-polygons"),
        group = grid_group_presence_extrap
      )
      overlay_groups <- c(overlay_groups, grid_group_presence_extrap)
    } else if (!is.null(data_store$regular_grid_presence_sf) && nrow(data_store$regular_grid_presence_sf) > 0) {
      m <- m %>% leaflet::addPolygons(
        data = data_store$regular_grid_presence_sf,
        color = "#1f9e89",
        weight = 1,
        fillColor = "#1f9e89",
        fillOpacity = 0.15,
        popup = ~paste("Occupied grid cell:", grid_id),
        options = leaflet::pathOptions(pane = "range-polygons"),
        group = grid_group_presence
      )
      overlay_groups <- c(overlay_groups, grid_group_presence)
    }

    if (!is.null(data_store$regular_grid_presence_points_sf) && nrow(data_store$regular_grid_presence_points_sf) > 0) {
      m <- m %>% leaflet::addPolygons(
        data = data_store$regular_grid_presence_points_sf,
        color = "#355C7D",
        weight = 1,
        fillColor = "#355C7D",
        fillOpacity = 0.10,
        popup = ~paste("Occurrence-point occupied grid cell:", grid_id),
        options = leaflet::pathOptions(pane = "range-polygons"),
        group = grid_group_presence_points
      )
      overlay_groups <- c(overlay_groups, grid_group_presence_points)
    }

    if (!is.null(data_store$regular_grid_presence_by_taxon) && length(data_store$regular_grid_presence_by_taxon) > 0) {
      taxa_grid_names <- names(data_store$regular_grid_presence_by_taxon)
      taxa_grid_colors <- get_species_color_map(taxa_grid_names)
      for (tx in taxa_grid_names) {
        tx_sf <- data_store$regular_grid_presence_by_taxon[[tx]]
        if (is.null(tx_sf) || nrow(tx_sf) == 0) next
        tx_col <- unname(taxa_grid_colors[tx])
        if (is.na(tx_col) || !nzchar(tx_col)) tx_col <- "#2C7FB8"
        tx_label <- if (!is.null(data_store$regular_grid_presence_extrap_by_taxon) && length(data_store$regular_grid_presence_extrap_by_taxon) > 0) {
          paste0(tx, " (", extrap_grid_label, ")")
        } else {
          paste0(tx, " (presence)")
        }
        tx_group <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - ", tx_label)
        m <- m %>% leaflet::addPolygons(
          data = tx_sf,
          color = tx_col,
          weight = 1,
          fillColor = tx_col,
          fillOpacity = 0.2,
          popup = ~paste0("Taxon: ", tx, "<br>Occupied grid cell: ", grid_id),
          options = leaflet::pathOptions(pane = "range-polygons"),
          group = tx_group
        )
        overlay_groups <- c(overlay_groups, tx_group)
      }
    }

    if (!is.null(data_store$regular_grid_presence_points_by_taxon) && length(data_store$regular_grid_presence_points_by_taxon) > 0) {
      taxa_grid_names <- names(data_store$regular_grid_presence_points_by_taxon)
      taxa_grid_colors <- get_species_color_map(taxa_grid_names)
      for (tx in taxa_grid_names) {
        tx_sf <- data_store$regular_grid_presence_points_by_taxon[[tx]]
        if (is.null(tx_sf) || nrow(tx_sf) == 0) next
        tx_col <- unname(taxa_grid_colors[tx])
        if (is.na(tx_col) || !nzchar(tx_col)) tx_col <- "#355C7D"
        tx_group <- paste0("Grid ", grid_res_txt, " x ", grid_res_txt, " deg - ", tx, " (occurrence points)")
        m <- m %>% leaflet::addPolygons(
          data = tx_sf,
          color = tx_col,
          weight = 1,
          fillColor = tx_col,
          fillOpacity = 0.12,
          popup = ~paste0("Taxon: ", tx, "<br>Occurrence-point grid cell: ", grid_id),
          options = leaflet::pathOptions(pane = "range-polygons"),
          group = tx_group
        )
        overlay_groups <- c(overlay_groups, tx_group)
      }
    }

    # Add extrapolation polygons with species names and dynamic opacity
    if (!is.null(data_store$loaded_shapefiles) && length(data_store$loaded_shapefiles) > 0) {
      for (i in seq_along(data_store$loaded_shapefiles)) {
        tryCatch({
          shp_sf <- normalize_to_wgs84_sf(data_store$loaded_shapefiles[[i]])
          empty_idx <- tryCatch(sf::st_is_empty(shp_sf), error = function(e) rep(FALSE, nrow(shp_sf)))
          keep_idx <- is.na(empty_idx) | !empty_idx
          if (!all(keep_idx)) {
            shp_sf <- shp_sf[keep_idx, , drop = FALSE]
          }
          if (nrow(shp_sf) == 0) {
            next
          }
          
          layer_key <- names(data_store$loaded_shapefiles)[i]
          layer_info <- data_store$shapefile_info[[layer_key]]
          if (!overlay_all_methods && !is.null(active_method_label) && !is.null(layer_info) && !is.null(layer_info$method) && !identical(layer_info$method, active_method_label)) {
            next
          }
          species_name <- layer_species_name(layer_info, layer_key)
          source_name <- if (!is.null(layer_info) && !is.null(layer_info$path)) {
            basename(layer_info$path)
          } else {
            basename(layer_key)
          }
          base_name <- tools::file_path_sans_ext(source_name)
          if (is_helper_points_layer(base_name)) {
            next
          }
          if (is_mst_internal_helper_layer(base_name, species_name)) {
            next
          }
          if (is_mst_onlyterminal_points_layer(base_name, species_name) && !identical(input$mst_level %||% "", "node_taxa")) {
            next
          }
          is_grid_layer <- grepl("^q[0-9]+$", base_name, ignore.case = TRUE) || grepl("^GRIDS_", base_name, ignore.case = TRUE)
          
          # Get color for this species from pre-calculated palette
          species_color <- unname(data_store$species_colors[species_name])
          if (is.na(species_color)) {
            species_color <- "#FF6B6B"  # Default color
          }
          
          # Get opacity from slider
          opacity <- poly_opacity
          if (is.null(opacity)) opacity <- 0.5
          
          is_summary_layer <- grepl("species_richness_bins", basename(layer_key), fixed = TRUE)
          if (length(selected_species) == 0 || species_name %in% selected_species || is_summary_layer) {
            geom_types <- as.character(sf::st_geometry_type(shp_sf, by_geometry = TRUE))

            poly_idx <- !is.na(geom_types) & geom_types %in% c("POLYGON", "MULTIPOLYGON")
            line_idx <- !is.na(geom_types) & geom_types %in% c("LINESTRING", "MULTILINESTRING")
            point_idx <- !is.na(geom_types) & geom_types %in% c("POINT", "MULTIPOINT")

            if (any(poly_idx)) {
              polygon_group <- if (is_grid_layer) {
                paste0(base_name, " - Grid")
              } else {
                paste0(species_name, " - Polygon")
              }
              m <- m %>%
                leaflet::addPolygons(
                  data = shp_sf[poly_idx, , drop = FALSE],
                  color = if (is_grid_layer) "#2b2b2b" else species_color,
                  weight = if (is_grid_layer) 0.8 else 1,
                  fillColor = if (is_grid_layer) "#FFFFFF" else species_color,
                  fillOpacity = if (is_grid_layer) 0 else opacity,
                  popup = if (is_grid_layer) paste("Grid layer:", base_name) else paste("Species:", species_name),
                  options = leaflet::pathOptions(pane = "range-polygons"),
                  group = polygon_group
                )
              overlay_groups <- c(overlay_groups, polygon_group)
            }

            if (any(line_idx)) {
              mst_group <- paste0(species_name, " - Tracks")
              m <- m %>%
                leaflet::addPolylines(
                  data = shp_sf[line_idx, , drop = FALSE],
                  color = species_color,
                  weight = 3,
                  opacity = pmax(opacity, 0.6),
                  popup = paste("MST Track:", species_name),
                  options = leaflet::pathOptions(pane = "mst-traces"),
                  group = mst_group
                )
              overlay_groups <- c(overlay_groups, mst_group)
            }

            if (any(point_idx)) {
              point_group <- paste0(species_name, " - Points")
              m <- m %>%
                leaflet::addCircleMarkers(
                  data = shp_sf[point_idx, , drop = FALSE],
                  radius = 4,
                  color = "black",
                  weight = 1,
                  fillColor = species_color,
                  fillOpacity = pmax(opacity, 0.5),
                  popup = paste("Point Layer:", species_name),
                  options = leaflet::pathOptions(pane = "occurrence-points"),
                  group = point_group
                )
              overlay_groups <- c(overlay_groups, point_group)
            }

            geom_collection_idx <- !is.na(geom_types) & geom_types %in% c("GEOMETRYCOLLECTION", "GEOMETRY")
            if (any(geom_collection_idx)) {
              geom_collection_sf <- shp_sf[geom_collection_idx, , drop = FALSE]

              extracted_poly <- suppressWarnings(sf::st_collection_extract(geom_collection_sf, "POLYGON"))
              if (!is.null(extracted_poly) && nrow(extracted_poly) > 0) {
                polygon_group <- if (is_grid_layer) {
                  paste0(base_name, " - Grid")
                } else {
                  paste0(species_name, " - Polygon")
                }
                m <- m %>%
                  leaflet::addPolygons(
                    data = extracted_poly,
                    color = if (is_grid_layer) "#2b2b2b" else species_color,
                    weight = if (is_grid_layer) 0.8 else 1,
                    fillColor = if (is_grid_layer) "#FFFFFF" else species_color,
                    fillOpacity = if (is_grid_layer) 0 else opacity,
                    popup = if (is_grid_layer) paste("Grid layer:", base_name) else paste("Species:", species_name),
                    options = leaflet::pathOptions(pane = "range-polygons"),
                    group = polygon_group
                  )
                overlay_groups <- c(overlay_groups, polygon_group)
              }

              extracted_line <- suppressWarnings(sf::st_collection_extract(geom_collection_sf, "LINESTRING"))
              if (!is.null(extracted_line) && nrow(extracted_line) > 0) {
                mst_group <- paste0(species_name, " - Tracks")
                m <- m %>%
                  leaflet::addPolylines(
                    data = extracted_line,
                    color = species_color,
                    weight = 3,
                    opacity = pmax(opacity, 0.6),
                    popup = paste("MST Track:", species_name),
                    options = leaflet::pathOptions(pane = "mst-traces"),
                    group = mst_group
                  )
                overlay_groups <- c(overlay_groups, mst_group)
              }

              extracted_point <- suppressWarnings(sf::st_collection_extract(geom_collection_sf, "POINT"))
              if (!is.null(extracted_point) && nrow(extracted_point) > 0) {
                point_group <- paste0(species_name, " - Points")
                m <- m %>%
                  leaflet::addCircleMarkers(
                    data = extracted_point,
                    radius = 4,
                    color = "black",
                    weight = 1,
                    fillColor = species_color,
                    fillOpacity = pmax(opacity, 0.5),
                    popup = paste("Point Layer:", species_name),
                    options = leaflet::pathOptions(pane = "occurrence-points"),
                    group = point_group
                  )
                overlay_groups <- c(overlay_groups, point_group)
              }
            }
          }
        }, error = function(e) {
          warning("Could not add map layer ", i, ": ", e$message)
        })
      }
    }

    if (isTRUE(input$viz_show_grid) && !is.null(data_store$loaded_rasters) && length(data_store$loaded_rasters) > 0) {
      for (i in seq_along(data_store$loaded_rasters)) {
        tryCatch({
          raster_key <- names(data_store$loaded_rasters)[i]
          rast_name <- basename(raster_key)

          # Ignore legacy helper rasters (q5, q10, etc.)
          if (grepl("_q[0-9]+", rast_name, ignore.case = TRUE)) {
            next
          }

          layer_info <- data_store$shapefile_info[[raster_key]]
          if (!overlay_all_methods && !is.null(active_method_label) && !is.null(layer_info) && !is.null(layer_info$method) && !identical(layer_info$method, active_method_label)) {
            next
          }
          raster_species <- if (!is.null(layer_info) && !is.null(layer_info$species_id) && !is.na(layer_info$species_id)) layer_info$species_id else NA_character_
          if (!is.na(raster_species) && length(selected_species) > 0 && !(raster_species %in% selected_species)) {
            next
          }
          rr <- raster::raster(data_store$loaded_rasters[[i]])
          vals <- raster::values(rr)
          vals <- vals[is.finite(vals)]
          if (length(vals) > 0) {
            pal <- leaflet::colorNumeric("YlOrRd", domain = vals, na.color = "transparent")
            raster_group <- if (is_mst_node_mode) {
              "mintreeall_ancterminal - Presence raster"
            } else if (!is.na(raster_species)) {
              paste0(raster_species, " - Presence raster")
            } else {
              paste0(rast_name, " - Raster")
            }
            m <- m %>% leaflet::addRasterImage(
              rr,
              colors = pal,
              opacity = 0.45,
              options = leaflet::gridOptions(zIndex = 1, pane = "presence-rasters"),
              group = raster_group
            )
            overlay_groups <- c(overlay_groups, raster_group)
          }
        }, error = function(e) {
          warning("Could not add raster ", i, ": ", e$message)
        })
      }
    }
    
    # Add occurrence points with colors by species
    if (!is.null(map_occurrence) && nrow(map_occurrence) > 0) {
      # Build a robust per-species color map for the current occurrence layer
      species_colors <- get_species_color_map(unique(map_occurrence$spp))
      
      # Add points by species to allow toggling them with the polygons
      for (sp in unique(map_occurrence$spp)) {
        if (length(selected_species) > 0 && !(sp %in% selected_species)) {
          next
        }
        sp_data <- map_occurrence[map_occurrence$spp == sp, ]
        sp_color <- unname(species_colors[sp])
        if (is.na(sp_color) || !nzchar(sp_color)) {
          sp_color <- "#2C7FB8"
        }
        point_group <- paste0(sp, " - Occurrence points")
        
        m <- m %>%
          leaflet::addCircleMarkers(
            data = sp_data,
            lng = ~long,
            lat = ~lat,
            radius = 5,
            color = "black",
            weight = 1,
            fillColor = sp_color,
            fillOpacity = 0.8,
            popup = ~paste("Species:", spp),
            options = leaflet::pathOptions(pane = "occurrence-points"),
            group = point_group
          )
        overlay_groups <- c(overlay_groups, point_group)
      }
    }

    overlay_groups <- unique(overlay_groups)
    
    # Add layer control
    m <- m %>%
      leaflet::addLayersControl(
        overlayGroups = overlay_groups,
        options = leaflet::layersControlOptions(collapsed = FALSE)
      )

    saved_groups <- isolate(distribution_visible_groups())
    if (!is.null(saved_groups) && length(saved_groups) > 0) {
      keep_groups <- intersect(overlay_groups, saved_groups)
      if (length(keep_groups) > 0) {
        hide_groups <- setdiff(overlay_groups, keep_groups)
        if (length(hide_groups) > 0) {
          m <- m %>% leaflet::hideGroup(hide_groups)
        }
      }
    }

    # Fit bounds to show all data
    if (!is.null(map_occurrence) && nrow(map_occurrence) > 0) {
      m <- m %>%
        leaflet::fitBounds(
          lng1 = min(map_occurrence$long, na.rm = TRUE),
          lat1 = min(map_occurrence$lat, na.rm = TRUE),
          lng2 = max(map_occurrence$long, na.rm = TRUE),
          lat2 = max(map_occurrence$lat, na.rm = TRUE)
        )
    }
    
    m
  })

  output$irregular_bins_map <- leaflet::renderLeaflet({
    req(!is.null(data_store$irregular_bins_richness))
    bins_sf <- normalize_to_wgs84_sf(data_store$irregular_bins_richness)
    vals <- bins_sf$n_species
    pal <- leaflet::colorNumeric("viridis", domain = vals, na.color = "transparent")

    leaflet::leaflet(bins_sf) %>%
      leaflet::addProviderTiles("CartoDB.Positron") %>%
      leaflet::addPolygons(
        fillColor = ~pal(n_species),
        fillOpacity = 0.8,
        color = "#333333",
        weight = 1,
        popup = ~paste0(
          "Richness (n_species): ", n_species,
          "<br>Species: ", ifelse(is.na(species_list) | species_list == "", "none", species_list)
        )
      ) %>%
      leaflet::addLegend(
        position = "bottomright",
        pal = pal,
        values = vals,
        title = "Species richness"
      )
  })

  output$irregular_bins_table <- DT::renderDataTable({
    req(!is.null(data_store$irregular_bins_species_table))
    DT::datatable(data_store$irregular_bins_species_table, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  
  output$phylo_tree_plot <- renderPlot({
    if (is.null(data_store$tree)) {
      plot(1, type = "n", main = "No tree loaded")
    } else {
      plot(data_store$tree, main = "Phylogenetic Tree")
    }
  })
  
  # ===== PAE-PCE ANALYSIS TAB =====
  
  # Reactive value to store PAE-PCE results
  pae_results <- reactiveVal(NULL)
  pae_input_matrix_snapshot <- reactiveVal(NULL)
  
  observeEvent(input$run_pae_pce, {
    # Check if using custom data
    use_custom <- input$pae_use_custom_data
    
    if (use_custom) {
      # Validate custom data inputs
      if (is.null(input$pae_custom_matrix) || is.null(input$pae_custom_shapefile)) {
        output$pae_pce_log <- renderText({
          "Error: Please upload both the presence-absence matrix and shapefile for custom data analysis."
        })
        return()
      }
    } else {
      # Check if MST extrapolation has been run
        if (is.null(data_store$mst_context) || is.null(data_store$extrap_method) || data_store$extrap_method != "mst") {
          output$pae_pce_log <- renderText({
            "Error: PAE-PCE analysis requires a presence-absence matrix generated by the Minimum Spanning Tree (MST) method.\nPlease go to Step 3, select 'Minimum Spanning Tree' and run the extrapolation first."
          })
          return()
        }
        
        # Check if shapefile is loaded
        if (is.null(data_store$mst_context$shapeFile)) {
          output$pae_pce_log <- renderText({
            "Error: No study area shapefile loaded. Please load a shapefile in Step 3."
          })
          return()
        }
    }
    
    # Show running status
    output$pae_pce_log <- renderText({
      "Running PAE-PCE analysis... This may take a while depending on the number of iterations and species."
    })
    
    # Run the analysis in a tryCatch block to handle errors gracefully
    tryCatch({
      # Get parameters - distinguish between custom and MST data
      if (use_custom) {
        # Load custom matrix
        mat_raw <- read.csv(input$pae_custom_matrix$datapath, row.names = 1)
        
        # Load custom shapefile
        shp_files <- input$pae_custom_shapefile
        shp_path <- shp_files$datapath[grep("\\.shp$", shp_files$name)]
        if (length(shp_path) == 0) {
          stop("Shapefile (.shp) not found in uploaded files")
        }
        shape <- rgdal::readOGR(shp_path, verbose = FALSE)
        res <- c(input$pae_custom_resol_x, input$pae_custom_resol_y)
      } else {
        # Use MST context data
        mat_raw <- data_store$mst_context$preabsMat
        shape <- data_store$mst_context$shapeFile
        res <- data_store$mst_context$resol
      }
      n_iter <- input$pae_n_iterations

      mat <- as.matrix(mat_raw)
      mat_num <- suppressWarnings(matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)))
      non_zero_rows <- rowSums(mat_num != 0, na.rm = TRUE) > 0
      is_root_row <- !is.null(rownames(mat)) & rownames(mat) == "ROOT"
      keep_rows <- non_zero_rows | is_root_row

      removed_rows <- rownames(mat)[!keep_rows]
      removed_rows <- removed_rows[!is.na(removed_rows)]
      mat <- mat[keep_rows, , drop = FALSE]

      filter_note <- paste0(
        "PAE input filter: removed ",
        sum(!keep_rows),
        " grid row(s) with all-zero species values",
        ifelse(length(removed_rows) > 0,
               paste0(" [", paste(head(removed_rows, 10), collapse = ", "), ifelse(length(removed_rows) > 10, ", ...", ""), "]"),
               "")
      )

      pae_input_matrix_snapshot(list(
        matrix = mat,
        method = data_store$extrap_method,
        timestamp = Sys.time(),
        n_rows_before = nrow(mat_raw),
        n_rows = nrow(mat),
        n_cols = ncol(mat),
        has_root = "ROOT" %in% rownames(mat),
        removed_zero_rows = sum(!keep_rows),
        removed_zero_row_names = removed_rows
      ))
      
      # Optional zoom parameters
      xmin <- if(!is.na(input$pae_xmin)) input$pae_xmin else NULL
      xmax <- if(!is.na(input$pae_xmax)) input$pae_xmax else NULL
      ymin <- if(!is.na(input$pae_ymin)) input$pae_ymin else NULL
      ymax <- if(!is.na(input$pae_ymax)) input$pae_ymax else NULL
      
      # Capture plot output
      output$pae_pce_plot <- renderPlot({
        # We need to call the function here to generate the plot
        # But we also want to capture the return value
        # So we'll use a local variable and then update the reactive value
        
        # Temporarily redirect output to capture log
        log_file <- tempfile()
        sink(log_file)
        
        result <- tryCatch({
          pae_fun <- get_pae_pce_function()
          pae_args <- list(
            preabsMat = mat,
            shapeFile = shape,
            resol = res,
            N = n_iter,
            gridView = input$pae_grid_view,
            labelGrid = input$pae_label_grid,
            nonHomoplasticSpeciesList = TRUE,  # Always generate species list
            random_seed = input$pae_seed,
            sobrepo = input$pae_sobrepo,
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax
          )

          do.call(pae_fun, pae_args)
        }, error = function(e) {
          sink()
          stop(e)
        })
        
        sink()
        
        # Read log
        log_content <- readLines(log_file)
        output$pae_pce_log <- renderText({
          paste(c(filter_note, log_content), collapse = "\n")
        })
        
        # Store result
        pae_results(result)
      })
      
    }, error = function(e) {
      output$pae_pce_log <- renderText({
        paste("Error during PAE-PCE analysis:", e$message)
      })
    })
  })
  
  # Render the results table
  output$pae_pce_table <- DT::renderDataTable({
    res <- pae_results()
    if (is.null(res)) return(NULL)

    # Check if it's a data frame (successful result) or a list (elegant termination)
    if (is.data.frame(res)) {
      DT::datatable(res, options = list(pageLength = 10, scrollX = TRUE))
    } else if (is.list(res) && !is.null(res$nonHomoplastic_species)) {
      assembled <- do.call(rbind, lapply(seq_along(res$nonHomoplastic_species), function(i) {
        x <- res$nonHomoplastic_species[[i]]
        if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NULL)
        x$iteration <- i
        x
      }))
      if (is.null(assembled) || nrow(assembled) == 0) {
        DT::datatable(data.frame(Message = "PAE-PCE completed but no non-homoplastic table was generated."), options = list(dom = 't'))
      } else {
        DT::datatable(assembled, options = list(pageLength = 10, scrollX = TRUE))
      }
    } else if (is.list(res) && !is.null(res$status)) {
      # It's an elegant termination object
      df <- data.frame(
        Status = res$status,
        Reason = res$reason,
        Iteration = res$iteration_parada
      )
      DT::datatable(df, options = list(dom = 't'))
    }
  })

  output$pae_input_matrix_info <- renderPrint({
    snap <- pae_input_matrix_snapshot()
    if (is.null(snap) || is.null(snap$matrix)) {
      cat("No PAE-PCE run matrix snapshot available yet.\n")
      cat("Run PAE-PCE once, then download the exact matrix used in that run.\n")
      return(invisible(NULL))
    }

    mat <- snap$matrix
    cat("PAE input matrix snapshot\n")
    cat("Timestamp:", format(snap$timestamp), "\n")
    cat("Source method:", snap$method %||% "unknown", "\n")
    if (!is.null(snap$n_rows_before)) {
      cat("Rows before filter:", snap$n_rows_before, "\n")
    }
    cat("Dimensions:", nrow(mat), "rows x", ncol(mat), "columns\n")
    if (!is.null(snap$removed_zero_rows)) {
      cat("Removed all-zero grid rows:", snap$removed_zero_rows, "\n")
    }
    cat("Contains ROOT row:", ifelse("ROOT" %in% rownames(mat), "yes", "no"), "\n")

    if (ncol(mat) > 0) {
      preview_cols <- head(colnames(mat), 8)
      cat("First species columns:", paste(preview_cols, collapse = ", "), "\n")
    }
  })

  output$download_pae_input_matrix <- downloadHandler(
    filename = function() {
      snap <- pae_input_matrix_snapshot()
      ts <- if (!is.null(snap) && !is.null(snap$timestamp)) {
        format(snap$timestamp, "%Y%m%d_%H%M%S")
      } else {
        format(Sys.time(), "%Y%m%d_%H%M%S")
      }
      paste0("pae_input_matrix_snapshot_", ts, ".txt")
    },
    content = function(file) {
      snap <- pae_input_matrix_snapshot()
      if (is.null(snap) || is.null(snap$matrix)) {
        writeLines("No PAE-PCE matrix snapshot available. Run PAE-PCE first.", con = file)
        return(invisible(NULL))
      }

      mat <- as.matrix(snap$matrix)
      write.table(mat, file = file, sep = "\t", quote = FALSE, col.names = NA)
    }
  )
  
  # Download handler for the table
  output$download_pae_table <- downloadHandler(
    filename = function() {
      "pae_pce_results.csv"
    },
    content = function(file) {
      res <- pae_results()
      if (!is.null(res) && is.data.frame(res)) {
        write.csv(res, file, row.names = FALSE)
      } else if (!is.null(res) && is.list(res) && !is.null(res$nonHomoplastic_species)) {
        assembled <- do.call(rbind, lapply(seq_along(res$nonHomoplastic_species), function(i) {
          x <- res$nonHomoplastic_species[[i]]
          if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NULL)
          x$iteration <- i
          x
        }))
        if (is.null(assembled)) {
          assembled <- data.frame(Message = "No non-homoplastic species table generated")
        }
        write.csv(assembled, file, row.names = FALSE)
      } else if (!is.null(res) && is.list(res) && !is.null(res$status)) {
        # Create a simple summary for elegant termination
        df <- data.frame(
          Status = res$status,
          Reason = res$reason,
          Iteration = res$iteration_parada
        )
        write.csv(df, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Error = "No results available"), file, row.names = FALSE)
      }
    }
  )
  
  # ===== EXPORT FILES TAB =====
  area_code_mapping_snapshot <- reactiveVal(NULL)

  resolve_matrix_for_export <- function(preference = NULL) {
    pref <- preference %||% (input$matrix_export_basis %||% "current")

    current_mat <- data_store$pres_abs
    regular_mat <- data_store$pres_abs_regular
    irregular_mat <- data_store$pres_abs_irregular

    selected <- current_mat
    selected_basis <- "current"
    notice <- NULL

    if (identical(pref, "regular")) {
      if (!is.null(regular_mat)) {
        selected <- regular_mat
        selected_basis <- "regular"
      } else {
        notice <- "Regular-grid matrix is not available for the last run; using current matrix instead."
      }
    } else if (identical(pref, "irregular")) {
      if (!is.null(irregular_mat)) {
        selected <- irregular_mat
        selected_basis <- "irregular"
      } else {
        notice <- "Irregular-polygon matrix is not available for the last run; using current matrix instead."
      }
    }

    list(
      matrix = selected,
      basis = selected_basis,
      notice = notice,
      is_irregular = if (identical(selected_basis, "irregular")) {
        TRUE
      } else if (identical(selected_basis, "regular")) {
        FALSE
      } else {
        isTRUE(data_store$matrix_is_irregular_aggregated) ||
          identical(data_store$extrap_method, "irregular_bins") ||
          (identical(data_store$extrap_method, "occurrence_only") &&
             ((input$points_occurrence_mode %||% "") %in% c("irregular_direct", "irregular_with_grid")))
      },
      available_regular = !is.null(regular_mat),
      available_irregular = !is.null(irregular_mat)
    )
  }

  output$matrix_export_status <- renderPrint({
    info <- resolve_matrix_for_export(input$matrix_export_basis)
    cat("Matrix export selector\n")
    cat("- Current matrix available:", !is.null(data_store$pres_abs), "\n")
    cat("- Regular-grid matrix available:", info$available_regular, "\n")
    cat("- Irregular-polygon matrix available:", info$available_irregular, "\n")
    cat("- Active basis for downloads:", info$basis, "\n")
    if (!is.null(info$notice) && nzchar(info$notice)) {
      cat("- Note:", info$notice, "\n")
    }
  })

  output$area_code_mapping_status <- renderPrint({
    snap <- area_code_mapping_snapshot()
    if (is.null(snap) || is.null(snap$mapping)) {
      cat("No area-code mapping generated yet.\n")
      cat("Export BioGeoBEARS/TNT/NEXUS to generate abbreviations.\n")
      return(invisible(NULL))
    }

    cat("Latest mapping format:", snap$format, "\n")
    cat("Generated at:", format(snap$timestamp), "\n")
    cat("Areas mapped:", nrow(snap$mapping), "\n")
    if (nrow(snap$mapping) > 0) {
      preview <- head(snap$mapping, 8)
      cat("Preview (code => area):\n")
      for (i in seq_len(nrow(preview))) {
        area_txt <- if ("label" %in% names(preview)) preview$label[i] else preview$original[i]
        cat("  ", preview$code[i], "=>", area_txt)
        if ("label" %in% names(preview) && !is.na(preview$original[i]) && nzchar(preview$original[i]) && preview$original[i] != area_txt) {
          cat(" [matrix id:", preview$original[i], "]")
        }
        cat("\n")
      }
      if (nrow(snap$mapping) > 8) {
        cat("  ...\n")
      }
    }
  })

  output$download_area_code_mapping <- downloadHandler(
    filename = function() {
      snap <- area_code_mapping_snapshot()
      ts <- if (!is.null(snap) && !is.null(snap$timestamp)) format(snap$timestamp, "%Y%m%d_%H%M%S") else format(Sys.time(), "%Y%m%d_%H%M%S")
      fmt <- if (!is.null(snap) && !is.null(snap$format)) snap$format else "export"
      paste0("area_code_mapping_", fmt, "_", ts, ".csv")
    },
    content = function(file) {
      snap <- area_code_mapping_snapshot()
      if (is.null(snap) || is.null(snap$mapping) || nrow(snap$mapping) == 0) {
        write.csv(data.frame(message = "No mapping available yet. Export BioGeoBEARS/TNT/NEXUS first."), file, row.names = FALSE)
        return(invisible(NULL))
      }
      write.csv(snap$mapping, file, row.names = FALSE)
    }
  )
  
  output$download_biogeobears <- downloadHandler(
    filename = function() { 
      method <- data_store$extrap_method
      if (is.null(method)) method <- "MST"
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      
      method_code <- switch(method,
        "occurrence_only" = "PTS",
        "buffer" = "BUFF",
        "convex_hull" = "MPC",
        "mst" = "MST",
        "irregular_bins" = "BINS",
        "MST"
      )
      
      paste0("pres_abs_", method_code, "_", matrix_sel$basis, "_geog.data") 
    },
    content = function(file) {
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      if (is.null(matrix_sel$matrix)) {
        writeLines("Error: No presence-absence matrix available", file)
      } else {
        tryCatch({
          method <- data_store$extrap_method
          if (is.null(method)) method <- "MST"
          
          method_code <- switch(method,
            "occurrence_only" = "PTS",
            "buffer" = "BUFF",
            "convex_hull" = "MPC",
            "mst" = "MST",
            "irregular_bins" = "BINS",
            "MST"
          )
          
          # We need to capture the output of range_BioGeoBears since it writes to a file
          # but we want to return it as a download
          
          # Create a temporary file
          temp_file <- tempfile()
          
          # Call the function (it will write to its own hardcoded path or the current dir)
          # We'll just use the logic from the function directly to write to our file
          
          mat <- matrix_sel$matrix
          
          # Remove ROOT (last row) if present
          if ("ROOT" %in% rownames(mat)) {
            mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
          } else if (nrow(mat) > 0) {
            # Sometimes ROOT is just the last row without name
            # Let's be safe and check if the last row is all zeros
            if (all(mat[nrow(mat), ] == 0)) {
              mat <- mat[-nrow(mat), , drop = FALSE]
            }
          }
          
          # Keep empty rows for regular grids; drop them for irregular areas
          if (isTRUE(matrix_sel$is_irregular)) {
            row_sums <- rowSums(mat)
            mat <- mat[row_sums > 0, , drop = FALSE]
          }
          
          # Abbreviate area names for compatibility (A, B, C, ...)
          area_codes <- abbreviate_area_labels(row.names(mat), label_lookup = get_area_label_lookup())
          row.names(mat) <- area_codes$labels
          taxa <- row.names(mat)
          n.cara <- ncol(mat)
          n.taxa <- nrow(mat)
          
          # Replace dots with underscores in species names
          colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))
          
          area_code_mapping_snapshot(list(
            mapping = area_codes$mapping,
            format = "BioGeoBEARS",
            timestamp = Sys.time()
          ))

          # Write to file
          cat(paste(n.cara, n.taxa, paste0('(', paste(taxa, collapse = ' '), ')'), sep = ' '), file = file, sep = '\n')
          
          for(i in 1:n.cara){
            cat(paste(colnames(mat)[i], paste(mat[,i], collapse = ''), sep = '\t'), file = file, sep = '\n', append = TRUE)
          }
          
        }, error = function(e) {
          writeLines(paste("Error generating BioGeoBEARS file:", e$message), file)
        })
      }
    }
  )
  
  output$download_tnt_matrix <- downloadHandler(
    filename = function() {
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      paste0("pres_abs_", matrix_sel$basis, ".tnt")
    },
    content = function(file) {
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      if (is.null(matrix_sel$matrix)) {
        writeLines("Error: No presence-absence matrix available", file)
      } else {
        tryCatch({
          # Ensure out_TNT directory exists
          dir.create(file.path(getwd(), "out_TNT"), showWarnings = FALSE)
          
          # Extract logic from tnt_matrix to write directly to file
          mat <- matrix_sel$matrix
          
          # TNT requires an outgroup row; preserve ROOT if present or create it
          if (!("ROOT" %in% rownames(mat))) {
            mat <- rbind(mat, ROOT = rep(0, ncol(mat)))
          }
          
          # Keep empty rows for regular grids; drop them for irregular areas
          if (isTRUE(matrix_sel$is_irregular)) {
            row_sums <- rowSums(mat)
            mat <- mat[row_sums > 0, , drop = FALSE]
          }
          
          # Abbreviate area names for compatibility (A, B, C, ...), preserving ROOT
          area_codes <- abbreviate_area_labels(row.names(mat), label_lookup = get_area_label_lookup())
          row.names(mat) <- area_codes$labels
          taxa <- row.names(mat)
          n.cara <- ncol(mat)
          n.taxa <- nrow(mat)
          
          # Replace dots and spaces with underscores in species names
          colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))
          colnames(mat) <- gsub(pattern = ' ', replacement = '_', x = colnames(mat))
          
          area_code_mapping_snapshot(list(
            mapping = area_codes$mapping,
            format = "TNT",
            timestamp = Sys.time()
          ))

          # Write TNT format to the download file
          cat('nstates num 32;', file = file, sep = '\n')
          cat('xread', file = file, sep = '\n', append = TRUE)
          cat(paste(n.cara, n.taxa, sep = ' '), file = file, sep = '\n', append = TRUE)
          
          for(i in 1:n.taxa){
            cat(paste(taxa[i], paste(mat[i,], collapse = ''), sep = ' '), file = file, sep = '\n', append = TRUE)
          }
          
          cat(';', file = file, sep = '\n', append = TRUE)
          cat('cnames', file = file, sep = '\n', append = TRUE)
          
          for(i in 1:n.cara){
            cat(paste0('{', i-1, ' ', colnames(mat)[i], ';'), file = file, sep = '\n', append = TRUE)
          }
          
          cat(';', file = file, sep = '\n', append = TRUE)
          cat('proc/;', file = file, sep = '\n', append = TRUE)
          
          # Also save a copy to out_TNT/pres_abs.tnt for the PAE-PCE script to use
          tnt_path <- file.path(getwd(), "out_TNT", "pres_abs.tnt")
          file.copy(file, tnt_path, overwrite = TRUE)
          
        }, error = function(e) {
          writeLines(paste("Error generating TNT matrix:", e$message), file)
        })
      }
    }
  )
  
  output$download_nexus <- downloadHandler(
    filename = function() { 
      method <- input$nexus_method
      if (is.null(method)) method <- "MST"
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      paste0("pres_abs_", method, "_", matrix_sel$basis, ".nex") 
    },
    content = function(file) {
      matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
      if (is.null(matrix_sel$matrix)) {
        writeLines("Error: No presence-absence matrix available", file)
      } else {
        tryCatch({
          method <- input$nexus_method
          if (is.null(method)) method <- "MST"
          
          # The range_nexus function writes to a file, so we need to capture it
          # First, we'll extract the logic from range_nexus to write directly to our file
          
          mat <- matrix_sel$matrix
          
          # Remove ROOT outgroup for NEXUS/BioGeoBEARS exports
          if ("ROOT" %in% rownames(mat)) {
            mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
          } else if (nrow(mat) > 0 && all(mat[nrow(mat), ] == 0)) {
            mat <- mat[-nrow(mat), , drop = FALSE]
          }
          
          # Keep empty rows for regular grids; drop them for irregular areas
          if (isTRUE(matrix_sel$is_irregular)) {
            row_sums <- rowSums(mat)
            mat <- mat[row_sums > 0, , drop = FALSE]
          }
          
          # Abbreviate area names for compatibility (A, B, C, ...)
          area_codes <- abbreviate_area_labels(row.names(mat), label_lookup = get_area_label_lookup())
          row.names(mat) <- area_codes$labels
          taxa <- row.names(mat)
          n.cara <- ncol(mat)
          n.taxa <- nrow(mat)
          
          # Replace dots and spaces with underscores in species names
          colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))
          colnames(mat) <- gsub(pattern = ' ', replacement = '_', x = colnames(mat))
          
          area_code_mapping_snapshot(list(
            mapping = area_codes$mapping,
            format = "NEXUS",
            timestamp = Sys.time()
          ))

          # Write NEXUS format
          cat('#NEXUS', file = file, sep = '\n')
          cat('begin data;', file = file, sep = '\n', append = TRUE)
          cat(paste0('dimensions ntax=', n.taxa, ' nchar=', n.cara, ';'), file = file, sep = '\n', append = TRUE)
          cat('format datatype=standard symbols="01" gap=-;', file = file, sep = '\n', append = TRUE)
          cat('', file = file, sep = '\n', append = TRUE)
          cat('CHARSTATELABELS', file = file, sep = '\n', append = TRUE)
          
          for(i in 1:n.cara){
            if(i == n.cara){
              cat(paste0(i, ' ', colnames(mat)[i], ';'), file = file, sep = '\n', append = TRUE)
            } else {
              cat(paste0(i, ' ', colnames(mat)[i], ','), file = file, sep = '\n', append = TRUE)
            }
          }
          
          cat('', file = file, sep = '\n', append = TRUE)
          cat('matrix', file = file, sep = '\n', append = TRUE)
          
          for(i in 1:n.taxa){
            cat(paste(taxa[i], paste(mat[i,], collapse = ''), sep = ' '), file = file, sep = '\n', append = TRUE)
          }
          
          cat(';', file = file, sep = '\n', append = TRUE)
          cat('end;', file = file, sep = '\n', append = TRUE)
          
        }, error = function(e) {
          writeLines(paste("Error generating NEXUS file:", e$message), file)
        })
      }
    }
  )
  
  output$download_tnt_script <- downloadHandler(
    filename = function() { "PAE_PCE_tntAUTO.run" },
    content = function(file) {
      # Create out_TNT directory if it doesn't exist
      dir.create(file.path(getwd(), "out_TNT"), showWarnings = FALSE)
      
      tryCatch({
        # First, we need to ensure the matrix file exists in the expected location
        # The generate_tnt_pae_pce function expects a matrix file to exist
        matrix_path <- "out_TNT/pres_abs.tnt"
        
        # If the matrix file doesn't exist but we have data, create it
        if (!file.exists(matrix_path) && !is.null(data_store$pres_abs)) {
          dir.create("out_TNT", showWarnings = FALSE)
          matrix_sel <- resolve_matrix_for_export(input$matrix_export_basis)
          
          # Create a temporary TNT matrix file
          mat <- matrix_sel$matrix %||% data_store$pres_abs
          
          # TNT bootstrap matrix requires ROOT outgroup row
          if (!("ROOT" %in% rownames(mat))) {
            mat <- rbind(mat, ROOT = rep(0, ncol(mat)))
          }
          
          # Keep empty rows for regular grids; drop them for irregular areas
          if (isTRUE(matrix_sel$is_irregular)) {
            row_sums <- rowSums(mat)
            mat <- mat[row_sums > 0, , drop = FALSE]
          }
          
          taxa <- row.names(mat)
          n.cara <- ncol(mat)
          n.taxa <- nrow(mat)
          
          colnames(mat) <- gsub(pattern = '\\.', replacement = '_', x = colnames(mat))
          colnames(mat) <- gsub(pattern = ' ', replacement = '_', x = colnames(mat))
          
          area_code_mapping_snapshot(list(
            mapping = area_codes$mapping,
            format = "TNT_script",
            timestamp = Sys.time()
          ))

          # Write TNT format to the expected path
          cat('nstates num 32;', file = matrix_path, sep = '\n')
          cat('xread', file = matrix_path, sep = '\n', append = TRUE)
          cat(paste(n.cara, n.taxa, sep = ' '), file = matrix_path, sep = '\n', append = TRUE)
          
          for(i in 1:n.taxa){
            cat(paste(taxa[i], paste(mat[i,], collapse = ''), sep = ' '), file = matrix_path, sep = '\n', append = TRUE)
          }
          
          cat(';', file = matrix_path, sep = '\n', append = TRUE)
          cat('cnames', file = matrix_path, sep = '\n', append = TRUE)
          
          for(i in 1:n.cara){
            cat(paste0('{', i-1, ' ', colnames(mat)[i], ';'), file = matrix_path, sep = '\n', append = TRUE)
          }
          
          cat(';', file = matrix_path, sep = '\n', append = TRUE)
          cat('proc/;', file = matrix_path, sep = '\n', append = TRUE)
        }
        
        # Now call the function with default parameters as requested
        script_path <- generate_tnt_pae_pce(
          matrix_file = matrix_path,
          output_file = "PAE_PCE_tntAUTO.run",
          max_iterations = 10,
          search_replicates = 4,
          search_method = "new_technology",
          output_dir = "out_TNT/"
        )
        
        # Read the generated script and write to the download file
        if (file.exists(script_path)) {
          script_content <- readLines(script_path)
          writeLines(script_content, file)
        } else {
          writeLines("Error: Script file was not generated properly.", file)
        }
      }, error = function(e) {
        writeLines(paste("Error generating TNT script:", e$message), file)
      })
    }
  )
  
  output$download_ndm <- downloadHandler(
    filename = function() { "data.xyd" },
    content = function(file) {
      if (is.null(data_store$occurrence)) {
        writeLines("Error: No occurrence data available", file)
      } else {
        tryCatch({
          # The toNDM function writes to a file, so we need to capture it
          # We'll use the occurrence data directly since it's the most reliable way
          
          data <- data_store$occurrence
          
          # Ensure data has correct columns
          if (!all(c("spp", "long", "lat") %in% colnames(data))) {
            if (all(c("species", "longitude", "latitude") %in% colnames(data))) {
              # Already correct
            } else {
              # Try to map columns
              colnames(data)[1:3] <- c("species", "longitude", "latitude")
            }
          } else {
            # Rename columns to match what toNDM expects internally
            colnames(data)[colnames(data) == "spp"] <- "species"
            colnames(data)[colnames(data) == "long"] <- "longitude"
            colnames(data)[colnames(data) == "lat"] <- "latitude"
          }
          
          # Remove NAs
          data <- data[complete.cases(data[, c("longitude", "latitude")]), ]
          
          # Get unique species
          species_list <- unique(data$species)
          species_list <- species_list[species_list != ""]
          
          # Write NDM format
          cat("xyddata", file = file, sep = "\n")
          
          for (sp in species_list) {
            sp_data <- data[data$species == sp, ]
            
            if (nrow(sp_data) > 0) {
              # Replace spaces with underscores in species name
              sp_name <- gsub(" ", "_", sp)
              
              cat(paste(sp_name, " ", nrow(sp_data)), file = file, sep = "\n", append = TRUE)
              
              for (i in 1:nrow(sp_data)) {
                cat(paste(sp_data$longitude[i], sp_data$latitude[i]), file = file, sep = "\n", append = TRUE)
              }
            }
          }
          
        }, error = function(e) {
          writeLines(paste("Error generating NDM file:", e$message), file)
        })
      }
    }
  )
  
  # Shapefiles Download Logic
  
  # Reactive value to store the list of files
  available_shapefiles <- reactiveVal(character(0))
  
  # Update file list when directory changes or refresh button is clicked
  observeEvent(c(input$shapefile_dir, input$refresh_shapefiles), {
    req(input$shapefile_dir)
    
    dir_path <- file.path(getwd(), input$shapefile_dir)
    
    if (dir.exists(dir_path)) {
      # Get all files in the directory
      all_files <- list.files(dir_path)
      
      # Filter for main downloadable outputs (.shp, .tif, .txt)
      # We don't list .dbf/.shx/.prj separately; they are bundled when downloading a .shp
      main_files <- all_files[grepl("\\.(shp|tif|txt)$", all_files, ignore.case = TRUE)]

      available_shapefiles(main_files)
    } else {
      available_shapefiles(character(0))
    }
  })
  
  # Generate UI for file downloads
  output$shapefile_download_ui <- renderUI({
    files <- available_shapefiles()
    
    if (length(files) == 0) {
      return(p("No files found in the selected directory. Run an extrapolation method first."))
    }
    
    # Create a list of download buttons
    download_buttons <- lapply(seq_along(files), function(i) {
      file_name <- files[i]
      # Create a unique ID for each button
      btn_id <- paste0("dl_shp_", i)
      
      div(
        style = "margin-bottom: 10px;",
        downloadButton(btn_id, label = paste("Download", file_name), class = "btn btn-sm btn-primary")
      )
    })
    
    do.call(tagList, download_buttons)
  })
  
  # Generate download handlers dynamically
  observe({
    files <- available_shapefiles()
    req(length(files) > 0)
    
    for (i in seq_along(files)) {
      local({
        my_i <- i
        my_file <- files[my_i]
        btn_id <- paste0("dl_shp_", my_i)
        
        output[[btn_id]] <- downloadHandler(
          filename = function() {
            # If it's a shapefile, we'll download a zip containing all components
            if (grepl("\\.shp$", my_file)) {
              paste0(tools::file_path_sans_ext(my_file), ".zip")
            } else {
              my_file
            }
          },
          content = function(file) {
            dir_path <- file.path(getwd(), input$shapefile_dir)
            source_path <- file.path(dir_path, my_file)
            
            if (grepl("\\.shp$", my_file)) {
              # For shapefiles, zip all related files (.shp, .shx, .dbf, .prj)
              base_name <- tools::file_path_sans_ext(my_file)
              related_files <- list.files(dir_path, pattern = paste0("^", base_name, "\\.(shp|shx|dbf|prj)$"), full.names = TRUE)
              
              # Create a temporary directory to store files before zipping
              temp_dir <- tempdir()
              zip_dir <- file.path(temp_dir, base_name)
              dir.create(zip_dir, showWarnings = FALSE)
              
              # Copy files to temp dir
              file.copy(related_files, zip_dir)
              
              # Zip them
              old_wd <- setwd(temp_dir)
              on.exit(setwd(old_wd))
              
              zip(file, files = base_name)
            } else {
              # For other files (.tif, .txt), just copy them
              file.copy(source_path, file)
            }
          }
        )
      })
    }
  })
  
  # ===== BIOGEOBEARS SETUP TAB =====
  bgb_model_results_table <- reactiveVal(NULL)
  bgb_model_results_raw <- reactiveVal(list())
  bgb_analysis_status <- reactiveVal("Waiting for BioGeoBEARS run...")
  bgb_last_inputs <- reactiveVal(NULL)
  bgb_last_plot_data <- reactiveVal(NULL)
  bgb_last_run_config <- reactiveVal(NULL)
  bgb_plot_error_text <- reactiveVal(NULL)
  bgb_pie_error_text <- reactiveVal(NULL)
  bgb_plot_diag_text <- reactiveVal(NULL)

  sanitize_taxon_names <- function(x) {
    x <- gsub("\\.", "_", x)
    x <- gsub("\\s+", "_", x)
    gsub("[^A-Za-z0-9_]", "_", x)
  }

  write_geog_from_current_matrix <- function(file_path) {
    mat <- data_store$pres_abs
    if (is.null(mat)) {
      stop("No matrix available from Step 3/5. Run an extrapolation first or upload .data file.")
    }

    if ("ROOT" %in% rownames(mat)) {
      mat <- mat[rownames(mat) != "ROOT", , drop = FALSE]
    }
    if (nrow(mat) > 0) {
      mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
    }
    if (nrow(mat) == 0 || ncol(mat) == 0) {
      stop("Current matrix has no informative data for BioGeoBEARS.")
    }

    area_codes <- abbreviate_area_labels(row.names(mat), label_lookup = get_area_label_lookup())
    row.names(mat) <- area_codes$labels
    area_code_mapping_snapshot(list(mapping = area_codes$mapping, format = "BioGeoBEARS_input", timestamp = Sys.time()))

    colnames(mat) <- sanitize_taxon_names(colnames(mat))
    mat[is.na(mat)] <- 0
    mat <- ifelse(mat > 0, 1, 0)

    # BioGeoBEARS::getranges_from_LagrangePHYLIP expects taxa in rows
    # (listed in header parentheses), and areas as line labels.
    # Internal workflows may produce either area x taxon or taxon x area.
    # Detect orientation from known species names and coerce to taxon x area.
    species_ref <- NULL
    if (!is.null(data_store$analysis_occurrence) && nrow(data_store$analysis_occurrence) > 0) {
      species_ref <- unique(as.character(data_store$analysis_occurrence$spp))
    } else if (!is.null(data_store$occurrence) && nrow(data_store$occurrence) > 0) {
      species_ref <- unique(as.character(data_store$occurrence$spp))
    }
    should_transpose <- FALSE
    if (!is.null(species_ref) && length(species_ref) > 0) {
      species_ref <- sanitize_taxon_names(species_ref)
      row_hits <- sum(rownames(mat) %in% species_ref)
      col_hits <- sum(colnames(mat) %in% species_ref)
      should_transpose <- (col_hits > row_hits)
    } else {
      tree_ref <- NULL
      if (!is.null(data_store$analysis_tree) && !is.null(data_store$analysis_tree$tip.label)) {
        tree_ref <- as.character(data_store$analysis_tree$tip.label)
      } else if (!is.null(data_store$tree) && !is.null(data_store$tree$tip.label)) {
        tree_ref <- as.character(data_store$tree$tip.label)
      }

      if (!is.null(tree_ref) && length(tree_ref) > 0) {
        tree_ref <- sanitize_taxon_names(tree_ref)
        row_hits <- sum(rownames(mat) %in% tree_ref)
        col_hits <- sum(colnames(mat) %in% tree_ref)
        should_transpose <- (col_hits > row_hits)
      } else {
        is_area_code <- function(x) grepl("^[A-Z]+$", x)
        row_area_score <- mean(is_area_code(rownames(mat)))
        col_area_score <- mean(is_area_code(colnames(mat)))
        should_transpose <- is.finite(row_area_score) && is.finite(col_area_score) && row_area_score > 0.8 && col_area_score < 0.5
      }
    }

    if (isTRUE(should_transpose)) {
      mat <- t(mat)
    }

    n_taxa <- nrow(mat)
    n_areas <- ncol(mat)
    taxa <- rownames(mat)
    areas <- colnames(mat)

    # Header: n_taxa n_areas (area names)
    # Body: one line per taxon with a 0/1 vector over areas.
    cat(paste(n_taxa, n_areas, paste0("(", paste(areas, collapse = " "), ")"), sep = " "), file = file_path, sep = "\n")
    for (i in seq_len(n_taxa)) {
      cat(paste(taxa[i], paste(mat[i, ], collapse = ""), sep = "\t"), file = file_path, sep = "\n", append = TRUE)
    }
  }

  write_geog_from_tipranges_df <- function(tip_df, file_path) {
    if (is.null(tip_df) || nrow(tip_df) == 0 || ncol(tip_df) == 0) {
      stop("Cannot write BioGeoBEARS geography file: empty taxa/area table.")
    }

    tip_df <- as.data.frame(tip_df, stringsAsFactors = FALSE)
    tip_df[is.na(tip_df)] <- 0
    tip_df <- ifelse(as.matrix(tip_df) > 0, 1, 0)
    rownames(tip_df) <- sanitize_taxon_names(rownames(tip_df))
    colnames(tip_df) <- sanitize_taxon_names(colnames(tip_df))

    n.taxa <- nrow(tip_df)
    n.cara <- ncol(tip_df)
    taxa <- rownames(tip_df)
    areas <- colnames(tip_df)

    cat(paste(n.taxa, n.cara, paste0("(", paste(areas, collapse = " "), ")"), sep = " "), file = file_path, sep = "\n")
    for (i in seq_len(n.taxa)) {
      cat(paste(taxa[i], paste(tip_df[i, ], collapse = ""), sep = "\t"), file = file_path, sep = "\n", append = TRUE)
    }
  }

  harmonize_tree_and_geog_files <- function(tree_file, geog_file) {
    tr <- read_bgb_tree_safe(tree_file)
    tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geog_file)
    geog_df <- as.data.frame(tipranges@df)

    tree_tips <- sanitize_taxon_names(as.character(tr$tip.label))
    geog_taxa <- sanitize_taxon_names(rownames(geog_df))
    rownames(geog_df) <- geog_taxa

    keep_taxa <- intersect(tree_tips, geog_taxa)
    missing_in_geog <- setdiff(tree_tips, geog_taxa)
    missing_in_tree <- setdiff(geog_taxa, tree_tips)

    if (length(keep_taxa) == 0) {
      tree_preview <- paste(utils::head(tree_tips, 5), collapse = ", ")
      geog_preview <- paste(utils::head(geog_taxa, 5), collapse = ", ")
      stop(
        paste0(
          "Tree and geography file share zero taxa after normalization. ",
          "Tree tips (n=", length(tree_tips), ") sample: [", tree_preview, "]; ",
          "Geography taxa (n=", length(geog_taxa), ") sample: [", geog_preview, "]"
        )
      )
    }

    changed <- (length(missing_in_geog) > 0 || length(missing_in_tree) > 0)
    if (changed) {
      if (length(missing_in_geog) >= length(tree_tips)) {
        stop("All tree tips are absent from the geography file after normalization.")
      }

      drop_from_tree <- setdiff(tree_tips, keep_taxa)
      if (length(drop_from_tree) > 0) {
        tr <- ape::drop.tip(tr, drop_from_tree)
      }
      tr$tip.label <- sanitize_taxon_names(tr$tip.label)
      ape::write.tree(tr, file = tree_file)

      keep_order <- tr$tip.label
      geog_df <- geog_df[keep_order, , drop = FALSE]
      write_geog_from_tipranges_df(geog_df, geog_file)
    }

    list(
      changed = changed,
      kept_n = length(keep_taxa),
      missing_in_geog = missing_in_geog,
      missing_in_tree = missing_in_tree
    )
  }

  get_bgb_areanames <- function(geogfn) {
    tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
    areas <- names(tipranges@df)
    areas <- as.character(areas)
    areas <- areas[!is.na(areas) & nzchar(areas)]
    if (length(areas) == 0) {
      stop("Could not extract area names from BioGeoBEARS geography file.")
    }
    areas
  }

  coerce_square_matrix_for_bgb <- function(mat, label, areanames, require_positive_offdiag = FALSE) {
    if (is.null(mat)) {
      stop(paste0(label, " matrix is NULL."))
    }

    m <- as.matrix(mat)
    suppressWarnings(storage.mode(m) <- "numeric")

    if (nrow(m) == 0 || ncol(m) == 0) {
      stop(paste0(label, " matrix is empty."))
    }
    if (nrow(m) != ncol(m)) {
      stop(paste0(label, " matrix must be square (nrow == ncol)."))
    }
    if (any(!is.finite(m), na.rm = TRUE)) {
      stop(paste0(label, " matrix contains non-finite values."))
    }

    rn <- rownames(m)
    cn <- colnames(m)
    if (is.null(rn) || any(!nzchar(rn))) {
      rn <- cn
    }
    if (is.null(cn) || any(!nzchar(cn))) {
      cn <- rn
    }
    if (is.null(rn) || is.null(cn) || any(!nzchar(rn)) || any(!nzchar(cn))) {
      stop(paste0(label, " matrix must have row and column names matching BioGeoBEARS area codes."))
    }
    rownames(m) <- as.character(rn)
    colnames(m) <- as.character(cn)

    if (!setequal(rownames(m), areanames) || !setequal(colnames(m), areanames)) {
      stop(
        paste0(
          label,
          " matrix area names do not match geography area names.\n",
          "Expected: ", paste(areanames, collapse = ", "), "\n",
          "Rows: ", paste(rownames(m), collapse = ", "), "\n",
          "Cols: ", paste(colnames(m), collapse = ", ")
        )
      )
    }

    m <- m[areanames, areanames, drop = FALSE]

    if (require_positive_offdiag) {
      diag(m) <- NA_real_
      if (any(m <= 0, na.rm = TRUE)) {
        stop(paste0(label, " matrix must have strictly positive off-diagonal values (>0)."))
      }
      diag(m) <- 0
    }

    m
  }

  write_bgb_matrix_list_file <- function(mats_list, out_file) {
    BioGeoBEARS::write_distances_to_fn(new_distmats_list = mats_list, outfn = out_file)
    out_file
  }

  file_looks_like_bgb_matrix_file <- function(file_info) {
    if (is.null(file_info)) {
      return(FALSE)
    }
    lines <- readLines(file_info$datapath, warn = FALSE)
    any(toupper(trimws(lines)) == "END")
  }

  count_bgb_matrix_blocks <- function(file_info) {
    if (is.null(file_info)) {
      return(NA_integer_)
    }
    lines <- readLines(file_info$datapath, warn = FALSE)
    lines <- gsub("\\r", "", lines)
    end_idx <- which(toupper(trimws(lines)) == "END")
    if (length(end_idx) > 0) {
      lines <- lines[seq_len(end_idx[1] - 1)]
    }
    if (length(lines) == 0) {
      return(0L)
    }
    is_blank <- trimws(lines) == ""
    blocks <- sum(!is_blank & c(TRUE, head(is_blank, -1)))
    as.integer(blocks)
  }

  summarize_upload_file <- function(file_info) {
    if (is.null(file_info)) {
      return(list(name = NA_character_, has_end = FALSE, blocks = NA_integer_))
    }
    list(
      name = file_info$name,
      has_end = file_looks_like_bgb_matrix_file(file_info),
      blocks = count_bgb_matrix_blocks(file_info)
    )
  }

  copy_uploaded_file_to_temp <- function(file_info, fileext = ".txt") {
    tmp <- tempfile(fileext = fileext)
    file.copy(file_info$datapath, tmp, overwrite = TRUE)
    tmp
  }

  prepare_bgb_constraints <- function(geogfn) {
    areanames <- get_bgb_areanames(geogfn)
    n_time <- 1L
    notes <- character(0)
    out <- list()

    if (isTRUE(input$use_time_slices)) {
      tp <- data_store$time_periods
      if (is.null(tp) || length(tp) == 0) {
        stop("'Use time slices' is enabled, but no time periods were loaded.")
      }
      n_time <- length(tp)
      out$timeperiods <- tp
      times_file <- tempfile(fileext = ".txt")
      writeLines(as.character(tp), con = times_file)
      out$timesfn <- times_file
      notes <- c(notes, paste0("time slices enabled (", n_time, " periods)"))
    }

    if (isTRUE(input$use_distance_matrix)) {
      if (!is.null(input$distance_matrix_file) && file_looks_like_bgb_matrix_file(input$distance_matrix_file)) {
        out$distsfn <- copy_uploaded_file_to_temp(input$distance_matrix_file)
        notes <- c(notes, "distance matrix applied from uploaded BioGeoBEARS file (x)")
      } else {
      if (is.null(data_store$distance_matrix) && !is.null(input$distance_matrix_file)) {
        data_store$distance_matrix <- load_named_square_matrix(
          file_info = input$distance_matrix_file,
          matrix_label = "Distance matrix",
          require_binary = FALSE,
          require_positive_offdiag = TRUE
        )
      }
      if (is.null(data_store$distance_matrix)) {
        has_file <- !is.null(input$distance_matrix_file)
        stop(
          paste0(
            "'Use distance matrix' is enabled, but no distance matrix was loaded. ",
            "File selected in UI: ", has_file, ". ",
            "Please reselect the file and click 'Load Matrix'."
          )
        )
      }
      dmat <- coerce_square_matrix_for_bgb(
        mat = data_store$distance_matrix,
        label = "Distance",
        areanames = areanames,
        require_positive_offdiag = TRUE
      )
      dlist <- rep(list(dmat), n_time)
      d_file <- tempfile(fileext = ".txt")
      write_bgb_matrix_list_file(dlist, d_file)
      out$distsfn <- d_file
      notes <- c(notes, "distance matrix applied (x)")
      }
    }

    if (!is.null(data_store$env_distance_matrix)) {
      envmat <- coerce_square_matrix_for_bgb(
        mat = data_store$env_distance_matrix,
        label = "Environmental distance",
        areanames = areanames,
        require_positive_offdiag = TRUE
      )
      elist <- rep(list(envmat), n_time)
      e_file <- tempfile(fileext = ".txt")
      write_bgb_matrix_list_file(elist, e_file)
      out$envdistsfn <- e_file
      notes <- c(notes, "environmental distance matrix applied (n)")
    } else if (!is.null(input$env_distance_matrix_file) && file_looks_like_bgb_matrix_file(input$env_distance_matrix_file)) {
      out$envdistsfn <- copy_uploaded_file_to_temp(input$env_distance_matrix_file)
      notes <- c(notes, "environmental distance matrix applied from uploaded BioGeoBEARS file (n)")
    }

    if (isTRUE(input$use_dispersal_multiplier)) {
      if (!is.null(input$dispersal_multipliers_file) && file_looks_like_bgb_matrix_file(input$dispersal_multipliers_file)) {
        out$dispersal_multipliers_fn <- copy_uploaded_file_to_temp(input$dispersal_multipliers_file)
        notes <- c(notes, "dispersal multipliers applied from uploaded BioGeoBEARS file (w)")
      } else {
      if (is.null(data_store$dispersal_multipliers) && !is.null(input$dispersal_multipliers_file)) {
        data_store$dispersal_multipliers <- load_named_square_matrix(
          file_info = input$dispersal_multipliers_file,
          matrix_label = "Dispersal multipliers matrix",
          require_binary = TRUE,
          require_positive_offdiag = FALSE
        )
      }
      if (is.null(data_store$dispersal_multipliers)) {
        stop("'Use dispersal multiplier' is enabled, but no dispersal multipliers matrix was loaded.")
      }
      dmm <- coerce_square_matrix_for_bgb(
        mat = data_store$dispersal_multipliers,
        label = "Dispersal multipliers",
        areanames = areanames,
        require_positive_offdiag = FALSE
      )
      dmlist <- rep(list(dmm), n_time)
      dm_file <- tempfile(fileext = ".txt")
      write_bgb_matrix_list_file(dmlist, dm_file)
      out$dispersal_multipliers_fn <- dm_file
      notes <- c(notes, "dispersal multipliers applied (w)")
      }
    }

    if (!is.null(data_store$areas_allowed_matrix)) {
      aam <- coerce_square_matrix_for_bgb(
        mat = data_store$areas_allowed_matrix,
        label = "Areas-allowed",
        areanames = areanames,
        require_positive_offdiag = FALSE
      )
      aalist <- rep(list(aam), n_time)
      aa_file <- tempfile(fileext = ".txt")
      write_bgb_matrix_list_file(aalist, aa_file)
      out$areas_allowed_fn <- aa_file
      notes <- c(notes, "areas-allowed matrix applied")
    } else if (!is.null(input$areas_allowed_file) && file_looks_like_bgb_matrix_file(input$areas_allowed_file)) {
      out$areas_allowed_fn <- copy_uploaded_file_to_temp(input$areas_allowed_file)
      notes <- c(notes, "areas-allowed matrix applied from uploaded BioGeoBEARS file")
    }

    if (!is.null(data_store$areas_adjacency_matrix)) {
      ajm <- coerce_square_matrix_for_bgb(
        mat = data_store$areas_adjacency_matrix,
        label = "Areas-adjacency",
        areanames = areanames,
        require_positive_offdiag = FALSE
      )
      ajlist <- rep(list(ajm), n_time)
      aj_file <- tempfile(fileext = ".txt")
      write_bgb_matrix_list_file(ajlist, aj_file)
      out$areas_adjacency_fn <- aj_file
      notes <- c(notes, "areas-adjacency matrix applied")
    } else if (!is.null(input$areas_adjacency_file) && file_looks_like_bgb_matrix_file(input$areas_adjacency_file)) {
      out$areas_adjacency_fn <- copy_uploaded_file_to_temp(input$areas_adjacency_file)
      notes <- c(notes, "areas-adjacency matrix applied from uploaded BioGeoBEARS file")
    }

    out$notes <- notes
    out
  }

  make_bgb_run_object <- function(trfn, geogfn, constraints = NULL) {
    run_obj <- BioGeoBEARS::define_BioGeoBEARS_run()
    run_obj$trfn <- trfn
    run_obj$geogfn <- geogfn
    run_obj$max_range_size <- input$bgb_max_range_size
    run_obj$min_branchlength <- input$bgb_min_branchlength
    run_obj$include_null_range <- isTRUE(input$bgb_include_null_range)
    run_obj$use_optimx <- identical(input$bgb_optimizer, "optimx")
    run_obj$num_cores_to_use <- input$bgb_num_cores
    run_obj$return_condlikes_table <- FALSE
    run_obj$calc_ancprobs <- TRUE

    # Compatibility guard for BioGeoBEARS builds where envdistsfn
    # is expected by check_BioGeoBEARS_run() but not pre-initialized.
    if (is.null(run_obj$envdistsfn)) {
      run_obj$envdistsfn <- NA
    }

    if (!is.null(constraints)) {
      if (!is.null(constraints$timesfn)) run_obj$timesfn <- constraints$timesfn
      if (!is.null(constraints$timeperiods)) run_obj$timeperiods <- constraints$timeperiods
      if (!is.null(constraints$distsfn)) run_obj$distsfn <- constraints$distsfn
      if (!is.null(constraints$envdistsfn)) run_obj$envdistsfn <- constraints$envdistsfn
      if (!is.null(constraints$dispersal_multipliers_fn)) run_obj$dispersal_multipliers_fn <- constraints$dispersal_multipliers_fn
      if (!is.null(constraints$areas_allowed_fn)) run_obj$areas_allowed_fn <- constraints$areas_allowed_fn
      if (!is.null(constraints$areas_adjacency_fn)) run_obj$areas_adjacency_fn <- constraints$areas_adjacency_fn
    }

    run_obj
  }

  detect_max_tipsize_from_geog <- function(geogfn) {
    tipranges <- BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
    tip_df <- as.data.frame(tipranges@df)
    if (nrow(tip_df) == 0 || ncol(tip_df) == 0) {
      return(0L)
    }
    rs <- rowSums(tip_df > 0, na.rm = TRUE)
    mx <- suppressWarnings(max(rs, na.rm = TRUE))
    if (!is.finite(mx)) {
      return(0L)
    }
    as.integer(mx)
  }

  extract_optimizer_diagnostics <- function(model_fit) {
    safe_get <- function(expr) {
      tryCatch(eval.parent(substitute(expr)), error = function(e) NULL)
    }

    conv_candidates <- list(
      safe_get(model_fit$optim_result$convergence),
      safe_get(model_fit$optim_result$opt_result$convergence),
      safe_get(model_fit$optim_result$details$convergence),
      safe_get(model_fit$outputs$optim_result$convergence),
      safe_get(model_fit$convergence)
    )

    msg_candidates <- list(
      safe_get(model_fit$optim_result$message),
      safe_get(model_fit$optim_result$opt_result$message),
      safe_get(model_fit$optim_result$details$message),
      safe_get(model_fit$outputs$optim_result$message),
      safe_get(model_fit$message)
    )

    conv <- NA_integer_
    for (candidate in conv_candidates) {
      val <- suppressWarnings(as.integer(candidate))[1]
      if (!is.na(val)) {
        conv <- val
        break
      }
    }

    msg <- NA_character_
    for (candidate in msg_candidates) {
      txt <- as.character(candidate)[1]
      if (!is.na(txt) && nzchar(txt)) {
        msg <- txt
        break
      }
    }

    list(convergence = conv, message = msg)
  }

  read_bgb_tree_safe <- function(tree_file) {
    if (is.null(tree_file) || !file.exists(tree_file)) {
      stop("Tree file used by BioGeoBEARS is not available in this session.")
    }

    raw_tree <- tryCatch(
      ape::read.tree(tree_file),
      error = function(e_tree) {
        tryCatch(
          ape::read.nexus(tree_file),
          error = function(e_nexus) {
            stop(
              paste0(
                "Could not read BioGeoBEARS tree file for plotting as Newick or Nexus. ",
                "Newick error: ", conditionMessage(e_tree), "; ",
                "Nexus error: ", conditionMessage(e_nexus)
              )
            )
          }
        )
      }
    )

    if (inherits(raw_tree, "multiPhylo")) {
      if (length(raw_tree) == 0) {
        stop("BioGeoBEARS tree file has no trees available for plotting.")
      }
      return(raw_tree[[1]])
    }

    raw_tree
  }

  get_bgb_results_params_table <- function(res_obj) {
    candidates <- list(
      tryCatch(res_obj$output@params_table, error = function(e) NULL),
      tryCatch(res_obj$outputs@params_table, error = function(e) NULL),
      tryCatch(res_obj$inputs$BioGeoBEARS_model_object@params_table, error = function(e) NULL)
    )

    for (candidate in candidates) {
      if (!is.null(candidate) && !is.null(dim(candidate)) && nrow(candidate) > 0) {
        return(candidate)
      }
    }

    NULL
  }

  get_bgb_plot_inputs <- function(res_obj, inputs) {
    merged <- inputs

    if (is.null(merged)) {
      merged <- list()
    }

    result_inputs <- tryCatch(res_obj$inputs, error = function(e) NULL)
    if (is.null(result_inputs)) {
      return(merged)
    }

    if (is.null(merged$trfn) && !is.null(result_inputs$trfn)) {
      merged$trfn <- result_inputs$trfn
    }
    if (is.null(merged$geogfn) && !is.null(result_inputs$geogfn)) {
      merged$geogfn <- result_inputs$geogfn
    }
    if (is.null(merged$include_null_range) && !is.null(result_inputs$include_null_range)) {
      merged$include_null_range <- result_inputs$include_null_range
    }

    merged
  }

  build_bgb_plot_data_cache <- function(inputs) {
    out <- list(tr = NULL, tipranges = NULL)
    if (is.null(inputs)) {
      return(out)
    }

    if (!is.null(inputs$trfn) && file.exists(inputs$trfn)) {
      out$tr <- tryCatch(read_bgb_tree_safe(inputs$trfn), error = function(e) NULL)
    }
    if (!is.null(inputs$geogfn) && file.exists(inputs$geogfn)) {
      out$tipranges <- tryCatch(
        BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = inputs$geogfn),
        error = function(e) NULL
      )
    }
    out
  }

  extract_plot_data_from_run_obj <- function(run_obj) {
    tr_obj <- tryCatch(run_obj$tr, error = function(e) NULL)
    if (inherits(tr_obj, "multiPhylo") && length(tr_obj) > 0) {
      tr_obj <- tr_obj[[1]]
    }

    tipranges_obj <- tryCatch(run_obj$tipranges, error = function(e) NULL)

    list(
      tr = tr_obj,
      tipranges = tipranges_obj
    )
  }

  plot_bgb_results_safe <- function(res_obj, model_name, inputs, plot_data = NULL, plotwhat = c("text", "pie"), title_suffix = NULL, label_offset_override = NULL, tipcex_override = NULL, statecex_override = NULL) {
    plotwhat <- match.arg(plotwhat)
    inputs <- get_bgb_plot_inputs(res_obj = res_obj, inputs = inputs)

    tr <- tryCatch(plot_data$tr, error = function(e) NULL)
    tipranges <- tryCatch(plot_data$tipranges, error = function(e) NULL)

    # Prefer re-reading tree/geography from files at plot time (as in the
    # standalone app workflow), which is the most reliable path for Shiny devices.
    read_plot_inputs_from_files <- function(inp) {
      out <- list(tr = NULL, tipranges = NULL)
      tr_path <- tryCatch(inp$trfn, error = function(e) NULL)
      geog_path <- tryCatch(inp$geogfn, error = function(e) NULL)
      if (!is.null(tr_path) && file.exists(tr_path)) {
        out$tr <- tryCatch(read_bgb_tree_safe(tr_path), error = function(e) NULL)
      }
      if (!is.null(geog_path) && file.exists(geog_path)) {
        out$tipranges <- tryCatch(BioGeoBEARS::getranges_from_LagrangePHYLIP(lgdata_fn = geog_path), error = function(e) NULL)
      }
      out
    }

    file_inputs <- read_plot_inputs_from_files(inputs)
    if (is.null(file_inputs$tr) || is.null(file_inputs$tipranges)) {
      res_inputs <- tryCatch(res_obj$inputs, error = function(e) NULL)
      if (!is.null(res_inputs)) {
        res_file_inputs <- read_plot_inputs_from_files(res_inputs)
        if (is.null(file_inputs$tr)) file_inputs$tr <- res_file_inputs$tr
        if (is.null(file_inputs$tipranges)) file_inputs$tipranges <- res_file_inputs$tipranges
      }
    }
    if (!is.null(file_inputs$tr)) tr <- file_inputs$tr
    if (!is.null(file_inputs$tipranges)) tipranges <- file_inputs$tipranges

    if (is.null(tr) || is.null(tipranges)) {
      cache_fallback <- build_bgb_plot_data_cache(inputs)
      if (is.null(tr)) {
        tr <- cache_fallback$tr
      }
      if (is.null(tipranges)) {
        tipranges <- cache_fallback$tipranges
      }
    }

    if (is.null(tr) && is.null(tipranges) && is.null(tryCatch(res_obj$inputs, error = function(e) NULL))) {
      stop("BioGeoBEARS plotting could not recover tree/geography inputs from session files or result object.")
    }

    params_table <- get_bgb_results_params_table(res_obj)

    addl_params <- list("j")
    if (!is.null(params_table) && "j" %in% rownames(params_table)) {
      j_type <- tryCatch(as.character(params_table["j", "type"]), error = function(e) NA_character_)
      j_est <- tryCatch(as.numeric(params_table["j", "est"]), error = function(e) NA_real_)

      if (isTRUE(identical(j_type, "fixed")) && (is.na(j_est) || isTRUE(all.equal(j_est, 0)))) {
        addl_params <- list()
      }
    }

    title_txt <- paste0("BioGeoBEARS ", model_name)
    if (!is.null(title_suffix) && nzchar(title_suffix)) {
      title_txt <- paste0(title_txt, " - ", title_suffix)
    }

    label_offset_default <- if (identical(plotwhat, "text")) 0.35 else 0.2
    label_offset_value <- suppressWarnings(as.numeric(label_offset_override)[1])
    if (is.na(label_offset_value) || !is.finite(label_offset_value)) {
      label_offset_value <- label_offset_default
    }

    tipcex_value <- suppressWarnings(as.numeric(tipcex_override)[1])
    if (is.na(tipcex_value) || !is.finite(tipcex_value) || tipcex_value <= 0) {
      tipcex_value <- 0.6
    }

    statecex_value <- suppressWarnings(as.numeric(statecex_override)[1])
    if (is.na(statecex_value) || !is.finite(statecex_value) || statecex_value <= 0) {
      statecex_value <- 0.7
    }

    args <- list(
      results_object = res_obj,
      analysis_titletxt = title_txt,
      addl_params = addl_params,
      plotwhat = plotwhat,
      label.offset = label_offset_value,
      tipcex = tipcex_value,
      statecex = statecex_value,
      splitcex = 0.6,
      titlecex = 0.8,
      plotsplits = FALSE,
      include_null_range = isTRUE(inputs$include_null_range),
      tr = tr,
      tipranges = tipranges
    )

    if (is.null(args$tr)) {
      args$tr <- NULL
    }
    if (is.null(args$tipranges)) {
      args$tipranges <- NULL
    }

    scriptdir <- system.file("extdata/a_scripts", package = "BioGeoBEARS")
    if (nzchar(scriptdir) && dir.exists(scriptdir)) {
      args$cornercoords_loc <- BioGeoBEARS::np(scriptdir)
    }

    call_plot_bgb <- function(plot_args, strict_text_warnings = FALSE) {
      tryCatch(
        withCallingHandlers(
          do.call(BioGeoBEARS::plot_BioGeoBEARS_results, plot_args),
          warning = function(w) {
            msg <- conditionMessage(w)
            if (isTRUE(strict_text_warnings) && grepl("mean\\.default|not numeric|não é numérico|retornando NA|returning NA", msg, ignore.case = TRUE)) {
              stop(simpleError(msg))
            }
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) e
      )
    }

    strict_warn <- identical(plotwhat, "text")
    out <- call_plot_bgb(args, strict_text_warnings = strict_warn)

    if (inherits(out, "error") && length(args$addl_params) > 0) {
      args$addl_params <- list()
      out <- call_plot_bgb(args, strict_text_warnings = strict_warn)
    }
    if (inherits(out, "error") && !is.null(args$cornercoords_loc)) {
      args$cornercoords_loc <- NULL
      out <- call_plot_bgb(args, strict_text_warnings = strict_warn)
    }
    res_include_null <- tryCatch(res_obj$inputs$include_null_range, error = function(e) NULL)
    if (inherits(out, "error") && !is.null(res_include_null) && !identical(args$include_null_range, isTRUE(res_include_null))) {
      args$include_null_range <- isTRUE(res_include_null)
      out <- call_plot_bgb(args, strict_text_warnings = strict_warn)
    }
    if (inherits(out, "error")) {
      fallback_args <- args
      fallback_args$tr <- NULL
      fallback_args$tipranges <- NULL
      out <- call_plot_bgb(fallback_args, strict_text_warnings = strict_warn)
    }

    if (inherits(out, "error")) {
      stop(conditionMessage(out))
    }

    if (identical(plotwhat, "pie") && !is.null(tr) && !is.null(tipranges)) {
      try({
        bgb_ns <- asNamespace("BioGeoBEARS")
        tips <- seq_along(tr$tip.label)
        areas <- bgb_ns$getareas_from_tipranges_object(tipranges)

        max_range_size <- tryCatch(as.numeric(res_obj$inputs$max_range_size)[1], error = function(e) NA_real_)
        if (is.na(max_range_size)) {
          max_range_size <- length(areas)
        }

        states_list_0based_index <- tryCatch(res_obj$inputs$states_list, error = function(e) NULL)
        if (is.null(states_list_0based_index)) {
          states_list_0based_index <- bgb_ns$rcpp_areas_list_to_states_list(
            areas,
            maxareas = max_range_size,
            include_null_range = isTRUE(res_obj$inputs$include_null_range)
          )
        }

        ranges_list <- tryCatch(
          bgb_ns$states_list_0based_to_ranges_txt_list(
            state_indices_0based = states_list_0based_index,
            areanames = areas
          ),
          error = function(e) NULL
        )
        if (is.null(ranges_list)) {
          ranges_list <- bgb_ns$areas_list_to_states_list_new(
            areas,
            maxareas = max_range_size,
            split_ABC = FALSE,
            include_null_range = isTRUE(res_obj$inputs$include_null_range)
          )
        }

        statenames <- unlist(ranges_list)
        relprobs_matrix <- res_obj$ML_marginal_prob_each_state_at_branch_top_AT_node
        MLstates <- bgb_ns$get_ML_states_from_relprobs(
          relprobs_matrix,
          statenames,
          returnwhat = "states",
          if_ties = "takefirst"
        )

        colors_matrix <- bgb_ns$get_colors_for_numareas(length(areas))
        colors_list_for_states <- bgb_ns$mix_colors_for_states(
          colors_matrix,
          states_list_0based_index,
          plot_null_range = isTRUE(res_obj$inputs$include_null_range)
        )
        cols_byNode <- bgb_ns$rangestxt_to_colors(ranges_list, colors_list_for_states, MLstates)

        ape::tiplabels(
          text = MLstates[tips],
          tip = tips,
          bg = cols_byNode[tips],
          col = "black",
          frame = "rect",
          cex = statecex_value,
          adj = c(0.5, 0.5)
        )
      }, silent = TRUE)
    }

    invisible(out)
  }

  draw_bgb_plot_with_replay <- function(plot_fun) {
    tmp_png <- tempfile(fileext = ".png")
    rec <- NULL

    tryCatch({
      grDevices::png(filename = tmp_png, width = 2200, height = 1400, res = 180)
      grDevices::dev.control(displaylist = "enable")
      plot_fun()
      rec <- grDevices::recordPlot()
      grDevices::dev.off()
    }, error = function(e) {
      if (!is.null(grDevices::dev.list())) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      stop(e)
    }, finally = {
      if (file.exists(tmp_png)) {
        unlink(tmp_png)
      }
    })

    if (is.null(rec)) {
      stop("BioGeoBEARS plot was generated without a replayable graphics record.")
    }

    grDevices::replayPlot(rec)
    invisible(TRUE)
  }

  should_use_bgb_replay <- function(inputs) {
    method <- tryCatch(inputs$extrap_method, error = function(e) NULL)
    isTRUE(method %in% c("mst", "buffer", "convex_hull"))
  }

  configure_bgb_model <- function(run_obj, model_code) {
    p <- run_obj$BioGeoBEARS_model_object@params_table

    set_fixed <- function(param, value) {
      p[param, "type"] <<- "fixed"
      p[param, "init"] <<- value
      p[param, "est"] <<- value
    }
    set_free <- function(param, value) {
      p[param, "type"] <<- "free"
      p[param, "init"] <<- value
      p[param, "est"] <<- value
    }

    set_fixed("j", 0)
    set_fixed("x", 0)
    set_fixed("n", 0)
    set_fixed("w", 1)

    # Numerical-stability bounds (following common BioGeoBEARS guidance,
    # especially useful for +J models like BAYAREALIKE+J on some platforms).
    if ("d" %in% rownames(p)) {
      p["d", "min"] <- 1e-7
      p["d", "max"] <- 4.9999999
    }
    if ("e" %in% rownames(p)) {
      p["e", "min"] <- 1e-7
      p["e", "max"] <- 4.9999999
    }
    if ("j" %in% rownames(p)) {
      p["j", "min"] <- 1e-5
      p["j", "max"] <- 0.99999
    }

    if (isTRUE(input$use_distance_matrix) && !is.null(data_store$distance_matrix)) {
      set_free("x", 0.1)
    }
    if (!is.null(data_store$env_distance_matrix)) {
      set_free("n", 1)
    }
    if (isTRUE(input$use_dispersal_multiplier) && !is.null(data_store$dispersal_multipliers)) {
      set_free("w", 1)
    }

    if (model_code %in% c("DECJ", "DIVALIKEJ", "BAYAREALIKEJ")) {
      set_free("j", 0.01)
    }

    if (model_code %in% c("DIVALIKE", "DIVALIKEJ")) {
      set_fixed("s", 0)
      p["ysv", "type"] <- "2-j"
      p["ys", "type"] <- "ysv*1/2"
      p["y", "type"] <- "ysv*1/2"
      p["v", "type"] <- "ysv*1/2"
      set_fixed("mx01v", 0.5)
    }

    if (model_code %in% c("BAYAREALIKE", "BAYAREALIKEJ")) {
      set_fixed("s", 0)
      set_fixed("v", 0)
      p["ysv", "type"] <- "1-j"
      p["ys", "type"] <- "ysv*1/1"
      p["y", "type"] <- "1-j"
      set_fixed("mx01y", 0.9999)
    }

    run_obj$BioGeoBEARS_model_object@params_table <- p
    BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = run_obj)
  }

  get_nested_base_model <- function(model_code) {
    switch(
      model_code,
      DECJ = "DEC",
      DIVALIKEJ = "DIVALIKE",
      BAYAREALIKEJ = "BAYAREALIKE",
      NULL
    )
  }

  get_model_param_estimate <- function(model_result, param_name) {
    if (is.null(model_result)) {
      return(NA_real_)
    }
    val <- tryCatch(as.numeric(model_result$output@params_table[param_name, "est"]), error = function(e) NA_real_)
    if (length(val) == 0 || !is.finite(val)) {
      return(NA_real_)
    }
    val
  }

  apply_nested_warm_start <- function(run_obj, model_code, nested_result, include_extra = FALSE) {
    if (is.null(nested_result)) {
      return(list(run_obj = run_obj, applied = FALSE, note = NULL))
    }

    p <- run_obj$BioGeoBEARS_model_object@params_table
    applied_params <- character(0)

    param_names <- c("d", "e")
    if (isTRUE(include_extra)) {
      param_names <- c(param_names, "x", "w", "n")
    }

    for (param_name in param_names) {
      if (!(param_name %in% rownames(p))) {
        next
      }
      if (!identical(p[param_name, "type"], "free")) {
        next
      }
      start_val <- get_model_param_estimate(nested_result, param_name)
      if (!is.finite(start_val)) {
        next
      }
      pmin <- suppressWarnings(as.numeric(p[param_name, "min"]))
      pmax <- suppressWarnings(as.numeric(p[param_name, "max"]))
      if (is.finite(pmin) && start_val <= pmin) {
        start_val <- pmin + 1e-8
      }
      if (is.finite(pmax) && start_val >= pmax) {
        start_val <- pmax - 1e-8
      }
      if (is.finite(start_val)) {
        p[param_name, "init"] <- start_val
        p[param_name, "est"] <- start_val
        applied_params <- c(applied_params, paste0(param_name, "=", signif(start_val, 5)))
      }
    }

    if (model_code %in% c("DECJ", "DIVALIKEJ", "BAYAREALIKEJ") && ("j" %in% rownames(p))) {
      j_start <- 1e-4
      jmin <- suppressWarnings(as.numeric(p["j", "min"]))
      jmax <- suppressWarnings(as.numeric(p["j", "max"]))
      if (is.finite(jmin) && j_start <= jmin) {
        j_start <- jmin + 1e-8
      }
      if (is.finite(jmax) && j_start >= jmax) {
        j_start <- jmax - 1e-8
      }
      if (is.finite(j_start)) {
        p["j", "init"] <- j_start
        p["j", "est"] <- j_start
      }
    }

    run_obj$BioGeoBEARS_model_object@params_table <- p
    run_obj <- BioGeoBEARS::fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = run_obj)

    if (length(applied_params) == 0) {
      return(list(run_obj = run_obj, applied = FALSE, note = NULL))
    }

    list(
      run_obj = run_obj,
      applied = TRUE,
      note = paste0("warm-start from nested model (", paste(applied_params, collapse = ", "), ")")
    )
  }

  output$bgb_input_status <- renderPrint({
    cat("BioGeoBEARS input status\n")
    cat("- Geographic source:", ifelse(input$bgb_geog_source %||% "current" == "current", "Current matrix (Step 3/5)", "Uploaded .data file"), "\n")
    cat("- Tree source:", ifelse(input$bgb_tree_source %||% "step1" == "step1", "Tree from Step 1", "Uploaded Newick file"), "\n")
    if ((input$bgb_geog_source %||% "current") == "current" && !is.null(data_store$pres_abs)) {
      mat_status <- data_store$pres_abs
      if ("ROOT" %in% rownames(mat_status)) {
        mat_status <- mat_status[rownames(mat_status) != "ROOT", , drop = FALSE]
      }
      if (nrow(mat_status) > 0) {
        cat("- Areas in current matrix:", nrow(mat_status), "\n")
        if (!isTRUE(data_store$matrix_is_irregular_aggregated) && nrow(mat_status) > 200) {
          cat("  Note: many areas can make BioGeoBEARS very slow/unstable. Consider aggregating matrix to irregular subdivisions in Step 3.\n")
        }
      }
    }
    if ((input$bgb_geog_source %||% "current") == "upload") {
      cat("- Uploaded .data:", ifelse(is.null(input$bgb_geog_file), "not loaded", input$bgb_geog_file$name), "\n")
    }
    if ((input$bgb_tree_source %||% "step1") == "upload") {
      cat("- Uploaded tree:", ifelse(is.null(input$bgb_tree_file), "not loaded", input$bgb_tree_file$name), "\n")
    }
  })

  observeEvent(input$run_analysis, {
    bgb_analysis_status("Running BioGeoBEARS models...")

    tryCatch({
      if (!requireNamespace("BioGeoBEARS", quietly = TRUE)) {
        stop(
          paste(
            "Package 'BioGeoBEARS' is not installed in this R environment.",
            "Install in the SAME R session used by Shiny with:",
            "install.packages('devtools', repos='https://cloud.r-project.org')",
            "devtools::install_github('nmatzke/BioGeoBEARS', INSTALL_opts='--byte-compile', upgrade='never')",
            sep = "\n"
          )
        )
      }

      models_selected <- input$bgb_models
      if (is.null(models_selected) || length(models_selected) == 0) {
        stop("Select at least one BioGeoBEARS model.")
      }

      model_priority <- c("DEC", "DIVALIKE", "BAYAREALIKE", "DECJ", "DIVALIKEJ", "BAYAREALIKEJ")
      models_selected <- c(intersect(model_priority, models_selected), setdiff(models_selected, model_priority))

      geog_file <- tempfile(fileext = ".data")
      if ((input$bgb_geog_source %||% "current") == "upload") {
        if (is.null(input$bgb_geog_file)) stop("Please upload a BioGeoBEARS .data file.")
        file.copy(input$bgb_geog_file$datapath, geog_file, overwrite = TRUE)
      } else {
        write_geog_from_current_matrix(geog_file)
      }

      tree_file <- tempfile(fileext = ".newick")
      if ((input$bgb_tree_source %||% "step1") == "upload") {
        if (is.null(input$bgb_tree_file)) stop("Please upload a Newick tree file.")
        file.copy(input$bgb_tree_file$datapath, tree_file, overwrite = TRUE)
      } else {
        if (is.null(data_store$tree)) stop("No tree loaded in Step 1.")
        ape::write.tree(data_store$tree, file = tree_file)
      }

      taxa_sync <- harmonize_tree_and_geog_files(tree_file = tree_file, geog_file = geog_file)
      if (isTRUE(taxa_sync$changed)) {
        sync_msg <- paste0(
          "Auto-harmonized tree and geography taxa (kept ", taxa_sync$kept_n,
          "; tree-only removed: ", length(taxa_sync$missing_in_geog),
          "; geog-only removed: ", length(taxa_sync$missing_in_tree), ")."
        )
        bgb_analysis_status(paste("Running BioGeoBEARS models...", sync_msg))
      }

      ntaxa <- length(read_bgb_tree_safe(tree_file)$tip.label)
      observed_max_tipsize <- detect_max_tipsize_from_geog(geog_file)
      requested_max_range <- as.integer(input$bgb_max_range_size)
      effective_max_range <- max(requested_max_range, observed_max_tipsize)

      if (effective_max_range > requested_max_range) {
        bgb_analysis_status(
          paste0(
            "Running BioGeoBEARS models... (auto-adjusted max range size from ",
            requested_max_range,
            " to ",
            effective_max_range,
            " because at least one tip occupies ",
            observed_max_tipsize,
            " areas)"
          )
        )
      }

      bgb_constraints <- prepare_bgb_constraints(geog_file)
      previous_model_raw <- bgb_model_results_raw()

      run_config <- list(
        timestamp = Sys.time(),
        geog_source = input$bgb_geog_source %||% "current",
        tree_source = input$bgb_tree_source %||% "step1",
        selected_models = models_selected,
        max_range_size = effective_max_range,
        max_range_size_requested = requested_max_range,
        max_tipsize_observed = observed_max_tipsize,
        min_branchlength = input$bgb_min_branchlength,
        include_null_range = isTRUE(input$bgb_include_null_range),
        optimizer = input$bgb_optimizer,
        num_cores = input$bgb_num_cores,
        use_time_slices = isTRUE(input$use_time_slices),
        use_distance_matrix = isTRUE(input$use_distance_matrix),
        use_dispersal_multiplier = isTRUE(input$use_dispersal_multiplier),
        warmstart_de = isTRUE(input$bgb_use_nested_starts),
        warmstart_xwn = isTRUE(input$bgb_use_nested_starts_xwn),
        time_periods = data_store$time_periods,
        constraints_notes = bgb_constraints$notes,
        taxa_sync = taxa_sync,
        uploaded_files = list(
          distance = summarize_upload_file(input$distance_matrix_file),
          env_distance = summarize_upload_file(input$env_distance_matrix_file),
          dispersal = summarize_upload_file(input$dispersal_multipliers_file),
          areas_allowed = summarize_upload_file(input$areas_allowed_file),
          areas_adjacency = summarize_upload_file(input$areas_adjacency_file),
          time_periods = if (is.null(input$time_periods_file)) NA_character_ else input$time_periods_file$name
        )
      )

      model_rows <- list()
      model_raw <- list()
      plot_data_acc <- list(tr = NULL, tipranges = NULL)
      messages <- c()
      failed_models <- list()

      for (m in models_selected) {
        run_obj <- make_bgb_run_object(trfn = tree_file, geogfn = geog_file, constraints = bgb_constraints)
        run_obj$max_range_size <- effective_max_range
        run_obj <- configure_bgb_model(run_obj, m)

        warm_note <- NULL
        if (isTRUE(input$bgb_use_nested_starts)) {
          base_model <- get_nested_base_model(m)
          if (!is.null(base_model)) {
            nested_result <- model_raw[[base_model]]
            if (is.null(nested_result) && !is.null(previous_model_raw) && length(previous_model_raw) > 0) {
              nested_result <- previous_model_raw[[base_model]]
            }
            warm_applied <- apply_nested_warm_start(
              run_obj = run_obj,
              model_code = m,
              nested_result = nested_result,
              include_extra = isTRUE(input$bgb_use_nested_starts_xwn)
            )
            run_obj <- warm_applied$run_obj
            warm_note <- warm_applied$note
          }
        }

        run_obj <- BioGeoBEARS::readfiles_BioGeoBEARS_run(run_obj)
        if (is.character(run_obj$timesfn)) {
          run_obj <- BioGeoBEARS::section_the_tree(inputs = run_obj, make_master_table = TRUE, plot_pieces = FALSE)
        }

        run_plot_data <- extract_plot_data_from_run_obj(run_obj)
        if (is.null(plot_data_acc$tr) && !is.null(run_plot_data$tr)) {
          plot_data_acc$tr <- run_plot_data$tr
        }
        if (is.null(plot_data_acc$tipranges) && !is.null(run_plot_data$tipranges)) {
          plot_data_acc$tipranges <- run_plot_data$tipranges
        }

        model_fit <- tryCatch({
          BioGeoBEARS::check_BioGeoBEARS_run(run_obj)
          BioGeoBEARS::bears_optim_run(run_obj)
        }, error = function(e) {
          failed_models[[m]] <<- conditionMessage(e)
          NULL
        })

        used_fallback_optimizer <- FALSE
        fit_diag <- if (!is.null(model_fit)) extract_optimizer_diagnostics(model_fit) else list(convergence = NA_integer_, message = NA_character_)
        needs_fallback <- is.null(model_fit) || (!is.na(fit_diag$convergence) && fit_diag$convergence != 0)

        if (needs_fallback) {
          fallback_run_obj <- make_bgb_run_object(trfn = tree_file, geogfn = geog_file, constraints = bgb_constraints)
          fallback_run_obj$max_range_size <- effective_max_range
          fallback_run_obj <- configure_bgb_model(fallback_run_obj, m)
          fallback_run_obj <- BioGeoBEARS::readfiles_BioGeoBEARS_run(fallback_run_obj)
          if (is.character(fallback_run_obj$timesfn)) {
            fallback_run_obj <- BioGeoBEARS::section_the_tree(inputs = fallback_run_obj, make_master_table = TRUE, plot_pieces = FALSE)
          }

          fallback_plot_data <- extract_plot_data_from_run_obj(fallback_run_obj)
          if (is.null(plot_data_acc$tr) && !is.null(fallback_plot_data$tr)) {
            plot_data_acc$tr <- fallback_plot_data$tr
          }
          if (is.null(plot_data_acc$tipranges) && !is.null(fallback_plot_data$tipranges)) {
            plot_data_acc$tipranges <- fallback_plot_data$tipranges
          }

          fallback_run_obj$use_optimx <- !isTRUE(run_obj$use_optimx)

          fallback_fit <- tryCatch({
            BioGeoBEARS::check_BioGeoBEARS_run(fallback_run_obj)
            BioGeoBEARS::bears_optim_run(fallback_run_obj)
          }, error = function(e2) {
            first_msg <- failed_models[[m]] %||% "primary optimizer failed"
            failed_models[[m]] <<- paste0(first_msg, " | fallback optimizer failed: ", conditionMessage(e2))
            NULL
          })

          if (!is.null(fallback_fit)) {
            fallback_diag <- extract_optimizer_diagnostics(fallback_fit)
            prefer_fallback <- is.null(model_fit) ||
              (is.na(fit_diag$convergence) || fit_diag$convergence != 0) ||
              (!is.na(fallback_diag$convergence) && fallback_diag$convergence == 0)

            if (prefer_fallback) {
              model_fit <- fallback_fit
              fit_diag <- fallback_diag
              used_fallback_optimizer <- TRUE
            }
          }
        }

        if (is.null(model_fit)) {
          messages <- c(messages, paste0("- ", m, " failed"))
          next
        }

        lnL <- as.numeric(model_fit$total_loglikelihood)
        k <- sum(model_fit$output@params_table$type == "free")
        aic <- 2 * k - 2 * lnL
        aicc <- if ((ntaxa - k - 1) > 0) aic + (2 * k * (k + 1)) / (ntaxa - k - 1) else NA_real_

        row <- data.frame(
          Model = m,
          LnL = lnL,
          nPar = k,
          AIC = aic,
          AICc = aicc,
          d = as.numeric(model_fit$output@params_table["d", "est"]),
          e = as.numeric(model_fit$output@params_table["e", "est"]),
          j = as.numeric(model_fit$output@params_table["j", "est"]),
          Convergence = fit_diag$convergence,
          Optimizer_message = ifelse(is.na(fit_diag$message), "", fit_diag$message),
          stringsAsFactors = FALSE
        )

        model_rows[[m]] <- row
        model_raw[[m]] <- model_fit
        msg_suffix <- if (used_fallback_optimizer) " [fallback optimizer]" else ""
        warm_suffix <- if (!is.null(warm_note)) paste0(" [", warm_note, "]") else ""
        conv_suffix <- if (!is.na(fit_diag$convergence) && fit_diag$convergence != 0) {
          paste0(" [optimizer convergence code ", fit_diag$convergence, "]")
        } else {
          ""
        }
        messages <- c(messages, paste0("- ", m, " completed (LnL=", round(lnL, 3), ")", msg_suffix, warm_suffix, conv_suffix))
      }

      if (length(model_rows) == 0) {
        fail_lines <- vapply(names(failed_models), function(mm) {
          paste0("- ", mm, ": ", failed_models[[mm]])
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
        stop(paste(c("All selected BioGeoBEARS models failed.", fail_lines), collapse = "\n"))
      }

      results_df <- do.call(rbind, model_rows)
      results_df <- results_df[order(results_df$AICc, na.last = TRUE), , drop = FALSE]
      best_aicc <- suppressWarnings(min(results_df$AICc, na.rm = TRUE))
      if (is.finite(best_aicc)) {
        results_df$DeltaAICc <- results_df$AICc - best_aicc
        rel_like <- exp(-0.5 * results_df$DeltaAICc)
        denom <- sum(rel_like, na.rm = TRUE)
        results_df$AICcWt <- if (is.finite(denom) && denom > 0) rel_like / denom else NA_real_
      } else {
        results_df$DeltaAICc <- NA_real_
        results_df$AICcWt <- NA_real_
      }

      bgb_model_results_table(results_df)
      bgb_model_results_raw(model_raw)
      last_inputs_obj <- list(
        trfn = tree_file,
        geogfn = geog_file,
        include_null_range = isTRUE(input$bgb_include_null_range),
        extrap_method = data_store$extrap_method,
        tr = plot_data_acc$tr,
        tipranges = plot_data_acc$tipranges
      )
      bgb_last_inputs(last_inputs_obj)
      file_plot_data <- build_bgb_plot_data_cache(last_inputs_obj)
      bgb_last_plot_data(list(
        tr = plot_data_acc$tr %||% file_plot_data$tr,
        tipranges = plot_data_acc$tipranges %||% file_plot_data$tipranges
      ))
      run_config$model_comparison <- results_df
      bgb_last_run_config(run_config)
      if (!is.null(bgb_constraints$notes) && length(bgb_constraints$notes) > 0) {
        messages <- c(messages, paste0("Constraints: ", paste(bgb_constraints$notes, collapse = "; ")))
      }
      fail_lines <- character(0)
      if (length(failed_models) > 0) {
        fail_lines <- vapply(names(failed_models), function(mm) {
          paste0("- ", mm, ": ", failed_models[[mm]])
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
      }
      status_header <- if (length(failed_models) > 0) {
        "BioGeoBEARS completed with partial failures."
      } else {
        "BioGeoBEARS completed."
      }
      bgb_analysis_status(paste(c(status_header, messages, fail_lines), collapse = "\n"))
    }, error = function(e) {
      bgb_analysis_status(paste0("BioGeoBEARS run failed: ", e$message))
      bgb_model_results_table(NULL)
      bgb_model_results_raw(list())
      bgb_last_inputs(NULL)
      bgb_last_plot_data(NULL)
      bgb_last_run_config(NULL)
    })
  })

  output$analysis_status <- renderText({
    bgb_analysis_status()
  })

  output$model_comparison_table <- DT::renderDataTable({
    res <- bgb_model_results_table()
    if (is.null(res)) return(NULL)
    DT::datatable(res, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$download_model_comparison <- downloadHandler(
    filename = function() {
      paste0("biogeobears_model_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      res <- bgb_model_results_table()
      if (is.null(res) || nrow(res) == 0) {
        write.csv(data.frame(message = "No model-comparison results available. Run BioGeoBEARS first."), file, row.names = FALSE)
        return(invisible(NULL))
      }
      write.csv(res, file, row.names = FALSE)
    }
  )

  output$download_bgb_models_rdata <- downloadHandler(
    filename = function() {
      paste0("biogeobears_fitted_models_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData")
    },
    content = function(file) {
      model_results <- bgb_model_results_raw()
      model_comparison <- bgb_model_results_table()
      bgb_inputs <- bgb_last_inputs()

      if (is.null(model_results) || length(model_results) == 0) {
        tmp_msg <- data.frame(message = "No fitted BioGeoBEARS models available. Run analysis first.", stringsAsFactors = FALSE)
        save(tmp_msg, file = file)
        return(invisible(NULL))
      }

      export_timestamp <- Sys.time()
      area_code_mapping <- area_code_mapping_snapshot()

      save(
        model_results,
        model_comparison,
        bgb_inputs,
        area_code_mapping,
        export_timestamp,
        file = file
      )
    }
  )

  output$download_bgb_run_report <- downloadHandler(
    filename = function() {
      paste0("biogeobears_run_config_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
    },
    content = function(file) {
      cfg <- bgb_last_run_config()
      if (is.null(cfg)) {
        writeLines("No BioGeoBEARS run configuration available. Run analysis first.", con = file)
        return(invisible(NULL))
      }

      lines <- c(
        "BioGeoBEARS Run Configuration Report",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("Run timestamp: ", format(cfg$timestamp, "%Y-%m-%d %H:%M:%S")),
        "",
        "[Inputs]",
        paste0("Geography source: ", cfg$geog_source),
        paste0("Tree source: ", cfg$tree_source),
        paste0("Models: ", paste(cfg$selected_models, collapse = ", ")),
        paste0("Max range size: ", cfg$max_range_size),
        paste0("Min branch length: ", cfg$min_branchlength),
        paste0("Include null range: ", cfg$include_null_range),
        paste0("Optimizer: ", cfg$optimizer),
        paste0("Cores: ", cfg$num_cores),
        "",
        "[Options]",
        paste0("Use time slices: ", cfg$use_time_slices),
        paste0("Use distance matrix (x): ", cfg$use_distance_matrix),
        paste0("Use dispersal multiplier (w): ", cfg$use_dispersal_multiplier),
        paste0("Warm-start d/e: ", cfg$warmstart_de),
        paste0("Warm-start x/w/n: ", cfg$warmstart_xwn),
        ""
      )

      if (!is.null(cfg$taxa_sync)) {
        lines <- c(
          lines,
          "[Taxa harmonization]",
          paste0("Auto-harmonized: ", isTRUE(cfg$taxa_sync$changed)),
          paste0("Kept taxa: ", cfg$taxa_sync$kept_n %||% NA_integer_),
          paste0("Tree-only removed: ", length(cfg$taxa_sync$missing_in_geog %||% character(0))),
          paste0("Geog-only removed: ", length(cfg$taxa_sync$missing_in_tree %||% character(0))),
          ""
        )
      }

      if (!is.null(cfg$time_periods) && length(cfg$time_periods) > 0) {
        lines <- c(lines, "[Time periods]", paste(cfg$time_periods, collapse = ", "), "")
      }

      if (!is.null(cfg$constraints_notes) && length(cfg$constraints_notes) > 0) {
        lines <- c(lines, "[Applied constraints]", paste0("- ", cfg$constraints_notes), "")
      }

      uf <- cfg$uploaded_files
      lines <- c(
        lines,
        "[Uploaded files]",
        paste0("Distance: ", uf$distance$name %||% "NA", " | blocks=", uf$distance$blocks %||% NA_integer_, " | END=", uf$distance$has_end %||% FALSE),
        paste0("Env distance: ", uf$env_distance$name %||% "NA", " | blocks=", uf$env_distance$blocks %||% NA_integer_, " | END=", uf$env_distance$has_end %||% FALSE),
        paste0("Dispersal: ", uf$dispersal$name %||% "NA", " | blocks=", uf$dispersal$blocks %||% NA_integer_, " | END=", uf$dispersal$has_end %||% FALSE),
        paste0("Areas allowed: ", uf$areas_allowed$name %||% "NA", " | blocks=", uf$areas_allowed$blocks %||% NA_integer_, " | END=", uf$areas_allowed$has_end %||% FALSE),
        paste0("Areas adjacency: ", uf$areas_adjacency$name %||% "NA", " | blocks=", uf$areas_adjacency$blocks %||% NA_integer_, " | END=", uf$areas_adjacency$has_end %||% FALSE),
        paste0("Time periods: ", uf$time_periods %||% "NA"),
        ""
      )

      if (!is.null(cfg$model_comparison) && nrow(cfg$model_comparison) > 0) {
        lines <- c(lines, "[Model comparison]", utils::capture.output(print(cfg$model_comparison, row.names = FALSE)))
      }

      writeLines(lines, con = file)
    }
  )

  observe({
    res <- bgb_model_results_table()
    choices <- if (!is.null(res) && nrow(res) > 0) res$Model else c("DEC", "DECJ")
    updateSelectInput(session, "lrt_model1", choices = choices, selected = choices[1])
    updateSelectInput(session, "lrt_model2", choices = choices, selected = choices[min(2, length(choices))])
    viz_choices <- if (!is.null(res) && nrow(res) > 0) res$Model else c("(Run BioGeoBEARS first)" = "")
    viz_selected <- if (!is.null(res) && nrow(res) > 0) res$Model[1] else ""
    updateSelectInput(session, "bgb_visual_model", choices = viz_choices, selected = viz_selected)
  })

  output$ancestral_ranges_plot <- renderPlot({
    req(nzchar(input$bgb_visual_model))
    raw <- bgb_model_results_raw()
    inputs <- bgb_last_inputs() %||% list()
    req(!is.null(raw), !is.null(raw[[input$bgb_visual_model]]))
    bgb_plot_error_text(NULL)
    bgb_plot_diag_text(NULL)

    res_obj <- raw[[input$bgb_visual_model]]

    tryCatch({
      plot_data <- bgb_last_plot_data()
      plot_call <- function() {
        plot_bgb_results_safe(
          res_obj = res_obj,
          model_name = input$bgb_visual_model,
          inputs = inputs,
          plot_data = plot_data,
          plotwhat = "text",
          title_suffix = "Ancestral ranges",
          label_offset_override = input$ancestral_ranges_label_offset,
          tipcex_override = input$ancestral_ranges_tip_cex_app,
          statecex_override = input$ancestral_ranges_state_cex_app
        )
      }
      if (isTRUE(should_use_bgb_replay(inputs))) {
        draw_bgb_plot_with_replay(plot_call)
      } else {
        plot_call()
      }
    }, error = function(e) {
      bgb_plot_error_text(conditionMessage(e))
      bgb_plot_diag_text(NULL)
      plot.new()
      title(main = paste0("Could not plot ancestral ranges (", input$bgb_visual_model, ")"))
      text(0.5, 0.5, labels = conditionMessage(e), cex = 0.9)
    })
  })

  output$range_uncertainty_plot <- renderPlot({
    req(nzchar(input$bgb_visual_model))
    raw <- bgb_model_results_raw()
    inputs <- bgb_last_inputs() %||% list()
    req(!is.null(raw), !is.null(raw[[input$bgb_visual_model]]))
    bgb_pie_error_text(NULL)

    res_obj <- raw[[input$bgb_visual_model]]

    tryCatch({
      plot_data <- bgb_last_plot_data()
      plot_call <- function() {
        plot_bgb_results_safe(
          res_obj = res_obj,
          model_name = input$bgb_visual_model,
          inputs = inputs,
          plot_data = plot_data,
          plotwhat = "pie",
          title_suffix = "Range uncertainty",
          label_offset_override = input$range_uncertainty_label_offset,
          tipcex_override = input$range_uncertainty_tip_cex_app,
          statecex_override = input$range_uncertainty_state_cex_app
        )
      }
      if (isTRUE(should_use_bgb_replay(inputs))) {
        draw_bgb_plot_with_replay(plot_call)
      } else {
        plot_call()
      }
    }, error = function(e) {
      bgb_pie_error_text(conditionMessage(e))
      plot.new()
      title(main = paste0("Could not plot uncertainty pies (", input$bgb_visual_model, ")"))
      text(0.5, 0.5, labels = conditionMessage(e), cex = 0.9)
    })
  })

  output$ancestral_ranges_status <- renderText({
    if (!nzchar(input$bgb_visual_model %||% "")) {
      return("Run BioGeoBEARS in Step 7 and select a model to visualize ancestral ranges.")
    }
    res <- bgb_model_results_table()
    warn <- ""
    if (!is.null(res) && "Convergence" %in% names(res)) {
      idx <- which(res$Model == input$bgb_visual_model)[1]
      if (!is.na(idx) && !is.na(res$Convergence[idx]) && res$Convergence[idx] != 0) {
        warn <- paste0(" (warning: optimizer convergence code ", res$Convergence[idx], ")")
      }
    }
    err <- bgb_plot_error_text()
    diag_txt <- bgb_plot_diag_text()
    if (!is.null(err) && nzchar(err)) {
      return(paste0("Ancestral ranges plot error: ", err))
    }
    base_txt <- paste0("Displaying ancestral ranges for model: ", input$bgb_visual_model, warn)
    if (!is.null(diag_txt) && nzchar(diag_txt)) {
      paste0(base_txt, " | ", diag_txt)
    } else {
      base_txt
    }
  })

  output$range_uncertainty_status <- renderText({
    if (!nzchar(input$bgb_visual_model %||% "")) {
      return("Run BioGeoBEARS in Step 7 and select a model to visualize uncertainty pies.")
    }
    err <- bgb_pie_error_text()
    if (!is.null(err) && nzchar(err)) {
      return(paste0("Range uncertainty plot error: ", err))
    }
    paste0("Displaying range uncertainty (pie charts) for model: ", input$bgb_visual_model)
  })

  output$download_bgb_ancestral_pdf <- downloadHandler(
    filename = function() {
      mdl <- input$bgb_visual_model %||% "model"
      mdl <- gsub("[^A-Za-z0-9_]+", "_", mdl)
      paste0("BioGeoBEARS_ancestral_", mdl, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      raw <- bgb_model_results_raw()
      inputs <- bgb_last_inputs() %||% list()
      plot_data <- bgb_last_plot_data()
      model_name <- input$bgb_visual_model %||% ""

      grDevices::pdf(file = file, width = 12, height = 8)
      on.exit(grDevices::dev.off(), add = TRUE)

      if (!nzchar(model_name) || is.null(raw) || is.null(raw[[model_name]])) {
        plot.new()
        text(0.5, 0.6, "No BioGeoBEARS model available for export", cex = 1.2)
        text(0.5, 0.5, "Run Step 7 and select a model in Step 4.", cex = 1)
        return(invisible(NULL))
      }

      res_obj <- raw[[model_name]]
      plot_bgb_results_safe(
        res_obj = res_obj,
        model_name = model_name,
        inputs = inputs,
        plot_data = plot_data,
        plotwhat = "text",
        title_suffix = "Ancestral ranges",
        label_offset_override = input$ancestral_ranges_label_offset
      )

      plot.new()
      plot_bgb_results_safe(
        res_obj = res_obj,
        model_name = model_name,
        inputs = inputs,
        plot_data = plot_data,
        plotwhat = "pie",
        title_suffix = "Range uncertainty",
        label_offset_override = input$range_uncertainty_label_offset,
        tipcex_override = input$range_uncertainty_tip_cex_app,
        statecex_override = input$range_uncertainty_state_cex_app
      )
    }
  )

  output$model_comparison_plot <- plotly::renderPlotly({
    res <- bgb_model_results_table()
    req(!is.null(res), nrow(res) > 0)

    plotly::plot_ly(
      data = res,
      x = ~Model,
      y = ~AICc,
      type = "bar",
      hovertemplate = "Model: %{x}<br>AICc: %{y:.3f}<extra></extra>",
      marker = list(color = "#2C7FB8")
    ) %>%
      plotly::layout(title = "Model Comparison (AICc)", xaxis = list(title = "Model"), yaxis = list(title = "AICc"))
  })

  output$lrt_results <- renderPrint({
    req(input$run_lrt)
    res <- bgb_model_results_table()
    if (is.null(res) || nrow(res) == 0) {
      cat("No model results available. Run BioGeoBEARS first.\n")
      return(invisible(NULL))
    }

    m1 <- input$lrt_model1
    m2 <- input$lrt_model2
    r1 <- res[res$Model == m1, , drop = FALSE]
    r2 <- res[res$Model == m2, , drop = FALSE]

    if (nrow(r1) == 0 || nrow(r2) == 0) {
      cat("Selected model(s) not found in current results.\n")
      return(invisible(NULL))
    }

    base <- if (r1$nPar <= r2$nPar) r1 else r2
    full <- if (r1$nPar <= r2$nPar) r2 else r1
    df <- full$nPar - base$nPar
    if (df <= 0) {
      cat("LRT requires nested models with different number of parameters.\n")
      return(invisible(NULL))
    }

    LR <- 2 * (full$LnL - base$LnL)
    pval <- stats::pchisq(LR, df = df, lower.tail = FALSE)

    cat("LRT Results:\n")
    cat("Base model:", base$Model, "(k=", base$nPar, ")\n", sep = "")
    cat("Full model:", full$Model, "(k=", full$nPar, ")\n", sep = "")
    cat("LR statistic:", round(LR, 4), "\n")
    cat("df:", df, "\n")
    cat("p-value:", signif(pval, 4), "\n")
  })
  

  # ===== MATRIX LOADING HANDLERS =====
  
  # Helper function to parse TNT format
  parse_tnt_matrix <- function(file_path) {
    lines <- readLines(file_path)
    header <- lines[1]
    header_parts <- strsplit(header, "\\s+")[[1]]
    n_taxa <- as.numeric(header_parts[1])
    n_chars <- as.numeric(header_parts[2])
    
    data_lines <- lines[2:(n_taxa + 1)]
    taxa_names <- character(n_taxa)
    matrix_data <- matrix(nrow = n_taxa, ncol = n_chars)
    
    for (i in seq_along(data_lines)) {
      parts <- strsplit(data_lines[i], "\\s+")[[1]]
      taxa_names[i] <- parts[1]
      matrix_data[i, ] <- as.numeric(strsplit(parts[2], "")[[1]])
    }
    
    rownames(matrix_data) <- taxa_names
    colnames(matrix_data) <- paste0("Char", 1:n_chars)
    
    return(list(matrix = matrix_data, n_taxa = n_taxa, n_chars = n_chars))
  }

  load_named_square_matrix <- function(file_info, matrix_label = "Matrix", require_binary = FALSE, require_positive_offdiag = FALSE) {
    if (is.null(file_info)) {
      stop("Please select a file")
    }

    split_tokens <- function(line) {
      x <- trimws(gsub("\\r", "", line))
      if (!nzchar(x)) {
        return(character(0))
      }
      if (grepl("\t", x)) {
        tok <- strsplit(x, "\t", fixed = TRUE)[[1]]
      } else if (grepl(",", x, fixed = TRUE)) {
        tok <- strsplit(x, ",", fixed = TRUE)[[1]]
      } else {
        tok <- strsplit(x, "\\s+")[[1]]
      }
      tok <- trimws(tok)
      tok[nzchar(tok)]
    }

    raw_lines <- readLines(file_info$datapath, warn = FALSE)
    clean_lines <- trimws(gsub("\\r", "", raw_lines))
    clean_lines <- clean_lines[nzchar(clean_lines)]
    clean_lines <- clean_lines[toupper(clean_lines) != "END"]

    if (length(clean_lines) < 2) {
      stop(paste0(matrix_label, " file is empty or incomplete."))
    }

    header_tokens <- split_tokens(clean_lines[1])
    n_areas <- length(header_tokens)
    if (n_areas < 2) {
      stop(paste0(matrix_label, " header must contain at least 2 area names."))
    }

    rows_numeric <- list()
    for (i in seq_along(clean_lines[-1])) {
      line_idx <- i + 1
      tok <- split_tokens(clean_lines[line_idx])
      if (length(tok) == 0) {
        next
      }

      if (length(tok) == (n_areas + 1)) {
        first_num <- suppressWarnings(as.numeric(tok[1]))
        if (is.na(first_num)) {
          tok <- tok[-1]
        }
      }

      if (length(tok) != n_areas) {
        next
      }

      nums <- suppressWarnings(as.numeric(tok))
      if (any(is.na(nums))) {
        next
      }

      rows_numeric[[length(rows_numeric) + 1]] <- nums
      if (length(rows_numeric) == n_areas) {
        break
      }
    }

    if (length(rows_numeric) != n_areas) {
      stop(
        paste0(
          matrix_label,
          " could not be parsed as a square matrix. Expected ",
          n_areas,
          " numeric rows matching the header."
        )
      )
    }

    mat <- do.call(rbind, rows_numeric)
    rownames(mat) <- header_tokens
    colnames(mat) <- header_tokens

    if (nrow(mat) != ncol(mat)) {
      stop(paste0(matrix_label, " must be square (same number of rows and columns)."))
    }
    if (any(!is.finite(mat), na.rm = TRUE)) {
      stop(paste0(matrix_label, " contains non-finite values."))
    }

    if (require_binary && !all(mat %in% c(0, 1, NA))) {
      stop(paste0(matrix_label, " should contain only 0 and 1 values."))
    }

    if (require_positive_offdiag) {
      tmp <- mat
      diag(tmp) <- NA_real_
      if (any(tmp <= 0, na.rm = TRUE)) {
        stop(paste0(matrix_label, " requires off-diagonal values > 0."))
      }
    }

    mat
  }
  
  load_distance_matrix_handler <- function(file_info) {
    tryCatch({
      if (is.null(file_info)) {
        output$distance_matrix_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        matrix_data <- load_named_square_matrix(
          file_info = file_info,
          matrix_label = "Distance matrix",
          require_binary = FALSE,
          require_positive_offdiag = TRUE
        )
        data_store$distance_matrix <- matrix_data
        
        output$distance_matrix_status <- renderPrint({
          cat("✓ Distance matrix loaded!\n")
          cat("Dimensions:", nrow(matrix_data), "x", ncol(matrix_data), "\n")
          cat("Areas:", paste(rownames(matrix_data), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$distance_matrix_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  }

  # Load distance matrix
  observeEvent(input$load_distance_matrix, {
    load_distance_matrix_handler(input$distance_matrix_file)
  })

  observeEvent(input$distance_matrix_file, {
    load_distance_matrix_handler(input$distance_matrix_file)
  })

  load_env_distance_matrix_handler <- function(file_info) {
    tryCatch({
      if (is.null(file_info)) {
        output$env_distance_matrix_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        matrix_data <- load_named_square_matrix(
          file_info = file_info,
          matrix_label = "Environmental distance matrix",
          require_binary = FALSE,
          require_positive_offdiag = TRUE
        )
        data_store$env_distance_matrix <- matrix_data

        output$env_distance_matrix_status <- renderPrint({
          cat("✓ Environmental distance matrix loaded!\n")
          cat("Dimensions:", nrow(matrix_data), "x", ncol(matrix_data), "\n")
          cat("Areas:", paste(rownames(matrix_data), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$env_distance_matrix_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  }

  # Load environmental distance matrix
  observeEvent(input$load_env_distance_matrix, {
    load_env_distance_matrix_handler(input$env_distance_matrix_file)
  })

  observeEvent(input$env_distance_matrix_file, {
    load_env_distance_matrix_handler(input$env_distance_matrix_file)
  })
  
  load_dispersal_multipliers_handler <- function(file_info) {
    tryCatch({
      if (is.null(file_info)) {
        output$dispersal_multipliers_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        matrix_data <- load_named_square_matrix(
          file_info = file_info,
          matrix_label = "Dispersal multipliers matrix",
          require_binary = TRUE,
          require_positive_offdiag = FALSE
        )
        
        data_store$dispersal_multipliers <- matrix_data
        
        output$dispersal_multipliers_status <- renderPrint({
          cat("✓ Dispersal multipliers loaded!\n")
          cat("Dimensions:", nrow(matrix_data), "x", ncol(matrix_data), "\n")
          cat("Areas:", paste(rownames(matrix_data), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$dispersal_multipliers_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  }

  # Load dispersal multipliers
  observeEvent(input$load_dispersal_multipliers, {
    load_dispersal_multipliers_handler(input$dispersal_multipliers_file)
  })

  observeEvent(input$dispersal_multipliers_file, {
    load_dispersal_multipliers_handler(input$dispersal_multipliers_file)
  })

  load_areas_allowed_handler <- function(file_info) {
    tryCatch({
      if (is.null(file_info)) {
        output$areas_allowed_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        matrix_data <- load_named_square_matrix(
          file_info = file_info,
          matrix_label = "Areas-allowed matrix",
          require_binary = TRUE,
          require_positive_offdiag = FALSE
        )
        data_store$areas_allowed_matrix <- matrix_data

        output$areas_allowed_status <- renderPrint({
          cat("✓ Areas-allowed matrix loaded!\n")
          cat("Dimensions:", nrow(matrix_data), "x", ncol(matrix_data), "\n")
          cat("Areas:", paste(rownames(matrix_data), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$areas_allowed_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  }

  # Load areas-allowed matrix
  observeEvent(input$load_areas_allowed, {
    load_areas_allowed_handler(input$areas_allowed_file)
  })

  observeEvent(input$areas_allowed_file, {
    load_areas_allowed_handler(input$areas_allowed_file)
  })

  load_areas_adjacency_handler <- function(file_info) {
    tryCatch({
      if (is.null(file_info)) {
        output$areas_adjacency_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        matrix_data <- load_named_square_matrix(
          file_info = file_info,
          matrix_label = "Areas-adjacency matrix",
          require_binary = TRUE,
          require_positive_offdiag = FALSE
        )
        data_store$areas_adjacency_matrix <- matrix_data

        output$areas_adjacency_status <- renderPrint({
          cat("✓ Areas-adjacency matrix loaded!\n")
          cat("Dimensions:", nrow(matrix_data), "x", ncol(matrix_data), "\n")
          cat("Areas:", paste(rownames(matrix_data), collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$areas_adjacency_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  }

  # Load areas-adjacency matrix
  observeEvent(input$load_areas_adjacency, {
    load_areas_adjacency_handler(input$areas_adjacency_file)
  })

  observeEvent(input$areas_adjacency_file, {
    load_areas_adjacency_handler(input$areas_adjacency_file)
  })
  
  # Load traits matrix (TNT format)
  observeEvent(input$load_traits_matrix, {
    tryCatch({
      if (is.null(input$traits_matrix_file)) {
        output$traits_matrix_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        tnt_data <- parse_tnt_matrix(input$traits_matrix_file$datapath)
        data_store$traits_matrix <- tnt_data$matrix
        
        output$traits_matrix_status <- renderPrint({
          cat("✓ Traits matrix loaded!\n")
          cat("Taxa:", tnt_data$n_taxa, "\n")
          cat("Characters:", tnt_data$n_chars, "\n")
          cat("Species:", paste(rownames(tnt_data$matrix)[1:min(5, nrow(tnt_data$matrix))], collapse = ", "))
          if (tnt_data$n_taxa > 5) cat(", ...")
          cat("\n")
        })
      }
    }, error = function(e) {
      output$traits_matrix_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  })
  
  # Load constraint matrix (TNT format)
  observeEvent(input$load_constraint_matrix, {
    tryCatch({
      if (is.null(input$constraint_matrix_file)) {
        output$constraint_matrix_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        tnt_data <- parse_tnt_matrix(input$constraint_matrix_file$datapath)
        data_store$constraint_matrix <- tnt_data$matrix
        
        output$constraint_matrix_status <- renderPrint({
          cat("✓ Constraint matrix loaded!\n")
          cat("Taxa:", tnt_data$n_taxa, "\n")
          cat("Characters:", tnt_data$n_chars, "\n")
          cat("Species:", paste(rownames(tnt_data$matrix)[1:min(5, nrow(tnt_data$matrix))], collapse = ", "))
          if (tnt_data$n_taxa > 5) cat(", ...")
          cat("\n")
        })
      }
    }, error = function(e) {
      output$constraint_matrix_status <- renderPrint({
        cat("Error loading matrix:\n", e$message, "\n")
      })
    })
  })
  
  # Load time periods
  observeEvent(input$load_time_periods, {
    tryCatch({
      if (is.null(input$time_periods_file)) {
        output$time_periods_status <- renderPrint({
          cat("Please select a file\n")
        })
      } else {
        time_data <- as.numeric(readLines(input$time_periods_file$datapath))
        time_data <- time_data[!is.na(time_data)]
        time_data <- time_data[time_data > 0]
        if (length(time_data) == 0) {
          stop("Time periods file must contain positive values (> 0).")
        }
        data_store$time_periods <- sort(unique(time_data), decreasing = FALSE)
        
        output$time_periods_status <- renderPrint({
          cat("✓ Time periods loaded!\n")
          cat("Number of periods:", length(data_store$time_periods), "\n")
          cat("Range:", min(data_store$time_periods), "-", max(data_store$time_periods), "Ma\n")
          cat("Periods (Ma, youngest to oldest):", paste(data_store$time_periods, collapse = ", "), "\n")
        })
      }
    }, error = function(e) {
      output$time_periods_status <- renderPrint({
        cat("Error loading time periods:\n", e$message, "\n")
      })
    })
  })

  observeEvent(input$time_periods_file, {
    if (!is.null(input$time_periods_file)) {
      tryCatch({
        time_data <- as.numeric(readLines(input$time_periods_file$datapath))
        time_data <- time_data[!is.na(time_data)]
        time_data <- time_data[time_data > 0]
        if (length(time_data) == 0) {
          stop("Time periods file must contain positive values (> 0).")
        }
        data_store$time_periods <- sort(unique(time_data), decreasing = FALSE)
        output$time_periods_status <- renderPrint({
          cat("✓ Time periods loaded!\n")
          cat("Number of periods:", length(data_store$time_periods), "\n")
          cat("Range:", min(data_store$time_periods), "-", max(data_store$time_periods), "Ma\n")
          cat("Periods (Ma, youngest to oldest):", paste(data_store$time_periods, collapse = ", "), "\n")
        })
      }, error = function(e) {
        output$time_periods_status <- renderPrint({
          cat("Error loading time periods:\n", e$message, "\n")
        })
      })
    }
  })
  
  # Summary of loaded matrices
  output$matrices_summary <- renderPrint({
    cat("=== LOADED MATRICES SUMMARY ===\n\n")
    
    if (!is.null(data_store$distance_matrix)) {
      cat("✓ Distance Matrix:", nrow(data_store$distance_matrix), "x", ncol(data_store$distance_matrix), "\n")
    } else {
      cat("✗ Distance Matrix: Not loaded\n")
    }

    if (!is.null(data_store$env_distance_matrix)) {
      cat("✓ Environmental Distance Matrix:", nrow(data_store$env_distance_matrix), "x", ncol(data_store$env_distance_matrix), "\n")
    } else {
      cat("✗ Environmental Distance Matrix: Not loaded\n")
    }
    
    if (!is.null(data_store$dispersal_multipliers)) {
      cat("✓ Dispersal Multipliers:", nrow(data_store$dispersal_multipliers), "x", ncol(data_store$dispersal_multipliers), "\n")
    } else {
      cat("✗ Dispersal Multipliers: Not loaded\n")
    }

    if (!is.null(data_store$areas_allowed_matrix)) {
      cat("✓ Areas-Allowed Matrix:", nrow(data_store$areas_allowed_matrix), "x", ncol(data_store$areas_allowed_matrix), "\n")
    } else {
      cat("✗ Areas-Allowed Matrix: Not loaded\n")
    }

    if (!is.null(data_store$areas_adjacency_matrix)) {
      cat("✓ Areas-Adjacency Matrix:", nrow(data_store$areas_adjacency_matrix), "x", ncol(data_store$areas_adjacency_matrix), "\n")
    } else {
      cat("✗ Areas-Adjacency Matrix: Not loaded\n")
    }
    
    if (!is.null(data_store$traits_matrix)) {
      cat("✓ Traits Matrix:", nrow(data_store$traits_matrix), "x", ncol(data_store$traits_matrix), "\n")
    } else {
      cat("✗ Traits Matrix: Not loaded\n")
    }
    
    if (!is.null(data_store$constraint_matrix)) {
      cat("✓ Constraint Matrix:", nrow(data_store$constraint_matrix), "x", ncol(data_store$constraint_matrix), "\n")
    } else {
      cat("✗ Constraint Matrix: Not loaded\n")
    }
    
    if (!is.null(data_store$time_periods)) {
      cat("✓ Time Periods:", length(data_store$time_periods), "periods\n")
    } else {
      cat("✗ Time Periods: Not loaded\n")
    }
  })

  output$bgb_strat_validation <- renderPrint({
    cat("=== TIME-STRATIFIED VALIDATION ===\n\n")

    use_ts <- isTRUE(input$use_time_slices)
    cat("Use time slices:", use_ts, "\n")

    if (!use_ts) {
      cat("- Status: 🟡 time stratification is disabled\n")
      return(invisible(NULL))
    }

    n_time <- if (!is.null(data_store$time_periods)) length(data_store$time_periods) else 0L
    cat("Expected strata (from time periods):", n_time, "\n")

    if (n_time == 0) {
      cat("- Status: 🔴 no time periods loaded\n")
      return(invisible(NULL))
    }

    status_counts <- c(green = 0L, yellow = 0L, red = 0L)

    check_file <- function(label, file_info, required = FALSE) {
      if (is.null(file_info)) {
        if (required) {
          cat("🔴 ", label, ": not provided (required)\n", sep = "")
          status_counts[["red"]] <<- status_counts[["red"]] + 1L
        } else {
          cat("🟡 ", label, ": not provided (optional)\n", sep = "")
          status_counts[["yellow"]] <<- status_counts[["yellow"]] + 1L
        }
        return(invisible(NULL))
      }

      blocks <- count_bgb_matrix_blocks(file_info)
      has_end <- file_looks_like_bgb_matrix_file(file_info)
      ok <- !is.na(blocks) && blocks >= n_time

      if (ok) {
        cat("🟢 ", label, ": ", file_info$name,
            " | blocks=", blocks, " | END=", has_end,
            " | needs >=", n_time, "\n", sep = "")
        status_counts[["green"]] <<- status_counts[["green"]] + 1L
      } else {
        cat("🔴 ", label, ": ", file_info$name,
            " | blocks=", blocks, " | END=", has_end,
            " | needs >=", n_time, "\n", sep = "")
        status_counts[["red"]] <<- status_counts[["red"]] + 1L
      }
    }

    check_file("Dispersal multipliers", input$dispersal_multipliers_file, required = isTRUE(input$use_dispersal_multiplier))
    check_file("Areas allowed", input$areas_allowed_file, required = FALSE)
    check_file("Areas adjacency", input$areas_adjacency_file, required = FALSE)
    if (isTRUE(input$use_distance_matrix)) {
      check_file("Distance matrix", input$distance_matrix_file, required = TRUE)
    }
    if (!is.null(input$env_distance_matrix_file)) {
      check_file("Environmental distance matrix", input$env_distance_matrix_file, required = FALSE)
    }

    cat("\nSummary: ",
        "🟢=", status_counts[["green"]], " ",
        "🟡=", status_counts[["yellow"]], " ",
        "🔴=", status_counts[["red"]], "\n", sep = "")

    if (status_counts[["red"]] > 0) {
      cat("Overall: 🔴 inconsistent/missing required stratified inputs\n")
    } else if (status_counts[["yellow"]] > 0) {
      cat("Overall: 🟡 runnable, with optional files missing\n")
    } else {
      cat("Overall: 🟢 ready for time-stratified run\n")
    }
  })

}
