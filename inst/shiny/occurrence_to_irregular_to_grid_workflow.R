# ==============================================================================
# ROTINA PASSO A PASSO:
# Occurrence Points -> Irregular Polygons -> Regular Grids
# ==============================================================================
#
# LÓGICA CENTRAL:
# Se uma espécie tem um ponto de ocorrência dentro de um polígono irregular
# (ex: uma província biogeográfica), TODAS as células de grid regular que
# intersectam esse polígono recebem presença = 1 para essa espécie.
#
# Resultado: Matriz de presença/ausência (linhas = grid cells, colunas = espécies)
# pronta para PAE-PCE.
#
# Execute linha por linha no RStudio.
# ==============================================================================

# ==============================================================================
# PASSO 0: Carregar pacotes
# ==============================================================================

library(terra)
library(sf)
library(raster)
library(dplyr)

# ==============================================================================
# PASSO 1: Carregar seus dados
# ==============================================================================

# --- 1A: Dados de ocorrência (CSV/TXT com colunas: spp, long, lat) ---
# Substitua pelo caminho do seu arquivo:
occ_data <- read.table("caminho/para/seu/arquivo_ocorrencias.txt", sep = "\t", dec = ".", header = TRUE)

# Verificar estrutura
head(occ_data)
# Deve ter pelo menos: spp, long, lat (ou SPECIES, LONG, LAT)
# Se suas colunas têm nomes diferentes, renomeie:
# names(occ_data)[names(occ_data) == "SPECIES"] <- "spp"
# names(occ_data)[names(occ_data) == "LONG"] <- "long"
# names(occ_data)[names(occ_data) == "LAT"] <- "lat"

cat("Pontos de ocorrência carregados:", nrow(occ_data), "\n")
cat("Espécies únicas:", length(unique(occ_data$spp)), "\n")
cat("Espécies:", paste(unique(occ_data$spp), collapse = ", "), "\n\n")

# --- 1B: Shapefile de polígonos irregulares (ex: províncias Morrone) ---
# Carregamento via terra::vect()
study_area_vect <- terra::vect("caminho/para/seu/shapefile.shp", crs = "+proj=longlat +datum=WGS84")

# Verificar CRS e colunas
print(terra::crs(study_area_vect))
print(names(study_area_vect))

# --- 1C: Definir qual coluna do shapefile identifica os polígonos ---
# Olhe os nomes das colunas acima e escolha a coluna de ID:
polygon_id_column <- "Provincias"  # <-- MUDE para o nome da sua coluna

# Verificar os polígonos disponíveis
cat("\nPolígonos encontrados:", length(unique(study_area_vect[[polygon_id_column]])), "\n")
print(unique(study_area_vect[[polygon_id_column]]))

# ==============================================================================
# PASSO 2: Converter terra::vect para sf e validar geometrias
# ==============================================================================
# Convertemos para sf porque st_join é mais eficiente para spatial joins massivos.
# O shapefile original permanece intacto como SpatVector.

study_area <- sf::st_as_sf(study_area_vect)
study_area <- sf::st_make_valid(study_area)
study_area <- study_area[!sf::st_is_empty(study_area), ]

cat("\nGeometrias validadas. Polígonos restantes:", nrow(study_area), "\n")

# Garantir CRS WGS84
if (is.na(sf::st_crs(study_area))) {
  sf::st_crs(study_area) <- 4326
}

# ==============================================================================
# PASSO 3: Converter pontos de ocorrência para objeto sf
# ==============================================================================

# Remover linhas com NA nas coordenadas
occ_data <- occ_data[complete.cases(occ_data[, c("long", "lat")]), ]

# Converter para sf
occ_sf <- sf::st_as_sf(occ_data, coords = c("long", "lat"), crs = 4326)

# Transformar para o mesmo CRS do shapefile (caso sejam diferentes)
occ_sf <- sf::st_transform(occ_sf, sf::st_crs(study_area))

cat("Pontos convertidos para sf:", nrow(occ_sf), "\n\n")

# ==============================================================================
# PASSO 4: Encontrar em qual polígono cada ponto cai (spatial join)
# ==============================================================================
# Para cada ponto, identificar qual polígono irregular o contém

occ_with_polygons <- sf::st_join(
  occ_sf, 
  study_area[, c(polygon_id_column, "geometry")], 
  join = sf::st_intersects,
  left = FALSE  # Remove pontos que não caem em nenhum polígono
)

# Quantos pontos foram atribuídos?
cat("Pontos dentro de polígonos:", nrow(occ_with_polygons), "\n")
cat("Pontos fora de polígonos (removidos):", nrow(occ_sf) - nrow(occ_with_polygons), "\n\n")

# Criar tabela de espécies por polígono (pares únicos)
species_in_polygons <- occ_with_polygons %>%
  sf::st_drop_geometry() %>%
  dplyr::select(spp, dplyr::all_of(polygon_id_column)) %>%
  dplyr::distinct()

cat("Pares espécie-polígono encontrados:", nrow(species_in_polygons), "\n")
print(species_in_polygons)
cat("\n")

# ==============================================================================
# PASSO 5: Criar o grid regular sobre a área de estudo
# ==============================================================================

# Definir resolução do grid (em graus, se WGS84)
grid_resolution <- 5  # <-- MUDE para a resolução desejada (ex: 1, 2, 5, 10)

cat("Criando grid regular com resolução:", grid_resolution, "graus\n")

# Criar raster cobrindo a extensão do shapefile
bb <- sf::st_bbox(study_area)
grid_raster <- raster::raster(
  xmn = bb["xmin"],
  xmx = bb["xmax"],
  ymn = bb["ymin"],
  ymx = bb["ymax"],
  resolution = c(grid_resolution, grid_resolution),
  crs = sp::CRS("+proj=longlat +datum=WGS84")
)

# Estender um pouco para pegar bordas
grid_raster <- raster::extend(grid_raster, c(1, 1))

# Converter raster para polígonos (cada célula vira um polígono)
grid_polygons <- raster::rasterToPolygons(grid_raster)

# Converter para sf
grid_sf <- sf::st_as_sf(grid_polygons)
grid_sf$grid_id <- seq_len(nrow(grid_sf))

cat("Grid total criado:", nrow(grid_sf), "células\n")

# Manter apenas células que intersectam a área de estudo
grid_keep <- lengths(sf::st_intersects(grid_sf, study_area)) > 0
grid_sf <- grid_sf[grid_keep, ]

cat("Grid após filtro (intersecta shapefile):", nrow(grid_sf), "células\n\n")

# ==============================================================================
# PASSO 6: Mapear grid cells -> polígonos irregulares
# ==============================================================================
# Para cada célula de grid, descobrir quais polígonos ela intersecta.
# Uma célula pode intersectar mais de um polígono (nas bordas).

grid_polygon_intersection <- sf::st_join(
  grid_sf, 
  study_area[, c(polygon_id_column, "geometry")], 
  join = sf::st_intersects
)

# Remover células que não intersectam nenhum polígono
grid_polygon_intersection <- grid_polygon_intersection[!is.na(grid_polygon_intersection[[polygon_id_column]]), ]

cat("Mapeamento grid->polígono criado:", nrow(grid_polygon_intersection), "associações\n")
cat("(Uma célula pode aparecer mais de uma vez se intersecta múltiplos polígonos)\n\n")

# ==============================================================================
# PASSO 7: Construir a matriz de presença/ausência
# ==============================================================================
# LÓGICA: Para cada espécie:
#   1. Encontrar em quais polígonos ela ocorre (do PASSO 4)
#   2. Encontrar TODAS as células de grid que intersectam esses polígonos (do PASSO 6)
#   3. Marcar essas células como presença = 1

species_list <- sort(unique(occ_data$spp))

# Inicializar matriz vazia (linhas = grid_id, colunas = espécies)
pa_matrix <- matrix(
  0, 
  nrow = nrow(grid_sf), 
  ncol = length(species_list),
  dimnames = list(as.character(grid_sf$grid_id), species_list)
)

cat("Preenchendo matriz de presença/ausência...\n")

for (sp in species_list) {
  # 1. Quais polígonos têm essa espécie?
  polys_with_sp <- species_in_polygons[[polygon_id_column]][species_in_polygons$spp == sp]
  polys_with_sp <- unique(polys_with_sp)
  
  if (length(polys_with_sp) > 0) {
    # 2. Quais grid cells intersectam esses polígonos?
    grids_with_sp <- grid_polygon_intersection$grid_id[
      grid_polygon_intersection[[polygon_id_column]] %in% polys_with_sp
    ]
    grids_with_sp <- unique(grids_with_sp)
    
    # 3. Marcar presença = 1
    pa_matrix[as.character(grids_with_sp), sp] <- 1
    
    cat("  ", sp, "-> presente em", length(polys_with_sp), "polígonos ->",
        length(grids_with_sp), "grid cells marcadas\n")
  }
}

cat("\nMatriz completa:", nrow(pa_matrix), "linhas x", ncol(pa_matrix), "colunas\n")

# ==============================================================================
# PASSO 8: Filtrar linhas vazias e adicionar ROOT (para PAE-PCE)
# ==============================================================================

# Remover grid cells sem nenhuma espécie
pa_matrix_filtered <- pa_matrix[rowSums(pa_matrix) > 0, , drop = FALSE]

cat("Após remover linhas vazias:", nrow(pa_matrix_filtered), "grid cells com espécies\n")

# Adicionar linha ROOT (necessária para PAE-PCE)
pa_matrix_final <- rbind(pa_matrix_filtered, ROOT = rep(0, ncol(pa_matrix_filtered)))

cat("Matriz final (com ROOT):", nrow(pa_matrix_final), "x", ncol(pa_matrix_final), "\n\n")

# Visualizar as primeiras linhas
cat("Primeiras linhas da matriz:\n")
print(head(pa_matrix_final, 10))

# ==============================================================================
# PASSO 9: Calcular riqueza por polígono irregular (para Leaflet)
# ==============================================================================
# Contar quantas espécies ocorrem em cada polígono

richness_by_polygon <- species_in_polygons %>%
  dplyr::group_by(.data[[polygon_id_column]]) %>%
  dplyr::summarise(
    n_species = dplyr::n_distinct(spp),
    species_list = paste(sort(unique(spp)), collapse = ", "),
    .groups = "drop"
  )

cat("Riqueza por polígono irregular:\n")
print(as.data.frame(richness_by_polygon))
cat("\n")

# ==============================================================================
# PASSO 10: Visualização
# ==============================================================================

# 10A: Plotar o mapa base com polígonos + grid + pontos
par(mfrow = c(1, 1))
plot(sf::st_geometry(study_area), 
     main = "Occurrence -> Irregular Polygons -> Regular Grids",
     border = "darkgray", lwd = 1.5, col = "lightyellow")

# Adicionar grid cells que têm pelo menos 1 espécie
occupied_grids <- grid_sf[grid_sf$grid_id %in% as.integer(rownames(pa_matrix_filtered)), ]
plot(sf::st_geometry(occupied_grids), add = TRUE, 
     col = rgb(0.2, 0.6, 1, 0.3), border = "blue", lwd = 0.5)

# Adicionar pontos de ocorrência
plot(sf::st_geometry(occ_sf), add = TRUE, col = "red", pch = 16, cex = 1.2)

legend("topright", 
       legend = c("Polígonos irregulares", "Grid cells com presença", "Pontos de ocorrência"),
       fill = c("lightyellow", rgb(0.2, 0.6, 1, 0.3), NA),
       border = c("darkgray", "blue", NA),
       pch = c(NA, NA, 16),
       col = c(NA, NA, "red"),
       cex = 0.8)

# 10B: Plotar uma espécie específica para verificar
sp_to_plot <- species_list[1]  # <-- Mude o índice para ver outra espécie
cat("\nPlotando espécie:", sp_to_plot, "\n")

# Polígonos onde a espécie ocorre
polys_sp <- species_in_polygons[[polygon_id_column]][species_in_polygons$spp == sp_to_plot]
polys_sp_sf <- study_area[study_area[[polygon_id_column]] %in% polys_sp, ]

# Grid cells onde a espécie tem presença
grids_sp <- rownames(pa_matrix_filtered)[pa_matrix_filtered[, sp_to_plot] == 1]
grids_sp_sf <- grid_sf[grid_sf$grid_id %in% as.integer(grids_sp), ]

plot(sf::st_geometry(study_area), 
     main = paste("Presença de", sp_to_plot),
     border = "darkgray", lwd = 1, col = "white")
plot(sf::st_geometry(polys_sp_sf), add = TRUE, col = "lightgreen", border = "darkgreen", lwd = 2)
plot(sf::st_geometry(grids_sp_sf), add = TRUE, col = rgb(1, 0, 0, 0.4), border = "red", lwd = 0.5)
plot(sf::st_geometry(occ_sf[occ_sf$spp == sp_to_plot, ]), add = TRUE, col = "black", pch = 16, cex = 1.5)

legend("topright",
       legend = c("Polígono com espécie", "Grid cells = presença", "Pontos de ocorrência"),
       fill = c("lightgreen", rgb(1, 0, 0, 0.4), NA),
       border = c("darkgreen", "red", NA),
       pch = c(NA, NA, 16),
       col = c(NA, NA, "black"),
       cex = 0.8)

# ==============================================================================
# PASSO 11: Exportar resultados (opcional)
# ==============================================================================

# Salvar matriz como CSV
write.csv(pa_matrix_final, "matriz_presenca_ausencia_grids.csv")
cat("\nMatriz salva em: matriz_presenca_ausencia_grids.csv\n")

# Salvar riqueza por polígono
write.csv(as.data.frame(richness_by_polygon), "riqueza_por_poligono.csv", row.names = FALSE)
cat("Riqueza salva em: riqueza_por_poligono.csv\n")

# ==============================================================================
# FIM DA ROTINA
# ==============================================================================
cat("\n")
cat("================================================================\n")
cat("RESUMO:\n")
cat("================================================================\n")
cat("  Espécies:", length(species_list), "\n")
cat("  Polígonos irregulares:", length(unique(study_area[[polygon_id_column]])), "\n")
cat("  Grid cells totais:", nrow(grid_sf), "\n")
cat("  Grid cells com presença:", nrow(pa_matrix_filtered), "\n")
cat("  Resolução do grid:", grid_resolution, "graus\n")
cat("  Matriz final:", nrow(pa_matrix_final), "x", ncol(pa_matrix_final), "\n")
cat("================================================================\n")
