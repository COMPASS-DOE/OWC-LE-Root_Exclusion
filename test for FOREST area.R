library(terra)
library(dplyr)


#####################################################################
##  Old Woman Creek — floresta inundada vs. terra firme (C-CAP 2006) bacia do estuario
#####################################################################

# install.packages(c("nhdplusTools", "terra", "sf", "dplyr", "tibble"))
library(nhdplusTools)
library(terra)
library(sf)
library(dplyr)

# ---- 1. Download do raster C-CAP -----------------------------------------
url_img  <- "https://chs.coast.noaa.gov/htdata/raster1/landcover/bulkdownload/hires/oh/oh_owc_nerr_2006_ccap_high_res_land_cover.img"
dest_img <- "oh_owc_nerr_2006_ccap_high_res_land_cover.img"

if (!file.exists(dest_img)) {
  download.file(url_img, destfile = dest_img, mode = "wb")
}

ccap <- rast(dest_img)
print(ccap)
print(crs(ccap, describe = TRUE))   # NAD83 / Ohio North (ftUS) esperado
print(res(ccap))                    # ~6.56 ft = 2 m esperado
plot(ccap)

# ---- 2. Delinear a bacia via NLDI (pour point = gauge USGS) --------------
gauge_id <- "USGS-04199155"   # Old Woman Creek at Berlin Rd near Huron OH

basin <- get_nldi_basin(list(featureSource = "nwissite", featureID = gauge_id))

basin_area_km2 <- st_area(basin) |> units::set_units("km^2") |> as.numeric()
cat(sprintf(
  "Área da bacia delineada via NLDI: %.2f km²\n(esperado ~57.2 km² pelo drainage area do gauge USGS)\n",
  basin_area_km2
))

st_write(basin, "owc_watershed_nldi.shp", delete_dsn = TRUE)

watershed <- vect(basin)
watershed <- project(watershed, crs(ccap))

# ---- 3. Checagem de cobertura ---------------------------------------------
watershed_area_km2 <- expanse(watershed, unit = "km") |> sum()
raster_area_km2    <- expanse(as.polygons(ext(ccap), crs = crs(ccap)), unit = "km") |> sum()

cat(sprintf(
  "Bacia (NLDI): %.2f km² | Raster C-CAP OWC 2006 (extent): %.2f km²\n",
  watershed_area_km2, raster_area_km2
))

if (raster_area_km2 < watershed_area_km2) {
  warning(sprintf(
    paste0(
      "O raster (%.2f km²) NÃO cobre a bacia inteira (%.2f km²) — cobre só a ",
      "reserva/entorno imediato. A tabela abaixo vai refletir só a fração ",
      "coberta, não a bacia toda. Pra cobrir o resto, seria preciso mosaicar ",
      "com oh_erie_co_ccap_hr_landcover20190528.img (ano 2019, não 2006)."
    ),
    raster_area_km2, watershed_area_km2
  ))
}

# ---- 4. Crop + mask --------------------------------------------------------
ccap_crop <- crop(ccap, watershed)
ccap_mask <- mask(ccap_crop, watershed)

# ---- 5. Área por classe via expanse() (CRS-aware, sem bug de unidade) -----
class_lookup <- tibble::tibble(
  value = c(0, 1, 2, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23),
  class_name = c(
    "Background", "Unclassified", "Impervious", "Developed Open Space",
    "Cultivated Crops", "Pasture/Hay", "Grassland/Herbaceous",
    "Deciduous Forest", "Evergreen Forest", "Scrub/Shrub",
    "Palustrine Forested Wetland", "Palustrine Scrub/Shrub Wetland",
    "Palustrine Emergent Wetland", "Estuarine Forested Wetland",
    "Estuarine Scrub/Shrub Wetland", "Estuarine Emergent Wetland",
    "Unconsolidated Shore", "Bare Land", "Open Water",
    "Palustrine Aquatic Bed", "Estuarine Aquatic Bed"
  )
)

area_tab <- expanse(ccap_mask, byValue = TRUE, unit = "km") %>%
  as_tibble() %>%
  rename(area_km2 = area) %>%
  left_join(class_lookup, by = "value") %>%
  mutate(area_ha = area_km2 * 100)

print(area_tab)

cat(sprintf(
  "Soma das classes: %.2f km² | Área da bacia (NLDI): %.2f km²\n",
  sum(area_tab$area_km2), watershed_area_km2
))

# ---- 6. Separar floresta inundada vs. terra firme --------------------------
flooded_forest_classes <- c("Palustrine Forested Wetland", "Estuarine Forested Wetland")
upland_forest_classes  <- c("Deciduous Forest", "Evergreen Forest")

summary_tab <- area_tab %>%
  mutate(
    forest_type = case_when(
      class_name %in% flooded_forest_classes ~ "Flooded forest (Palustrine/Estuarine)",
      class_name %in% upland_forest_classes  ~ "Upland forest (Deciduous/Evergreen)",
      TRUE ~ "Other (non-forest)"
    )
  ) %>%
  group_by(forest_type) %>%
  summarise(
    area_ha  = sum(area_ha, na.rm = TRUE),
    area_km2 = sum(area_km2, na.rm = TRUE),
    .groups  = "drop"
  )

print(summary_tab)
############################
library(terra)
library(tidyterra)
library(ggplot2)

# ---- 1. Reclassificar (matriz corrigida) ------------------------------------
flooded_vals <- c(13, 16)    # Palustrine / Estuarine Forested Wetland
upland_vals  <- c(9, 10)     # Deciduous / Evergreen Forest

rcl_matrix <- rbind(
  cbind(flooded_vals, 1),
  cbind(upland_vals, 2)
)

print(rcl_matrix)   # confere ANTES de rodar: precisa mostrar
#      13 1
#      16 1
#       9 2
#      10 2

ccap_forest <- classify(ccap_mask, rcl = rcl_matrix, others = NA)

# ---- Checagem antes de seguir -------------------------------------------------
print(unique(values(ccap_forest)))   # agora deve mostrar só NA, 1, 2 — nada mais

coltab(ccap_forest) <- NULL
ccap_forest <- as.factor(ccap_forest)
levels(ccap_forest) <- data.frame(
  value = c(1, 2),
  forest_type = c(
    "Flooded forest (Palustrine)",
    "Upland forest (Deciduous/Evergreen)"
  )
)

print(levels(ccap_forest))
plot(ccap_forest)   # checkpoint em base R

# ---- 2. Mapa ggplot ------------------------------------------------------------
p_forest_map <- ggplot() +
  geom_spatvector(data = watershed, fill = "grey95", color = "grey40", linewidth = 0.6) +
  geom_spatraster(data = ccap_forest, aes(fill = forest_type)) +
  scale_fill_manual(
    values = c(
      "Flooded forest (Palustrine)"         = "#3182BD",
      "Upland forest (Deciduous/Evergreen)" = "#31A354"
    ),
    na.translate = FALSE,
    na.value = "transparent",
    name = "Forest type"
  ) +
  labs(
    title = "Old Woman Creek watershed — forest cover by hydrologic condition",
    subtitle = "Flooded forest: 4.44 km² (7.7%)  |  Upland forest: 6.15 km² (10.6%)  |  Total: 10.59 km² (18.3%)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_blank()
  )

p_forest_map

#ggsave("owc_flooded_vs_upland_forest.png", p_forest_map, width = 8, height = 6, dpi = 300)


library(sf)

# ---- Ponto da reserva (WGS84) e reprojeção pro CRS do raster ----------------
owc_reserve_wgs84 <- st_sfc(
  #st_point(c(-82.5116, 41.3761)), #visitor center
  st_point(c(-82.50635, 41.37576)),#Transition
  crs = 4326
)

owc_reserve_proj <- st_transform(owc_reserve_wgs84, crs = crs(ccap_forest))
owc_reserve_coords <- st_coordinates(owc_reserve_proj)

# ---- Mapa com a estrela ------------------------------------------------------
p_forest_map <- ggplot() +
  geom_spatvector(data = watershed, fill = "grey95", color = "grey40", linewidth = 0.6) +
  geom_spatraster(data = ccap_forest, aes(fill = forest_type)) +
  scale_fill_manual(
    values = c(
      "Flooded forest (Palustrine)"         = "#3182BD",
      "Upland forest (Deciduous/Evergreen)" = "#31A354"
    ),
    na.translate = FALSE,
    na.value = "transparent",
    name = "Forest type"
  ) +
  geom_point(
    aes(x = owc_reserve_coords[1, "X"], y = owc_reserve_coords[1, "Y"]),
    shape = 8, color = "red", size = 4, stroke = 1.3
  ) +
  annotate(
    "text",
    x = owc_reserve_coords[1, "X"], y = owc_reserve_coords[1, "Y"],
    label = "Old Woman Creek NERR",
    color = "red", size = 3, fontface = "bold",
    vjust = -1.2
  ) +
  labs(
    title = "Old Woman Creek watershed — forest cover by hydrologic condition",
    subtitle = "Flooded forest: 4.44 km² (7.7%)  |  Upland forest: 6.15 km² (10.6%)  |  Total: 10.59 km² (18.3%)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_blank()
  )

p_forest_map


#####################################
######################################             BACIA da TRANSICAO

library(nhdplusTools)
library(sf)
library(terra)

# ---- 1. Ponto de interesse (WGS84) -----------------------------------------
new_point <- st_sfc(
  st_point(c(-82.50635, 41.37576)),   # ordem: (lon, lat)
  crs = 4326
)

# ---- 2. Achar o COMID do trecho de rio mais próximo -------------------------
new_comid <- discover_nhdplus_id(new_point)
print(new_comid)   # confere que retornou um valor válido antes de seguir

# ---- 3. Delinear a bacia a partir desse COMID --------------------------------
basin_new <- get_nldi_basin(
  nldi_feature = list(featureSource = "comid", featureID = as.character(new_comid))
)

basin_new_area_km2 <- st_area(basin_new) |> units::set_units("km^2") |> as.numeric()
cat(sprintf("Área da nova bacia (via ponto): %.2f km²\n", basin_new_area_km2))

# ---- 4. Confere se o ponto agora cai dentro do polígono ----------------------
inside_check <- st_within(new_point, basin_new)
print(inside_check)   # deve retornar o índice 1 na lista, não list()

# ---- 5. Salvar e substituir a bacia usada no resto do pipeline ---------------
st_write(basin_new, "owc_watershed_newpoint.shp", delete_dsn = TRUE)

watershed <- vect(basin_new)
watershed <- project(watershed, crs(ccap))   # crs(ccap) já deve estar carregado da sessão anterior

# ---- 6. Recalcula cobertura raster x bacia nova -------------------------------
watershed_area_km2 <- expanse(watershed, unit = "km") |> sum()
raster_area_km2    <- expanse(as.polygons(ext(ccap), crs = crs(ccap)), unit = "km") |> sum()

cat(sprintf(
  "Bacia (novo ponto): %.2f km² | Raster C-CAP OWC 2006 (extent): %.2f km²\n",
  watershed_area_km2, raster_area_km2
))

library(dplyr)
library(tidyterra)
library(ggplot2)

# ---- 1. Crop + mask com a bacia nova ------------------------------------------
ccap_crop <- crop(ccap, watershed)
ccap_mask <- mask(ccap_crop, watershed)

# ---- 2. Área por classe --------------------------------------------------------
class_lookup <- tibble::tibble(
  value = c(0, 1, 2, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23),
  class_name = c(
    "Background", "Unclassified", "Impervious", "Developed Open Space",
    "Cultivated Crops", "Pasture/Hay", "Grassland/Herbaceous",
    "Deciduous Forest", "Evergreen Forest", "Scrub/Shrub",
    "Palustrine Forested Wetland", "Palustrine Scrub/Shrub Wetland",
    "Palustrine Emergent Wetland", "Estuarine Forested Wetland",
    "Estuarine Scrub/Shrub Wetland", "Estuarine Emergent Wetland",
    "Unconsolidated Shore", "Bare Land", "Open Water",
    "Palustrine Aquatic Bed", "Estuarine Aquatic Bed"
  )
)

area_tab <- expanse(ccap_mask, byValue = TRUE, unit = "km") %>%
  as_tibble() %>%
  rename(area_km2 = area) %>%
  left_join(class_lookup, by = "value") %>%
  mutate(area_ha = area_km2 * 100)

print(area_tab)
cat(sprintf("Soma das classes: %.2f km² | Área da bacia (novo ponto): %.2f km²\n",
            sum(area_tab$area_km2), watershed_area_km2))

# ---- 3. Separar floresta inundada vs. terra firme ------------------------------
flooded_forest_classes <- c("Palustrine Forested Wetland", "Estuarine Forested Wetland")
upland_forest_classes  <- c("Deciduous Forest", "Evergreen Forest")

summary_tab <- area_tab %>%
  mutate(
    forest_type = case_when(
      class_name %in% flooded_forest_classes ~ "Flooded forest (Palustrine/Estuarine)",
      class_name %in% upland_forest_classes  ~ "Upland forest (Deciduous/Evergreen)",
      TRUE ~ "Other (non-forest)"
    )
  ) %>%
  group_by(forest_type) %>%
  summarise(area_ha = sum(area_ha, na.rm = TRUE),
            area_km2 = sum(area_km2, na.rm = TRUE),
            .groups = "drop")

print(summary_tab)   # esses números vão substituir os do subtítulo do mapa abaixo

# ---- 4. Reclassificar (matriz corrigida) ---------------------------------------
flooded_vals <- c(13, 16)
upland_vals  <- c(9, 10)

rcl_matrix <- rbind(
  cbind(flooded_vals, 1),
  cbind(upland_vals, 2)
)

ccap_forest <- classify(ccap_mask, rcl = rcl_matrix, others = NA)
print(unique(values(ccap_forest)))   # checkpoint: só NA, 1, 2

coltab(ccap_forest) <- NULL
ccap_forest <- as.factor(ccap_forest)
levels(ccap_forest) <- data.frame(
  value = c(1, 2),
  forest_type = c(
    "Flooded forest (Palustrine)",
    "Upland forest (Deciduous/Evergreen)"
  )
)

# ---- 5. Coordenadas do ponto no CRS do raster (já calculado antes) -------------
new_point_proj <- st_transform(new_point, crs = crs(ccap_forest))
new_point_coords <- st_coordinates(new_point_proj)

# ---- 6. Mapa final ---------------------------------------------------------------
p_forest_map <- ggplot() +
  geom_spatvector(data = watershed, fill = "grey95", color = "grey40", linewidth = 0.6) +
  geom_spatraster(data = ccap_forest, aes(fill = forest_type)) +
  scale_fill_manual(
    values = c(
      "Flooded forest (Palustrine)"         = "#3182BD",
      "Upland forest (Deciduous/Evergreen)" = "#31A354"
    ),
    na.translate = FALSE,
    na.value = "transparent",
    name = "Forest type"
  ) +
  geom_point(
    aes(x = new_point_coords[1, "X"], y = new_point_coords[1, "Y"]),
    shape = 8, color = "red", size = 4, stroke = 1.3
  ) +
  annotate(
    "text",
    x = new_point_coords[1, "X"], y = new_point_coords[1, "Y"],
    label = "Old Woman Creek NERR",
    color = "red", size = 3, fontface = "bold",
    vjust = -1.2
  ) +
  labs(
    title = "Old Woman Creek watershed — forest cover by hydrologic condition",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid = element_blank(),
    axis.text = element_blank()
  )

p_forest_map


# ---- Verificação: a bacia antiga é subconjunto da bacia nova? -----------------
old_gauge <- list(featureSource = "nwissite", featureID = "USGS-04199155")
basin_old <- get_nldi_basin(nldi_feature = old_gauge)

# reprojeta os dois pro mesmo CRS antes de comparar
basin_old_proj <- st_transform(basin_old, crs = st_crs(basin_new))

# a antiga deveria estar quase inteiramente DENTRO da nova
overlap_area_km2 <- st_intersection(basin_old_proj, basin_new) |>
  st_area() |> units::set_units("km^2") |> as.numeric()

old_area_km2 <- st_area(basin_old_proj) |> units::set_units("km^2") |> as.numeric()
new_area_km2 <- st_area(basin_new) |> units::set_units("km^2") |> as.numeric()

cat(sprintf(
  "Bacia antiga (gauge): %.2f km²\nBacia nova (seu ponto): %.2f km²\nSobreposição: %.2f km² (%.1f%% da bacia antiga está contida na nova)\nÁrea extra na nova (não coberta pela antiga): %.2f km²\n",
  old_area_km2, new_area_km2, overlap_area_km2,
  100 * overlap_area_km2 / old_area_km2,
  new_area_km2 - overlap_area_km2
))

# ---- Mapa comparativo dos dois contornos, pra visualizar a diferença -----------
ggplot() +
  geom_sf(data = basin_new, fill = "#3182BD", alpha = 0.3, color = "#3182BD", linewidth = 0.8) +
  geom_sf(data = basin_old_proj, fill = "#DE2D26", alpha = 0.3, color = "#DE2D26", linewidth = 0.8) +
  geom_sf(data = new_point, color = "black", size = 3, shape = 8) +
  labs(title = "Bacia antiga (vermelho, gauge USGS) vs. nova (azul, seu ponto)") +
  theme_minimal()



# ---- Baixa a rede de drenagem (flowlines) dentro da bacia nova ----------------
flowlines <- get_nhdplus(AOI = basin_new, realization = "flowline")

cat(sprintf("Número de trechos de rio (flowlines) dentro da bacia: %d\n", nrow(flowlines)))

# ---- Mapa: bacia + rede de drenagem completa -----------------------------------
ggplot() +
  geom_sf(data = basin_new, fill = "grey90", color = "black", linewidth = 0.6) +
  geom_sf(data = flowlines, color = "#3182BD", linewidth = 0.5) +
  geom_sf(data = new_point, color = "red", size = 3, shape = 8) +
  labs(
    title = "Rede de drenagem dentro da bacia delineada",
    subtitle = sprintf("%d trechos de rio/afluentes", nrow(flowlines))
  ) +
  theme_minimal()