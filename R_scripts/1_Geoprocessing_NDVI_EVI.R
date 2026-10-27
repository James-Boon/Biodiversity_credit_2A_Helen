# =========
# PACKAGES
# =========

libs <- c("dplyr", "sf", "terra", "corrplot")

installed_libs <- libs %in% rownames(installed.packages())

if (any(installed_libs == F)) {install.packages(
  libs[!installed_libs])}

invisible(lapply(libs, library, character.only = T))
rm(list = ls())

# ======================
# LOAD PLANETSCOPE DATA
# ======================

sat_image_knepp <- rast("./Data/Remote_sensing_data/Knepp_24/20240913_111340_90_24ed_3B_AnalyticMS_SR_8b_clip.tif")
sat_image_attenborough <- rast("./Data/Remote_sensing_data/Attenborough_24/20240811_111515_48_24bd_3B_AnalyticMS_SR_8b_clip.tif")

# ============
# INSPECT DATA
# ============
sat_image_knepp
names(sat_image_knepp)
crs(sat_image_knepp)

sat_image_attenborough
names(sat_image_attenborough)
crs(sat_image_attenborough)

# =================
# CALCULATE INDICES
# =================
source("./R_scripts/Functions/calc_indices.R")
knepp_stack <- calc_indices(sat_image_knepp)
attenborough_stack <- calc_indices(sat_image_attenborough)

# =======================
# SHAPEFILE GEOPROCESSING
# =======================

# KNEPP
knepp_dir <- "./Data/Knepp_edited_shapefiles/"

knepp_files <- list.files(
  path = knepp_dir,
  pattern = "\\.shp$",
  full.names = TRUE
)
knepp_files

knepp_category_lookup <- c(
  "Grasslands"  = "Grassland",
  "Ponds_Lakes" = "Lakes",
  "Scrub"       = "Scrub",
  "Woodlands"   = "Woodland"
)

habitats_knepp <- knepp_files |>
  lapply(function(f) {
    x <- st_read(f, quiet = TRUE)
    
    file_stub <- tools::file_path_sans_ext(basename(f))
    x$Category <- knepp_category_lookup[file_stub]
    
    x
  }) |>
  bind_rows()

print(habitats_knepp)
names(habitats_knepp)
table(habitats_knepp$Habitat)

clean_condition <- function(x) {
  gsub("^[0-9]+\\.\\s*", "", trimws(x))
}

habitats_knepp <- habitats_knepp |>
  transmute(
    poly_uid = dplyr::row_number(),
    ID = as.character(ID),
    Habitat = as.character(Habitat),
    Category = as.character(Category),
    Area = as.numeric(Area),
    Condition = clean_condition(as.character(Condition)),
    geometry = geometry
  )

st_crs(habitats_knepp)
# ATTENBOROUGH
Att_dir <- "./Data/Attenborough_edited_shapefiles/"

Att_files <- list.files(
  path = Att_dir,
  pattern = "\\.shp$",
  full.names = TRUE
)
Att_files

habitats_att <- Att_files |>
  lapply(function(f) {
    x <- st_read(f, quiet = TRUE)
    
    # convert all non-geometry columns to character
    x <- x |>
      mutate(across(-geometry, as.character))
    
    x
  }) |>
  bind_rows()

print(habitats_att)
names(habitats_att)
table(habitats_att$Habitat)

habitats_att <- habitats_att |>
  transmute(
    poly_uid = dplyr::row_number(),
    ID = as.character(polygon_nu),
    Habitat = as.character(UK_hab),
    Category = as.character(Category),
    Area = as.numeric(Area),
    Condition = clean_condition(as.character(Condition)),
    geometry = geometry
  )

print(habitats_att)
names(habitats_att)

table(habitats_att$Habitat, useNA = "ifany")
table(habitats_att$Condition, useNA = "ifany")
table(habitats_att$Category, useNA = "ifany")

st_crs(habitats_att)
habitats_att <- st_transform(habitats_att, 27700)

# Write polygons
dir.create("./Data/Processed_data", showWarnings = FALSE)
dir.create("./Data/Processed_data/Knepp/", showWarnings = FALSE)
dir.create("./Data/Processed_data/Attenborough/", showWarnings = FALSE)

st_write(
  habitats_knepp,
  "./Data/Processed_data/Knepp/Polygons_habitats_knepp.gpkg",
  delete_dsn = TRUE
)

st_write(
  habitats_att,
  "./Data/Processed_data/Attenborough/Polygons_habitats_att.gpkg",
  delete_dsn = TRUE
)

# ====================
# REPROJECTION TO BNG
# ====================

knepp_stack_bng <- project(
  knepp_stack,
  "EPSG:27700",
  method = "bilinear"
)

attenborough_stack_bng <- project(
  attenborough_stack,
  "EPSG:27700",
  method = "bilinear"
)

knepp_stack_bng
crs(knepp_stack_bng)

knepp_stack_bng
crs(attenborough_stack_bng)

# ====================
# WRITE RASTER STACKS
# ====================

writeRaster(
  knepp_stack_bng,
  "./Data/Processed_data/Knepp/Knepp_stack_bng.tif",
  overwrite = TRUE
)

writeRaster(
  attenborough_stack_bng,
  "./Data/Processed_data/Attenborough/Attenborough_stack_bng.tif",
  overwrite = TRUE
)

# ==================
# POLYGON EXTRACTION
# ==================

# Knepp
habitats_knepp_vect <- vect(habitats_knepp)
knepp_extract <- terra::extract(knepp_stack_bng, habitats_knepp_vect)
head(knepp_extract)

knepp_lookup <- as.data.frame(habitats_knepp_vect) |>
  mutate(extract_row = row_number())

knepp_extract <- terra::extract(knepp_stack_bng, habitats_knepp_vect) |>
  rename(extract_row = ID) |>
  left_join(knepp_lookup, by = "extract_row")

# Attenborough
habitats_att_vect <- vect(habitats_att)
att_extract <- terra::extract(attenborough_stack_bng, habitats_att_vect)
head(att_extract)

att_lookup <- as.data.frame(habitats_att_vect) |>
  mutate(extract_row = row_number())

att_extract <- terra::extract(attenborough_stack_bng, habitats_att_vect) |>
  rename(extract_row = ID) |>
  left_join(att_lookup, by = "extract_row")

# =====================
# SUMMARISE IN POLYGONS
# =====================

att_polygon_summary <- att_extract |>
  group_by(poly_uid) |>
  summarise(
    ID        = first(ID),
    Habitat   = first(Habitat),
    Category  = first(Category),
    Area      = first(Area),
    Condition = first(Condition),
    across(
      c(coastal_blue, blue, green_i, green, yellow, red, rededge, nir,
        ndvi, evi, ndwi),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(Site = "Attenborough")
write.csv(att_polygon_summary, "./Data/Processed_data/Attenborough/Att_habitat_summary.csv",row.names = F)

knepp_polygon_summary <- knepp_extract |>
  group_by(poly_uid) |>
  summarise(
    ID        = first(ID),
    Habitat   = first(Habitat),
    Category  = first(Category),
    Area      = first(Area),
    Condition = first(Condition),
    across(
      c(coastal_blue, blue, green_i, green, yellow, red, rededge, nir,
        ndvi, evi, ndwi),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(Site = "Knepp")
write.csv(knepp_polygon_summary, "./Data/Processed_data/Knepp/Knepp_habitat_summary.csv",row.names = F)

combined_df <- bind_rows(
  att_polygon_summary,
  knepp_polygon_summary
)

clean_df <- combined_df |>
  select(
    Site, poly_uid, ID, Habitat, Category, Condition, Area,
    ndvi_mean, ndvi_sd,
    evi_mean, evi_sd,
    ndwi_mean, ndwi_sd
  )

write.csv(
  clean_df,
  "./Data/Processed_data/combined_df.csv",
  row.names = FALSE
)

# ====================
# CORRELATION ANALYSIS
# ====================

num_vars <- clean_df |>
  select(ndvi_mean, ndvi_sd, evi_mean, evi_sd, ndwi_mean, ndwi_sd)

cor_matrix <- cor(num_vars, use = "complete.obs")

dir.create("./Plots/Correlation_plot", showWarnings = FALSE)

png(
  filename = "./Plots/Correlation_plot/correlation_matrix.png",
  width = 2000,
  height = 2000,
  res = 300
)

corrplot(
  cor_matrix,
  method = "square",
  type = "lower",
  addCoef.col = "black",
  tl.cex = 0.9,
  number.cex = 0.7
)

dev.off()
