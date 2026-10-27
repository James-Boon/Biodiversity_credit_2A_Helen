# =========
# CLEAN ENV
# =========
rm(list = ls())

# =========
# PACKAGES
# =========
libs <- c("dplyr", "sf", "terra", "corrplot")

installed_libs <- libs %in% rownames(installed.packages())
if (any(!installed_libs)) {
  install.packages(libs[!installed_libs])
}

invisible(lapply(libs, library, character.only = TRUE))

source("./R_scripts/Functions/summarise_Li_raster_to_polygons.R")
# =================
# LOAD LiDAR DATA
# =================
lidar_dir <- "./Data/Remote_sensing_data/Knepp_24/LiDAR/"

lidar_files <- c(
  canopy_cover = "Canopy_Cover.tif",
  lai          = "LAI.tif",
  maxch        = "Maxchm.tif",
  meanch       = "mean_ch.tif",
  mean_lai     = "mean_LAI.tif",
  rugosity     = "Rugosity.tif",
  vch          = "VCH.tif"
)

# Read LiDAR rasters
lidar_rasters <- lapply(names(lidar_files), function(nm) {
  r <- rast(file.path(lidar_dir, lidar_files[[nm]]))
  names(r) <- nm
  crs(r) <- "EPSG:27700"
  r
})
names(lidar_rasters) <- names(lidar_files)

# ===================
# EXTRACT TO POLYGONS
# ===================
summarise_Li_raster_to_polygons <- function(r, polygons_sf, poly_id = "ID", var_name = names(r)) {
  polygons_proj <- st_transform(polygons_sf, crs(r))
  polygons_vect <- vect(polygons_proj)
  
  lookup <- as.data.frame(polygons_vect) |>
    mutate(extract_row = dplyr::row_number())
  
  ext_df <- terra::extract(r, polygons_vect) |>
    rename(extract_row = ID) |>
    left_join(lookup, by = "extract_row")
  
  mean_name <- paste0(var_name, "_mean")
  sd_name   <- paste0(var_name, "_sd")
  
  ext_df |>
    group_by(.data[[poly_id]]) |>
    summarise(
      !!mean_name := mean(.data[[var_name]], na.rm = TRUE),
      !!sd_name   := sd(.data[[var_name]], na.rm = TRUE),
      .groups = "drop"
    )
}



habitats_knepp <- st_read("./Data/Processed_data/Knepp/Polygons_habitats_knepp.gpkg", quiet = TRUE)
knepp_spectral_summary <- read.csv("./Data/Processed_data/Knepp/Knepp_habitat_summary.csv")

lidar_summaries <- lapply(names(lidar_rasters), function(nm) {
  summarise_Li_raster_to_polygons(
    r = lidar_rasters[[nm]],
    polygons_sf = habitats_knepp,
    poly_id = "poly_uid",
    var_name = nm
  )
})

knepp_lidar_summary <- Reduce(function(x, y) left_join(x, y, by = "poly_uid"), lidar_summaries)


knepp_full_summary <- left_join(
  knepp_spectral_summary,
  knepp_lidar_summary,
  by = "poly_uid"
)

Full_df <- knepp_full_summary |>
  rename(
    # Spectral
    coastal_blue_mean = coastal_blue_mean,
    coastal_blue_sd   = coastal_blue_sd,
    
    green1_mean  = green_i_mean,
    green1_sd    = green_i_sd,
    
    # LiDAR
    Canopy_cover_mean   = canopy_cover_mean,
    Canopy_cover_sd     = canopy_cover_sd,
    
    Canopy_height_max_mean = maxch_mean,
    Canopy_height_max_sd   = maxch_sd,
    
    Canopy_height_mean = meanch_mean,
    Canopy_height_sd   = meanch_sd,
    
    Canopy_height_var_mean  = vch_mean,
    Canopy_height_var_sd    = vch_sd
  ) %>%
  dplyr::select(poly_uid,Site, Habitat, Category, Area, Condition, ndvi_mean, ndvi_sd, evi_mean, evi_sd,
                Canopy_cover_mean, Canopy_cover_sd, lai_mean, lai_sd, Canopy_height_max_mean, Canopy_height_max_sd,
                Canopy_height_mean, Canopy_height_sd, rugosity_mean, rugosity_sd, Canopy_height_var_mean, Canopy_height_var_sd)

write.csv(Full_df, "./Data/Processed_data/Knepp/Full_df_Knepp_LiDAR.csv", row.names = F)

# ==================
# CORRELATION MATRIX
# ==================

num_vars <- knepp_full_summary |>
  select(
    ndvi_mean, ndvi_sd,
    evi_mean, evi_sd,
    ndwi_mean, ndwi_sd,
    canopy_cover_mean, canopy_cover_sd,
    lai_mean, lai_sd,
    maxch_mean, maxch_sd,
    meanch_mean, meanch_sd,
    rugosity_mean, rugosity_sd,
    vch_mean, vch_sd
  )

nice_names <- c(
  ndvi_mean = "NDVI (mean)",
  ndvi_sd   = "NDVI (SD)",
  evi_mean  = "EVI (mean)",
  evi_sd    = "EVI (SD)",
  ndwi_mean = "NDWI (mean)",
  ndwi_sd   = "NDWI (SD)",
  
  canopy_cover_mean = "Canopy cover (mean)",
  canopy_cover_sd   = "Canopy cover (SD)",
  
  lai_mean = "Leaf Area Index (mean)",
  lai_sd   = "Leaf Area Index (SD)",
  
  maxch_mean = "Maximum canopy height (mean)",
  maxch_sd   = "Maximum canopy height (SD)",
  
  meanch_mean = "Mean canopy height (mean)",
  meanch_sd   = "Mean canopy height (SD)",
  
  vch_mean = "Canopy height variability (mean)",
  vch_sd   = "Canopy height variability (SD)",
  
  rugosity_mean = "Canopy rugosity (mean)",
  rugosity_sd   = "Canopy rugosity (SD)"
)

cor_matrix <- cor(num_vars, use = "complete.obs")
# Copy matrix so you don’t overwrite original
cor_matrix_named <- cor_matrix

# Replace row/column names
colnames(cor_matrix_named) <- nice_names[colnames(cor_matrix)]
rownames(cor_matrix_named) <- nice_names[rownames(cor_matrix)]
corrplot(
  cor_matrix_named,
  method = "square",
  type = "lower",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.8,
  cl.cex = 0.8
)

png(
  "./Plots/Correlation_plot/knepp_Li_correlation_matrix_.png",
  width = 4000,   # increase size
  height = 4000,
  res = 400       # higher DPI
)

corrplot(
  cor_matrix_named,
  method = "square",
  type = "lower",
  addCoef.col = "black",
  order = "hclust",   # improves layout
  
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.2,       # 🔥 increase label size
  
  number.cex = 0.6,   # 🔥 increase correlation numbers
  
  cl.cex = 1.2        # color legend size
)

dev.off()
