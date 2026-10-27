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
