calc_indices <- function(img) {
  # Check required bands exist
  required_bands <- c("blue", "green", "red", "nir")
  missing_bands <- setdiff(required_bands, names(img))
  
  if (length(missing_bands) > 0) {
    stop("Missing required bands: ", paste(missing_bands, collapse = ", "))
  }
  
  # Scale reflectance if image is stored as integer SR values
  blue  <- img[["blue"]]  / 10000
  green <- img[["green"]] / 10000
  red   <- img[["red"]]   / 10000
  nir   <- img[["nir"]]   / 10000
  
  # Indices
  ndvi <- (nir - red) / (nir + red)
  names(ndvi) <- "ndvi"
  
  evi <- 2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))
  names(evi) <- "evi"
  
  ndwi <- (green - nir) / (green + nir)
  names(ndwi) <- "ndwi"
  
  # Return original image + indices
  c(img, ndvi, evi, ndwi)
}