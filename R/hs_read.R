#' Read hyperspectral imagery from the NEON AOP
#'
#' The `hs_read` function reads hyperspectral imagery from NEON's Aerial
#' Observation Platform, and returns a Raster* object.
#'
#' You can learn more the NEON AOP hyperspectral imagery at:
#' https://www.neonscience.org/data-collection/airborne-remote-sensing
#'
#' @param filename Path to an .h5 file containing L3 hyperspectral data (char)
#' @param bands Indices of bands to read (integer)
#' @param crop Optional extent object to use in cropping
#' @return Raster* object containing hyperspectral data
#'
#' @export
hs_read <- function(filename, bands, crop = FALSE){
  stopifnot(all(bands > 0))
  h5_extent <- hs_extent(filename)
  use_h5_extent <- isFALSE(crop)
  if (use_h5_extent) {
    out_extent <- h5_extent
    dims <- get_dims(filename)
    index <- list(bands, seq_len(dims[2]), seq_len(dims[3]))
  } else {
    index_extent <- calculate_index_extent(clip_extent = crop,
                                           h5_extent = h5_extent)

    index <- list(bands,
                  index_extent$xmin_idx:(index_extent$xmax_idx),
                  index_extent$ymin_idx:(index_extent$ymax_idx))
    out_extent <- raster::intersect(h5_extent, crop)
  }
  r <- read_hs_values(filename, index)
  raster::extent(r) <- out_extent
  r
}

# library(raster)
# library(viridis)
# r <- hs_read('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5',
#              bands = 1:3)
# r
# plot(r, col = cividis(100))
# ce <- raster::extent(c(257500, 259000, 4111500, 4112100))
# cr <- hs_read('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5',
#               bands = c(1, 150, 300),
#               crop = ce)
# plot(cr, add = TRUE, col = viridis(100))



get_dims <- function(filename) {
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]
  dims <- reflectance$dims
  file_h5$close_all()
  dims
}

# get_dims('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5')



hs_extent <- function(filename){
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  path <- paste0(site, '/Reflectance/Metadata/Coordinate_System/Map_Info')
  map_info <- unlist(strsplit(file_h5[[path]]$read(), ','))
  file_h5$close_all()
  dims <- get_dims(filename)
  xy_resolution <- as.numeric(c(map_info[2], map_info[3]))
  xmin <- as.numeric(map_info[4])
  xmax <- xmin + dims[2] * xy_resolution[1]
  ymax <- as.numeric(map_info[5])
  ymin <- ymax - dims[3] * xy_resolution[2]
  raster::extent(xmin, xmax, ymin, ymax)
}

# e <- hs_extent('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5')


get_epsg <- function(filename) {
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  epsg_path <- paste0(site, '/Reflectance/Metadata/Coordinate_System/EPSG Code')
  epsg_code <- file_h5[[epsg_path]]$read()
  file_h5$close_all()
  epsg_code
}


get_wavelength <- function(filename, bands) {
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  w <- file_h5[[paste0(site, '/Reflectance/Metadata/Spectral_Data/Wavelength')]]
  wavelength_nm <- w$read()[bands]
  file_h5$close_all()
  band_indices <- paste0('band', bands)
  wavelength_names <- paste0(round(wavelength_nm), 'nm')
  layer_names <- paste(band_indices, wavelength_names, sep = '_')
  layer_names
}

# get_wavelength('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5',
#                bands = 1:10)


read_hs_values <- function(filename, index){
  bands <- index[[1]]
  stopifnot(all(bands > 0))
  n_band <- length(bands)

  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name

  reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]

  arrays <- lapply(bands, function(x) {
    raster::raster(reflectance$read(args = list(x, index[[2]], index[[3]])))
  })

  r <- raster::t(raster::stack(arrays))

  na_value <- reflectance$attr_open('Data_Ignore_Value')$read()
  r[r == na_value] <- NA
  scale_factor <- reflectance$attr_open('Scale_Factor')$read()
  r <- r / scale_factor
  file_h5$close_all()

  raster::projection(r) <- paste0("+init=epsg:", get_epsg(filename))
  names(r) <- get_wavelength(filename, bands)
  r
}

# library(raster)
# library(viridis)
# b <- read_hs_values('data/NEON_D17_SJER_DP3_257000_4111000_reflectance.h5',
#              list(c(1, 30, 50), 1:100, 1:100))
# plot(b, col = viridis(100))


calculate_index_extent <- function(clip_extent, h5_extent){
  dummy_raster <- raster::raster(h5_extent, resolution = 1)
  cells <- raster::cellsFromExtent(dummy_raster, extent = clip_extent)
  colrange <- range(raster::colFromCell(dummy_raster, cells)) # x
  rowrange <- range(raster::rowFromCell(dummy_raster, cells)) # y
  list(xmin_idx = min(colrange), xmax_idx = max(colrange),
       ymin_idx = min(rowrange), ymax_idx = max(rowrange))
}
