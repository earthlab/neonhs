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
#' @examples
#' library(raster)
#' path_to_file <- system.file('extdata', 'ex.h5', package = 'neonaop')
#'
#' # read the full spatial extent of a file
#' r <- hs_read(path_to_file, bands = 1:4)
#'
#' # read a subset of a file based on a cropping extent
#' extent_to_read <- extent(c(257000, 257020, 4111980, 4112100))
#' subset_r <- hs_read(path_to_file, bands = 1:4, crop = extent_to_read)
#'
#' @export
hs_read <- function(filename, bands, crop = FALSE){
  stopifnot(all(bands > 0))
  h5_extent <- hs_extent(filename)
  use_h5_extent <- isFALSE(crop)
  if (use_h5_extent) {
    out_extent <- h5_extent
    dims <- hs_dims(filename)
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


#' Get the dimensions of a hyperspectral reflectance HDF5 file
#' 
#' The `hs_dims` function returns the dimensions of reflectance data contained
#' within an HDF5 file for NEON's L3 hyperspectral reflectance data. 
#' In most cases, these dimensions will be 426 X 1000 X 1000: 426 bands and 
#' images that are 1000m by 1000m in their spatial extent at 1 meter resolution.
#' 
#' @param filename Path to an .h5 file containing hyperspectral data (char)
#' 
#' @return an integer vector of length 3 containing the number of bands, 
#' number of x pixels, and number of y pixels. 
#' 
#' @examples
#' path_to_file <- system.file('extdata', 'ex.h5', package = 'neonaop')
#' hs_dims(path_to_file)
#' 
#' @export
hs_dims <- function(filename) {
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]
  dims <- reflectance$dims
  file_h5$close_all()
  dims
}


#' Get the spatial extent of a hyperspectral image
#' 
#' @param filename Path to an .h5 file containing hyperspectral data (char)
#' 
#' @return a raster::extent object that contains the min and max x and y 
#' coordinates of a hyperspectral image
#' 
#' @examples
#' path_to_file <- system.file('extdata', 'ex.h5', package = 'neonaop')
#' hs_extent(path_to_file)
#' 
#' @export
hs_extent <- function(filename){
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  path <- paste0(site, '/Reflectance/Metadata/Coordinate_System/Map_Info')
  map_info <- unlist(strsplit(file_h5[[path]]$read(), ','))
  file_h5$close_all()
  dims <- hs_dims(filename)
  xy_resolution <- as.numeric(c(map_info[2], map_info[3]))
  xmin <- as.numeric(map_info[4])
  xmax <- xmin + dims[2] * xy_resolution[1]
  ymax <- as.numeric(map_info[5])
  ymin <- ymax - dims[3] * xy_resolution[2]
  raster::extent(xmin, xmax, ymin, ymax)
}


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


calculate_index_extent <- function(clip_extent, h5_extent){
  dummy_raster <- raster::raster(h5_extent, resolution = 1)
  cells <- raster::cellsFromExtent(dummy_raster, extent = clip_extent)
  colrange <- range(raster::colFromCell(dummy_raster, cells)) # x
  rowrange <- range(raster::rowFromCell(dummy_raster, cells)) # y
  list(xmin_idx = min(colrange), xmax_idx = max(colrange),
       ymin_idx = min(rowrange), ymax_idx = max(rowrange))
}
