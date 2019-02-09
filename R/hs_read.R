#' Read hyperspectral imagery from the NEON AOP
#'
#' The `hs_read` function reads hyperspectral imagery from NEON's Airborne
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
  bands <- sort(bands)
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

#' Efficiently extract spectra using a SpatialPointsDataFrame
#' 
#' The `hs_extract_pts` function efficiently extracts hyperspectral reflectance 
#' data at spatial points. This is more efficient than reading the entire 
#' raster using `hs_read`, then using `raster::extract`.
#' 
#' Points that are not contained within the extent of the raster will return NA
#' values.
#' 
#' @inheritParams hs_read
#' @param pts a `SpatialPointsDataFrame` object defining points where data will
#' be extracted
#' 
#' @export
hs_extract_pts <- function(filename, pts, bands) {
  stopifnot(all(bands > 0))
  bands <- sort(bands)
  h5_dims <- hs_dims(filename)
  h5_extent <- hs_extent(filename)
  r <- raster::raster(h5_extent, crs = hs_proj4string(filename), 
                      nrows = h5_dims[2], ncols = h5_dims[3])
  
  pts[['row_identifier']] <- seq_len(nrow(pts))
  pts_in_raster <- raster::intersect(pts, r)
  if (nrow(pts_in_raster) == 0) {
    stop('None of the points are within the hyperspectral raster extent.')
  }
  
  # assuming pts is a spatialpointsdataframe, find row/col indices
  cells <- raster::cellFromXY(r, pts_in_raster)
  rowcols <- raster::rowColFromCell(r, cells)
  
  # read the values from the raster
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]
  vals <- matrix(nrow = length(cells), ncol = length(bands))
  for (i in seq_along(cells)) {
    vals[i, ] <- reflectance[bands, rowcols[i, 'row'], rowcols[i, 'col']]
  }
  vals <- hs_clean(vals, reflectance)
  file_h5$close_all()
  
  colnames(vals) <- hs_wavelength(filename, bands)
  extracted_vals <- cbind(pts_in_raster['row_identifier'], vals)
  res <- sp::merge(pts, extracted_vals, by = 'row_identifier')
  cols_to_keep <- !(names(res) %in% c('row_identifier', 
                                      'coords.x1', 
                                      'coords.x2'))
  res[, cols_to_keep]
}


#' Get the dimensions of a hyperspectral reflectance HDF5 file
#' 
#' The `hs_dims` function returns the dimensions of reflectance data contained
#' within an HDF5 file for NEON's L3 hyperspectral reflectance data. 
#' In most cases, these dimensions will be 426 X 1000 X 1000: 426 bands and 
#' images that are 1000m by 1000m in their spatial extent at 1 meter resolution.
#' 
#' @inheritParams hs_read
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
#' @inheritParams hs_read
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


#' Get the EPSG code of a hyperspectral image
#' 
#' This function returns the [EPSG](http://www.epsg.org/) code that defines
#' the spatial coordinate reference system of an image.
#' 
#' @inheritParams hs_read
#' @return character representation of an epsg code, e.g., '4326'
#' 
#' @export
hs_epsg <- function(filename) {
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
  site <- file_h5$ls()$name
  epsg_path <- paste0(site, '/Reflectance/Metadata/Coordinate_System/EPSG Code')
  epsg_code <- file_h5[[epsg_path]]$read()
  file_h5$close_all()
  epsg_code
}

#' Get the proj4string for a hyperspectral image
#' 
#' This returns the proj4string representation of the coordinate reference 
#' system of an image.
#' 
#' @inheritParams hs_read
#' @return character proj4string representation
#' 
#' @export
hs_proj4string <- function(filename) {
  epsg <- hs_epsg(filename)
  paste0("+init=epsg:", epsg)
}

#' Get wavelengths from a hyperspectral image
#' 
#' The `hs_wavelength` returns a character vector containing band indices and
#' wavelengths rounded to the nearest nanometer, e.g., band1_384nm, band2_389nm, 
#' where band... gives the band index (first band, second band, etc.) and 
#' ...nm gives the wavelength of that band in nanometers. 
#' 
#' @inheritParams hs_read
#' @return character vector containing band indices and wavelengths
#' @export
hs_wavelength <- function(filename, bands) {
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
  r <- hs_clean(r, reflectance)
  file_h5$close_all()

  raster::projection(r) <- hs_proj4string(filename)
  names(r) <- hs_wavelength(filename, bands)
  r
}


hs_clean <- function(vals, reflectance) {
  na_value <- reflectance$attr_open('Data_Ignore_Value')$read()
  vals[vals == na_value] <- NA
  scale_factor <- reflectance$attr_open('Scale_Factor')$read()
  vals <- vals / scale_factor
  vals
}


calculate_index_extent <- function(clip_extent, h5_extent){
  dummy_raster <- raster::raster(h5_extent, resolution = 1)
  cells <- raster::cellsFromExtent(dummy_raster, extent = clip_extent)
  colrange <- range(raster::colFromCell(dummy_raster, cells)) # x
  rowrange <- range(raster::rowFromCell(dummy_raster, cells)) # y
  list(xmin_idx = min(colrange), xmax_idx = max(colrange),
       ymin_idx = min(rowrange), ymax_idx = max(rowrange))
}
