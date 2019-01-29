# Generate fake h5 file containing a smaller subset of NEON AOP data
#
# This is an internal script that is used to create inst/extdata/ex.h5,
# which is used for examples and also for testing. This is just one band from
# the NEON AOP hyperspectral data.
#
# First, get a real NEON AOP file for SJER locally via
# neonUtilities::byTileAOP(dpID="DP3.30006.001", site="SJER", year="2017",
#                          easting = 257179, northing = 4111285,
#                          check.size = FALSE)

create_new_example_file <- FALSE

if (create_new_example_file) {
  library(hdf5r)
  path <- list.files(pattern = 'NEON_D17_SJER_DP3_257000_4111000_reflectance.h5',
                     recursive = TRUE, full.names = TRUE)
  real_file <- H5File$new(path, mode = 'r+')

  # Create the example file and some relevant groups ------------------------
  ex_path <- 'inst/extdata/ex.h5'
  unlink(ex_path)
  ex_file <- H5File$new(ex_path, mode = 'w')
  ex_file$create_group('SJER')
  ex_file$create_group('SJER/Reflectance')
  ex_file$create_group('SJER/Reflectance/Metadata')
  ex_file$create_group('SJER/Reflectance/Metadata/Coordinate_System')
  ex_file$create_group('SJER/Reflectance/Metadata/Spectral_Data')



  # Save a subset of bands with relevant attributes -------------------------
  refl_path <- 'SJER/Reflectance/Reflectance_Data'
  real_reflectance <- real_file[[refl_path]]

  n <- real_reflectance$dims[1] # take all 426 wavelengths
  nxy <- 30
  first_n_bands <- real_reflectance$read(args = list(1:n,
                                                     1:nxy,
                                                     1:nxy))
  ex_file[[refl_path]] <- first_n_bands

  na_value <- real_reflectance$attr_open('Data_Ignore_Value')$read()
  h5attr(ex_file[[refl_path]], 'Data_Ignore_Value') <- na_value

  scale_factor <- real_reflectance$attr_open('Scale_Factor')$read()
  h5attr(ex_file[[refl_path]], 'Scale_Factor') <- scale_factor



  # Save map info metadata --------------------------------------------------

  info_path <- 'SJER/Reflectance/Metadata/Coordinate_System/Map_Info'
  real_map_info <- real_file[[info_path]]
  ex_file[[info_path]] <- real_map_info$read()

  epsg_path <- 'SJER/Reflectance/Metadata/Coordinate_System/EPSG Code'
  real_epsg <- real_file[[epsg_path]]
  ex_file[[epsg_path]] <- real_epsg$read()


  # Save band wavelengths in nanometers ------------------------------------
  wv_path <- 'SJER/Reflectance/Metadata/Spectral_Data/Wavelength'
  real_wv <- real_file[[wv_path]]
  ex_file[[wv_path]] <- real_wv$read()[1:n] # only save first n bands

  # Close files ------------------------------------------------------------
  real_file$close_all()
  ex_file$close_all()
  unlink('DP3.30006.001/', recursive = TRUE)
}
