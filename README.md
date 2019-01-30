
<!-- README.md is generated from README.Rmd. Please edit that file -->

# neonaop

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)

The goal of neonaop is to make data from the NEON AOP easier to use. The
National Ecological Observatory Network (NEON) collects hyperspectral
imagery via its Aerial Observation Platform (AOP) at a 1 meter spatial
resolution for 426 different wavelengths.

## Installation

You can install the development version of neonaop via:

``` r
#install.packages('devtools')
devtools::install_github('earthlab/neonaop')
```

## Example

This is a basic example which shows you how to read some bands from L3
hyperspectral reflectance data as a multi-layer raster:

``` r
library(neonaop)
library(raster)
library(viridis)

path_to_file <- system.file('extdata', 'ex.h5', package = 'neonaop')
r <- hs_read(path_to_file, bands = c(1, 50, 100, 400))
r
#> class       : RasterBrick 
#> dimensions  : 50, 50, 2500, 4  (nrow, ncol, ncell, nlayers)
#> resolution  : 1, 1  (x, y)
#> extent      : 257000, 257050, 4111950, 4112000  (xmin, xmax, ymin, ymax)
#> coord. ref. : +init=epsg:32611 +proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 
#> data source : in memory
#> names       : band1_384nm, band50_629nm, band100_879nm, band400_2382nm 
#> min values  :      0.0101,       0.0042,        0.0615,         0.0033 
#> max values  :      0.1488,       0.1183,        0.7033,         0.2072
```

``` r
plot(r, col = cividis(100), axes = FALSE, box = FALSE)
```

<img src="man/figures/README-plot-1.png" width="100%" />
