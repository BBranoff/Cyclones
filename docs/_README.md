Cyclones
================
July, 2026

# Introduction

For workflow vignettes, see the [documention
website](https://bbranoff.github.io/Cyclones/)

Cyclones is a collection of R utilities for producing variable temporal
and spatial resolution gridded and time-stacked representations of
historic tropical cyclone wind, precipitation, and storm surge.
Processing can be done sequentially for low and moderate resolution
products, or in parallel for high resolution.

For now, Cyclones is hosted only here on Github and can be installed as
follows:

``` r
## install.packages("remotes")
remotes::install_github("BBranoff/Cyclones")

## for the development version (unstable)
remotes::install_github("BBranoff/Cyclones@development")
```

# Recent ’Development” Branch Updates (Aug 2026):

- Added ‘aggregate_product()’ function and ‘agg’ argument to get_winds()
  and get_precip() for whole storm aggregated outputs
- Changed ‘get_winds()’ processing to only compute native timesteps and
  interpolate everything in between for faster processing
- Added pre-loaded msw interval model coefficients to ‘cons_stormdat()’.
  This allows the msw interval translation to occur when limited data is
  supplied. \*\* Adding ‘.new’ to the “msw_int” argument will override
  these coefficients and compute new models on the supplied data.
- Added Google Cloud as a source for ecmwf precipitation. Using the
  source=“emcwf_GC” argument will retrieve the data from Google Cloud
  instead of the ECMWF service. \*\* This was done to avoid getting
  queued in the ECMWF service. \*\* The GC products are .tif, whereas
  the ECMWF products are .grib, and the extents are slightly different.
  Thus, the ‘get_ecmwf()’ function will delete any storm data that is
  incomplete if one or the other is requested, ensuring a storm is
  sourced from either ‘ecmwf’ or ‘ecmef_GC’, but not both.
- Storm surge products from storms spanning the 180 meridian were very
  large, encompassing all longitudes. This has been resolved by shifting
  and/or rotating these storm extents to retreive the necessary data.
  Results are thus much smaller in size for these storms.
- Added partial ERDDAP functionality to ‘get_storms()’, allowing single
  storm downloads instead of entire IBTrACS subsets. NOTE: There are
  issues with the ERDDAP server and this functionality is not yet
  complete.
- added ‘todir’ argument to ‘make_extents()’, allowing indidivual storm
  extents to be saved upon completion.

# Highlights

- Produce raster images of cyclone winds, rains, and storm surge, either
  at specific times or as aggregate
  - Choose data sources and temporal and spatial resolutions
- Ingest storm tracks and wind and pressure extents from IBTrACS
  - Consolidate data across meteorological agencies
- For winds:
  - Choose among Boose, Holland, Willoughby, or Thin Plate Spline
    methods
  - Model missing extents from similar storms
  - Calculate maximum sustained winds, wind direction, and power
  - Compare calculated winds against hurricane hunter drop sondes
- For rainfall:
  - Choose among ERA5, MSWEP, or GPM sources
  - Calculate PRE-storm and storm precipitation
- For storm surge:
  - Conglomerate across measured (NOAA & USGS), and modeled (CMIP6)
    water levels
  - Compare against SRTM elevation to approximate water depth

# Future Implementation Ideas

- Integrate NISAR data for:
  - Soil moisture
  - Flooding extent
  - Biomass loss
- Integrate STAR NESDIS SAR for windspeed comparisons
- Integrate LANDSAT for NDVI
- Integrate available LiDAR or NAIP clouds for canopy structure
- Integrate SLOSH for storm surge
