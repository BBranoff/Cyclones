Cyclones
================
Ben Branoff
January 23, 2026

# Coming Soon

This is the future site of the Cyclones R package currently in
development. A stable repo will be uploaded as soon as its available.
May 1, 2026

# Highlights

- Produce raster images of cyclone winds, rains, and storm surge, either
  at specific times or as aggregate
  - Choose data sources and temporal and spatial resolutions

# Future Implementation Ideas

- Integrate NISAR data for:
  - Soil moisture
  - Flooding extent
  - Biomass loss
- Integrate STAR NESDIS SAR for windspeed comparisons
- Integrate LANDSAT for NDVI
- Integrate available LiDAR or NAIP clouds for canopy structure
- Integrate SLOSH for storm surge

# Introduction

R utilities for producing variable temporal and spatial resolution
gridded and time-stacked representations of tropical cyclone wind,
precipitation, and storm surge.

A note on lists, rasters, wrapping, and parallelization. Most of the
functions in ‘Cyclones’ produce and operate over lists, primarily
because they work nicely with lapply() and its parallel equivalents.
This makes it relatively easy to process multiple storms in parallel or
to process a single storm in parallel. However, raster data via ‘terra’
can not be passed across parallel workers unless they are wrapped. Thus,
for now, the raster products in Cyclones will always be wrapped and must
first be unwrapped before they can be utilized in traditional workflows.

Parallel processing is mostly only necessary for multiple storms or for
high resolution outputs and we recommend users familiarize themselves
with the sequential default functions and only then utilize parallel
when ready for batch processing. For internal functions, parallelism is
done via snowfall at the timestep level. This is best for high
resolution (both temporal (\<15 mins) and/or spatial (\<10000 m))
outputs and especially when producing Thin Plate Spline (TPS) outputs,
and is probably best utilized if at least 5-10 cpus could be set aside
for the effort, but there hasnt been any in-depth test of performance
yet. Parallelizing at the storm level is demonstrated below in snowfall,
but could be done with other libraries at the user’s discretion and also
could be done with only a few (\<5) cpus. The functions will not allow
you to do both internal and external parallelization, but otherwise, the
user is expected to know how many cpus can be safely utilized on their
system.

# Get Storm Tabular Data

Use the ‘get_storms()’ function to gather the available time series,
location, and meteorological information for a particular storm. If a
local, pre-downloaded IBTrACS .csv or .nc data source is available, it
can be used, either as the file location or as a pre-loaded dataset. If
not, the data will be downloaded from the web.

``` r
###  install from github
```
