Cyclones
================
Ben Branoff
January 23, 2026

# Coming Soon

This is the future site of the Cyclones, or CycDstr R package currently
in development. A stable repo will be uploaded as soon as its available.
Feb. 26, 2026

# Introduction

R utilities for producing variable temporal and spatial resolution
gridded and time-stacked representations of tropical cyclone wind,
precipitation, and storm surge.

A note on parallelization. Most of the primary functions in this library
can be parallelized. This is mostly only necessary for multiple storms
or for high resolution outputs and we recommend users familiarize
themselves with the sequential default functions and then utilizing
parallel when ready for batch processing. For internal functions,
parallelism is done via snowfall at the timestep level. This is best for
high resolution (both temporal (\<15 mins) and/or spatial (\<10000 m))
outputs and especially when producing Thin Plate Spline (TPS) outputs,
and is probably best utilized if at least 5-10 cpus could be set aside
for the effort, but there hasnt been any in-depth test of performance
yet. Parallelizing at the storm level is demonstrated below in snowfall,
but could be done with other libraries at the user’s discretion and also
could be done with only a few (\<5) cpus. The functions will not allow
you to do both internal and external parallelization, but otherwise, the
user is expected to know how many cpus can be safely utilized on their
system.

# Highlights
* Produce rasters of cyclone winds, rains, and storm surge
  * Choose data sources and temporal and spatial resolutions

# Coming soon
* Integrate NISAR data for:
  * Soil moisture
  * Flooding extent
  * Biomass loss

# Get Storm Tabular Data

Use the ‘get_storms()’ function to gather the available time series,
location, and meteorological information for a particular storm. If a local, 
pre-downloaded IBTrACS .csv or .nc data source is available, it can be used, 
either as the file location or as a pre-loaded dataset. If not, the data will be
downloaded from the web. The IBTrACS data comes with data from multiple 
meteorological agencies across the world, each often following different 
protocols for data reporting. By default, data will be consolidated to represent
USA variables, meaning those variables will be selected if multiple are 
available, or if only non-USA variables are available, those will be translated 
to represent USA reporting standards. USA variables are preferred because they 
are most ubiquitous, but these preferences can be changed with the 
'cons_stormdat()' function. 

``` r
library(Cyclones)

##  storms can be singular or plural and can be identified specifically, or not
storms <- get_storms(source="hurdat",name="Maria",basin="NA",season=2017)
storms <- get_storms(source="hurdat",name=c("Maria","MICHAEL"),basin="NA",season=2017)
names(storms)

##  alternatively, all of the data can be downloaded and filtered later
##  we recommend this if wind models need to be constructed in later steps
##  downloading all storms will likely require extending the timeout time
##  options(timeout = 300)
allstorms <- get_storms(ib_filt="ALL")
storms <- allstorms[grep("MARIA_2017|MICHAEL_2018",names(allstorms))]
```

# Get Wind Products

With the tabular data loaded, they can now be used to build rasters of
wind, precipitation, and/or storm surge, as well as a combined Tropical
Cyclone Severity Scale (TCSS) (Bloemendale et al. (2021)).

## Build Extent Models

If the Thin Plate Spline (TPS) wind method is desired, wind extent
models are built into a lookup table (The empirical TPS method is shown
below to be more accurate than the theoretical models). These lookup
tables should be based on as full or limited a set of storms as
necessary for the objective. For generalized modeling, its best to use a
full set of storms.

``` r
## build the models based on the entire storm dataset.  
## The models are built into lookup tables, 
## with a different model for every quadrant of the storm, 
## for every storm size and in every basin. 
mods <- build_models(allstorms)
```

## Make Wind Extent Spatial Features

The above models will allow the reconstruction of certain wind extents
when missing, which is often the case for older (pre 2018) storms.
Again, if generalized wind fields are desired, the other non-TPS methods
can be used as explained further below. But if the objective is to
create wind field predictions specific to each quadrant and each basin
and each storm category based on the other storms in the data, we use
these models and the TPS interpolation.

``` r
##  generate linestring (default) and polygon (optionally) wind extent features for each timestep 
##  in each storm. When used as shown below with multiple storms as an input, the functions will 
##  only perform the action on the first storm
windextents <- make_extent(storms,mods=mods,type=c("linestrings","polygons"))

##  To apply the function to all storms, use lapply or a parallel equivalent (snowfall::sfLapply)
windextents <- lapply(storms,make_extent, mods=mods)
library(snowfall)
sfInit(parallel=TRUE, cpus=2) 
sfLibrary(sf)
sfLibrary(terra)
sfLibrary(dplyr)
windextents <- sfLapply(storms,make_extent, mods=mods)

##  If both linestrings and polygons were returned, they can be pulled out accordingly
linestrings <- lapply(windextents,function(x) x[x$extent_type=="linestrings",])
polygons <- lapply(windextents,function(x) x[x$extent_type=="polygons",])
##  the individual storms are projected in their own crs, centered on the storm. To get them all in the same CRS, if desired
linestrings <- lapply(linestrings,st_transform,crs=4326)
```

These spatial features can also be good for certain storm visualizations
in which a fully gridded raster is unnecessary.

![](README_files/figure-gfm/plot-extents-1.png)<!-- -->![](README_files/figure-gfm/plot-extents-2.png)<!-- -->

## Build Wind Rasters

Next, we interpolate and/or model the winds across space and time,
creating a more complete gridded representation that can be useful in
any number of applications. This is accomplished through either one of
three theoretical models, or through a more empirical Thin Plate Spline,
which utilizes the above linestrings, which were generated from a set of
models. The ‘make_winds’ function accepts both spatial and temporal
resolution specifications. Although accuracy and smoothness are
optimized at smaller resolutions, computation times suffer. For initial
testing, a spatial resolution of 20000 m and a temporal resolution of 60
mins is reasonable and default. Because terra::raster objects are saved
on disc rather than memory, they can not be passed and collected as part
of parallel processes. Therefore, results from get_wind are wrapped
rasters in case a parallel operation is desired. They must be separated
and unwrapped with a simple ‘unwrap()’ call before they can be used for
most operations. All storm layers are projected onto a Lambert Azimuthal
Equal Area reference system that is centered on the storm’s extent,
which preserves area across large geographic regions and avoids
geographic coordinate issues at the 180th meridian.

``` r
winds <- lapply(windextents,get_wind,methods=c("all"))
### unwrap files in memory
winds <- lapply(winds, function(x) lapply(x,unwrap))
names(winds)
## load files saved to disk from the above
winds <- lapply(winds, function(x) rast(unlist(x)))
## load files saved to disk previously, providing the directory names as saved
##  in this case, if files were saved other than the working directory, provide that as the 'todir' argument
winds <- lapply(list("MARIA2017_5000_10","MICHAEL2018_5000_10"),get_wind,loadrasts=FALSE)
winds <- unlist(winds)

###  the layers are stored in nested lists in which individual time steps are nested within each storm
winds$MARIA_2017_NA_2017260N12310

plot(winds$MARIA_2017_NA_2017260N12310$MSW_MARIA_2017_TPS)

####################
#  parallel times via snowfall
#  best for high resolution output of one storm
####################

##  take only one storm at a time
winds <- get_wind(windextents[[1]],methods=c("all"),t_res=5,s_res=15000,parallel=TRUE,cpus=15)
##  take only one storm at a time
winds <- get_wind(windextents[[1]],methods=c("all"),t_res=5,s_res=15000,parallel=TRUE,cpus=15)

### unwrap files in memory
winds <- lapply(winds, function(x) lapply(x,unwrap))
## load files saved to disk
winds <- lapply(winds, function(x) rast(unlist(x)))

####################
#  parallel storms via snowfall
#  best for low to medium resolution for multiple storms
####################
library(snowfall)
sfStop()
logtmp <- tempfile(fileext=".txt")  ## if debugging is of interest
sfInit(parallel=TRUE, cpus=2,slaveOutfile = logtmp) ## WARNING: Do not overdo the number of CPUs. Verify your machine's capacity beforehand.
##  the parallel cpus will each do one storm, so no point in giving it more cpus than there are storms to process.
sfLibrary(sf)
sfLibrary(terra)
sfLibrary(dplyr)
sfLibrary(fields) # for the TPS models
sfLibrary(geosphere)  # for the bearing function
sfLibrary(Cyclones)# for the current library
#sfLibrary(snowfall) # to be able to write to log file
winds <- sfLapply(windextents,Cyclones:::get_wind,methods="all",s_res=20000,t_res=60,parallel=TRUE)
sfStop()
## unwrap
winds <- lapply(winds, function(x) lapply(x,unwrap))
```

![](README_files/figure-gfm/getwinds-1.png)<!-- -->

    ## Warning: Removed 264 rows containing missing values or values outside the scale range
    ## (`geom_raster()`).

![](README_files/figure-gfm/getwinds-2.png)<!-- -->

Another look at one time step with the published wind extents overlayed
from the National Hurricane Center.
![](README_files/figure-gfm/plotwinds-1.png)<!-- -->

The rasters produced via the above get_winds routine contain a layer for
each time step as specified in t_res. They are listed and wrapped, and
may need to be separated and unwrapped before they can be used. Another
way to represent the rasters, however is by aggregating across the
storm’s lifetime. Getting the maximum sustained wind for each pixel, for
example, would be a popular product. The agg_winds function will provide
this as well as other common metrics, such as the total power exerted on
a 1 sq. meter surface during the duration of the storm, the direction of
the maximum winds, or the duration of category 1+ or category 3+ winds.

## Compare With Any Corresponding Drop Sondes

The compare_winds function will compare the computed raster data with
drop sonde data from the National Hurricane Center, if available. Drop
sondes are acquired from [NOAA’s Hurricane Research
Division](https://www.aoml.noaa.gov/hrd/data_sub/hurr.html "NOAA's Hurricane Research Division").
They are high accuracy gps beacons, whose trajectory and sensors provide
high frequency wind information (i.e. velocity, direction). Data are
available for select storms dating back to 1960. If available, data is
ingested, parsed, aggregated, and saved to either a temp drive or a
permanent repository, which makes processing more stable if the same
storms will be reanalyzed in the future (no need to download). Once the
sonde data is ready, they are cross referenced with the wind rasters and
matched by location (accuracy dependent on raster resolution) and time
(within 5 minutes default). Although drop sondes collect information
along the entire vertical profile of the storm, only surface winds (\<30
m in elevation) are compared here, as those are most comparable to the
raster layers.

The compare_winds function will output a list of spatial feature data
frames in a geographic coordinate system, with each entry corresponding
to one set of sonde measurements and the corresponding extracted raster
value. It can be used to assess and examine the accuracy of the wind
products as they related to direct measurements from the sondes.

``` r
##  Either the rasters (as filenames or loaded into memory), the linestrings, or both can be compared against the sondes
##  indicating the ddir will save and load the sonde data from that directory
comps  <- compare_winds(rasts=winds,ddir="./Drop Sondes/")
comps_r <- do.call(rbind,comps$rsamps)
###  get the value closest to the ground for each sonde in each storm
comps_r_filt <- comps_r |> group_by(sonde, storm, method,var) |>
  filter(ZW.m==min(ZW.m,na.rm=TRUE)) |>
  ungroup()
```

    ## Warning: Removed 3507 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](README_files/figure-gfm/compare_winds-1.png)<!-- -->

    ## Warning: Removed 2751 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](README_files/figure-gfm/compare_winds-2.png)<!-- --> The above
graphic represents the drop sonde comparisons for the raster products of
about 40 storms since 2016, each at 5 km spatial resolution and 10 min
temporal resolution. Although each storm will be different, overall the
TPS method outperforms the other models in terms of the root mean square
error (RMSE). Most of the error for all methods stems from sonde
locations in the eye-wall, where wind dynamics can change drastically
over relatively small distances.

# Get Precipitation Products

The Cyclones package can also ingest and process precipitation data for
all storms after 1940. There are three potential sources for the data:
[the European Space Agency’s ERA5
(ecmwf)](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels),
[the National Aeronautics and Space Administration’s (NASA) Global
Precipitation Measurement Mission
(gpm)](https://gpm.nasa.gov/missions/GPM), and [the Multi-Source
Weighted-Ensemble Precipitation (mswep)](https://www.gloh2o.org/). All
data sources require some sort initial registration before they can be
accessed via the api, as well as additional packages for extracting and
assembling the data. They are global, hourly or three-hourly data at ~25
and 10 km resolution. The get_precip function will request data from the
source (default is ecmwf), download it into a given directory (if
specified), and then create storm specific gridded precipitation
outputs. As with wind, the outputs can remain stacked at different
times, or they can be returned aggregated to produced summary values for
each pixel across the duration of the storm.

``` r
##  to download the ERA5 website, you must first register with either ECMWF or NASA, or request access to MSWEP:
##  [register with ECMWF for ERA5.](https://www.ecmwf.int/user/login)
##  [register with NASA Earth Data for GPM IMERG.](https://urs.earthdata.nasa.gov/)
##  [request access to MSWEP](https://www.gloh2o.org/)

##  once registered with ECMWF, you will recieve an API token that looks like this:  abcd1234-foo-bar-98765431-XXXXXXXXXX
##  This will give you access to the data repositories.
##  To set your api key in R and make it accessible for future use this:
##  (This only needs to be run once per R installation and setup and should not be run every time) 
#   library(ecmwfr)
#   wf_set_key(key = "abcd1234-foo-bar-98765431-XXXXXXXXXX")

##  or, once registered at nasa earthdata, run the following with your username and password 
##  library(earthdatalogin)
#   edl_netrc(username = "XXXX",pasword="XXXX")

##  or, once the MSWEP data has been granted, access is through Google Drive and can be granted via (follow the pop-up prompts):
#   library(googledrive)
#   drive_auth()

##  The system is now configured to access the data.
##  To not have to re-download storm data again, store the files in a local repository
##  otherwise, if dpath is not provided, they will be stored in a temporary directory
##  either way, the function will provide the resulting rasters, which can be processed and stored as desired
##  you may also choose to aggregate the rasters into single (non-timestamped) whole storm layers by setting agg=TRUE
##  the default source is ecmwf
precip <- lapply(storms,get_precip,dpath="./Precips")
precip_agg <- lapply(precip,agg_precip)
```
