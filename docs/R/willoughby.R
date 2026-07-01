## adapted from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
#' @importFrom terra time
willoughby <- function(L,extents,tmpRas,todir,overwrite=FALSE,smooth=FALSE,eye_option=NULL) {
  id =unique(extents$ID[!is.na(extents$ID)])
  cens <- extents |> filter(location=="track points")
  cens <- shift_track(cens)
  cent <- cens |> slice(L)
  if (!is.null(todir)){
    tofile <- paste0(todir,"/willoughby_",id,"_",unique(format(cent$date,"%Y%m%d%H%M")),".tif")
    if (file.exists(tofile)&overwrite==FALSE) return(tofile)
  }
  tmpRasP <- tmpRas$rTempP
  tmpRas <- tmpRas$rTempG
  if (class(tmpRas)[1]=="PackedSpatRaster"){tmpRas=suppressWarnings(rast(unwrap(tmpRas)));tmpRasP=suppressWarnings(rast(unwrap(tmpRasP)))}
  #tmpRas <- crop(tmpRas,ext(extents))
  msw = cent$maxwind.ms#round(cent$maxwind.ms/1.94384)
  lat=cent$centerY_shifted_geo
  rmw <- round(cent$rmw.m/1000)
  if (is.na(rmw)) rmw <- round(46.4 * exp(-0.0155 * msw + 0.0169 * abs(lat))) #in km
  maxpress= cent$maxpress.mb/0.01
  minpress=cent$minpress.mb/0.01
  deltP <- maxpress - minpress
  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001
  x1 <- 287.6 - 1.942 * msw + 7.799 * log(rmw) + 1.819 * abs(lat)
  x2 <- 25
  a <- 0.5913 + 0.0029 * msw - 0.1361 * log(rmw) - 0.0042 * abs(lat)
  n <- 2.1340 + 0.0077 * msw - 0.4522 * log(rmw) - 0.0038 * abs(lat)
  vr <- distEye
  vr[distEye >= rmw] <- msw * ((1 - a) * exp(-abs((distEye[distEye >= rmw] - rmw) / x1)) + a * exp(-abs(distEye[distEye >= rmw] - rmw) / x2))
  vr[distEye < rmw] <- msw * abs((distEye[distEye < rmw] / rmw)^n)
  if (any(vr>100)) browser()
  tmpRasA <- computeAsymmetry("Chen",vr,x,y,cent$vxDeg_geo, cent$vyDeg_geo,
                             cent$stormSpeed_geo,
                             distEye, rmw, lat)
  dist <- sqrt(x * x + y * y)
  tmpRasA$wind[dist > 2.5] <- NA
  tmpRasA$wind <- round(tmpRasA$wind,3)
  terra::values(tmpRas) <- tmpRasA$wind
  direction <- computeDirection(x, y, st_coordinates(st_transform(cent,4326))[,2])
  direction <- rast(t(matrix(direction,ncol=dim(tmpRas)[1])),crs=crs(tmpRas),extent=ext(tmpRas))
  ##  to get power, must first compute pressure field
  ###  the pressure field should just be a linear interpolation from the center (minpress) to the roci (maxpress)
  ##  start with the distance to the center
  distr <- rast(matrix(distEye,ncol=dim(tmpRas)[1]),crs=crs(tmpRas),extent=ext(tmpRas))
  distr <- resample(distr,tmpRas)
  ##  the pressure shouldnt increase after the roci, so set all values greater than roci to roci
  distr[distr>(cent$roci.m/1000)] <- cent$roci.m/1000
  P <- minpress+(distr/0.001)*(deltP/(cent$roci.m))
  dens <- P*100/(287.058*298)
  kW <- as.numeric(cent$dt)*(0.5*dens*1*(tmpRas)^3)/1000
  kW[is.na(kW)] <- 0
  terra::units(tmpRas) <- "m/s"
  terra::units(kW) <- "kW"
  terra::units(direction) <- "deg"
  tmpRas <- rast(list(power.kW=kW,msw.ms=tmpRas,windir.deg=direction))
  tmpRas <- mask(tmpRas,tmpRas$msw.ms)
  names(tmpRas) <- paste0("willoughby_",names(tmpRas))
  terra::time(tmpRas) <- rep(cent$date,each=3)
  if (smooth){
    tmpRas <- stormRsmooth(tmpRas,s_res,tmpRasP)
  }else{
    tmpRas <- project(tmpRas,tmpRasP)#rast(lapply(tmpRas,function(x,y) {project(x,y)},y=tmpRasP))
  }

  message(paste0("\rCalculating Willoughby wind field for ",paste(unique(cens$name),unique(format(cens$date,"%Y")),sep="_"),
                 ": %",round(100*L/nrow(cens),1)),appendLF = FALSE)
  wrap(tmpRas)
}
