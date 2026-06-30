## originally from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
#' @importFrom terra ext
#' @importFrom dplyr nth
#' @export
boose <- function(L,extents,tmpRas,todir,overwrite=FALSE,smooth=FALSE,eye_option=NULL) {

  id =unique(extents$ID[!is.na(extents$ID)])
  cens <- extents |> filter(location=="track points")
  cens <- shift_track(cens)
  cent <- cens |> slice(L)
  if (!is.null(todir)){
    tofile <- paste0(todir,"/boose_",id,"_",unique(format(cent$date,"%Y%m%d%H%M")),".tif")
    if (file.exists(tofile)&overwrite==FALSE) return(tofile)
  }
  tmpRasP <- tmpRas$rTempP
  tmpRas <- tmpRas$rTempG
  if (class(tmpRas)[1]=="PackedSpatRaster"){tmpRas=suppressWarnings(rast(unwrap(tmpRas)));tmpRasP=suppressWarnings(rast(unwrap(tmpRasP)))}
  vx <- cent |> pull(vxDeg_geo)
  vy <- cent |> pull(vyDeg_geo)
  vh <- cent |> pull(stormSpeed_geo)
  #tmpRas <- crop(tmpRas,ext(extents))#unlist(st_drop_geometry(cent[,c('xmin','xmax','ymin','ymax')]))))
  rho <- 1 # air density
  msw = cent$maxwind.ms#/1.94384 ## convert to m/s ?
  maxpress= cent$maxpress.mb/0.01  ##  convert to atm?
  minpress=cent$minpress.mb/0.01
  rmw <- cent$rmw.m/1000
  ###  using higher resolution (l=50) reduces speed significantly
  land <- rnaturalearth::ne_countries() |> st_cast("POLYGON")|>st_transform(crs(tmpRas))|>rasterize(y=tmpRas,field=1,background=0)
  # Computing coordinates of raster
  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001
  deltP <- maxpress - minpress
  if (deltP==0) deltP=1
  b <- rho * exp(1) * msw**2 / deltP
  vr <- distEye
  vr <- sqrt((rmw / vr)**b * exp(1 - (rmw / vr)**b))
  if (cent$centerY_shifted_geo >= 0) {
    # Northern Hemisphere, t is clockwise
    angle <- atan2(vy, vx) - atan2(y, x)
  } else {
    # Southern Hemisphere, t is counterclockwise
    angle <- atan2(y, x) - atan2(vy, vx)
  }
  landIntersect <- terra::extract(land,crds)
  vr[landIntersect == 1] <- 0.8 * (msw - (1 - sin(angle[landIntersect == 1])) * vh[landIntersect == 1] / 2) * vr[landIntersect == 1]
  vr[landIntersect == 0] <- (msw - (1 - sin(angle[landIntersect == 0])) * vh / 2) * vr[landIntersect == 0]
  vr <- round(vr,3)
  dist <- sqrt(x * x + y * y)
  vr[dist > 2.5] <- NA
  terra::values(tmpRas) <- vr

  direction <- computeDirectionBoose(x, y, st_coordinates(st_transform(cent,4326))[,2], landIntersect)
  direction <- rast(t(matrix(direction,ncol=dim(tmpRas)[1])),crs=crs(tmpRas),extent=ext(tmpRas))
  ##  to get power, must first compute pressure field
  ###  the pressure field should just be a linear interpolation from the center to the roci
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
  names(tmpRas) <- paste0("boose_",names(tmpRas))
  terra::time(tmpRas) <- rep(cent$date,each=3)
  if (smooth){
    tmpRas <- stormRsmooth(tmpRas,s_res,tmpRasP)
  }else{
    tmpRas <- rast(lapply(tmpRas,function(x,y) {project(x,y)},y=tmpRasP))
  }
  message(paste0("\rCalculating Boose wind field for ",paste(unique(cens$name),unique(format(cens$date,"%Y")),sep="_"),
                 ": %",round(100*L/(nrow(cens)+1),1)),appendLF = FALSE)


  return(wrap(tmpRas))
}

