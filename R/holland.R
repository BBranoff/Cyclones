## originally from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
holland <- function(L, extents,tmpRas,todir,overwrite=FALSE,smooth=FALSE,quiet=FALSE,eye_option=NULL) {
  id =unique(extents$ID[!is.na(extents$ID)])
  cens <- extents |> filter(location=="track points")
  cens <- shift_track(cens)
  cent <- cens |>slice(L)
  if (!is.null(todir)){
    tofile <- paste0(todir,"/holland_",id,"_",unique(format(cent$date,"%Y%m%d%H%M")),".tif")
    if (file.exists(tofile)&overwrite==FALSE) return(tofile)
  }
  tmpRasP <- tmpRas$rTempP
  tmpRas <- tmpRas$rTempG
  cens <- cens |> filter(location=="track points")
  if (class(tmpRas)[1]=="PackedSpatRaster"){tmpRas=suppressWarnings(rast(unwrap(tmpRas)));tmpRasP=suppressWarnings(rast(unwrap(tmpRasP)))}
  #tmpRas <- crop(tmpRas,ext(extents))
  msw = cent$maxwind.ms#round(cent$maxwind.ms/1.94384)
  rmw <- round(cent$rmw.m/1000)
  maxpress= cent$maxpress.mb/0.01
  minpress=cent$minpress.mb/0.01
  lat=cent$centerY_shifted_geo
  rho <- 1.15 # air densiy
  f <- 2 * 7.29 * 10**(-5) * sin(lat) # Coriolis parameter
  deltP <- maxpress - minpress
  if (deltP==0) deltP=1
  ##  if deltP is too small, the below exponent utilizing b will be too big to compute and will result in NA
  ##  a below step corrects for this by setting those value to the minimum
  b <- rho * exp(1) * msw**2 / deltP

  crds <- terra::crds(tmpRas, na.rm = FALSE)
  x <- crds[, 1] - cent$centerX_shifted_geo
  y <- crds[, 2] - cent$centerY_shifted_geo
  # Computing distances to the eye of the storm in m
  distEye <- terra::distance(
    x = crds,
    y = st_coordinates(cent|>st_transform(4326)|>st_shift_longitude()),
    lonlat = TRUE
  )* 0.001


  vr <- distEye
  vr <- sqrt(b / rho * (rmw / distEye)**b * deltP * exp(-(rmw / distEye)**b) + (distEye * f / 2)**2) - distEye* f / 2
  ####  the missing values are at the minimum
  ###   the missing values create issues in the asymmetry calc
  #plot(vr,ylim=c(0,.1))
  ##  find the missing values
  #vr2 <- vr
  #vr2[is.na(vr2)] <- min(vr,na.rm=TRUE,NA)-.1*min(vr,na.rm=TRUE,NA)
  #vr2[vr2!=min(vr2)] <- NA
  ### plot the missing values
  #points(vr2,col="red")

  ## set the missing values to the minimum
  vr[is.na(vr)] <- min(vr,na.rm=TRUE)

  if (any(is.na(vr))) browser()
  tmpRasA <- computeAsymmetry("Chen",vr,x,y,cent$vxDeg_geo, cent$vyDeg_geo,
                              cent$stormSpeed_geo,
                              distEye, rmw, lat)
  #if (L>=79) browser()
  dist <- sqrt(x * x + y * y)
  tmpRasA$wind[dist > 2.5] <- NA
  tmpRasA$wind <- round(tmpRasA$wind,3)
  terra::values(tmpRas) <- tmpRasA$wind
  direction <- computeDirection(x, y, st_coordinates(st_transform(cent,4326))[,2])
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
  names( tmpRas) <- paste0("holland_",names(tmpRas))
  terra::time(tmpRas) <- rep(cent$date,each=3)

  if (smooth){
    tmpRas <- stormRsmooth(tmpRas,s_res,tmpRasP)
  }else{
    tmpRas <- rast(lapply(tmpRas,function(x,y) {project(x,y)},y=tmpRasP))
  }

  #comps <- compare_winds(rasts=tmpRas,shape=data.frame(name=cent$name,date=cent$date))
  #comps$rsource="Holland"
  if (!quiet) message(paste0("\rCalculating Holland wind field for ",paste(unique(cens$name),unique(format(cens$date,"%Y")),sep="_"),
                 ": %",round(100*L/(nrow(cens)+1),1)),appendLF = FALSE)
  if (!is.null(todir)){
    if (!(file.exists(tofile)&overwrite==FALSE)) writeRaster(tmpRas,tofile,overwrite=overwrite)
    return(tofile)
  }else{
    return(  wrap(tmpRas))
  }
}
