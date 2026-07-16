## originally from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
holland <- function(L, extents,tmpRas,smooth=FALSE,quiet=FALSE,eye_option=NULL) {
  dat <- prep_theoretical_data(L,extents,tmpRas)
  cent <- dat$cent
  tmpRasP <- tmpRas$rTempP
  tmpRas <- dat$tmpRas
  rho <- 1.15 # air density
  lat=cent$centerY_shifted_geo
  b <- rho * exp(1) * dat$msw**2 / dat$deltP
  f <- 2 * 7.29 * 10**(-5) * sin(lat) # Coriolis parameter
  vr <- sqrt(b / rho * (dat$rmw / dat$vr)**b * dat$deltP * exp(-(dat$rmw / dat$vr)**b) + (dat$vr * f / 2)**2) - dat$vr* f / 2
  ## set the missing values to the minimum
  vr[is.na(vr)] <- min(vr,na.rm=TRUE)
  tmpRasA <- computeAsymmetry("Chen",vr,dat$x,dat$y,cent$vxDeg_geo, cent$vyDeg_geo,
                              cent$stormSpeed_geo,
                              dat$vr, dat$rmw, lat)
  #dist <- sqrt(dat$x * dat$x + dat$y * dat$y)
  tmpRasA$wind <- round(tmpRasA$wind,3)
  wind <- direction <- dat$tmpRas
  terra::values(wind) <- tmpRasA$wind
  terra::values(direction) <- tmpRasA$direction
  message(paste0("\rCalculating Holland wind field for ",paste(unique(cent$name),unique(format(cent$date,"%Y")),sep="_"),
                 ": %",round(100*L/(nrow(extents|>filter(location=="track points"))+1),1)),appendLF = FALSE)
  package_theoretical_data (dat,wind,tmpRasP,direction,meth="holland")

}
