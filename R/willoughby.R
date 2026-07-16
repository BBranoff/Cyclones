## adapted from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
#' @importFrom terra time
willoughby <- function(L,extents,tmpRas,smooth=FALSE,eye_option=NULL) {
  dat <- prep_theoretical_data(L,extents,tmpRas)
  tmpRasP <- tmpRas$rTempP
  tmpRas <- dat$tmpRas
  cent <- dat$cent
  lat=cent$centerY_shifted_geo
  msw <- dat$msw
  ##  willoughby has a method to estimate the RMW if missing
  rmw <- dat$rmw
  if (is.na(rmw)) rmw <- round(46.4 * exp(-0.0155 * dat$msw + 0.0169 * abs(lat))) #in km
  x1 <- 287.6 - 1.942 *msw + 7.799 * log(rmw) + 1.819 * abs(lat)
  x2 <- 25
  a <- 0.5913 + 0.0029 * msw - 0.1361 * log(rmw) - 0.0042 * abs(lat)
  n <- 2.1340 + 0.0077 * msw - 0.4522 * log(rmw) - 0.0038 * abs(lat)
  vr <- dat$vr
  vr[dat$vr >= rmw] <- msw * ((1 - a) * exp(-abs((vr[dat$vr >= rmw] - rmw) / x1)) + a * exp(-abs(vr[dat$vr >= rmw] - rmw) / x2))
  vr[dat$vr < rmw] <- msw * abs((dat$vr[dat$vr < rmw] / rmw)^n)
  tmpRasA <- computeAsymmetry("Chen",vr,dat$x,dat$y,cent$vxDeg_geo, cent$vyDeg_geo,
                             cent$stormSpeed_geo,
                             dat$vr, rmw, lat)
  tmpRasA$wind <- round(tmpRasA$wind,3)
  wind <- direction <- dat$tmpRas
  terra::values(wind) <- tmpRasA$wind
  terra::values(direction) <- tmpRasA$direction
  message(paste0("\rCalculating Willoughby wind field for ",paste(unique(cent$name),unique(format(cent$date,"%Y")),sep="_"),
                 ": %",round(100*L/(nrow(extents|>filter(location=="track points"))+1),1)),appendLF = FALSE)
  package_theoretical_data (dat,wind,tmpRasP,direction,meth="willoughby")
}
