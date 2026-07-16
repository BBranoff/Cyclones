## originally from stormR (Delaporte B, Ibanez T, Keppel G, Jullien S, Menkes C, Arsouze T (2024))
## adapted here for compatability with other functions
#' @importFrom terra ext
#' @importFrom dplyr nth
#' @export
boose <- function(L,extents,tmpRas,todir,smooth=FALSE,eye_option=NULL) {

  dat <- prep_theoretical_data(L,extents,tmpRas)
  cent <- dat$cent
  msw <- dat$msw
  tmpRasP <- tmpRas$rTempP
  tmpRas <- dat$tmpRas
  rho <- 1 # air density
  b <- rho * exp(1) * msw**2 / dat$deltP

  ###  using higher resolution (l=50) reduces speed significantly
  land <- rnaturalearth::ne_countries() |> st_cast("POLYGON")|>
    st_transform(crs(tmpRas))|>
    rasterize(y=tmpRas,field=1,background=0)
  vr <- sqrt((dat$rmw / dat$vr)**b * exp(1 - (dat$rmw / dat$vr)**b))
  vx <- cent |> pull(vxDeg_geo)
  vy <- cent |> pull(vyDeg_geo)
  vh <- cent |> pull(stormSpeed_geo)
  if (cent$centerY_shifted_geo >= 0) {
    # Northern Hemisphere, t is clockwise
    angle <- atan2(vy, vx) - atan2(dat$y, dat$x)
  } else {
    # Southern Hemisphere, t is counterclockwise
    angle <- atan2(dat$y, dat$x) - atan2(vy, vx)
  }
  landIntersect <- terra::extract(land,dat$crds)
  vr[landIntersect == 1] <- 0.8 * (msw - (1 - sin(angle[landIntersect == 1])) * vh[landIntersect == 1] / 2) * vr[landIntersect == 1]
  vr[landIntersect == 0] <- (msw - (1 - sin(angle[landIntersect == 0])) * vh / 2) * vr[landIntersect == 0]
  vr <- round(vr,3)
  ##  distance only needed if masking by distance from center
  ##  really this should be theoretically derived as the asymptote of the curve
  #dist <- sqrt(dat$x * dat$x + dat$y * dat$y)
  #vr[dist > 2.5] <- NA

  terra::values(tmpRas) <- vr
  direction <- computeDirectionBoose(dat$crds[,1], dat$crds[,2], st_coordinates(st_transform(dat$cent,4326)),landIntersect,tmpRas)
  message(paste0("\rCalculating Boose wind field for ",paste(unique(cent$name),unique(format(cent$date,"%Y")),sep="_"),
                 ": %",round(100*L/(nrow(extents|>filter(location=="track points"))+1),1)),appendLF = FALSE)
  package_theoretical_data (dat,tmpRas,tmpRasP,direction,meth="boose")
}

