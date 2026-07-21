#' @importFrom terra vect project crs
#' @keywords internal

get_dir <- function(x, y, lonlat, landIntersect,DirRas,meth) {
  site_bear <- calculate_bearing(y,x,lonlat[2],lonlat[1],DirRas)
  terra::values(DirRas) <- site_bear
  if (lonlat[2] >= 0) {
    DirRas <- DirRas-90
  } else {
    DirRas <- DirRas+90
  }
  if (meth=="boose"){
    if (lonlat[2] >= 0) {
      DirRas[landIntersect == 1] <- DirRas[landIntersect == 1] - 40
      DirRas[landIntersect == 0] <- DirRas[landIntersect == 0] - 20
    } else {
      DirRas[landIntersect == 1] <- DirRas[landIntersect == 1] + 40
      DirRas[landIntersect == 0] <- DirRas[landIntersect == 0] + 20
    }
  }
  DirRas[DirRas< 0] <- DirRas[DirRas< 0] + 360
  DirRas[DirRas> 360] <- DirRas[DirRas> 360] - 360
  return(DirRas)
}

calculate_bearing <- function(lat1, lon1, lat2, lon2,tmpras){
  ##  from Boose HURRECON model
  d2r = 0.017453292519943295  # pi / 180
  r2d = 57.29577951308232  # 180 / pi
  # nearly same point
  bear=rep(NA,length(lon1))
  if (any(abs(lat2 - lat1) < 0.000001 & (abs(lon2 - lon1) < 0.000001))){
    bear[is.na(bear)] = rep(0,length(lon1))
  } else{
    ##  date line
    if (any(lon1) > 90 & lon2 < -90) lon2 <- lon2 + 360
    if (any(lon1) < -90 & lon2 > 90)   lon1[lon1< -90] <- lon1[lon1< -90] + 360
    #if (any(lon1<0)) lon1 <- lon1+360
    #if (lon2<0) lon2 <- lon2+360
    # same longitude
    if (any(lon1 == lon2)){
      if (any(lat1[lon1 == lon2] > lat2)) bear[lon1 == lon2&lat1 > lat2] = 180
      bear[lon1 == lon2&lat1 <= lat2] = 0
    }

    # different longitude
    # convert degrees to radians
    rlat1 = d2r*lat1
    rlat2 = d2r*lat2
    rlon1 = d2r*lon1
    rlon2 = d2r*lon2
    B2 <- atan2(sin(rlon2-rlon1)*cos(rlat2), cos(rlat1)*sin(rlat2) - sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1))
    #B2 =atan2(cos(rlat1)*sin(rlat2) -sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1),
    #          sin(rlon2-rlon1)*cos(rlat2))
    # convert radians to degrees
    B = r2d*B2
    #B[which(B<0)] <- B[which(B<0)]+360
    #bear[is.na(bear)] <- B[is.na(bear)]
    bear[which(lon1 < lon2)]= B[which(lon1 < lon2)] # quadrants I, IV
    bear[which(lon1 >= lon2)] = 360 + B[which(lon1 >= lon2)]  # quadrants II, III

  }
  ###  Hurrecon calculates the bearing from multiple sites to the storms center at different times
  ###  we have reversed it and are taking the bearing from the center at one time to every cell around the center
  ###  so must reverse the bearing
  return (bear)
}
get_dir_old <- function(r,C){
  ###
  ## find out if storm straddles antimeridian
  e <- ext(project(r,"epsg:4326"))
  if (e[1]>0&e[2]<0) {
    r_points <- as.data.frame(rotate(project(r,"epsg:4326")),xy=TRUE)
    r_points$bearing <- geosphere::bearing(r_points[,c("x","y")],C |> st_transform(4326)|>st_shift_longitude()|>st_coordinates())
    r_points$WD <- (r_points$bearing+ 90 - 180)%% 360
    WD <-rasterize(vect(r_points[,c("x","y","WD")],geom=c("x","y"),crs="epsg:4326"),rotate(project(r,"epsg:4326")),field="WD")
  }else{
    r_points <- as.data.frame(project(r,"epsg:4326"),xy=TRUE)
    r_points$bearing <- geosphere::bearing(r_points[,c("x","y")],C |> st_transform(4326)|>st_coordinates())
    r_points$WD <- (r_points$bearing+ 90 - 180)%% 360
    WD <-rasterize(vect(r_points[,c("x","y","WD")],geom=c("x","y"),crs="epsg:4326"),project(r,"epsg:4326"),field="WD")
  }
  WD <- project(WD,r)

  ## antimeridian storms
  #r_points <- as.data.frame(project(r,crs2),xy=TRUE) ## this is a 'shifted' crs, 4326 will not work at the antimeridian
  #r_points$bearing <- geosphere::bearing(r_points[,c("x","y")],C |> st_transform(4326)|>st_coordinates())
  #r_points$WD <- (r_points$bearing+ 90 - 180)%% 360
  #WD <-rasterize(vect(r_points[,c("x","y","WD")],geom=c("x","y"),crs="+proj=longlat +datum=WGS84 +lon_0=180"),project(r,"+proj=longlat +datum=WGS84 +lon_0=180"),field="WD")
  WD
}
