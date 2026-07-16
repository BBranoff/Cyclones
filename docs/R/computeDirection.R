computeDirection <- function(x, y, lat) {
  azimuth <- -(atan2(y, x) - pi / 2)

  azimuth[azimuth < 0] <- azimuth[azimuth < 0] + 2 * pi

  if (lat >= 0) {
    direction <- azimuth * 180 / pi - 90
  } else {
    direction <- azimuth * 180 / pi + 90
  }

  direction[direction < 0] <- direction[direction < 0] + 360
  direction[direction > 360] <- direction[direction > 360] - 360
  ###  flip direction to indicate where the wind is coming from, not going to
  direction <- ifelse(direction<180,direction+180,direction-180)
  return(direction)
}
computeDirectionBoose <- function(x, y, lonlat, landIntersect,DirRas) {
  ##  StormR methods
  # azimuth <- -(atan2(y, x) - pi / 2)
  # azimuth[azimuth < 0] <- azimuth[azimuth < 0] + 2 * pi
  #
  # if (lonlat[2] >= 0) {
  #   direction <- azimuth * 180 / pi - 90
  #   direction[landIntersect == 1] <- direction[landIntersect == 1] - 40
  #   direction[landIntersect == 0] <- direction[landIntersect == 0] - 20
  # } else {
  #   direction <- azimuth * 180 / pi + 90
  #   direction[landIntersect == 1] <- direction[landIntersect == 1] + 40
  #   direction[landIntersect == 0] <- direction[landIntersect == 0] + 20
  # }
  #
  # direction[direction < 0] <- direction[direction < 0] + 360
  # direction[direction > 360] <- direction[direction > 360] - 360
  # ###  flip direction to indicate where the wind is coming from, not going to
  # direction <- ifelse(direction<180,direction+180,direction-180)

  ## Boose methods from GitHub Hurrecon
  ##  https://github.com/hurrecon-model/HurreconPython/blob/master/python/hurrecon.py
  site_bear <- calculate_bearing(y,x,lonlat[2],lonlat[1],DirRas)
  values(DirRas) <- site_bear
   if (lonlat[2] >= 0) {
     DirRas[landIntersect == 1] <- DirRas[landIntersect == 1] - 40
     DirRas[landIntersect == 0] <- DirRas[landIntersect == 0] - 20
     DirRas[DirRas< 0] <- DirRas[DirRas< 0] + 360
   } else {
     DirRas[landIntersect == 1] <- DirRas[landIntersect == 1] + 40
     DirRas[landIntersect == 0] <- DirRas[landIntersect == 0] + 20
     DirRas[DirRas> 0] <- DirRas[DirRas> 0] - 360
   }
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
      B2 =atan2(cos(rlat1)*sin(rlat2) -sin(rlat1)*cos(rlat2)*cos(rlon2-rlon1),
                sin(rlon2-rlon1)*cos(rlat2))
      # convert radians to degrees
      B = r2d*B2
      B[which(B<0)] <- B[which(B<0)]+360
      bear[is.na(bear)] <- B[is.na(bear)]
      #bear[which(lon1 < lon2)]= B[which(lon1 < lon2)] # quadrants I, IV
      #bear[which(lon1 >= lon2)] = 360 + B[which(lon1 >= lon2)]  # quadrants II, III

  }
  return (bear)
}
