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
computeDirectionBoose <- function(x, y, lat, landIntersect) {
  azimuth <- -(atan2(y, x) - pi / 2)

  azimuth[azimuth < 0] <- azimuth[azimuth < 0] + 2 * pi

  if (lat >= 0) {
    direction <- azimuth * 180 / pi - 90
    direction[landIntersect == 1] <- direction[landIntersect == 1] - 40
    direction[landIntersect == 0] <- direction[landIntersect == 0] - 20
  } else {
    direction <- azimuth * 180 / pi + 90
    direction[landIntersect == 1] <- direction[landIntersect == 1] + 40
    direction[landIntersect == 0] <- direction[landIntersect == 0] + 20
  }

  direction[direction < 0] <- direction[direction < 0] + 360
  direction[direction > 360] <- direction[direction > 360] - 360
  ###  flip direction to indicate where the wind is coming from, not going to
  direction <- ifelse(direction<180,direction+180,direction-180)
  return(direction)
}
