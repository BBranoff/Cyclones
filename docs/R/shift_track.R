#' @importFrom dplyr filter
#' @importFrom sf st_length st_line_sample st_crs st_geometry
#' @keywords internal
shift_track <- function(points){
    ###  interpolate (linearly) missing values and create variables necessary for theoretical models
  points_geo <- points |>
    ###  get the shifted coords necessary for the theoretical equation methods
    st_transform(4326) |>
    mutate(dt=difftime(lead(date),date,units="hours"))|>
    st_shift_longitude()

  points_geo$centerX_shifted_geo <- st_coordinates(points_geo)[,1]
  points_geo$centerY_shifted_geo <- st_coordinates(points_geo)[,2]
  ###  storm speed in km/h ?
  points_geo$stormSpeed_geo <- terra::distance(
    x=vect(cbind(points_geo$centerX_shifted_geo, points_geo$centerY_shifted_geo),crs="epsg:4326"),
    y=vect(cbind(dplyr::lead(points_geo$centerX_shifted_geo), dplyr::lead(points_geo$centerY_shifted_geo)),crs="epsg:4326"),
    pairwise=TRUE) * (0.001 / as.numeric(points_geo$dt)) / 3.6
  points_geo <- points_geo |>
    group_by(date)|>
    mutate(rmw = min(dist.m[dist.m>0],na.rm=TRUE)) |>
    ungroup()|>
    mutate(vxDeg_geo=(dplyr::lead(as.numeric(centerX_shifted_geo)) - as.numeric(centerX_shifted_geo)) / as.numeric(dt),
           vyDeg_geo=(dplyr::lead(as.numeric(centerY_shifted_geo)) - as.numeric(centerY_shifted_geo)) / as.numeric(dt)) |>
    tidyr::fill(stormSpeed_geo:vyDeg_geo)
  points <- st_transform(points_geo,st_crs(points))

  points
}
