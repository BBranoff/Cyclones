#' @importFrom terra vect project crs
#' @keywords internal
get_dir <- function(r,C){
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
