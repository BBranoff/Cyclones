tps_fast <- function(tracks,r,todir=NULL,overwrite=FALSE,smooth=FALSE,trim=FALSE,eye_option="given"){
  #reye <- r
  r <- rast(unwrap(r$rTempP))
  centers=tracks |>filter(location=="track points",source=="native")|>mutate(dt=difftime(lead(date),date,units="hours"))
  lines <- tracks |> filter(date %in% centers$date)
  #if (L==length(unique(lines$date[!lines$location=="track points"]))){return(NULL)}
  lines$date <- as.POSIXct(lines$date,tz="UTC")
  ###  find the outer bounds of the storm
  lne <- lines |> filter(location!="track")
  dts <- unique(lne$date[!lne$location %in% c("track","track points")])
  ###  need to have all the timesteps first
  ###  get all of those
  rasts <- lapply(dts,function(d,lne,cents,r,eye_opt){
    line1 <- lne |>
      filter(date==d)
    rs1 <- tps_interpolate(line1,cents[cents$date==d,],r,trim=F,eye_opt=eye_opt)
    cat(which(dts==d))
    rs1
  },lne=lne,cents=centers,r=r,eye_opt="maxwind")
  ###  then interpolate between
  rasts_int <- lapply(seq_along(dts), function(d,cents){
    d2 <- dts[d+1]
    r1 <- rasts[[d]]
    r2 <- rasts[[d+1]]
    d <- dts[d]
    custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",
                                       st_coordinates(cents[cents$date==d,] |>st_transform(4326))[,"X"],
                                      " +lat_0=",
                                       st_coordinates(cents[cents$date==d,]|>st_transform(4326))[,"Y"])
    coord_dif <- st_coordinates(cents[cents$date==d,]) -
         st_coordinates(cents[cents$date==d2,])
    r2s <- shift(r2,dx=coord_dif[1],dy=coord_dif[2])
    r2s <- project(r2s,r1)
    dates <- tracks$date[tracks$date<=d2&tracks$date>=d] |> unique()
    rall <- c(rep(r1,length(dates)-1),r2s)
    ##  keep the first date and the last date but set all dates in between to NA
    rall <- rast(lapply(1:nlyr(rall),function(x,y){if(x>3&x<(nlyr(y)-2)) setValues(y[[x]],NA) else(y[[x]])},y=rall))
    ##  approximate the values in between the endpoints, band by band
    rall <- rast(lapply(unique(names(rall)), function(n,y){approximate(y[[names(y)==n]])},y=rall))
    #set the time and remove the last date if not the end point
    terra::time(rall) <- rep(dates,3)
    ### reorder by time
    rall <- rall[[order(terra::time(rall))]]
    ###  now shift the in between dates to their correct position
    rall <- lapply(1:nlyr(rall),function(x,y,t){
       coord_dif <- st_coordinates(t|>filter(location=="track points",date==terra::time(y[[x]])))-
         st_coordinates(t|>filter(location=="track points",date==terra::time(y[[1]])))
       if(x>3) shift(y[[x]],coord_dif[1], coord_dif[2])
       else(y[[x]])},y=rall,t=tracks)
    rall <- rast(lapply(rall, function(x) resample(x,r,method="bilinear")))
    rall <- extend(rall,r)

  })


  d1 <- unique(dts)[L]
  #d2 <- unique(dts)[which(unique(dts)==d1)+1]
  #dates <- centers$date[centers$date<=d2&centers$date>=d1] |> unique()
  tofiles <-paste0(todir,"/tps_",unique(tracks$ID[!is.na(tracks$ID)]),"_",unique(format(d1,"%Y%m%d%H%M")),".tif")
  if (!is.null(todir)){
    if (all(file.exists(tofiles))&overwrite==FALSE) return(tofiles)
  }
  #timestep <- as.numeric(difftime(d2,d1,units="hours"))

  ###  now mask out the outer bounds of the storm
  ###  use the interpolated roci values for each timestep
  ##  add the method to the names
  rall <- rs1
  names(rall) <- paste0("tps_",names(rall))
  rall <- mask(rall,st_buffer(line1|>filter(location=="track points"),line1|>filter(location=="track points")|>pull(roci.m)))
  rall <- rast(lapply(unique(dates),function(d,r,t){
    mask(rall[[terra::time(rall)==d]],st_buffer(c[c$date==d,],c$roci.m[c$date==d]))},
    r=rall,t=tracks))
  ## remove last date if not the end point
  #if (!end) {rall <- rall[[1:(nlyr(rall)-3)]];tofiles <- tofiles[-length(tofiles)]}


  #line2 <- lne |>
  #  filter(date==d2)
  #if(d1==unique(dts)[length(unique(dts))-1]){end=T}else{end=F}
  ###  create custom crs centered on current track segment
  ##  this reduces geometry calculation errors from using a general global crs
  #custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",
  #                  st_coordinates(centers[centers$date==d1,] |>st_transform(4326))[,"X"],
  #                 " +lat_0=",
  #                  st_coordinates(centers[centers$date==d1,]|>st_transform(4326))[,"Y"])
  ###  calculate the delta x and delta y between the two locations
  ###  we will use this to align the two dates on top of eachother so we can then interpolate between them
  # coord_dif <- st_coordinates(centers[centers$date==d2,] |> st_transform(custCRS))-
  #   st_coordinates(centers[centers$date==d1,]|> st_transform(custCRS))
  ## get all of the dates and calculate the time difference between each
  ##  this is important for calculating energy dissipation, which is a function of the time
  #dt <- unique(difftime(dates, lag(dates),unit="hours"))
  #dt <- as.numeric(dt[!is.na(dt)])
  ###  now shift the second line segment to be centered on the first
  ##  this will allow us to interpolate raster representations between the two times
  #line2.2 <- st_as_sfc(line2 |> st_transform(custCRS))-matrix(data=coord_dif,ncol=2)
  #line2.2 <- st_sf(geom=line2.2)
  # sf::st_crs(line2.2) <- custCRS
  #line2.2 <- line2.2 |> st_transform(st_crs(line1))  |>
  #  st_as_sf() |>
  #  rename(geometry=geom)
  #line2.2<- cbind(line2.2,st_drop_geometry(line2))
  ##  create the wind velocity and power fields for the two segments,
  ##  which are now centered over the same geographical space

  ###  include the time duration in the power calculation
  # rs2 <- tps_interpolate(line2.2,centers[centers$date==d2,],r,eye_opt=eye_option)
  # ###  stretch so they have the same spatial extent
  # ###  this is necessary to stack and interpolate between them
  # shared_extent <- union(ext(rs1), ext(rs2))
  # rs1 <- extend(rs1,shared_extent)
  # rs2 <- extend(rs2,shared_extent)
  # rall <- c(rep(rs1,length(dates)-1),rs2)
  # ##  keep the first date and the last date but set all dates in between to NA
  # rall <- rast(lapply(1:nlyr(rall),function(x,y){if(x>3&x<(nlyr(y)-2)) setValues(y[[x]],NA) else(y[[x]])},y=rall))
  # ##  approximate the values in between the endpoints, band by band
  # rall <- rast(lapply(unique(names(rall)), function(n,y){approximate(y[[names(y)==n]])},y=rall))
  # #set the time and remove the last date if not the end point
  # terra::time(rall) <- rep(dates,3)
  # ### reorder by time
  # rall <- rall[[order(terra::time(rall))]]
  # ###  now shift the in between dates to their correct position
  # rall <- lapply(1:nlyr(rall),function(x,y,c){
  #   coord_dif <- st_coordinates(c[c$date==terra::time(y[[x]]),])-
  #     st_coordinates(c[c$date==terra::time(y[[1]]),])
  #   if(x>3) shift(y[[x]],coord_dif[1], coord_dif[2])
  #   else(y[[x]])},y=rall,c=centers)
  # rall <- rast(lapply(rall, function(x) resample(x,r,method="bilinear")))
  #rall <- extend(rall,r)
  ###  now mask out the outer bounds of the storm
  ###  use the interpolated roci values for each timestep
  ##  add the method to the names
  rall <- rs1
  names(rall) <- paste0("tps_",names(rall))
  rall <- mask(rall,st_buffer(line1|>filter(location=="track points"),line1|>filter(location=="track points")|>pull(roci.m)))
  #rall <- rast(lapply(unique(dates),function(d,r,c){
  #  mask(rall[[terra::time(rall)==d]],st_buffer(c[c$date==d,],c$roci.m[c$date==d]))},
  #  r=rall,c=centers))
  ## remove last date if not the end point
  #if (!end) {rall <- rall[[1:(nlyr(rall)-3)]];tofiles <- tofiles[-length(tofiles)]}
  msg = paste("\rCalculating wind field via Thin Plate Spline for ",paste(unique(lines$name),unique(format(lines$date,"%Y")),sep="_"),
              ": %",round(100*L/ (length(unique(lines$date[!lines$location %in% c("track","track points")]))+1),1))
  #if (messagefun=='sf') sfClusterEval(cat(msg))
  #else
  message(msg,appendLF = FALSE)
  Sys.sleep(0.01)
  ##  wrap the results for parallel compatability
  if (!is.null(todir)){
    lapply(tofiles, function(x){
      if (!(file.exists(x)&overwrite==FALSE)) writeRaster(rall[[format(terra::time(rall),"%Y%m%d%H%M")==gsub(".tif","",strsplit(basename(x),"_")[[1]][6])]],x,overwrite=overwrite)
    })
    return(tofiles)
  }else{
    return(wrap(rall))
  }
}
