interp_rast <- function(r,track,parallel){
  #r <- rast(lapply(r,unwrap))
  #times <- sort(unique(time(r)))
  #times <- sort(unique(lapply()))
  rtrim <- lapply(r,function(x) trim(unwrap(x)))
  rtrim <- lapply(rtrim,wrap)
  if (parallel){
    cat("\rinterpolating to desired frequency.")
    #sfExport("rtrim","track")
    r_int <- sfClusterApplyLB(1:(length(rtrim)-1),interp_func,rtrim,track)
  } else {
    r_int <- lapply(1:(length(rtrim)-1),interp_func,rtrim,track)
  }
  r_int <- lapply(unlist(r_int),function(x) project(unwrap(x),unwrap(r[[1]])))
  r_int
}
interp_func <- function(x,rtrim,track){
  #if (parallel) r <- unwrap(r)
  r1 <- unwrap(rtrim[[x]])
  r2 <- unwrap(rtrim[[x+1]])
  t1 <- time(r1)
  t2 <- time(r2)
  track1 <- track|>filter(date==t1,location=="track points")
  track2 <- track|>filter(date==t2,location=="track points")
  ##  to get the rocis
  bounds1 <- st_as_sf(as.polygons(classify(r1[[1]], cbind(floor(minmax(r1[[1]])[1]), Inf, 1)), dissolve = TRUE))
  bounds2 <- st_as_sf(as.polygons(classify(r2[[1]], cbind(floor(minmax(r2[[1]])[1]), Inf, 1)), dissolve = TRUE))
  custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",
                    st_coordinates(track1 |>st_transform(4326))[,"X"],
                    " +lat_0=",
                    st_coordinates(track2|>st_transform(4326))[,"Y"])
  coord_dif <- st_coordinates(track1) - st_coordinates(track2)
  r2s <- shift(r2,dx=coord_dif[1],dy=coord_dif[2])
  r2s <- project(r2s,r1)
  dates <- track$date[track$date<=t2&track$date>=t1] |> unique()
  rall <- c(rep(r1,length(dates)-1),r2s)
  ####  This will be per date, not per layer. So should be 1/3 the number of layers
  areas <- approx(c(as.numeric(st_area(bounds1)),as.numeric(st_area(bounds2))),n=length(dates))$y
  ###  now need to interpolate
  ##  keep the first date and the last date
  rall <- rast(lapply(1:nlyr(rall),function(x,y){if(x>3&x<(nlyr(y)-2)) setValues(y[[x]],NA) else(y[[x]])},y=rall))
  # must first set the time, approximate uses it to determine distance
  terra::time(rall) <- rep(dates,each=nlyr(r1))
  ##  approximate the values in between the endpoints, band by band
  rall <- rast(lapply(unique(names(rall)), function(n,y){approximate(y[[names(y)==n]])},y=rall))
  ### reorder by time
  rall <- rall[[order(terra::time(rall))]]
  ###  now shift the in between dates to their correct position
  ###  the last date will be the first date in the next iteration, so no need to include it
  ###  unless this is the last iteration
  rall <- lapply(seq_along(dates),function(d,y,t,ar){
    ##  the first date is already correct, no need to to anything
    if (d==1){
      rp <- y[[time(y)==dates[d]]]
    } else if (d==length(dates)&x<(length(rtrim)-1)){
      ##  last date can be dropped unless this is the last iteration
      rp <- NULL
    }else{
      ##  the following dates need to be shifted and cropped
      coord_dif <- st_coordinates(t|>filter(location=="track points",date==dates[d]))-
        st_coordinates(t|>filter(location=="track points",date==dates[1]))
      y1 = y[[time(y)==dates[d]]]
      ###  everything needs to be shifted to its actual coordinates
      rp <- shift(y1,coord_dif[1], coord_dif[2])
      ##  all but the last need to be trimmed
      if(d<length(dates)){
        rp <- mask(rp,st_buffer(t|>filter(location=="track points",date==dates[d]),sqrt(ar[[d]]/pi)))
      }
    }
    rp
  },y=rall,t=track,ar=areas)
  cat(paste0("\rinterpolating to desired frequency: % ",round(100*x/length(rtrim))))
  rall <- Filter(Negate(is.null), rall)
  lapply(rall,wrap)
}

