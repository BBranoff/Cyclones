
#' @importFrom terra rast ext crop rasterize interpolate mask disagg trim aggregate app nlyr union unwrap units
#' @importFrom fields Tps
#' @importFrom sf st_agr
tps <- function(L,tracks,r,todir=NULL,overwrite=FALSE,smooth=FALSE,trim=FALSE,eye_option="given",cpus=FALSE){
  #reye <- r
  r <- rast(unwrap(r$rTempP))
  centers=tracks |>filter(location=="track points")|>mutate(dt=difftime(lead(date),date,units="hours"))
  lines <- tracks |> filter(date %in% centers$date)
  if (L==length(unique(lines$date[!lines$location=="track points"]))){return(NULL)}
  lines$date <- as.POSIXct(lines$date,tz="UTC")
  ###  find the outer bounds of the storm
  lne <- lines |> filter(location!="track")
  ### if the raster is coming from a parrallel call, unwrap it first
  #if (class(r)[1]=="PackedSpatRaster") r=suppressWarnings(rast(unwrap(r)));messagefun='sf'
  ###  get the date stamps of the lines
  dts <- unique(lne$date[!lne$location %in% c("track","track points")])
  d1 <- unique(dts)[L]
  #d2 <- unique(dts)[which(unique(dts)==d1)+1]
  #dates <- centers$date[centers$date<=d2&centers$date>=d1] |> unique()
  tofiles <-paste0(todir,"/tps_",unique(tracks$ID[!is.na(tracks$ID)]),"_",unique(format(d1,"%Y%m%d%H%M")),".tif")
  if (!is.null(todir)){
    if (all(file.exists(tofiles))&overwrite==FALSE) return(tofiles)
  }
  #timestep <- as.numeric(difftime(d2,d1,units="hours"))
  line1 <- lne |>
    filter(date==d1)
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
  rs1 <- tps_interpolate(line1,centers[centers$date==d1,],r,eye_opt=eye_option)
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
  terra::time(rall) <- rep(d1,3)
  rall <- mask(rall,st_buffer(line1|>filter(location=="track points"),line1|>filter(location=="track points")|>pull(roci.m)))
  rall <- project(rall,r)
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
 return(wrap(rall))
}
tps_interpolate <- function(line,center,r,trim=FALSE,eye_opt){
  ###  from Holland 1980
  Beta <- 1.0036+0.0173*as.numeric(max(line$kts,na.rm=TRUE)+
                                     0.0313*log(min(line$rmw.m,na.rm=TRUE))+
                                     0.0087*st_coordinates(center|>st_transform((4326)))[,"Y"])
  ##  now the pressure at each swath extent
  line <- line |> ungroup()|>
    mutate(P=minpress.mb+((maxpress.mb-minpress.mb)*exp(-(min(dist.m[kts>0],na.rm=TRUE)/dist.m)^Beta)))
  ### crop the line to the outer extent of the storm
  ###  if the roci was missing because the extents weren't modeled, set it to 500km
  if ("ROCI" %in% line$location) outer=line |> filter(location=="ROCI")
  else(outer=line |> filter(location=="track points")|>st_buffer(500000))
  ###   Using the original ROCI can create extreme and unnatural shifts in wind Velocity with the thin spline method
  ###   to avoid this, bump the 0 velocity line out a bit further
  ###   we will still zero out the original ROCI later, this is only for the thine spline interpolation
  sf::st_agr(outer) <- "constant"  ##  to avoid sf warning about repeating sub geometries
  line <- bind_rows(line,
                    st_buffer(outer|>st_cast("POLYGON"),outer$dist.m*.1) |>
                      mutate(location="ROCI2",dist.m=dist.m+dist.m*.1)|>st_cast("LINESTRING"))
  ##  cropping will chnage extent relative to other methods
  ##  but necessary to avoid way to ong processing times for interpolating empty space
  ##  results are reprojected back to original r afterwards
  r_cr <- crop(r,line|>filter(location=="ROCI2"))#rast(crs=crs(r),res=res(r),extent=ext(line))
  ####  remove the eye and the outer storm limits (note: currently keeping these in as they were important when comparing to sonde data. Removing them speeds up analysis significantly)
  ###  we can set those to zero later
  # n <- 2.1340 + 0.0077 * center$maxWind*0.51444 - 0.4522 * log(center$rmw/1000) - 0.0038 * abs(center$centerY_shifted_geo)
  # deltP <- (center$maxpress-center$minpress)/0.01
  # rho <- 1.15 # air density
  # b <- rho * exp(1) * (center$maxWind*0.51444)**2 / (deltP)
  # f <- 2 * 7.29 * 10**(-5) * sin(center$centerY_shifted_geo)
  # v_eye_will <- center$maxWind*0.51444 * abs((seq(1,center$rmw*0.001) / (center$rmw*0.001))^n)
  # v_eye_holl <- sqrt(b / rho * ((center$rmw/1000) / seq(1,center$rmw*0.001))**b * deltP * exp(-((center$rmw/1000)/ seq(1,center$rmw*0.001))**b) + (seq(1,center$rmw*0.001) * f / 2)**2) - seq(1,center$rmw*0.001)* f / 2
  # browser()
  # for (v in 1:length(v_eye_holl)){
  #   if ("eye" %in% line$location){
  #     if (line$dist_m_min[line$location=="eye"]*0.001<seq(1,center$rmw*0.001)[v]){
  #       eye <- st_buffer(center,1000*seq(1,center$rmw*0.001)[v]) |> st_cast("LINESTRING")
  #       eye$kts = v_eye_holl[v]
  #       eye$location = "eye_new"
  #       line <- bind_rows(line, eye)
  #     }
  #   }else{
  #   eye <- st_buffer(center,1000*seq(1,center$rmw*0.001)[v]) |> st_cast("LINESTRING")
  #   eye$kts = v_eye_holl[v]
  #   eye$location = "eye_new"
  #   line <- bind_rows(line, eye)}
  # }
  #inner <- line[line$dist_m_min==min(line$dist_m_min[line$dist_m_min>0],na.rm=TRUE),]
  #if ("eye" %in% inner$location) inner <- st_buffer(inner,inner$dist_m_min*1.1)
  #reye = list( rTempP=crop(reye$rTempP,inner,snap="out",touches=FALSE),
  #             rTempG=crop(reye$rTempG,st_shift_longitude(st_transform(inner,4326)),snap="out",touches=FALSE))
  #eye <- unwrap(holland(L,cens,reye,todir=NULL,quiet=TRUE))$holland_msw.ms
  #eye <- mask(eye,st_buffer(center,inner$dist_m_min))
  ###  rasterize the remaining swaths, one for velocity and one for pressure
  ###  the more features that are included, the longer the TPS interpolation takes.
  if ("eye" %in% line$location&&eye_opt=="maxwind") {
    line <- line|>mutate(kts=if_else(location=="eye",max(kts,na.rm=TRUE)*0.8,kts))
  }
  rv <-  rasterize(bind_rows(line|>filter(!location %in%c("ROCI","track","track points")),
                             st_buffer(line[line$location=="track points",],1)|>mutate(kts=10)|>st_cast("LINESTRING")), r_cr, "kts",touches=T,fun="mean")
  #rv2 <-  rasterize(
    ###  the normal extents
  #  line|>filter(!location %in%c("ROCI","track","track points")),
  #  r, "kts",touches=T,fun="mean")

  rP <-  rasterize(bind_rows(line|>filter(!location %in%c("ROCI","track","track points")),
                             st_buffer(line[line$location=="track points",],1)|>mutate(P=center$minpress.mb)|>st_cast("LINESTRING")), r_cr, "P",touches=T)
  ###  convert the xyz values to a dataframe for thin spline interpolation
  xyv <- as.data.frame(rv, xy=T,na.rm=F)
  xyP <- as.data.frame(rP, xy=T,na.rm=F)
  options(warn = -1)
  ###  fit the thin spline interpolation model
  tps_v <- Tps(xyv[,1:2], xyv[,3])
  #ftps_v <- fastTps(xyv[,1:2], xyv[,3],aRange=1)
  if (!any(is.na(xyP$P))) tps_P <-  Tps(xyP[,1:2], xyP[,3])
  else tps_P <- NULL
  ###  use the model to predict the unknown wind speeds
  p_v <- interpolate(r_cr, tps_v)
  if ("eye" %in% line$location&&eye_opt=="given") {
    p_v <- mask(p_v,line|>filter(location=="eye")|>st_cast("POLYGON"),updatevalue=0,inverse=TRUE,touches=FALSE)
  }
  ###
  ###     This was working decently, but does not take the actual 'eye' value
  ####
  ##  to be able to interpolate the wind changes at the eye wall, need a small res raster
  ##  regardless of what the output raster is
  if (max(line$kts,na.rm=TRUE)>=64&&eye_opt=="80-50"){
    r_eye <- crop(disagg(r_cr,res(r_cr)[1]/1000),st_buffer(line[line$location=="track points",],min(line$dist.m[line$kts==max(line$kts,na.rm=TRUE)],na.rm=TRUE)))
    ## how many steps?
    rv_eye <-rasterize(
      bind_rows(
        ##  this should be the maximum extent of the eye
        line[line$kts==max(line$kts,na.rm=TRUE),],
        ###  the below 'inner eyewalls'  will overlap the maximum wind if the rmw is too small
        ###  make sure its first so that the maximum wind gets priority if so
        st_buffer(line[line$location=="track points",],0.8*min(line$dist.m[line$kts==max(line$kts,na.rm=TRUE)],na.rm=TRUE))|>
          mutate(kts=max(line$kts,na.rm=TRUE)*.5)|>st_cast("LINESTRING"),
        ###  another
        st_buffer(line[line$location=="track points",],0.5*min(line$dist.m[line$kts==max(line$kts,na.rm=TRUE)],na.rm=TRUE))|>
          mutate(kts=max(line$kts,na.rm=TRUE)*.1)|>st_cast("LINESTRING"),
        ###  then the center
        st_buffer(line[line$location=="track points",],1)|>mutate(kts=0)|>st_cast("LINESTRING")),
      r_eye, "kts",touches=T,fun="mean")
    xyv_eye <- as.data.frame(rv_eye, xy=T,na.rm=F)
    tps_veye <- Tps(xyv_eye[,1:2], xyv_eye[,3])
    p_v_eye <-interpolate(r_eye,tps_veye)
    p_v_eye <- project(p_v_eye,p_v)
    ##  overlay the eye with the rest of the storm and take the minimum
    p_v <- app(c(p_v,p_v_eye),"min",na.rm=TRUE)
  }

  ###  in case the model predicts higher than recorded wind speeds, set them to the maximum known
  p_v <- terra::clamp(p_v,lower=0,upper=max(line$kts,na.rm=TRUE))
  #p_v <- terra::mosaic(terra::sprc(list(p_v,eye)), fun="min")
  #p_v[p_v>max(line$kts,na.rm=T)] <- max(line$kts,na.rm=T)


  ####  do the same for the minimum values

  #p_P[p_P<min(line$minpress,na.rm=T)] <- min(line$minpress,na.rm=T)
  #p_P[p_P>max(line$maxpress,na.rm=T)] <- max(line$maxpress,na.rm=T)
  #p_v[p_v<0] <- 0
  ###  now set the outer limits of the storm to missing or blank
  p_v <- mask(p_v,outer |> st_cast("POLYGON"),updatevalue=NA)
  ###  now model the winds inside of the eye based on the other methods
  # if ("eye" %in% unique(line$location)){
  #p_v <- mask(p_v,line |> filter(location=="eye")|>st_cast("POLYGON"),inverse=T,updatevalue=-999)
  #willoughby, distEYE is in km


  #  eye <- st_buffer(center,min(line$dist_m_min[line$dist_m_min>0],na.rm=TRUE))

  #p_v[p_v==-999] <- sample(seq(from=1.1,to=9.9,b=0.1), size=length(p_v[p_v==-999]), replace=TRUE)
  #  names(p_v) <- "lyr.1"
  #}

  ###  trim if desired
  if (trim) p_v <- trim(p_v)
  windir <- get_dir(p_v,center)

  ##  to calculate power:
  ###  the density of air is a conversion function of the pressure
  ##  assuming the temperature is around 25 C or 298K

  ##  power is then a function of the density and the pressure
  ##  https://www.e-education.psu.edu/emsc297/node/649
  if (!is.null(tps_P)){
    p_P <- interpolate(r_cr, tps_P)
    p_P <- terra::clamp(p_P,lower=min(line$minpress.mb,na.rm=T),upper=max(line$maxpress.mb,na.rm=T))
    dens <- p_P*100/(287.058*298)
    kW <- as.numeric(center$dt)*(0.5*dens*1*(p_v)^3)/1000
    kW[is.na(kW)] <- 0
  } else{
    kW <- setValues(p_v,NA)
  }
  ###  here, wind velocity in knots is converted to m/s
  p_v <- p_v/1.94384
  terra::units(p_v) <- "m/s"
  terra::units(kW) <- "kW"
  terra::units(windir) <- "deg"
  ##  and everything is converted from Watts to kWh by multiplying by the time duration
  ## we are calculating the power on 1 square meter of a surface for the duration of the current segment
  rast(list(power.kW=kW,msw.ms=p_v,windir.deg=windir))
}
