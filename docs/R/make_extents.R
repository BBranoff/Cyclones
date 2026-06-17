#' Create spatial simple features, linestrings or polygons, representing the storm track and wind extents of a storm at each timestep.
#'@description
#'`make_extents()` will accept the output of `get_storms()` and create simple features geometry collections and data frames for each storm, with geo-located features
#'representing the maximum wind extent of various wind speeds at each time step of the storm, as well as the outer pressure extent (ROCI) and tracks representing the center of the storm.
#'By default, only the native time steps are included (usually every 3 hours), but the 't_res' argument can be used to change the temporal resolution and create
#'interpolated extents at different, usually more frequent, time steps. Unlike other functions in Cyclones, the 'make_extents()' iterations are relatively faster and computation times
#'typically increase, rather than decrease, in parallel. Thus, parallel computation can only be performed at the storm level .
#'#'
#' @param storm a singular storm tabular data set as output by the `get_storms()` function. If a collection of storms is entered, only the first will be processed,
#' unless utilized in an lapply or parallel equivalent, in which case each will be processed individually.
#'
#' @param mods the result of the `build_models()` function. Used to model necessary information when missing or erroneous, usually for older storms.
#' If absent, only native extents are produced and other rare missing or likely erroneous data are not adjusted. Examples of erroneous data include radii
#' of maximum wind that are greater than extents of lower wind speeds, or minimum central pressures that are greater than the pressure at the last closed isobar.
#' @param type the desired class of output sf object, either 'linestrings','polygons', or 'all'.
#' @param t_res the desired temporal resolution of the output features, in minutes.
#' @returns Simple features geometry collection containing time stamped storm tracks and either linestrings, polygons, or both representing the maximum extent
#' of wind speeds, the eye wall (if present) and the Radius of the Last Closed Isobar (ROCI). All geometries are relative to a custom coordinate reference
#' system (crs) centered on the storm's centroid and in a Lambert azimuthal equal-area projection.
#' @importFrom purrr possibly
#' @importFrom dplyr last rowwise left_join join_by pull lead n all_of bind_cols
#' @importFrom sf st_sfc st_point st_transform st_polygon st_linestring st_as_sf st_buffer st_cast st_set_geometry st_set_crs st_drop_geometry st_union st_coordinates st_geometry st_geometry_type st_multilinestring st_line_merge st_line_sample
#' @keywords internal
make_extents <- function(storm,mods=NULL,type="linestrings",t_res=NULL,agency="CONS",cpus=NULL){
  if (!"linestrings"%in% type&(!"polygons" %in% type)&(!"all" %in% type)) stop("geometry return type unknown")
  storm <- checkstorm(storm,agency)
  #if (mods=="Cyclones") mods <- data(Cyclones::storm_mods)
  #mult <- storm$mult;storm <- storm$dta
  points <- storm |>
    mutate(ISO_TIME =as.POSIXct(ISO_TIME,format="%Y-%m-%d %H:%M:%S",tz="UTC"),
           MONTH  = format(ISO_TIME,"%m"),
           LAT =as.numeric(LAT),
           LON = as.numeric(LON),
           source="native",
           BASIN=if_else(is.na(BASIN),"NA",BASIN)) |>
    st_as_sf(coords=c("LON","LAT"),crs=4326,remove=FALSE)|>
    rename(date=ISO_TIME)
  if (!is.null(t_res)){
   points <- interp_track(points,t_res)|>
     st_as_sf(coords=c("LON","LAT"),crs=4326,remove=FALSE)
  }
  options(warn = 0)
  ###  parallel taking much longet than serial....
  if (!is.null(cpus)) {
    on.exit(sfStop())
    initiatepar(cpus,type="extents")
    dates_par <- split(seq_along(unique(points$date)),cut(seq_along(unique(points$date)),cpus))
    custCRSall <- paste0("+proj=laea +lat_0=",round(mean(points$LAT,na.rm=TRUE))," +lon_0=",round(mean(points$LON,na.rm=TRUE))," +lat_1=",ceiling((max(points$LAT,na.rm=TRUE)-min(points$LAT,na.rm=TRUE))/6),
                         " +lat_2=",ceiling(5*(max(points$LAT,na.rm=TRUE)-min(points$LAT,na.rm=TRUE))/6)," +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
    track <- sfLapply(dates_par,function(x,pnts,md,typ,cr){
      pnts <- pnts |> filter(date %in% unique(pnts$date)[x])
      tr <- lapply(1:nrow(pnts),apply_tracks,pnts,md,typ=type)
      tr <- lapply(tr,st_transform,cr)
      do.call(rbind,tr)
    },md=mods,pnts=points,typ=type,cr=custCRSall)
  }else{
    track <- lapply(1:nrow(points),apply_tracks,m=mods,pts=points,typ=type)
  }
  if (all(unique(unlist(lapply(track,names)))==c("swaths","linestrings"))) extents <- list(swaths=do.call(rbind,lapply(track,function(x) x$swaths)),linestrings=do.call(rbind,lapply(track,function(x) x$linestrings)))
  else extents <-do.call(rbind,track)
  cat("\n")
  #bundle <- list(polygons=swaths,linestrings=linestrings)
  #names(bundle) <- unique(linestrings$ID)
  return(extents)
}
apply_tracks <- function(L,pts,m,typ){
  tracks <- pts  |>
    mutate(geometry_lead = lead(geometry, default = NULL)) |>
    # drop the NA row created by lagging
    slice(-n()) |>
    mutate(line = st_sfc(purrr::map2(
      .x = geometry,
      .y = geometry_lead,
      .f = ~{st_union(c(.x, .y)) |> suppressWarnings(st_cast("LINESTRING"))}
    ))) |>
    st_drop_geometry()|>
    select(-geometry_lead) |>
    rename(geometry="line")|>
    st_set_geometry("geometry")|>
    st_set_crs(4326)|>
    st_cast("LINESTRING")
  pt <- pts |> slice(L)
  options(dplyr.summarise.inform = FALSE)
  id=unique(pt$ID)
  name=pt$NAME
  year=pt$SEASON
  dte = pt$date
  basin=pt$BASIN
  distmods <- m$distmods
  emod=m$emod
  rocimod=m$rocimod
  pocimod=m$pocimod
  minpresss=m$minpress
  ###  north american (NA) basins may return NA instead of  "NA"
  if(is.na(basin)){basin="NA"}
  X <- pt$LON
  Y <- pt$LAT
  ###  create a custom CRS centered on the track. This helps preserve true distances and reduces error
  ##  associated with crs warping
  custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",X," +lat_0=",Y," +datum=WGS84")  ##  adding the WGS84 datum to avoid reprojection issues when no internet connection is available for PROJ grib libraries
  ###  create another CRS for the whole set of tracks for this storm
  custCRSall <- paste0("+proj=laea +lat_0=",round(mean(pts$LAT,na.rm=TRUE))," +lon_0=",round(mean(pts$LON,na.rm=TRUE))," +lat_1=",ceiling((max(pts$LAT,na.rm=TRUE)-min(pts$LAT,na.rm=TRUE))/6),
                       " +lat_2=",ceiling(5*(max(pts$LAT,na.rm=TRUE)-min(pts$LAT,na.rm=TRUE))/6)," +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
  ##  adding the WGS84 datum to avoid reprojection issues when no internet connection is available for PROJ grib libraries
  ## get the center point of the current track
  center <- st_sfc(st_point(c(X,Y),dim="XY"),crs = 4326) |>
    st_transform(custCRS)
  ###  take relevant information from the IBTrACS data
  Cat <- as.numeric(pt$USA_SSHS)
  maxwind <- round(as.numeric(pt$WIND))
  if (length(maxwind)==0){maxwind <- NA}
  ###  in rare, mostly minor storms, the maxwind is missing
  ##   assign this based on the category
  if(is.na(maxwind)&!is.null(m)){
    maxwind=ifelse(Cat==5,137,ifelse(Cat==4,113,ifelse(Cat==3,96,ifelse(Cat==2,83,ifelse(Cat==1,64,ifelse(Cat==0,34,ifelse(Cat==-1,24,ifelse(Cat==-2,14,ifelse(Cat==-3,10,NA)))))))))
  }
  maxwinddist <- as.numeric(pt$RMW)
  if (length(maxwinddist)==0){maxwinddist <- NA}
  minpress <- as.numeric(pt$PRES)
  if (length(minpress)==0){minpress <- NA}
  #### in only a very few cases, the minimum pressure is missing,
  ####  take the average for the storms of the same size class in the same month and in the same basin
  if (is.na(minpress)&!is.null(m)){
    minpress=minpresss$PRES[minpresss$USA_SSHS==Cat&minpresss$MONTH==pt$MONTH&minpresss$BASIN==basin]
    if (length(minpress)==0){
      minpress=mean(minpresss$PRES[minpresss$MONTH==pt$MONTH&minpresss$BASIN==basin],na.rm=TRUE)
    }
  }
  maxpress <- as.numeric(pt$POCI)
  if (length(maxpress)==0){maxpress=NA}
  if (is.na(maxpress)) maxpress <- mean(pts$POCI,na.rm=TRUE)
  ### there are more cases of  missing outer pressure, so model this
  if (is.na(maxpress)&!is.null(m)){
    if (length(pocimod$model[pocimod$USA_SSHS==Cat&pocimod$BASIN==basin])>0){
      maxpress <- minpress+predict(pocimod$model[pocimod$USA_SSHS==Cat&pocimod$BASIN==basin][[1]],newdata=data.frame(WIND=maxwind))
    }else{
      maxpress <- minpress
    }
  }
  ## for some reason, IBTracs sometimes lists the POCI as lower than the central pressure. Likely these are no longer 'cyclones' but dissipated
  if (maxpress<=minpress&!is.null(m)) maxpress=minpress+1
  ROCId <- as.numeric(pt$ROCI)
  if (length( ROCId)==0){ ROCId=NA}
  EYE <- as.numeric(pt$EYE)
  if (length(EYE)==0){ EYE=NA}
  swathpoints <- c(NE34 = pt$R34_NE,SE34 = pt$R34_SE,NW34 = pt$R34_NW,SW34 = pt$R34_SW,
                   NE50 = pt$R50_NE,SE50 = pt$R50_SE,NW50 = pt$R50_NW,SW50 = pt$R50_SW,
                   NE64 = pt$R64_NE,SE64 = pt$R64_SE,NW64 = pt$R64_NW,SW64 = pt$R64_SW)
  swathpoints <- sapply(swathpoints,as.numeric)
  swathpoints[is.na(swathpoints)] <-0
  maxswath <- setNames(rep(maxwinddist,4),paste0(c("NE","SE","SW","NW"),maxwind))
  ###  store the native given information
  swathorig <- list(swathpoints=swathpoints,maxswath=maxswath[!names(maxswath)%in%names(swathpoints)])
  swathsource <- rep(pt$source,length(c(swathorig$swathpoints,swathorig$maxswath)))
  if (!is.null(m)){
    swathnew <- model_extent(swathorig,m,maxwind,maxwinddist,Cat,basin)
    ###  sometimes the modeled extents are greater than the ROCId, which doesnt really make sense
    swathnew$dist[swathnew$dist>ROCId] <- ROCId
    #swathcomp <- merge(data.frame(val_new=swathnew,name=names(swathnew)),
    #                   data.frame(val_orig=c(swathorig$swathpoints,swathorig$maxswath),source=swathsource ,
    #                              name=c(names(swathorig$swathpoints),names(swathorig$maxswath))),by="name",all.x=TRUE)
    #swathcomp$source[!is.na(swathcomp$source)] <- paste(swathcomp$source[!is.na(swathcomp$source)],"& modeled")
    swathnew$source <- paste(unique(swathsource),"&",swathnew$source)
    swathpoints <- setNames(swathnew$dist,paste0(swathnew$quad,swathnew$windspeed))
    #swathpoints <- setNames(swathcomp$val_new,swathcomp$name)
    swathsource <- setNames(swathnew$source,paste0(swathnew$quad,swathnew$windspeed))
    #if (!all(is.na(swathcomp$val_orig))&&any(swathcomp$val_new!=swathcomp$val_orig)) browser()
  }else{
    swathpoints <- c(swathorig$swathpoints,swathorig$maxswath)
  }
  ## sometimes the combination of swaths and maximum wind creates duplicates, remove those
  swathpoints <- swathpoints[!duplicated(data.frame(val=swathpoints,name=names(swathpoints)))]
  if (all(swathpoints[!is.na(swathpoints)]==0)){
    point <- pt |>
      mutate(kts=NA,location="track points",dist.m=0,rmw.m=maxwinddist*1852,ID=id,name=name,maxwind.ms=as.numeric(maxwind)*0.5144,
             minpress.mb=as.numeric(minpress),
             maxpress.mb=as.numeric(maxpress),roci.m=ROCId*1852,centerX=X,centerY=Y)|>
      select(kts,BASIN,location,source,ID,name,date,dist.m,rmw.m,maxwind.ms,minpress.mb,maxpress.mb,roci.m,centerX,centerY)|>
      rename(basin=BASIN)
    trck <- tracks |>
      filter(date==dte)|>
      mutate(kts=NA,location="track",dist.m=0,rmw.m=maxwinddist*1852,ID=id,name=name,maxwind.ms=as.numeric(maxwind)*0.5144,
             minpress.mb=as.numeric(minpress),
             maxpress.mb=as.numeric(maxpress),roci.m=ROCId*1852,centerX=X,centerY=Y)|>
      select(kts,BASIN,location,source,ID,name,date,dist.m,rmw.m,maxwind.ms,minpress.mb,maxpress.mb,roci.m,centerX,centerY)|>
      rename(basin=BASIN)
    linestrings=bind_rows(trck|>st_transform(custCRSall),point|>st_transform(custCRSall))
    cat(paste("\rBuilding wind extents for ",unique(linestrings$ID)," : %",round(100*L/nrow(pts),1)))
    return(linestrings)
  }
  if (is.na(maxwinddist)) maxwinddist <- min(swathpoints,na.rm=TRUE)
  for (i in which(swathpoints>0)){
    swath <- matrix(c(0,0,
                      0,swathpoints[i]*1852,
                      swathpoints[i]*cos(80*pi/180)*1852,swathpoints[i]*sin(80*pi/180)*1852,
                      swathpoints[i]*cos(70*pi/180)*1852,swathpoints[i]*sin(70*pi/180)*1852,
                      swathpoints[i]*cos(60*pi/180)*1852,swathpoints[i]*sin(60*pi/180)*1852,
                      swathpoints[i]*cos(50*pi/180)*1852,swathpoints[i]*sin(50*pi/180)*1852,
                      swathpoints[i]*cos(40*pi/180)*1852,swathpoints[i]*sin(40*pi/180)*1852,
                      swathpoints[i]*cos(30*pi/180)*1852,swathpoints[i]*sin(30*pi/180)*1852,
                      swathpoints[i]*cos(20*pi/180)*1852,swathpoints[i]*sin(20*pi/180)*1852,
                      swathpoints[i]*cos(10*pi/180)*1852,swathpoints[i]*sin(10*pi/180)*1852,
                      swathpoints[i]*1852,0,
                      0,0),ncol=2,byrow=T)
    swath <- st_sfc(st_polygon(list(swath)),crs=custCRS)
    swath <- rot(swath,names(swathpoints)[i])
    swath <- st_as_sf(swath,crs=custCRS)
    swath$location <- substr(names(swathpoints[i]),1,2)
    swath$kts <- as.numeric(gsub("NE|SE|SW|NW","",names(swathpoints)[i]))
    swath$dist_m <- swathpoints[i]*1852
    swath$source <- swathsource[i]
    if (exists("swaths", inherits=FALSE))
      swaths <- rbind(swaths,swath)
    else
      swaths <- swath
    # }
  }
  #######
  ####  Handle the eye wall
  #######
  if(!is.na(EYE)){
    eyerad <- EYE/2
    eyesource <- unique(pt$source)
  }else{
    eyerad <- NA
    if (!is.null(m)){
      ##  one approach is 80% of the minimum wind distance
      #eyerad=min(swathpoints[swathpoints>0])*.8
      ## another approach is to use the models
      if (Cat>=0){
        eyerad <- data.frame(minwinddist=min(swathpoints[swathpoints>0],na.rm=TRUE),USA_SSHS=Cat,BASIN=basin) |>
          group_nest(USA_SSHS,BASIN) |>
          left_join(emod|>select(USA_SSHS,BASIN,model),by=c("USA_SSHS","BASIN")) |>
          mutate(pred=map2(model,data,predict))|>
          select(-model) |>
          tidyr::unnest(c(data,pred))|>
          mutate(pred=ceiling(pred))|>pull(pred)
        ###  decrease the radius a little, need to adjust slightly for the gradual increase in wind speeds?
        eyerad <- eyerad*0.9
      }
      eyesource=paste0(unique(pt$source)," & modeled")
    }
    else{eyerad=NA;eyesource=unique(pt$source)}
  }
  if (!is.na(eyerad)){
    eye = st_buffer(center,eyerad*1852) |>
      st_as_sf() |>
      mutate(location="eye",kts=0,dist_m=eyerad*1852,source=eyesource)
  }
  #######
  ### Handle the ROCI
  #######
  if (is.na(ROCId)){
    if (!is.null(m)){
      rocisource=paste0(unique(pt$source)," & modeled")
      if (Cat>0){
        ###  the ROCI models are only good for well-formed storms
        ROCId <- data.frame(USA_SSHS=Cat,BASIN=basin,dist=max(swathpoints)) |>
          group_nest(USA_SSHS,BASIN) |>
          left_join(rocimod|>select(USA_SSHS,BASIN,model),by=c("USA_SSHS","BASIN")) |>
          mutate(pred=map2(model,data,predict))|>
          select(-model) |>
          tidyr::unnest(c(data,pred))|>
          mutate(pred=round(pred))|>pull(pred)
        if (ROCId<max(swathpoints,na.rm=TRUE))  ROCId=max(swathpoints,na.rm=TRUE)+10
        ###  otherwise, just assume they are a little more than the maximum wind extent
      }else{
        ROCId=max(swathpoints,na.rm=TRUE)*1.1
     }
      ####  Otherwise the ROCI is assumed to be the maximum wind extent
      ###  shouldnt assume ROCI if modelling is not desired
      ###  can set to 500km if necessary in later steps
    #}else{
      #rocisource=unique(pt$source)
    #  ROCId <- max(swathpoints)+10
    }
  }else{rocisource=unique(pt$source)}
  if (!is.na(ROCId)){
    ROCI = st_buffer(center,ROCId*1852) |>
      st_as_sf()|>
      mutate(location="ROCI",kts=0,dist_m=ROCId*1852,source=rocisource)
    swaths <- bind_rows(swaths,ROCI)
  }


  if (!is.na(eyerad)) swaths <- rbind(swaths,eye)
  swaths <- swaths |>
    group_by(kts,location) |>
    summarise(dist.m = mean(dist_m),
              source=unique(source)) |>
    ungroup() |>
    st_cast()#|>
  swaths$ID <- id;
  swaths$basin <- basin;
  swaths$name <- name;
  swaths$date <- dte;
  swaths$maxwind.ms <- as.numeric(maxwind)*0.5144;
  swaths$rmw.m <- maxwinddist*1852
  swaths$minpress.mb <- as.numeric(minpress);
  swaths$maxpress.mb <- as.numeric(maxpress);
  swaths$roci.m <- ROCId*1852;
  swaths$centerX <- X;
  swaths$centerY <- Y;
  swaths <- swaths |>
    ungroup() |>
    suppressWarnings(st_cast("POLYGON")) |>
    rename(geometry=x)|>
    st_set_geometry("geometry")

  point <- pt |>
    #st_transform(custCRSall)|>
    mutate(kts=NA,location="track points",dist.m=0,rmw.m=maxwinddist*1852,ID=id,name=name,maxwind.ms=as.numeric(maxwind)*0.5144,
           minpress.mb=as.numeric(minpress),
           maxpress.mb=as.numeric(maxpress),roci.m=ROCId*1852,centerX=X,centerY=Y)|>
    select(kts,BASIN,location,source,ID,name,date,dist.m,rmw.m,maxwind.ms,minpress.mb,maxpress.mb,roci.m,centerX,centerY)|>
    rename(basin=BASIN)
  trck <- tracks |>
    filter(date==dte)|>
    #st_transform(custCRSall)|>
    mutate(kts=NA,location="track",dist.m=0,rmw.m=maxwinddist*1852,ID=id,name=name,maxwind.ms=as.numeric(maxwind)*0.5144,
           minpress.mb=as.numeric(minpress),
           maxpress.mb=as.numeric(maxpress),roci.m=ROCId*1852,centerX=X,centerY=Y)|>
    select(kts,BASIN,location,source,ID,name,date,dist.m,rmw.m,maxwind.ms,minpress.mb,maxpress.mb,roci.m,centerX,centerY)|>
    rename(basin=BASIN)
  if ("linestrings" %in% typ|"all" %in% typ){
    linestrings <- stdh_cast_substring(swaths)
    linestrings <- bind_rows(linestrings|>st_transform(custCRSall),trck|>st_transform(custCRSall),point|>st_transform(custCRSall))
  }
  if ("polygons" %in% typ|"all" %in% typ){
    swaths <- bind_rows(swaths|>st_transform(custCRSall),trck|>st_transform(custCRSall),point|>st_transform(custCRSall))
  }
  cat(paste("\rBuilding wind extents for ",unique(linestrings$ID)," : %",round(100*L/nrow(pts),1)))
  if (("linestrings" %in% typ)&length(typ==1)) extents=linestrings
  else if (("polygons" %in% typ)&length(typ==1)) extents=swaths
  else if (("linestrings" %in% typ & "polygons" %in% typ)|typ=="all")  extents <- list(swaths=swaths,linestrings=linestrings)#|> tidyr::nest(.by =ID)
  return(extents)
}
interp_track <- function(pnts,tres,pre=NULL,wind=TRUE){
  if (unique(st_geometry_type(pnts))=="POINT"){
    lines <- pnts  |>
      mutate(geometry_lead = lead(geometry, default = NULL)) |>
      # drop the NA row created by lagging
      slice(-n()) |>
      mutate(line = st_sfc(purrr::map2(
        .x = geometry,
        .y = geometry_lead,
        .f = ~{st_union(c(.x, .y)) |> suppressWarnings(st_cast("LINESTRING"))}
      ))) |>
      st_drop_geometry()|>
      select(-geometry_lead) |>
      rename(geometry="line")|>
      st_set_geometry("geometry")|>
      st_set_crs(4326)|>
      st_cast("LINESTRING") |> rename(centerX=LON,centerY=LAT)
    interpt=TRUE
  }else{
    lines=pnts
    interpt=FALSE
  }
  ####  mark the quadrant extent start and end points, dont want to interpolate between them
 if (wind) pnts <- pnts |>
    mutate(across(R34_NE:R64_NW, ~if_else(.x==0,NA,.x)))|>
    mutate(across(R34_NE:R64_NW, function(x){
      res <- rle(is.na(x))
      ends <- cumsum(res$lengths)
      starts <- c(1, ends[-length(ends)] + 1)
      missing_runs <- data.frame(
        start = starts[res$values],
        end = ends[res$values]
      )
      missing_runs <- as.vector(t(missing_runs))
      replace(x,missing_runs,-999)
    }))
  t_res_times <- seq(min(pnts$date,na.rm=TRUE),max(pnts$date,na.rm=TRUE),by=tres*60)
  newpoints <- lapply(1:nrow(lines),function(n,tr,trests,pts){
    ## create custom CRS centered on segment to minimize geometry calculation errors
    custCRS <- paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=",tr$centerX[n]," +lat_0=",tr$centerY[n])
    t1 <- tr$date[n]
    if (n==nrow(tr)){
      t2 <- pts$date[n+1]
      t_res_ts <- trests[trests>=t1&trests<=t2]
    }else{
      t2 <- tr$date[n+1]
      t_res_ts <- trests[trests>=t1&trests<t2]
    }
    ###  get the time difference between the two adjacent segments
    dt_track = as.numeric(difftime(t2,t1,units="hours"))
    # }
    ###  only necessary to interpolate if they are more than the desired temporal resolution
    if (length(t_res_ts)>0){
      ###  get the time difference between the interpolated times
      dt <- as.numeric(difftime(t_res_ts,t1,units="hours"))
      ##  the intervals are the cumulative proportion of the segment as calculated by the
      ##  dt time inervals relative to the duration of the entire segment
      ints <- dt[!is.na(dt)]/dt_track
      ###  now sample along the track segment according to the intervals calculated
      ###  transform back to original crs
      t  <- tr[n,] |>
        st_transform(custCRS) |>
        st_line_sample(type="regular",sample=ints) |>
        st_as_sf(crs=custCRS) |>
        st_transform(st_crs(tr))|>
        st_cast("POINT")|>
        mutate(date=t_res_ts)
      sf::st_geometry(t)="geometry"
      t
    }
  },tr=lines,trests=t_res_times,pts=pnts)
  newpoints <- do.call(rbind,newpoints)
  if (!is.null(pre)){
    pre=max(pre)
    newpoints <- bind_rows(newpoints |> slice(rep(1,pre/(tres/60)))|>
                             mutate(date=seq(newpoints$date[1]-pre*60*60,newpoints$date[1]-tres*60,by=tres*60)),
                           newpoints)
    }
  newpoints$LON <- st_coordinates(newpoints)[,1]
  newpoints$LAT <- st_coordinates(newpoints)[,2]
  newpoints <- newpoints |> left_join(st_drop_geometry(pnts)|>select(-any_of(c("LON","LAT"))),by=join_by(date))
  if (interpt) newpoints <- newpoints |>
    mutate(across(c(SEASON,LAT,LON,USA_SSHS:R64_NW),as.numeric)) |>
    ###  zeros in the quadrant data shouldnt be interpolated as zero, they should be missing
    ###  similarly, if a storm loses a wind extent and regains it later, the values in between should not be interpolated
    mutate(across(STORM_SPEED:R64_NW,zoo::na.approx,na.rm=FALSE))|>
    mutate(across(R34_NE:R64_NW,~if_else(.x<0,NA,.x)),
      source=if_else(is.na(NAME),"interpolated","native"),
      MONTH=format(date,"%m"),
      USA_SSHS = if_else(is.na(USA_SSHS),
                         if_else(WIND<=33,-1,
                                 if_else(WIND<=63,-0,
                                         if_else(WIND<=82,1,
                                                 if_else(WIND<=95,2,
                                                         if_else(WIND<=112,3,
                                                                 if_else(WIND<=136,4,5)))))),USA_SSHS))|>
    tidyr::fill(c(SID:BASIN,ID))
  newpoints
}
model_extent <- function(swathpnts,m,mxwnd,mxwnd_dist,C,bas){
  maxswath <- swathpnts$maxswath
  swathpnts <- swathpnts$swathpoints
  mod <- m$distmods
  emod <- m$emod
  pocimod <- m$pocimod
  rocimod <- m$rocimod
  minpresss <- m$minpress
  ###  first, develop a vector of wind speeds to model for distance
  swathpnts_blank <- data.frame(windspeed=rep(c(34,50,64,mxwnd),4),quad=rep(c("NE","SE","SW","NW"),each=4),USA_SSHS=C,BASIN=bas,source="modeled")
  swathpnts_blank <- swathpnts_blank[swathpnts_blank$windspeed<=mxwnd,]
  swathpnts_blank <- swathpnts_blank[!duplicated(swathpnts_blank),]
  ###  join the available distances
  repl <- swathpnts[grep(max(as.numeric(gsub("NE|SE|SW|NW","",names(swathpnts[swathpnts>0]))),na.rm=T),names(swathpnts))]
  repl <- repl[repl>0]
  if (!is.na(mxwnd_dist)){
    ####  if the RMW is already determined,
    ####  insert the maximum wind and its distance
    diff <- repl-mxwnd_dist
    repl <- replace(repl,1:length(repl),mxwnd_dist)
    ###  sometimes the maximum wind distance is greater than that supplied by the quadrant data
    ##  this is either due to asymmetry or sometimes seemingly to error
    ##  in these cases, make the maximum wind distance proportionally shorter than the next known highest wind distance
    # if modelling is desired
    repl[diff<=0] <- as.numeric(gsub('\\D+', "", names(repl[diff<=0])))/mxwnd*repl[diff<=0]
    names(repl) <- suppressWarnings(gsub(paste(gsub("NE|SE|SW|NW","",names(repl)),collapse="|"),mxwnd,names(repl)))
    swathpnts <- c(swathpnts,c(repl))
  }
  swathpnts <- swathpnts[swathpnts>0]
  swathpnts <- data.frame(windspeed=names(swathpnts),dist=swathpnts)|>
    mutate(quad=substr(windspeed,1,2),
           windspeed=as.numeric(gsub(paste(quad,collapse="|"),"",windspeed)))
  swathpnts <- swathpnts_blank |>
    left_join(swathpnts,by=join_by(quad,windspeed))
  ###  model the missing distances
  ###  models could only be fit for category 0 storms and above
  #if (C>0){
    swathpnts <- swathpnts |>
      group_nest(quad,USA_SSHS,BASIN) |>
      left_join(mod|>select(quad,USA_SSHS,BASIN,model_asymp),by=c("USA_SSHS","quad","BASIN")) |>
      left_join(mod|>filter(quad=="all")|>select(USA_SSHS,BASIN,model_asymp),by=c("USA_SSHS","BASIN"))|>
      left_join(mod|>filter(quad=="all",BASIN=="all")|>select(USA_SSHS,model_asymp)|>rename(model_asymp.z=model_asymp),by=c("USA_SSHS"))|>
      bind_cols(mod|>filter(quad=="all",BASIN=="all",USA_SSHS==999)|>select(model_asymp))|>
      mutate(pred.x = map2(model_asymp.x,data,possibly(predict,otherwise=NA_real_)),
             pred.y = map2(model_asymp.y,data,possibly(predict,otherwise=NA_real_)),
             pred.z = map2(model_asymp.z,data,possibly(predict,otherwise=NA_real_)),
             pred = map2(model_asymp,data,possibly(predict,otherwise=NA_real_))) |>
      tidyr::unnest(c(data,pred.x,pred.y,pred.z,pred))|>
      mutate(dist_pred = round(coalesce(pred.x,pred.y,pred.z,pred)))
    ##  for smaller storms, use the cat 1 storm model
  # }else{
  #   swathpnts <- swathpnts |>
  #     group_nest(quad,BASIN) |>
  #     left_join(mod|>filter(USA_SSHS==1)|>select(quad,model_asymp,BASIN),by=c("quad","BASIN")) |>
  #     mutate(dist=map2(model_asymp,data,predict))|>
  #     select(-model_asymp) |>
  #     tidyr::unnest(c(data,dist))|>
  #     mutate(dist=round(dist))
  # }
  swathpnts$source[!is.na(swathpnts$dist)] <- "native"
  swathpnts$dist[is.na(swathpnts$dist)] <- swathpnts$dist_pred[is.na(swathpnts$dist)]
  ###  make sure a lower windspeed in the same quadrant is not the lower than a higher windspeed
  swathpnts <- swathpnts |>
    arrange(rev(windspeed))|>
    group_by(quad)|>
    mutate(dist = if_else(windspeed==max(windspeed,na.rm=TRUE),dist,if_else(dist<=lag(dist),lag(dist)*1.1,dist)))

  #swathpnts[swathpnts<10] <-10
  swathpnts
}
rot = function(ls,q){
  if (grepl("SE",q)){
    a=pi/2
  }else if (grepl("SW",q)){
    a=pi
  }else if (grepl("NW",q)){
    a=1.5*pi
  }else{return(ls)}
  mat <- matrix(c(round(cos(a)), round(sin(a)), -round(sin(a)), round(cos(a))), 2, 2)
  ls*mat
}

stdh_cast_substring <- function(x, to = "MULTILINESTRING") {
  ####  if the geometry was sent after the extents are already produced
  ###  separate the polygons first
  ## https://rpubs.com/dieghernan/Cast-line-subsegments-R
  if (!all(unique(st_geometry_type(x)) %in% c("POLYGON", "LINESTRING"))) {
    yg <- x[!st_geometry_type(x)=="POLYGON",]
    x <- x[st_geometry_type(x)=="POLYGON",]
    #stop("Input should be  LINESTRING or POLYGON")
  }
  ggg <- st_geometry(x)
  for (k in 1:length(st_geometry(ggg))) {
    sub <- ggg[k]
    geom <- lapply(
      1:(length(st_coordinates(sub)[, 1]) - 1),
      function(i)
        rbind(
          as.numeric(st_coordinates(sub)[i, 1:2]),
          as.numeric(st_coordinates(sub)[i + 1, 1:2])
        )
    ) |>
      st_multilinestring() |>
      st_sfc()

    if (k == 1) {
      endgeom <- geom
    }
    else {
      endgeom <- rbind(endgeom, geom)
    }
  }
  endgeom <- endgeom |> st_sfc(crs = st_crs(x))
  if (class(x)[1] == "sf") {
    endgeom <- st_set_geometry(x, endgeom)
  }
  if (to == "LINESTRING") {
    endgeom <- endgeom |> st_cast("LINESTRING")
  }
  endgeom <- suppressWarnings(endgeom |> group_by(location,kts)|> st_cast("LINESTRING"))
  coords <- st_coordinates(endgeom)
  origin__seg <- unique(coords[rowSums(coords[,1:2]==0)==2,"L1"])
  endgeom <- endgeom |>
    ungroup()|>
    slice(-origin__seg)|>
    group_by(across(-geometry)) |>
    summarise()|>
    st_line_merge()
  if (exists("yg",inherits = FALSE)) endgeom <- rbind(endgeom,yg)
  return(endgeom)
}
