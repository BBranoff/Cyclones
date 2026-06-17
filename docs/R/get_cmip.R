#' @importFrom sf st_crop st_as_sf st_geometry st_intersection
#' @importFrom ecmwfr wf_request
#' @importFrom dplyr c_across
get_cmip <- function(storm,dpath=NULL,var="waterlevel"){
  if(is.null(dpath)) dpath=tempdir()
  if (!dir.exists(paste0(dpath,"/CMIP6/"))) dir.create(paste0(dpath,"/CMIP6/"))
  dpath <- paste0(dpath,"/CMIP6/")
  YRMNTH <- unique(format(seq(min(storm$ISO_TIME),max(storm$ISO_TIME),by="day"),"%Y %m"))
  roci <- storm |> pull(contains("CONS_ROCI"))|>mean(na.rm=TRUE)
  if (is.na(roci)) roci <- 500000/1852
  if (!class(storm)[1]=="sf") storm <- st_as_sf(storm,coords=c("LON","LAT"),crs=4326)
  e <- st_bbox(st_buffer(storm,as.numeric(roci)*1852))
  ###  first get previous year's averages
  YR <- format(max(storm$ISO_TIME),"%Y")
  fnames <- paste0(dpath,ifelse(YR>2014,"future","historical"),"_tide_actual-value_",YR,"_",c("MSL","MLLW","MHHW"),"_v1.nc")
  ###  get tides
  tides <- cmip_request(fnames,dpath,vars=c("mean_sea_level","annual_mean_of_highest_high_water","annual_mean_of_lowest_low_water"),YR=YR,m=NULL)
  ### now get daily maximums
 if (as.numeric(YR)<=2024){
   WLsurges <- lapply(YRMNTH, function(ym) {
     YR = strsplit(ym," ")[[1]][1]
     m = strsplit(ym," ")[[1]][2]
     fnames = paste0(dpath,"/reanalysis_",c("waterlevel","surge"),"_dailymax_",YR,"_",m,"_v3.nc")
     WLsurge <- cmip_request(fnames,dpath,vars=c("storm_surge_residual","total_water_level"),YR=YR,m=m)
     WLsurge <- WLsurge |>tidyr::pivot_longer(cols=-c("lon","lat","variable"),names_to = "Date")
     WLsurge$Date <- as.Date(as.numeric(WLsurge$Date))
     WLsurge
   })
   WLsurges <- do.call(rbind,WLsurges)|>
     tidyr::pivot_wider(id_cols=c("lon","lat","Date"),names_from = "variable")|>
     left_join(tides |> tidyr::pivot_wider(id_cols=c("lon","lat"),names_from="variable"),by=join_by(lon==lon,lat==lat))
 } else{
   options(warn=1)
   warning("CMIP daily maximum water levels not available past 2024. Returning only tidal average. Surge may not be available if other water level sources are not.")
   options(warn=1)
   WLsurges <- tides |> tidyr::pivot_wider(id_cols=c("lon","lat"),names_from="variable")
 }
  geom <- st_as_sf(WLsurges|>distinct(lon,lat),coords=c("lon","lat"),crs=4326,remove=FALSE)
  e <- st_union(st_geometry(st_buffer(storm,as.numeric(roci)*1852)))
  suppressMessages({sf_use_s2(FALSE)
  geoms <- geom[which(lengths(st_intersects(geom,st_transform(e,4326)))>0),]
  sf_use_s2(TRUE)
  })
  WLsurges <- WLsurges |> filter(paste0(lon,lat) %in% paste0(geoms$lon,geoms$lat))|>
    st_as_sf(coords=c("lon","lat"),crs=4326)
  if (as.numeric(YR)<=2024){
    WLsurges <- WLsurges |>
    group_by(geometry)|>
    filter(surge==max(surge,na.rm=TRUE))|>
    ungroup()
  }

  cat("\r", paste(rep(" ", 50), collapse=""), "\r")
  WLsurges
}
cmip_request <- function(fns,dp,vars,YR,m){
  if (!any(file.exists(fns))){
    ###  for tidal averages
    if ("annual_mean_of_highest_high_water" %in% vars){
      request <- list(
        dataset_short_name = "sis-water-level-change-indicators-cmip6",
        variable =vars,
        derived_variable = "absolute_value",
        #product_type= "reanalysis",
        statistic= "1_year",
        experiment= ifelse(YR<=2014,"historical","future"),
        period= YR,
        target = "TMPFILE")
    }else{
      ###  for daily maximums
      # if (YR>2024){
      #   request <-list(
      #     dataset_short_name = "sis-water-level-change-timeseries-cmip6",
      #     variable = vars,
      #     experiment = "future",
      #     model = "cmcc_cm2_vhr4",
      #     ###  daily max not available after 2024
      #     temporal_aggregation = "10_min",
      #     year = YR,
      #     month = m,
      #     version = "v1",
      #     target = "TMPFILE")
      # }else{
        request <-list(
          dataset_short_name = "sis-water-level-change-timeseries-cmip6",
          variable = vars,
          experiment = "reanalysis",
          temporal_aggregation = "daily_maximum",
          year = YR,
          month = m,
          version = "v3",
          target = "TMPFILE")
      #}

    }
    max_attempts =5;attempt=0
    while(attempt < max_attempts ) {
      attempt <- attempt + 1
      tryCatch({
        wf_request(request)
      }, error = function(e){
        if (attempt >= max_attempts) {
          stop(e)
        }else{
          Sys.sleep(0.5)
        }
      })
    }
    unzip(paste0(tempdir(),"/",request$target,".zip"),exdir = dp)
  }
  valsdf <-lapply(fns,function(fn) {
    vals <- ncdf4::nc_open(fn)
    on.exit(ncdf4::nc_close(vals))
    lat <- ncdf4::ncvar_get(vals,"station_y_coordinate")
    lon <- ncdf4::ncvar_get(vals,"station_x_coordinate")
    var <- names(vals$var)[!names(vals$var) %in% c("station_x_coordinate","station_y_coordinate")]
    val <-  ncdf4::ncvar_get(vals,var)
    valsdf <- data.frame(cbind(lon,lat,val))
    if (!grepl("surge|waterlevel",fn)) names(valsdf) <- c("lon","lat","value")
    else  names(valsdf) <- c("lon","lat",seq.Date(as.Date(paste(YR,m,01,sep="-")),seq(as.Date(paste(YR,as.numeric(m),01,sep="-")),by="month",length.out=2)[2]-1))
    valsdf$variable <- var
    valsdf
  })
  valsdf <- do.call(rbind,valsdf)
  valsdf
}

