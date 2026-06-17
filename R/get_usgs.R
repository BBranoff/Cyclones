get_usgs <- function(storm,key=NULL,dpath=NULL){
  if(is.null(dpath)) dpath=tempdir()
  if (!dir.exists(paste0(dpath,"/USGS/"))) dir.create(paste0(dpath,"/USGS"))
  dpath <- paste0(dpath,"/USGS/")
  fname = paste0(dpath,unique(storm$ID),"_waterlevels.csv")
  if (file.exists(fname)) return(st_as_sf(read.csv(fname),coords=c("lon","lat"),crs=4326))
  stormdates <- as.POSIXct(storm$ISO_TIME)
  storm <- st_as_sf(storm,coords=c("LON","LAT"),crs=4326)
  roci <- storm |> pull(contains("CONS_ROCI"))|>mean(na.rm=TRUE)
  if (is.na(roci)) roci <- 500000/1852
  extnt <- st_transform(st_union(st_geometry(st_buffer(storm,as.numeric(roci)*1852))),4326)
  bbox <- st_bbox(extnt)#st_bbox(c(-67.34398,17.88728,-65.3330,18.52592))#PR
  ##  using dataRetrieval pacakge
  statdets <- jsonlite::fromJSON(paste0("https://api.waterdata.usgs.gov/ogcapi/v0/collections/monitoring-locations/items?bbox=",
                                        paste(bbox,collapse=","),"&limit=10000&f=json&site_type=Estuary&api_key=",key))
  statdets <- statdets$features
  if (length(statdets)==0){warning("No USGS data available for this storm."); return(NULL)}
  #statdets_missingdatum <- statdets$properties$id[is.na(statdets$properties$vertical_datum)]
  if (length(unique(statdets$id))>75) ids <- split(unique(statdets$id),cut(seq_along(unique(statdets$id)), ceiling(length(unique(statdets$id))/75), labels = FALSE))
  else ids <- unique(statdets$id)
  statmet <- lapply(seq_along(ids), function(i) {
    u <- paste0("https://api.waterdata.usgs.gov/ogcapi/v0/collections/time-series-metadata/items?monitoring_location_id=",
                paste(unlist(ids[i]),collapse=","),"&parameter_code=00065&limit=50000&f=json&api_key=",key)
    max_attempts =5;attempt=0;success=FALSE
    met <- while(attempt < max_attempts && !success) {
      attempt <- attempt + 1
      tryCatch({
        mt <- jsonlite::fromJSON(u)
        success <- TRUE
        return(mt)
      }, error = function(e){
        if (attempt >= max_attempts) {
          warning(paste0("data for at least one USGS station could not be downloaded: ",e))
          return(NULL)
        }else{
          Sys.sleep(0.5)
        }
      })
    }
    cat("\rfetching USGS station meta data: %",round(100*i/length(ids),1))
    met
  })
  statmet <-  Filter(Negate(is.null), statmet)
  if (length(statmet)==0) {warning("No USGS data available for this storm."); return(NULL)}
  statmetfeat <- lapply(statmet, function(x) f= x$features)
  statmetfeat <- statmetfeat|>bind_rows()
  statmetfeat <- statmetfeat|>filter(as.POSIXct(statmetfeat$properties$begin)<min(stormdates)&as.POSIXct(statmetfeat$properties$end)>max(stormdates))
  statmetfeat <- statmetfeat[!is.na(statmetfeat$geometry$type),]
  if (nrow(statmetfeat)==0) {warning("No USGS data available for this storm."); return(NULL)}
  statmetfeat$lon <- unlist(sapply(statmetfeat$geometry$coordinates,"[[",1))
  statmetfeat$lat <- unlist(sapply(statmetfeat$geometry$coordinates,"[[",2))
  statmetfeat$monitoring_location_id <- statmetfeat$properties$monitoring_location_id
  ids <- unique(statmetfeat$monitoring_location_id)
  dailies<- lapply(seq_along(ids), function(x){
    stat_ids = unique(statmetfeat$properties$statistic_id[statmetfeat$properties$monitoring_location_id==ids[x]])
    if (any(c("00001","00002","00003","00021","00024","00022","00023") %in% stat_ids)){
      ur <- paste0("https://api.waterdata.usgs.gov/ogcapi/v0/collections/daily/items?monitoring_location_id=",
                   #"USGS-02326550&time=",
                   ids[x],"&time=",
                   paste(format(min(stormdates)-365*60*60*24,"%Y-%m-%dT%H:%M:%SZ"),format(max(stormdates),"%Y-%m-%dT%H:%M:%SZ"),sep="/"),
                   "&parameter_code=00065&limit=50000&f=json&api_key=",key)
    }else{
      ur <- paste0("https://api.waterdata.usgs.gov/ogcapi/v0/collections/continuous/items?monitoring_location_id=",
                   #"USGS-02326550&time=",
                   ids[x],"&time=",
                   paste(format(min(stormdates)-365*60*60*24,"%Y-%m-%dT%H:%M:%SZ"),format(max(stormdates),"%Y-%m-%dT%H:%M:%SZ"),sep="/"),
                   "&parameter_code=00065&limit=50000&f=json&api_key=",key)
    }
    max_attempts =3;attempt=0;success=FALSE
    dat=while(attempt < max_attempts && !success) {
      attempt <- attempt + 1
      tryCatch({
        dt <- jsonlite::fromJSON(ur)
        success <- TRUE
        return(dt)
      }, error = function(e){
        if (attempt >= max_attempts) {
          return(NULL)
        }else{
          Sys.sleep(0.5)
        }
      })
    }
    Sys.sleep(0.1)
    cat("\rDownloading/loading USGS water level data: %",round(100*x/length(ids),1))
    dat
  })
  cat("\r", paste(rep(" ", 50), collapse=""), "\r")
 # }
  dailies <-  Filter(Negate(is.null), dailies)
  dailies_f <- lapply(dailies,function(x) x$features)
  dailies_f <- dailies_f |> bind_rows()
  if (nrow(dailies_f)==0) {warning("No USGS data available for this storm."); return(NULL)}
  dailies_f <- dailies_f$properties |>
    mutate(time =as.POSIXct(time),value = as.numeric(value),day=format(time,"%j"),
           value=if_else(unit_of_measure=="ft",value*.3048,value))|>
    ## get mean high tide
    group_by(monitoring_location_id,day) |>
    mutate(HT = max(value,na.rm=TRUE)) |>
    group_by(monitoring_location_id) |>
    mutate(MSL = mean(value[time<min(stormdates,na.rm=TRUE)],na.rm=TRUE),
           MHT = mean(HT[time<min(stormdates,na.rm=TRUE)],na.rm=TRUE)) |>
    ungroup()|>
    filter(time>=min(stormdates))|>
    group_by(monitoring_location_id,time)|>
    filter(value==max(value,na.rm=TRUE))

  #dailies_f <- dailies_f |> left_join(statdets$properties |> select(id,vertical_datum,altitude),by=join_by(monitoring_location_id==id))
  dailies_f <- dailies_f |> left_join(statmetfeat|> select(monitoring_location_id,lon,lat)|>distinct(),by=join_by(monitoring_location_id))|>
    select(-qualifier) |> ungroup()
  write.csv(data.frame(dailies_f),fname)
  return(st_as_sf(read.csv(fname),coords=c("lon","lat"),crs=4326))

  #ggplot(dailies_f,aes(x=time,y=surge*0.3048,col=monitoring_location_id))+geom_line()+
  #  facet_wrap(~monitoring_location_id)
}
