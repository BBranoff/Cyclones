get_noaaWL <- function(storm,dpath=NULL){
  if(is.null(dpath)) dpath=tempdir()
  if (!dir.exists(paste0(dpath,"/NOAA/"))) dir.create(paste0(dpath,"/NOAA"))
  dpath <- paste0(dpath,"/NOAA/")
  fname = paste0(dpath,unique(storm$ID),"_waterlevels.csv")
  if (file.exists(fname)) return(st_as_sf(read.csv(fname),coords=c("Lon","Lat"),crs=4326))
  stormdates <- as.POSIXct(storm$ISO_TIME)
  ###  get station list
  max_attempts =5;attempt=0;success=FALSE
  if (!any(grepl("NOAAstats",list.files(tempdir())))){
    while(attempt < max_attempts && !success) {
      attempt <- attempt + 1
      tryCatch({
        statdts <- jsonlite::fromJSON("https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json?expand=details&type=historicwl")
        saveRDS(statdts,tempfile(pattern="NOAAstats",fileext = ".rds"))
        success <- TRUE
      }, error = function(e){
        if (attempt >= max_attempts) {
          stop(paste0("could not get NOAA data: ",e))
        }else{
          Sys.sleep(0.5)
        }
      })
    }
  }
  NOAAstat_temp <- grep("NOAAstats",list.files(tempdir(),full.names=TRUE),value=TRUE)
  statdets <- readRDS(NOAAstat_temp)
  statdets <- statdets$stations
  statdets <- st_as_sf(statdets,coords=c("lng","lat"),crs=4326)
  storm <- st_as_sf(storm,coords=c("LON","LAT"),crs=4326)
  roci <- storm |> pull(contains("CONS_ROCI"))|>mean(na.rm=TRUE)
  if (is.na(roci)) roci <- 500000/1852
  extnt <- st_transform(st_union(st_geometry(st_buffer(storm,as.numeric(roci)*1852))),4326)
  suppressMessages({sf_use_s2(FALSE)
  statdets_storm <- statdets[which(lengths(st_intersects(statdets,extnt))>0),]
  sf_use_s2(TRUE)})
  statdets_storm$details$removed[statdets_storm$details$removed==""] <- format(as.POSIXct(Sys.time(),tz="UTC"),"%Y-%m-%d %H:%M:%S")
  statdets_storm$details$removed <- as.POSIXct(statdets_storm$details$removed,tz="UTC")
  statdets_storm <- statdets_storm[!is.na(statdets_storm$details$established),]
  statdets_storm$details$established[statdets_storm$details$established==""] <- statdets_storm$details$origyear[statdets_storm$details$established==""]
  statdets_storm$details$established<- as.POSIXct(statdets_storm$details$established,tz="UTC")
  statdets_storm <- statdets_storm[statdets_storm$details$established<min(stormdates),]
  statdets_storm <- statdets_storm[statdets_storm$details$removed>max(stormdates),]
  if (length(statdets_storm$details$id)==0) {warning("No NOAA data available for this storm."); return(NULL)}
  urls <- paste0("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?begin_date=",format(min(stormdates)-365*24*60*60,"%Y%m%d"),"&end_date=",format(max(stormdates),"%Y%m%d"),"&station=",statdets_storm$details$id,
                 "&product=daily_max_min&datum=MSL&time_zone=gmt&units=metric&application=USFShurricaneResearch&format=csv")
  stationdat <- lapply(1:length(urls), function(i) {
    statid <- gsub(".*station=(.+)&product.*", "\\1", urls[i])
    #chck <- readLines(paste0("https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/",statid,"/details.json"))
    #startdt <- chck[grep("established",chck)]
    #startdt <- gsub('  \"established\": \"|\",',"",startdt)
    #startdt <- as.POSIXct(startdt,format="%Y-%m-%d %H:%M:%S")
    #enddt <- chck[grep("removed",chck)]
    #enddt <- gsub('  \"removed\": \"|\",',"",enddt)
    #enddt <- as.POSIXct(enddt,format="%Y-%m-%d %H:%M:%S")
    #if (min(stormdates)>startdt&max(stormdates)<enddt){
      dat <- tryCatch(read.csv(urls[i]),error=function(e)e)
      if (!"message" %in% names(dat)&&nrow(dat)>1){
        #dat <- st_as_sf(dat,geometry=statdets_storm$geometry[statdets_storm$details$id==statid])
        dat$Lat = unique(st_coordinates(statdets_storm[statdets_storm$details$id==statid,])[,2])
        dat$Lon = unique(st_coordinates(statdets_storm[statdets_storm$details$id==statid,])[,1])
        dat$storm=unique(storm$ID)
        #dat$station=sts$id[i]
        return(dat)
      }
    #}
    if (i%%10==0){ Sys.sleep(1);cat(paste0("\rDownloading NOAA Tides & Currents Water Levels: %",round(100*i/length(urls))))}
  })
  cat("\r", paste(rep(" ", 50), collapse=""), "\r")
  stationdat <- do.call(rbind,stationdat)
  if (is.null(stationdat)) {warning("No NOAA data available for this storm."); return(NULL)}
  stationdat <- stationdat |>
    group_by(stationId,date)|>
    filter(value==max(value,na.rm=TRUE)) |>
    group_by(stationId) |>
    mutate(date=as.POSIXct(paste(date,time),format="%Y-%m-%d %H:%M", tz="UTC"),
           meanHT=mean(value[date<min(stormdates,na.rm=TRUE)],na.rm=TRUE))|>
    filter(date>=min(stormdates))
  ###  find mean high tide and then filter to only storm dates
  write.csv(stationdat,fname)
  return(st_as_sf(read.csv(fname),coords=c("Lon","Lat"),crs=4326))
}

