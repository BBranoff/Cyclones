#' @importFrom dplyr row_number filter slice_min full_join join_by
#' @importFrom sf st_transform
#' @importFrom terra extract crs distance minmax classify values time ext is.lonlat rotate
#######  function to compare winds stored as the original and/or enhanced IBTrACS linestrings data or the rasters, with the modeled high winds both stored as
########   vectors, with the interpolated raster values
compare_winds <- function(rasts,checksondes=TRUE,ddir=NULL,dt=5){
  cat("Checking input rasters...")
  rasts <- checkrasts(rasts)
  if (!is.null(rasts)){
    storm <- tolower(lapply(strsplit(names(rasts),"_"),function(x)paste(x[1:2],collapse="")))
  }
  cat(paste0("\n",length(storm)," storms found in input: ", paste(unique(storm),collapse=", ")))
  cat("\nChecking for sonde data...")
  sondedata <- lapply(storm, get_sondes, ddir=ddir,checksondes=checksondes)
  names(sondedata) <- storm
  sondedata <- Filter(Negate(isFALSE),sondedata)
  raststorm <- names(sondedata)
  cat(paste0("\n",length(raststorm)," storms matched with sonde data: ", paste(unique(raststorm),collapse=", ")))
  r1 <- rasts[[1]]
  if (is.list(rasts)){
    ###  double subsetting bracket here does not work, unless only one storm
    #rasts <- rasts[grep(paste(raststorm,collapse="|"),gsub("_","",names(rasts)),ignore.case = TRUE)]
    rsamps <- lapply(raststorm, function(st) {
      ###  check for temporal overalap
      ##  whittle down the rasters and sondedata by time
      ##  this greatly reduces the processing time
      sondedata_st <- sondedata[[st]]
      #sondedata_st$time <- as.POSIXct(sondedata_st$time)
      rasts_st <- rasts[[ which(tolower(unlist(lapply(strsplit(names(rasts),"_"),function(x) paste(x[1:2],collapse=""))))==st)]]
      #  cannot call rast on an already loaded raster, this should be taken care of in checkrasts()
      #rasts_st <- rast(rasts_st)
      if (!any(time(rasts_st[[grep("msw|windir",names(rasts_st))]])>=(min(sondedata_st$time,na.rm=TRUE)-dt*60)&
               time(rasts_st[[grep("msw|windir",names(rasts_st))]])<=(max(sondedata_st$time,na.rm=TRUE)+dt*60))){
        warning(paste0("No temporal overalap within specified 'dt' value for storm: ",st))
        return(NULL)
      }
      rasts_st <- rasts_st[[grep("msw|windir",names(rasts_st))]]
      rasts_st <- rasts_st[[time(rasts_st)>=(min(sondedata_st$time,na.rm=TRUE)-dt*60)&
                              time(rasts_st)<=(max(sondedata_st$time,na.rm=TRUE)+dt*60)]]
      cat(paste0("\nSampling at drop sonde locations (<30 m elevation) for ",st))
      rtimes <- time(rasts_st)
      rtimes <- lapply(rtimes,function(x,ts,dt) {
        diffts <- as.numeric(abs(difftime(x,ts,units="mins")))
        sondes <- which(diffts<=dt)
        sondes
      },ts=sondedata_st$time,dt=dt)
      rasts_st <- rasts_st[[which(unlist(lapply(rtimes,length))>0)]]
      sondedata_st <- sondedata_st |> filter(time %in% unique(time)[unique(unlist(rtimes))],
                                             ZW.m<=30)
      ###  if the rasters come from StormR, they are rotated and need to be corrected to overaly with the sonde coordinates
      if (is.lonlat(rasts_st)&&(ext(rasts_st)[1]>=0&ext(rasts_st)[2]>180)) rasts_st <- rotate(rasts_st)
      rast_ex <- extract(rasts_st,st_transform(sondedata_st,crs(rasts_st)))
      names(rast_ex)[2:ncol(rast_ex)] <- paste(time(rasts_st),names(rast_ex)[2:ncol(rast_ex)],gsub("_",",",rep(sources(rasts_st),times=sources(rasts_st,nlyr=TRUE)$nlyr)),sep="_")
      rast_ex <- suppressMessages(sondedata_st |> dplyr::bind_cols(rast_ex)|>
        tidyr::pivot_longer(!c(IX:ID,geometry),names_sep = "_",names_to = c("time_r","method","var","rastfile")) |>
        mutate(time_r = as.POSIXct(time_r,format="%Y-%m-%d %H:%M:%S",tz="UTC"),
               rastfile=gsub(",","_",rastfile),
               percerror = if_else(var=="msw.ms",100*abs((value-WS.ms)/(WS.ms)),100*abs((value-WD)/(WD))))|>
        filter(abs(time_r-time)<=5*60))
      rast_ex$storm <- st
      rast_ex
    })
    names(rsamps) <- raststorm
  }else{
    browser()
    #rsamps <- lapply(rasts,rsampf_f,s=sondedata)
  }
  return(rsamps)
}
rsamp_f <- function(x,s,dt){
  m <- strsplit(names(x),"_")[[1]][1]
  var <- strsplit(names(x),"_")[[1]][2]
  ###  filter sondedata and rasts now so that the time difference is achieved
  sonde <- s |>  mutate(difft=difftime(time(x),time,units="mins")) |>
    filter(abs(difft)<dt,ZW.m<=30) |>
    st_transform(crs(x))
  if (nrow(sonde)>0){
    sonde$rsamp <- extract(x,sonde)[,2]
    #if (any(!is.na(sonde$rsamp))&&any(sonde$rsamp[!is.na(sonde$rsamp)]==0)&m=="tps"&var=="msw.ms") browser()
    #sonde$minDistWithin10kts.m <-  sonde$rsource <- sonde$var <- NA
    #sonde_30ft <- which(sonde$ZW.m<10)
    #if (length(sonde_30ft)>0){
    #  dst<- lapply(sonde_30ft,function(i,y,so) {
    #    #y[y<(&y>(as.numeric(so[i,]$WS.ms)*1.94-10)] <- -9999
    #    target_y <- classify(y, matrix(c(-Inf, as.numeric(so[i,]$WS.ms)*1.94-10, NA, as.numeric(so[i,]$WS.ms)*1.94+10, Inf, NA), ncol=3, byrow=TRUE), others=-9999)
    #    if (all(is.na(values(target_y)))) return(Inf)
    #    d <- distance(target_y)
    #    d <- extract(d, so[i,])[,2]
    #    #d <- distance(y,so[i,],rasterize=TRUE,target=-9999)
    #    #d <- mask(d,y)
    #    #d <- minmax(d)[1]
    #    d
    #  },y=x,so=sonde)
    #  dst <- do.call(c,dst)
    #}else{dst <- NA}
    #sonde$minDistWithin10kts.m[sonde_30ft] <- dst
    sonde$rsource <- m
    sonde$var <- var
    sonde <- st_transform(sonde, 4326)
    coords <- st_coordinates(sonde)
    sonde$LAT <- coords[,2]
    sonde$LON <- coords[,1]
    sonde <- sonde[!is.na(sonde$rsamp),]
    sonde$sourcefile = sources(x)
    sonde
  }
}
