get_ecmwf <- function(cent,dpath=tempdir(),dfiles=NULL,GC=FALSE){
  ###  times in ecmwfr are for the previous hour, so no need to change for consistency with the other sources (which are for the next hour or three hours)
  times = sort(seq(min(cent$date,na.rm=TRUE),max(cent$date,na.rm=TRUE),by=unique(as.numeric(cent$dt)[!is.na(as.numeric(cent$dt))])*60*60))
  id = unique(cent$ID[!is.na(cent$ID)])
  ext <- st_bbox(st_shift_longitude(st_transform(st_buffer(cent,500000),4326)))
  ##  get the rounded min max lat and lon
  ext[1:2] <- floor(ext[1:2])
  ext[3:4] <- ceiling(ext[3:4])
  ###  get the blocks of full days. These can be requested all at once
  days = unique(format(times, "%j"))
  requests <- lapply(days,function(d,t,e,id){
    ts <- t[format(t,"%j")==d]
    request <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type = "reanalysis",
      variable = "total_precipitation",#,"10m_v_component_of_wind"),
      year = unique(format(ts,"%Y")),
      month =  unique(format(ts,"%m")),
      day =  unique(format(ts,"%d")),
      time = format(seq(min(ts),max(ts),"hours"),"%H:%M"),
      #time = c("00:00", "01:00", "02:00", "03:00","04:00","05:00",
       #        "06:00","07:00", "08:00", "09:00", "10:00","11:00",
       #        "12:00","13:00", "14:00", "15:00", "16:00","17:00",
       #        "18:00","19:00","20:00", "21:00", "22:00", "23:00"),
      data_format = "grib",
      download_format = "unarchived",
      area = c(e$ymax,e$xmin,e$ymin,e$xmax),
      target = paste0("precip_ecmwf_",id,"_",unique(format(ts,"%Y-%m%-%d")))
    )},t=times,e=ext,id=id)
  ###  look in the temp directory if already there
  ###  this saves time and download resources when the same file can be used to recalculate a new time step frequency for example
  if (is.null(dfiles)){
    ##  if dfiles is null, it was called outside of parallel, so no need for recursive
    dfiles <- list.files(dpath,recursive=TRUE,pattern="ecmwf",full.names=TRUE)
    dfiles <- dfiles[grep(id,dfiles)]
  }
  if (length(dfiles)>0&&any(lengths(strsplit(basename(dfiles),"_"))>=7))  dfiles <- dfiles[which(lengths(strsplit(basename(dfiles),"_"))>=7)];dfiles <- dfiles[format(as.Date(gsub(".grib|.tif","",sapply(strsplit(basename(dfiles),"_"),"[[",7))),"%j") %in% format(cent$date,"%j")]
  ###  cant mix ecmwf and ecmwf_GC, extents are slightly different
  ##  so remove one if the other is present
  if (GC&any(grepl(".grib",dfiles))){ file.remove(dfiles[grep(".grib",dfiles)]);dfiles <- character(0)}
  if (!GC&any(grepl(".tif",dfiles))){ file.remove(dfiles[grep(".tif",dfiles)]);dfiles <- character(0)}

  ###  if all of the files are already downloaded, just use those
  if (all(paste0("precip_ecmwf_",id,"_",unique(format(times,"%Y-%m%-%d")),".grib") %in% basename(dfiles))|
      all(paste0("precip_ecmwf_",id,"_",unique(format(times,"%Y-%m%-%d")),".tif") %in% basename(dfiles))){
    requests <- list(done=dfiles)
    ##  if only some of them, download the missing ones
  }else if(any(paste0("precip_ecmwf_",id,"_",unique(format(times,"%Y-%m%-%d")),rep(c(".grib",".tif"),each=length(unique(format(times,"%Y-%m%-%d"))))) %in% basename(dfiles))){
    requests <- list(done=dfiles,request=requests[-grep(paste0(gsub(".tif|.grib","",basename(dfiles)),collapse="|"),sapply(requests,"[[",11))])
  }else{
    requests <- list(request=requests)
  }##  or download them all
  if ("request" %in% names(requests)){
    getreqs <-requests$request
    wrkrs = length(getreqs)
    max_attempts =5;attempt=0
    if (!GC){
      while(attempt < max_attempts ) {
        attempt <- attempt + 1
        get_reqs <- tryCatch({
          ecmwfr::wf_request_batch(getreqs,path=dpath,workers=wrkrs)
          paste0(dpath,lapply(getreqs,function(x) x$target),".grib")
        }, error = function(e){
          if (attempt >= max_attempts) {
            stop(e)
          }else{
            Sys.sleep(0.5)
            return(NULL)
          }
        })
        if (!is.null(get_reqs)) break
      }
    } else{
      get_reqs <- get_ecmwf_google(getreqs,dpath)
    }
    if ("done" %in% names(requests)){
      reqs <- c(get_reqs,requests$done)
    }else{
      reqs <- get_reqs
    }
  }else{
    reqs <- requests$done
  }
  cat("\n")
  wrap(rast(reqs)*1000)
}
get_ecmwf_google <- function(req,df){
  lapply(req, function(x){
    url <- paste0("https://storage.googleapis.com/gcp-public-data-arco-era5/raw/date-variable-single_level/",x$year,"/",x$month,"/",x$day,
                  "/total_precipitation/surface.nc")
    tf <- tempfile(fileext=".nc")
    attempt = 0
    while(attempt < 5 ) {
      attempt <- attempt + 1
      tryCatch({
        download.file(url,tf,mode="wb")
      }, error = function(e){
        if (attempt >= 5) {
          if (file.exists(tf)) file.remove(tf)
          stop(e)
        }else{
          if (file.exists(tf)) file.remove(tf)
          Sys.sleep(0.5)
          return(NULL)
        }
      })
      if (file.exists(tf)) break
    }
    r <- rast(tf)
    r <- crop(r,ext(c(x$area[2],x$area[4],x$area[3],x$area[1])))#st_shift_longitude(st_as_sf(st_as_sfc(st_bbox(c(x$area[2],x$area[4],x$area[3],x$area[1]))),crs=4326)))
    r <- r[[format(time(r),"%H:%M") %in% x$time]]
    ### in at least one GC layer, some hours are duplicated and empty
    r <- r[[!duplicated(time(r))]]
    targ <- paste0(df,x$target,".tif")
    writeRaster(r,targ)
    file.remove(tf)
    targ
  })|>unlist()
}

