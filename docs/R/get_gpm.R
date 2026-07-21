#' @importFrom earthdatalogin edl_download
#' @importFrom terra crs ext
get_gpm <- function(cent,dpath){
  times = sort(unique(cent$date[!is.na(cent$date)]))
  ###  IMERG data re half hourly
  dstart <- first(times[as.numeric(format(times,"%H")) %in% seq(1,23,by=1)])
  ###  the file names are centered on the start times, but we want to capture the previous 30 minutes, not the next 30 minutes
  ###  so subtract 30 minutes from the times so they correspond to the times we want
  dtimes <- seq.POSIXt(dstart,max(times),by=30*60)-30*60
  id = unique(cent$ID[!is.na(cent$ID)])
  ext <- st_bbox(st_transform(st_buffer(cent,500000),4326))
  ##  get the rounded min max lat and lon
  ext[1:2] <- floor(ext[1:2])
  ext[3:4] <- ceiling(ext[3:4])
  jday <- format(dtimes,"%j")
  day <- format(dtimes,"%d")
  year <- format(dtimes,"%Y")
  month <- format(dtimes,"%m")
  HH <- format(dtimes,"%H")
  MM <- format(dtimes,"%M")
  tstep <- difftime(dtimes,as.Date(dtimes),units="mins")
  tstep <- sprintf("%04d",tstep)
  ## need to check for different early/late products
  ## default to late product
  files <- list(late=list( url=paste0("https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHH.07/",year,"/",jday,"/3B-HHR.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.HDF5"),
                          destfiles = paste0(dpath,"/3B-HHR.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.h5")),
                early_late=list(url=paste0("https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHL.07/",year,"/",jday,"/3B-HHR-L.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.HDF5"),
                                destfiles = paste0(dpath,"/3B-HHR-L.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.h5")),
                early=list(url=paste0("https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHE.07/",year,"/",jday,"/3B-HHR-E.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.HDF5"),
                           destfiles = paste0(dpath,"/3B-HHR-E.MS.MRG.3IMERG.",year,month,day,"-S",HH,MM,"00-E",HH,as.numeric(MM)+29,"59.",tstep,".V07B.h5")))
  precips <- lapply(1:length(files$late$url),function(f,fles,e,t){
    cat("\rDownloading/loading NASA GPM precipitation: %",round(100*f/length(fles$late$url),1))
    for (FILE in fles){
      if (file.exists(FILE$destfiles[f])){
        r=rast(FILE$destfiles[f])
        break
      } else if (file.exists(gsub("V07B","V07C",FILE$destfiles[f]))){
        r=rast(gsub("V07B","V07C",FILE$destfiles[f]))
        break
      }else{
        r <- NULL
        if (rvest::session(FILE$url[f])$response$status_code==404){
          ###  if there is not an unauthrized access code, it doesnt exist, try the other versions
          next
        }else {
          ##  otherwise try to download it
          r <- downr(FILE$url[f],FILE$destfiles[f])
          ###  sometimes the V07B is not available, so try V07C
          if (is.null(r)){
            r <- downr(u=gsub("V07B","V07C",FILE$url[f]),
                       d=gsub("V07B","V07C",FILE$destfiles[f]))
          }
          if (!is.null(r)){
            break
          }
        }
      }
    }
    if (is.null(r)){
      return(NULL)
    }else{
      ####  the netcdf data are reversed and transposed relative to terra's read format
      r <- terra::t(terra::rev(r$precipitation))
      terra::crs(r)="EPSG:4326"
      terra::ext(r) =c(-180, 180, -90, 90)
      ###  time should reflect the previous 30 minutes, add back in the 30 minutes that we subtracted to get the naming convention
      terra::time(r) <-  t[f]+30*60
      return(crop(r,e))
    }
  },fles=files,e=ext,t=dtimes)
  precips <- Filter(Negate(is.null),precips)
  cat("\n")
  wrap(rast(precips))
}
downr <- function(u,d,max_attempts =5,attempt=0,success=FALSE) {
  while(attempt < max_attempts && !success) {
    attempt <- attempt + 1
    tryCatch({
      ####  edl_download will get an empty file that does not exist
      edl_download(u,d,quiet=FALSE,overwrite=TRUE)
      return(rast(d))
    }, error = function(e){
      file.remove(d)
      if (attempt >= max_attempts) {
        return(NULL)
      }else{
        Sys.sleep(2)
      }
    })
  }
}
