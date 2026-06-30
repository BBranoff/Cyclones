#' @importFrom googledrive drive_ls as_id drive_download
get_mswep <- function(cent,dpath=NULL){
  times = sort(unique(cent$date[!is.na(cent$date)]))
  dstart <- first(times[as.numeric(format(times,"%H")) %in% seq(0,23,by=3)])
  ###  for MSWEP, the file names are the starting times
  ###  so we need to subtract 3 hours to capture the first timestep and so on
  dtimes <- seq.POSIXt(dstart,max(times),by=3*60*60)-3*60*60
  id = unique(cent$ID[!is.na(cent$ID)])
  ext <- st_bbox(st_transform(st_buffer(cent,500000),4326))
  ##  get the rounded min max lat and lon
  ext[1:2] <- floor(ext[1:2])
  ext[3:4] <- ceiling(ext[3:4])
  precips <- lapply(dtimes, function(t,e,dp){
    fname1 <- format(t,"%Y%j")
    fname2 <- sprintf("%02d",as.numeric(format(t,"%H")))
    #fname2 <- ceiling(as.numeric(format(t,"%H"))/3)*3
    fname=paste(fname1,fname2,"nc",sep=".")
    ###  for past files
    drivefiles_past <- drive_ls(as_id("1DVR90Ud1C444bTOPeENX-3I7tgqPLnao"), q = paste0("name = '",fname,"'"),recursive=FALSE)
    ##  for more recent files
    drivefiles_recent <- drive_ls(as_id("1XBonFS0_t3aSM_C4CYobwsN-spmj29vo"), q = paste0("name = '",fname,"'"),recursive=FALSE)
    drivefiles <- list(drivefiles_past,drivefiles_recent)[which(c(nrow(drivefiles_past),nrow(drivefiles_recent))>0)][[1]]
    ##  terra will return a raster with missing crs if filename is not correctly formatted
    fp <-  normalizePath(paste0(dp,fname))
    if (file.exists(fp)){
      r=rast(fp)$precipitation
    }else if(nrow(drivefiles)==1){
      downr <- function(df,max_attempts =5,attempt=0,success=FALSE) {
        while(attempt < max_attempts && !success) {
          attempt <- attempt + 1
          tryCatch({
            drive_download(df$id[1],path=fp,overwrite=TRUE)
            success <- TRUE
            return(rast(fp))
          }, error = function(e){
            file.remove(fp)
            if (attempt >= max_attempts) {
              return(NULL)
            }else{
              Sys.sleep(2)
            }
          })
        }
      }
      r <- downr(drivefiles)
      if (!is.null(r))  r <- r$precipitation
      #   ###  some files have missing extents and layer names
      #   if (crs(r)!=""){ r <- r$precipitation
      #   ###  if thats the case, get the extent of a 'similar' file
      #   }else{
      #     fname <- gsub(format(t,"%Y%j"),paste0(as.numeric(format(t,"%Y"))+1,format(t,"%j")),fname)
      #     drivefiles_past <- drive_ls(as_id("1DVR90Ud1C444bTOPeENX-3I7tgqPLnao"), q = paste0("name = '",fname,"'"),recursive=FALSE)
      #     drivefiles_recent <- drive_ls(as_id("1XBonFS0_t3aSM_C4CYobwsN-spmj29vo"), q = paste0("name = '",fname,"'"),recursive=FALSE)
      #     drivefiles <- list(drivefiles_past,drivefiles_recent)[which(c(nrow(drivefiles_past),nrow(drivefiles_recent))>0)][[1]]
      #     if (file.exists(paste0(dp,fname)))  rtmp=rast(paste0(dp,"/",fname))$precipitation
      #     else rtmp <- downr(drivefiles)
      #     if (crs(rtmp)!=""){
      #       crs(r) <- crs(rtmp)
      #       ext(r) <- ext(rtmp)
      #       names(r) <- "precipitation"
      #       writeRaster(r,paste0(dp,fname))
      #     }
      #   }
      # }
    }else{
      times <- times[-which(times==t)]
      return(NULL)
    }
    cat("\rDownloading/loading MSWEP precipitation: %",round(100*which(dtimes==t)/length(dtimes),1))
    if (!is.null(r)){
      ###  for MSWEP, the file names are the starting times, but we want the end times
      ###  so add in the 3 hours
      terra::time(r) <- rep(as.POSIXct(gsub(".nc","",fname),format="%Y%j.%H",tz="UTC")+3*60*60,nlyr(r))
      if (crs(r)=="") terra::crs(r) <- "epsg:4326";terra::ext(r) <- c(-180,180,-90,90)
      crop(r,e)
    }else{
      NULL
    }
  },e=ext,dp=dpath)
  if (all(is.null(precips))) stop("Download or loading MSWEP precipitation failed.")
  cat("\n")
  wrap(rast(precips))
}
