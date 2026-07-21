#' @importFrom sf sf_use_s2 st_intersects
#' @importFrom terra mosaic
#' @importFrom earthdatalogin edl_download
get_srtm <- function(storm,dpath=NULL,coastonly=TRUE){
  if (is.null(dpath)) dpath=tempdir()
  ## start with storm footprint
  if (class(storm)[1]=="sf"){
    ###  need some way to manage missing rocis
    ##  for now, just take the mean
    roci <- storm |> pull(contains("CONS_ROCI"))|>mean(na.rm=TRUE)
    ###  if not present, take the 500km radius
    if (is.na(roci)) roci <- 500000/1852
    extnt <- st_transform(st_union(st_geometry(st_buffer(storm,as.numeric(roci)*1852))),4326)
    aoi = st_bbox(extnt)
  }else if (class(storm)[1]=="SpatRaster"){
    #storm <- app(storm,mean,na.rm=TRUE)
    #storm <- lapply(storm,function(x){
     # x <- trim(storm)
    #  x[x==0] <-NA
    #  nonNA_mask <- !is.na(x)
    #  nonNA_mask[nonNA_mask == 0] <- NA
    #  r_sf <- st_geometry(st_as_sf(terra::as.polygons(nonNA_mask,dissolve=TRUE)))
    #} )
    extnt <- st_transform(storm,4326)
    aoi = st_bbox(extnt)
  }else{
    stop("Unknown storm class in get_srtm function")
  }
 # aoi <- st_transform(st_bbox(trim(Kat)),4326)
  #aoi <- st_transform(aoi,4326)
  ##  get the rounded min max lat and lon
  aoi[1:2] <- floor(aoi[1:2])
  aoi[3:4] <- ceiling(aoi[3:4])
  aoi <- st_as_sfc(aoi)
  ## need to tile the aoi to ensure proper download
  aoi_t <- sf::st_make_grid(aoi, cellsize=1)
  suppressMessages({
    sf_use_s2(FALSE)
    aoi_t <- aoi_t[which(lengths(st_intersects(aoi_t,extnt))>0),]|>st_as_sf()
    aoi_t_land <- aoi_t[which(lengths(st_intersects(aoi_t,rnaturalearth::ne_countries(scale=50)))>0),]|>st_as_sf()
    coast <- rnaturalearth::ne_coastline(scale=50)
    aoi_t_coast <- aoi_t[which(lengths(st_intersects(aoi_t,st_transform(st_buffer(st_transform(coast,st_crs(storm)),.75),4326)))>0),]|>st_as_sf()
    aoi_t_coast <- aoi_t_coast[aoi_t_coast$x %in% aoi_t_land$x,]
    sf_use_s2(TRUE)
  })

  if(coastonly){
    lats <- unlist(lapply(1:nrow(aoi_t_coast),function(x) st_bbox(aoi_t_coast[x,])[2]))#lats <- seq(floor(bbox[2]),ceiling(bbox[4]))
    lons <-unlist(lapply(1:nrow(aoi_t_coast),function(x) st_bbox(aoi_t_coast[x,])[1]))#lons <- seq(floor(bbox[1]),ceiling(bbox[3]))
    tiles <- aoi_t_coast
  }else{
    lats <- unlist(lapply(1:nrow(aoi_t_land),function(x) st_bbox(aoi_t_land[x,])[2]))#lats <- seq(floor(bbox[2]),ceiling(bbox[4]))
    lons <-unlist(lapply(1:nrow(aoi_t_land),function(x) st_bbox(aoi_t_land[x,])[1]))#lons <- seq(floor(bbox[1]),ceiling(bbox[3]))
    tiles <- aoi_t_land
  }
  if(nrow(tiles)==0) {warning("No elevation tiles in storm extent");return(NULL)}
  lons <- ifelse(lons<0,paste0("W",sprintf("%03d",abs(lons))),paste0("E", sprintf("%03d",abs(lons))))
  lats <-  ifelse(lats<0,paste0("S",sprintf("%02d",abs(lats))),paste0("N", sprintf("%02d",abs(lats))))
  latlons <- paste0(lats,lons)
  if (!dir.exists(paste0(dpath,"/SRTM/"))) dir.create(paste0(dpath,"/SRTM"))
  destfiles <- paste0(dpath,"/SRTM/",latlons,".SRTMGL3S.hgt.zip")
  url <- paste0("https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/SRTMGL3S.003/",gsub(".zip","",basename(destfiles)),"/",basename(destfiles))
  elevs <- lapply(1:length(url),function(f,us,dfiles,e,t){
    #if (grepl("N34E140",gsub(".SRTMGL3.|.zip","",dfiles[f])))browser()
    if (file.exists(gsub(".SRTMGL3.|.zip","",dfiles[f]))){
      r=rast(gsub(".SRTMGL3.|.zip","",dfiles[f]))
    }else{
      r <- downr(us[f],dfiles[f])
      if (!is.null(r)) r=rast(r)
    }
    cat("\rDownloading/loading NASA SRTM elevation: %",round(100*f/length(us),1))
    if (!is.null(r)){
      return(list(rast=r,tile=t[f,]))
    }else{
      return(NULL)
    }
  },us=url,dfiles=destfiles,e=bbox,t=tiles)
  cat("\r", paste(rep(" ", 50), collapse=""), "\r")
  elevs <- Filter(Negate(is.null), elevs)
  tiles <- lapply(elevs,function(x) x$tile)
  tiles <- do.call(rbind,tiles)
  elevs <- lapply(elevs,function(x) x$rast)
  #aoi_t_coast <- aoi_t_coast[unlist(elevs),]
  #elevs <- rast(destfiles[unlist(elevs)])
  #if (coastonly) elevs <- list(elevs=elevs,aois=aoi_t_coast)
  #cat("\nMosaicing SRTM elevation\n")
  #elev <- do.call(mosaic, c(elevs, fun = mean))
  return(list(elevs=elevs,tiles=tiles))
}
# downr <- function(u,df,max_attempts =5,attempt=0,success=FALSE) {
#   while(attempt < max_attempts && !success) {
#     attempt <- attempt + 1
#     tryCatch({
#       edl_download(u,dest =df,quiet=FALSE,overwrite=TRUE)
#       unzip(df,exdir=dirname(df))
#       file.remove(df)
#       if (file.exists(file.exists(gsub(".SRTMGL3.|.zip","",df)))){
#         success <- TRUE
#         return(rast(gsub(".SRTMGL3.|.zip","",df)))
#       }
#     }, error = function(e){
#       if (attempt >= max_attempts) {
#         return(NULL)
#       }else{
#         Sys.sleep(0.5)
#       }
#     })
#   }
# }
