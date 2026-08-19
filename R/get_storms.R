#' Gather storm(s) tabular data.
#'
#'`get_storms()` will accept an already loaded data frame, a filename, or a target web source to download cyclone tabular data and
#' format it for processing in the Cyclones functions.
#'
#' @param source string ('ncei' or 'hurdat'), filename, or pre-loaded data frame. If string, should be a target web source (either 'ncei' or 'hurdat') for downloading,
#' or a filename of a .csv or .nc file containing tabular cyclone data. Can also be pre-loaded data.frame. In the case of an existing data.frame, it should be consistent with
#' IBTrACS or HURDAT format.
#' @param id Either the SID (from IBTrACS) or the USA_ATCF_ID (from both hurdat and IBTrACS) used to identify individual storms. Results will be filtered to these values.
#' @param name The official name of the storm(s). Results will be filtered to these values.
#' @param season The year associated with the date of the storm(s). Results will be filtered to these values.
#' @param basin The oceanic basin of the storm(s). Results will be filtered to these values. In case a storm traverses multiple basins, its origin basin is used.
#' @param ib_filt A filtering value for downloading specific IBTrACS data and to avoid downloading the entire dataset if not required. Valid values are:
#' 'ACTIVE', 'ALL', 'last3years', 'since1980', or the basins: "EP", "NA", "NI", "SA", "SI", "SP", "WP".
#' @param consolidate To run 'cons_stormdat()' or not, to consolidate values across agencies. When TRUE, 'CONS_' columns will be added to represent the consolidated values and, unless
#' the 'cols' parameter requests otherwise, the original columns will be dropped. If consolidation is not run, downstream functions may fail or return unexpected results,
#' as they may not be able to distinguish which columns to use. Also, for the maximum sustained wind columns, the interval period is inconsistent among agencies. Consolidation
#' translates these values to a common interval. See 'cons_stormdat()'.
#' @param cols Which columns to retain in the output if consolidate is TRUE. Default is that only consolidated values are retained, in addition to storm identification, location,
#' and time information. Input values can be specific column names or agency prefixes, in which case all column from that agency are retained. If consolidation is FALSE,
#' all original columns are retained when 'cols' is empty.
#' @param erddap Exploratory, waiting on ERDDAP host fixes. Whether to retrieve filtered IBTRACS data from ERDDAP, which negates the need to download a full subset if only
#' one storm is needed.
#' @param returndf if the results shuold be returned as a single data frame (TRUE) or as a list of single storm data frames (FALSE, default).
#' @param ... Additional arguments passed to cons_stormdat(), e.g. 'pref', 'msw_int'.
#' @returns A list of data.frames, one for each storm, containing tabular information required for Cyclone processing. The names for each data.frame in the list are
#' unique identifiers for storms and used throughout Cyclones as filenames for various functions.
#' @export
#' @examples
#' # default will get the last 3 years of storms from IBTrACS (ncei)
#' storms <- get_storms()
#'
#' # get a specific storm
#' # older years may require extending the timeout for download
#' options(timeout=300)
#' Maria_2017 <- get_storms(name="MARIA",season=2017)
#' Maria_2017 <- get_storms(id="2017260N12310")
#'
#' # get all storms for a basin
#' NorAtl <- get_storms(basin="NA")
#'
#' @importFrom dplyr filter mutate any_of select pull slice group_split syms coalesce group_keys first if_else group_by
get_storms <- function(source="ncei",id=NULL,name=NULL,season=NULL,basin=NULL,ib_filt=NULL,consolidate=TRUE,cols=NULL,erddap=FALSE,returndf=FALSE,...){
  tmf <- tempfile(pattern=paste0(c(source,id,name,season,basin,ib_filt),collapse="_"))
  on.exit(unlink(tmf), add = TRUE)
  if (is.data.frame(source)){
    dat=source
  }else if (is.character(source)){
    if(source=="ncei") {
      cat("Downloading IBTrACS from: https://www.ncei.noaa.gov/products/international-best-track-archive\n")
      if (erddap){#&!any(is.null(c(id,name,season,basin)))) {
        dat = get_ibERDDAP(id,name,season,basin)
      }else{dat=NULL}
      if (is.null(dat)){
        if (is.null(season)&is.null(ib_filt)) season = (as.numeric(format(Sys.time(),"%Y")))
        if (is.null(ib_filt)&&isTRUE(as.numeric(season)>=(as.numeric(format(Sys.time(),"%Y"))-3))) ib_filt = "last3years" else if (is.null(ib_filt)&&isTRUE(as.numeric(season)>1980)) ib_filt= "since1980"
        ## get appropriate url
        if (!is.null(ib_filt)){
          if (!ib_filt %in% c("ACTIVE","ALL","EP","NA","NI","SA","SI","SP","WP","last3years","since1980")) stop("ib_filt value not found in https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/") else
            URL <- paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.",ib_filt,".list.v04r01.csv")
        }else{
          if (!is.null(basin)){
            if (!basin %in% c("EP","NA","NI","SA","SI","SP","WP")) stop("Provided basin value not found in IBTrACS") else
              URL =  paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.",basin,".list.v04r01.csv")
          }else{
            URL <- paste0("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.ALL.list.v04r01.csv")
          }
        }
        ###  download with url
        tryCatch({
          # Attempt to download the file
          download.file(URL,paste0(tmf,".csv"))
          dat = read.csv(paste0(tmf,".csv"))
          message("Download successful!")
        },
        error = function(e) {
          # Check if the error message indicates a timeout
          if (grepl("Timeout", e$message, ignore.case = TRUE)) {
            warning("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
          } else {
            # Handle other potential errors (e.g., file not found, permission issues)
            warning(paste("Download failed with error:", e$message))
          }
          # Return NULL or an indicator of failure
          invisible(NULL)
        },
        warning = function(w) {
          if (grepl("Timeout|!= reported length", w$message, ignore.case = TRUE)) {
            stop("Download timeout reached. Please check your internet connection or increase the R timeout option. Alternatively, download and source data directly from browser.")
          }else{
            # Handle other warnings, such as "downloaded length != reported length"
            message(paste("A warning occurred during download:", w$message))
          }
        })
      }
      ##  HURDAT
    }else if (source=="hurdat"){
      #consolidate=FALSE
      cat("HURDAT data missing ROCI column. This value will be modelled/estimated when needed to provide maximum storm extent.\n")
      if (is.null(basin)){
        url <- c("https://www.aoml.noaa.gov/hrd/hurdat/hurdat2.html","https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac.html")
        basins <- c("NA","EP")
      }else if (basin=="NA"){
        cat("Downloading North Atlantic HURDAT2 from: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html\n")
        url = "https://www.aoml.noaa.gov/hrd/hurdat/hurdat2.html"
      }else if(basin=="EP"){
        cat("Downloading Northeast Pacific HURDAT2 from: https://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html\n")
        url = "https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac.html"
      }else{
        stop("Basin not available in HURDAT. HURDAT available for North Atlantic 'NA' and 'EP' basins.")
      }
      if (length(url)==2){
        for (u in 1:2){
          dt <- get_hurdat(url[u],tmf,id,basins[u])
          if (u!=1) dat <- rbind(dat,dt)
          else dat <- dt
        }
      }else{
        dat <- get_hurdat(url,tmf,id,basin)
      }
    }else if(grepl(".nc",source)){
      if (!requireNamespace("ncdf4", quietly = TRUE)) {
        stop("Package \"ncdf4\" must be installed to use this function (and some others). \n Use `install.packages(\"ncdf4\")` to install it.", call. = FALSE)
      }
     dat <- get_nc(source)
    } else {
      if (!grepl(".csv",source)) stop("file must be sourced from 'ncei', 'hurdat' or from a .nc or .csv file")
      dat=read.csv(source)
    }
  }else{
    stop("unknown data source")
  }
  if (any(grepl("nmile",dat|>slice(1)))) dat <- dat |> slice(-1)
  dat <- dat |>
    ## get the best available unique ID, SID is preferred, but USA_ATCFID if not
    mutate(ID=coalesce(!!!syms(intersect(c("SID","USA_ATCF_ID"),names(dat)))),
           ##  longitudes are inconsistent for storms that cross the dateline, make consistent
           LON=if_else(as.numeric(LON)>180,as.numeric(LON)-360,as.numeric(LON))) |>
    ## some storms traverse basins so must take only the original basin in the name
    group_by(ID) |>
    ## cant have colons in name, mucks up filenames in downstream processes
    mutate(ID=paste(gsub(":","-",NAME),SEASON,first(BASIN),ID, sep="_"))|>ungroup()
  ###  for now, consolidation only works for ibtracs data
  if (consolidate==TRUE&&(source=="ncei"|any(grepl("_WIND",names(dat))))){
    dat <- cons_stormdat(dat)
    dat=dat |>
      select(any_of(c("ID","SID","USA_ATCF_ID","SEASON","NAME","BASIN","ISO_TIME","LAT","LON","USA_SSHS","STORM_SPEED")),
             contains(c("CONS_",cols)))
  }
  dat <- dat |>
    mutate(ISO_TIME=as.POSIXct(ISO_TIME,tz="UTC"),
           BASIN=if_else(is.na(BASIN),"NA",BASIN))
  if (!is.null(basin)){
    if (toupper(basin) %in% c("EP","NA","NI","SA","SI","SP","WP"))
      dat <- dat |> filter(BASIN %in% toupper(basin))
    else
      stop("Provided basin value not found in data. Check source and other filters")
  }
  if (!is.null(season)){
    if (as.numeric(season) %in% seq(1851,format(Sys.time(),"%Y")))
      dat <- dat |> filter(as.numeric(SEASON) %in% as.numeric(season))
    else
      stop("Provided season value not found in data. Check source and other filters.")
  }
  if (!is.null(name)){
    names=dat |> pull(NAME)|>unique()
    if (sum(grepl(paste0(toupper(name),collapse="|"),names),na.rm=TRUE)>0){
      if (sum(grepl(paste0(toupper(name),collapse="|"),names),na.rm=TRUE)<length(name)){warning("Not all names found in data. Check other filters.")}
      dat <- dat |> filter(NAME %in% toupper(name))
    }else{ stop("Provided season value not found in data. Check source and other filters.")}
  }
  if (source=="hurdat"){
    ##  add in the STORMSPEED column after filtering to save time
    sfdat <- st_as_sf(dat,coords=c("LON","LAT"),crs=4326)|>
      group_by(USA_ATCF_ID)|>
      mutate(distance=sf::st_distance(geometry, lead(geometry), by_element = TRUE),
             timestep_hrs=difftime(lead(ISO_TIME),ISO_TIME,units="hours"),
             STORMSPEED=round(as.numeric(distance)/1852/as.numeric(timestep_hrs),1))
    dat$STORM_SPEED <- sfdat$STORMSPEED
  }
  if (!is.null(id)){
    if ((id %in% unique(dat$SID)|id %in% unique(dat$USA_ATCF_ID)))
      dat <- dat |> filter(SID %in% id|USA_ATCF_ID %in% id)
    else
      stop("Provided id value not found in data. Check source and other filters. HURDAT uses the USA_ATCFID while IBTrACS uses either USA_ATCFID or SID.")
  }
  if (!returndf){
    dat <- dat |> group_by(ID)
    IDs <- group_keys(dat)|>pull(ID)
    dat <- dat |>group_split()
    names(dat) <- IDs
  }
  return(dat)
}
get_hurdat <- function(u,tf,ID,basin){
  tryCatch({
    # Attempt to download the file
    if (length(list.files(dirname(tf),pattern=paste((strsplit(basename(tf),"_"))[[1]][1:3],collapse="_")))>0){
      cat(paste0("Loading previously downloaded file from temp folder: ",tempdir()))
      tf <-gsub(".txt","",list.files(dirname(tf),pattern=paste((strsplit(basename(tf),"_"))[[1]][1:3],collapse="_"),full.names = TRUE))
    }
    else (download.file(u,paste0(tf,".txt")))
    hurdat <- readLines(paste0(tf,".txt"))
    message("Download successful!")
  },
  error = function(e) {
    # Check if the error message indicates a timeout
    if (grepl("Timeout", e$message, ignore.case = TRUE)) {
      stop("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
    } else {
      # Handle other potential errors (e.g., file not found, permission issues)
      stop(paste("Download failed with error:", e$message))
    }
    # Return NULL or an indicator of failure
    invisible(NULL)
  },
  warning = function(w) {
    if (grepl("Timeout", w$message, ignore.case = TRUE)) {
      stop("Download timeout reached. Please check your internet connection or increase the R timeout option. Alernatively, download and source data directly from browser.")
    }else{
      # Handle other warnings, such as "downloaded length != reported length"
      message(paste("A warning occurred during download:", w$message))
    }
  })
  start <- grep("body",head(hurdat)) +2
  hurdat <- hurdat[start:length(hurdat)]
  ls <- lengths(strsplit(hurdat,","))
  hdlines <- which(ls==ls[1])
  idlines <- hurdat[hdlines]
  stormids <- sapply(strsplit(idlines,","),"[[",1)
  stormlengths <- as.numeric(sapply(strsplit(idlines,","),"[[",3))
  RLES <- if (is.null(ID)) cbind(hdlines, stormlengths) else cbind(hdlines[grep(paste(ID,collapse="|"),stormids)],stormlengths[grep(paste(ID,collapse="|"),stormids)])
  RLES <- split(RLES,seq(nrow(RLES)))
  storms <- lapply(RLES,function(x,y){
    read.table(text=y[(x[1]+1):(x[1]+x[2])],sep=",",header=FALSE,
               col.names=c("date","time","recID","Status","LAT","LON","USA_MSW","USA_PRES",
                           "USA_R34_NE", "USA_R34_SE","USA_R34_SW","USA_R34_NW",
                           "USA_R50_NE", "USA_R50_SE","USA_R50_SW","USA_R50_NW",
                           "USA_R64_NE", "USA_R64_SE","USA_R64_SW","USA_R64_NW","USA_RMW")) |>
      mutate(USA_ATCF_ID=strsplit(y[x[1]],",")[[1]][1],
             NAME=gsub("^\\s+|\\s+$","",strsplit(y[x[1]],",")[[1]][2]),
             ISO_TIME=as.POSIXct(paste0(substr(date,1,4),"-",substr(date,5,6),"-",substr(date,7,8)," ",sprintf("%04d",time)),format="%Y-%m-%d %H",tz="UTC"),
             LAT = as.numeric(if_else(grepl("N",LAT),gsub("N| ","",LAT),
                                      if_else(grepl("S",LAT),paste0("-",gsub("S| ","",LAT)),LAT))),
             LON = as.numeric(if_else(grepl("E",LON),gsub("E| ","",LON),
                                      if_else(grepl("W",LON),paste0("-",gsub("W| ","",LON)),LON))),
             USA_SSHS = if_else(is.na(USA_MSW),NA,
                                if_else(USA_MSW>=135,5,if_else(USA_MSW>=114,4,if_else(USA_MSW>=96,3,if_else(USA_MSW>=84,2,if_else(USA_MSW>=65,1,if_else(USA_MSW>=34,0,-1)))))))
             )|>
      rename(USA_WIND=USA_MSW)
  },y=hurdat)
  dat <- do.call(rbind,storms) |>
    mutate(BASIN=basin,SEASON=substr(date,1,4))

  dat
}
get_nc <-  function(srce){
  dat_nc <- ncdf4::nc_open(srce)
  #vars <- tolower(c("SID","USA_ATCF_ID","SEASON","NAME","BASIN","ISO_TIME","LAT","LON","USA_WIND","USA_PRES","USA_SSHS",
  #                  "USA_R34_NE", "USA_R34_SE","USA_R34_SW","USA_R34_NW",
  #                  "USA_R50_NE", "USA_R50_SE","USA_R50_SW","USA_R50_NW",
  #                  "USA_R64_NE", "USA_R64_SE","USA_R64_SW","USA_R64_NW",
  #                  "USA_ROCI","USA_POCI","USA_RMW","REUNION_RMW","BOM_RMW","USA_EYE","STORM_SPEED"))
  vars <- lapply(dat_nc$var,function(x) ncdf4::ncvar_get(dat_nc,x$name))
  varlns <- lengths(vars)
  valvars <- vars[varlns==varlns[names(varlns)=="iso_time"]]
  quadvars <- vars[varlns==varlns[names(varlns)=="usa_r34"]]
  idvars <- vars[varlns==varlns[names(varlns)=="sid"]]
  nstorms <-varlns[names(varlns)=="sid"]
  dat <- lapply(1:nstorms, function(x){
    ids <- lapply(idvars,function(i){id = rep(i[x],idvars$numobs[x])})
    ids<- do.call(cbind,ids)
    vals <- lapply(valvars,function(v){ val = v[1:idvars$numobs[x],x]})
    vals <- do.call(cbind,vals)
    quads <- lapply(seq_along(quadvars), function(q){
      quad=quadvars[[q]][,1:idvars$numobs[x],x]
      if (idvars$numobs[x]==1) quad <- data.frame(cbind(quad[1],quad[2],quad[3],quad[4]))
      else quad <- data.frame(cbind(quad[1,],quad[2,],quad[3,],quad[4,]))
      names(quad)<-paste(names(quadvars[q]),c("NE","SE","SW","NW"),sep="_")
      quad
      })
    quads <- do.call(cbind,quads)
    cbind(ids,vals,quads)
    })
  dat <- do.call(rbind,dat)
  names(dat) <- toupper(names(dat))
  dat
}
get_ibERDDAP <- function(id,name,season,basin){
  baseurl <- "https://erddap.aoml.noaa.gov/hdb/erddap/tabledap/IBTRACS"
  baseurl <- paste0("https://erddap.aoml.noaa.gov/hdb/erddap/tabledap/IBTRACS_","last3years",".csv")
  colurl = "?numobs%2Csid%2Cseason%2Cnumber%2Cbasin%2Csubbasin%2Cname%2Csource_usa%2Csource_jma%2Csource_cma%2Csource_hko%2Csource_new%2Csource_reu%2Csource_bom%2Csource_nad%2Csource_wel%2Csource_td5%2Csource_td6%2Csource_ds8%2Csource_neu%2Csource_mlc%2Ctime%2Ciso_time%2Cnature%2Clatitude%2Clongitude%2Cwmo_wind%2Cwmo_pres%2Cwmo_agency%2Ctrack_type%2Cmain_track_sid%2Cdist2land%2Clandfall%2Ciflag%2Cusa_agency%2Cusa_atcf_id%2Cusa_lat%2Cusa_lon%2Cusa_record%2Cusa_status%2Cusa_wind%2Cusa_pres%2Cusa_sshs%2Cusa_r34%2Cusa_r50%2Cusa_r64%2Cusa_poci%2Cusa_roci%2Cusa_rmw%2Cusa_eye%2Ctokyo_lat%2Ctokyo_lon%2Ctokyo_grade%2Ctokyo_wind%2Ctokyo_pres%2Ctokyo_r50_dir%2Ctokyo_r50_long%2Ctokyo_r50_short%2Ctokyo_r30_dir%2Ctokyo_r30_long%2Ctokyo_r30_short%2Ctokyo_land%2Ccma_lat%2Ccma_lon%2Ccma_cat%2Ccma_wind%2Ccma_pres%2Chko_lat%2Chko_lon%2Chko_cat%2Chko_wind%2Chko_pres%2Cnewdelhi_lat%2Cnewdelhi_lon%2Cnewdelhi_grade%2Cnewdelhi_wind%2Cnewdelhi_pres%2Cnewdelhi_ci%2Cnewdelhi_dp%2Cnewdelhi_poci%2Creunion_lat%2Creunion_lon%2Creunion_type%2Creunion_wind%2Creunion_pres%2Creunion_tnum%2Creunion_ci%2Creunion_rmw%2Creunion_r34%2Creunion_r50%2Creunion_r64%2Cbom_lat%2Cbom_lon%2Cbom_type%2Cbom_wind%2Cbom_pres%2Cbom_tnum%2Cbom_ci%2Cbom_rmw%2Cbom_r34%2Cbom_r50%2Cbom_r64%2Cbom_roci%2Cbom_poci%2Cbom_eye%2Cbom_pos_method%2Cbom_pres_method%2Cnadi_lat%2Cnadi_lon%2Cnadi_cat%2Cnadi_wind%2Cnadi_pres%2Cwellington_lat%2Cwellington_lon%2Cwellington_wind%2Cwellington_pres%2Cds824_lat%2Cds824_lon%2Cds824_stage%2Cds824_wind%2Cds824_pres%2Ctd9636_lat%2Ctd9636_lon%2Ctd9636_stage%2Ctd9636_wind%2Ctd9636_pres%2Ctd9635_lat%2Ctd9635_lon%2Ctd9635_wind%2Ctd9635_pres%2Ctd9635_roci%2Cneumann_lat%2Cneumann_lon%2Cneumann_class%2Cneumann_wind%2Cneumann_pres%2Cmlc_lat%2Cmlc_lon%2Cmlc_class%2Cmlc_wind%2Cmlc_pres%2Cusa_gust%2Cbom_gust%2Cbom_gust_per%2Creunion_gust%2Creunion_gust_per%2Cusa_seahgt%2Cusa_searad%2Cstorm_speed%2Cstorm_dir"
  url <- paste0(baseurl,colurl)
  if (!is.null(season)&&(season>=2022|season<=2025)){
    #yearfilt=paste0("&season=~%22",paste0(toupper(season),collapse="|"),"%22")
    yearfilt=paste0("&season%3E=",min(as.numeric(season)),"&season%3C=",max(as.numeric(season)))
    #timefilt=paste0("&time%3E=",min(season),"-01-01T00%3A00%3A00Z&time%3C=",max(season),"-12-31T00%3A00%3A00Z")
    url <- paste0(url,yearfilt)
  }else{yearfilt=NULL}
  if (!is.null(name)){
    namefilt=paste0("&name=~%22",paste0(toupper(name),collapse="|"),"%22")
    url <- paste0(url,namefilt)
  }else{namefilt=NULL}
  if (!is.null(basin)){
    basinfilt=paste0("&basin=~%22",paste0(toupper(basin),collapse="|"),"%22")
    url <- paste0(url,basinfilt)
  }else{basinfilt=NULL}
  dat <- tryCatch({
    read.csv(paste0(baseurl,colurl,namefilt)) |> slice(-1)
    },
    error = function(x){
      warning(paste("No Storms Found at: ", baseurl))
      return(NULL)}
    )
  names(dat) <- toupper(names(dat))
  dat <- dat |>
    rename(LON=LONGITUDE,LAT=LATITUDE)
  return(dat)
}
