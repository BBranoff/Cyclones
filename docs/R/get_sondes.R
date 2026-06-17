get_sondes <-  function(storm,ddir=NULL,checksondes=TRUE){
  storm <- tolower(storm)
  ###  read the parsed data from the temp dir if created this session
  tmdir= tempdir()
  donefile  = paste0(tmdir,"/",storm,".rds")
  if (file.exists(donefile)) {return(readRDS(donefile))}else {created=TRUE}
  ##  otherwise, attempt to download or read the local files if already downloaded
  if (!is.null(ddir)){
    if (!dir.exists(ddir)) stop("\nSonde download directory not found")
    fls <- list.files(ddir)
    flstorms <- sapply(strsplit(fls,"_"),"[[",1)
    fls <- fls[flstorms==storm]
  } else{cat("\nNo directory supplied, downloading from www.aoml.noaa.gov/hrd/Storm_pages/");fls=NULL}
  tryURL <- function(storm){
    address=paste0("https://www.aoml.noaa.gov/hrd/Storm_pages/",storm,"/sonde.html")
    # Try to open the URL
    con <- try(suppressWarnings(url(address)), silent = TRUE)
    # If open() failed, or if the URL can't be read, return FALSE
    if (inherits(con, "try-error")) {
      return(FALSE)
    }
    # Attempt to read the first line to confirm it's actually readable
    res <- try(suppressWarnings(readLines(con, n = 1)), silent = TRUE)
    close(con)
    # If readLines failed (e.g. 404), return FALSE
    if (inherits(res, "try-error")) {
      return(FALSE)
    }
    # If we made it here, the page exists
    page = rvest::read_html(paste0("https://www.aoml.noaa.gov/hrd/Storm_pages/",storm,"/sonde.html"))
    lnks <- page |> rvest::html_nodes("a") |>rvest::html_attr("href")
    lnks <- lnks[grep("FRD",lnks)]
    lnks
  }
  try_retry <- function(url, destfile, ..., maxcount = 5) {
    cat("\nDownloading dropsonde data from: https://www.aoml.noaa.gov/hrd/data_sub/dropsonde.html")
    count <- 0
    repeat{
      Sys.sleep(0.5)
      try(download.file(url, destfile, ...),silent=TRUE)
      count <- count + 1
      if (file.exists(destfile) || count >= maxcount)
        break
    }
  }
  if (length(fls)==0) lnks <- tryURL(storm) else lnks<-gsub(paste0(storm,"_"),"",fls)
  if (length(lnks)==1){ cat(paste0("\nNo Sonde data detected for storm: ",storm));return(FALSE)}
  sondedata <- lapply(seq_along(lnks),function(x,lnks,trtry,dd,st){
    if (is.null(dd)) tf <- tempfile() else tf=paste0(dd,"/",st,"_",basename(lnks[x]))
    if (!file.exists(tf)&checksondes) trtry(lnks[x],tf,maxcount=3)
    if (!file.exists(tf)) return(NULL)
    dat <- readLines(tf)
    sondes <- grep("DROPWINDSONDE PROCESSING RECORD",dat)
    ends <- which(dat=="")
    headers <- grep("IX",dat)
    times <- grep("Aircraft",dat)
    dats=lapply(seq_along(sondes),function(i,s,h,d,t,e){
      sonde <- strsplit(d[s[i]],"Sonde: ")[[1]]
      if (length(sonde)==2) sonde <- sonde[[2]] else sonde=NA
      datetime = scan(text=d[t[i]],what="",quiet=TRUE)
      datetime = as.POSIXct(paste0(datetime[grep("Date:",datetime)+1],datetime[grep("Time:",datetime)+1]),format="%y%m%d %H%M%S",tz="UTC")
        if(is.na(datetime)) browser()
      if (i==length(s)){
        if (length(which(e>h[i]))>0) end=e[which(e>h[i])[1]]-1
        else end = length(d)
      }else{
        end =s[i+1]-2
      }
      dat <- data.frame(read.table(text=d[(h[i]+1):(end)],header=FALSE,
                                   col.names=c("IX","t.s","P.mb","T.C","RH.pc","Z.m","WD","WS.ms","U.ms","V.ms","NS ","WZ.ms","ZW.m","FP","FT","FH","FW","LAT.N","LON.E")))
      dat$sonde=s[i]
      dat$time = datetime+as.numeric(dat$t.s)
      dat
    },s=sondes,h=headers,d=dat,t=times,e=ends)
    dats <- do.call(rbind,dats)
    dats
  },lnks=lnks,trtry=try_retry,dd=ddir,st=storm)
  sondedata <- Filter(Negate(is.null), sondedata)
  if (length(sondedata)>0){
    sondedata <- do.call(rbind,sondedata)
    sondedata <- sondedata|>
      mutate(storm=storm) |>
      filter(LAT.N != -999,WS.ms!=-999,Z.m!=-999)|>
      sf::st_as_sf(coords=c("LON.E","LAT.N"),crs=4326)
    if (created) saveRDS(sondedata, donefile)
    return(sondedata)
  }else{
    cat("\nNo sonde data returned for this storm")
  }
}
