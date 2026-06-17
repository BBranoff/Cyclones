checkrasts <- function(rasts){
  if (!unique(unlist((lapply(rasts,class)))) %in% c("SpatRaster","character")) stop("Unknown raster inputs. Supply either storm raster bricks or filenames.")
  nmes <- names(rasts)
  if (is.null(nmes)){
      nmes <- unlist(lapply(rasts,names))
      if (is.null(nmes)){
          nmes <- unique(unlist(lapply(rasts,function(x){paste(strsplit(basename(x),"_")[[1]][2:5],collapse="_")})))
          meths <- lapply(rasts,function(x){strsplit(basename(x),"_")[[1]][1]})
          rasts <- unlist(rasts)
          rasts <- setNames(lapply(nmes, function(x) rasts[grep(x,rasts)]),unique(nmes))
      }else{
        rasts <- lapply(rasts,function(x) x[[grep("msw|windir",names(x))]])
        nmes <- unique(unlist(lapply(rasts,function(x){paste(strsplit(basename(sources(x)),"_")[[1]][2:5],collapse="_")})))
        meths <- lapply(rasts,function(x){strsplit(basename(x),"_")[[1]][1]})
        rasts <- setNames(rasts,unique(nmes))#setNames(sapply(rasts,"[[",1),unique(nmes))
      }
  }else{
    ID <- sapply(strsplit(nmes,"_"),"[[",4)
    ##  if the filenames were unlisted somehow
    ##  re group them
    #if (all(substr(nmes[1:9],nchar(nmes[1:9]),nchar(nmes[1:9]))==as.character(c(1:9))))
    if (nchar(ID[1])>13){
      ID <- substr(ID,1,13)
      nmes <- strsplit(names(rasts),"_")
      nmes <- unlist(lapply(seq_along(nmes),function(x,n,id){paste(c(n[[x]][1:3],id[x]),collapse="_")},n=nmes,id=ID))
      rasts <- setNames(lapply(unique(nmes), function(x) rasts[nmes==x]),unique(nmes))
    }
  }
  ####  now load rasters if only the filenames provided
  if (unique(unlist((lapply(rasts,class))))=="character"){
    cat("\nloading rasters...")
    rasts <- lapply(rasts,rast)
  }
  rasts
}
