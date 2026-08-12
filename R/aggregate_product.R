aggregate_product <- function(rasts,type="wind",vars=NULL){
  aggs <- list()
  if (type=="wind"){
    if (is.null(vars)) vars=c("msw","pow","dur")
    for (v in vars){
      if (v=="msw"){
        aggs <- append(aggs,setNames(list(rast(lapply(rasts,function(x) app(x[[grepl("msw",names(x))]],fun="max",na.rm=TRUE)))),"max_msw"))
      }else if (v=="pow"){
        aggs <- append(aggs,setNames(list(rast(lapply(rasts, function(x) app(x[[grepl("msw",names(x))]],fun="sum",na.rm=TRUE)))),"tot_pow"))
      }else if(v=="dur"){
        dur <- lapply(rasts, function(x){
          msw <- x[[grep("msw",names(x))]]
          ###  we need to know the timestep, 0.5 hrs
          dt <- unique(difftime(time(msw),dplyr::lag(time(msw)),units="hours"))
          dt <- as.numeric(dt[!is.na(dt)])
          ### remove layers that never reach hurricane force winds
          msw <- msw[[minmax(msw)[2, ]>=33]]
          ###  reset the msw values less than Cat 1 hurricane  speed (33 m/s) to NA and the other to the time step
          dur1 <- lapply(msw,classify,matrix(c(0,33,NA,33,Inf,dt),ncol=3,byrow=TRUE), include.lowest=TRUE)
          dur1 <- app(rast(dur1),fun="sum",na.rm=TRUE)
          dur3 <- lapply(msw,classify,matrix(c(0,50,NA,50,Inf,dt),ncol=3,byrow=TRUE), include.lowest=TRUE)
          dur3 <- app(rast(dur3),fun="sum",na.rm=TRUE)
          list(dur1,dur3)
        })
        aggs <- append(aggs,setNames(unlist(dur),c("Cat1hrs","Cat3hrs")))
      }
    }
  }else if(type=="precip"){
    if (is.null(vars)) vars=c("tot","totpre")
    for (v in vars){
      if (v=="tot"){
        stormagg <- lapply(rasts,function(x) app(x$storm,fun="sum",na.rm=TRUE))
        if (length(stormagg)>1){ stormagg <- rast(stormagg)
        }else{ names(stormagg[[1]]) <- names(stormagg)[1];stormagg<-stormagg[[1]]}
        aggs <- append(aggs,setNames(list(stormagg),"tot_storm"))
      }else if (v=="totpre"){
        prehrs <- sapply(strsplit(unique(unlist(lapply(rasts, function(x) names(x$prestorm)))),"_"),"[[",3)
        if (length(prehrs)==1) prehrs=list(prehrs)
        preagg <- unlist(lapply(rasts, function(x) lapply(x$prestorm, app,fun="sum",na.rm=TRUE)))
        if (!grepl("tot_storm",names(aggs))){
          stormagg <- lapply(rasts,function(x) app(x$storm,fun="sum",na.rm=TRUE))
          if (length(stormagg)>1){ stormagg <- rast(stormagg)
          }else{ names(stormagg[[1]]) <- names(stormagg)[1];stormagg<-stormagg[[1]]}
        }
        aggs <- append(aggs,setNames(lapply(prehrs,function(h){
          rast(lapply(1:nlyr(stormagg),function(x){
            sum(stormagg[[x]],preagg[grepl(names(stormagg)[x],names(preagg))&grepl(h,names(preagg))][[1]],na.rm=TRUE)
          }))
        }),paste0("tot_stormprestorm-",prehrs)))
      }
    }
  }else{stop("non wind or rain aggregation no available yet.")}
  aggs
}
