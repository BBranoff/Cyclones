#' Consolidate storm tabular data.
#'
#'`cons_stormdat()` will accept a singular data frame sourced from IBTrACS or a listed data set produced by 'get_storms()' and will consolidate
#'the desired variables across the various agencies into one column. In the case of multiple values for a given variable, either the preferred agency
#'is used or a function is applied to summarise the values.
#' @param dat A singular dataframe or listed data set, sourced from IBTrACS and containing multiple analogous values from different agencies.
#' @param vars Short names for the variables to be consolidated. Values of 'wind', 'pres', and 'rmw' correspond to the maximum sustained
#' wind, the minimum central pressure, and the radius of maximum wind, respectively. A value of 'quads' refers to the quadrant specific
#' wind extents of 34, 50, and 64 knot winds in the NE, SE, SW, and NW quadrants.
#' @param msw_int The target interval for the maximum sustained winds, which is not consistent among the agencies. The USA uses a one-minute
#' interval, most other agencies use a 10-minute interval, but 2-minute and 3-minute intervals are also used. When consolidating the wind columns,
#' linear models are used to translate any values whose interval is not the target interval. The models are computed from all rows in which multiple
#' interval values are provided.
#' @param pref Agency preferences for all values, listed in order of preference. If these values are present, they will be used (after conversion to the
#' target interval). These can be any of the agency or model abbreviations in the IBTrACS data (USA, WMO, TOKYO, CMA, HKO, KMA, NEWDELHI, REUNION, BOM,
#' NADI, WELLINGTON, DS824, TD9636, TD9635, NEUMANN, MLC). However, because storms sometimes traverse multiple agency jurisdictions, WMOfirst or WMOlast
#' can also be used to avoid sudden changes in values.
#' @param fun The summary function to be applied to any instances of multiple values.
#' @returns A data frame or listed data set, whatever the input was, with new columns whose names are prefixed with 'CONS_' representing the consolidated values.
#' @importFrom dplyr if_any contains arrange
#' @importFrom tidyr fill
cons_stormdat <- function(dat,vars=c("wind","pres","rmw","quads","roci","poci","eye"),msw_int="1min",pref=c("USA","WMO"),fun="mean"){

  if ('vctrs_list_of' %in% class(dat)){listed=TRUE;dat <- lapply(dat,data.frame)|> bind_rows()}else{listed=FALSE}
  dat <- dat|>mutate(tempn = seq(1:n()))
  if (sum(grepl("WMO",pref))>1) stop("Too many WMO preference inputs. Choose either 'WMO', 'WMOfirst', or 'WMOlast'.")
  if ("quads"%in%vars) vars <- c(vars[-which(vars=="quads")],as.vector(outer(c("R34_","R50_","R64_"),c("NE","SE","SW","NW"),paste0)))
  for (v in vars){
    varstring <- paste0("_",toupper(v))
    varcols <-  dat|>
      select(contains(c("tempn","ID","BASIN","WMO_AGENCY",varstring)))|>
      ##  set ID to ensure we get everything and reassemble
      mutate(across(everything(),~if_else(.x==" ",NA,.x)),
             across(contains(varstring),~if_else(as.numeric(.x)<0,NA,as.numeric(.x))))
    ## storm observations with no values
    ## cant do anything with these for now
    novals <- varcols |>
      slice(which(rowSums(is.na(varcols|>select(contains(varstring))))==sum(grepl(varstring,names(varcols)))))
    ### remaining are either single or multiple vals
    ### exclude WMO because that will always be a repeat of another value
    multvals <- varcols |>
      select(-contains("WMO"))|>
      slice(which(rowSums(!is.na(varcols|>select(contains(varstring))))>1))
    ### use these to build the inter-interval models for wind
    ### if a preferred source is given, that is the target value for the model
    ###  otherwise all of the desired interval columns are modelled together
    if (grepl("wind",v,ignore.case = TRUE)&&!is.null(msw_int)) interval_mods <- do_interval_mods(multvals,msw_int,pref)
    else interval_mods <- NULL
    vals <- do_pref(varcols,varstring,msw_int,pref,fun, interval_mods)
    vals <- bind_rows(novals,vals)|>select(tempn,contains(varstring))
    dat <- dat |> left_join(vals |> select(tempn,contains("CONS")),by=join_by(tempn))
  }
  if (listed){
    dat <- dat |> group_by(ID)
    IDs <- group_keys(dat)|>pull(ID)
    dat <- dat |>group_split()
    names(dat) <- IDs
  }
  dat
}
do_pref <- function(df,varstr,target_int,prf,fn,mod){
  if (!is.null(prf)&&any(grepl("WMO",prf))){
    if (!"WMO_AGENCY" %in% names(df)) df$WMO_AGENCY <- "USA"
    else df$WMO_AGENCY <- toupper(df$WMO_AGENCY)
    df <- df |>
      mutate(WMO_AGENCY= if_else(WMO_AGENCY %in% c("HURDAT_ATL","ATCF","HURDAT_EPA","CPHC"),"USA",
                                 if_else(WMO_AGENCY==" ",NA,WMO_AGENCY)))|>
      group_by(ID)
    ###  fill in the agency
    if (any(prf=="WMO")){ df <- df|> fill(WMO_AGENCY,.direction="down");df <- df|> fill(WMO_AGENCY,.direction="up")|>ungroup()}
    else if (any(prf=="WMOfirst")) df <- df|> mutate(WMO_AGENCY=first(WMO_AGENCY[!is.na(WMO_AGENCY)]))|>ungroup()
    else if (any(prf=="WMOlast")) df <- df|> mutate(WMO_AGENCY=last(WMO_AGENCY[!is.na(WMO_AGENCY)]))|>ungroup()
  }
  if (varstr=="_WIND"){
    cols10min <- c("TOKYO_WIND","HKO_WIND","KMA_WIND","REUNION_WIND","BOM_WIND","NADI_WIND","WELLINGTON_WIND")
    cols1min <- c("USA_WIND","DS824_WIND","TD9636_WIND","TD9635_WIND","NEUMANN_WIND","MLC_WIND")
    cols3min <- "NEWDELHI_WIND"
    cols2min <- "CMA_WIND"
    ###  convert each column to the target wind interval
    for (df_int in c("10min","1min","2min","3min")) {
      if(df_int=="10min") cols <-cols10min
      else if(df_int=="1min") cols <- cols1min
      else if(df_int=="2min") cols=cols2min
      else if (df_int=="3min") cols= cols3min
      ###  only do the columns that exists in the dataset
      cols <- grep(paste(names(df),collapse="|"),cols,value=TRUE)
      for (C in cols){
        ###  if the column is in another time interval group, model the target interval groups
        ###  !is.null(mod[[df_int]]) should only be true if the columns exist, checked in do_interval_mods
        if (!target_int==df_int&&!is.null(mod[[df_int]])){
          m <- mod[[df_int]]
          df[!is.na(df[,C]),C] <- round(predict(m,newdata=data.frame(x=df[!is.na(df[,C]),]|>pull({{C}}))))
        }
        ###  sometimes the WMO_WIND is missing, even though the wind from the corresponding WMO agency is not. fill these in.
        if (!is.null(prf)&&any(grepl("WMO",prf))){
          df$WMO_WIND[!is.na(df$WMO_AGENCY)&!is.na(df[,C])&df$WMO_AGENCY==gsub("_WIND","",C)] <- df[!is.na(df$WMO_AGENCY)&!is.na(df[,C])&df$WMO_AGENCY==gsub("_WIND","",C),]|>pull({{C}})
        }
      }
    }
  }
  if (!is.null(prf)){
    prefcol <- grep(varstr,names(df),ignore.case=TRUE,value=TRUE)
    ##  to maintain the original order of the preferences, must loop through them individually
    prefcol <- lapply(prf,function(p) grep(p,prefcol,ignore.case=TRUE,value=TRUE))|>unlist()
  }
  ###  for multivals
  multvals <- df |>
    slice(which(rowSums(!is.na(df|>select(-contains("WMO"))|>select(contains(varstr))))>1))
  consvar <- paste0("CONS",varstr)
  if (nrow(multvals)>0){
    if (fn=="mean"){
      multvals[,consvar] <- round(rowMeans(multvals|>select(-contains("WMO"))|>select(contains(varstr)),na.rm=TRUE))
    }else if(fn=="max"){
      multvals[,consvar] <- round(do.call(pmax,c(multvals|>select(-contains("WMO"))|>select(contains(varstr)),na.rm=TRUE)))
    }else if(fn=="min"){
      multvals[,consvar] <- round(do.call(pmin,c(multvals|>select(-contains("WMO"))|>select(contains(varstr)),na.rm=TRUE)))
    }else{
      cat(paste0("'",fn,"' not recognized, returning NA for non preferential values."))
      multvals[,consvar] <- NA
    }
    ###  if a preference was given replace the mean with that value
    ###  if multiple were given, start with the first and work down the list
    if (!is.null(prf)&&length(prefcol)>0){
      multvals <- multvals |>
        mutate({{consvar}} := coalesce(!!!syms(c(prefcol,consvar))))
    }
  }
  ### for single vals
  if (any(grepl(varstr,names(df)))) {
    singlevals <- df |>
      select(-contains("WMO"))|>
      slice(which(rowSums(!is.na(across(contains(varstr))))==1))|>
      mutate("{consvar}" := do.call(coalesce,across(contains(varstr))))
    ###  a handful of WMO vals are the only vals for some reason, add those in too
    ###  only if WMO for the current variable is available
    if (paste0("WMO",varstr) %in% names(df)){
      singlevals <- bind_rows(singlevals,
                              df |>
                                slice(which(rowSums(!is.na(across(contains(varstr))))==1&!is.na(across(contains(paste0("WMO",varstr))))))|>
                                mutate("{consvar}" := !!!syms(paste0("WMO",varstr))))
    }
    vals <- bind_rows(multvals,singlevals)
  }else{
    vals <- multvals
  }
  vals
}
do_interval_mods <- function(df,target_int,prf){
  cols10 <- c("TOKYO_WIND","HKO_WIND","KMA_WIND","REUNION_WIND","BOM_WIND","NADI_WIND","WELLINGTON_WIND")
  cols1 <- c("USA_WIND","DS824_WIND","TD9636_WIND","TD9635_WIND","NEUMANN_WIND","MLC_WIND")
  cols2 <- "CMA_WIND"
  cols3 <- "NEWDELHI_WIND"
  cols <- list("10"=cols10,"1"=cols1,"2"=cols2,"3"=cols3)
  predcols <- cols[-which(gsub("min","",target_int)==c("10","1","2","3"))]
  predcols <- lapply(predcols,function(x) grep(paste(names(df),collapse="|"),x,value=TRUE))
  predcols <- Filter(Negate(function(x) length(x)==0),predcols)
  if (length(predcols)==0) return(NULL)
  respcols <- cols[which(gsub("min","",target_int)==c("10","1","2","3"))]
  mods <- list()
  for (p in 1:length(predcols)){
    ###  if a preference is given, just model that value
    ###  Note, if WMO is the preference, this cannot be modeled as they often represent multiple time intervals
    if(!is.null(prf)&&sum(grepl(paste(prf,collapse="|"),respcols[[1]]))>0){
      ###  even if multiple preference were given, just use the first as the target for the models
      r=respcols[[1]][grep(paste(prf,collapse="|"),respcols[[1]])]
      df_long <- df |> select(r,predcols[[p]])|>
        pivot_longer(cols=predcols[[p]])|>
        rename(y:={{r}},x=value)|>
        filter(!if_any(c(y,x),~is.na(.)))
      ##  otherwise model all the values in the time interval
    }else{
      df_long <- lapply(respcols[[1]], function(r){
        df |> select(r,predcols[[p]])|>
          pivot_longer(cols=predcols[[p]])|>
          rename(y:={{r}},x=value)|>
          filter(!if_any(c(y,x),~is.na(.)))
      })
      df_long <- do.call(rbind,df_long)
    }
    mod <- tryCatch({
      lm(y~x,data=df_long)
     }, error = function(e){
      NULL})
    mods <- append(mods,list(mod))
    names(mods)[length(mods)] <- paste0( names(predcols)[[p]],"min")
  }
  mods
}

