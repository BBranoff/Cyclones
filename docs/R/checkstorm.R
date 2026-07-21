#' Internal check to verify if input data is properly formatted
#'
#' @param dta to check formatting on. Should be the result of `get_storms()`
#'
#' @return a single data frame that is the same as the entered into parent function, or one of the listed data frames entered into the parent function.
#' @importFrom dplyr rename_with
checkstorm <- function(dta,agncy){
  options(warn = 1)
  if (paste0(class(dta),collapse="") %in% c("vctrs_list_ofvctrs_vctrlist","list")){
    if (length(dta)>1){
      #mult=TRUE
      classes <- lapply(dta, function(x) paste0(class(x),collapse=""))
      if (!all(unlist(classes)=="tbl_dftbldata.frame")) stop("Unknown data input for one or more storms. Use get_storms to retrieve and/or format data for further processing.")
      warning("More than one storm in data. Only first will be used. Use lapply or a parallel equivalent to repeat for multiple storms.")
    }#else{mult=FALSE}
    dta = dta[[1]]
  } else if(paste0(class(dta),collapse="")!="tbl_dftbldata.frame") {
    stop("Unknown data input. Use get_storms to retrieve and/or format data for further processing.")
  }
  if (!is.null(agncy)){
    if (!all(paste0(agncy,"_",c("WIND","PRES","RMW","R34_NE","R34_SE","R34_SW","R34_NW","R50_NE","R50_SE","R50_SW","R50_NW","R64_NE","R64_SE","R64_SW","R64_NW")) %in% names(dta))){
      stop(paste("The following required columns are missing: ",
                 paste(grep(paste(names(dta),collapse="|"),paste0(agncy,"_",c("WIND","PRES","RMW","R34_NE","R34_SE","R34_SW","R34_NW","R50_NE","R50_SE","R50_SW","R50_NW","R64_NE","R64_SE","R64_SW","R64_NW")),value=TRUE,invert = TRUE),collapse=","),
                 ". Change the 'agency' parameter or use 'cons_stormdat()' to consolidate and/or retain agency columns before running 'make_extents()'.",collapse=""))
    }else{
      dta <- dta |>
        rename_with(~ gsub(paste0(agncy,"_"), "", .x), any_of(paste0(agncy,"_",c("WIND","PRES","RMW","ROCI","POCI","EYE","R34_NE","R34_SE","R34_SW","R34_NW","R50_NE","R50_SE","R50_SW","R50_NW","R64_NE","R64_SE","R64_SW","R64_NW"))))
    }
  }
  return(dta)
}
