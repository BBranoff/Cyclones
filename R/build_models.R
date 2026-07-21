#' Build models of storm characteristics from existing data.
#'@description
#'`build_models()` will accept the output of `get_storms()` and create a series of statistical models relating storm characteristics.
#'This includes the relationship between a given wind speed and its distance from the center, the eyewall radius as a function of the radius of maximum wind,
#'the pressure differential in a storm as a function of the maximum wind speed, & the radius of the last closed isobar as a function of the maximum radius of
#'given wind speeds. Models are nested by basin, Saffir-Simpson intensity category, and by storm quadrant. Models will only reflect the variability of input
#'storms. If input storms were previously consolidated, those values will be modeled, otherwise all values are modeled.
#'
#'These models are currently only used in the `make_extents()` function to fill in missing data when necessary, primarily to build wind extents for older storms.
#'
#' @param storms A list of storm tabular data as output by the `get_storms()` function. Ideally, this should include the largest possible pool of storms relevant for
#' subsequent analyses, as the models will reflect variability in these input storms only.
#' @param keep_data Whether to retain the training data with a prediction column for each set of models.
#'
#' @returns A list of nested tibbles, one for each of the model groups:
#' * distmods:  asymptotic and linear models of distance from the storm center as a function of wind speed. Models nested by basin, quadrant, and Saffir-Simpson category.
#' * emod:      linear models of eye radius as a function of the radius of maximum winds. Models nested by basin and Saffir-Simpson category.
#' * rocimod:   linear models of the radius of the last closed isobar (ROCI: the storm's maximum extent) as a function of the maximum radius of windspeeds.
#' Models nested by basin and Saffir-Simpson category.
#' * pocimod:   linear models of the pressure differential (the minimum pressure at the center and the atmospheric pressure at the ROCI) as a function of the maximum wind speed.
#' Models nested by basin and Saffir-Simpson category.
#'
#' For all models, variable units are native to those of USA IBTrACS variables: speed is windspeed in knots, dist is distance from center in nautical miles, and pressure is millibars.
#' @export
#' @examples
#'
#' #  build models off all IBTrACS storms
#' allstorms <- get_storms(ib_filt="ALL")
#' mods <- build_models(allstorms)
#'
#'# plot training and prediction data
#'distmoddat <- mods$distmods |> unnest(data)
#'samp <- sample(seq(1:nrow(distmoddat)),1000)
#'library(ggplot2)
#'ggplot(data=distmoddat[samp,],aes(x=windspeed,y=dist,col=factor(USA_SSHS),shape=BASIN))+geom_point()+geom_line(aes(y=pred_asymp.dist))+geom_line(aes(y=pred_lm.dist),lty=2)+facet_wrap(~BASIN)
#'
#' @importFrom dplyr select mutate rename filter bind_rows if_else across group_by ungroup summarise nest_by reframe distinct group_nest left_join cur_data
#' @importFrom tidyr pivot_longer unnest
#' @importFrom purrr map map2 safely
#' @importFrom stats dist lm nls predict setNames
build_models <- function(storms=NULL,strip=FALSE){
  if (is.null(storms)) stop("No data supplied")#;storms <- get_storms(source="ncei",ib_filt=NULL)
  if ('vctrs_list_of' %in% class(storms)){
    storms <- lapply(storms,data.frame) |>
      bind_rows()
  }
  storms <- storms|>
    mutate(across(contains(c("USA_SSHS","EYE","WIND","RMW","PRES","ROCI","POCI","R34","R50","R64")),as.numeric))

  storms_long <- storms |>
     pivot_longer(cols=contains(c("WIND","PRES","RMW","EYE","ROCI","POCI","R34","R50","R64")),names_to = c("source","var"),names_pattern="([^_]+)_(.*)")

  ############################################################
  #######   create models needed later
  ######    We need these to reconstruct maximum wind extents from given data
  ########################################################
  ###  a non-linear model of the wind speed distance as a function of wind speed and storm size by quadrant
  ###  in general, the SW quadrant has higher winds closer to the center, and the NE quadrant has the same wind speeds farther out
  ###  storm size also seems to be important
  ###  Asymptotic models by storm size and quadrant
  ###  general form of the models: In the following steps, these are done for specific subsets of data
  #modquad <- lm(data=quaddist|>filter(USA_SSHS>=1),dist~log(speed):(quad*BASIN))
  #modquad_asymp <- nls(dist ~ SSasymp(speed, Asym,R0, lrc), data = quaddist|> filter(USA_SSHS>=1),control = list(maxiter = 500))
  ###  try the model grouped by storm size (saffir simpson)
  ###  create a new dataset to predict from the model, with all possibilities of wind speed
  ###  create a 'safe' version of nls that will return NA instead of aborting
  nls_safe <- safely(nls,otherwise=NA)
  lm_safe <- safely(lm,otherwise=NA)
  ###  catch the missing models before predicting
  predict_safe <- function(mod,dat) if(length(mod)==1) return(NA) else(predict(mod,newdata=dat))
  #quaddist_new <- quaddist |>
  #  filter(USA_SSHS>=1)|>
  #  distinct(USA_SSHS,quad,BASIN)|>
  #  group_by(USA_SSHS,quad,BASIN) |>
  #  reframe(windspeed=1:185)|>
  #  nest_by(USA_SSHS,quad,BASIN)
  ##  build models on the original data, nested by both size, quadrant, and BASIN
  distdat <- storms_long|> filter(USA_SSHS>=0,grepl("WIND|RMW|NE|SE|SW|NW",var)) |>
    filter(if(any(grepl("CONS",source))) grepl("CONS",source),!is.na(value))
  distdat <- bind_rows(
    distdat |> filter(var %in% c("WIND","RMW")) |> tidyr::pivot_wider(names_from=var) |> filter(!is.na(RMW)) |> rename(windspeed=WIND,dist=RMW)|>mutate(quad=list(c("NE","SE","SW","NW")))|>
      unnest(quad),
    distdat |> filter(!var %in% c("WIND","RMW")) |>
      mutate(windspeed=as.numeric(gsub("R","",sapply(strsplit(var, "_"), `[`, 1))),
             quad=sapply(strsplit(var, "_"), `[`, 2)) |>
      rename(dist=value)
    )
   distmods <- distdat |> group_nest(USA_SSHS,quad,BASIN) |>
      mutate(#data2=quaddist_new$data,  ###  predict fitted models onto the new dataset with the full range of windspeeds. use the $result component of the _safe methods
           model_asymp = map(data,~nls_safe(dist ~ SSasymp(windspeed, Asym,R0, lrc), data = .,control = list(maxiter = 500))$result),
           model_lm = map(data,~lm_safe(dist ~log(windspeed), data = .)$result),
           ## add conditionals if the safe models return NA
           data=map2(model_asymp,data,~mutate(.y,pred_asymp.dist=predict_safe(.x,.y))),
           data=map2(model_lm,data,~mutate(.y,pred_lm.dist=predict_safe(.x,.y))))
   ##  combine with all quadrant model
   distmods <- bind_rows(distmods,
                     distdat |> group_nest(USA_SSHS,BASIN) |>
                       mutate(quad="all",
                              model_asymp = map(data,~nls_safe(dist ~ SSasymp(windspeed, Asym,R0, lrc), data = .,control = list(maxiter = 500))$result),
                              model_lm = map(data,~lm_safe(dist ~log(windspeed), data = .)$result),
                              ## add conditionals if the safe models return NA
                              data=map2(model_asymp,data,~mutate(.y,pred_asymp.dist=predict_safe(.x,.y))),
                              data=map2(model_lm,data,~mutate(.y,pred_lm.dist=predict_safe(.x,.y)))))
   ## combine with all quadrant, all basin model
   distmods <- bind_rows(distmods,
                     distdat |> group_nest(USA_SSHS) |>
                       mutate(quad="all",BASIN="all",
                              model_asymp = map(data,~nls_safe(dist ~ SSasymp(windspeed, Asym,R0, lrc), data = .,control = list(maxiter = 500))$result),
                              model_lm = map(data,~lm_safe(dist ~log(windspeed), data = .)$result),
                              ## add conditionals if the safe models return NA
                              data=map2(model_asymp,data,~mutate(.y,pred_asymp.dist=predict_safe(.x,.y))),
                              data=map2(model_lm,data,~mutate(.y,pred_lm.dist=predict_safe(.x,.y)))))
   ## combine with non nested data
   distmods <- bind_rows(distmods,
                     distdat |> group_nest()|>
                       mutate(quad="all",BASIN="all",USA_SSHS=999,
                              model_asymp = map(data,~nls_safe(dist ~ SSasymp(windspeed, Asym,R0, lrc), data = .,control = list(maxiter = 500))$result),
                              model_lm = map(data,~lm_safe(dist ~log(windspeed), data = .)$result),
                              ## add conditionals if the safe models return NA
                              data=map2(model_asymp,data,~mutate(.y,pred_asymp.dist=predict_safe(.x,.y))),
                              data=map2(model_lm,data,~mutate(.y,pred_lm.dist=predict_safe(.x,.y)))))
  ### take the predicted and original values to plot with
  if (any(is.na(distmods$model_asymp))) warning("Empty models for some Basin:Category:Quadrant groups. Insuffucient data the likely cause for those groups.")
   #####  visualize
   # ggplot(distmods %>% unnest(data) %>% group_by(USA_SSHS,quad,BASIN,windspeed) %>%
   #          summarise(dist=mean(dist,na.rm=T),
   #                    fit_lm=mean(pred_lm.dist),
   #                    fit_asymp=mean(pred_asymp.dist)),
   #        aes(x=windspeed,y=dist,col=factor(quad),shape=BASIN))+
   #   geom_point()+
   #   facet_wrap(~USA_SSHS)+
   #   geom_line(aes(x=windspeed,y=fit_asymp,lty="asympt."))+
   #   #geom_line(aes(x=speed,y=fit_lm,lty="linear"))+
   #   scale_linetype_manual(values=c(1,2))+
   #   labs(lty="model type",col="Quadrant")+
   #   ggtitle("North Atlantic Tropical Cyclone Wind Speed Distance by Storm Saffir-Simpson Scale\nGrouped models")+
   #   xlab("Wind Speed (kts)")+ylab("Max. distance from center (nmi)")

  ####  the non-linear model is a better fit for all combinations of quad and saffir simpson
  #####
  ##  also need to model the eyewall radius
  ##  this is crucial for maximum wind speeds and capturing the intense transition from eye to eye wall
  ####
   eyedat <- distdat |>
     dplyr::right_join(storms_long|>filter(var=="EYE",!is.na(value))|>select(SID,ISO_TIME,value),by=join_by(SID,ISO_TIME))|>
     group_by(SID,ISO_TIME) |>
     filter(dist==min(dist,na.rm=T))|>
     rename(minwinddist = dist)|>
     mutate(eyedist=as.numeric(value)/2) |>
     ungroup()
  eyemod <- eyedat|>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(eyedist~0 + minwinddist,data=.)),
           data=map2(model,data,~mutate(.y,pred.eyedist=predict(.x,.y))))
  # ggplot(data=eyemod %>% tidyr::unnest(data),
  #        aes(x=minwinddist,y=eyedist,col=BASIN))+
  #   geom_point()+
  #   facet_wrap(~USA_SSHS)+
  #   #stat_smooth(method="lm",formula="y~x+0")+
  #   geom_line(aes(y=pred.eyedist),lty=2)+
  #   geom_abline(slope=1,intercept=0)+
  #   geom_vline(xintercept=0)+
  #   geom_hline(yintercept=0)+
  #   ggtitle("North Atlantic Tropical Cyclone Eye Wall Distance and Distance of Maximum Wind by Storm Saffir-Simpson Scale\nGrouped models")+
  #   xlab("Minimum Distance of Maximum Wind (nmi)")+ylab("Eyewall Radius (nmi)")+
  #   theme_bw()


  ##  and storms size
  ####
  ###  for larger storms, the size can be reasonably predicted from the largest known wind swath
  rocidat <- storms_long|>
    filter(var=="ROCI",source=="CONS",!is.na(value))|>
    rename(ROCI=value)|>
    left_join(storms_long |> filter(grepl("_NE|_SE|_SW|_NW",var))|>
                group_by(var)|>
                filter(if("CONS" %in% source) source=="CONS" else TRUE)|>
                ungroup()|>
               tidyr::pivot_wider(id_cols=c(ID,ISO_TIME),names_from=var),by=join_by(ID,ISO_TIME))|>
    tidyr::pivot_longer(cols=contains(c("_NE","_SE","_SW","_NW")),values_to="dist")|>
    filter(!is.na(dist))|>
    group_by(ID,ISO_TIME) |>
    filter(dist==max(dist,na.rm=T))|>
    ungroup()
  ROCImod <- rocidat |>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(ROCI~0+dist,data=.)),
           data=map2(model,data,~mutate(.y,pred.roci=predict(.x,.y))))
  # ggplot(data=ROCImod %>% tidyr::unnest(data),
  #        aes(x=dist,y=ROCI,col=BASIN))+
  #   geom_point(alpha=0.5)+
  #   facet_wrap(~USA_SSHS)+
  #   #stat_smooth(method="lm",formula="y~x")+
  #   geom_line(aes(y=pred.roci),lty=2,size=1.5)+
  #   geom_abline(slope=1,intercept=0)+
  #   ggtitle("North Atlantic Tropical Cyclone Outer Isobar Distance and Maximum Swath Distance\nGrouped models")+
  #   xlab("Maximum Distance of Minimum Wind (nmi)")+ylab("Radius of Last Closed Isobar (nmi)")+
  #   xlim(0,1000)+ylim(0,1000)+
  #   theme_bw()

  pocidat <-  storms_long|> filter(grepl("POCI|PRES|WIND",var))|>
    group_by(var)|>
    filter(if("CONS" %in% source) source=="CONS" else TRUE)|>
    ungroup()|>
    tidyr::pivot_wider(id_cols=c(ID,ISO_TIME,USA_SSHS,BASIN),names_from=var)|>
    group_by(ID,ISO_TIME) |>
    mutate(Pressdif = POCI-PRES)|>
    filter(!is.na(Pressdif),Pressdif>0)|>
    ungroup()

  ###  also need to predict POCI when it is absent
  POCImod <- pocidat |>
    group_nest(USA_SSHS,BASIN) |>
    mutate(model=map(data,~lm(Pressdif~WIND,data=.)),
           data=map2(model,data,~mutate(.y,pred.pressdiff=predict(.x, .y))))

  # ggplot(POCImod %>% tidyr::unnest(data),
  #        aes(x=WIND,y=Pressdif,col=BASIN))+
  #   geom_point()+
  #   facet_wrap(~USA_SSHS)+
  #   #stat_smooth(method="lm",formula=y~x)+
  #   geom_line(aes(y=pred.pressdiff),lty=2)+
  #   ggtitle("North Atlantic Tropical Cyclone Pressure Difference and Maximum Wind Speed\nGrouped models")+
  #   xlab("Maximum Wind Speed (kts)")+ylab("Storm Pressure Difference (mbar)")
  if (strip){
    distmods <- distmods |> select(-data)|>
      mutate(model_asymp=map(model_asymp,stripmod),
             model_lm = map(model_lm,stripmod))
    distmods <- list(model_asymp = distmods$model_asymp,model_lm=distmods$model_lm,dat=distmods|>select(-c(model_asymp,model_lm)))
    eyemod <- eyemod |> select(-data)|>
      mutate(model=map(model,stripmod))
    eyemod <- list(model = eyemod$model,dat=eyemod|>select(-model))
    ROCImod <- ROCImod|> select(-data)|>
      mutate(model=map(model,stripmod))
    ROCImod  <- list(model = ROCImod$model,dat=ROCImod|>select(-model))
    POCImod <- POCImod |>select(-data)|>
      mutate(model=map(model,stripmod))
    POCImod  <- list(model = POCImod$model,dat=POCImod|>select(-model))
  }
  ## also need the average minimum pressure by month, basin, and storm size for rare cases when missing
  minpress <- storms_long |> filter(var=="PRES")|>
    mutate(MONTH=format.Date(ISO_TIME,"%m")) |>
    filter(if("CONS" %in% source) source=="CONS" else TRUE)|>
    group_by(BASIN,USA_SSHS,MONTH) |>
    summarise(PRES=mean(value,na.rm=TRUE))
  mods <- list(distmods=distmods,#|> select(-c(data,pred_asymp.dist,pred_lm.dist)),
               emod=eyemod,#|> select(-c(data,pred)),
               rocimod=ROCImod,#|> select(-c(data,pred)),
               pocimod=POCImod,#|> select(-c(data,pred)),
               minpress = minpress)

 mods
}
stripmod = function(cm) {
  cm$y = c()
  cm$model = c()

  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()


  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()

  cm
}
