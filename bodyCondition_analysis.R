library(tidyverse)
library(readr)
library(mgcv)
library(MuMIn)
library(ggplot2)

custsum <- function(x){
  mu <- mean(x)
  quant <- quantile(x,c(.025,0.5,0.975))
  std <- sd(x)
  out <- c("mean"=mu,"sd"=std,quant)
  return(out)
}

adjusted <- FALSE # TRUE or FALSE for adjusted or unadjusted measurements
observers <- c("ELR","GMB","HLZ") # observers to include
means <- TRUE # use observers means? TRUE or FALSE
title_var <- paste0("UAS Body Condition ",ifelse(adjusted,"adjusted","unadjusted"),ifelse(length(observers)<3,paste0(" ",paste0(observers,collapse=" & ")),""),ifelse(means & length(observers)>1," mean",""))

rawData <- read_csv("data/BodyCondition_ModelDevelopmentExport_20241017_SMK.csv",na=c("","NULL"))
modData <- rawData %>% select(target_identifier,field_efforts_id,measured_by,adjusted_image_dt,scale_bar_cm,species,sex,age_class,age,adjusted_range_m,mass_kg_during_mean,
                              paste0("total_length_cm",ifelse(adjusted,"_adj","")),
                              paste0("total_length_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_01_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_01_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_02_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_02_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_03_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_03_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_04_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_04_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_05_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_05_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_06_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_06_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_07_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_07_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_08_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_08_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_09_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_09_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("width_10_cm",ifelse(adjusted,"_adj","")),
                              paste0("width_10_cm",ifelse(adjusted,"_adj",""),"_s"),
                              paste0("ax_width_cm",ifelse(adjusted,"_adj","")),
                              paste0("max_width_cm",ifelse(adjusted,"_adj",""))) %>% 
  rename(total_length=paste0("total_length_cm",ifelse(adjusted,"_adj","")),
         total_length_s=paste0("total_length_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_01=paste0("width_01_cm",ifelse(adjusted,"_adj","")),
         width_01_s=paste0("width_01_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_02=paste0("width_02_cm",ifelse(adjusted,"_adj","")),
         width_02_s=paste0("width_02_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_03=paste0("width_03_cm",ifelse(adjusted,"_adj","")),
         width_03_s=paste0("width_03_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_04=paste0("width_04_cm",ifelse(adjusted,"_adj","")),
         width_04_s=paste0("width_04_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_05=paste0("width_05_cm",ifelse(adjusted,"_adj","")),
         width_05_s=paste0("width_05_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_06=paste0("width_06_cm",ifelse(adjusted,"_adj","")),
         width_06_s=paste0("width_06_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_07=paste0("width_07_cm",ifelse(adjusted,"_adj","")),
         width_07_s=paste0("width_07_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_08=paste0("width_08_cm",ifelse(adjusted,"_adj","")),
         width_08_s=paste0("width_08_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_09=paste0("width_09_cm",ifelse(adjusted,"_adj","")),
         width_09_s=paste0("width_09_cm",ifelse(adjusted,"_adj",""),"_s"),
         width_10=paste0("width_10_cm",ifelse(adjusted,"_adj","")),
         width_10_s=paste0("width_10_cm",ifelse(adjusted,"_adj",""),"_s"),
         ax_width=paste0("ax_width_cm",ifelse(adjusted,"_adj","")),
         max_width=paste0("max_width_cm",ifelse(adjusted,"_adj","")))

modData$field_efforts_id <- factor(modData$field_efforts_id)
modData$measured_by <- factor(modData$measured_by)
modData$species <- relevel(factor(modData$species),"Ringed seal")
modData$age_class <- factor(modData$age_class)
modData$sex <- factor(modData$sex)
modData$season <- factor(ifelse(month(modData$adjusted_image_dt)<=8,"low","high")) # low reserves before September, high reserves after
modData$monthNum <- month(modData$adjusted_image_dt)
modData$month <- factor(modData$monthNum)
modData$year <- factor(year(modData$adjusted_image_dt))
modData$ID <- modData$target_identifier
for(i in unique(modData$target_identifier))
  modData$ID[which(modData$target_identifier==i)] <- paste0(i,modData$year[which(modData$target_identifier==i)],modData$season[which(modData$target_identifier==i)])
modData$target_identifier <- factor(modData$target_identifier)
modData$replicate <- factor(as.numeric(factor(modData$adjusted_image_dt)))
if(isTRUE(observers=="random")){
  modData <- modData %>% group_by(replicate) %>% sample_n(1) %>% ungroup()
} else {
  if(length(observers)<length(unique(modData$measured_by))) modData <- modData %>% filter(measured_by %in% observers)
  if(means & length(observers)>1){
    modData <- modData %>% group_by(replicate) %>% mutate(total_length=mean(total_length)) %>% 
      mutate(ax_width=mean(ax_width)) %>% 
      mutate(max_width=mean(max_width)) %>% 
      mutate(width_01=mean(width_01)) %>% 
      mutate(width_02=mean(width_02)) %>% 
      mutate(width_03=mean(width_03)) %>% 
      mutate(width_04=mean(width_04)) %>% 
      mutate(width_05=mean(width_05)) %>% 
      mutate(width_06=mean(width_06)) %>% 
      mutate(width_07=mean(width_07)) %>% 
      mutate(width_08=mean(width_08)) %>% 
      mutate(width_09=mean(width_09)) %>% 
      mutate(width_10=mean(width_10)) %>% 
      mutate(total_length_s=mean(total_length_s)) %>% 
      mutate(width_01_s=mean(width_01_s)) %>% 
      mutate(width_02_s=mean(width_02_s)) %>% 
      mutate(width_03_s=mean(width_03_s)) %>% 
      mutate(width_04_s=mean(width_04_s)) %>% 
      mutate(width_05_s=mean(width_05_s)) %>% 
      mutate(width_06_s=mean(width_06_s)) %>% 
      mutate(width_07_s=mean(width_07_s)) %>% 
      mutate(width_08_s=mean(width_08_s)) %>% 
      mutate(width_09_s=mean(width_09_s)) %>% 
      mutate(width_10_s=mean(width_10_s)) %>% 
      slice_head(n = 1) %>% ungroup()
  } 
}

id <- gam(mass_kg_during_mean ~ s(target_identifier, bs = "re"),family=gaussian(link="log"),data=modData,method="ML")

alt <- gam(mass_kg_during_mean ~ s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

id.alt <- gam(mass_kg_during_mean ~ s(target_identifier, bs = "re")+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

spec <- gam(mass_kg_during_mean ~ species,family=gaussian(link="log"),data=modData,method="ML")

length_s.spec <- gam(mass_kg_during_mean ~ species + total_length_s,family=gaussian(link="log"),data=modData,method="ML")

length.spec <- gam(mass_kg_during_mean ~ species + total_length,family=gaussian(link="log"),data=modData,method="ML")

lengthXspec <- gam(mass_kg_during_mean ~ total_length * species,family=gaussian(link="log"),data=modData,method="ML")

length_sXspec <- gam(mass_kg_during_mean ~ total_length_s * species,family=gaussian(link="log"),data=modData,method="ML")

length2.spec <- gam(mass_kg_during_mean ~ total_length + I(total_length^2) + species,family=gaussian(link="log"),data=modData,method="ML")

length2_s.spec <- gam(mass_kg_during_mean ~ total_length_s + I(total_length_s^2) + species,family=gaussian(link="log"),data=modData,method="ML")

lengthXspec.length2Xspec <- gam(mass_kg_during_mean ~ total_length * species + I(total_length^2) * species,family=gaussian(link="log"),data=modData,method="ML")

length_sXspec.length2_sXspec <- gam(mass_kg_during_mean ~ total_length_s * species + I(total_length_s^2) * species,family=gaussian(link="log"),data=modData,method="ML")

length.spec.alt <- gam(mass_kg_during_mean ~ total_length + species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

length_s.spec.alt <- gam(mass_kg_during_mean ~ total_length_s + species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

length.spec.alt <- gam(mass_kg_during_mean ~ total_length + species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

length_s.spec.alt <- gam(mass_kg_during_mean ~ total_length_s + species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

lengthXspec.alt <- gam(mass_kg_during_mean ~ total_length * species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

length_sXspec.alt <- gam(mass_kg_during_mean ~ total_length_s * species+s(adjusted_range_m),family=gaussian(link="log"),data=modData,method="ML")

ax.max.spec <- gam(mass_kg_during_mean ~ ax_width + max_width + species,family=gaussian(link="log"),data=modData,method="ML")

axXspec.maxXspec <- gam(mass_kg_during_mean ~ ax_width * species + max_width * species,family=gaussian(link="log"),data=modData,method="ML")

width.spec <- gam(mass_kg_during_mean ~ width_01+width_02+width_03+width_04+width_05+width_06+width_07+width_08+width_09+width_10+species,family=gaussian(link="log"),data=modData,method="ML")

widthXspec <- gam(mass_kg_during_mean ~ width_01* species+width_02* species+width_03* species+width_04* species+width_05* species+width_06* species+width_07* species+width_08* species+width_09* species+width_10* species,family=gaussian(link="log"),data=modData,method="ML")

width_s.spec <- gam(mass_kg_during_mean ~ width_01_s+width_02_s+width_03_s+width_04_s+width_05_s+width_06_s+width_07_s+width_08_s+width_09_s+width_10_s + species,family=gaussian(link="log"),data=modData,method="ML")

width_sXspec <- gam(mass_kg_during_mean ~ width_01_s* species+width_02_s* species+width_03_s* species+width_04_s* species+width_05_s* species+width_06_s* species+width_07_s* species+width_08_s* species+width_09_s* species+width_10_s* species,family=gaussian(link="log"),data=modData,method="ML")

length.ax.max.spec <- gam(mass_kg_during_mean ~ total_length +ax_width + max_width + species,family=gaussian(link="log"),data=modData,method="ML")

length_s.ax.max.spec <- gam(mass_kg_during_mean ~ total_length_s +ax_width + max_width + species,family=gaussian(link="log"),data=modData,method="ML")

lengthXspec.axXspec.maxXspec <- gam(mass_kg_during_mean ~ total_length * species+ax_width * species + max_width * species,family=gaussian(link="log"),data=modData,method="ML")

length_sXspec.axXspec.maxXspec <- gam(mass_kg_during_mean ~ total_length_s * species+ax_width * species + max_width * species,family=gaussian(link="log"),data=modData,method="ML")

length.width.spec <- gam(mass_kg_during_mean ~ total_length+width_01+width_02+width_03+width_04+width_05+width_06+width_07+width_08+width_09+width_10+species,family=gaussian(link="log"),data=modData,method="ML")

length.width.spec.id <- gam(mass_kg_during_mean ~ total_length+width_01+width_02+width_03+width_04+width_05+width_06+width_07+width_08+width_09+width_10+species+s(target_identifier, bs = "re"),family=gaussian(link="log"),data=modData,method="ML")

length_s.width_s.spec <- gam(mass_kg_during_mean ~ total_length_s+width_01_s+width_02_s+width_03_s+width_04_s+width_05_s+width_06_s+width_07_s+width_08_s+width_09_s+width_10_s+species,family=gaussian(link="log"),data=modData,method="ML")

length_s.width_s.spec.id <- gam(mass_kg_during_mean ~ total_length_s+width_01_s+width_02_s+width_03_s+width_04_s+width_05_s+width_06_s+width_07_s+width_08_s+width_09_s+width_10_s+species+s(target_identifier, bs = "re"),family=gaussian(link="log"),data=modData,method="ML")

lengthXspec.widthXspec <- gam(mass_kg_during_mean ~ total_length* species+width_01* species+width_02* species+width_03* species+width_04* species+width_05* species+width_06* species+width_07* species+width_08* species+width_09* species+width_10* species,family=gaussian(link="log"),data=modData,method="ML")

length_sXspec.width_sXspec <- gam(mass_kg_during_mean ~ total_length_s* species+width_01_s* species+width_02_s* species+width_03_s* species+width_04_s* species+width_05_s* species+width_06_s* species+width_07_s* species+width_08_s* species+width_09_s* species+width_10_s* species,family=gaussian(link="log"),data=modData,method="ML")

if(length(observers)>1 & means==FALSE){
  lengthXspec.widthXspec.obs <- gam(mass_kg_during_mean ~ total_length* species+width_01* species+width_02* species+width_03* species+width_04* species+width_05* species+width_06* 
                                          species+width_07* species+width_08* species+width_09* species+width_10* species+s(measured_by, bs = "re"),family=gaussian(link="log"),data=modData,method="ML")
  length_sXspec.width_sXspec.obs <- gam(mass_kg_during_mean ~ total_length_s* species+width_01_s* species+width_02_s* species+width_03_s* species+width_04_s* species+width_05_s* species+width_06_s* species+width_07_s* species+width_08_s* species+width_09_s* species+width_10_s* species+s(measured_by, bs = "re"),family=gaussian(link="log"),data=modData,method="ML")
}

if(length(observers)>1 & means==FALSE){
  aic <- AICc(alt,ax.max.spec,axXspec.maxXspec,id,id.alt,length_s.ax.max.spec,length_s.spec,length_s.spec.alt,length_s.width_s.spec,length_s.width_s.spec.id,length_sXspec,length_sXspec.alt,length_sXspec.axXspec.maxXspec,length_sXspec.length2_sXspec,length_sXspec.width_sXspec,length.ax.max.spec,length.spec,length.spec.alt,length.width.spec,length.width.spec.id,length2_s.spec,length2.spec,lengthXspec,lengthXspec.alt,lengthXspec.axXspec.maxXspec,lengthXspec.length2Xspec,lengthXspec.widthXspec,spec,width_s.spec,width_sXspec,width.spec,widthXspec,lengthXspec.widthXspec.obs,length_sXspec.width_sXspec.obs)
  bic <- BIC(alt,ax.max.spec,axXspec.maxXspec,id,id.alt,length_s.ax.max.spec,length_s.spec,length_s.spec.alt,length_s.width_s.spec,length_s.width_s.spec.id,length_sXspec,length_sXspec.alt,length_sXspec.axXspec.maxXspec,length_sXspec.length2_sXspec,length_sXspec.width_sXspec,length.ax.max.spec,length.spec,length.spec.alt,length.width.spec,length.width.spec.id,length2_s.spec,length2.spec,lengthXspec,lengthXspec.alt,lengthXspec.axXspec.maxXspec,lengthXspec.length2Xspec,lengthXspec.widthXspec,spec,width_s.spec,width_sXspec,width.spec,widthXspec,lengthXspec.widthXspec.obs,length_sXspec.width_sXspec.obs)
  
} else {
  aic <- AICc(alt,ax.max.spec,axXspec.maxXspec,id,id.alt,length_s.ax.max.spec,length_s.spec,length_s.spec.alt,length_s.width_s.spec,length_s.width_s.spec.id,length_sXspec,length_sXspec.alt,length_sXspec.axXspec.maxXspec,length_sXspec.length2_sXspec,length_sXspec.width_sXspec,length.ax.max.spec,length.spec,length.spec.alt,length.width.spec,length.width.spec.id,length2_s.spec,length2.spec,lengthXspec,lengthXspec.alt,lengthXspec.axXspec.maxXspec,lengthXspec.length2Xspec,lengthXspec.widthXspec,spec,width_s.spec,width_sXspec,width.spec,widthXspec)
  bic <- BIC(alt,ax.max.spec,axXspec.maxXspec,id,id.alt,length_s.ax.max.spec,length_s.spec,length_s.spec.alt,length_s.width_s.spec,length_s.width_s.spec.id,length_sXspec,length_sXspec.alt,length_sXspec.axXspec.maxXspec,length_sXspec.length2_sXspec,length_sXspec.width_sXspec,length.ax.max.spec,length.spec,length.spec.alt,length.width.spec,length.width.spec.id,length2_s.spec,length2.spec,lengthXspec,lengthXspec.alt,lengthXspec.axXspec.maxXspec,lengthXspec.length2Xspec,lengthXspec.widthXspec,spec,width_s.spec,width_sXspec,width.spec,widthXspec)
}

bicOrd <- order(bic$BIC)
aic <- aic[bicOrd,]
bic <- bic[bicOrd,]

exclude <- c('s(target_identifier)',"s(year)","s(season)","s(measured_by)","s(age)","s(sex)","s(sex):speciesBearded seal","s(sex):speciesRinged seal","s(sex):speciesSpotted seal", "s(age):speciesBeardedseal","s(age):speciesRinged seal","s(age):speciesSpotted seal","s(season):speciesBearded seal","s(season):speciesRinged seal","s(season):speciesSpotted seal")

crossval <- list()
predci <- list()
for(j in 1:nrow(bic)){
  
  modName <- rownames(bic)[j]
  mod <- get(modName)
  
  qt <- qt(0.975,df=mod$df.null)
  
  pred <- predict(mod,type="response",se.fit=TRUE,exclude=exclude)
  plot(modData$mass_kg_during_mean[which(!is.na(modData$mass_kg_during_mean))],pred$fit,xlab="mass_kg_during_mean",ylab="prediction")
  predci[[modName]] <- data.frame(modData %>% 
                         filter(!is.na(modData$mass_kg_during_mean)) %>% 
                         select(target_identifier,adjusted_image_dt,species,season,total_length,mass_kg_during_mean) %>%
                         rename(ID=target_identifier,date=adjusted_image_dt,mass=mass_kg_during_mean),
                       pred=pred$fit,
                       se=pred$se.fit,
                       lci=pred$fit - qt*pred$se.fit,
                       uci=pred$fit + qt*pred$se.fit)
  lines(modData$mass_kg_during_mean[which(!is.na(modData$mass_kg_during_mean))],modData$mass_kg_during_mean[which(!is.na(modData$mass_kg_during_mean))])
  predci[[modName]]$cv <- predci[[modName]]$se/predci[[modName]]$pred * 100
  predci[[modName]]$cover <- (predci[[modName]]$lci <= predci[[modName]]$mass & predci[[modName]]$uci >= predci[[modName]]$mass)
  predci[[modName]]$bias <- predci[[modName]]$pred - predci[[modName]]$mass
  predci[[modName]]$mse <- (predci[[modName]]$mass - predci[[modName]]$pred)^2
  predci[[modName]] <- predci[[modName]][order(predci[[modName]]$ID,predci[[modName]]$date),]
  plot(mod,residuals = TRUE,all.terms=TRUE,seWithMean=TRUE,shift=coef(mod)[1],pages=1)
  qq.gam(mod)
  
  mean(predci[[modName]]$cover)
  summary(predci[[modName]]$cv)
  summary(predci[[modName]]$bias/predci[[modName]]$mass*100)
  
  mean((predci[[modName]] %>% filter(species!="Bearded seal"))$cover)
  summary((predci[[modName]] %>% filter(species!="Bearded seal"))$cv)
  summary((predci[[modName]] %>% filter(species!="Bearded seal"))$bias/(predci[[modName]] %>% filter(species!="Bearded seal"))$obs*100)
  
  IDs <- unique((modData %>% filter(!is.na(mass_kg_during_mean)))$target_identifier)
  crossval[[modName]] <- vector('list',length(IDs))
  names(crossval[[modName]]) <- IDs
  for(i in IDs){
    fitsub <- tryCatch(gam(mod$formula,family=gaussian(link="log"),data=modData,subset=which(modData$target_identifier!=i),drop.unused.levels = FALSE,method="ML"),warning=function(w) w)
    if(inherits(fitsub,"warning")){
      if(grepl("gam.fit3 algorithm did not converge",fitsub)){
        fitsub <- tryCatch(gam(mod$formula,family=gaussian(link="log"),data=modData,subset=which(modData$target_identifier!=i),drop.unused.levels = FALSE),warning=function(w) w)
        if(inherits(fitsub,"warning")) stop(j," failed to converge")
      }
    }
    qtsub <- qt(0.975,df=fitsub$df.null)
    predsub <- tryCatch(predict(fitsub,newdata=modData[which(modData$target_identifier==i),],type="response",se.fit=TRUE,exclude=exclude),error=function(e) e)
    predcisub <- data.frame(modData %>% 
                              filter(modData$target_identifier==i) %>% 
                              select(target_identifier,measured_by,adjusted_range_m,species,sex,season,total_length,mass_kg_during_mean) %>% 
                              rename(ID=target_identifier,obs=measured_by,mass=mass_kg_during_mean),
                            pred=predsub$fit,
                            se=predsub$se.fit,
                            lci=predsub$fit - qtsub*predsub$se.fit,
                            uci=predsub$fit + qtsub*predsub$se.fit)
    predcisub$cv <- predcisub$se/predcisub$pred * 100
    predcisub$cover <- (predcisub$lci <= predcisub$mass & predcisub$uci >= predcisub$mass)
    predcisub$bias <- predcisub$pred-predcisub$mass
    predcisub$mse <- (predcisub$mass-predcisub$pred)^2
    crossval[[modName]][[i]] <- list(fit=fitsub,pred=predsub,predci=predcisub)
  }
  
  mean(unlist(lapply(crossval[[modName]][1:length(IDs)],function(x) x$predci$cover)))
  summary(unlist(lapply(crossval[[modName]][1:length(IDs)],function(x) x$predci$cv)))
  summary(unlist(lapply(crossval[[modName]][1:length(IDs)],function(x) x$predci$bias/x$predci$mass*100)))
  summary(unlist(lapply(crossval[[modName]][1:length(IDs)],function(x) x$predci$mse)))
  
  mean(unlist(lapply(crossval[[modName]][which(!(IDs %in% c("Noatak","Tuliq")))],function(x) x$predci$cover)))
  summary(unlist(lapply(crossval[[modName]][which(!(IDs %in% c("Noatak","Tuliq")))],function(x) x$predci$cv)))
  summary(unlist(lapply(crossval[[modName]][which(!(IDs %in% c("Noatak","Tuliq")))],function(x) x$predci$bias/x$predci$mass*100)))
  summary(unlist(lapply(crossval[[modName]][which(!(IDs %in% c("Noatak","Tuliq")))],function(x) x$predci$mse)))

}

order(unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$cover))))),decreasing=TRUE)
order(unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$cv))))),decreasing=FALSE)
order(unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) abs(y$predci$bias/y$predci$mass*100)))))),decreasing=FALSE)
order(unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$mse))))),decreasing=FALSE)

tab <- cbind(BIC=bic[,2],AICc=aic[,2],Cov=unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$cover))))),cv=unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$cv))))),PRB=unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$bias/y$predci$mass*100))))),mse=unlist(lapply(crossval,function(x) mean(unlist(lapply(x,function(y) y$predci$mse))))))
rownames(tab) <- gsub("X"," * ",gsub("."," + ",rownames(tab),fixed=TRUE),fixed=TRUE)
tab

overallCov <- c(mean=mean(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$cover))),se=sqrt(mean(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$cover)))*(1-mean(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$cover))))/mod$df.null),NA,NA,NA) # prediction coverage
overallCV <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$cv))) # prediction cv
overallPRB <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$bias/x$predci$mass*100))) # prediction PRB
overallMSE <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[1:length(IDs)],function(x) x$predci$mse))) # prediction MSE

# bearded
BeardedCov <- c(mean=mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$cover))),
  se=sqrt(mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$cover)))*(1-mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$cover))))/mod$df.null),NA,NA,NA)# coverage  
BeardedCV <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$cv))) # cv  
BeardedPRB <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$bias/x$predci$mass*100))) # PRB  
BeardedMSE <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Bearded seal")])))],function(x) x$predci$mse))) # MSE  

# ringed
RingedCov <- c(mean=mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$cover))),
                se=sqrt(mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$cover)))*(1-mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$cover))))/mod$df.null),NA,NA,NA)# coverage  
RingedCV <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$cv))) # cv  
RingedPRB <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$bias/x$predci$mass*100))) # PRB  
RingedMSE <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Ringed seal")])))],function(x) x$predci$mse))) # MSE  

# spotted
SpottedCov <- c(mean=mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$cover))),
                se=sqrt(mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$cover)))*(1-mean(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$cover))))/mod$df.null),NA,NA,NA)# coverage  
SpottedCV <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$cv))) # cv  
SpottedPRB <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$bias/x$predci$mass*100))) # PRB  
SpottedMSE <- custsum(unlist(lapply(crossval$length_sXspec.width_sXspec[which((IDs %in% unique(modData$target_identifier[which(modData$species=="Spotted seal")])))],function(x) x$predci$mse))) # MSE  

rbind(rbind(overallCov,overallCV,overallPRB,overallMSE),
      rbind(BeardedCov,BeardedCV,BeardedPRB,BeardedMSE),
      rbind(RingedCov,RingedCV,RingedPRB,RingedMSE),
      rbind(SpottedCov,SpottedCV,SpottedPRB,SpottedMSE))

dfpred <- predci$length_sXspec.width_sXspec
dfpred <- dfpred[order(dfpred$mass),]
pdf("results/Pred.pdf",width=8,height=8)
ggplot(dfpred,aes(x=mass,y=pred)) + geom_point(aes(col=species ), size=2) + xlab("true mass (kg)") + ylab("predicted mass (kg)") + geom_abline(intercept=0,slope=1,col=1) + geom_errorbar(aes(ymin = lci, ymax = uci, col=species), alpha=0.4, width = 1, size=1)
dev.off()
png("results/Pred.png",width=8,height=8,units="in",res=600)
ggplot(dfpred,aes(x=mass,y=pred)) + geom_point(aes(col=species ), size=2) + xlab("true mass (kg)") + ylab("predicted mass (kg)") + geom_abline(intercept=0,slope=1,col=1) + geom_errorbar(aes(ymin = lci, ymax = uci, col=species), alpha=0.4, width = 1, size=1)
dev.off()

dfcv <- do.call(rbind,lapply(crossval$length_sXspec.width_sXspec,function(x) x$predci))
dfcv <- dfcv[order(dfcv$mass),]
pdf("results/cvPred.pdf",width=8,height=8)
ggplot(dfcv,aes(x=mass,y=pred)) + geom_point(aes(col=species ), size=3) + xlab("true mass (kg)") + ylab("predicted mass (kg)") + geom_abline(intercept=0,slope=1,col=1) + geom_errorbar(aes(ymin = lci, ymax = uci, col=species), alpha=0.4, width = 2, size=1)
dev.off()
png("results/cvPred.png",width=8,height=8,units="in",res=600)
ggplot(dfcv,aes(x=mass,y=pred)) + geom_point(aes(col=species ), size=3) + xlab("true mass (kg)") + ylab("predicted mass (kg)") + geom_abline(intercept=0,slope=1,col=1) + geom_errorbar(aes(ymin = lci, ymax = uci, col=species), alpha=0.4, width = 2, size=1)
dev.off()

save.image("results/bodyCondition_analysis.RData")