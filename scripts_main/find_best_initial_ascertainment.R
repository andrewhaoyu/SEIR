#find the best initial ascertainment through best likelihood
loglh_func <- function(pars){
  
  
  ypred <- SEIRpred(pars, init_settings = init_sets_list)
  onset_obs <- init_sets_list$daily_new_case
  ypred <- ypred[, "Onset_expect"]
  phi = pars[length(pars)]
  #p = phi/(phi+as.numeric(ypred))
  # meant to suppress warnings when ypred is negative
  suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
                                #size = phi,
                                size = 1/phi,
                                mu = ypred,log=T))
  
  #dnbinom(x = as.numeric(onset_obs), 
  #          size = phi,
  #         mu = ypred,log=T))
  
  #if(any(p == 0) || any(is.nan(p))){
  if(any(is.nan(p))){  
    logL <- -Inf
  }else{
    logL <- sum(p)
  }
  return(logL)
}

statename = c("NY","MA",
              "FL","MI")
#downloaded from https://covidtracking.com/data/download
allData <- read.csv("../data/all-states-history.csv")
#keep date to 08/31/2020
library(lubridate)
date_in_model <- as.Date(allData$date,format="%m/%d/%Y")
idx <- which(date_in_model<="20-08-31")
allData <- allData[idx,]
#population number (downloaded from https://www.census.gov/data/datasets/time-series/demo/popest/2010s-state-total.html)
stateName = c("New York","Massachusetts",
              "Florida","Michigan")


library(readr)
library(coda)
library(cairoDevice)
file.dir <- "/data/zhangh24/SEIR/output/"
setwd(file.dir)

files <- dir(path=file.dir,pattern = "pars_est_run_")
library(data.table)
library(reshape2)
library(ggplot2)
logLikelihood_mat <- matrix(0,4,9)
for(i1 in 1:4){
  for(i3 in 1:9){
    
    result.list <- list()
    temp = 1
    
    for(i2 in 1:3){
      file <- paste0(paste0("pars_est_run_",i1,"_",i2,"_",i3,".txt"))
      if(file%in%files){
        pars_estimate_main = read.table(file=file,header=T)
        result.list[[temp]] <- pars_estimate_main
        temp = temp+1
      }
      
      
    }
    result <- rbindlist(result.list)
    pars <- colMeans(result)
    
    #plug in the population number
    population <- read.csv("../data/state_population.csv")
    idx <- which(population$State==stateName[i1])
    N = population$Population[idx]
    
    #idx <- which(allData$stateName==statename[i1])
    idx <- which(allData$state==statename[i1])
    
    print(statename[i1])
    stateData <- allData[idx,]
    #order the data by date
    stateData$date = as.Date(stateData$date,format="%m/%d/%Y")
    stateData = stateData[order(stateData$date),]
    
    #find first date with positive cases more than 20
    jdx <- which(stateData$positiveIncrease>50)
    #start analysis date
    jan1_idx = min(jdx)
    
    stateDataClean = stateData[jan1_idx:nrow(stateData),]
    all.date <- stateDataClean$date
    #leave 10 days for prediction
    n.days <- nrow(stateDataClean)-10
    n.days.all <- nrow(stateDataClean)
    days_to_fit <- 1:n.days
    #install.packages("lubridate")
    
    date_in_model <- stateDataClean$date
    start.date <- date_in_model[1]
    end.date <- date_in_model[n.days]
    all.cut.date <- c(floor_date(seq(start.date,end.date,by="month"),unit="month")+14,
                      ceiling_date(seq(start.date, end.date, by = 'month'), unit = "month")-1)
    all.cut.date <- all.cut.date[order(all.cut.date)]
    #remove the first cut date smaller than the start date
    if(as.numeric(all.cut.date[1]<=start.date)){
      all.cut.date <- all.cut.date[-1]
    }
    
    #remove the first cut date if it's too close to start.date
    if(as.numeric(all.cut.date[1]-start.date)<=7){
      all.cut.date <- all.cut.date[-1]
    }
    
    #remove the last cut date bigger than the end date
    if(end.date<=all.cut.date[length(all.cut.date)]){
      all.cut.date <- all.cut.date[-length(all.cut.date)]
    }
    
    #remove the end cut date if it's too close to end.date
    if(as.numeric(end.date-all.cut.date[length(all.cut.date)]<=7)){
      all.cut.date <- all.cut.date[-length(all.cut.date)]
    }
    idx <- which(date_in_model%in%all.cut.date)
    days.to.fit <- 1:length(date_in_model)
    n.stage <- length(idx)+1
    stage_intervals <- list()
    for(l in 1:(n.stage)){
      if(l<=(n.stage-1)&l>=2){
        stage_intervals[[l]] = c(start=idx[l-1]+1,end = idx[l])  
      }else if(l==1){
        stage_intervals[[l]] = c(start=1,end = idx[l])  
      }else{
        stage_intervals[[l]] = c(start=idx[l-1]+1,end = n.days)  
      }
      
    }
    flowN <- rep(0,n.stage)
    Di = 3.5
    Dp = 2.75
    De = 2.45
    
    alpha = 0.55
    Dh = 30
    
    
    Dq <- rep(0,n.stage)
    
    GenerateDq <- function(cut.date){
      if(cut.date<="20-04-01"){
        return(10)
      }else if(cut.date<="20-04-15"){
        return(6)
      } else{
        return(3)
      }
    }
    
    for(i in 1:(n.stage-1)){
      Dq[i] <- GenerateDq(all.cut.date[i])
    }
    Dq[length(Dq)] = 3
    
    r0_vec = c(0.05,0.10,0.15,0.20,0.23,0.30,0.35,0.40,0.5)
    r0 = r0_vec[i3]
    init_sets_list=get_init_sets_list(r0=r0,
                                      Di = Di,
                                      Dp = Dp,
                                      De = De,
                                      Dq =Dq,
                                      alpha = alpha,
                                      Dh = Dh,
                                      N=N,
                                      flowN = flowN,
                                      jan1_idx = jan1_idx,
                                      stateDataClean = stateDataClean,
                                      stateInput = statename[i1],
                                      stage_intervals = stage_intervals,
                                      stateData=stateData,
                                      method = method)
    
    logLikelihood_mat[i1,i3] <- loglh_func(pars)
  }
}

which.max(logLikelihood_mat[1,])
which.max(logLikelihood_mat[2,])
which.max(logLikelihood_mat[3,])
which.max(logLikelihood_mat[4,])
