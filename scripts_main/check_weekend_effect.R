i1 = as.numeric(args[[1]])
code_root = "/Users/zhangh24/GoogleDrive/covid_project/SEIR/"
#code_root = "/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"
#code_root = "/dcl01/chatterj/data/hzhang1/temp/SEIR/"
setwd(paste0(code_root, "scripts_main"))
#install.packages("BayesianTools")
library(BayesianTools)
#install.packages("vioplot")
library(vioplot)
#install.packages("corrplot")
library(corrplot)
library(readr)
#install.packages("cairoDevice")
#library(cairoDevice)
library(dplyr)
##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu.R"))

#source(paste0(code_root, "R/init_cond_update.R"))

  source(paste0(code_root, "R/fun_SEIRfitting.R"))  
  source(paste0(code_root, "R/init_cond_update.R"))

source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/correlationPlot_modified.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##


#use covidtracing data to analyze
#downloaded from https://covidtracking.com/data/download
# statename = c("NY","MA",
#                             "FL","MI")
# #
# allData <- read.csv("../data/all-states-history.csv")
# #keep date to 08/31/2020
# library(lubridate)
# date_in_model <- as.Date(allData$date,format="%m/%d/%Y")
# idx <- which(date_in_model<="20-08-31")
# allData <- allData[idx,]
# #population number (downloaded from https://www.census.gov/data/datasets/time-series/demo/popest/2010s-state-total.html)
# stateName = c("New York","Massachusetts",
#               "Florida","Michigan")
# #
# #plug in the population number
# population <- read.csv("../data/state_population.csv")
# idx <- which(population$State==stateName[i1])
# N = population$Population[idx]
# #
# 
# idx <- which(allData$state==statename[i1])
# 
# print(statename[i1])
# stateData <- allData[idx,]
# #order the data by date
# stateData$date = as.Date(stateData$date,format="%m/%d/%Y")
# stateData = stateData[order(stateData$date),]



#use JHU data to analyze
#download data from https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_State.csv
for(i1 in 1:4){
  statename = c("New York","Massachusetts",
                "Florida","Michigan")
  
  allData <- read.csv("../data/US_State_data.csv")
  
  library(lubridate)
  
  date_in_model <- as.Date(allData$date,format="%Y-%m-%d")
  allData$date <- date_in_model
  idx <- which(date_in_model<="2020-08-31")
  
  allData <- allData[idx,]
  idx <- which(allData$stateName==statename[i1])
  N = allData$population[idx][1]
  idx <- which(allData$stateName==statename[i1])
  stateData <- allData[idx,]
  
  #find first date with positive cases more than 20
  jdx <- which(stateData$positiveIncrease>50)
  #start analysis date
  jan1_idx = min(jdx)
  
  stateDataClean = stateData[jan1_idx:nrow(stateData),]
  all.date <- stateDataClean$date
  weekday <- weekdays(all.date)
  
  weekend <- weekday
  idx <- which(weekday=="Saturday"|weekday=="Sunday")
  weekend[idx] <- "weekend"
  idx <- which(weekday!="Saturday"&weekday!="Sunday")
  weekend[idx] <- "weekday"
  stateDataClean$weekend <- weekend
  
  library(ggplot2)
  
  p <- ggplot(stateDataClean) + 
    geom_point(aes(date,positiveIncrease,color=weekend))+
    theme_Publication()
  png(filename = paste0("../output/result_101220/weekday_effect_",i1,".png"),
      width = 10, height = 8, res =300, units = "in")
  print(p)
  dev.off()
}


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



alpha_vec <- c(0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85)
alpha = alpha_vec[i6]
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
r0_vec = c(0.10,0.15,0.20,0.23,0.30,0.35,0.40,0.5)
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

# good initial conditions
# c(1.284, 0.384, 0.174, 0.096, 0.161, -0.046, -0.379, 0.569)
if(i4==1){
  beta_shape1 <- 1
  beta_shape2 <- 1
  
}else if(i4==2){
  beta_shape1 <- 7.3
  beta_shape2 <- 24.6
  
}


library(invgamma)
SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = paste0("101220_",i1,"_",i2,"_",i6), output_ret = T, skip_MCMC=F,
            all.date = all.date,
            n_burn_in=170000,
            n_iterations=2000000,
            method = method)

## to evaluate convergence, we run another two rounds of this program
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)
