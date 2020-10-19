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
for(i1 in 1:8){
  statename = c("New York","Massachusetts",
                "Florida","Michigan",
                "Connecticut",
                "Illinois",
                "Indiana",
                "Louisiana")
  
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
    theme_Publication()+
    ggtitle(paste0(statename[i1]," daily reported cases by JHU"))
  png(filename = paste0("../output/result_101220/weekday_effect_",i1,".png"),
      width = 10, height = 8, res =300, units = "in")
  print(p)
  dev.off()
}


