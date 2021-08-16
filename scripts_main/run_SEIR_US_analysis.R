## IMPORTANT: set code_root properly!
# code_root="/home/dgwu/covid19/NatSEIR_Rcode/"
#code_root="~/jianguoyun/Nutstore/covid19/NatSEIR_Rcode/"
#code_root="/Users/zhangh24/Desktop/codes_covid19_v2/"
# code_root="C:/Users/xingj/Documents/WangLabAdmin/COVID-19/NatSEIR_Rcode/"
#load state number
args = commandArgs(trailingOnly = T)
#i1 represent state
#i2 represent replicate
#i3 represent initial ascertainment
#i4 represent prior
#i5 represent reparametrization
#i6 represent asympotomatic infection rate
i1 = as.numeric(args[[1]])
i2 = 1
i3 = as.numeric(args[[2]])
i4 = 1
i5 = 1
#i6 = as.numeric(args[[3]])
set.seed(i1*1000+i2*100+i3)
# ind = as.numeric(args[[1]])
# #number of statess
# #number of replicates
# n_rep <- 3
# 
# i1_opt = c(1:4)
# i2_opt = c(1:3)
# i3_opt = c(1:2)
# total = length(i1_opt)*length(i2_opt)*length(i3_opt)
# i1_vec = i2_vec = i3_vec = rep(0,total)
# temp = 1
# for(i1 in i1_opt){
#   for(i2 in i2_opt){
#     for(i3 in i3_opt){
#       i1_vec[temp] = i1_opt[i1]
#       i2_vec[temp] = i2_opt[i2]
#       i3_vec[temp] = i3_opt[i3]
#       temp = temp+1
#     }
#   }
#   
# }
# i1 = i1_vec[ind]
# i2 = i2_vec[ind]
#i3 = i3_vec[ind]
method_vec = c("poisson","nb")
method = method_vec[2]
code_root = "/data/zhangh24/SEIR/"
#code_root = "/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"
#code_root = "/dcl01/chatterj/data/hzhang1/temp/SEIR/"
setwd(paste0(code_root, "scripts_main"))
#install.packages("BayesianTools")
library(BayesianTools)
#install.packages("vioplot")
#library(vioplot)
#install.packages("corrplot")
library(corrplot)
library(readr)
#install.packages("cairoDevice")
#library(cairoDevice)
library(dplyr)
##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu_update.R"))

#source(paste0(code_root, "R/init_cond_update.R"))
#if(i5 ==1){
  source(paste0(code_root, "R/fun_SEIRfitting.R"))  
  source(paste0(code_root, "R/init_cond_update.R"))
# }else if(i5 ==2){
#   source(paste0(code_root, "R/fun_SEIRfitting_new.R"))  
#   source(paste0(code_root, "R/init_cond_new.R"))
# }

source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/correlationPlot_modified.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##


#use covidtracing data to analyze
#downloaded from https://covidtracking.com/data/download
statename = c("FL",
              "MA")

# statename = c("NY",
#               "MA",
#               "FL",
#               "MI",
#               "CT",
#               "LA",
#               "MO","UT","IN","OH")
#downloaded from https://covidtracking.com/data/download
allData <- read.csv("../data/all-states-history.csv")
#keep date to 08/31/2020
library(lubridate)
#leave off days for prediction
leave_days = 11

date_in_model <- as.Date(allData$date,format="%Y-%m-%d")
# idx <- which(date_in_model<="2021-1-11"&
#                date_in_model>="2020-11-01")
#allData <- allData[idx,]
#population number (downloaded from https://www.census.gov/data/datasets/time-series/demo/popest/2010s-state-total.html)
stateName = c("Florida",
              "Massachusetts")

# stateName = c("New York","Massachusetts",
#               "Florida","Michigan",
#               "Connecticut",
#               "Louisiana",
#               "Missouri",
#               "Utah","Indiana","Ohio")
#
#plug in the population number
population <- read.csv("../data/state_population.csv")
idx <- which(population$State==stateName[i1])
N = population$Population[idx]
#

idx <- which(allData$state==statename[i1])

print(statename[i1])
stateData <- allData[idx,]
#order the data by date
stateData$date = as.Date(stateData$date,format="%Y-%m-%d")
stateData = stateData[order(stateData$date),]


cbind(weekdays(as.Date(stateData$date))
  ,
      stateData$positiveIncrease)
#use JHU data to analyze
#download data from https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_State.csv
# statename = c("New York","Massachusetts",
#              "Florida","Michigan")
# 
# allData <- read.csv("../data/US_State_data.csv")
# 
# library(lubridate)
# 
# date_in_model <- as.Date(allData$date,format="%Y-%m-%d")
# allData$date <- date_in_model
# idx <- which(date_in_model<="2020-08-31")
# 
# allData <- allData[idx,]
# idx <- which(allData$stateName==statename[i1])
# N = allData$population[idx][1]
# idx <- which(allData$stateName==statename[i1])
# stateData <- allData[idx,]

#find first date with positive cases more than 20
jdx <- which(stateData$positiveIncrease>50)
#start analysis date
jan1_idx = min(jdx)

stateDataClean = stateData[jan1_idx:nrow(stateData),]

#clean the outlier data for MA and CT
#MA added antibody tests results into the data
#the data suddenly increased a lot
#to avoid 
# if(i1 ==2){
#   idx <- which(stateDataClean$date=="2020-06-01")
#   stateDataClean$positiveIncrease[idx]  =  as.integer((stateDataClean$positiveIncrease[idx-1]+stateDataClean$positiveIncrease[idx+1])/2)
# }
# #CT data is usually 0 during weekend
# #it's due to no reporting during weekend
# #to aviod this, we dropped weekend CT data from the loglikelihood after 2020-07-04
# if(i1==5){
#   idx <- which(stateDataClean$date>="2020-07-04"&
#                  (weekdays(stateDataClean$date)=="Saturday"|
#                     weekdays(stateDataClean$date)=="Sunday"))
#   subset.id = which(c(1:(nrow(stateDataClean)-leave_days))%in%idx==F)
# }
all.date <- stateDataClean$date
#leave 10 days for prediction
n.days <- nrow(stateDataClean)-leave_days
n.days.all <- nrow(stateDataClean)
days_to_fit <- 1:n.days
#install.packages("lubridate")

date_in_model <- stateDataClean$date
start.date <- date_in_model[1]
end.date <- date_in_model[n.days]
# all.cut.date <- c(floor_date(seq(start.date,end.date,by="month"),unit="month")+14,
# ceiling_date(seq(start.date, end.date, by = 'month'), unit = "month")-1)
# all.cut.date <- all.cut.date[order(all.cut.date)]
all.cut.date <- c(ceiling_date(seq(start.date, end.date, by = 'month'), unit = "month")-1)
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
#add addtional cut for NY data due to dramatic change
# if(i1==1){
#   all.cut.date<- c(as.Date("2020-03-20"),all.cut.date)
# }
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



alpha <- 0.55

Dh = 30




Dq <- rep(0,n.stage)

GenerateDq <- function(cut.date){
  if(cut.date<="2020-04-01"){
    return(10)
  }else if(cut.date<="2020-04-30"){
    return(6)
  } else{
    return(3)
  }
}

for(i in 1:(n.stage-1)){
  Dq[i] <- GenerateDq(all.cut.date[i])
}
Dq[length(Dq)] = 3
library(dplyr)

r0_vec = c(0.05,0.075,0.10,0.125,0.15,0.20,0.225,0.25,0.275,0.30)+0.15
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
#if(i4==1){
  beta_shape1 <- 1
  beta_shape2 <- 1
  
# }else if(i4==2){
#   beta_shape1 <- 7.3
#   beta_shape2 <- 24.6
#   
# }


library(invgamma)
#update the outlier
idx <- which(init_sets_list$daily_new_case<0)
init_sets_list$daily_new_case[idx]= 0
idx <- which(init_sets_list$daily_new_case_all<0)
init_sets_list$daily_new_case_all[idx]= 0
delta_mean <- 0
delta_sd <- 0.5
#beta_shape1 <- 7.3
#beta_shape2 <- 24.6
beta_shape1 <- 1
beta_shape2 <- 1
gamma_shape = 3
gamma_rate = 1


init_sets_list$init_states = c(18875899,34200,35931,10451,21189,78039,2422035)
names(init_sets_list$init_states) = c("S","E","P","I","A","H","R")





SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = paste0("080221_",i1,"_",i2,"_",i3), output_ret = T, skip_MCMC=F,
            all.date = all.date,
           # n_burn_in=2800,
           # n_iterations=30000,
            n_burn_in=50000,
            n_iterations=1000000,
            method = method)
#run_id = paste0("040521_",i1,"_",i2,"_",i3)
#quantile(mcmc_pars_estimate_original[,"r2"],probs = c(0.025,0.975))
## to evaluate convergence, we run another two rounds of this program
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)
