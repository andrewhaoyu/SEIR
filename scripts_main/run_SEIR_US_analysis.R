## IMPORTANT: set code_root properly!
# code_root="/home/dgwu/covid19/NatSEIR_Rcode/"
#code_root="~/jianguoyun/Nutstore/covid19/NatSEIR_Rcode/"
#code_root="/Users/zhangh24/Desktop/codes_covid19_v2/"
# code_root="C:/Users/xingj/Documents/WangLabAdmin/COVID-19/NatSEIR_Rcode/"
#load state number
args = commandArgs(trailingOnly = T)
ind = as.numeric(args[[1]])
#number of statess
#number of replicates
n_rep <- 3

i1_opt = c(1:4)
i2_opt = c(1:3)
i3_opt = c(1:2)
total = length(i1_opt)*length(i2_opt)*length(i3_opt)
i1_vec = i2_vec = i3_vec = rep(0,total)
temp = 1
for(i1 in i1_opt){
  for(i2 in i2_opt){
    for(i3 in i3_opt){
      i1_vec[temp] = i1_opt[i1]
      i2_vec[temp] = i2_opt[i2]
      i3_vec[temp] = i3_opt[i3]
      temp = temp+1
    }
  }
  
}
i1 = i1_vec[ind]
i2 = i2_vec[ind]
i3 = i3_vec[ind]
method_vec = c("possion","nb")
method = method_vec[i3]
#code_root = "/data/zhangh24/SEIR/"
code_root = "/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"
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
source(paste0(code_root, "R/fun_SEIRfitting.R"))
source(paste0(code_root, "R/init_cond_update.R"))
source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/correlationPlot_modified.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##

statename = c("New York","Massachusetts",
              "Florida","Michigan")
allData <- read.csv("../data/US_State_data.csv")
idx <- which(allData$stateName==statename[i1])
N = allData$population[idx[1]]
print(statename[i1])
stateData <- allData[idx,]
#find first date with positive cases more than 20
jdx <- which(stateData$positiveIncrease>50)
#start analysis date
jan1_idx = min(jdx)

stateDataClean = stateData[jan1_idx:nrow(stateData),]
all.date <- as.Date(stateDataClean$date)
#leave 10 days for prediction
n.days <- nrow(stateDataClean)-10
n.days.all <- nrow(stateDataClean)
days_to_fit <- 1:n.days
#install.packages("lubridate")
library(lubridate)
date_in_model <- as.Date(stateDataClean$date)
start.date <- as.Date(stateDataClean$date[1])
end.date <- as.Date(stateDataClean$date[n.days])
all.cut.date <- c(floor_date(seq(start.date,end.date,by="month"),unit="month")+14,
ceiling_date(seq(start.date, end.date, by = 'month'), unit = "month")-1)
all.cut.date <- all.cut.date[order(all.cut.date)]
#remove the first cut date if it's too close to start.date
if(as.numeric(all.cut.date[1]-start.date)<=7){
  all.cut.date <- all.cut.date[-1]
}
#remove the end cut date if it's too close to end.date
if(as.numeric(end.date-all.cut.date[length(all.cut.date)]<=8)){
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
Di = 2.9
Dp = 4
De = 2.9

alpha = 0.55
Dh = 30
#force Apr 1, Apr 15, 
N <- stateDataClean$population[1]

Dq <- rep(0,n.stage)

GenerateDq <- function(cut.date){
  if(cut.date<="2020-04-01"){
    return(10)
  }else if(cut.date<="2020-04-15"){
    return(6)
  } else{
    return(3)
  }
}

for(i in 1:(n.stage-1)){
  Dq[i] <- GenerateDq(all.cut.date[i])
}
Dq[length(Dq)] = 3
init_sets_list=get_init_sets_list(r0 = 0.23,
                                  N = N,
                                  Dq = Dq,
                                  flowN = flowN,
                                  jan1_idx = jan1_idx,
                                  stateDataClean = stateDataClean,
                                  stateInput = statename[i1],
                                  stage_intervals = stage_intervals,
                                  stateData=stateData,
                                  method = method)

# good initial conditions
# c(1.284, 0.384, 0.174, 0.096, 0.161, -0.046, -0.379, 0.569)

SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = paste0(i1,"_",i2,"_",i3), output_ret = T, skip_MCMC=F,
            all.date = all.date,
            n_burn_in=40000,
            n_iterations=500000,
            method = method)

## to evaluate convergence, we run another two rounds of this program
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)
