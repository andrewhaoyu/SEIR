rm(list = ls())
## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory where the repository README file is located. 
## For more information, please read the repository README file
code_root="/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"

setwd(paste0(code_root, "scripts_main"))





#confirm number of stage
stage.list <- rep(0,4)
i2 = 1
for(i1 in 1:4){
  
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
  stage.list[i1] = n.stage
}







library(readr)
library(coda)
library(cairoDevice)

# if (!file.exists("../output/pars_est_run_main_analysis.txt") | 
#     !file.exists("../output/pars_est_run_main_analysis_rep1.txt") |
#     !file.exists("../output/pars_est_run_main_analysis_rep2.txt")) {
#   stop("Outputs from main analysis cannot be found.\n
#       Probably scripts_main/Run_SEIR_main_analysis.R has not yet been run or code_root is not set correctly.")
# }

#i2 = 1
for(i1 in 1:4){
  n.stage = stage.list[i1]
  pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)))
  pars_estimate_main=read.table(paste0("../output/pars_est_run_",i1,"_",1,".txt"), header=T)
  pars_estimate_main_rep1=read.table(paste0("../output/pars_est_run_",i1,"_",2,".txt"), header=T)
  pars_estimate_main_rep2=read.table(paste0("../output/pars_est_run_",i1,"_",3,".txt"), header=T)
  
  mcmc_main=mcmc(data=pars_estimate_main)
  mcmc_rep1=mcmc(data=pars_estimate_main_rep1)
  mcmc_rep2=mcmc(data=pars_estimate_main_rep2)
  
  mcmc_3traj=mcmc.list(mcmc_main,mcmc_rep1,mcmc_rep2)
  
  gelman.diag(mcmc_3traj)
  # our results: multivariate psrf = 1
  p_vec <- c(10000:20000)
  
  plot_par_3traj = function(par_name, plotmath_name) {
    # red-like: #BC3C29, blue-like: #0072B5, orange-like: #E18727
    plot(p_vec, pars_estimate_main[p_vec, par_name], type="l", col="#0072B5", main=plotmath_name, xlab="", ylab="")
    points(p_vec, pars_estimate_main_rep1[p_vec, par_name], type="l", col="#BC3C29")
    points(p_vec, pars_estimate_main_rep2[p_vec, par_name], type="l", col="#E18727")
  }
  
  cairo_pdf(paste0("../output/mcmc_convergence",i1,".pdf"), width=10, height=10)
  par(mfrow=c(3, ceiling(ncol(pars_estimate_main)/3)))
  for(k in 1:ncol(pars_estimate_main)){
    plot_par_3traj(colnames(pars_estimate_main)[k], colnames(pars_estimate_main)[k])
  }
  dev.off()
  
  
}
