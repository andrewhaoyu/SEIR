args = commandArgs(trailingOnly = T)
i3 = as.numeric(args[[1]])
#check the ascertainment rate for different states
#code_root="/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"
code_root="/Users/zhangh24/GoogleDrive/covid_project/SEIR/"
setwd(paste0(code_root, "scripts_main"))
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
  library(lubridate)
  
  date_in_model <- as.Date(allData$date,format="%Y-%m-%d")
  allData$date <- date_in_model
  idx <- which(date_in_model<="2020-08-31")
  
  allData <- allData[idx,]
  
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
  n.stage <- length(idx)
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


transform_var_main_stage=function(pars) {
  n.stage <- (length(pars)-1)/2
  #b_vec <- pars[1:n.stage]
  b_vec <- pars[1:n.stage]
  for(l in 2:n.stage){
    b_temp = b_vec[l-1]*exp(b_vec[l])
    b_vec[l] = b_temp
  }
  r_vec <- pars[(n.stage+1):(2*n.stage)]
  r1 = r_vec[1]
  for(l in 2:n.stage){
    rtemp = 1 / (1 + (1 - r_vec[l-1]) / (r_vec[l-1] * exp(r_vec[l])))
    r_vec[l] = rtemp
  }
  
  return(list(b_vec, r_vec))
}
#for(i3 in 1:2){
  analysis_result_table <- NULL
  for(i1 in 1:4){
    n.stage = stage.list[i1]
    pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)))
    mcmc_pars_estimate=read.table(paste0("../output/result_101220/output/pars_est_run_101220_",i1,"_",2,"_",5,".txt"), header=T)
    r_mat <- mcmc_pars_estimate[,(n.stage+1):(2*n.stage)]
    colnames(r_mat) <- paste0("r",c(1:n.stage))
    
    # for(l in 1:nrow(mcmc_pars_estimate)){
    #   r_mat[l,] <- transform_var_main_stage(mcmc_pars_estimate[l,])[[2]]
    # }
    # 
    r_est = transform_var_main_stage(colMeans(mcmc_pars_estimate))[[2]]
    # r_est_low = apply(r_mat,2,function(x){quantile(x,0.025)})
    # r_est_high = apply(r_mat,2,function(x){quantile(x,0.975)})
    # 
    analysis_result_table_temp <- r_est
    # pla = 2
    # for(l in 1:length(r_est)){
    #   analysis_result_table_temp[l] =  paste0(round(r_est[l],pla)," (",
    #                                      round(r_est_low[l],pla),
    #                                      "-",
    #                                      round(r_est_high[l],pla)
    #                                      ,")")
    # }
    # 
    analysis_result_table <- rbind(analysis_result_table,
                                   analysis_result_table_temp)
    
  }
  write.csv(analysis_result_table,file=
              paste0("../output/ascertainment_rate_101320.csv"))
  
#}
