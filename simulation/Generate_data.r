#Generate six stages analysis
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(BayesianTools)
#install.packages("vioplot")
#library(vioplot)
#install.packages("corrplot")
library(corrplot)
library(readr)
#install.packages("cairoDevice")
#library(cairoDevice)
library(dplyr)
library(invgamma)
setwd("/data/zhangh24/SEIR/simulation")
source("fun_SEIRfitting.R")
source("fun_SEIRpred.R")
source("fun_SEIRsimu.R")
set.seed(i1)
n_stage = 6
n.stage = n_stage
initial.ascertainment = 0.20
N = 21477737
stage_intervals = list()
total= 0
for(i in 1:n_stage){
  stage_intervals[[i]] <- c(start = total+1,end = total+15)
  total = total+15
}
b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
r_vec = c(0.2,0.25,0.3,0.35,0.4,0.30)
Di = 3.5
Dp = 2.75
De = 2.45
Dq_vec = c(10,6,3,3,3,3)
alpha = 0.55
Dh = 30
flowN_vec = c(0,0,0,0,0,0)
days_to_fit = c(1:total)
init_states <- c(21470907,2295,3215,264,1056,0,0)
names(init_states) <- c("S","E","P","I","A","H","R")
SEIR_mat =SEIRsimu(stage_intervals,b_vec,r_vec,
         Di,Dp,
         De,Dq_vec,
         alpha,Dh,N,flowN_vec,init_states,
         days_to_fit)





transform_var_fun=function(pars) {
  n.stage <- length(pars)/2
  b_vec <- pars[1:n.stage]
  r_vec <- pars[(n.stage+1):(2*n.stage)]
  r1 = r_vec[1]
  for(l in 2:n.stage){
    rtemp = 1 / (1 + (1 - r_vec[l-1]) / (r_vec[l-1] * exp(r_vec[l])))
    r_vec[l] = rtemp
  }
  
  return(list(b_vec, r_vec))
}


onset_obs = SEIR_mat[,"Onset_expect"]
par_lower = c(0,rep(-10,n_stage-1),0,rep(-10,n_stage-1),0)
par_upper = c(3,rep(10,n_stage-1),1,rep(10,n_stage-1),1000)
n_iterations = 800000
n_burn_in = 72000
delta_mean <- 0
delta_sd <- 1
#beta_shape1 <- 7.3
#beta_shape2 <- 24.6
beta_shape1 <- 1
beta_shape2 <- 1
gamma_shape = 1
gamma_rate = 1

est_result = SEIRfitting(
  n_burn_in=n_burn_in,
  n_iterations=n_iterations,
  all.date = all.date,
  onset_obs,
  init_states,
  n_stage,
  par_lower,
  par_upper)

save(est_result,file = paste0("/data/zhangh24/SEIR/result/simulation/seir_result_",i1,".rdata"))
