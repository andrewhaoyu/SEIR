#Generate six stages analysis
args = commandArgs(trailingOnly = T)
#i1 replicates
i1 = as.numeric(args[[1]])
#days in each stage
#i2 = as.numeric(args[[2]])
#prior distribution
#i3 = as.numeric(args[[3]])
#isolation time
#i4 = as.numeric(args[[4]])

i2 = 2
i3 = 2
i4 = 1

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
source("fun_SEIRsimu_update.R")
source("estimate_R.R")
logit <- function(x){
  log(x/(1-x))
}

set.seed(i1)
n_stage = 6
n.stage = n_stage
initial.ascertainment = 0.20
N = 21477737
stage_intervals = list()
total= 0
days_in_each_stage_vec = c(15,30)
days_in_each_stage = days_in_each_stage_vec[i2]
for(i in 1:n_stage){
  stage_intervals[[i]] <- c(start = total+1,end = total+days_in_each_stage)
  total = total+days_in_each_stage
}
b_vec = c(0.81,0.22,0.19,0.28,0.31,0.48)
#b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
#b_vec = c(0.7,0.3)
#phi_vec = c(0,0.2)
# ascertainment_mat = matrix(c(c(0.05,0.07),
#                              (c(0.05,0.07)+0.2)),ncol=2)

# ascertainment_mat = matrix(c(c(0.05,0.07,0.09,0.10,0.11,0.12),
#                              (c(0.05,0.07,0.09,0.10,0.11,0.12)+0.2)),ncol=2)
phi = 0.2
r_vec = c(0.23,0.24,0.38,0.57,0.69,0.75)

Di = 3.5
Dp = 2.75
De = 2.45
#Dq_vec = c(10,6,3,3,3,3)

Dq_vec = c(10,6,3,3,3,3)
if(i4 ==2){
  Dq_vec = c(5,3)
}
alpha = 0.55
Dh = 30
flowN_vec = rep(0,n.stage)
days_to_fit = c(1:total)
init_states <- c(21470907,2295,3215,264,1056,0,0)
names(init_states) <- c("S","E","P","I","A","H","R")
SEIR_mat= SEIRpred(stage_intervals,b_vec,r_vec,
         Di,Dp,
         De,Dq_vec,
         alpha,Dh,N,flowN_vec,init_states,
         days_to_fit)
 



onset_expect = SEIR_mat[,"Onset_expect"]
S_mat = SEIR_mat[,"S"]
true_prevalence <- (N-S_mat)/N

onset_obs <- rnbinom(n = length(onset_expect),
                     size = 1/phi,
                     mu = onset_expect)
par_lower = c(0,rep(-10,n_stage-1),0,rep(-10,n_stage-1),0)
par_upper = c(3,rep(10,n_stage-1),1,rep(10,n_stage-1),1000)
n_iterations = 100000
n_burn_in = 9090
delta_mean <- 0
delta_sd <- 1
#beta_shape1 <- 7.3
#beta_shape2 <- 24.6
beta_shape1 <- 1
beta_shape2 <- 1
gamma_shape = 1
gamma_rate = 1

if(i3 ==2){
  delta_mean <- 0
  delta_sd <- 0.5
  #beta_shape1 <- 7.3
  #beta_shape2 <- 24.6
  beta_shape1 <- 1
  beta_shape2 <- 1
  gamma_shape = 3
  gamma_rate = 1
}

logit = function(x){
  log(x/(1-x))
}

transformed_b = b_vec
transformed_r = r_vec
for(l in 2:length(b_vec)){
  transformed_b[l] = log(b_vec[l])-log(b_vec[l-1])
  transformed_r[l] = logit(r_vec[l])-logit(r_vec[l-1])
}

# true_pars = c(transformed_b,transformed_r,phi)
# 
# loglh_func(true_pars)
# tmp_est = var_trans_fun(pars)
# b_vec = tmp_est[[1]]
# r_vec = tmp_est[[2]]
# 
# ypred <- SEIRpred(stage_intervals,b_vec,r_vec,
#                   Di,Dp,
#                   De,Dq_vec,
#                   alpha,Dh,N,flowN_vec,init_states,
#                   days_to_fit)
# onset_obs <- onset_obs
# ypred <- ypred[, "Onset_expect"]
# phi = pars[length(pars)]
# #p = phi/(phi+as.numeric(ypred))
# # meant to suppress warnings when ypred is negative
# 
# suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
#                               #size = phi,
#                               size = 1/phi,
#                               mu = ypred,log=T))
# 
# 



est_result = SEIRfitting(
  n_burn_in=n_burn_in,
  n_iterations=n_iterations,
  onset_obs,
  init_states,
  n_stage,
  par_lower,
  par_upper)

est_result$true_prevalence = true_prevalence

#save(est_result,file = paste0("/data/zhangh24/SEIR/result/simulation/seir_result_",i1,"_",i2,"_",i3,".rdata"))
save(est_result,file = paste0("/data/zhangh24/SEIR/result/simulation/fourstage_two_stage_",i1,"_",i2,"_",i3,"_",i4,".rdata"))