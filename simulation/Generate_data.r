#Generate six stages analysis
args = commandArgs(trailingOnly = T)
#i1 replicates
i1 = as.numeric(args[[1]])
#i2 phi
i2 = as.numeric(args[[2]])
#i3 ascertainment
i3 = as.numeric(args[[3]])
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
source("estimate_R.R")
logit <- function(x){
  log(x/(1-x))
}

set.seed(i1)
n_stage = 2
n.stage = n_stage
initial.ascertainment = 0.20
N = 21477737
stage_intervals = list()
total= 0
for(i in 1:n_stage){
  stage_intervals[[i]] <- c(start = total+1,end = total+15)
  total = total+30
}
#b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
b_vec = c(0.7,0.3)
phi_vec = c(0,0.2)
ascertainment_mat = matrix(c(c(0.05,0.07),
                             (c(0.05,0.07)+0.2)),ncol=2)

# ascertainment_mat = matrix(c(c(0.05,0.07,0.09,0.10,0.11,0.12),
#                              (c(0.05,0.07,0.09,0.10,0.11,0.12)+0.2)),ncol=2)
phi = phi_vec[i2]
r_vec = ascertainment_mat[,i3]

Di = 3.5
Dp = 2.75
De = 2.45
#Dq_vec = c(10,6,3,3,3,3)
Dq_vec = c(10,6)
alpha = 0.55
Dh = 30
flowN_vec = c(0,0)
days_to_fit = c(1:total)
init_states <- c(21470907,2295,3215,264,1056,0,0)
names(init_states) <- c("S","E","P","I","A","H","R")
SEIR_mat= SEIRpred(stage_intervals,b_vec,r_vec,
         Di,Dp,
         De,Dq_vec,
         alpha,Dh,N,flowN_vec,init_states,
         days_to_fit)
 



onset_expect = SEIR_mat[,"Onset_expect"]

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

est_result = SEIRfitting(
  n_burn_in=n_burn_in,
  n_iterations=n_iterations,
  all.date = all.date,
  onset_obs,
  init_states,
  n_stage,
  par_lower,
  par_upper)

#save(est_result,file = paste0("/data/zhangh24/SEIR/result/simulation/seir_result_",i1,"_",i2,"_",i3,".rdata"))
save(est_result,file = paste0("/data/zhangh24/SEIR/result/simulation/two_stage_",i1,"_",i2,"_",i3,".rdata"))