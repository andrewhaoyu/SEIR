mcmc_pars_estimate  = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header=T)
startValue = colMeans(mcmc_pars_estimate)
print(startValue)
loglh_func(startValue)


n.stage = 10
delta = rep(0,n_stage-1)
load("../output/r_vec_1.rdata")
logit <- function(x){log(x/(1-x))}
for(k in 1:(n.stage-1)){
  delta[k] = logit(r_vec[k+1])-logit(r_vec[k])
  
}










n.stage = 10
c1 = 0.38585
c0 = -0.52274
startValue[n.stage+1] = -0.52274
startValue[n.stage+2] = 0.38585
loglh_func(startValue)

save(r_vec,file = "../output/r_vec_1.rdata")
mcmc_pars_estimate  = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header=T)
startValue = colMeans(mcmc_pars_estimate)
c1 = startValue[n.stage+2]
c0 = startValue[n.stage+1]
r_vec = logit_inv(c0+c1*test_stage)
save(r_vec,file = "../output/r_vec_2.rdata")
