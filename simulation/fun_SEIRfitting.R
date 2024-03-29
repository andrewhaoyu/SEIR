default_pars_density <- function(pars) {
  n.stage = (length(pars)-1)/2
  d_vec <- rep(NA, length(pars))
  ##b12, b3, b4, b5
  #bvec
  # for(i in c(1:n.stage)) {
  #   d_vec[i] <- dunif(pars[i], 0, 2, log = T)
  # }
  d_vec[1] =  log(1/2)
  d_vec[2:n.stage] =  dnorm(pars[2:(n.stage)],delta_mean, delta_sd, log = T)
  #rvec1
  d_vec[n.stage+1] = dbeta(pars[n.stage+1],beta_shape1, beta_shape2, log = T)
  #rvec
  d_vec[(n.stage+2):(2*n.stage)] = 
    dnorm(pars[(n.stage+2):(2*n.stage)],delta_mean, delta_sd, log = T)
  phi = pars[2*n.stage+1]
  #d_vec[2*n.stage+1] = dgamma(phi,gamma_shape,gamma_rate,log =T)
  d_vec[2*n.stage+1] = dinvgamma(phi,shape = gamma_shape,
                                 rate = gamma_rate,log =T)
  return(sum(d_vec))
  #rvec
  # for(l in (n.stage+1):(2*n.stage)){
  #   r_temp = pars[l]
  #   d_vec[l] = dbeta(r_temp, beta_shape1, beta_shape2, log = T)
  # }
  
}

default_pars_sampler <- function(n.stage=n.stage) {
  s_vec <- matrix(NA, 1, 2*n.stage+1)
  
  
  ## b1
  s_vec[, 1] <- runif(1, 0, 2) 
  ## b2...
  s_vec[, 2:n.stage] <- rnorm(n.stage-1, delta_mean, delta_sd) 
  #r1
  s_vec[, n.stage+1] <- rbeta(1, beta_shape1, beta_shape2)
  #r2...
  s_vec[, (n.stage+2):(2*n.stage)] <- rnorm(n.stage-1, delta_mean, delta_sd)
  
  # ## b12, b3, b4, b5
  # s_vec[, 1:n.stage] <- runif(n.stage, 0, 2) 
  # s_vec[, n.stage+1] <- rbeta(1, beta_shape1, beta_shape2)
  # s_vec[, (n.stage+2):(2*n.stage)] <- rnorm(n.stage-1, delta_mean, delta_sd)
  s_vec[, 2*n.stage+1] <- rinvgamma(1,
                                    shape = gamma_shape,
                                    rate = gamma_rate)
  
  #rgamma(1,gamma_shape,gamma_rate)
  return(s_vec)
}

loglh_func <- function(pars){
  tmp_est = var_trans_fun(pars)
  b_vec = tmp_est[[1]]
  r_vec = tmp_est[[2]]
  
  ypred <- SEIRpred(stage_intervals,b_vec,r_vec,
                    Di,Dp,
                    De,Dq_vec,
                    alpha,Dh,N,flowN_vec,init_states,
                    days_to_fit)
  onset_obs <- onset_obs
  ypred <- ypred[, "Onset_expect"]
  phi = pars[length(pars)]
  #p = phi/(phi+as.numeric(ypred))
  # meant to suppress warnings when ypred is negative
  
  suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
                                #size = phi,
                                size = 1/phi,
                                mu = ypred,log=T))
  
  
  
  #dnbinom(x = as.numeric(onset_obs), 
  #          size = phi,
  #         mu = ypred,log=T))
  
  #if(any(p == 0) || any(is.nan(p))){
  if(any(is.nan(p))){  
    logL <- -Inf
  }else{
    logL <- sum(p)
  }
  return(logL)
}

var_trans_fun= function(pars) {
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


trans_ori_to_delta= function(pars) {
  n.stage <- (length(pars)-1)/2
  #b_vec <- pars[1:n.stage]
  delta_b_vec <- pars[1:n.stage]
  for(l in 2:n.stage){
    delta_b_temp = log(b_vec[l])-log(b_vec[l-1])
    delta_b_vec[l] = delta_b_temp
  }
  delta_r_vec <- pars[(n.stage+1):(2*n.stage)]
  logit = function(x){log(x/(1-x))}
  r1 = r_vec[1]
  for(l in 2:n.stage){
    deltartemp = logit(r_vec[l])-logit(r_vec[l-1])
    delta_r_vec[l] = deltartemp
  }
  phi = pars[length(pars)]
  return(c(delta_b_vec, delta_r_vec,phi))
}


SEIRfitting=function(
  n_burn_in=n_burn_in,
  n_iterations=n_iterations,
  onset_obs,
  init_states,
  n_stage,
  par_lower,
  par_upper) {
  
  n.stage = n_stage
  pars_density = default_pars_density
  pars_sampler = default_pars_sampler
  loglh_func = loglh_func
  pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)),paste0("phi"))
  
  ################################################################################################################
  ## take a try
  # loglh_func(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0))
  
  ## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
  
  
  ## take a try 
  
  ## take a try 
  # pars_sampler(n = 1)
  pars_prior <- createPrior(density = pars_density, sampler = function(){pars_sampler(n.stage = n_stage)}, 
                            lower = par_lower, upper = par_upper)
  
  
  bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
  
  
  
  ## DRAM: Adaptive MCMC, prior optimization, delayed rejection
  # startValue = c(b12=1.2, b3=0.4, b4=0.2, b5=0.1, r12=0.5, delta3=-1, delta4=0, delta5=0)
  # startValue = c(b12 = 1.359, b3 = 0.537, b4 = 0.203, b5 = 0.196, r12 = 0.305, delta3 = -0.964, delta4 = -0.593, delta5 = -0.309)
  # mh_settings = list(startValue = startValue,
  #                    adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10,
  #                    message = T)
  #startValue = trans_ori_to_delta(c(b_vec,r_vec,phi))
  mh_settings = list(
    #startValue = matrix(startValue,1,length(startValue)),
    iterations = n_iterations,
    startValues = true_pars,
    thin = 10,
    message = T
  )
  # mh_settings = list(startValue = startValue,
  #                    adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10,
  #                    message = T)
  #mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
  mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "DEzs", settings = mh_settings)
  #plot(mh_out)
  #plot(mh_out)
  #plot(mh_out)
  
  mcmc_pars_estimate <- getSample(mh_out)
  mcmc_pars_estimate <- mcmc_pars_estimate[(n_burn_in+2):nrow(mcmc_pars_estimate),]
 
   #loglh_func(colMeans(mcmc_pars_estimate))
  # colMeans(mcmc_pars_estimate_original)
  # loglh_func(trans_ori_to_delta(c(b_vec,r_vec,phi)))
  transform_delta_to_orginal=function(pars) {
    #n.stage <- (length(pars)-1)/2
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
    result <- pars
    result[1:(2*n.stage)] = c(b_vec,r_vec)
    return(result)
  }
  
  
  mcmc_pars_estimate_original = 
    t(apply(mcmc_pars_estimate,1,transform_delta_to_orginal))
  colnames(mcmc_pars_estimate_original) = c(paste0("b",1:n.stage),paste0("r",1:n.stage),"phi")
  
  
  # estSEAIP_mat <- apply(mcmc_pars_estimate, 1, function(x) SEIRsimu(pars = x, stage_intervals,
  #                                                                   Di,Dp,
  #                                                                   De,Dq_vec,
  #                                                                   alpha,Dh,N,flowN_vec,init_states,
  #                                                                   days_to_fit,num_periods = n.stage)[, c("S","E","I", "A", "P","Onset_expect")])
  # ptime <- 1:length(onset_obs)
  # estS_mat <- estSEAIP_mat[ptime,]
  # estE_mat <- estSEAIP_mat[ptime + length(ptime),]
  # estI_mat <- estSEAIP_mat[ptime + length(ptime)*2, ]
  # estA_mat <- estSEAIP_mat[ptime + length(ptime)*3, ]
  # estP_mat <- estSEAIP_mat[ptime + length(ptime) * 4, ]
  # estN_mat <- estSEAIP_mat[ptime + length(ptime) * 5, ]
  # prevalence <- rowMeans((N-estS_mat)/N)
  # prevalence_low <- apply(estS_mat,1,function(x){quantile((N-x)/N,0.025)})
  # prevalence_high <- apply(estS_mat,1,function(x){quantile((N-x)/N,0.975)})

  
  estRt_mat = apply(mcmc_pars_estimate_original,1,function(x){estimate_R(x,
                                                                         Di,
                                                                         Dp,
                                                                         Dq_vec,
                                                                         N,
                                                                         flowN_vec,
                                                                         n_stage)}) 
  
  
  est = colMeans(mcmc_pars_estimate_original)
  est_low = apply(mcmc_pars_estimate_original,2,function(x){quantile(x,0.025)})
  est_high = apply(mcmc_pars_estimate_original,2,function(x){quantile(x,0.975)})
  
  
  estRt = rowMeans(estRt_mat)
  Rt_low = apply(estRt_mat,1,function(x){quantile(x,0.025)})
  Rt_high = apply(estRt_mat,1,function(x){quantile(x,0.975)})
  
  # pars_est = colMeans(mcmc_pars_estimate)
  # pars_low = apply(mcmc_pars_estimate,2,function(x){quantile(x,0.025)})
  # pars_high = apply(mcmc_pars_estimate,2,function(x){quantile(x,0.975)})
  return(list(est,est_low,est_high,
              estRt,
              Rt_low,
              Rt_high,
              diagnosis = gelmanDiagnostics(mh_out)))
              # prevalence,
              # prevalence_low,
              # prevalence_high))
  
  
  #par_str=rep("c",n_pars)
  
  
}









