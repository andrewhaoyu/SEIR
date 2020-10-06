Func_list = function(method){
  if(method=="poisson"){
    default_pars_density <- function(pars) {
      n.stage = length(pars)/2
      d_vec <- rep(NA, length(pars))
      ##b12, b3, b4, b5
      #bvec
      # for(i in c(1:n.stage)) {
      #   d_vec[i] <- dunif(pars[i], 0, 2, log = T)
      # }
      d_vec[1:n.stage] =  log(1/2)
      
      #rvec1
      d_vec[n.stage+1] = dbeta(pars[n.stage+1],beta_shape1, beta_shape2, log = T)
      #rvec
      d_vec[(n.stage+2):(2*n.stage)] = 
        dnorm(pars[(n.stage+2):(2*n.stage)],delta_mean, delta_sd, log = T)
      #rvec
      # for(l in (n.stage+1):(2*n.stage)){
      #   r_temp = pars[l]
      #   d_vec[l] = dbeta(r_temp, beta_shape1, beta_shape2, log = T)
      # }
      return(sum(d_vec))
    }
    
    default_pars_sampler <- function(n.stage=n.stage) {
      s_vec <- matrix(NA, 1, 2*n.stage)
      
      ## b12, b3, b4, b5
      s_vec[, 1:n.stage] <- runif(n.stage, 0, 2) 
      s_vec[, n.stage+1] <- rbeta(1, beta_shape1, beta_shape2)
      s_vec[, (n.stage+2):(2*n.stage)] <- rnorm(n.stage-1, delta_mean, delta_sd)
      return(s_vec)
    }
    ## likelihood function
    loglh_func <- function(pars){
      ypred <- SEIRpred(pars, init_settings = init_sets_list)
      ypred <- ypred[, "Onset_expect"]
      onset_obs <- init_sets_list$daily_new_case
      # meant to suppress warnings when ypred is negative
      suppressWarnings(p <- dpois(as.numeric(onset_obs), ypred,log=T))
      
      #if(any(p == 0) || any(is.nan(p))){
      if(any(is.nan(p))){  
        logL <- -Inf
      }else{
        logL <- sum(p)
      }
      return(logL)
    }
    pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)))
    return(list(default_pars_density,default_pars_sampler,loglh_func,pars_name))
  }else if(method=="nb"){
    default_pars_density <- function(pars) {
      n.stage = (length(pars)-1)/2
      d_vec <- rep(NA, length(pars))
      ##b12, b3, b4, b5
      #bvec
      # for(i in c(1:n.stage)) {
      #   d_vec[i] <- dunif(pars[i], 0, 2, log = T)
      # }
      d_vec[1:n.stage] =  log(1/2)
      
      #rvec1
      d_vec[n.stage+1] = dbeta(pars[n.stage+1],beta_shape1, beta_shape2, log = T)
      #rvec
      d_vec[(n.stage+2):(2*n.stage)] = 
        dnorm(pars[(n.stage+2):(2*n.stage)],delta_mean, delta_sd, log = T)
      phi = pars[2*n.stage+1]
      d_vec[2*n.stage+1] = dgamma(phi,gamma_shape,gamma_rate,log =T)
      return(sum(d_vec))
      #rvec
      # for(l in (n.stage+1):(2*n.stage)){
      #   r_temp = pars[l]
      #   d_vec[l] = dbeta(r_temp, beta_shape1, beta_shape2, log = T)
      # }
      
    }
    
    default_pars_sampler <- function(n.stage=n.stage) {
      s_vec <- matrix(NA, 1, 2*n.stage+1)
      
      ## b12, b3, b4, b5
      s_vec[, 1:n.stage] <- runif(n.stage, 0, 2) 
      s_vec[, n.stage+1] <- rbeta(1, beta_shape1, beta_shape2)
      s_vec[, (n.stage+2):(2*n.stage)] <- rnorm(n.stage-1, delta_mean, delta_sd)
      s_vec[, 2*n.stage+1] <- rgamma(1,gamma_shape,gamma_rate)
      return(s_vec)
    }
    ## likelihood function
    loglh_func <- function(pars){
      
      
      ypred <- SEIRpred(pars, init_settings = init_sets_list)
      onset_obs <- init_sets_list$daily_new_case
      ypred <- ypred[, "Onset_expect"]
      phi = pars[length(pars)]
      #p = phi/(phi+as.numeric(ypred))
      # meant to suppress warnings when ypred is negative
      suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
                                    size = phi,
                                    mu = ypred,log=T))
      
      #if(any(p == 0) || any(is.nan(p))){
      if(any(is.nan(p))){  
        logL <- -Inf
      }else{
        logL <- sum(p)
      }
      return(logL)
    }
    pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)),paste0("phi"))
    return(list(default_pars_density,default_pars_sampler,loglh_func,pars_name))
  }
  
}

#debug the results
init_sets_list;
randomize_startValue=T;
startValue=NA;
output_ret=T;
run_id = paste0(i1,"_",i2,"_",i3);
skip_MCMC=F;
panel_B_R_ylim=4;
plot_combined_fig=T;
calc_clearance=T;
n_burn_in=130000;
n_iterations=1500000;
all.date = all.date;
method = method;
mcmc_result <- read.table(paste0("../output/pars_est_run_",run_id,".txt"),header=T)
startValue = colMeans(mcmc_result)
init_settings = init_sets_list
pars = startValue 
ypred <- SEIRpred(pars, init_settings = init_sets_list)
onset_obs <- init_sets_list$daily_new_case
ypred <- ypred[, "Onset_expect"]
phi = pars[length(pars)]
#p = phi/(phi+as.numeric(ypred))
# meant to suppress warnings when ypred is negative
p_fit <- dnbinom(x = as.numeric(onset_obs), 
                 size = phi,
                 mu = ypred,log=T)

p_max <- dnbinom(x = as.numeric(onset_obs), 
                 size = phi,
                 mu = as.numeric(onset_obs),log=T)
D = 2*(p_max-p_fit)
D_residual = ifelse(onset_obs-ypred>=0,1,-1)*sqrt(D)
time = c(1:length(D_residual))
plot(time,D_residual)

