#library(cairoDevice)

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
      #b
      s_vec[, 1:n.stage] <- runif(n.stage, 0, 2) 
      #r1
      s_vec[, n.stage+1] <- rbeta(1, beta_shape1, beta_shape2)
      #r2..
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
    ## likelihood function
    loglh_func <- function(pars){
      
      
      ypred <- SEIRpred(pars, init_settings = init_sets_list)
      onset_obs <- init_sets_list$daily_new_case
      ypred <- ypred[, "Onset_expect"]
      phi = pars[length(pars)]
      #p = phi/(phi+as.numeric(ypred))
      # meant to suppress warnings when ypred is negative
      if(i1!=5){
        suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
                                      #size = phi,
                                      size = 1/phi,
                                      mu = ypred,log=T))
        
      }else{
        #special setting for CT due to the weekend effect
        suppressWarnings(p <- dnbinom(x = as.numeric(onset_obs), 
                                      #size = phi,
                                      size = 1/phi,
                                      mu = ypred,log=T))
        
      }
        
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
    pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)),paste0("phi"))
    return(list(default_pars_density,default_pars_sampler,loglh_func,pars_name))
  }
  
}

## wrapper for the analysis run
## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
#' @param init_set_list   initial settings produced by init_cond.R
#' @param randomize_startValue    this function will randomly generate an initial condtion if this argument is set to T; If you want to specify your own initial condition, set to F
#' @param startValue  If randomize_startValue is set to T, you can your own initial condition using this argument; If randomize_startValue is set to T, this argument will be ignored
#' @param output_ret  Whether to output parameter estimates output by MCMC
#' @param run_id this ID is meant to distinguish different runs (different parameters, random seeds, etc.). Run_ID will be included in the file names of all outputs
#' @param skip_MCMC This is meant for redrawing all results without rerunning MCMC
#' @param panel_B_R_ylim the y limit in panel B of the main result plot
SEIRfitting=function(init_sets_list, 
                     randomize_startValue=F, 
                     startValue=NA, 
                     output_ret=T, 
                     run_id, 
                     skip_MCMC=F, 
                     panel_B_R_ylim=4,
                     plot_combined_fig=T,
                     calc_clearance=T,
                     n_burn_in=20000,
                     n_iterations=320000,
                     all.date = all.date,
                     method = method) {
  if (randomize_startValue & !is.na(startValue)) {
    print("startValue will be ignored since randomize_startValue is set to TRUE!")
  } else if (!randomize_startValue & is.na(startValue)) {
    print("Please specify a startValue since you have set randomize_startValue to FALSE! Exiting!")
    q(save="no")
  }
  Fun_list_result = Func_list(method)
  pars_density=Fun_list_result[[1]]
  pars_sampler=Fun_list_result[[2]]
  loglh_func = Fun_list_result[[3]]
  pars_name=Fun_list_result[[4]]
  
  
  
  onset_obs <- init_sets_list$daily_new_case
  init_states <- init_sets_list$init_states
  n_pars = length(pars_name)
  n.stage = length(init_sets_list$stage_intervals)
  n_stage = n.stage
  ## take a try: pars = c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
  # SEIRpred(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0), init_settings = init_sets_list)[, "Onset_expect"]
  ################################################################################################################
  ## take a try
  # loglh_func(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0))
  
  ## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
  
  
  ## take a try 
  
  ## take a try 
  # pars_sampler(n = 1)
  pars_prior <- createPrior(density = pars_density, sampler = function(){pars_sampler(n.stage = n.stage)}, 
                            lower = init_sets_list$par_lower, upper = init_sets_list$par_upper)
  
  if (!skip_MCMC) {
    bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
  
    if (randomize_startValue) {  
      startValue=pars_sampler(n.stage = n.stage)
      best_logl = loglh_func(startValue)
      for(l in 1:1000){
        temp = pars_sampler(n.stage = n.stage)
        temp_logl = loglh_func(temp)
        if(best_logl<temp_logl){
          startValue = temp
        }
      }
      
      print(loglh_func(startValue))
      while (is.infinite(loglh_func(startValue))) {
        startValue=pars_sampler(n.stage = n.stage)
      }
    }
    
    ## DRAM: Adaptive MCMC, prior optimization, delayed rejection
    # startValue = c(b12=1.2, b3=0.4, b4=0.2, b5=0.1, r12=0.5, delta3=-1, delta4=0, delta5=0)
    # startValue = c(b12 = 1.359, b3 = 0.537, b4 = 0.203, b5 = 0.196, r12 = 0.305, delta3 = -0.964, delta4 = -0.593, delta5 = -0.309)
    # mh_settings = list(startValue = startValue,
    #                    adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10,
    #                    message = T)
    mh_settings = list(
      #startValue = startValue,
                       iterations = n_iterations,
                       
                       thin = 10,
                       message = T)
    #mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
    mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "DEzs", settings = mh_settings)
    #plot(mh_out)
    #plot(mh_out)
    mcmc_pars_estimate <- getSample(mh_out)
    mcmc_pars_estimate <- mcmc_pars_estimate[(n_burn_in+2):nrow(mcmc_pars_estimate),]
                                    #start = n_burn_in+2) 
    #mcmc_pars_estimate <- getSample(mh_out,                
                                    #start = n_burn_in+2, 
                                    #thin = 1)  ## set start = 2002 as the burn in period
    mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
    print("get sample finished")
    colnames(mcmc_pars_estimate) <- pars_name 
    
    if (output_ret) {
      if(i2==1){
        save(mh_out,file = paste0("../output/mcmc_out",run_id,".rdata"))  
      }
      write.table(mcmc_pars_estimate, paste0("../output/pars_est_run_",run_id,".txt"), quote = F, row.names = F, sep = "\t")
    }
  } else {
    mcmc_pars_estimate = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header = T)
    pars_name = names(mcmc_pars_estimate)
  }
  # n.par = ncol(mcmc_pars_estimate)
  # summary_string=paste0(paste(pars_name, collapse = ","), "\n")
  # 
  # par_str=list()
  # for (i_par in 1:n_pars) {
  #   par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate[,i_par]),2), " (",
  #                           round(quantile(mcmc_pars_estimate[,i_par],0.025),2)," - " , 
  #                           round(quantile(mcmc_pars_estimate[,i_par],0.975),2), ")")
  # }
  # print("summary string finished")
  # summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
  # 
  # estRt_mat <- apply(mcmc_pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_sets_list))
  # 
  # summary_string = paste0(summary_string, paste0("stage",1:n_stage,collapse=","), "\n")
  # 
  # r_str=list()
  # 
  # if (n_stage>1) {
  #   for (i_stage in 1:n_stage) {
  #     r_str[[i_stage]]=paste0(round(mean(estRt_mat[i_stage,]),2), " (",
  #                             round(quantile(estRt_mat[i_stage,],0.025),2)," - " , 
  #                             round(quantile(estRt_mat[i_stage,],0.975),2), ")")
  #   }
  # } else {
  #   r_str[[1]]=paste0(round(mean(estRt_mat),2), " (",
  #                           round(quantile(estRt_mat,0.025),2)," - " , 
  #                           round(quantile(estRt_mat,0.975),2), ")")
  # }
  # 
  # summary_string = paste0(summary_string, paste(r_str,collapse = ", "),"\n\n")
  # print("estimate Rt finished")
  # # if (calc_clearance) {
  # #   clearance_date = Findzero(mcmc_pars_estimate, init_sets_list)
  # #   
  # #   summary_string = paste0(summary_string, paste(names(clearance_date), collapse = ", "))
  # #   
  # #   summary_string = paste0(summary_string, "\n", paste(clearance_date, collapse = ", "), "\n")
  # # }
  # 
  # write_file(summary_string, paste0("../output/summary_run_",run_id,".txt"))
  # 
  # cairo_pdf(paste0("../output/par_cor_run_",run_id,".pdf"),width=10,height=10)
  # correlationPlot_modified(mcmc_pars_estimate, scaleCorText = F)
  # dev.off()
  # print("plot correlation plot finished")
  # #png(paste0("../output/par_hist_run_",run_id,".png"))
  # pdf(paste0("../output/par_hist_run_",run_id,".pdf"),width = 9, height = 10)
  # par(mfrow = c(4, ceiling(n.par/4)))
  # for(i in 1:n_pars) {
  #   hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
  #   rm(i)
  # }
  # dev.off()
  # print("plot hist plot finished")
  # #png(paste0("../output/par_traj_run_",run_id,".png"), width=1000, height=500)
  # pdf(paste0("../output/par_traj_run_",run_id,".pdf"),width = 9, height = 10)
  # par(mfrow = c(4, ceiling(n.par/4)))
  # for(i in 1:n_pars) {
  #   plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
  #   rm(i)
  # }
  # dev.off()
  # print("plot tracj finished")
  # if (plot_combined_fig) {
  #   SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim,
  #            stage_intervals=stage_intervals,all.date = all.date)
  # }
  
  #par(mfrow = c(1, 1))
  # corrplot(cor(mcmc_pars_estimate))
  # pairs(mcmc_pars_estimate)
  
}









