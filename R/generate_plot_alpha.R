#library(cairoDevice)
## wrapper for the analysis run
## Create BayesianSetup and settings, lower/upper for parameters: b12, b3, b3, b5, r12, delta3, delta4, delta5
#' @param init_set_list   initial settings produced by init_cond.R
#' @param randomize_startValue    this function will randomly generate an initial condtion if this argument is set to T; If you want to specify your own initial condition, set to F
#' @param startValue  If randomize_startValue is set to T, you can your own initial condition using this argument; If randomize_startValue is set to T, this argument will be ignored
#' @param output_ret  Whether to output parameter estimates output by MCMC
#' @param run_id this ID is meant to distinguish different runs (different parameters, random seeds, etc.). Run_ID will be included in the file names of all outputs
#' @param skip_MCMC This is meant for redrawing all results without rerunning MCMC
#' @param panel_B_R_ylim the y limit in panel B of the main result plot
GeneratePlot=function(init_sets_list, 
                    run_id = runid,
                     panel_B_R_ylim=6,
                     all.date = all.date) {
  
  plot_combined_fig=T
  calc_clearance=T
  
  onset_obs <- init_sets_list$daily_new_case
  init_states <- init_sets_list$init_states
  n.stage = length(init_sets_list$stage_intervals)
  stage_intervals = init_sets_list$stage_intervals
  n_stage = n.stage
  mcmc_pars_estimate  = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header=T)
  pars_name = colnames(mcmc_pars_estimate)
  n.par = ncol(mcmc_pars_estimate)
  n_pars = length(pars_name)
  #summary_string=paste0(paste(pars_name, collapse = ","), "\n")
  
  
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  #get the b and r under orignal scale
  mcmc_pars_estimate_original <- mcmc_pars_estimate
  
  transform_delta_to_orginal=function(pars) {
    n.stage <- (length(pars)-2)/2
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
  colnames(mcmc_pars_estimate_original) = c(paste0("b",1:n.stage),paste0("r",1:n.stage),"phi","alpha")
  par_str=rep("c",n_pars)
  for (i_par in 1:n_pars) {
    par_str[i_par]=paste0(round(mean(mcmc_pars_estimate_original[,i_par]),2), " (",
                            round(quantile(mcmc_pars_estimate_original[,i_par],0.025),2)," - " , 
                            round(quantile(mcmc_pars_estimate_original[,i_par],0.975),2), ")")
  }
  names(par_str) = colnames(mcmc_pars_estimate_original)
  print("summary string finished")

  
  
  estRt_mat <- apply(mcmc_pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_sets_list))
  
  
  
  r_str=rep("c",n_stage)
  
  if (n_stage>1) {
    for (i_stage in 1:n_stage) {
      r_str[i_stage]=paste0(round(mean(estRt_mat[i_stage,]),2), " (",
                              round(quantile(estRt_mat[i_stage,],0.025),2)," - " , 
                              round(quantile(estRt_mat[i_stage,],0.975),2), ")")
    }
  } else {
    r_str[1]=paste0(round(mean(estRt_mat),2), " (",
                      round(quantile(estRt_mat,0.025),2)," - " , 
                      round(quantile(estRt_mat,0.975),2), ")")
  }
  
  names(r_str) <- paste0("Rt",c(1:n_stage))
  
  summary_string <- c(par_str,r_str)
  
  write.csv(summary_string, paste0("../output/summary_run_",run_id,".csv"))

  cairo_pdf(paste0("../output/par_cor_run_",run_id,".pdf"),width=10,height=10)
  correlationPlot_modified(mcmc_pars_estimate, scaleCorText = F)
  dev.off()
  print("plot correlation plot finished")
  #png(paste0("../output/par_hist_run_",run_id,".png"))
  pdf(paste0("../output/par_hist_run_",run_id,".pdf"),width = 9, height = 10)
  par(mfrow = c(4, ceiling(n.par/4)))
  for(i in 1:n_pars) {
    hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
    rm(i)
  }
  dev.off()
  print("plot hist plot finished")
  #png(paste0("../output/par_traj_run_",run_id,".png"), width=1000, height=500)
  pdf(paste0("../output/par_traj_run_",run_id,".pdf"),width = 9, height = 10)
  par(mfrow = c(4, ceiling(n.par/4)))
  for(i in 1:n_pars) {
    plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
    rm(i)
  }
  dev.off()
  print("plot tracj finished")
  if (plot_combined_fig) {
    SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim,
             stage_intervals=stage_intervals,all.date = all.date)
  }

  par(mfrow = c(1, 1))
  corrplot(cor(mcmc_pars_estimate))
  pairs(mcmc_pars_estimate)

}









