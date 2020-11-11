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
                      panel_B_R_ylim=7,
                      all.date = all.date) {
  init_settings = init_sets_list
  plot_combined_fig=T
  calc_clearance=T
  
  onset_obs <- init_sets_list$daily_new_case
  init_states <- init_sets_list$init_states
  n.stage = length(init_sets_list$stage_intervals)
  stage_intervals = init_sets_list$stage_intervals
  n_stage = n.stage
  #library(coda)
  run_id_split = strsplit(run_id,split="_")
  date.input <- run_id_split[[1]][1]
  id1 = run_id_split[[1]][2]
  id3 = run_id_split[[1]][4]
  pars_estimate_main=read.table(paste0("../output/pars_est_run_",date.input,"_",id1,"_",1,"_",id3,".txt"), header=T)
  pars_estimate_main_rep1=read.table(paste0("../output/pars_est_run_",date.input,"_",id1,"_",2,"_",id3,".txt"), header=T)
  pars_estimate_main_rep2=read.table(paste0("../output/pars_est_run_",date.input,"_",id1,"_",3,"_",id3,".txt"), header=T)
  # mcmc_main=mcmc(data=pars_estimate_main)
  # mcmc_rep1=mcmc(data=pars_estimate_main_rep1)
  # mcmc_rep2=mcmc(data=pars_estimate_main_rep2)
  # mcmc_3traj=mcmc.list(mcmc_main,mcmc_rep1,mcmc_rep2)
  # mcmc_3traj=mcmc.list(mcmc_main,mcmc_rep1,mcmc_rep2)
  # gelman.diag(mcmc_3traj)
  p_vec <- c(1:nrow(pars_estimate_main))
  plot_par_3traj = function(par_name, plotmath_name) {
    # red-like: #BC3C29, blue-like: #0072B5, orange-like: #E18727
    plot(p_vec, pars_estimate_main[p_vec, par_name], type="l", col="#0072B5", main=plotmath_name, xlab="", ylab="")
    points(p_vec, pars_estimate_main_rep1[p_vec, par_name], type="l", col="#BC3C29")
    points(p_vec, pars_estimate_main_rep2[p_vec, par_name], type="l", col="#E18727")
  }
  png(paste0("../output/mcmc_convergence",run_id,".png"), width=15, height=10,res = 300, units = "in")
  par(mfrow=c(3, ceiling(ncol(pars_estimate_main)/3)))
  for(k in 1:ncol(pars_estimate_main)){
    plot_par_3traj(colnames(pars_estimate_main)[k], colnames(pars_estimate_main)[k])
  }
  dev.off()
  mcmc_pars_estimate  = read.table(paste0("../output/pars_est_run_",run_id,".txt"), header=T)
  pars_name = colnames(mcmc_pars_estimate)
  n.par = ncol(mcmc_pars_estimate)
  n_pars = length(pars_name)
  #summary_string=paste0(paste(pars_name, collapse = ","), "\n")
  
  
  #get the b and r under orignal scale
  
  
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
  par_str=rep("c",n_pars)
  

  
  for (i_par in 1:n_pars) {
    par_str[i_par]=paste0(round(mean(mcmc_pars_estimate_original[,i_par]),2), " (",
                          round(quantile(mcmc_pars_estimate_original[,i_par],0.025),2)," - " , 
                          round(quantile(mcmc_pars_estimate_original[,i_par],0.975),2), ")")
  }
  names(par_str) = colnames(mcmc_pars_estimate_original)
  print("summary string finished")
  ascertainment = rep(0,n_stage)
  ascertainment_low = rep(0,n_stage)
  ascertainment_high = rep(0,n_stage)
  stage_date = all.date[1:n_stage]
  for(i in 1:n_stage){
    ascertainment[i] <- round(mean(mcmc_pars_estimate_original[,n_stage+i]),2)
    ascertainment_low[i] <- round(quantile(mcmc_pars_estimate_original[,n_stage+i],0.025),2)
    ascertainment_high[i] <- round(quantile(mcmc_pars_estimate_original[,n_stage+i],0.975),2)
    stage_date[i] = all.date[as.integer((stage_intervals[[i]][1]+stage_intervals[[i]][2])/2)]
  }
  plot.data <- data.frame(stage_date,ascertainment,
                          ascertainment_low,
                          ascertainment_high)
  
  
  p = ggplot(plot.data,aes(x= stage_date,y =ascertainment ))+
    geom_line()+
    geom_ribbon(aes(ymin=ascertainment_low,ymax=ascertainment_high),alpha = 0.2)+
    xlab("Time-period")+
    ylab("Ascertainment (95%CI)")+
    theme_Publication()
  png(file = paste0("ascertainment_",run_id,".png"),width = 10,height =8, res = 300, units = "in")
  print(p)
  dev.off()
  
  n.days <- length(init_settings$daily_new_case_all)
  stages <- length(stage_intervals)
  start.vec <- rep(0,n_stage)
  end.vec <- rep(0,n_stage)
  for(k in 1:stages){
    start.vec[k] <- stage_intervals[[k]][1]
    end.vec[k] <- stage_intervals[[k]][2]
  }
  #days cut = total days/stages
  init_settings$days_to_fit <- 1:n.days
  onset_obs_all <- as.numeric(init_settings$daily_new_case_all)
  ptime <- 1:length(onset_obs_all)
  mydate <- format(all.date,"%b %d")
  stage.label <- rep("c",n_stage)
  for(l in 1:n_stage){
    stage.label[l] = paste0(format(all.date[start.vec[l]],"%b %d")," - ",format(all.date[end.vec[l]],"%d"))
  }
  
  rt.time <- data.frame(stage.label,ascertainment_vec)
  
  
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
  #pdf(paste0("../output/par_hist_run_",run_id,".pdf"),width = 9, height = 10)
  # par(mfrow = c(4, ceiling(n.par/4)))
  # for(i in 1:n_pars) {
  #   hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
  #   rm(i)
  # }
  # dev.off()
  # print("plot hist plot finished")
  #png(paste0("../output/par_traj_run_",run_id,".png"), width=1000, height=500)
  # pdf(paste0("../output/par_traj_run_",run_id,".pdf"),width = 9, height = 10)
  # par(mfrow = c(4, ceiling(n.par/4)))
  # for(i in 1:n_pars) {
  #   plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
  #   rm(i)
  # }
  # dev.off()
 # print("plot tracj finished")
  if (plot_combined_fig) {
    SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim,
             stage_intervals=stage_intervals,all.date = all.date)
  }
  
  par(mfrow = c(1, 1))
  corrplot(cor(mcmc_pars_estimate))
  pairs(mcmc_pars_estimate)
  
}









