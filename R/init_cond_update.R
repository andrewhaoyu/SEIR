generate_init_condi <- function(r0,
                                Di = 2.9,
                                Dp = 4,
                                De = 2.9,
                                Dq ,
                                alpha = 0.55,
                                Dh = 30,
                                N ,
                                flowN ,
                                stateInput,
                                stateDataClean,
                                jan1_idx,
                                stage_intervals,
                                stateData,
                                method
) {
  
  stopifnot(r0>=0 & r0<=1 & Di>=0 & Dp>=0 & De>=0 & all(Dp>=0) & alpha>=0 & alpha<=1 & Dh>=0 & N>=0 & all(flowN>=0))
  
  ## N            : population size
  ## H0           : initial number of hospitalized cases based on the reports
  ## R0           : initial number of removed individuals
  ## De           : latent period
  ## r0           : initial ascertainment rate
  ## realData     : real data from the CDC
  realData_all <- stateData %>% select(positiveIncrease)
  R0 <- 0
  H0 <- 0
  n.stage <- length(stage_intervals)
  E0 <- sum(realData_all[(jan1_idx+round(Dp)):(jan1_idx+round(Dp)+round(De)-1),1]) / r0 ## Jan 3-5 for De=2.9 and Dp=2.3
  # E0 <- (40 + 23 + 47) / r0                                              
  P0 <- sum(realData_all[jan1_idx:(jan1_idx+round(Dp)-1),1]) / r0                     ## Jan 1-2 for Dp=2.3
  # P0 <- (41 + 34) / r0
  I0 <- sum(realData_all[(jan1_idx-round(Di)):(jan1_idx-1),1])                             ## Dec 29-31 for Di=2.9
  # I0 <- 11 + 13 + 10                                     
  A0 <- I0 * (1 - r0) / r0
  S0 <- N - E0 - P0 - I0 - A0 - H0 - R0
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0), 0)
  
  #leave 10 days for prediction
  realData_all <- stateDataClean
  realData <- realData_all[1:(nrow(realData_all)-10),]
  daily_new_case <- realData$positiveIncrease
  daily_new_case_all <- realData_all$positiveIncrease
  realData_all <- realData_all %>% select(positiveIncrease)
  realData <- realData %>% select(positiveIncrease)
  ##
  
  ## helper function
  # transform variables to a form that SEIRpred can use
  # so that SEIRpred can be re-used as much as possible
 
  
  if(method=="possion"){
    par_lower = c(rep(0,n.stage),0,rep(-10,n.stage-1))
    par_upper =  c(rep(2,n.stage),1,rep(10,n.stage-1))
    transform_var_main_stage=function(pars) {
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
  }else if(method=="nb"){
    par_lower = c(rep(0,n.stage),0,rep(-10,n.stage-1),0)
    par_upper =  c(rep(2,n.stage),1,rep(10,n.stage-1),500)
    transform_var_main_stage=function(pars) {
      n.stage <- (length(pars)-1)/2
      b_vec <- pars[1:n.stage]
      r_vec <- pars[(n.stage+1):(2*n.stage)]
      r1 = r_vec[1]
      for(l in 2:n.stage){
        rtemp = 1 / (1 + (1 - r_vec[l-1]) / (r_vec[l-1] * exp(r_vec[l])))
        r_vec[l] = rtemp
      }
      
      return(list(b_vec, r_vec))
    }
  }
  
  return(list(Di=Di,
              Dp=Dp,
              De=De,
              Dq=Dq,
              alpha=alpha,
              Dh=Dh,
              N=N,
              flowN=flowN,
              daily_new_case = daily_new_case, 
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              days_to_fit=1:length(daily_new_case),
              stage_intervals=stage_intervals,
              var_trans_fun=transform_var_main_stage,
              par_lower = par_lower,
              par_upper =  par_upper)
  )
  # TODO: please confirm the following:
  # boundaries for delta3-5 will not be used, they are here merely to meet the formality imposed by runMCMC
}

# get_init_sets_list is an alias of generate_init_condi in order not to break exsiting code
get_init_sets_list = generate_init_condi

delta_mean <- 0
delta_sd <- 1
beta_shape1 <- 7.3
beta_shape2 <- 24.6
gamma_shape = 5
gamma_rate = 1
