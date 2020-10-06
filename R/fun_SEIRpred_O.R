## Deterministic SEIR model
## five periods: Jan 1-9 (index 1-9), Jan 10-22 (index 10-22), Jan 23-Feb 1 (index 23-32), Feb 2-16 (index 33-47), Feb 17- (index 48-60)
#' @param pars                    a vetor of parameters: c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
#' @param init_settings           a list of initial values and known parameters
#' @param Dp                      presymptomatic infectious period
#' @param Di                      symptomatic infectious period
#' @param De                      latent period
#' @param Dq                      duration from illness onset to hospitalization
#' @param Dh                      hospitalization period                   
#' @param alpha                   ratio of the transmission rate of unascertained over ascertained case
#' @param N                       population size
#' @param flowN_vec               daily inbound and outbound size during five periods (n)
#' @param init_states             initial c(S, E, P, Is, A, H, R)
#' @param days_to_fit             the days to fit the model
#' @param b                       transmission rate of ascertained cases
#' @param r                       ascertainment rate  
#' @param stage_intervals         default corresponding staget intervals (Jan1, Jan9), (Jan10, Jan22), (Jan23, Feb1), (Feb2, Feb16), (Feb17, last day (varying))
#' 
#################################################################################################
SEIRpred <- function(pars, 
                     init_settings) {
  
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  tmp_ret=init_settings$var_trans_fun(pars)
  #transmission rate
  beta_vec=tmp_ret[[1]]
  #asceratinment rate from presymptomatic
  pi_p_vec=tmp_ret[[2]]
  #asceratinment rate from symptom onset
  pi_o_vec=tmp_ret[[3]]
  
  delta_vec = tmp_ret[[4]]
  
  #Ke latent period
  Ke <- init_settings$Ke
  #Kp presympotomatic period
  Kp <- init_settings$Kp
  #Ks sympotomatic infectious period
  Ks <- init_settings$Ks
  #Ko from sympotom onset to ascertained period
  Ko_vec <- init_settings$Ko_vec
  #Ka from ascertained to isloation or hospitalization
  Ka <- init_settings$Ka
  #Kh hospitalization period
  Kh <- init_settings$Kh
  #Ki infectious period after isolation
  Ki <- init_settings$Ki
  #total population
  N <- init_settings$N
  
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  ## ODE function based on deterministic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    beta <- stage_pars[1]
    pi_p <- stage_pars[2]
    pi_o <- stage_pars[3]
    delta <- stage_pars[4]
    Ko <- stage_pars[5]
    
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    O <- states_old[4]
    U <- states_old[5]
    A <- states_old[6]
    I <- states_old[7]
    H <- stages_old[8]
    R <- stages_old[9]
    ## new values
    S_new <- S - beta * S * (alpha * P + alpha* U + O + A + 0.5 * alpha * I) / N 
    E_new <- E + beta * S * (alpha * P + alpha* U + O + A + 0.5 * alpha * I) / N  - E / Ke 
    P_new <- P +  E / Ke  - P / Kp 
    U_new <- U + (1 - pi_p - pi_o) * P / Kp - U / Ks
    O_new <- O + pi_p * P / Kp - O / Ko
    A_new <- A + pi_o * P / Kp + O / Ko - A / Ka
    I_new <- I + delta * A / Ka - I / Ki
    H_new <- H + (1 - delta) *A / Ka - H / Kh
  
    R_new <- R + H / Kh + I / Ki + U / Ks
    Ascertain_expect <- pi_p * P / Kp + pi_o * O / Ko
    ##
    return(c(S_new, E_new, P_new, U_new, O_new, A_new, I_new, H_new, R_new, Ascertain_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "U", "O", "A", "I", "H", "R", "Ascertain_expect")
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    stage_pars_setings <- c(beta = beta_vec[i_stage], pi_p = pi_p_vec[i_stage], pi_o = pi_o_vec[i_stage], delta = delta_vec[i_stage])
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  return(states_mat)
}



