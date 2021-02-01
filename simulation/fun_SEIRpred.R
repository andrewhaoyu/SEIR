SEIRpred <- function(stage_intervals,
                     b_vec,
                     r_vec,
                     Di,
                     Dp,
                     De,
                     Dq_vec,
                     alpha,
                     Dh,
                     N,
                     flowN_vec,
                     init_states,
                     days_to_fit) {
  
  n_stage = length(stage_intervals)
  num_periods = n_stage
  ## ODE function based on deterministic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    Dq <- stage_pars[3]
    n <- stage_pars[4]
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    H <- states_old[6]
    R <- states_old[7]
    ## new values
    S_new <- S - b * S * (alpha * P + I + alpha * A) / N + n - n * S / N
    E_new <- E + b * S * (alpha * P + I + alpha * A) / N - E / De - n * E / N
    P_new <- P +  E / De  - P / Dp - n * P / N
    I_new <- I + r * P / Dp - I / Di - I / Dq
    A_new <- A + (1 - r) * P / Dp - A / Di - n * A / N
    H_new <- H + I / Dq - H / Dh
    R_new <- R + H / Dh + (A + I) / Di - n * R / N
    Onset_expect <- r * P / Dp
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  return(states_mat)
}

