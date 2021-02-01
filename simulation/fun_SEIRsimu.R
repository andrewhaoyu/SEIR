SEIRsimu <- function(stage_intervals,
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
  ## ODE function based on stochastic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    Dq <- stage_pars[3]
    n <- stage_pars[4]
    #n <- rpois(1, lambda = as.numeric(stage_pars[4]))     ## stochastic, Poisson Distribution
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    H <- states_old[6]
    R <- states_old[7]
    ## S
    ## meaning S->E, S->, S->S
    pS_vec <- c(min(1,b * (alpha * P + I + alpha * A) / N), n / N, max(0,1 - b * (alpha * P + I + alpha * A) / N - n / N))
    sample_S <- rmultinom(1, size = S, prob = pS_vec)
    ## E
    ## meaning E->P, E->, E->E
    pE_vec <- c(1 / De, n / N, 1 - 1 / De - n / N)
    sample_E <- rmultinom(1, size = E, prob = pE_vec)
    ## P
    ## meaning P->I, P->A, P->, P->P
    pP_vec <- c(r / Dp, (1 - r) / Dp, n/N, 1 - 1 / Dp - n/N)
    sample_P <- rmultinom(1, size = P, prob = pP_vec)
    ## I
    ## meaning I->H, I->R, I->I
    pI_vec <- c(1 / Dq, 1 / Di, 1 - 1 / Dq - 1 / Di)
    sample_I <- rmultinom(1, size = I, prob = pI_vec)
    ## A
    ## meaning A->R, A->, A->A
    pA_vec <- c(1 / Di, n / N, 1 - 1 / Di - n / N)
    sample_A <- rmultinom(1, size = A, prob = pA_vec)
    ## H
    ## meaning H->R, H->H
    pH_vec <- c(1 / Dh, 1 - 1 / Dh)
    sample_H <- rmultinom(1, size = H, prob = pH_vec)
    ## R
    ## meaning R->, R->R
    pR_vec <- c(n / N, 1 - n / N)
    sample_R <- rmultinom(1, size = R, prob = pR_vec)
    ## new values
    S_new <- sample_S[3] + n
    E_new <- sample_E[3] + sample_S[1]
    P_new <- sample_P[4] + sample_E[1]
    I_new <- sample_I[3] + sample_P[1]
    A_new <- sample_A[3] + sample_P[2]
    H_new <- sample_H[2] + sample_I[1]
    R_new <- sample_R[2] + sample_I[2] + sample_A[1] + sample_H[1]
    Onset_expect <- sample_P[1]
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
  ## evovle the system according to the discretized ODEs
  stage_start <- stage_end <- rep(0,n_stage)
  for(l in 1:n_stage){
    stage_start[l] = stage_intervals[[l]][1]
    stage_end[l] = stage_intervals[[l]][2]
  }
  stage_end[l]  = length(days_to_fit)
  myold_states <- init_states
  for (i_stage in 1:n_stage) {
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_start[i_stage]:stage_end[i_stage]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  if(num_periods == n_stage) {  ## total 5 periods: Jan 1-9, Jan 10-22, Jan 23-Feb 1, Feb 2-16, Feb 17-
    states_mat <- states_mat
  }
  ## num_periods=4: only 4 periods: Jan 1-9, Jan 10-22, Jan 23-Feb 1, Feb 2-
  ## num_periods=3: only 3 periods: Jan 1-9, Jan 10-22, Jan 23-
  ## num_periods=2: only 2 periods: Jan 1-9, Jan 10
  else if (num_periods %in% c(2:(n_stage-1))) {  
    i_stage <- num_periods
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_start[i_stage]:length(days_to_fit)) {
      myold_states <- states_mat[d - 1, -1]
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
    }
  }
  else {
    print("num_periods has to be within 2 to total number of stage!")
    #q(save="no")
  }
  
  return(states_mat)
}
