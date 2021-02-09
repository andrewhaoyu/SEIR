estimate_R <- function(pars,
                       Di,
                       Dp,
                       Dq_vec,
                       N,
                       flowN_vec,
                       n_stage) {
  b_vec = pars[1:n_stage]
  r_vec = pars[(n_stage+1):(2*n_stage)]
  #
  R0_est <- rep(NA, n_stage)
  for(i in 1:n_stage) {
    b <- b_vec[i]
    r <- r_vec[i]
    Dq <- Dq_vec[i]
    n <- flowN_vec[i]
    R0_est[i] <- alpha * b / (1 / Dp + n / N) + (1 - r) * alpha * b / (1 / Di + n / N) + r * b / (1 / Di + 1 / Dq)
    rm(i, b, r, Dq, n)
  }
  return(R0_est)
}

