#Generate six stages analysis
n_stage = 6
n.stage = n_stage
initial.ascertainment = 0.20
N = 21477737
stage_intervals = list()
total= 0
for(i in 1:n_stage){
  stage_intervals[[i]] <- c(start = total+1,end = total+15)
  total = total+15
}
b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
r_vec = c(0.2,0.25,0.3,0.35,0.4,0.30)
Di = 3.5
Dp = 2.75
De = 2.45
Dq_vec = c(10,6,3,3,3,3)
alpha = 0.55
Dh = 30
flowN_vec = c(0,0,0,0,0,0)
days_to_fit = c(1:total)
init_states <- c(21470907,2295,3215,264,1056,0,0)
names(init_states) <- c("S","E","P","I","A","H","R")
SEIR_mat =SEIRsimu(stage_intervals,b_vec,r_vec,
         Di,Dp,
         De,Dq_vec,
         alpha,Dh,N,flowN_vec,init_states,
         days_to_fit)

onset_obs = SEIR_mat[,"Onset_expect"]
par_lower = c(0,rep(-10,n_stage-1),0,rep(-10,n_stage-1),0)
par_upper = c(3,rep(10,n_stage-1),1,rep(10,n_stage-1),1000)
n_iterations = 150000
n_burn_in = 13600
delta_mean <- 0
delta_sd <- 1
#beta_shape1 <- 7.3
#beta_shape2 <- 24.6
beta_shape1 <- 1
beta_shape2 <- 1
gamma_shape = 1
gamma_rate = 1



