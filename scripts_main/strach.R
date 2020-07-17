r0=0.23;
Di = 2.9;
Dp = 2.3;
De = 2.9;
Dq = c(21, 15, 10, 6, 2);
alpha = 0.55;
Dh = 30;
N = 10000000;
flowN = c(500000, 800000, 0, 0, 0)
R0 <- 0
H0 <- 27
realData <- read.csv("../data/Covid19CasesWH_simulated_Jan1-Feb29.csv", row.names = 1)[-c(1:24), ]
realData_all <- read.csv("../data/Covid19CasesWH_simulated_Jan1-Feb29.csv", row.names = 1)  # the 25th row correspond to 1 Jan

jan1_idx = 25






randomize_startValue = T;
run_id = "main_analysis"; 
output_ret = T;
skip_MCMC=F


init_sets_list;
randomize_startValue=T;
startValue=NA;
output_ret=T;
run_id=0;
skip_MCMC=F; 
panel_B_R_ylim=4;
plot_combined_fig=T;
pars_density=default_pars_density;
pars_sampler=default_pars_sampler;
pars_name=c(paste0("b",c(1:n.stage)),"r1",paste0("delta",c(2:n.stage)));
calc_clearance=T;
n_burn_in=12000;
n_iterations=180000


pars = startValue














Di=Di
Dp=Dp
De=De
Dq=Dq
alpha=alpha
Dh=Dh
N=N
flowN=flowN
daily_new_case = daily_new_case
daily_new_case_all = daily_new_case_all
init_states = init_states
days_to_fit=1:length(daily_new_case)
stage_intervals=stage_intervals
var_trans_fun=transform_var_main_stage
par_lower = c(rep(0,stages),0,rep(-10,stages-1))
par_upper =  c(rep(2,stages),1,rep(10,stages-1))
