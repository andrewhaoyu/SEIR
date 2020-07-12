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
pars_name=c("b12", "b3", "b4", "b5", "r12", "delta3", "delta4", "delta5");
calc_clearance=T;
n_burn_in=4000;
n_iterations=180000


pars = startValue