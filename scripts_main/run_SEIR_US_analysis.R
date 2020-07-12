## IMPORTANT: set code_root properly!
# code_root="/home/dgwu/covid19/NatSEIR_Rcode/"
#code_root="~/jianguoyun/Nutstore/covid19/NatSEIR_Rcode/"
#code_root="/Users/zhangh24/Desktop/codes_covid19_v2/"
code_root="/data/zhangh24/SEIR/"
# code_root="C:/Users/xingj/Documents/WangLabAdmin/COVID-19/NatSEIR_Rcode/"

setwd(paste0(code_root, "scripts_main"))
#install.packages("BayesianTools")

library(BayesianTools)
#install.packages("vioplot")
library(vioplot)
#install.packages("corrplot")
library(corrplot)
library(readr)
#install.packages("cairoDevice")
library(cairoDevice)

##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu.R"))
source(paste0(code_root, "R/fun_SEIRfitting.R"))
source(paste0(code_root, "R/init_cond_update.R"))
source(paste0(code_root, "R/fun_R0estimate.R"))
source(paste0(code_root, "R/correlationPlot_modified.R"))
source(paste0(code_root, "R/fun_SEIRplot.R"))
source(paste0(code_root, "R/fun_Findzero.R"))
##
library(dplyr)
init_sets_list=get_init_sets_list(r0 = 0.23)

# good initial conditions
# c(1.284, 0.384, 0.174, 0.096, 0.161, -0.046, -0.379, 0.569)

SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = "main_analysis", output_ret = T, skip_MCMC=F)

## to evaluate convergence, we run another two rounds of this program
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
# SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)
