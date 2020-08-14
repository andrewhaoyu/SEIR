mcmc_pars_estimate = read.table(paste0("../output/pars_est_run_",run_id,".txt"),header=T)
new.data <- cbind(estN_mean,onset_obs_all)
head(new.data)
colnames(new.data) <- c("estN","onset")
estN <- new.data[,1]
onset <- new.data[,2]
plot(new.data[,1],new.data[,2])
outcome <- ((onset-estN)^2-onset)/estN
model <- lm(outcome~estN)
summary(model)
