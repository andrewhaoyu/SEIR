init_sets_list;
randomize_startValue=T;
startValue=NA;
output_ret=T;
run_id = paste0(i1,"_",i2,"_",i3);
skip_MCMC=F;
panel_B_R_ylim=4;
plot_combined_fig=T;
calc_clearance=T;
n_burn_in=20000;
n_iterations=500000;
all.date = all.date;
method = method;
pars = startValue


library(data.table)
s.t <- list(data.frame(a=rnorm(100),b=rnorm(100),chain=rep("chain1",100),x=c(1:100)),
          data.frame(a=rnorm(100),b=rnorm(100),chain=rep("chain2",100),x=c(1:100)))
s.t.data <- rbindlist(s.t)

s.t.data.long <- melt(s.t.data,id.var=c("chain","x"))

ggplot(s.t.data.long)+
  geom_line(aes(x= x,y = value,color=chain))+
  facet_wrap(~variable)+
  theme_Publication()


library(ggmcmc)
data(radon)
s.radon.short <- radon$s.radon.short
S <- ggs(s.radon.short)
ggs_traceplot(S)
