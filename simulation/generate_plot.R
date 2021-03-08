setwd("/Users/zhangh24/GoogleDrive/covid_project/SEIR")

load("./output/simulation/sum_two_stage_2.rdata")
library(ggplot2)
n.stage = 6
method_vec = c(rep("nb",n.stage),rep("poisson",n.stage))
pars = c(0.70,0.30,0.25,0.32,0.20)


bias = result$est_mean-pars
data.plot = data.frame(result,var_name = factor(row.names(result),
                                                levels = c("b1","b2","r1","r2","phi")),
                       Bias = result$est_mean-pars)
library(ggplot2)
p1 = ggplot(data.plot)+
  geom_bar(aes(var_name,Bias),
           stat = "identity",
           position = position_dodge())+
  theme_Publication()+
  xlab(NULL)+
  ggtitle("Bias of parameters")

p2 = ggplot(data.plot)+
  geom_bar(aes(var_name,est_cover),
           stat = "identity",
           position = position_dodge())+
  theme_Publication()+
  xlab(NULL)+
  ylab("Coverage")+
  ggtitle("Coverage of parameters")+
  coord_cartesian(ylim=c(0.90,1))+
  geom_hline(yintercept = 0.95, col="red",linetype ="dashed")
library(dplyr)
data.plot.sub = data.plot %>% 
  mutate(var_name = var_name) %>% 
  filter(var_name %in%c("r1","r2")) %>% 
  mutate(var_name = factor(var_name,levels = c("r1","r2")))
p3 = ggplot(data.plot,aes(x = var_name,y = est_mean))+
  geom_bar(
               stat = "identity",
               position=position_dodge())+
  geom_errorbar(aes(ymin=est_low,ymax=est_high))+
  theme_Publication()+
  xlab(NULL)+
  ylab("Ascertainment (95% CI)")+
  ggtitle("Ascertainment (95% CI)")+
  ylim(c(0,1))



for(i in 1:4){
  b_vec = c(bias.table[i,1:n.stage],bias.table[i+4,1:n.stage])
  r_vec = c(bias.table[i,(n.stage+1):(2*n.stage)],bias.table[i+4,(n.stage+1):(2*n.stage)])
  stages = rep(paste0("Stage ",1:n.stage),2)
  library(ggplot2)
  data.plot = data.frame(stages,method_vec,
                         transmission_rate = as.numeric(b_vec),
                         ascertainment_rate = as.numeric(r_vec))
  
}
ggplot(data.plot,aes(x = stages, y = transmission_rate))+
  geom_bar(stat = "identity",fill = method_vec)+
  theme_Publication()
