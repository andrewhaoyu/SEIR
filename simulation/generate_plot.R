setwd("/Users/zhangh24/GoogleDrive/covid_project/SEIR")

bias.table = read.csv("./output/simulation/bias_table.csv")
library(ggplot2)
n.stage = 6
method_vec = c(rep("nb",n.stage),rep("poisson",n.stage))
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
