setwd("/Users/zhangh24/GoogleDrive/covid_project/SEIR")
#i2: 1 means 15 days in each stage
#i2: 2 means 30 days in each stage
#i3: 1 means phi as inverse gamma (1,1), delta as normal (0,1)
#i3: 2 means phi as inverse gamma (1,3), delta as normal (0,0.5)
pars = c(0.70,0.30,0.25,0.32,0.20)
for(i2 in 1:2){
  load(paste0("./output/simulation/sum_two_stage_",i2,"_",1,".rdata"))    
  data.plot1 = data.frame(result,var_name = factor(row.names(result),
                                                  levels = c("b1","b2","r1","r2","phi")),
                         Bias = result$est_mean-pars,
                         method = "Noninformative prior")
  load(paste0("./output/simulation/sum_two_stage_",i2,"_",2,".rdata"))    
  data.plot2 = data.frame(result,var_name = factor(row.names(result),
                                                  levels = c("b1","b2","r1","r2","phi")),
                         Bias = result$est_mean-pars,
                         method = "Informative prior")
  data.plot <- rbind(data.plot1,data.plot2)
    p1 = ggplot(data.plot)+
      geom_bar(aes(var_name,Bias,fill =method),
               stat = "identity",
               position = position_dodge())+
      theme_Publication()+
      xlab(NULL)+
      ggtitle("Bias of parameters")+
      scale_fill_Publication()
    png(filename = paste0("./output/simulation/bias_compare_",i2,".png"),width = 8, height =6, res=300,units = "in")
    print(p1)
    dev.off()
    p2 = ggplot(data.plot)+
      geom_bar(aes(var_name,est_cover,fill =method),
               stat = "identity",
               position = position_dodge())+
      theme_Publication()+
      xlab(NULL)+
      ylab("Coverage")+
      ggtitle("Coverage of parameters")+
      coord_cartesian(ylim=c(0.90,1))+
      geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
      scale_fill_Publication()
    png(filename = paste0("./output/simulation/coverage_compare_",i2,".png"),width = 8, height =6, res=300,units = "in")
    print(p2)
    dev.off()
    library(dplyr)
    
    p3 = ggplot(data.plot,aes(x = var_name,y = est_mean,ymin=est_low,ymax=est_high,fill=method))+
      geom_bar(
        stat = "identity",
        position=position_dodge())+
      geom_errorbar(position=position_dodge())+
      theme_Publication()+
      xlab(NULL)+
      ylab("Ascertainment (95% CI)")+
      ggtitle("Ascertainment (95% CI)")+
      ylim(c(0,1))+
      scale_fill_Publication()
    png(filename = paste0("./output/simulation/est_compare_",i2,".png"),width = 8, height =6, res=300,units = "in")
    print(p3)
    dev.off()
    
  
}

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
