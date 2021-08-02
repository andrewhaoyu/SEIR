b_vec = c(0.81,0.22,0.19,0.28,0.31,0.48)
r_vec = c(0.23,0.24,0.38,0.57,0.69,0.75)
n_stage = 6
phi_vec = 0.2
ascertainment_mat = c(0.23,0.24,0.38,0.57,0.69,0.75)
phi = 0.2
r_vec = c(0.23,0.24,0.38,0.57,0.69,0.75)

Di = 3.5
Dp = 2.75
De = 2.45
#Dq_vec = c(10,6,3,3,3,3)
n_stage  = n.stage = length(r_vec)
Dq_vec = c(10,6,3,3,3,3)
alpha = 0.55
Dh = 30
flowN_vec = rep(0,n.stage)
days_to_fit = c(1:total)
N = 21477737
i2 = 2
i3 = 1
i4 = 1
Rt = estimate_R(pars = c(b_vec,r_vec),
                Di = Di,
                Dp = Dp,
                Dq_vec = Dq_vec,
                N = N,
                flowN_vec = flowN_vec,
                n_stage = n_stage)

pars = c(b_vec,r_vec,phi)



true_Rt = Rt




n.rep =2000
#2 ascertainment settings
#2 overdispersion settings
#2 methods
total = n.rep
n.stage= 6
pars_est = matrix(0,total,n.stage*2+1)
pars_est_cover = matrix(0,total,n.stage*2+1)
Rt_est = matrix(0,total,n.stage)
Rt_cover = matrix(0,total,n.stage)
pre_est = matrix(0,total,n_stage*30)
pre_cover = matrix(0,total,n_stage*30)
temp = 1
files = dir(path = "/data/zhangh24/SEIR/result/simulation",pattern = "sixstage_two_stage_",full.names = T)
    for(i1 in 1:n.rep){
      file = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_two_stage_",i1,"_",i2,"_",i3,"_",i4,".rdata")
      if(file %in% files){
        load(file)
        if(est_result$diagnosis[2]<1.1){
          #true_prevalence = est_result[[10]]
          pars_est[temp,] = est_result[[1]]
          est_low = est_result[[2]]
          est_high = est_result[[3]]
          pars_est_cover[temp,] = (pars>=est_low)*(pars<=est_high)
          
          Rt_est[temp,] = est_result[[4]]
          Rt_low = est_result[[5]]
          Rt_high = est_result[[6]]
          Rt_cover[temp,] = (true_Rt>=Rt_low)*(true_Rt<=Rt_high)
          
          temp = temp+1
        }
        
        
        
      }
      
      #load(paste0("/data/zhangh24/SEIR/result/simulation/seir_pos_result_",i1,"_",i2,"_",i3,".rdata"))
    }
 
pars_est = pars_est[1:(temp-1),]
pars_est_cover = pars_est_cover[1:(temp-1),]
Rt_est = Rt_est[1:(temp-1),]
Rt_cover = Rt_cover[1:(temp-1),]
#bias plot
bias_est = colMeans(pars_est)-pars


data.plot = data.frame(x = paste0("b",c(1:n_stage)),
                       bias = bias_est[1:n_stage])
library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,bias),
           stat = "identity",
           )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Bias")+
  ggtitle("Bias of b")+
  coord_cartesian(ylim=c(-0.3,0.3))+
  #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_bias_b.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()

data.plot = data.frame(x = paste0("r",c(1:n_stage)),
                       bias = bias_est[(1:n_stage)+n.stage])
#library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,bias),
           stat = "identity",
  )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Bias")+
  ggtitle("Bias of r")+
  #coord_cartesian(ylim=c(0.90,1))+
  #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_bias_r.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()

cover_est = colMeans(pars_est_cover)

data.plot = data.frame(x = paste0("b",c(1:n_stage)),
                       cover = cover_est[1:n_stage])
#library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,cover),
           stat = "identity",
  )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Coverage")+
  ggtitle("Coverage of b")+
  coord_cartesian(ylim=c(0.70,1))+
  geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_cover_b.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()

data.plot = data.frame(x = paste0("r",c(1:n_stage)),
                       cover = cover_est[(1:n_stage)+n.stage])
#library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,cover),
           stat = "identity",
  )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Coverage")+
  ggtitle("Coverage of r")+
  coord_cartesian(ylim=c(0.65,1))+
  geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_cover_r.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()



bias_Rt = colMeans(Rt_est)-Rt
cover_Rt = colMeans(Rt_cover)



data.plot = data.frame(x = paste0("Rt",c(1:n_stage)),
                       bias = bias_Rt[1:n_stage])
library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,bias),
           stat = "identity",
  )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Bias")+
  ggtitle("Bias of Rt")+
  #coord_cartesian(ylim=c(0.90,1))+
  #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_bias_Rt.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()


data.plot = data.frame(x = paste0("Rt",c(1:n_stage)),
                       cover = cover_Rt[1:n_stage])
library(ggplot2)
p = ggplot(data.plot)+
  geom_bar(aes(x,cover),
           stat = "identity",
  )+
  theme_Publication()+
  xlab(NULL)+
  ylab("Coverage")+
  ggtitle("Coverage of Rt")+
  coord_cartesian(ylim=c(0.80,1))+
  geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  scale_fill_Publication()


png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_cover_Rt.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()




ascert_est = pars_est[,(1:n.stage)+1]
aecert_low = apply(ascert_est,2,function(x){quantile(x,0.05)})
aecert_high = apply(ascert_est,2,function(x){quantile(x,0.95)})
ascert_avg = colMeans(ascert_est)

data.plot = data.frame(x = c(1:n_stage),
                                  ascert_avg,
                                  aecert_low,
                                  aecert_high)

p = ggplot(data.plot,aes(x = x))+
  geom_ribbon(aes(ymin = aecert_low,
                  ymax = aecert_high), fill = "grey70")+
  geom_line(aes(y = ascert_avg))+
  theme_Publication()+
  xlab("Stage")+
  ylab("Ascertainment")+
  ggtitle("Ascertainment")

png(filename = paste0("/data/zhangh24/SEIR/result/simulation/sixstage_ascert_trend.png"),width = 10, height = 8, unit = "in",res = 300)
print(p)
dev.off()

  # coord_cartesian(ylim=c(0.80,1))+
  # geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
  #scale_fill_Publication()

bias_pre = colMeans(pre_est)-true_prevalence
cover_pre = colMeans(pre_cover)

result = list(bias_est,cover_est,
              bias_Rt,cover_Rt,
              bias_pre,
              cover_pre)







n_stage = 6


n.rep =1000
#2 ascertainment settings
#2 overdispersion settings
#2 methods
total = n.rep
n.stage= 6
pars_est = matrix(0,total,n.stage*2+1)
pars_est_cover = matrix(0,total,n.stage*2+1)
i2_vec = rep(0,total)
i3_vec = rep(0,total)
method_vec = rep(0,total)
temp = 1

    pars = c(b_vec,r_vec,phi)
    
    for(i1 in 1:n.rep){
      load(paste0("/data/zhangh24/SEIR/result/simulation/sixstage_two_stage_",i1,"_",i2,"_",i3,"_",i4,".rdata"))
      pars_est[temp,] = est_result[[1]]
      est_low = est_result[[2]]
      est_high = est_result[[3]]
      pars_est_cover[temp,] = (pars>=est_low)*(pars<=est_high)
      i2_vec[temp] = i2
      i3_vec[temp] = i3
      method_vec[temp] = "nb"
      temp = temp+1
      #load(paste0("/data/zhangh24/SEIR/result/simulation/seir_pos_result_",i1,"_",i2,"_",i3,".rdata"))
    }
 

var.names = names(est_result[[1]])

for(i2 in 1:2){
  for(i3 in 1:2){
    
    
    phi = phi_vec[i2]
    r_vec = ascertainment_mat[,i3]
    
    pars = c(b_vec,r_vec,phi)
    
    
    for(i1 in 1:n.rep){
      load(paste0("/data/zhangh24/SEIR/result/simulation/seir_pos_result_",i1,"_",i2,"_",i3,".rdata"))
      pars_est[temp,] = c(est_result[[1]],0)
      est_low = c(est_result[[2]],0)
      est_high = c(est_result[[3]],0)
      pars_est_cover[temp,] = (pars>=est_low)*(pars<=est_high)
      i2_vec[temp] = i2
      i3_vec[temp] = i3
      method_vec[temp] = "poisson"
      temp = temp+1
      #load(paste0("/data/zhangh24/SEIR/result/simulation/seir_pos_result_",i1,"_",i2,"_",i3,".rdata"))
    }
  }
  
}


result = data.frame(pars_est,pars_est_cover,i2_vec,i3_vec,method_vec)
colnames(result)[1:(n.stage*2+1)] = var.names
colnames(result)[(1:(n.stage*2+1))+(n.stage*2+1)] = paste0("cover_",var.names)
case_option = c("No overdispersion-Low ascertainment",
                "No overdispersion-High ascertainment",
                "Overdispersion-Low ascertainment",
                "Overdispersion-High ascertainment")

case = 4
bias_table = matrix(0,2*case,2*n.stage)
cover_table = matrix(0,2*case,2*n.stage)

temp =1 
for(method in c("nb","poisson")){
  for(i2 in 1:2){
    for(i3 in 1:2){
      
      result.sub = result %>% 
        filter(i2_vec == i2&
                 i3_vec==i3&
                 method_vec == method)
      phi = phi_vec[i2]
      r_vec = ascertainment_mat[,i3]
      
      pars = c(b_vec,r_vec)
      
      bias_table[temp,] = round(colMeans(result.sub[,1:(2*n.stage)])-pars,3)
      cover_table[temp,] = round(colMeans(result.sub[,(2*n.stage+1)+(1:(2*n.stage))]),3)
      temp = temp+1
    }
    
  }
}
colnames(bias_table) = var.names[1:(2*n.stage)]
colnames(cover_table) = var.names[1:(2*n.stage)]
write.csv(bias_table,file = "/data/zhangh24/SEIR/result/simulation/bias_table.csv")
write.csv(cover_table,file = "/data/zhangh24/SEIR/result/simulation/cover_table.csv")
