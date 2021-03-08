b_vec = c(0.7,0.3)
phi_vec = c(0,0.2)
ascertainment_mat = matrix(c(c(0.05,0.07),
                             (c(0.05,0.12)+0.2)),ncol=2)




n.rep =1000
#2 ascertainment settings
#2 overdispersion settings
#2 methods
total = n.rep
n.stage= 2
est = matrix(0,total,n.stage*2+1)
est_cover_mat = matrix(0,total,n.stage*2+1)
est_low_mat = matrix(0,total,n.stage*2+1)
est_high_mat = matrix(0,total,n.stage*2+1)
pars_est = matrix(0,total,n.stage*2+1)
pars_est_cover = matrix(0,total,n.stage*2+1)
i2_vec = rep(0,total)
method_vec = rep(0,total)
temp = 1
i2 = 1
i3 = 2

phi = 0.2
r_vec = ascertainment_mat[,i3]

pars = c(b_vec,r_vec,phi)
for(i2 in 1:2){
  temp = 1
  for(i1 in 1:n.rep){
    load(paste0("/data/zhangh24/SEIR/result/simulation/two_stage_",i1,"_",i2,"_",i3," .rdata"))
    est[temp,] = est_result[[1]]
    est_low_mat[temp,] = est_result[[2]]
    est_high_mat[temp,] = est_result[[3]]
    est_low = est_result[[2]]
    est_high = est_result[[3]]
    est_cover_mat[temp,] = (pars>=est_low)*(pars<=est_high)
    
    pars_est[temp,] = est_result[[7]]
    # pars_est_low = est_result[[8]]
    # pars_est_high = est_result[[9]]
    # pars_est_cover[temp,] = (pars>=pars_est_low)*(pars<=pars_est_high)
    i2_vec[temp] = i2
    method_vec[temp] = "nb"
    temp = temp+1
    #load(paste0("/data/zhangh24/SEIR/result/simulation/seir_pos_result_",i1,"_",i2,"_",i3,".rdata"))
    
  }
  est_mean = colMeans(est)
  est_cover = colMeans(est_cover_mat)
  est_low = colMeans(est_low_mat)
  est_high = colMeans(est_high_mat)  
  result = data.frame(est_mean,
                      est_cover,
                      est_low,
                      est_high)
  rownames(result) = c("b1","b2","r1","r2","phi")
  save(result,file = paste0("/data/zhangh24/SEIR/result/simulation/sum_two_stage_",i2,"_",i3,".rdata"))
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

