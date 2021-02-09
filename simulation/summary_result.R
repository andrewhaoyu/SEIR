b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
r_vec = c(0.05,0.07,0.09,0.10,0.11,0.12)

n_stage = 6
n.rep = 1000
b_est = matrix(0,n.rep,6)
b_est_low = matrix(0,n.rep,6)
b_est_high = matrix(0,n.rep,6)
b_cover = matrix(0,n.rep,6)
Rt_est = matrix(0,n.rep,6)
Rt_cover = matrix(0,n.rep,6)
r_est = matrix(0,n.rep,6)
r_cover = matrix(0,n.rep,6)
pars = c(b_vec,r_vec)

Rt_true = estimate_R(pars,
                     Di,
                     Dp,
                     Dq_vec,
                     N,
                     flowN_vec,
                     n_stage)

for(i1 in 1:n.rep){
  load(paste0("/data/zhangh24/SEIR/result/simulation/seir_result_",i1,".rdata"))
  b_est[i1,] = est_result[[1]][1:n_stage]
  b_est_low[i1,] = est_result[[2]][1:n_stage]
  b_est_high[i1,] = est_result[[3]][1:n_stage]
  r_est[i1,] = est_result[[1]][(n_stage+1):(2*n_stage)]
  b_cover[i1,] = (est_result[[2]][1:n_stage]<=b_vec)*
    (est_result[[3]][1:n_stage]>=b_vec)
  r_cover[i1,] =  (est_result[[2]][(n_stage+1):(2*n_stage)]<=r_vec)*
    (est_result[[3]][(n_stage+1):(2*n_stage)]>=r_vec)
  Rt_est[i1,] = rowMeans(est_result[[4]])
  Rt_low = apply(est_result[[4]],1,function(x){quantile(x,0.025)})
  Rt_high = apply(est_result[[4]],1,function(x){quantile(x,0.975)})
  Rt_cover[i1,] = (Rt_low<=Rt_true)*(Rt_high>=Rt_true)
}
round(colMeans(b_est)-b_vec,3)
colMeans(r_est)
colMeans(b_cover)
colMeans(r_cover)
round(colMeans(r_est)-r_vec,3)
round(colMeans(Rt_est)-Rt_true,3)
colMeans(Rt_cover)
round(Rt_true,2)
  
  
  
  
  
  
  
  
  
  

