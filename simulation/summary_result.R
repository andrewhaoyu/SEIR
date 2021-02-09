b_vec = c(0.7,0.3,0.3,0.3,0.48,0.48)
r_vec = c(0.2,0.25,0.3,0.35,0.4,0.30)

n_stage = 6
n.rep = 1000
b_est = matrix(0,n.rep,6)
b_cover = matrix(0,n.rep,6)
r_est = matrix(0,n.rep,6)
r_cover = matrix(0,n.rep,6)
for(i1 in 1:n.rep){
  load(paste0("/data/zhangh24/SEIR/result/simulation/seir_result_",i1,".rdata"))
  b_est[i1,] = est_result[[1]][1:n_stage]
  r_est[i1,] = est_result[[1]][(n_stage+1):(2*n_stage)]
  b_cover[i1,] = (est_result[[2]][1:n_stage]<=b_vec)*
    (est_result[[3]][1:n_stage]>=b_vec)
  r_cover[i1,] =  (est_result[[2]][(n_stage+1):(2*n_stage)]<=r_vec)*
    (est_result[[3]][(n_stage+1):(2*n_stage)]>=r_vec)
  
}
colMeans(b_est)
colMeans(r_est)
colMeans(b_cover)
colMeans(r_cover)
colMeans(r_est)-r_vec
