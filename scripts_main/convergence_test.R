run_id = paste0("080221_",i1,"_",i2,"_",i3)
code_root = "/data/zhangh24/SEIR/"
#code_root = "/n/holystore01/LABS/xlin/Lab/hzhang/SEIR/"
#code_root = "/dcl01/chatterj/data/hzhang1/temp/SEIR/"
setwd(paste0(code_root, "scripts_main"))


i1_vec = rep(0,1000)

i3_vec = rep(0,1000)
df = rep(0,1000)
temp = 1
for(i1 in 1:2){
  for(i3 in 1:10){
    run_id = paste0("080221_",i1,"_",i2,"_",i3)
    load(paste0("../output/mcmc_out",run_id,".rdata"))
    diagnosis = gelmanDiagnostics(mh_out)
    df[temp] = diagnosis$mpsrf
    i1_vec[temp] = i1
    i3_vec[temp] = i3
    temp = temp + 1
  }
}
i1_vec = i1_vec[1:(temp-1)]
i3_vec = i3_vec[1:(temp-1)]
df = df[1:(temp-1)]
result = cbind(i1_vec,i3_vec,df)
result

