library(readr)
library(coda)
library(cairoDevice)
file.dir <- "/data/zhangh24/SEIR/output/"
setwd(file.dir)

files <- dir(path=file.dir,pattern = "pars_est_run_")
library(data.table)
library(reshape2)
library(ggplot2)
for(i1 in 1:4){
  for(i3 in 1:9){
    
    result.list <- list()
    temp = 1
    
    for(i2 in 1:3){
      file <- paste0(paste0("pars_est_run_",i1,"_",i2,"_",i3,".txt"))
      if(file%in%files){
        pars_estimate_main = read.table(file=file,header=T)
        n <- nrow(pars_estimate_main)
        pars_estimate_main$x = c(1:n)
        pars_estimate_main$chain = rep(paste0("chain",i2),n)
        result.list[[temp]] <- pars_estimate_main
        temp = temp+1
      }
      
      
    }
    result <- rbindlist(result.list)
    data.long <- melt(result,id.vars=c("x","chain"))
    p <- ggplot(data.long)+
      geom_line(aes(x= x,y = value,color=chain))+
      facet_wrap(~variable)+
      theme_Publication()
    png(filename=paste0("./check_convergence/check_convergence_",i1,"_",i3,".png"),height=8,width = 10,units = "in",res = 300)
    print(p)
    dev.off()
  }
}







