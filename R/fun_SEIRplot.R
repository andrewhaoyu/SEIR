## SEIR model plots for the six panels
## five periods: Jan 1-9 (index 1-9), Jan 10-22 (index 10-22), Jan 25-Feb 1 (index 23-32), Feb 2-16 (index 33-47), Feb 17- (index 48-60)
#' @param pars_estimate           a vetor of parameters: c(b12, b3, b3, b5, r12, delta3, delta4, delta5)
#' @param init_settings           a list of initial values and known parameters
#' @param Dp                      presymptomatic infectious period
#' @param Di                      symptomatic infectious period
#' @param De                      latent period
#' @param Dq                      duration from illness onset to hospitalization
#' @param Dh                      hospitalization period                   
#' @param alpha                   ratio of the transmission rate of unascertained over ascertained case
#' @param N                       population size
#' @param flowN_vec               daily inbound and outbound size during five periods (n)
#' @param init_states             initial c(S, E, P, Is, A, H, R)
#' @param days_to_fit             the days to fit the model
#' @param b                       transmission rate of ascertained cases
#' @param r                       ascertainment rate 
#' @param num_periods             number of periods to simulate
#################################################################################################
SEIRplot <- function(pars_estimate, file_name, init_settings, panel_B_R_ylim=4,
                     stage_intervals=stage_intervals,
                     all.date = all.date) {
  #total days to predict
  n.days <- length(init_settings$daily_new_case_all)
  stages <- length(stage_intervals)
  start.vec <- rep(0,stages)
  end.vec <- rep(0,stages)
  for(k in 1:stages){
    start.vec[k] <- stage_intervals[[k]][1]
    end.vec[k] <- stage_intervals[[k]][2]
  }
  #days cut = total days/stages
  init_settings$days_to_fit <- 1:n.days
  
  library(vioplot)
 my.color = c("#BC3C29FF","#0072B5FF", "#E18727FF", "#7876B1FF", "#FFDC91FF","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
  ##
  onset_obs_all <- as.numeric(init_settings$daily_new_case_all)
  ptime <- 1:length(onset_obs_all)
  mydate <- format(all.date,"%b %d")
  stage.label <- rep("c",stages)
  for(l in 1:stages){
    stage.label[l] = paste0(format(all.date[start.vec[l]],"%b %d")," - ",format(all.date[end.vec[l]],"%d"))
  }
  #
  pdf(paste0("../output/Figure_", file_name, ".pdf"), width = 9, height = 10)
  par(mar = c(4, 5, 2.5, 1))
  layout(matrix(c(1:3), byrow = T, nrow = 3))
  
  l <- stages
  estSEAIP_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = n.stage)[, c("S","E","I", "A", "P","Onset_expect")])
   
  estS_mat <- estSEAIP_mat[ptime,]
  estE_mat <- estSEAIP_mat[ptime + length(ptime),]
  estI_mat <- estSEAIP_mat[ptime + length(ptime)*2, ]
  estA_mat <- estSEAIP_mat[ptime + length(ptime)*3, ]
  estP_mat <- estSEAIP_mat[ptime + length(ptime) * 4, ]
  estN_mat <- estSEAIP_mat[ptime + length(ptime) * 5, ]
  
  #estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = l)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start A
  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1,end.vec), labels = mydate[c(1,end.vec)])
  #
  
  abline(v = end.vec, lty = 3, lwd = 2, col = "darkgrey")
  text(end.vec, par()$usr[4], labels = mydate[end.vec], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:end.vec[l]], rev(ptime[1:end.vec[l]])), c(estN_up[1:end.vec[l]], rev(estN_low[1:end.vec[l]])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:end.vec[l])], rev(ptime[-c(1:end.vec[l])])), c(estN_up[-c(1:end.vec[l])], rev(estN_low[-c(1:end.vec[l])])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:end.vec[l]], estN_mean[1:end.vec[l]], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:end.vec[l])], estN_mean[-c(1:end.vec[l])], col = "#0072B5FF", pch = 17, cex = 0.8)
  points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  #text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)
  
  
  ##########################################   Panel B  ##########################################################
  estRt_mat <- apply(pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_settings))
  estRt_mat <- t(estRt_mat)
  ##
  rt_mean <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) mean(x)), 2))
  rt_low <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.025)), 2))
  rt_up <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.975)), 2))
  #
  vioplot(estRt_mat,names = NA, ylim = c(0, panel_B_R_ylim), col = my.color[1:stages], ylab = "", xlab = "")
  #vioplot(estRt_mat[, 1], estRt_mat[, 2], estRt_mat[, 3], estRt_mat[, 4], estRt_mat[, 5], names = NA, ylim = c(0, panel_B_R_ylim), col = c("#BC3C29FF","#0072B5FF", "#E18727FF", "#7876B1FF", "#FFDC91FF"), ylab = "", xlab = "")
  mtext("Outbreak period (2020)", side = 1, line  = 3, cex = 1.01)
  mtext(expression("R"["0"]), side = 2, line = 3, cex = 1.01)
  axis(1, at = 1:n.stage, tick = F, labels = stage.label)
  abline(h = 1, lwd = 2, lty = 3, col = "red")
  #
  for(l in 1:stages){
    if(min(estRt_mat[, l]>2)){
      text(l, min(estRt_mat[, l]) - 0.3, labels = rt_mean[l])
      text(l, min(estRt_mat[, l]) - 0.7, labels = paste("(", rt_low[l], "-", rt_up[l], ")", sep = ""))
    }else{
      text(l, max(estRt_mat[, l]) + 0.3, labels = rt_mean[l])
      text(l, max(estRt_mat[, l]) + 0.7, labels = paste("(", rt_low[l], "-", rt_up[l], ")", sep = ""))
      }
    
  }
  # text(1, min(estRt_mat[, 1]) - 0.2, labels = rt_mean[1])
  # text(1, min(estRt_mat[, 1]) - 0.45, labels = paste("(", rt_low[1], "-", rt_up[1], ")", sep = ""))
  # text(2, min(estRt_mat[, 2]) - 0.2, labels = rt_mean[2])
  # text(2, min(estRt_mat[, 2]) - 0.45, labels = paste("(", rt_low[2], "-", rt_up[2], ")", sep = ""))
  # text(3, max(estRt_mat[, 3]) + 0.4, labels = rt_mean[3])
  # text(3, max(estRt_mat[, 3]) + 0.15, labels = paste("(", rt_low[3], "-", rt_up[3], ")", sep = ""))
  # text(4, max(estRt_mat[, 4]) + 0.4, labels = rt_mean[4])
  # text(4, max(estRt_mat[, 4]) + 0.15, labels = paste("(", rt_low[4], "-", rt_up[4], ")", sep = ""))
  # text(5, max(estRt_mat[, 5]) + 0.4, labels = rt_mean[5])
  # text(5, max(estRt_mat[, 5]) + 0.15, labels = paste("(", rt_low[5], "-", rt_up[5], ")", sep = ""))
  # #text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "B", xpd = T, cex = 2)
  
  
  
  ##########################################   Panel C-E  ##########################################################
  # letters <- c("C","D","E")
  # for(l in (stages-1):2){
  #   estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = l)[, "Onset_expect"])
  #   estN_mean <- round(apply(estN_mat, 1, mean), 0)
  #   estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  #   estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  #   # start A
  #   plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  #   mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  #   mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  #   axis(1, at = seq(1, n.days, cut.days), labels = mydate[seq(1, n.days, cut.days)])
  #   #
  #   
  #   abline(v = end.vec, lty = 3, lwd = 2, col = "darkgrey")
  #   text(end.vec, par()$usr[4], labels = mydate[end.vec], col = "darkgrey", pos = 3, xpd = T)
  #   #
  #   polygon(c(ptime[1:end.vec[l]], rev(ptime[1:end.vec[l]])), c(estN_up[1:end.vec[l]], rev(estN_low[1:end.vec[l]])), col = "#F39B7FB2", border = NA)
  #   polygon(c(ptime[-c(1:end.vec[l])], rev(ptime[-c(1:end.vec[l])])), c(estN_up[-c(1:end.vec[l])], rev(estN_low[-c(1:end.vec[l])])), col = "#4DBBD5B2", border = NA)
  #   #
  #   points(ptime[1:end.vec[l]], estN_mean[1:end.vec[l]], col = "#BC3C29FF", pch = 16, cex = 0.8)
  #   points(ptime[-c(1:end.vec[l])], estN_mean[-c(1:end.vec[l])], col = "#0072B5FF", pch = 17, cex = 0.8)
  #   points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  #   #
  #   legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  #   #text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)
  #   
  #   
  #   
  # }
  
  
  
  ##########################################   Panel F  ##########################################################
  
  
  estI_mean <- apply(estI_mat, 1, mean)
  estA_mean <- apply(estA_mat, 1, mean)
  estP_mean <- apply(estP_mat, 1, mean)
  estAIP_dat <- rbind(estI_mean, estA_mean, estP_mean)
  barpos <- barplot(estAIP_dat, col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"), xlab = "", ylab = "", border = "NA")
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of active infectious cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = barpos[ c(1,end.vec,length(all.date))], labels = mydate[ c(1,end.vec,length(all.date))])
  legend("topleft", legend = c("Presymptomatic (P)", "Unascertained (A)", "Ascertained (I)"), fill = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF"), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "F", xpd = T, cex = 2)
  ##figure_F finished
  dev.off()
  prevalence <- rowMeans((N-estS_mat)/N)
  prevalence_low <- apply(estS_mat,1,function(x){quantile((N-x)/N,0.025)})
  prevalence_high <- apply(estS_mat,1,function(x){quantile((N-x)/N,0.975)})
  
  library(ggplot2)
  if(statename[i1]%in%CDC$Statename){
    data = data.frame(date= all.date,Prevalance=100*prevalence,Prevalence_low = 100*prevalence_low,Prevalence_high= 100*prevalence_high)
    
    CDC_filter <- CDC %>% filter(Statename==statename[i1]) %>% 
      mutate(date = as.Date(Infection_date),format="%y-%m-%d") %>% 
      select(Statename,Prevalance,Prevalance_low,Prevalance_high,date)
    idx <- which(data$date%in%CDC_filter$date)
    data_select = data %>% select(date,Prevalance) %>% 
      rename(Prevalance_pred = Prevalance)
    
    CDC_match = left_join(CDC_filter,data_select,by="date")
    
    least_square = sum(CDC_match$Prevalance-CDC_match$Prevalance_pred)^2/(nrow(CDC_match))
    p <- ggplot(data,aes(x=date))+geom_line(aes(y = Prevalance))+
      geom_ribbon(aes(ymin=Prevalence_low,ymax=Prevalence_high),alpha = 0.2)+
      geom_point(data= CDC_filter,aes(x=date,y = Prevalance))+
      geom_errorbar(data=CDC_filter,aes(ymin = Prevalance_low,ymax=Prevalance_high))+
      theme_Publication()+
      ggtitle(paste0("Prevalance estimate in ",paste0(statename[i1])," (least square = ",least_square))
    png(paste0("../output/Prevalance_", file_name, ".png"), width = 10, height = 10,res=300,units="in")
    print(p)
    dev.off()
    
  }else{
    data = data.frame(date= all.date,Prevalance=100*prevalence,Prevalence_low = 100*prevalence_low,Prevalence_high= 100*prevalence_high)
    
    CDC_filter <- CDC %>% filter(Statename==statename[i1]) %>% 
      mutate(date = as.Date(Infection_date),format="%y-%m-%d") %>% 
      select(Statename,Prevalance,Prevalance_low,Prevalance_high,date)
    p <- ggplot(data,aes(x=date))+geom_line(aes(y = Prevalance))+
      geom_ribbon(aes(ymin=Prevalence_low,ymax=Prevalence_high),alpha = 0.2)+
      #geom_point(data= CDC_filter,aes(x=date,y = Prevalance))+
      #geom_errorbar(data=CDC_filter,aes(ymin = Prevalance_low,ymax=Prevalance_high))+
      theme_Publication()+
      ylim(0,100)+
      ggtitle(paste0("Prevalance estimate in ",paste0(statename[i1])))
    png(paste0("../output/Prevalance_", file_name, ".png"), width = 10, height = 10,res=300,units="in")
    print(p)
    dev.off()
    
  }
  
}