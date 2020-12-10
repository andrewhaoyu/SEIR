#generate plots for ascertainment
library(ggplot2)
library(data.table)
source("./R/theme_Publication.R")
#FL
data <- as.data.frame(fread("./output/result_102520/summary_run_102520_3_1_3.csv",header=T))
n_stage = 10
test = c(0.09846624,0.10552756,0.06723429,0.04204461,0.03397413,0.0464,0.14795901,0.12427998,0.09779907,0.06216735)
logit <- function(x){log(x/(1-x))}
ascertainment = c(0.15,0.13,0.14,0.21,0.38,0.46,0.54, 0.54, 0.58, 0.56)
ascertainment_low = c(0.1,0.05,0.05,0.06,0.12,0.12,0.14,0.12,0.12,0.1)
ascertainment_high = c(0.23,0.24,0.29,0.43,0.75,0.84,0.9,0.94,0.96,0.97)
stage_label <- as.Date(c("2020-03-24", "2020-04-08", "2020-04-23",
                 "2020-05-08","2020-05-23","2020-06-08","2020-06-24",
                 "2020-07-08","2020-07-24","2020-08-11"),format="%Y-%m-%d")

plot.data <- data.frame(stage_label,ascertainment,ascertainment_low,
                        ascertainment_high,test)
p = ggplot(plot.data,aes(x = log(test),y = logit(ascertainment)))+
  geom_point()+
  geom_smooth(method='lm')+
  xlab("Time-period")+
  ylab("Ascertainment (95%CI)")+
  theme_Publication()
  
  #geom_ribbon(aes(ymin=ascertainment_low,ymax=ascertainment_high),alpha = 0.2)+
  print(p)
model <- lm(logit(ascertainment)~test)
summary(model)
  png(filename = "./output/result_102520/FL_ascertainment.png",width = 12,height = 8, units = "in",res=300)
print(p)
dev.off()











test = c(1288.083,2265.733,3132.667,3918.333,5883.625,7481.800,7739.000,9153.333,11146.125,8980.875)/10000
logit <- function(x){log(x/(1-x))}
ascertainment = c(0.37,0.38,0.41,0.42,0.43,0.47,0.46, 0.46, 0.47, 0.42)
#ascertainment_low = c(0.1,0.05,0.05,0.06,0.12,0.12,0.14,0.12,0.12,0.1)
#ascertainment_high = c(0.23,0.24,0.29,0.43,0.75,0.84,0.9,0.94,0.96,0.97)
#test = c(0.397,0.709,1.2223,1.4229,1.139,1.115,1.2455,1.46,1.79581,1.912369)
plot.data <- data.frame(test,ascertainment)
                        #ascertainment_low,
                        #ascertainment_high,test)
p = ggplot(plot.data,aes(x = log(test),y = logit(ascertainment)))+
  geom_point()+
  geom_smooth(method='lm')+
  xlab("Time-period")+
  ylab("Ascertainment (95%CI)")+
  theme_Publication()

#geom_ribbon(aes(ymin=ascertainment_low,ymax=ascertainment_high),alpha = 0.2)+
print(p)
model <- lm(logit(ascertainment)~test)
summary(model)

