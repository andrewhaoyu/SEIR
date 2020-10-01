#process HWF data
#difference between PCR test result and sympotom onset
mean.diff <- rep(0,9)

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Feb_1_15.csv",header=T)
mean.diff[1] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Feb_16_29.csv",header=T)
mean.diff[2] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Mar_1_15.csv",header=T)
mean.diff[3] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Mar_16_31.csv",header=T)
mean.diff[4] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])


HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Apr_1_15.csv",header=T)
mean.diff[5] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])


HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_Apr_16_30.csv",header=T)
mean.diff[6] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])


HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_May_1_15.csv",header=T)
mean.diff[7] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_May_16_31.csv",header=T)
mean.diff[8] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])

HWF <- read.csv("./data/HWF_data/HWF_PCR_Swab_Table_June_1_18.csv",header=T)
mean.diff[9] <- sum(HWF[,1]*HWF[,2])/sum(HWF[,2])


write.csv(mean.diff,file = "./data/HWF_data/mean_diff.csv")
