#Load survival table
df1 <- read.csv("~/R-tables/KM_events8.csv", header = TRUE, check.names = TRUE)

#Enter cohort sizes and treatment groups
cohorts=c(5,5,5,0)

#Subsets survival table
#Extracts event data
library(reshape2)
#Arm1
arm1 = cbind(df1$Time,df1[,2:(cohorts[1]+1)])
arm1_melt <- melt(arm1, id=c("df1$Time"))
colnames(arm1_melt) <- c("Time", "ID", "Event")
arm1_melt <- arm1_melt[with(arm1_melt, order(-Event)), ]
arm1_event <- subset(arm1_melt, Event ==0 | Event ==1) 
Tx <- rep(1, cohorts[1])
arm1_event <- cbind(arm1_event, Tx)
#Arm2
arm2 = cbind(df1$Time,df1[,(cohorts[1]+2):(cohorts[1] + cohorts[2]+1)])
arm2_melt <- melt(arm2, id=c("df1$Time"))
colnames(arm2_melt) <- c("Time", "ID", "Event")
arm2_melt <- arm2_melt[with(arm2_melt, order(-Event)), ]
arm2_event <- subset(arm2_melt, Event ==0 | Event ==1)
Tx <- rep(2, cohorts[2])
arm2_event <- cbind(arm2_event, Tx)
#Arm3
if (cohorts[3] != 0) {
  arm3 = cbind(df1$Time,df1[,(cohorts[1] + cohorts[1]+2):(cohorts[1] + cohorts[2] + cohorts[3] + 1)])
  arm3_melt <- melt(arm3, id=c("df1$Time"))
  colnames(arm3_melt) <- c("Time", "ID", "Event")
  arm3_melt <- arm3_melt[with(arm3_melt, order(-Event)), ]
  arm3_event <- subset(arm3_melt, Event ==0 | Event ==1)
  Tx <- rep(3, cohorts[3])
  arm3_event <- cbind(arm3_event, Tx)
}
#Arm4
if (cohorts[4] != 0) {
  arm4 = cbind(df1$Time,df1[,(cohorts[1] + cohorts[2] + cohorts[3] + 2):(cohorts[1] + cohorts[2] + cohorts[3] + cohorts[4] + 1)])
  arm4_melt <- melt(arm4, id=c("df1$Time"))
  colnames(arm3_melt) <- c("Time", "ID", "Event")
  arm4_melt <- arm4_melt[with(arm4_melt, order(-Event)), ]
  arm4_event <- subset(arm4_melt, Event ==0 | Event ==1)
  Tx <- rep(4, cohorts[4])
  arm4_event <- cbind(arm4_event, Tx)
}
#Concatenates event data
if (cohorts[3]==0){
  SurvivalData <- rbind(arm1_event, arm2_event)
}
if (cohorts[3]!=0 & cohorts[4]==0){
  SurvivalData <- rbind(arm1_event, arm2_event, arm3_event)
}
if (cohorts[4]!=0){
  SurvivalData <- rbind(arm1_event, arm2_event, arm3_event, arm4_event)
}

colnames(SurvivalData) <- c("Time", "Number", "Event", "Treatment")
#####################################################################################################################################################################################################################################################
#Make survival plot
library(OIsurv)
library(ggplot2)
library(GGally)

attach(SurvivalData)
mySurv1 <- Surv(Time, Event)
myFit1 <- survfit(mySurv1 ~ Treatment)
ggsurv(myFit1, back.white=TRUE) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 1.1)) + #Adjusts axes scales
  guides(linetype = FALSE) +
  scale_colour_discrete( #Labels legends
    name   = 'Treatment',
    breaks = c(1,2,3),
    labels = c('Mock', 'LPS', 'CpG')
)
detach(SurvivalData)
##########################################################################################################################################
#Survival analysis
summary(myFit)
myFit$surv # outputs the Kaplan-Meier estimate at each t_i
myFit$time # t_i
myFit$n.risk # Y_i
myFit$n.event # d_i
myFit$std.err # standard error of the K-M estimate at t_i
myFit$lower # lower pointwise estimates (alternatively, $upper)

#Calculates KM-statistics
attach(SurvivalData)
survdiff(formula = Surv(Time, Event) ~ Treatment)
detach(SurvivalData)

#Output matrix of p-values for KM analysis

##########################################################################################################################################

library(OIsurv)
#Create survival objects
  #Arm1
    attach(arm1_event)
    mySurv1 <- Surv(Time, Event)
    (myFit1 <- survfit(mySurv1 ~ 1))
    detach(arm1_event)
  #Arm2
    attach(arm2_event)
    mySurv2 <- Surv(Time, Event)
    (myFit2 <- survfit(mySurv2 ~ 1))
    detach(arm2_event)
  #Arm3
    if (cohorts[3] != 0) {
      attach(arm3_event)
      mySurv3<- Surv(Time, Event)
      (myFit3 <- survfit(mySurv3 ~ 1))
      detach(arm3_event)
    }
  #Arm4
    if (cohorts[4] != 0) {
      attach(arm4_event)
      mySurv4<- Surv(Time, Event)
      (myFit4 <- survfit(mySurv4 ~ 1))
      detach(arm4_event)
    } 

    