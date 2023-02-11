rm(list=ls())

#############################
# Source Files
#############################

library(did)
source("SPC_Ft.R")

#############################
# Dataset
# Zika Virus in Brazil
# The dataset is used in Universal Difference-in-Differences for Causal Inference in Epidemiology 
# (Tchetgen Tchetgen, Park, Richardson, 2023) https://arxiv.org/abs/2302.00840
# The source of the dataset are given below:
# Pre- and Post-treatment Outcomes, Treatment, and log population: 
#       zika_Table2.tab in  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY
# log population density and proportion of female
#       https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados
#############################

LongData <- read.csv("Zika_Brazil.csv")
# Y    = Outcome = birth rate
# Trt  = Indicator of whether a municipality belongs to Pernambuco (PE) (Trt=1) or belongs to Rio Grande do Sul (RS) (Trt=0)
# Time = 2014 (Time=0) or 2016 (Time=1)
# X_1  = log population
# X_2  = log population density
# X_3  = proportion of female
# ID   = municipality ID

# We only use Y,Trt,Time,ID in Tchetgen Tchetgen, Park, Richardson (2023) Single Proxy Control
# LongData: used for DiD

LongData <- LongData[,c("Y","Trt","Time","ID")]
Y0     <- LongData$Y[LongData$Time==0]
Y1     <- LongData$Y[LongData$Time==1]
Trt    <- LongData$Trt[LongData$Time==0]

# (Short) Data: used for Crude & Single Proxy Control
Data <- data.frame(cbind(Y0,Y1,Trt))
mY0  <- mean(Data$Y0)   # Pre-treatment outcome mean 
Data$Y0 <- Data$Y0-mY0  # Shift-mean for stability
Data$Y1 <- Data$Y1-mY0  # Shift-mean for stability


#############################
# DiD - Parallel Trends
#############################

PT.OR.Mar <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~1,
                    data=LongData,
                    bstrap = F,
                    est_method = "reg")
PT.PS.Mar <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~1,
                    data=LongData,
                    bstrap = F,
                    est_method = "ipw")
PT.DR.Mar <- att_gt(yname="Y",
                    tname="Time",
                    idname="ID",
                    gname="Trt",
                    xformla=~1,
                    data=LongData,
                    bstrap = F,
                    est_method = "dr")

PP <- function(ll){
  c(ll$att, ll$att-qnorm(0.975)*ll$se, ll$att+qnorm(0.975)*ll$se)
}

Result.PT <- cbind( c(PP(PT.OR.Mar)), 
                    c(PP(PT.PS.Mar)),
                    c(PP(PT.DR.Mar))) 

Result.PT 
# DID-Parallel Trends
# Est -1.1912000
# LB  -1.5067090
# UB  -0.8756909

#############################
# Crude 
#############################

mean(Data$Y1[Data$Trt==1]) - mean(Data$Y1[Data$Trt==0])
t.test(Data$Y1[Data$Trt==1],
       Data$Y1[Data$Trt==0])
# Crude
# Est 3.384296
# LB  2.952956
# UB  3.815636

#############################
# Single Proxy Control
#############################

RESULT <- matrix(0,2,6)  # Result summary
W.Origin <- Data$Y0      # Original W
Y.Origin <- Data$Y1      # Original Y 

for(ii in c(1,2)){
  
  num.bin    <- as.numeric(ii==1)*5 + as.numeric(ii==2)*10
  Bin0.Cut   <- quantile(W.Origin,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
  Bin1.Cut   <- quantile(Y.Origin,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
  Disc.Y0    <- apply(matrix(W.Origin,length(W.Origin),1),1,function(v){sum(v>Bin0.Cut)})
  Disc.Y1    <- apply(matrix(Y.Origin,length(Y.Origin),1),1,function(v){sum(v>Bin1.Cut)})
  
  Ymat.D <- Wmat.D <- matrix(0,length(Disc.Y0),num.bin)
  for(jj in 1:num.bin){
    Wmat.D[Disc.Y0==jj-1,jj] <- 1
    Ymat.D[Disc.Y1==jj-1,jj] <- 1
  }
  
  Ymat.D <- cbind(Ymat.D)  # Discretized Y ; used in EPS estimation 
  Wmat.D <- cbind(Wmat.D)  # Discretized W ; used in ORB estimation
  Ymat.D1 <- Ymat.D        
  Ymat.D1[,num.bin] <- 1   # make the last bin as the reference for EPS estimation
  
  FACTOR <- 2
  
  Bin0.2.Cut   <- quantile(W.Origin,seq(0,1,length=FACTOR*num.bin+1))[-c(1,1+FACTOR*num.bin)]
  Bin1.2.Cut   <- quantile(Y.Origin,seq(0,1,length=FACTOR*num.bin+1))[-c(1,1+FACTOR*num.bin)]
  Disc.2.Y0    <- apply(matrix(W.Origin,length(W.Origin),1),1,function(v){sum(v>Bin0.2.Cut)})
  Disc.2.Y1    <- apply(matrix(Y.Origin,length(Y.Origin),1),1,function(v){sum(v>Bin1.2.Cut)})
  
  Ymat.2.D <- Wmat.2.D <- matrix(0,length(Disc.2.Y0),FACTOR*num.bin)
  for(jj in 1:(FACTOR*num.bin)){
    Wmat.2.D[Disc.2.Y0==jj-1,jj] <- 1
    Ymat.2.D[Disc.2.Y1==jj-1,jj] <- 1
  }
  
  Ymat.2.D <- cbind(Ymat.2.D)  # Discretized W ; used in EPS estimation 
  Wmat.2.D <- cbind(Wmat.2.D)  # Discretized Y ; used in ORB estimation
  
  Result.PS.D <- SPC.PS(Y.Origin,Trt,
                        Ymat=Ymat.D1,
                        Wmat=Wmat.2.D,
                        penalty=10^(-6))
  
  Result.OR.D <- SPC.OR(Y.Origin,Trt,
                        Ymat=Ymat.2.D,
                        Wmat=Wmat.D)
  
  Result.DR.D <- SPC.DR(Y.Origin,Trt,
                        Ymat.PS=Ymat.D1,
                        Wmat.PS=Wmat.2.D,
                        Ymat.OR=Ymat.2.D,
                        Wmat.OR=Wmat.D,
                        penalty=10^(-6)) 
  
  RESULT[ii,] <- c( Result.PS.D$Est,Result.PS.D$SE,
                    Result.OR.D$Est,Result.OR.D$SE,
                    Result.DR.D$Est,Result.DR.D$SE )
  
  
}


# Single Proxy Control
#            EPS-Est    EPS-SE   ORB-Est   ORB-SE    DR-Est     DR-SE
# 5-bin    -1.694627 0.4542382 -2.470063 0.701895 -1.857519 0.5133388
# 10-bin   -1.941184 0.4160650 -2.601973 1.028506 -1.816989 0.5760705

C1 <- c("\\multicolumn{1}{|c|}{\\multirow{2}{*}{5}}",
        "\\multicolumn{1}{|c|}{}",
        "\\multicolumn{1}{|c|}{\\multirow{2}{*}{10}}",
        "\\multicolumn{1}{|c|}{}")
C2 <- c("\\multicolumn{1}{c|}{Estimate}",
        "\\multicolumn{1}{c|}{95$\\%$ CI}",
        "\\multicolumn{1}{c|}{Estimate}",
        "\\multicolumn{1}{c|}{95$\\%$ CI}")
FP <- function(tt){
  c( sprintf("\\multicolumn{1}{c|}{%0.3f}",RESULT[1,tt]),
     sprintf("\\multicolumn{1}{c|}{(%0.3f, %0.3f)}", 
             RESULT[1,tt]-RESULT[1,tt+1]*qnorm(0.975),
             RESULT[1,tt]+RESULT[1,tt+1]*qnorm(0.975) ),
     sprintf("\\multicolumn{1}{c|}{%0.3f}",RESULT[2,tt]),
     sprintf("\\multicolumn{1}{c|}{(%0.3f, %0.3f)}", 
             RESULT[2,tt]-RESULT[2,tt+1]*qnorm(0.975),
             RESULT[2,tt]+RESULT[2,tt+1]*qnorm(0.975) ) )
}

C3 <- FP(1)
C4 <- FP(3)
C5 <- FP(5)
C6 <- c( "\\\\ \\cline{2-5}" ,
         "\\\\ \\hline",
         "\\\\ \\cline{2-5}" ,
         "\\\\ \\hline")
DDD <- data.frame( cbind(C1," & ",C2, " & ", C3, " & ", C4, " & ", C5, C6) )
print(DDD,row.names=F)

#############################
# Sensitivity Analysis
# Implemented in SPC_Sensitivity_Cluster.R
# Recommend to use clusters
#############################

W.Origin <- Data$Y0
Y.Origin <- Data$Y1

num.bin    <- 5
Bin0.Cut   <- quantile(W.Origin,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
Bin1.Cut   <- quantile(Y.Origin,seq(0,1,length=num.bin+1))[-c(1,1+num.bin)]
Disc.Y0    <- apply(matrix(W.Origin,length(W.Origin),1),1,function(v){sum(v>Bin0.Cut)})
Disc.Y1    <- apply(matrix(Y.Origin,length(Y.Origin),1),1,function(v){sum(v>Bin1.Cut)})

Ymat.D <- Wmat.D <- matrix(0,length(Disc.Y0),num.bin)
for(jj in 1:num.bin){
  Wmat.D[Disc.Y0==jj-1,jj] <- 1
  Ymat.D[Disc.Y1==jj-1,jj] <- 1
}

Ymat.D <- cbind(Ymat.D)  # Discretized Y ; used in EPS estimation 
Wmat.D <- cbind(Wmat.D)  # Discretized W ; used in ORB estimation
Ymat.D1 <- Ymat.D        
Ymat.D1[,num.bin] <- 1   # make the last bin as the reference for EPS estimation

FACTOR <- 2

Bin0.2.Cut   <- quantile(W.Origin,seq(0,1,length=FACTOR*num.bin+1))[-c(1,1+FACTOR*num.bin)]
Bin1.2.Cut   <- quantile(Y.Origin,seq(0,1,length=FACTOR*num.bin+1))[-c(1,1+FACTOR*num.bin)]
Disc.2.Y0    <- apply(matrix(W.Origin,length(W.Origin),1),1,function(v){sum(v>Bin0.2.Cut)})
Disc.2.Y1    <- apply(matrix(Y.Origin,length(Y.Origin),1),1,function(v){sum(v>Bin1.2.Cut)})

Ymat.2.D <- Wmat.2.D <- matrix(0,length(Disc.2.Y0),FACTOR*num.bin)
for(jj in 1:(FACTOR*num.bin)){
  Wmat.2.D[Disc.2.Y0==jj-1,jj] <- 1
  Ymat.2.D[Disc.2.Y1==jj-1,jj] <- 1
}

Ymat.2.D <- cbind(Ymat.2.D)  # Discretized W ; used in EPS estimation 
Wmat.2.D <- cbind(Wmat.2.D)  # Discretized Y ; used in ORB estimation

sens.para <- c(0,0.439,0.614,0.8043)   # the values are found from grid-search

Sens.Result <- list()

Sens.Result[[1]] <- SPC.PS.Sensitivity(Y.Origin = Y.Origin,
                                       Trt = Trt,
                                       Ymat = Ymat.D1,
                                       Wmat = Wmat.2.D,
                                       Purturb = W.Origin*sens.para[1],
                                       penalty=10^(-6))
Sens.Result[[2]] <- SPC.PS.Sensitivity(Y.Origin = Y.Origin,
                                       Trt = Trt,
                                       Ymat = Ymat.D1,
                                       Wmat = Wmat.2.D,
                                       Purturb = W.Origin*sens.para[2],
                                       penalty=10^(-6))
Sens.Result[[3]] <- SPC.PS.Sensitivity(Y.Origin = Y.Origin,
                                       Trt = Trt,
                                       Ymat = Ymat.D1,
                                       Wmat = Wmat.2.D,
                                       Purturb = W.Origin*sens.para[3],
                                       penalty=10^(-6))
Sens.Result[[4]] <- SPC.PS.Sensitivity(Y.Origin = Y.Origin,
                                       Trt = Trt,
                                       Ymat = Ymat.D1,
                                       Wmat = Wmat.2.D,
                                       Purturb = W.Origin*sens.para[4],
                                       penalty=10^(-6))

RR <- rbind( sens.para,
             sapply(1:4,function(tt){Sens.Result[[tt]]$Est}),
             sapply(1:4,function(tt){Sens.Result[[tt]]$Est-qnorm(0.975)*Sens.Result[[tt]]$SE}),
             sapply(1:4,function(tt){Sens.Result[[tt]]$Est+qnorm(0.975)*Sens.Result[[tt]]$SE}))

RR
#                   COCA      COCA(NULL)   COCA(Positive)     COCA(Crude)
# sens.para    0.0000000     0.439000000      0.614000000      0.8043000
# Estimate    -1.6946269    -0.661963017      0.007796038      0.9267923
# LB          -2.5849174    -1.328821835     -0.962489263     -1.1183557
# UB          -0.8043364     0.004895802      0.978081338      2.9719403

print(data.frame( cbind( c("Sensitivity Parameter & ","Estimate &",
                           "95\\% CI & "),
                         c( sprintf("%0.3f &",RR[1,1]),
                            sprintf("%0.3f &",RR[2,1]),
                            sprintf("(%0.3f,  %0.3f) &",RR[3,1],RR[4,1]) ),
                         c( sprintf("%0.3f &",RR[1,2]),
                            sprintf("%0.3f &",RR[2,2]),
                            sprintf("(%0.3f,  %0.3f) &",RR[3,2],RR[4,2]) ),
                         c( sprintf("%0.3f &",RR[1,3]),
                            sprintf("%0.3f &",RR[2,3]),
                            sprintf("(%0.3f,  %0.3f) &",RR[3,3],RR[4,3]) ),
                         c( sprintf("%0.3f \\\\ \\hline",RR[1,4]),
                            sprintf("%0.3f \\\\ \\hline",RR[2,4]),
                            sprintf("(%0.3f,  %0.3f) \\\\ \\hline",RR[3,4],RR[4,4]) ) )),row.names=F)


 