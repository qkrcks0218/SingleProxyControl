rm(list=ls())

################################################################################
# Source File
################################################################################

source("0.SPC_Ft.R")
library(did)

################################################################################
# Dataset
# Zika Virus in Brazil
# The dataset is used in 
# Universal Difference-in-Differences for Causal Inference in Epidemiology 
# (Tchetgen Tchetgen, Park, Richardson, 2023) https://arxiv.org/abs/2302.00840
# The source of the dataset are given below:
# birth rates in 2014/2016, Treatment, and log population in 2014: 
#     zika_Table2.tab in
#     https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ENG0IY
# log population density and proportion of female in 2014:
#     https://www.ibge.gov.br/en/statistics/social/income-expenditure-and-consumption/18391-2010-population-census.html?=&t=resultados
# 
# birth rates in 2013:
#     estimated population:
#     https://www.ibge.gov.br/estatisticas/sociais/populacao/9103-estimativas-de-populacao.html?edicao=17283&t=downloads
#     number of new births:
#     https://opendatasus.saude.gov.br/dataset/sistema-de-informacao-sobre-nascidos-vivos-sinasc/resource/2b5225fd-520d-4e5b-92e8-d7c801d16090
################################################################################

N <- 673
LongData <- read.csv("Zika_Brazil_2013.csv")
LongData <- LongData[,c("Y","Trt","Time","X_1","X_2","X_3","ID")]
Yn1      <- LongData$Y[LongData$Time==-1]    # birth rate in 2013
Y0       <- LongData$Y[LongData$Time==0]     # birth rate in 2014
Y1       <- LongData$Y[LongData$Time==1]     # birth rate in 2016
X1       <- LongData$X_1[LongData$Time==0]   # log-population in 2014
X2       <- LongData$X_2[LongData$Time==0]   # log-population density in 2014
X3       <- LongData$X_3[LongData$Time==0]   # proportion of females in 2014
Trt      <- LongData$Trt[LongData$Time==0]   # Ind(municipality \in PE)

################################################################################
# Re-scale variables for numerical stability
################################################################################

Data     <- data.frame(cbind(Yn1,Y0,Y1,Trt,X1,X2,X3))
mYn1     <- mean(Data$Yn1[Data$Trt==0])    
mY0      <- mean(Data$Y0[Data$Trt==0])
mY1      <- mean(Data$Y1[Data$Trt==0])
sYn1     <- sd(Data$Yn1[Data$Trt==0])
sY0      <- sd(Data$Y0[Data$Trt==0])
sY1      <- sd(Data$Y1[Data$Trt==0])
mX1      <- mean(Data$X1[Data$Trt==0])
mX2      <- mean(Data$X2[Data$Trt==0])
mX3      <- mean(Data$X3[Data$Trt==0])
sX1      <- sd(Data$X1[Data$Trt==0])
sX2      <- sd(Data$X2[Data$Trt==0])
sX3      <- sd(Data$X3[Data$Trt==0])

Data$Yn1 <- (Data$Yn1-mYn1)/sYn1
Data$Y0  <- (Data$Y0-mY0)/sY0
Data$Y1  <- (Data$Y1-mY1)/sY1
Data$X1  <- (Data$X1-mX1)/sX1
Data$X2  <- (Data$X2-mX2)/sX2
Data$X3  <- (Data$X3-mX3)/sX3
Yn1      <- Data$Yn1
Y1       <- Data$Y1
Y0       <- Data$Y0
X1       <- as.numeric(Data$X1)
X2       <- as.numeric(Data$X2)
X3       <- as.numeric(Data$X3)

################################################################################
# Standard DiD under parallel trends
################################################################################

LongData.W2014  <- LongData[LongData$Time!=-1,]
LongData.W2013  <- LongData[LongData$Time!=0,]
LongData.W20134 <- LongData

LongData.W2013$Time  <- as.integer((LongData.W2013$Time+1)/2)

PT.OR.W2014 <- att_gt(yname="Y",
                      tname="Time",
                      idname="ID",
                      gname="Trt",
                      xformla=~1+X_1+X_2+X_3,
                      data=LongData.W2014,
                      bstrap = F,
                      est_method = "reg") 
PT.OR.W2013 <- att_gt(yname="Y",
                      tname="Time",
                      idname="ID",
                      gname="Trt",
                      xformla=~1+X_1+X_2+X_3,
                      data=LongData.W2013,
                      bstrap = F,
                      est_method = "reg") 
PT.OR.W20134 <- att_gt(yname="Y",
                       tname="Time",
                       idname="ID",
                       gname="Trt",
                       xformla=~1+X_1+X_2+X_3,
                       data=LongData.W20134,
                       bstrap = F,
                       est_method = "reg") 

PP <- function(ll){
  LL <- length(ll$att)
  c(ll$att[LL],
    ll$se[LL], 
    ll$att[LL]-qnorm(0.975)*ll$se[LL], 
    ll$att[LL]+qnorm(0.975)*ll$se[LL])
}

Result.PT <- rbind( c(PP(PT.OR.W2014)),
                    c(PP(PT.OR.W2013)),
                    c(PP(PT.OR.W20134)))

Result.PT
# DID-Parallel Trends
# [1,] -1.156146 -1.545715 -0.7665766
# [2,] -1.041166 -1.424079 -0.6582523
# [3,] -1.041166 -1.424079 -0.6582523

################################################################################
# Crude estimatnd E(Y|A=1) - E(Y|A=0)
################################################################################

mean(Data$Y1[Data$Trt==1]*sY1) - mean(Data$Y1[Data$Trt==0]*sY1)
t.test(Data$Y1[Data$Trt==1]*sY1, Data$Y1[Data$Trt==0]*sY1)

# Crude
# Est 3.384296
# LB  2.952956
# UB  3.815636

################################################################################
# Single proxy control: find optimal regularization parameters
################################################################################

Result.PS.D <- Result.OR.D <- Result.DR.D <- list()
CV.PS <- CV.OR <- CV.DR <- list()
PSCV <- ORCV <- DRCV <- list()
Folder     <- c("1NCO1","1NCO2","2NCO")
Y.Origin   <- Data$Y1

for(TYPE in 1:3){
  
  PSCV[[TYPE]] <- as.matrix(read.csv(sprintf("CV/Merge_ResultPS_I_%s.csv",
                                             Folder[TYPE])))
  ORCV[[TYPE]] <- as.matrix(read.csv(sprintf("CV/Merge_ResultOR_I_%s.csv",
                                             Folder[TYPE])))
  DRCV[[TYPE]] <- as.matrix(read.csv(sprintf("CV/Merge_ResultDR_I_%s.csv",
                                             Folder[TYPE])))
  PSCV[[TYPE]] <- PSCV[[TYPE]][order(PSCV[[TYPE]][,1]),]
  ORCV[[TYPE]] <- ORCV[[TYPE]][order(ORCV[[TYPE]][,1]),]
  DRCV[[TYPE]] <- DRCV[[TYPE]][order(DRCV[[TYPE]][,1]),]
  
  Summary <- function(tt){
    
    MMM <- function(VV){ 
      WEIGHT <- diag(rep(1,dim(VV)[2]))
      t(apply(VV,2,mean))%*% WEIGHT %*%(apply(VV,2,mean))
    }
    
    c( MMM(PSCV[[TYPE]][PSCV[[TYPE]][,1]==tt,-c(1:2)]),
       MMM(ORCV[[TYPE]][ORCV[[TYPE]][,1]==tt,-c(1:2)]),
       MMM(DRCV[[TYPE]][DRCV[[TYPE]][,1]==tt,-c(1:2)]))
  }
  
  GRID <- seq(-5.5,-1,by=0.5)
  LOO <- sapply(GRID, Summary)
  
  PS.Summary <- cbind(c(GRID),
                      rep(c(1),each=length(GRID)),
                      as.vector(t(LOO[c(1),])))
  OR.Summary <- cbind(c(GRID),
                      rep(c(1),each=length(GRID)),
                      as.vector(t(LOO[c(2),])))
  DR.Summary <- cbind(c(GRID),
                      rep(c(1),each=length(GRID)),
                      rep(c(1),each=length(GRID)),
                      as.vector(t(LOO[c(3),])))
  
  PS.Summary.Smooth <- PS.Summary
  PS.Summary.Smooth[1:length(GRID),3] <- 
    sapply(GRID,function(vv){median( PS.Summary[ PS.Summary[,2]==1 & 
                                                   abs(GRID-vv)<=0.25,3 ] )})
  
  OR.Summary.Smooth <- OR.Summary
  OR.Summary.Smooth[1:length(GRID),3] <- 
    sapply(GRID,function(vv){median( OR.Summary[ OR.Summary[,2]==1 & 
                                                   abs(GRID-vv)<=0.25,3 ] )})
  
  DR.Summary.Smooth <- DR.Summary
  DR.Summary.Smooth[1:length(GRID),4] <- 
    sapply(GRID,function(vv){median( DR.Summary[ DR.Summary[,3]==1 & 
                                                   abs(GRID-vv)<=0.25,4 ] )})
  
  layout(matrix(1:6,2,3,byrow=T))
  par(mar=c(2,2,1,1))
  
  plot(PS.Summary.Smooth[,c(1,3)],col=PS.Summary.Smooth[,2],pch=19)
  plot(OR.Summary.Smooth[,c(1,3)],col=OR.Summary.Smooth[,2],pch=19)
  plot(DR.Summary.Smooth[,c(1,4)],
       col=DR.Summary.Smooth[,2],pch=-15+4+DR.Summary.Smooth[,3]*15)
  
  FF <- function(DDD){ 
    RRR <- unique(DDD[,1])
    RRR2 <- RRR*0
    for(bb in 1:length(RRR)){
      RRR2[bb] <- sd(DDD[ abs(DDD[,1]-RRR[bb])<=0.05 ,2])
    }
    cbind( RRR,RRR2 )
  }
  
  DDD <- data.frame(PSCV[[1]])
  plot( FF(DDD) ,pch=19)
  DDD <- data.frame(ORCV[[1]])
  plot( FF(DDD) ,pch=19)
  DDD <- data.frame(DRCV[[1]])
  plot( FF(DDD) ,pch=19)
  
  CV.PS[[TYPE]] <- PS.Summary.Smooth[which.min(PS.Summary.Smooth[,3]),]
  CV.OR[[TYPE]] <- OR.Summary.Smooth[which.min(OR.Summary.Smooth[,3]),]
  CV.DR[[TYPE]] <- DR.Summary.Smooth[which.min(DR.Summary.Smooth[,4]),] 
  
}

CV.PS
CV.OR
CV.DR

################################################################################
# Single proxy control: parametric estimators
################################################################################

NUMNN <- 0
Wmat.OR <- Wmat.2.D <- STARTING.PS <- STARTING.OR <- list()

Ymat.PS    <- cbind(1, Y1, X1, X2, X3 )
Ymat.2.D   <- cbind(1, Y1, X1, X2, X3, Y1*X1, Y1*X2, Y1*X3 )

STARTING.PS <- rep(0,dim(Ymat.PS)[2])
STARTING.PS[1] <- -2

SUMMARY <- list()

for(TYPE in 1:3){
  
  if(TYPE==1){
    Wmat.OR[[1]]    <- cbind(1, Y0, X1, X2, X3  )
    Wmat.2.D[[1]]   <- cbind(1, Y0, X1, X2, X3, Y0*X1, Y0*X2, Y0*X3 )
  } else if(TYPE==2){
    Wmat.OR[[2]]    <- cbind(1, Yn1, X1, X2, X3  )
    Wmat.2.D[[2]]   <- cbind(1, Yn1, X1, X2, X3, Yn1*X1, Yn1*X2, Yn1*X3 )
  } else if(TYPE==3){
    Wmat.OR[[3]]    <- cbind(1, Yn1, Y0, X1, X2, X3  )
    Wmat.2.D[[3]]   <- cbind(1, Yn1, Y0, X1, X2, X3, Yn1*X1, 
                             Yn1*X2, Yn1*X3, Y0*X1, Y0*X2, Y0*X3 )
  }
  
  STARTING.OR[[TYPE]] <- rep(0,dim(Wmat.OR[[TYPE]])[2])
  
  Result.PS.D[[TYPE]] <- 
    SPC.PS(Y.Origin,Trt,
           Ymat=Ymat.PS,
           Wmat=Wmat.2.D[[TYPE]],
           penalty=10^(CV.PS[[TYPE]][1]),
           starting=STARTING.PS,
           Var.Type = "Eff",
           Optim.Method="Nelder-Mead",
           radius=1,
           Num.NN=NUMNN) 
  
  Result.OR.D[[TYPE]] <- 
    SPC.OR(Y.Origin,Trt,
           Ymat=Ymat.2.D,
           Wmat=Wmat.OR[[TYPE]],
           penalty=10^(CV.OR[[TYPE]][1]),
           starting=STARTING.OR[[TYPE]],
           Var.Type = "Eff",
           Optim.Method="Nelder-Mead",
           radius=1,
           Num.NN=NUMNN)
  
  Result.DR.D[[TYPE]] <- 
    SPC.DR(Y.Origin,Trt,
           Ymat.PS=Ymat.PS,
           Wmat.PS=Wmat.2.D[[TYPE]],
           Ymat.OR=Ymat.2.D,
           Wmat.OR=Wmat.OR[[TYPE]],
           penalty=10^(CV.DR[[TYPE]][1]),
           Var.Type = "Eff",
           Optim.Method="Nelder-Mead",
           radius=1,
           COEF.PS=STARTING.PS,
           COEF.OR=STARTING.OR[[TYPE]], 
           Num.NN=NUMNN)
  
  
  RESULT2 <- rbind( c(Result.PS.D[[TYPE]]$Est,Result.PS.D[[TYPE]]$SE),
                    c(Result.OR.D[[TYPE]]$Est,Result.OR.D[[TYPE]]$SE),
                    c(Result.DR.D[[TYPE]]$Est,Result.DR.D[[TYPE]]$SE) )
  
  RESULT3 <- RESULT2*sY1
  
  SUMMARY[[TYPE]] <- cbind(RESULT3,
                           RESULT3[,1]-qnorm(0.975)*RESULT3[,2],
                           RESULT3[,1]+qnorm(0.975)*RESULT3[,2]) 
  
}

SUMMARY

# [[1]]
# [,1]      [,2]      [,3]       [,4]
# [1,] -2.297552 0.4915874 -3.261046 -1.3340589
# [2,] -3.446455 0.5454687 -4.515554 -2.3773556
# [3,] -1.832753 0.5187626 -2.849509 -0.8159971
# 
# [[2]]
# [,1]      [,2]      [,3]      [,4]
# [1,] -2.333947 0.3652781 -3.049879 -1.618015
# [2,] -3.560421 0.5198239 -4.579257 -2.541585
# [3,] -2.235067 0.5023708 -3.219695 -1.250438
# 
# [[3]]
# [,1]      [,2]      [,3]      [,4]
# [1,] -2.461953 0.4567892 -3.357243 -1.566663
# [2,] -3.598942 0.7027659 -4.976338 -2.221547
# [3,] -2.181911 0.4151927 -2.995674 -1.368149

################################################################################
# Over-identification
################################################################################

Var.Type <- "Eff"
Num.NN <- 0
radius <- 1

data.PS <- data.OR <- data.DR <- list()
for(TYPE in 1:3){
  data.PS[[TYPE]] <- data.frame(cbind(Y1, Trt, Ymat.PS, Wmat.2.D[[TYPE]]))
  colnames(data.PS[[TYPE]]) <- c( "Outcome",
                                  "Trt",
                                  sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                                  sprintf("W_%0.3d",1:dim(Wmat.2.D[[TYPE]])[2]))
  data.OR[[TYPE]] <- data.frame(cbind(Y1, Trt, Ymat.2.D, Wmat.OR[[TYPE]]))
  colnames(data.OR[[TYPE]]) <- c( "Outcome",
                                  "Trt",
                                  sprintf("Y_%0.3d",1:dim(Ymat.2.D)[2]),
                                  sprintf("W_%0.3d",1:dim(Wmat.OR[[TYPE]])[2]))
  
  data.DR[[TYPE]] <- data.frame(cbind(Y1, Trt, Ymat.PS, Wmat.2.D[[TYPE]], 
                                      Ymat.2.D, Wmat.OR[[TYPE]]))
  
  colnames(data.DR[[TYPE]]) <- c( "Outcome",
                                  "Trt",
                                  sprintf("PY_%0.3d",1:dim(Ymat.PS)[2]),
                                  sprintf("PW_%0.3d",1:dim(Wmat.2.D[[TYPE]])[2]),
                                  sprintf("OY_%0.3d",1:dim(Ymat.2.D)[2]),
                                  sprintf("OW_%0.3d",1:dim(Wmat.OR[[TYPE]])[2]))
  
}

theta.PS <- theta.OR <- theta.DR <- list()
for(TYPE in 1:3){
  theta.PS[[TYPE]] <- Result.PS.D[[TYPE]]$Full.Est[-c(1)]
  theta.OR[[TYPE]] <- Result.OR.D[[TYPE]]$Full.Est[-c(1)]
  theta.DR[[TYPE]] <- Result.DR.D[[TYPE]]$Full.Est[-c(1)]
}

Res.PS <- Weight.PS <- 
  Res.OR <- Weight.OR <- 
  Res.DR <- Weight.DR <- list()
for(TYPE in 1:3){
  Res.PS[[TYPE]] <- GMM.Moment.PS(theta.PS[[TYPE]], as.matrix(data.PS[[TYPE]]))
  Weight.PS[[TYPE]] <- ginv( t(Res.PS[[TYPE]])%*%(Res.PS[[TYPE]])/N )
  
  Res.OR[[TYPE]] <- GMM.Moment.OR(theta.OR[[TYPE]], as.matrix(data.OR[[TYPE]]))
  Weight.OR[[TYPE]] <- ginv( t(Res.OR[[TYPE]])%*%(Res.OR[[TYPE]])/N )
  
  Res.DR[[TYPE]] <- GMM.Moment.DR(theta.DR[[TYPE]], as.matrix(data.DR[[TYPE]]))
  Weight.DR[[TYPE]] <- ginv( t(Res.DR[[TYPE]])%*%(Res.DR[[TYPE]])/N )
}

Res.PS.Full <- cbind(Res.PS[[1]],Res.PS[[2]],Res.PS[[3]])
Res.OR.Full <- cbind(Res.OR[[1]],Res.OR[[2]],Res.OR[[3]])
Res.DR.Full <- cbind(Res.DR[[1]],Res.DR[[2]],Res.DR[[3]])

BLK.DIAG.MAT <- function(VV){
  VV.Full <- matrix(0,
                    dim(VV[[1]])[1]+dim(VV[[2]])[1]+dim(VV[[3]])[1],
                    dim(VV[[1]])[2]+dim(VV[[2]])[2]+dim(VV[[3]])[2])
  VV.Full[1:dim(VV[[1]])[1],1:dim(VV[[1]])[2]] <- VV[[1]]
  VV.Full[dim(VV[[1]])[1]+1:dim(VV[[2]])[1],
          dim(VV[[1]])[2]+1:dim(VV[[2]])[2]] <- VV[[2]]
  VV.Full[dim(VV[[1]])[1]+dim(VV[[2]])[1]+1:dim(VV[[3]])[1],
          dim(VV[[1]])[2]+dim(VV[[2]])[2]+1:dim(VV[[3]])[2]] <- VV[[3]]
  return(VV.Full)
}

Weight.PS.Full <- BLK.DIAG.MAT(Weight.PS)
Weight.OR.Full <- BLK.DIAG.MAT(Weight.OR)
Weight.DR.Full <- BLK.DIAG.MAT(Weight.DR)

Grad.PS <- GMM.Moment.Grad.Merge(theta.PS,data.PS,type="PS")
Grad.OR <- GMM.Moment.Grad.Merge(theta.OR,data.OR,type="OR")
Grad.DR <- GMM.Moment.Grad.Merge(theta.DR,data.DR,type="DR")

Sand.PS <- 
  solve( (Grad.PS)%*%(Weight.PS.Full)%*%t(Grad.PS) )%*%(Grad.PS)%*%(Weight.PS.Full)
Sand.OR <- 
  solve( (Grad.OR)%*%(Weight.OR.Full)%*%t(Grad.OR) )%*%(Grad.OR)%*%(Weight.OR.Full)
Sand.DR <- 
  solve( (Grad.DR)%*%(Weight.DR.Full)%*%t(Grad.DR) )%*%(Grad.DR)%*%(Weight.DR.Full)

VAR.PS  <- Sand.PS%*%(t(Res.PS.Full)%*%(Res.PS.Full)/N)%*%(t(Sand.PS))/N
VAR.OR  <- Sand.OR%*%(t(Res.OR.Full)%*%(Res.OR.Full)/N)%*%(t(Sand.OR))/N
VAR.DR  <- Sand.DR%*%(t(Res.DR.Full)%*%(Res.DR.Full)/N)%*%(t(Sand.DR))/N

colnames(VAR.PS) <- c("ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                      "ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                      "ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]))

colnames(VAR.OR) <- c("ATT1","ATT0",sprintf("W_%0.3d",1:dim(Wmat.OR[[1]])[2]),
                      "ATT1","ATT0",sprintf("W_%0.3d",1:dim(Wmat.OR[[2]])[2]),
                      "ATT1","ATT0",sprintf("W_%0.3d",1:dim(Wmat.OR[[3]])[2]))

colnames(VAR.DR) <- c("ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                      sprintf("W_%0.3d",1:dim(Wmat.OR[[1]])[2]),
                      "ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                      sprintf("W_%0.3d",1:dim(Wmat.OR[[2]])[2]),
                      "ATT1","ATT0",sprintf("Y_%0.3d",1:dim(Ymat.PS)[2]),
                      sprintf("W_%0.3d",1:dim(Wmat.OR[[3]])[2]))

COMPARE <- function(W1=1,W2=3,type="PS"){
  
  if(type=="PS"){
    MM <- colnames(VAR.PS)
    VV <- VAR.PS*sY1
    EE <- c(Result.PS.D[[1]]$Est,Result.PS.D[[2]]$Est,Result.PS.D[[3]]$Est)*sY1
  } else if (type=="OR"){
    MM <- colnames(VAR.OR)
    VV <- VAR.OR*sY1
    EE <- c(Result.OR.D[[1]]$Est,Result.OR.D[[2]]$Est,Result.OR.D[[3]]$Est)*sY1
  } else if (type=="DR"){
    MM <- colnames(VAR.DR)
    VV <- VAR.DR*sY1
    EE <- c(Result.DR.D[[1]]$Est,Result.DR.D[[2]]$Est,Result.DR.D[[3]]$Est)*sY1
  }
  contrast <- rep(0,dim(VV)[1])
  
  contrast[ which(MM=="ATT1")[W1] ] <- 1
  contrast[ which(MM=="ATT1")[W2] ] <- -1
  contrast[ which(MM=="ATT0")[W1] ] <- -1
  contrast[ which(MM=="ATT0")[W2] ] <- 1
  
  Result <- list()
  Result$Est <- EE[c(W1,W2)]
  Result$Diff <- EE[c(W1)] - EE[c(W2)]
  Result$SE   <- sqrt(as.numeric(t(contrast)%*%VV%*%(contrast)))
  return(Result)
}


COMPARE(1,3,"PS")
COMPARE(2,3,"PS")

COMPARE(1,3,"OR")
COMPARE(2,3,"OR")

COMPARE(1,3,"DR")
COMPARE(2,3,"DR")

################################################################################
# Save Result
################################################################################

Result.PT        <- Result.PT
Result.COCA.Para <- SUMMARY
Result.COCA.Para.PS.Compare.13 <- COMPARE(1,3,"PS")
Result.COCA.Para.PS.Compare.23 <- COMPARE(2,3,"PS")
Result.COCA.Para.OR.Compare.13 <- COMPARE(1,3,"OR")
Result.COCA.Para.OR.Compare.23 <- COMPARE(2,3,"OR")
Result.COCA.Para.DR.Compare.13 <- COMPARE(1,3,"DR")
Result.COCA.Para.DR.Compare.23 <- COMPARE(2,3,"DR")

save(Result.PT,
     Result.COCA.Para,
     Result.COCA.Para.PS.Compare.13,
     Result.COCA.Para.PS.Compare.23,
     Result.COCA.Para.OR.Compare.13,
     Result.COCA.Para.OR.Compare.23,
     Result.COCA.Para.DR.Compare.13,
     Result.COCA.Para.DR.Compare.23,
     file="Result_Parametric.RData")
