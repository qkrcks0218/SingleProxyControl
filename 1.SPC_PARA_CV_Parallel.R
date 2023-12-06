################################################################################
# Parallel Computing
################################################################################

args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

################################################################################
# Source File
################################################################################

source("0.SPC_Ft.R")

CV.PS <- function(cv,iter,typeOR){
  PS.Penalty <- SPC.PS(Y.Origin[-cv],
                       Trt[-cv],
                       Ymat=Ymat.PS[-cv,],
                       Wmat=Wmat.2.D[[typeOR]][-cv,],
                       penalty=10^iter,
                       starting=STARTING.PS,
                       radius=1,
                       Var.Type = "Eff",
                       Optim.Method="Nelder-Mead",
                       Num.NN=NUMNN)
  InvPr0 <- (1+exp( Ymat.PS[cv,]%*%PS.Penalty$Full.Est[-c(1,2,3)] ))
  
  PSW <- link(Ymat.PS[cv,]%*%PS.Penalty$Full.Est[-c(1,2,3)])
  
  ATT1 <- Trt[cv]*(Y.Origin[cv]-PS.Penalty$Full.Est[2])
  ATT0 <- (1-Trt[cv])*(Y.Origin[cv])*(PSW) - 
    (1-Trt[cv])*(PSW)*PS.Penalty$Full.Est[3]
  
  Residual <- c(PS.Penalty$Full.Est[1],
                ATT1,ATT0,as.numeric( (1-Trt[cv])*InvPr0-1 )*Wmat.2.D[[typeOR]][cv,])
  
  return(Residual)
}

CV.OR <- function(cv,iter,typeOR){
  OR.Penalty <- SPC.OR(Y.Origin[-cv],
                       Trt[-cv],
                       Ymat=Ymat.2.D[-cv,],
                       Wmat=Wmat.OR[[typeOR]][-cv,],
                       penalty=10^iter,
                       starting=STARTING.OR[[typeOR]],
                       radius=1,
                       Var.Type = "Eff",
                       Optim.Method="Nelder-Mead",
                       Num.NN=NUMNN)
  
  BR <- as.numeric( Wmat.OR[[typeOR]][cv,]%*%OR.Penalty$Full.Est[-c(1,2,3)] )
  
  ATT1 <- Trt[cv]*(Y.Origin[cv]-OR.Penalty$Full.Est[2])
  ATT0 <- Trt[cv]*(BR-OR.Penalty$Full.Est[3])
  
  Residual <- c(OR.Penalty$Full.Est[1],
                ATT1,ATT0,Ymat.2.D[cv,]*(1-Trt[cv])*(Y.Origin[cv] - BR))
  
  return(Residual)
}

CV.DR <- function(cv,iter,typeOR){
  
  DR.Penalty <- SPC.DR(Y.Origin[-cv],
                       Trt[-cv],
                       Ymat.PS=Ymat.PS[-cv,],
                       Wmat.PS=Wmat.2.D[[typeOR]][-cv,],
                       Ymat.OR=Ymat.2.D[-cv,],
                       Wmat.OR=Wmat.OR[[typeOR]][-cv,],
                       penalty=10^iter,
                       radius=1,
                       Var.Type = "Eff",
                       Optim.Method="Nelder-Mead",
                       COEF.PS=STARTING.PS,
                       COEF.OR=STARTING.OR[[typeOR]],
                       Num.NN=NUMNN)
  
  PI.para <- 
    DR.Penalty$Full.Est[which(substr(names(DR.Penalty$Full.Est),1,2)=="PY")]
  BR.para <- 
    DR.Penalty$Full.Est[which(substr(names(DR.Penalty$Full.Est),1,2)=="OW")]
  
  InvPr0 <- (1+exp( Ymat.PS[cv,]%*%PI.para ))
  Residual.PS <- as.numeric( (1-Trt[cv])*InvPr0-1 )*Wmat.2.D[[typeOR]][cv,]
  
  PSW <- link(Ymat.PS[cv,]%*%PI.para)
  
  BR <- as.numeric( Wmat.OR[[typeOR]][cv,]%*%BR.para )
  Residual.OR <- Ymat.2.D[cv,]*(1-Trt[cv])*(Y.Origin[cv] - BR)
  
  ATT1 <- Trt[cv]*(Y.Origin[cv]-DR.Penalty$Full.Est[2])
  ATT0 <- ((1-Trt[cv])*(Y.Origin[cv] - BR)*PSW+(Trt[cv])*BR) - 
    Trt[cv]*DR.Penalty$Full.Est[3]
  
  Residual <- c(DR.Penalty$Full.Est[1],
                ATT1,ATT0,Residual.PS,Residual.OR)
  return(Residual) 
}

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
# Single proxy control: cross-validation
################################################################################

RESULT     <- matrix(0,2,6)
Y.Origin   <- Data$Y1
Folder     <- c("1NCO1","1NCO2","2NCO")
GRID       <- sort( unique( c( seq(-5.5,-1,by=0.05) ) ) )

NUMNN      <- 0
Wmat.OR    <- Wmat.2.D <- STARTING.PS <- STARTING.OR <- list()

Ymat.PS    <- cbind(1, Y1, X1, X2, X3 )
Ymat.2.D   <- cbind(1, Y1, X1, X2, X3, Y1*X1, Y1*X2, Y1*X3 )

STARTING.PS <- rep(0,dim(Ymat.PS)[2])
STARTING.PS[1] <- -2

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
  
  RRR.PS.I     <- matrix(0,length(GRID),dim(Wmat.2.D[[TYPE]])[2]+4)
  RRR.OR.I     <- matrix(0,length(GRID),dim(Ymat.2.D)[2]+4)
  RRR.DR.I   <- matrix(0,length(GRID),
                       dim(Ymat.2.D)[2]+dim(Wmat.2.D[[TYPE]])[2]+4)
  
  for(tt in 1:length(GRID)){
    
    Lambda <- GRID[tt]
    
    RRR.PS.I[tt,] <- c(Lambda,CV.PS(BATCH,Lambda,TYPE))
    RRR.OR.I[tt,] <- c(Lambda,CV.OR(BATCH,Lambda,TYPE))
    RRR.DR.I[tt,] <- c(Lambda,CV.DR(BATCH,Lambda,TYPE))
    
  }
  
  colnames(RRR.PS.I) <- c("Para", "ATT", "Moment_ATT1", "Moment_ATT0", 
                          sprintf("Moment_X_%0.3d",1:(dim(RRR.PS.I)[2]-4)))
  colnames(RRR.OR.I) <- c("Para", "ATT", "Moment_ATT1", "Moment_ATT0", 
                          sprintf("Moment_X_%0.3d",1:(dim(RRR.OR.I)[2]-4)))
  colnames(RRR.DR.I) <- c("Para", "ATT", "Moment_ATT1", "Moment_ATT0", 
                          sprintf("Moment_X_%0.3d",1:(dim(RRR.DR.I)[2]-4)))
  
  write.csv(RRR.PS.I,sprintf("CV/%s/SingleStart_ResultPS_I_B%0.5d.csv",
                             Folder[TYPE],BATCH+2000),row.names=F)
  write.csv(RRR.OR.I,sprintf("CV/%s/SingleStart_ResultOR_I_B%0.5d.csv",
                             Folder[TYPE],BATCH+2000),row.names=F)
  write.csv(RRR.DR.I,sprintf("CV/%s/SingleStart_ResultDR_I_B%0.6d.csv",
                             Folder[TYPE],BATCH+2000),row.names=F)
  
}

