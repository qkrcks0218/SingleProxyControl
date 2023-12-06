rm(list=ls())

#############################
# Source Files
#############################

source("0.SPC_Ft.R")

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
# Single proxy control: sensitivity analysis
################################################################################

Wmat.OR <- Wmat.2.D <- STARTING.PS <- STARTING.OR <- list()
Ymat.PS    <- cbind(1, Y1, X1, X2, X3 )
Ymat.2.D   <- cbind(1, Y1, X1, X2, X3, Y1*X1, Y1*X2, Y1*X3 )
STARTING.PS <- rep(0,dim(Ymat.PS)[2])

CV.PS         <- list()
CV.PS[[1]]    <- c(-4,1)
CV.PS[[2]]    <- c(-3.5,1)
CV.PS[[3]]    <- c(-4,1)

NNNUM         <- 4
RR            <- 0.1
Intercept.Num <- 11

alpha.W.Vec <- matrix(c(6,0.69,0.98,4.37,
                        6,0.69,1.49,2.60,
                        6,1.04,1.32,7.03),3,4,byrow=T)

alpha.intercept <- matrix(c(6,1,8,7,
                            6,1,4,3,
                            6,5,3,6),3,4,byrow=T)

Result.PS.W.Sen <- list()

for(TYPE in 1:3){
  
  Result.PS.W.Sen[[TYPE]] <- list()
  
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
  
  STARTING.PS[1] <- -2
  
  Result.PS.W.Sen[[TYPE]][[1]] <- 
    SPC.PS(Y1,
           Trt,
           Ymat=Ymat.PS,
           Wmat=Wmat.2.D[[TYPE]],
           penalty=10^(CV.PS[[TYPE]][1]),
           starting=STARTING.PS,
           Var.Type = "Eff",
           Optim.Method="Nelder-Mead",
           radius=0,
           Num.NN=0) 
  
  for(bb in 2:4){
    
    alpha.W <- alpha.W.Vec[TYPE,bb]
    
    STARTING.PS[1] <- -2-alpha.W-seq(-0.5,0.5,length=Intercept.Num)[alpha.intercept[TYPE,bb]]
    
    
    if(TYPE==1){
      
      Result.PS.W.Sen[[TYPE]][[bb]] <- 
        try(SPC.PS.Sensitivity(Y.Origin = Y1,
                               Trt = Trt,
                               Ymat = Ymat.PS,
                               Wmat = Wmat.2.D[[TYPE]],
                               Purturb = (Y0)*alpha.W,
                               penalty=10^(CV.PS[[TYPE]][1]),
                               starting=STARTING.PS,
                               Var.Type = "Eff",
                               Optim.Method="Nelder-Mead",
                               radius=RR,
                               Num.NN=NNNUM),silent = T)
      
    } else if (TYPE==2){
      Result.PS.W.Sen[[TYPE]][[bb]] <- 
        try(SPC.PS.Sensitivity(Y.Origin = Y1,
                               Trt = Trt,
                               Ymat = Ymat.PS,
                               Wmat = Wmat.2.D[[TYPE]],
                               Purturb = (Yn1)*alpha.W,
                               penalty=10^(CV.PS[[TYPE]][1]),
                               starting=STARTING.PS,
                               Var.Type = "Eff",
                               Optim.Method="Nelder-Mead",
                               radius=RR,
                               Num.NN=NNNUM),silent = T)
    } else if (TYPE==3){
      Result.PS.W.Sen[[TYPE]][[bb]] <- 
        try(SPC.PS.Sensitivity(Y.Origin = Y1,
                               Trt = Trt,
                               Ymat = Ymat.PS,
                               Wmat = Wmat.2.D[[TYPE]],
                               Purturb = (Yn1+Y0)*alpha.W/2,
                               penalty=10^(CV.PS[[TYPE]][1]),
                               starting=STARTING.PS,
                               Var.Type = "Eff",
                               Optim.Method="Nelder-Mead",
                               radius=RR,
                               Num.NN=NNNUM),silent = T)
    }
    
  }
  
  
}

################################################################################
# Single proxy control: sensitivity analysis summary
################################################################################

PRINT <- function(SS,v1,v2,v3,v4,ay){
  
  v1 <- as.numeric(v1)
  v2 <- as.numeric(v2)
  v3 <- as.numeric(v3)
  v4 <- as.numeric(v4)
  
  RR <- 
    rbind( sprintf("\\multirow{4}{*}{$%s$} & $\\alpha_w$ & %0.3f & %0.3f & %0.3f & %0.3f \\\\ \\cline{2-6}", 
                   SS,v1[1],v2[1],v3[1],v4[1]), 
           sprintf(" & $\\alpha_y$ & %0.3f & %0.3f & %0.3f & %0.3f \\\\ \\cline{2-6}", 
                   ay[1],ay[2],ay[3],ay[4]), 
           sprintf(" & Estimate & %0.3f & %0.3f & %0.3f & %0.3f \\\\ \\cline{2-6}", 
                   v1[2],v2[2],v3[2],v4[2]), 
           sprintf(" & 95\\%% CI & (%0.3f,%0.3f) & (%0.3f,%0.3f) & (%0.3f,%0.3f) & (%0.3f,%0.3f) \\\\ \\hline", 
                   v1[4],v1[5],v2[4],v2[5],v3[4],v3[5],v4[4],v4[5]) )
  RR <- data.frame(RR)
  print(RR,row.names=F)
  
}

Result.Raw <- as.matrix( read.csv("Sensitivity/Merge_Sensitivity_1NCO2.csv") )
Result.Merge <- matrix(0,dim(Result.Raw)[1],6)
Result.Merge[,1] <- Result.Raw[,1]
Result.Merge[,2] <- apply(Result.Raw[,1+1:11],1,median)
Result.Merge[,3] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( Result.Raw[ii,1+11+which(Result.Raw[ii,1+1:11]-
                                                  Result.Merge[ii,2]==0)] )[1]
         })
Result.Merge[,4:5] <- cbind(Result.Merge[,2]-qnorm(0.975)*Result.Merge[,3],
                            Result.Merge[,2]+qnorm(0.975)*Result.Merge[,3])
Result.Merge[,6] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( which(Result.Raw[ii,1+1:11]-
                               Result.Merge[ii,2]==0) )[1]
         })

UB  <- Result.Merge[,5]
Est <- Result.Merge[,2]
plot(Result.Merge[,1], UB, col=as.numeric(UB>0)+1, pch=19, ylim=c(-0.5,3.5),type='b')
head(Result.Merge[UB>0,1]); head(Result.Merge[UB>0,])
plot(Result.Merge[,1], Est, col=as.numeric(Est>0)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[Est>0,1]); head(Result.Merge[Est>0,])
plot(Result.Merge[,1], UB, col=as.numeric(UB>2.953)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[UB>2.953,1]); head(Result.Merge[UB>2.953,])

ay1 <- c(Result.PS.W.Sen[[2]][[1]]$Full.Est[5],
         Result.PS.W.Sen[[2]][[2]]$Full.Est[5],
         Result.PS.W.Sen[[2]][[3]]$Full.Est[5],
         Result.PS.W.Sen[[2]][[4]]$Full.Est[5])

RR1 <- PRINT( "W_1",
              Result.Merge[1,], 
              Result.Merge[alpha.W.Vec[2,2]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[2,3]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[2,4]==Result.Merge[,1],],
              ay1)

Result.Raw <- as.matrix( read.csv("Sensitivity/Merge_Sensitivity_1NCO1.csv") )
Result.Merge <- matrix(0,dim(Result.Raw)[1],6)
Result.Merge[,1] <- Result.Raw[,1]
Result.Merge[,2] <- apply(Result.Raw[,1+1:11],1,median)
Result.Merge[,3] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( Result.Raw[ii,1+11+which(Result.Raw[ii,1+1:11]-
                                                  Result.Merge[ii,2]==0)] )[1]
         })
Result.Merge[,4:5] <- cbind(Result.Merge[,2]-qnorm(0.975)*Result.Merge[,3],
                            Result.Merge[,2]+qnorm(0.975)*Result.Merge[,3])
Result.Merge[,6] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( which(Result.Raw[ii,1+1:11]-
                               Result.Merge[ii,2]==0) )[1]
         })

UB  <- Result.Merge[,5]
Est <- Result.Merge[,2]
plot(Result.Merge[,1], UB, col=as.numeric(UB>0)+1, pch=19, ylim=c(-0.5,3.5),type='b')
head(Result.Merge[UB>0,1]); head(Result.Merge[UB>0,])
plot(Result.Merge[,1], Est, col=as.numeric(Est>0)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[Est>0,1]); head(Result.Merge[Est>0,])
plot(Result.Merge[,1], UB, col=as.numeric(UB>2.953)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[UB>2.953,1]); head(Result.Merge[UB>2.953,])

ay2 <- c(Result.PS.W.Sen[[1]][[1]]$Full.Est[5],
         Result.PS.W.Sen[[1]][[2]]$Full.Est[5],
         Result.PS.W.Sen[[1]][[3]]$Full.Est[5],
         Result.PS.W.Sen[[1]][[4]]$Full.Est[5])


RR2 <- PRINT( "W_2",
              Result.Merge[1,], 
              Result.Merge[alpha.W.Vec[1,2]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[1,3]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[1,4]==Result.Merge[,1],],
              ay2)

Result.Raw <- as.matrix( read.csv("Sensitivity/Merge_Sensitivity_2NCO.csv") )
Result.Merge <- matrix(0,dim(Result.Raw)[1],6)
Result.Merge[,1] <- Result.Raw[,1]
Result.Merge[,2] <- apply(Result.Raw[,1+1:11],1,median)
Result.Merge[,3] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( Result.Raw[ii,1+11+which(Result.Raw[ii,1+1:11]-
                                                  Result.Merge[ii,2]==0)] )[1]
         })
Result.Merge[,4:5] <- cbind(Result.Merge[,2]-qnorm(0.975)*Result.Merge[,3],
                            Result.Merge[,2]+qnorm(0.975)*Result.Merge[,3])
Result.Merge[,6] <- 
  sapply(1:dim(Result.Raw)[1],
         function(ii){
           as.numeric( which(Result.Raw[ii,1+1:11]-
                               Result.Merge[ii,2]==0) )[1]
         })

UB  <- Result.Merge[,5]
Est <- Result.Merge[,2]
plot(Result.Merge[,1], UB, col=as.numeric(UB>0)+1, pch=19, ylim=c(-0.5,3.5),type='b')
head(Result.Merge[UB>0,1]); head(Result.Merge[UB>0,])
plot(Result.Merge[,1], Est, col=as.numeric(Est>0)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[Est>0,1]); head(Result.Merge[Est>0,])
plot(Result.Merge[,1], UB, col=as.numeric(UB>2.953)+1, pch=19, ylim=c(-0.5,3.5))
head(Result.Merge[UB>2.953,1]); head(Result.Merge[UB>2.953,])

ay3 <- c(Result.PS.W.Sen[[3]][[1]]$Full.Est[5],
         Result.PS.W.Sen[[3]][[2]]$Full.Est[5],
         Result.PS.W.Sen[[3]][[3]]$Full.Est[5],
         Result.PS.W.Sen[[3]][[4]]$Full.Est[5])


RR3 <- PRINT( "(W_1,W_2)",
              Result.Merge[1,], 
              Result.Merge[alpha.W.Vec[3,2]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[3,3]==Result.Merge[,1],],
              Result.Merge[alpha.W.Vec[3,4]==Result.Merge[,1],],
              ay3)

print(rbind(RR1,RR2,RR3),row.names=F)
