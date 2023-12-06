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

alpha.W.Grid  <- seq(0,9,length=901)
alpha.W       <- alpha.W.Grid[BATCH]
NNNUM         <- 4
RR            <- 0.1
Intercept.Num <- 11

if(alpha.W==0){
  
  Result <- matrix(c(0,
                     rep(-2.297552,Intercept.Num),
                     rep(0.4915874,Intercept.Num)),
                   1,2*Intercept.Num+1)
  colnames(Result) <- c("alphaW",
                        sprintf("Original.Est_B%0.3d",1:Intercept.Num),
                        sprintf("Original.SE_B%0.3d",1:Intercept.Num))
  write.csv( Result,
             sprintf("Sensitivity/1NCO1/Result_%0.5d.csv",BATCH),
             row.names=F)
  
  Result <- 
    matrix(c(0,
             rep(-2.333947,Intercept.Num),
             rep(0.3652781,Intercept.Num)),
           1,2*Intercept.Num+1)
  colnames(Result) <- c("alphaW",
                        sprintf("Original.Est_B%0.3d",1:Intercept.Num),
                        sprintf("Original.SE_B%0.3d",1:Intercept.Num))
  write.csv( Result,
             sprintf("Sensitivity/1NCO2/Result_%0.5d.csv",BATCH),
             row.names=F)
  
  Result <- 
    matrix(c(0,
             rep(-2.461953,Intercept.Num),
             rep(0.4567892,Intercept.Num)),
           1,2*Intercept.Num+1)
  colnames(Result) <- c("alphaW",
                        sprintf("Original.Est_B%0.3d",1:Intercept.Num),
                        sprintf("Original.SE_B%0.3d",1:Intercept.Num))
  write.csv( Result,
             sprintf("Sensitivity/2NCO/Result_%0.5d.csv",BATCH),
             row.names=F)
  
} else {
  
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
    
    Result.PS.W.Sen <- list()
    
    for(ii in 1:Intercept.Num){
      STARTING.PS[1] <- -2-alpha.W-seq(-0.5,0.5,length=Intercept.Num)[ii]
      
      if(TYPE==1){
        
        Result.PS.W.Sen[[ii]] <- 
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
        Result.PS.W.Sen[[ii]] <- 
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
        Result.PS.W.Sen[[ii]] <- 
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
      
      if(class(Result.PS.W.Sen[[ii]])=="try-error"){
        Result.PS.W.Sen[[ii]]     <- list()
        Result.PS.W.Sen[[ii]]$Est <- -1000
        Result.PS.W.Sen[[ii]]$SE  <- 0
      }
    }
    
    FE <- function(ll){
      c(Result.PS.W.Sen[[ll]]$Est,
        Result.PS.W.Sen[[ll]]$SE)*sY1
    }
    
    Result <- matrix(c(alpha.W,
                       as.vector(t(sapply(1:Intercept.Num,FE)))),
                     1,2*Intercept.Num+1)
    
    colnames(Result) <- c("alphaW",
                          sprintf("Original.Est_B%0.3d",1:Intercept.Num),
                          sprintf("Original.SE_B%0.3d",1:Intercept.Num))
    
    if(TYPE==1){
      write.csv( Result,
                 sprintf("Sensitivity/1NCO1/Result_%0.5d.csv",BATCH),
                 row.names=F)
      
    } else if (TYPE==2){
      write.csv( Result,
                 sprintf("Sensitivity/1NCO2/Result_%0.5d.csv",BATCH),
                 row.names=F)
      
    } else if (TYPE==3){
      write.csv( Result,
                 sprintf("Sensitivity/2NCO/Result_%0.5d.csv",BATCH),
                 row.names=F)
      
    }
    
  }
}


