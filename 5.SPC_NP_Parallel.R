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

source("0.Functions_MM.R")
source("0.SPC_Ft.R")
library(MASS)
library(splines)

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
Wn1      <- LongData$Y[LongData$Time==-1]    # birth rate in 2013
W0       <- LongData$Y[LongData$Time==0]     # birth rate in 2014
Y1       <- LongData$Y[LongData$Time==1]     # birth rate in 2016
X1       <- LongData$X_1[LongData$Time==0]   # log-population in 2014
X2       <- LongData$X_2[LongData$Time==0]   # log-population density in 2014
X3       <- LongData$X_3[LongData$Time==0]   # proportion of females in 2014
A        <- LongData$Trt[LongData$Time==0]   # Ind(municipality \in PE)

################################################################################
# Re-scale variables for numerical stability
################################################################################

Data     <- data.frame(cbind(W0,Wn1,Y1,A,X1,X2,X3))
mWn1     <- mean(Data$Wn1[Data$A==0])    
mW0      <- mean(Data$W0[Data$A==0])
mY1      <- mean(Data$Y1[Data$A==0])
sWn1     <- sd(Data$Wn1[Data$A==0])
sW0      <- sd(Data$W0[Data$A==0])
sY1      <- sd(Data$Y1[Data$A==0])
mX1      <- mean(Data$X1[Data$A==0])
mX2      <- mean(Data$X2[Data$A==0])
mX3      <- mean(Data$X3[Data$A==0])
sX1      <- sd(Data$X1[Data$A==0])
sX2      <- sd(Data$X2[Data$A==0])
sX3      <- sd(Data$X3[Data$A==0])

Data$Wn1 <- (Data$Wn1-mWn1)/sWn1
Data$W0  <- (Data$W0-mW0)/sW0
Data$Y1  <- (Data$Y1-mY1)/sY1
Data$X1  <- (Data$X1-mX1)/sX1
Data$X2  <- (Data$X2-mX2)/sX2
Data$X3  <- (Data$X3-mX3)/sX3
Wn1      <- Data$Wn1
Y1       <- Data$Y1
W0       <- Data$W0
X1       <- as.numeric(Data$X1)
X2       <- as.numeric(Data$X2)
X3       <- as.numeric(Data$X3)

################################################################################
# Minimax Estimation
################################################################################

CF <- 2
NumCV <- 5

set.seed(BATCH)

Pos11 <- sample(which(Data$A==1), 92)
Pos01 <- sample(which(Data$A==0), 244)
Pos12 <- setdiff(which(Data$A==1),Pos11)
Pos02 <- setdiff(which(Data$A==0),Pos01)

Data <- Data[c( c(Pos11,Pos01)[sample(336,336)],
                c(Pos12,Pos02)[sample(337,337)]),]

SS.Index <- list()
SS.Index[[1]] <-   1:336
SS.Index[[2]] <- 337:673

CV.Index <- list()
CV.Index[[1]] <- CV.Index[[2]] <- list()
CV.CUT <- seq(0,336,length=NumCV+1)
CV.Index[[1]][[1]] <- 0*67+1:67
CV.Index[[1]][[2]] <- 1*67+1:67
CV.Index[[1]][[3]] <- 2*67+1:67
CV.Index[[1]][[4]] <- 3*67+1:67
CV.Index[[1]][[5]] <- 4*67+1:68

CV.Index[[2]][[1]] <- 0*67+1:67
CV.Index[[2]][[2]] <- 1*67+1:67
CV.Index[[2]][[3]] <- 2*67+1:67
CV.Index[[2]][[4]] <- 3*67+1:68
CV.Index[[2]][[5]] <- 4*67+1+1:68

ncIF <- list()

for(TYPE in 1:3){
  
  X.pos <- which( substr(colnames(Data),1,1)=="X" )
  A.pos <- which( substr(colnames(Data),1,1)=="A" )
  Y.pos <- which( substr(colnames(Data),1,1)=="Y" )
  if(TYPE==1){
    W.pos <- which( substr(colnames(Data),1,1)=="W" )[1]
  } else if (TYPE==2){
    W.pos <- which( substr(colnames(Data),1,1)=="W" )[2]
  } else if (TYPE==3){
    W.pos <- which( substr(colnames(Data),1,1)=="W" )
  }
  
  Y  <- as.matrix(Data[,Y.pos])
  A  <- as.matrix(Data[,A.pos])
  W  <- as.matrix(Data[,W.pos])
  X  <- as.matrix(Data[,X.pos])
  YX <- as.matrix(cbind(Y,X))
  WX <- as.matrix(cbind(W,X))
  
  N.MM <- Y.MM <- W.MM <- A.MM <- X.MM <- YX.MM <- WX.MM <- list()
  for(ss in 1:CF){
    N.MM [[ss]] <- length(SS.Index[[ss]])
    Y.MM [[ss]] <- Y[SS.Index[[ss]],]
    W.MM [[ss]] <- W[SS.Index[[ss]],]
    A.MM [[ss]] <- A[SS.Index[[ss]],]
    X.MM [[ss]] <- X[SS.Index[[ss]],]
    YX.MM[[ss]] <- YX[SS.Index[[ss]],]
    WX.MM[[ss]] <- WX[SS.Index[[ss]],]
  }
  
  ParaGrid.PerWX.TarYX <- expand.grid(as.numeric(FT_BWHeuristic(WX))*c(3,3.25),
                                      as.numeric(FT_BWHeuristic(YX))*c(3,3.25),
                                      10^c(-4,-3.75,-3.5),
                                      10^c(-4,-3.75,-3.5))
  
  ParaGrid.PerYX.TarWX <- expand.grid(as.numeric(FT_BWHeuristic(YX))*c(3,3.25),
                                      as.numeric(FT_BWHeuristic(WX))*c(3,3.25),
                                      10^c(-4,-3.75,-3.5),
                                      10^c(-4,-3.75,-3.5))
  
  ParaGrid.pi  <- as.matrix(ParaGrid.PerWX.TarYX)
  ParaGrid.br  <- as.matrix(ParaGrid.PerYX.TarWX)
  
  RISK.pi <- RISK.br <- list() 
  for(ss in 1:CF){
    RISK.pi[[ss]] <- RISK.br[[ss]] <- matrix(10^5,dim(ParaGrid.pi)[1],4*NumCV) 
    
    colnames(RISK.pi[[ss]]) <-  colnames(RISK.br[[ss]]) <-  
      c(sprintf("Pr%0.1d",1:NumCV),
        sprintf("P0%0.1d",1:NumCV),
        sprintf("Sq%0.1d",1:NumCV),
        sprintf("PMMR%0.1d",1:NumCV))
  }
  
  Bound.pi  <- c(-100,100)
  Bound.br  <- c(-100,100)
  
  ##############################################################################
  # Practical Regularization: estimate of EPS = GMM Estimate * MM Estimate
  ##############################################################################
  
  WX.Spline <- cbind( 1,W,X1,X2,X3 )
  YX.Spline <- cbind( 1,Y,X1,X2,X3 )
  
  WX.Spline.MM <- list()
  YX.Spline.MM <- list()
  WX.Spline.MM[[1]] <- WX.Spline[ SS.Index[[1]], ]
  WX.Spline.MM[[2]] <- WX.Spline[ SS.Index[[2]], ]
  YX.Spline.MM[[1]] <- YX.Spline[ SS.Index[[1]], ]
  YX.Spline.MM[[2]] <- YX.Spline[ SS.Index[[2]], ] 
  
  Lambda <- rep(-3,3)
  
  PI.GMM <- list()
  PI.GMM[[1]] <- SPC.PS(Y[SS.Index[[1]]],
                        A[SS.Index[[1]]],
                        Ymat=YX.Spline[SS.Index[[1]],],
                        Wmat=WX.Spline[SS.Index[[1]],],
                        penalty=10^(Lambda[TYPE]),
                        starting=c(-2,rep(0,4)),
                        Var.Type = "Eff",
                        Optim.Method="Nelder-Mead",
                        radius=1,
                        Num.NN=0)
  
  PI.GMM[[2]] <- SPC.PS(Y[SS.Index[[2]]],
                        A[SS.Index[[2]]],
                        Ymat=YX.Spline[SS.Index[[2]],],
                        Wmat=WX.Spline[SS.Index[[2]],],
                        penalty=10^(Lambda[TYPE]),
                        starting=c(-2,rep(0,4)),
                        Var.Type = "Eff",
                        Optim.Method="Nelder-Mead",
                        radius=1,
                        Num.NN=0)
  
  Pi.Factor.Coef <- list()
  Pi.Factor.Coef[[1]] <- PI.GMM[[1]]$Full.Est[-c(1,2,3)]
  Pi.Factor.Coef[[2]] <- PI.GMM[[2]]$Full.Est[-c(1,2,3)]
  
  ##############################################################################
  # EPS Minimax Estimation: Cross-validation
  ##############################################################################
  
  for(Para.Iter in 1:dim(ParaGrid.pi)[1]){
    
    for(ss in 1:CF){
      for(cv in 1:NumCV){
        
        CV.Split.Index <- list()
        CV.Split.Index[[1]] <- SS.Index[[ss]][ -CV.Index[[ss]][[cv]] ]
        CV.Split.Index[[2]] <- SS.Index[[ss]][ CV.Index[[ss]][[cv]] ]
        
        Factor <- exp(-YX.Spline[SS.Index[[ss]],]%*%Pi.Factor.Coef[[ss]])
        
        Factor.CV <- list()
        Factor.CV[[1]] <- Factor[-CV.Index[[ss]][[cv]]]
        Factor.CV[[2]] <- Factor[CV.Index[[ss]][[cv]]]
        
        Y.CV  <- list()
        A.CV  <- list()
        X.CV  <- list()
        W.CV  <- list()
        WX.CV <- list()
        YX.CV <- list()
        
        pihat.ML.CV <- list()
        
        for(cvest in 1:2){
          Y.CV [[cvest]] <- Y[CV.Split.Index[[cvest]],]
          A.CV [[cvest]] <- A[CV.Split.Index[[cvest]],]
          X.CV [[cvest]] <- X[CV.Split.Index[[cvest]],]
          W.CV [[cvest]] <- W[CV.Split.Index[[cvest]],]
          WX.CV[[cvest]] <- cbind(W.CV[[cvest]],X.CV[[cvest]])
          YX.CV[[cvest]] <- cbind(Y.CV[[cvest]],X.CV[[cvest]])
        }
        
        ## pi
        
        pos.pi.CV <- list()
        Coef.pi.CV <- Intercept.pi.CV <- list()
        
        for(cvest in 1:2){
          pos.pi.CV[[cvest]]  <- 1:length(Y.CV[[cvest]])
          Coef.pi.CV[[cvest]] <- (1-A.CV[[cvest]])*Factor.CV[[cvest]]
          Intercept.pi.CV[[cvest]] <- -A.CV[[cvest]]
        }
        
        CV.pi <- 
          FT_CV_Minimax( Coef.Train      = Coef.pi.CV[[1]][pos.pi.CV[[1]] ],
                         Intercept.Train = Intercept.pi.CV[[1]][pos.pi.CV[[1]] ],
                         X.Perturb.Train = WX.CV[[1]][pos.pi.CV[[1]],],
                         X.Target.Train  = YX.CV[[1]][pos.pi.CV[[1]],],
                         Coef.Valid      = Coef.pi.CV[[2]][pos.pi.CV[[2]] ],
                         Intercept.Valid = Intercept.pi.CV[[2]][pos.pi.CV[[2]] ],
                         X.Perturb.Valid = WX.CV[[2]][pos.pi.CV[[2]],],
                         X.Target.Valid  = YX.CV[[2]][pos.pi.CV[[2]],],
                         bw.Perturb      = ParaGrid.pi[Para.Iter,1],
                         bw.Target       = ParaGrid.pi[Para.Iter,2],
                         lambda.Perturb  = ParaGrid.pi[Para.Iter,3],
                         lambda.Target   = ParaGrid.pi[Para.Iter,4],
                         bound = Bound.pi,
                         SVM.bias = TRUE,
                         bias.input = 1)
        
        RISK.pi[[ss]][Para.Iter,c(cv,cv+NumCV,cv+2*NumCV,cv+3*NumCV)] <- 
          c(CV.pi$Error.Proj,
            CV.pi$Error.Proj.0,
            CV.pi$Error.Sq,
            CV.pi$Error.PMMR)
      }
      print(c(Para.Iter,ss))
    }
    print(c(Para.Iter))
  }
  
  Opt.Para.pi <- list()
  if( dim(RISK.pi[[ss]])[1] > 1){
    for(ss in 1:CF){
      pos.p0 <- which( substr(colnames(RISK.pi[[ss]]),1,4)=="PMMR" )
      Opt.Para.pi[[ss]] <- 
        as.numeric(ParaGrid.pi[which.min(apply(RISK.pi[[ss]][,pos.p0],1,mean)),])
    }
  } else {
    for(ss in 1:CF){
      Opt.Para.pi[[ss]] <- ParaGrid.pi[1,]
    }
  }
  
  ##############################################################################
  # EPS Minimax Estimation: Evaluation
  ##############################################################################
  
  pos.pi.MM <- Coef.pi.MM <- Intercept.pi.MM <- list()
  pi.MM <- list()
  pi.predict <- list()
  
  for(ss in 1:CF){
    
    Factor <- exp(-YX.Spline[SS.Index[[ss]],]%*%Pi.Factor.Coef[[ss]])
    
    pos.pi.MM[[ss]]  <- 1:length(A.MM[[ss]])
    Coef.pi.MM[[ss]] <- (1-A.MM[[ss]])*Factor
    Intercept.pi.MM[[ss]] <- -A.MM[[ss]]
    
    pi.MM[[ss]] <- 
      FT_Minimax( Coef            = Coef.pi.MM[[ss]][pos.pi.MM[[ss]] ],
                  Intercept       = Intercept.pi.MM[[ss]][pos.pi.MM[[ss]] ],
                  X.Perturb       = WX.MM[[ss]][pos.pi.MM[[ss]],],
                  X.Target        = YX.MM[[ss]][pos.pi.MM[[ss]],],
                  bw.Perturb      = Opt.Para.pi[[ss]][1],
                  bw.Target       = Opt.Para.pi[[ss]][2],
                  lambda.Perturb  = Opt.Para.pi[[ss]][3],
                  lambda.Target   = Opt.Para.pi[[ss]][4],
                  SVM.bias =TRUE,
                  bias.input = 1)
    
    pi.predict[[ss]] <- function(YX.New.Input,
                                 YX.New.Spline.Input){
      
      RKHS.Fit <- 
        (FT_RBF(X    = YX.MM[[ss]][pos.pi.MM[[ss]],],
                X.new = YX.New.Input,
                bw    = Opt.Para.pi[[ss]][2])%*%pi.MM[[ss]]$gamma + 
           pi.MM[[ss]]$bias)*(exp(YX.New.Spline.Input%*%Pi.Factor.Coef[[ss]]))
      
      RKHS.Fit
    }
  }
  
  pihat.NoCF <- rep(0,N)
  for(ss in 1:CF){
    pihat.NoCF[SS.Index[[ss]]] <- pi.predict[[ss]](YX.MM[[ss]],
                                                   YX.Spline.MM[[ss]])
  }
  
  pihat.CF <- rep(0,N)
  for(ss in 1:CF){
    pihat.CF[SS.Index[[3-ss]]] <- pi.predict[[ss]](YX.MM[[3-ss]],
                                                   YX.Spline.MM[[3-ss]])
  }
  
  pihat.NoCF.MM <- list()
  for(ss in 1:CF){
    pihat.NoCF.MM[[ss]]  <- pihat.NoCF[SS.Index[[ss]] ]
  } 
  
  ##############################################################################
  # Bridge Function Minimax Estimation: Cross-validation
  ##############################################################################
  
  for(Para.Iter in 1:dim(ParaGrid.br)[1]){
    
    for(ss in 1:CF){
      for(cv in 1:NumCV){
        
        
        CV.Split.Index <- list()
        CV.Split.Index[[1]] <- SS.Index[[ss]][ -CV.Index[[ss]][[cv]] ]
        CV.Split.Index[[2]] <- SS.Index[[ss]][ CV.Index[[ss]][[cv]] ]
        
        Y.CV  <- list()
        A.CV  <- list()
        X.CV  <- list()
        W.CV  <- list()
        WX.CV <- list()
        YX.CV <- list()
        
        brhat.ML.CV <- list()
        
        for(cvest in 1:2){
          Y.CV [[cvest]] <- Y[CV.Split.Index[[cvest]],]
          A.CV [[cvest]] <- A[CV.Split.Index[[cvest]],]
          X.CV [[cvest]] <- X[CV.Split.Index[[cvest]],]
          W.CV [[cvest]] <- W[CV.Split.Index[[cvest]],]
          WX.CV[[cvest]] <- cbind(W.CV[[cvest]],X.CV[[cvest]])
          YX.CV[[cvest]] <- cbind(Y.CV[[cvest]],X.CV[[cvest]])
        }
        
        ## br
        
        pos.br.CV <- list()
        Coef.br.CV <- Intercept.br.CV <- list()
        
        for(cvest in 1:2){
          pos.br.CV[[cvest]]       <- 1:length(A.CV[[cvest]])
          Coef.br.CV[[cvest]]      <- (1-A.CV[[cvest]])
          Intercept.br.CV[[cvest]] <- -(1-A.CV[[cvest]])*Y.CV[[cvest]]
        }
        
        CV.br <- 
          FT_CV_Minimax( Coef.Train      = Coef.br.CV[[1]][pos.br.CV[[1]] ],
                         Intercept.Train = Intercept.br.CV[[1]][pos.br.CV[[1]] ],
                         X.Perturb.Train = YX.CV[[1]][pos.br.CV[[1]],],
                         X.Target.Train  = WX.CV[[1]][pos.br.CV[[1]],],
                         Coef.Valid      = Coef.br.CV[[2]][pos.br.CV[[2]] ],
                         Intercept.Valid = Intercept.br.CV[[2]][pos.br.CV[[2]] ],
                         X.Perturb.Valid = YX.CV[[2]][pos.br.CV[[2]],],
                         X.Target.Valid  = WX.CV[[2]][pos.br.CV[[2]],],
                         bw.Perturb      = ParaGrid.br[Para.Iter,1],
                         bw.Target       = ParaGrid.br[Para.Iter,2],
                         lambda.Perturb  = ParaGrid.br[Para.Iter,3],
                         lambda.Target   = ParaGrid.br[Para.Iter,4],
                         bound = Bound.br,
                         SVM.bias = FALSE,
                         bias.input = NULL)
        
        RISK.br[[ss]][Para.Iter,c(cv,cv+NumCV,cv+2*NumCV,cv+3*NumCV)] <- 
          c(CV.br$Error.Proj,
            CV.br$Error.Proj.0,
            CV.br$Error.Sq,
            CV.br$Error.PMMR)
      }
      print(c(Para.Iter,ss))
    }
    print(c(Para.Iter))
  }
  
  Opt.Para.br <- list()
  if( dim(RISK.br[[ss]])[1] > 1){
    for(ss in 1:CF){
      pos.p0 <- which( substr(colnames(RISK.br[[ss]]),1,4)=="PMMR" )
      Opt.Para.br[[ss]] <- 
        as.numeric(ParaGrid.br[which.min(apply(RISK.br[[ss]][,pos.p0],1,mean)),])
    }
  } else {
    for(ss in 1:CF){
      Opt.Para.br[[ss]] <- ParaGrid.br[1,]
    }
  }
  
  ##############################################################################
  # Bridge Function Minimax Estimation: Evaluation
  ##############################################################################
  
  pos.br.MM <- Coef.br.MM <- Intercept.br.MM <- list()
  br.MM <- list()
  br.predict <- list()
  
  for(ss in 1:CF){
    
    pos.br.MM[[ss]]       <- 1:length(A.MM[[ss]])
    Coef.br.MM[[ss]]      <- (1-A.MM[[ss]])
    Intercept.br.MM[[ss]] <- -(1-A.MM[[ss]])*Y.MM[[ss]]
    
    br.MM[[ss]] <- 
      FT_Minimax( Coef            = Coef.br.MM[[ss]][pos.br.MM[[ss]] ],
                  Intercept       = Intercept.br.MM[[ss]][pos.br.MM[[ss]] ],
                  X.Perturb       = YX.MM[[ss]][pos.br.MM[[ss]],],
                  X.Target        = WX.MM[[ss]][pos.br.MM[[ss]],],
                  bw.Perturb      = Opt.Para.br[[ss]][1],
                  bw.Target       = Opt.Para.br[[ss]][2],
                  lambda.Perturb  = Opt.Para.br[[ss]][3],
                  lambda.Target   = Opt.Para.br[[ss]][4],
                  SVM.bias = FALSE,
                  bias.input = NULL)
    
    br.predict[[ss]] <- function(WX.New.Input){
      FT_RBF(X     = WX.MM[[ss]][pos.br.MM[[ss]],],
             X.new = WX.New.Input,
             bw    = Opt.Para.br[[ss]][2])%*%br.MM[[ss]]$gamma + 
        br.MM[[ss]]$bias
    }
  }
  
  brhat.NoCF <- rep(0,N)
  for(ss in 1:CF){
    brhat.NoCF[SS.Index[[ss]]] <- br.predict[[ss]](WX.MM[[ss]])
  }
  
  brhat.CF <- rep(0,N)
  for(ss in 1:CF){
    brhat.CF[SS.Index[[3-ss]]] <- br.predict[[ss]](WX.MM[[3-ss]])
  }
  
  brhat.NoCF.MM <- list()
  for(ss in 1:CF){
    brhat.NoCF.MM[[ss]]  <- brhat.NoCF[SS.Index[[ss]] ]
  }
  
  Weight.CF  <- pihat.CF
  BFhat.CF  <- brhat.CF
  
  ncIF[[TYPE]] <- ((Y)*A - ( A*BFhat.CF + (1-A)*Weight.CF*(Y-BFhat.CF) ))*sY1
  
}

NCIF <- cbind(ncIF[[1]],ncIF[[2]],ncIF[[3]])
EFF  <- apply(NCIF,2,mean)/mean(A)
VAR  <- var((NCIF - matrix(A,N,3)*matrix(EFF,N,3,byrow=T))/mean(A))/N

RRR  <- c(BATCH, EFF, diag(VAR)^(0.5), as.vector(VAR))

RESULT <- matrix(RRR ,1, length(RRR))
colnames(RESULT) <- c("BATCH",
                      "Eff1","Eff2","Eff3",
                      "SE1","SE2","SE3",
                      "VAR11","VAR12","VAR13",
                      "VAR21","VAR22","VAR23",
                      "VAR31","VAR32","VAR33")
 
write.csv(RESULT,
          sprintf("NP/NP_CrossFitting/Result_BATCH%0.4d.csv",BATCH),
          row.names=F)
