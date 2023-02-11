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

source("SPC_Ft.R")

#############################
# Dataset
# Zika Virus in Brazil
# The dataset is used in Tchetgen Tchetgen, Park, Richardson (2023) Universal Difference-in-Differences for Causal Inference in Epidemiology
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
# Single Proxy Control
# Sensitivity
#############################

W.Origin <- Data$Y0
Y.Origin <- Data$Y1

alpha.W.Grid <- seq(-0.1,0.9,length=1001)
alpha.W <- alpha.W.Grid[BATCH]

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


Result.PS.W.Sen <- try(SPC.PS.Sensitivity(Y.Origin = Y.Origin,
                                           Trt = Trt,
                                           Ymat = Ymat.D1,
                                           Wmat = Wmat.2.D,
                                           Purturb = W.Origin*alpha.W,
                                           penalty=10^(-6)))

if(class(Result.PS.W.Sen)=="try-error"){
  Result.PS.W.Sen <- list()
  Result.PS.W.Sen$Est <- NA
  Result.PS.W.Sen$SE <- NA
}

Result <- matrix(c(num.bin,
                   alpha.W,
                   Result.PS.W.Sen$Est,
                   Result.PS.W.Sen$SE),1,4)

colnames(Result) <- c("Bin","alphaW","Original.Est","Original.SE")


write.csv( Result,
           sprintf("Result_Reg/Result_%0.5d.csv",BATCH),
           row.names=F)
