library(MASS)

expit <- function(v){exp(v)/(1+exp(v))}
link <- expit

########################################
# Moment Equation
########################################

GMM.Moment.PS <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  A <- data[,d.A.pos]
  Y <- data[,d.Y.orig.pos]
  Pi  <- link( as.matrix(data[,d.Y.pos])%*%theta[-c(1,2)] )
  IPW <- 1/(1-Pi)
  ResPS <- ((1-A)*IPW - 1)
  Res <- matrix(ResPS,length(ResPS),length(d.W.pos))*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- (1-A)*(Y)*(PSW) - (1-A)*(PSW)*theta[2]
  return( cbind(ATT1,ATT0,Res) )
}


GMM.Moment.PS.Sensitivity <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  d.Ptb.pos <- which( colnames(data)=="Ptb" )
  
  A <- data[,d.A.pos]
  Y <- data[,d.Y.orig.pos]
  Ptb <- data[,d.Ptb.pos]
  Pi  <- link( as.matrix(data[,d.Y.pos])%*%theta[-c(1,2)] + Ptb )
  IPW <- 1/(1-Pi)
  ResPS <- ((1-A)*IPW - 1)
  Res <- matrix(ResPS,length(ResPS),length(d.W.pos))*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- (1-A)*(Y)*(PSW) - (1-A)*(PSW)*theta[2]
  
  
  return( cbind(ATT1,ATT0,Res) )
}


GMM.Moment.OR <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  A <- as.numeric( data[,d.A.pos] )
  Y <- as.numeric( data[,d.Y.orig.pos] )
  Ym <- as.matrix( data[,d.Y.pos] )
  Wm <- as.matrix( data[,d.W.pos] )
  Res <- matrix( (1-A)*( Wm%*%theta[-c(1,2)] - Y ), length(Y), dim(Ym)[2] )*Ym
  
  Bridge <- Wm%*%theta[-c(1,2)]
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- A*(Bridge-theta[2])
  
  return( cbind(ATT1,ATT0,Res) )
}


GMM.Moment.DR <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.PY.pos <- which( substr(colnames(data),1,2)=="PY" )
  d.PW.pos <- which( substr(colnames(data),1,2)=="PW" )
  d.OY.pos <- which( substr(colnames(data),1,2)=="OY" )
  d.OW.pos <- which( substr(colnames(data),1,2)=="OW" )
  
  A <- as.numeric( data[,d.A.pos] )
  Y <- as.numeric( data[,d.Y.orig.pos] )
  PYm <- as.matrix( data[,d.PY.pos] )
  PWm <- as.matrix( data[,d.PW.pos] )
  OYm <- as.matrix( data[,d.OY.pos] )
  OWm <- as.matrix( data[,d.OW.pos] )
  
  theta.PS <- theta[2+1:dim(PYm)[2]]
  theta.OR <- theta[2+dim(PYm)[2]+1:dim(OWm)[2]]
  
  Pi  <- link( PYm%*%theta.PS )
  IPW <- 1/(1-Pi)
  PRes <- matrix(((1-A)*IPW - 1),length(Y),dim(PWm)[2])*(PWm)
  PSW <- Pi/(1-Pi)
  
  ORes <- matrix( (1-A)*( OWm%*%theta.OR - Y ), length(Y), dim(OYm)[2] )*(OYm)
  
  Bridge <- OWm%*%theta.OR
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- ((1-A)*(Y - Bridge)*PSW+(A)*Bridge) - A*theta[2]
  
  return( cbind(ATT1,ATT0,PRes,ORes) )
}


########################################
# Gradient
########################################

GMM.Moment.Grad.PS <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  A <- data[,d.A.pos]
  Y <- data[,d.Y.orig.pos]
  Pi  <- link( as.matrix(data[,d.Y.pos])%*%theta[-c(1,2)])
  IPW <- 1/(1-Pi)
  ResPS <- ((1-A)*IPW - 1)
  Res <- matrix(ResPS,length(ResPS),length(d.W.pos))*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- (1-A)*(PSW)*(Y-theta[2])
  
  
  GRAD.True <- matrix(0,length(theta),2+length(d.W.pos))
  GRAD.True[1,1]  <- -mean(A)
  GRAD.True[2,2]  <- -mean((1-A)*(PSW))
  
  GRAD.True[-c(1,2),2] <- as.numeric( apply(matrix((1-A)*PSW*(Y-theta[2]),dim(data)[1],length(d.Y.pos))*
                                              as.matrix(data[,d.Y.pos]),2,mean) )
  
  for(tt in 1:length(d.W.pos)){
    GRAD.True[-c(1,2),2+tt] <- 
      as.numeric( apply(matrix( (1-A)*PSW ,dim(data)[1],length(d.Y.pos))*
                          as.matrix(data[,d.Y.pos])*
                          matrix(data[,d.W.pos[tt]],dim(data)[1],length(d.Y.pos)),2,mean) )
  }
  
  
  return( GRAD.True )
  
  # dg_1(O,t)/dt_1 dg_2(O,t)/dt_1 dg_3(O,t)/dt_1 ...
  # dg_1(O,t)/dt_2 dg_2(O,t)/dt_2 dg_3(O,t)/dt_2 ...
  # dg_1(O,t)/dt_3 dg_2(O,t)/dt_3 dg_3(O,t)/dt_3 ...
  # (2+num.inner)*(2+dim.outer)
}

GMM.Moment.Grad.PS.Sensitivity <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  d.Ptb.pos <- which( colnames(data)=="Ptb" )
  
  A <- data[,d.A.pos]
  Y <- data[,d.Y.orig.pos]
  Ptb <- data[,d.Ptb.pos]
  Pi  <- link( as.matrix(data[,d.Y.pos])%*%theta[-c(1,2)] + Ptb )
  IPW <- 1/(1-Pi)
  ResPS <- ((1-A)*IPW - 1)
  Res <- matrix(ResPS,length(ResPS),length(d.W.pos))*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- (1-A)*(PSW)*(Y-theta[2])
  
  
  GRAD.True <- matrix(0,length(theta),2+length(d.W.pos))
  GRAD.True[1,1]  <- -mean(A)
  GRAD.True[2,2]  <- -mean((1-A)*(PSW))
  
  GRAD.True[-c(1,2),2] <- as.numeric( apply(matrix((1-A)*PSW*(Y-theta[2]),dim(data)[1],length(d.Y.pos))*
                                              as.matrix(data[,d.Y.pos]),2,mean) )
  
  for(tt in 1:length(d.W.pos)){
    GRAD.True[-c(1,2),2+tt] <- 
      as.numeric( apply(matrix( (1-A)*PSW ,dim(data)[1],length(d.Y.pos))*
                          as.matrix(data[,d.Y.pos])*
                          matrix(data[,d.W.pos[tt]],dim(data)[1],length(d.Y.pos)),2,mean) )
  }
  
  
  return( GRAD.True )
  
  # dg_1(O,t)/dt_1 dg_2(O,t)/dt_1 dg_3(O,t)/dt_1 ...
  # dg_1(O,t)/dt_2 dg_2(O,t)/dt_2 dg_3(O,t)/dt_2 ...
  # dg_1(O,t)/dt_3 dg_2(O,t)/dt_3 dg_3(O,t)/dt_3 ...
  # (2+num.inner)*(2+dim.outer)
}


GMM.Moment.Grad.OR <-  function(theta,data){
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  
  A <- as.numeric( data[,d.A.pos] )
  Y <- as.numeric( data[,d.Y.orig.pos] )
  Ym <- as.matrix( data[,d.Y.pos] )
  Wm <- as.matrix( data[,d.W.pos] )
  
  Res <- matrix( (1-A)*( Wm%*%theta[-c(1,2)] - Y ), length(Y), dim(Ym)[2] )*Ym
  
  Bridge <- Wm%*%theta[-c(1,2)]
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- A*(Bridge-theta[2])
  
  GRAD.True <- matrix(0,length(theta),2+length(d.Y.pos))
  GRAD.True[1,1]  <- -mean(A)
  GRAD.True[2,2]  <- -mean(A)
  
  GRAD.True[-c(1,2),2] <- as.numeric( apply(A*Wm,2,mean) )
  
  for(tt in 1:length(d.Y.pos)){
    GRAD.True[-c(1,2),2+tt] <- apply((1-A)*Wm*Ym[,tt],2,mean)
  }
  
  
  return( GRAD.True )
  
  # dg_1(O,t)/dt_1 dg_2(O,t)/dt_1 dg_3(O,t)/dt_1 ...
  # dg_1(O,t)/dt_2 dg_2(O,t)/dt_2 dg_3(O,t)/dt_2 ...
  # dg_1(O,t)/dt_3 dg_2(O,t)/dt_3 dg_3(O,t)/dt_3 ...
  # (2+num.inner)*(2+dim.outer)
}


GMM.Moment.Grad.DR <-  function(theta,data){
  
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.PY.pos <- which( substr(colnames(data),1,2)=="PY" )
  d.PW.pos <- which( substr(colnames(data),1,2)=="PW" )
  d.OY.pos <- which( substr(colnames(data),1,2)=="OY" )
  d.OW.pos <- which( substr(colnames(data),1,2)=="OW" )
  
  A <- as.numeric( data[,d.A.pos] )
  Y <- as.numeric( data[,d.Y.orig.pos] )
  PYm <- as.matrix( data[,d.PY.pos] )
  PWm <- as.matrix( data[,d.PW.pos] )
  OYm <- as.matrix( data[,d.OY.pos] )
  OWm <- as.matrix( data[,d.OW.pos] )
  
  theta.PS <- theta[2+1:dim(PYm)[2]]
  theta.OR <- theta[2+dim(PYm)[2]+1:dim(OWm)[2]]
  
  Pi  <- link( PYm%*%theta.PS )
  IPW <- 1/(1-Pi)
  PRes <- matrix(((1-A)*IPW - 1),length(Y),dim(PWm)[2])*(PWm)
  PSW <- Pi/(1-Pi)
  
  ORes <- matrix( (1-A)*( OWm%*%theta.OR - Y ), length(Y), dim(OYm)[2] )*(OYm)
  
  Bridge <- OWm%*%theta.OR
  
  ATT1 <- A*(Y-theta[1])
  ATT0 <- ((1-A)*(Y - Bridge)*PSW+(A)*Bridge) - A*theta[2]
  
  GRAD.True <- matrix(0,length(theta),2+length(d.PW.pos)+length(d.OY.pos))
  GRAD.True[1,1]  <- -mean(A)
  GRAD.True[2,2]  <- -mean(A)
  
  GRAD.True[2+1:length(d.PY.pos),2] <- 
    as.numeric( apply(matrix((1-A)*PSW*(Y-Bridge),dim(data)[1],length(d.PY.pos))*
                        as.matrix(data[,d.PY.pos]),2,mean) )
  GRAD.True[2+length(d.PY.pos)+1:length(d.OW.pos),2] <- apply( matrix( (1-A)*(-PSW)+A ,dim(data)[1],dim(OWm)[2])*OWm,2,mean)
  
  
  for(tt in 1:length(d.PW.pos)){
    GRAD.True[2+1:length(d.PY.pos),2+tt] <- 
      as.numeric( apply(matrix( (1-A)*PSW ,dim(data)[1],length(d.PY.pos))*
                          as.matrix(data[,d.PY.pos])*
                          matrix(data[,d.PW.pos[tt]],dim(data)[1],length(d.PY.pos)),2,mean) )
  }
  
  
  for(tt in 1:length(d.OY.pos)){
    GRAD.True[2+length(d.PY.pos)+1:length(d.OW.pos),2+length(d.PW.pos)+tt] <- apply((1-A)*OWm*OYm[,tt],2,mean)
  }
  
  
  
  return( GRAD.True )
  
  # dg_1(O,t)/dt_1 dg_2(O,t)/dt_1 dg_3(O,t)/dt_1 ...
  # dg_1(O,t)/dt_2 dg_2(O,t)/dt_2 dg_3(O,t)/dt_2 ...
  # dg_1(O,t)/dt_3 dg_2(O,t)/dt_3 dg_3(O,t)/dt_3 ...
  # (2+num.inner)*(2+dim.outer)
}


########################################
# Variance Matrix
########################################

MEAT <- function(theta,data,TT="PS"){
  
  if(TT=="PS"){
    RESIDUALS <- GMM.Moment.PS(theta,as.matrix(data))
    Sigma <- t(RESIDUALS)%*%(RESIDUALS)/dim(RESIDUALS)[1]
  } else if (TT=="PS.Sensitivity"){
    RESIDUALS <- GMM.Moment.PS.Sensitivity(theta,as.matrix(data))
    Sigma <- t(RESIDUALS)%*%(RESIDUALS)/dim(RESIDUALS)[1]
  } else if (TT=="OR"){
    RESIDUALS <- GMM.Moment.OR(theta,as.matrix(data))
    Sigma <- t(RESIDUALS)%*%(RESIDUALS)/dim(RESIDUALS)[1]
  } else if (TT=="DR"){
    RESIDUALS <- GMM.Moment.DR(theta,as.matrix(data))
    Sigma <- t(RESIDUALS)%*%(RESIDUALS)/dim(RESIDUALS)[1]
  }
  
  
  return(Sigma)
}

VMAT.GMM.Exact <- function(theta,data,TT="PS",Var.Type="Eff",penalty.vec){
  
  if(TT=="PS"){
    GMM.Moment <- GMM.Moment.PS
    GMM.Grad   <- GMM.Moment.Grad.PS
  } else if (TT=="PS.Sensitivity"){
    GMM.Moment <- GMM.Moment.PS.Sensitivity
    GMM.Grad   <- GMM.Moment.Grad.PS.Sensitivity
  } else if (TT=="OR"){
    GMM.Moment <- GMM.Moment.OR
    GMM.Grad   <- GMM.Moment.Grad.OR
  } else if (TT=="DR"){
    GMM.Moment <- GMM.Moment.DR
    GMM.Grad   <- GMM.Moment.Grad.DR
  }
  
  RESIDUALS <- GMM.Moment(theta,as.matrix(data))
  Sigma <- t(RESIDUALS)%*%(RESIDUALS)/dim(RESIDUALS)[1]
  if(Var.Type=="Eff"){
    WEIGHT <- ginv(Sigma)
  } else {
    WEIGHT <- diag(rep(1,dim(RESIDUALS)[2]))
  }
  
  Grad.Mat <- GMM.Grad(theta,data)
  
  Stable <- diag(penalty.vec)
  StableSq <- diag(penalty.vec^2)
  
  if(Var.Type=="Eff"){
    VV <- try(ginv(Grad.Mat%*%WEIGHT%*%t(Grad.Mat)),silent=T)
    VV2 <- VV%*%(Grad.Mat%*%WEIGHT)
    OPT <- VV2%*%(Sigma)%*%t(VV2)
  } else {
    # OPT <- try(ginv(Grad.Mat%*%WEIGHT%*%t(Grad.Mat)) %*% (
    #   Grad.Mat%*%(WEIGHT)%*%(Sigma)%*%(WEIGHT)%*%t(Grad.Mat)
    # ) %*% ginv(Grad.Mat%*%WEIGHT%*%t(Grad.Mat)),silent=T)
    OPT <- try(ginv(Grad.Mat%*%WEIGHT%*%t(Grad.Mat)),silent=T)
  }
  return(OPT)
}

########################################
# Optim
########################################

SOLVE.GMM <- function(FN.List,
                      PS,
                      data,
                      Num.NN,
                      radius,
                      Var.Type,
                      Optim.Method,
                      penalty.vec){
  GMM.Moment.Current <- FN.List$GMM.Moment.Current
  GMM.Grad.Current <- FN.List$GMM.Grad.Current
  VMAT.Current <- FN.List$VMAT.Current
  MEAT.Current <- FN.List$MEAT.Current
  
  if(Var.Type=="Eff"){
    ########### 1st stage opt.
    
    if(Num.NN > 0){
      Candidates <- matrix(0,Num.NN+1,length(PS)+1)
      for(bb in 1:Num.NN){
        
        set.seed(bb)
        
        CPS <- Candidates[bb,1+1:length(PS)] <- PS+runif(length(PS))*radius-0.5*radius
        
        Candidates[bb,1] <- optim(par=CPS ,fn=function(theta){
          t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
        },
        gr=function(theta){
          2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
        },method=Optim.Method)$value
        
      }
      bb <- Num.NN+1
      Candidates[bb,1+1:length(PS)] <- PS
      Candidates[bb,1] <- optim(par=CPS ,fn=function(theta){
        t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
      },
      gr=function(theta){
        2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
      },method=Optim.Method)$value
      
      PS <- Candidates[which.min(Candidates[,1]),-1]
    } else {
      PS <- optim(par=PS ,fn=function(theta){
        t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
      },
      gr=function(theta){
        2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
      },method=Optim.Method)$par
    }
    
    ########### 1st stage opt. (end)
    
    ########### 2nd stage opt.
    
    OPT <- VMAT.Current(PS, data)
    PS <- optim(par=PS ,fn=function(theta){
      t(apply(GMM.Moment.Current(theta,data),2,mean))%*%
        (ginv( MEAT.Current(PS,data) ))%*%
        (apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
    },
    gr=function(theta){
      2*GMM.Grad.Current(theta,data)%*%
        (ginv( MEAT.Current(PS,data) ))%*%
        apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
    },method=Optim.Method)$par
    
  } else {
    ## Inefficient, but more stable
    ########### 1st stage opt.
    
    if(Num.NN > 0){
      
      Candidates <- matrix(0,Num.NN+1,length(PS)+1)
      for(bb in 1:Num.NN){
        
        set.seed(bb)
        
        CPS <- Candidates[bb,1+1:length(PS)] <- PS+runif(length(PS))*radius-0.5*radius
        
        Candidates[bb,1] <- optim(par=CPS ,fn=function(theta){
          t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
        },
        gr=function(theta){
          2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
        },method=Optim.Method)$value
        
      }
      bb <- Num.NN+1
      Candidates[bb,1+1:length(PS)] <- PS
      Candidates[bb,1] <- optim(par=CPS ,fn=function(theta){
        t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
      },
      gr=function(theta){
        2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
      },method=Optim.Method)$value
      
      PS <- Candidates[which.min(Candidates[,1]),-1]
    } else {
      
      PS <- optim(par=PS ,fn=function(theta){
        t(apply(GMM.Moment.Current(theta,data),2,mean))%*%(apply(GMM.Moment.Current(theta,data),2,mean)) + sum((theta*penalty.vec)^2)
      },
      gr=function(theta){
        2*GMM.Grad.Current(theta,data)%*%apply(GMM.Moment.Current(theta,data),2,mean) + 2*penalty.vec*theta
      },method=Optim.Method)$par
    }
  }
  
  
  return(PS)
  
}



########################################
# Estimations
########################################

# Ymat=Ymat.PS[[1]]
# Wmat=Wmat.2.D
# penalty=10^(-5)
# starting=STARTING.PS[[ 1 ]]
# Var.Type = "Eff"
# Num.NN=0
# radius=1
# Optim.Method="Nelder-Mead"
# 

SPC.PS <- function(Y.Origin,Trt,
                   Ymat,Wmat,
                   Num.NN=0,
                   radius=1,
                   Var.Type="Eff",
                   starting=NULL,
                   Optim.Method="Nelder-Mead",
                   penalty=10^(-5)){
  
  if(is.null(dim(Ymat))){
    Y <- matrix(Ymat,length(Ymat),1)
  } else {
    Y <- Ymat
  }
  if(is.null(dim(Wmat))){
    W <- matrix(Wmat,length(Wmat),1)
  } else {
    W <- Wmat
  }
  
  data <- data.frame(cbind(Y.Origin,Trt,
                           Y,W))
  colnames(data) <- c( "Outcome",
                       "Trt",
                       sprintf("Y_%0.3d",1:dim(Y)[2]),
                       sprintf("W_%0.3d",1:dim(W)[2]))
  
  
  # COEF <- glm(Trt~0+Y,data=data,family="binomial")$coefficient
  COEF <- starting
  if(is.null(COEF)){
    COEF <- rep(0,dim(Y)[2])
  } 
  penalty.vec <- rep(0,2+length(COEF))
  penalty.vec[2+which(apply(Y,2,mean)!=1)] <- sqrt(penalty)
  
  VMAT.Current <- function(theta,data){
    VMAT.GMM.Exact(theta, data , TT="PS", Var.Type=Var.Type, penalty.vec )
  }
  MEAT.Current <- function(theta,data){
    MEAT(theta, data , TT="PS" )
  }
  GMM.Moment.Current <- function(theta,data){
    GMM.Moment.PS(theta,data)
  }
  GMM.Grad.Current <- function(theta,data){
    GMM.Moment.Grad.PS(theta,data)
  }
  
  GRID <- function(vv){
    d.Y.orig.pos <- which(colnames(data)=="Outcome")
    d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
    d.W.pos <- which( substr(colnames(data),1,1)=="W" )
    Pi  <- link( as.matrix(data[,d.Y.pos])%*%vv )
    IPW <- 1/(1-Pi)
    Res <- ((1-data$Trt)*IPW - 1)*(data[,d.W.pos])
    sum( ( apply(Res,2,mean) )^2 + sum((vv*penalty.vec[2+length(COEF)])^2) )
  }
  
  OPTIM1 <- matrix(0,Num.NN+1,length(COEF)+1)
  bb <- Num.NN+1
  OPT <- optim(par=COEF,
               fn=GRID,
               method=Optim.Method)
  OPTIM1[bb,] <- c(OPT$value,OPT$par)
  
  if(Num.NN>0){
    for(bb in 1:Num.NN){
      set.seed(bb)
      OPT <- optim(par=COEF+runif(length(COEF))*radius-0.5*radius,
                   fn=GRID,
                   method=Optim.Method)
      OPTIM1[bb,] <- c(OPT$value,OPT$par)
    }
  } else {
    OPTIM1 <- rbind(OPTIM1,
                    OPTIM1)
  }
  PS <- OPTIM1[which.min(OPTIM1[,1]),-1]
  
  
  

  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  Pi  <- expit( as.matrix(data[,d.Y.pos])%*%PS )
  IPW <- 1/(1-Pi)
  Res <- ((1-data$Trt)*IPW - 1)*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  ATT1 <- mean( data$Outcome[data$Trt==1] )
  ATT0 <- mean((1-data$Trt)*(data$Outcome)*(PSW))/mean((1-data$Trt)*(PSW))
  ATT <- ATT1-ATT0
  
  ParaStart <- c(ATT,ATT1,ATT0,PS)
  PS <- as.numeric(ParaStart[-1])
  OPT <- VMAT.Current(PS, data)
  Sigma <- MEAT.Current(PS, data) 
  
  
  FN.List <- list()
  FN.List$GMM.Moment.Current <- GMM.Moment.Current
  FN.List$GMM.Grad.Current <- GMM.Grad.Current
  FN.List$VMAT.Current <- VMAT.Current
  FN.List$MEAT.Current <- MEAT.Current
  
  
  GMM <- list()
  GMM$coefficients <- SOLVE.GMM(FN.List,PS,data,Num.NN,radius,
                                Var.Type=Var.Type,Optim.Method=Optim.Method,
                                penalty.vec = penalty.vec)
  GMM$vcov <- VMAT.Current(as.numeric(GMM$coefficients), data)/dim(data)[1]
  
  CONT <- rbind( c(1,-1,rep(0,length(GMM$coefficients)-2)), 
                 diag(rep(1,length(GMM$coefficients))) )
  
  CEE.Full <- c(CONT%*%GMM$coefficients)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("Y_%0.3d",1:dim(Y)[2]))
  
  CVar.Full <- CONT%*%as.matrix(GMM$vcov)%*%t(CONT)
  
  Result <- list()
  Result$Est <- as.numeric( CEE.Full[1] )
  Result$SE  <- as.numeric( sqrt(CVar.Full[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}

SPC.PS.Sensitivity <- function(Y.Origin,Trt,
                               Ymat,Wmat,
                               Purturb,
                               Num.NN=0,
                               radius=1,
                               Var.Type="Eff",
                               starting=NULL,
                               Optim.Method="Nelder-Mead",
                               penalty=10^(-5)){
  
  if(is.null(dim(Ymat))){
    Y <- matrix(Ymat,length(Ymat),1)
  } else {
    Y <- Ymat
  }
  if(is.null(dim(Wmat))){
    W <- matrix(Wmat,length(Wmat),1)
  } else {
    W <- Wmat
  }
  
  data <- data.frame(cbind(Y.Origin,Trt,
                           Y,W,
                           Purturb))
  colnames(data) <- c( "Outcome",
                       "Trt",
                       sprintf("Y_%0.3d",1:dim(Y)[2]),
                       sprintf("W_%0.3d",1:dim(W)[2]),
                       "Ptb")
  
  COEF <- starting
  if(is.null(COEF)){
    COEF <- rep(0,dim(Y)[2])
  } 
  penalty.vec <- rep(0,2+length(COEF))
  penalty.vec[2+which(apply(Y,2,mean)!=1)] <- sqrt(penalty) 
  
  VMAT.Current <- function(theta,data){
    VMAT.GMM.Exact(theta, data , TT="PS.Sensitivity", Var.Type=Var.Type, penalty.vec )
  }
  MEAT.Current <- function(theta,data){
    MEAT(theta, data , TT="PS.Sensitivity" )
  }
  GMM.Moment.Current <- function(theta,data){
    GMM.Moment.PS.Sensitivity(theta,data)
  }
  GMM.Grad.Current <- function(theta,data){
    GMM.Moment.Grad.PS.Sensitivity(theta,data)
  }
  
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  d.Ptb.pos <- which( colnames(data)=="Ptb" )
  
  GRID <- function(vv){
    d.Y.orig.pos <- which(colnames(data)=="Outcome")
    d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
    d.W.pos <- which( substr(colnames(data),1,1)=="W" )
    d.Ptb.pos <- which( colnames(data)=="Ptb" )
    Pi  <- link( as.matrix(data[,d.Y.pos])%*%vv + data[,d.Ptb.pos] )
    IPW <- 1/(1-Pi)
    Res <- ((1-data$Trt)*IPW - 1)*(data[,d.W.pos])
    sum( ( apply(Res,2,mean) )^2 + sum((vv*penalty.vec[2+length(COEF)])^2) )
  }
  
  OPTIM1 <- matrix(0,Num.NN+1,length(COEF)+1)
  bb <- Num.NN+1
  OPT <- optim(par=COEF,
               fn=GRID,
               method=Optim.Method)
  OPTIM1[bb,] <- c(OPT$value,OPT$par)
  
  if(Num.NN>0){
    for(bb in 1:Num.NN){
      set.seed(bb)
      OPT <- optim(par=COEF+runif(length(COEF))*radius-0.5*radius,
                   fn=GRID,
                   method=Optim.Method)
      OPTIM1[bb,] <- c(OPT$value,OPT$par)
    }
  } else {
    OPTIM1 <- rbind(OPTIM1,
                    OPTIM1)
  }
  PS <- OPTIM1[which.min(OPTIM1[,1]),-1]
  
  Pi  <- expit( as.matrix(data[,d.Y.pos])%*%PS + data[,d.Ptb.pos] )
  IPW <- 1/(1-Pi)
  Res <- ((1-data$Trt)*IPW - 1)*(data[,d.W.pos])
  PSW <- Pi/(1-Pi)
  ATT1 <- mean( data$Outcome[data$Trt==1] )
  ATT0 <- mean((1-data$Trt)*(data$Outcome)*(PSW))/mean((1-data$Trt)*(PSW))
  ATT <- ATT1-ATT0
  
  ParaStart <- c(ATT,ATT1,ATT0,PS)
  PS <- as.numeric(ParaStart[-1])
  OPT <- VMAT.Current(PS, data)
  Sigma <- MEAT.Current(PS, data) 
  
  
  FN.List <- list()
  FN.List$GMM.Moment.Current <- GMM.Moment.Current
  FN.List$GMM.Grad.Current <- GMM.Grad.Current
  FN.List$VMAT.Current <- VMAT.Current
  FN.List$MEAT.Current <- MEAT.Current
  
  
  GMM <- list()
  GMM$coefficients <- SOLVE.GMM(FN.List,PS,data,Num.NN,radius,Var.Type=Var.Type,Optim.Method=Optim.Method,penalty.vec = penalty.vec)
  GMM$vcov <- VMAT.Current(as.numeric(GMM$coefficients), data)/dim(data)[1]
  
  CONT <- rbind( c(1,-1,rep(0,length(GMM$coefficients)-2)), 
                 diag(rep(1,length(GMM$coefficients))) )
  
  CEE.Full <- c(CONT%*%GMM$coefficients)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("Y_%0.3d",1:dim(Y)[2]))
  
  CVar.Full <- CONT%*%as.matrix(GMM$vcov)%*%t(CONT)
  
  
  
  Result <- list()
  Result$Est <- as.numeric( CEE.Full[1] )
  Result$SE  <- as.numeric( sqrt(CVar.Full[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}


# Ymat=Ymat.2.D
# Wmat=Wmat.OR[[ CV.OR[2] ]]
# penalty=10^(CV.OR[1])
# starting=STARTING.OR[[ CV.OR[2] ]]
# Var.Type = "Eff"
# Optim.Method="Nelder-Mead"
# Num.NN=NUMNN


SPC.OR <- function(Y.Origin,Trt,
                   Ymat,Wmat,
                   Num.NN=0,
                   radius=1,
                   Var.Type="Eff",
                   starting=NULL,
                   Optim.Method="Nelder-Mead",
                   penalty=0){
  
  if(is.null(dim(Ymat))){
    Y <- matrix(Ymat,length(Ymat),1)
  } else {
    Y <- Ymat
  }
  if(is.null(dim(Wmat))){
    W <- matrix(Wmat,length(Wmat),1)
  } else {
    W <- Wmat
  }
  
  
  data <- data.frame(cbind(Y.Origin,Trt,
                           Y,W))
  colnames(data) <- c( "Outcome",
                       "Trt",
                       sprintf("Y_%0.3d",1:dim(Y)[2]),
                       sprintf("W_%0.3d",1:dim(W)[2]))
  
  
  COEF <- starting
  if(is.null(COEF)){
    COEF <- lm(Y.Origin~0+W)$coefficients
  } 
  
  penalty.vec <- rep(0,2+length(COEF))
  penalty.vec <- rep(0,2+length(COEF))
  penalty.vec[2+which(apply(W,2,mean)!=1)] <- sqrt(penalty)
  
  VMAT.Current <- function(theta,data){
    VMAT.GMM.Exact(theta, data , TT="OR", Var.Type=Var.Type, penalty.vec )
  }
  MEAT.Current <- function(theta,data){
    MEAT(theta, data , TT="OR" )
  }
  GMM.Moment.Current <- function(theta,data){
    GMM.Moment.OR(theta,data)
  }
  GMM.Grad.Current <- function(theta,data){
    GMM.Moment.Grad.OR(theta,data)
  }
  
  GRID <- function(vv){
    d.Y.orig.pos <- which(colnames(data)=="Outcome")
    d.A.pos <- which(colnames(data)=="Trt")
    d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
    d.W.pos <- which( substr(colnames(data),1,1)=="W" )
    
    A <- as.numeric( data[,d.A.pos] )
    Y <- as.numeric( data[,d.Y.orig.pos] )
    Ym <- as.matrix( data[,d.Y.pos] )
    Wm <- as.matrix( data[,d.W.pos] )
    Res <- matrix( (1-A)*( Wm%*%vv - Y ), length(Y), dim(Ym)[2] )*Ym
    sum( ( apply(Res,2,mean) )^2 )
  }
  
  OPTIM1 <- matrix(0,Num.NN+1,length(COEF)+1)
  bb <- Num.NN+1
  OPT <- optim(par=COEF,
               fn=GRID,
               method=Optim.Method)
  OPTIM1[bb,] <- c(OPT$value,OPT$par)
  
  if(Num.NN>0){
    for(bb in 1:Num.NN){
      set.seed(bb)
      OPT <- optim(par=COEF+runif(length(COEF))*radius-0.5*radius,
                   fn=GRID,
                   method=Optim.Method)
      OPTIM1[bb,] <- c(OPT$value,OPT$par)
    }
  } else {
    OPTIM1 <- rbind(OPTIM1,
                    OPTIM1)
  }
  PS <- OPTIM1[which.min(OPTIM1[,1]),-1]
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
  d.W.pos <- which( substr(colnames(data),1,1)=="W" )
  
  Bridge <- Wmat%*%PS
  
  ATT1 <- mean( data$Outcome[data$Trt==1] )
  ATT0 <- mean( Bridge[data$Trt==1] )
  ATT <- ATT1-ATT0
  
  ParaStart <- c(ATT,ATT1,ATT0,PS)
  PS <- as.numeric(ParaStart[-1])
  OPT <- VMAT.Current(PS, data)
  Sigma <- MEAT.Current(PS, data) 
  
  FN.List <- list()
  FN.List$GMM.Moment.Current <- GMM.Moment.Current
  FN.List$GMM.Grad.Current <- GMM.Grad.Current
  FN.List$VMAT.Current <- VMAT.Current
  FN.List$MEAT.Current <- MEAT.Current
  
  GMM <- list()
  GMM$coefficients <- SOLVE.GMM(FN.List,PS,data,Num.NN,radius,Var.Type=Var.Type,Optim.Method=Optim.Method,penalty.vec = penalty.vec)
  GMM$vcov <- VMAT.Current(as.numeric(GMM$coefficients), data)/dim(data)[1]
  
  CONT <- rbind( c(1,-1,rep(0,length(GMM$coefficients)-2)), 
                 diag(rep(1,length(GMM$coefficients))) )
  
  CEE.Full <- c(CONT%*%GMM$coefficients)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("W_%0.3d",1:dim(W)[2]))
  
  CVar.Full <- CONT%*%as.matrix(GMM$vcov)%*%t(CONT)
  
  Result <- list()
  Result$Est <- as.numeric( CEE.Full[1] )
  Result$SE  <- as.numeric( sqrt(CVar.Full[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
}






# Ymat.PS <- Ymat
# Wmat.PS <- Wmat
# Ymat.OR <- Ymat
# Wmat.OR <- Wmat
# Num.NN=0
# radius=1
# Var.Type="Ineff"


SPC.DR <- function(Y.Origin,Trt,
                   Ymat.PS,Wmat.PS,
                   Ymat.OR,Wmat.OR,
                   Num.NN=0,
                   radius=1,
                   Var.Type="Eff",
                   Optim.Method="Nelder-Mead",
                   penalty=10^(-5),
                   COEF.PS=NULL,
                   COEF.OR=NULL){
  
  if(is.null(dim(Ymat.PS))){
    Y.PS <- matrix(Ymat.PS,length(Ymat.PS),1)
  } else {
    Y.PS <- Ymat.PS
  }
  if(is.null(dim(Ymat.OR))){
    Y.OR <- matrix(Ymat.OR,length(Ymat.OR),1)
  } else {
    Y.OR <- Ymat.OR
  }
  if(is.null(dim(Wmat.PS))){
    W.PS <- matrix(Wmat.PS,length(Wmat.PS),1)
  } else {
    W.PS <- Wmat.PS
  }
  if(is.null(dim(Wmat.OR))){
    W.OR <- matrix(Wmat.OR,length(Wmat.OR),1)
  } else {
    W.OR <- Wmat.OR
  }
  data <- data.frame(cbind(Y.Origin,Trt,
                           Y.PS,W.PS,
                           Y.OR,W.OR))
  colnames(data) <- c( "Outcome",
                       "Trt",
                       sprintf("PY_%0.3d",1:dim(Y.PS)[2]),
                       sprintf("PW_%0.3d",1:dim(W.PS)[2]),
                       sprintf("OY_%0.3d",1:dim(Y.OR)[2]),
                       sprintf("OW_%0.3d",1:dim(W.OR)[2]))
  
  if(is.null(COEF.PS)){
    # COEF.PS <- glm(Trt~0+Y.PS,data=data,family="binomial")$coefficient
    COEF.PS <- rep(0,dim(Y.PS)[2])
  }
  if(is.null(COEF.OR)){
    COEF.OR <- lm(Y.Origin~0+W.OR)$coefficients
  }
  
  
  penalty.vec <- rep(0,2+length(COEF.PS)+length(COEF.OR))
  penalty.vec[2+which(apply(cbind(Ymat.PS,Wmat.OR),2,mean)!=1)] <- sqrt(penalty) 
  
  VMAT.Current <- function(theta,data){
    VMAT.GMM.Exact(theta, data , TT="DR", Var.Type=Var.Type, penalty.vec )
  }
  MEAT.Current <- function(theta,data){
    MEAT(theta, data , TT="DR" )
  }
  GMM.Moment.Current <- function(theta,data){
    GMM.Moment.DR(theta,data)
  }
  GMM.Grad.Current <- function(theta,data){
    GMM.Moment.Grad.DR(theta,data)
  }
  
  GRID.PS <- function(vv){
    d.Y.orig.pos <- which(colnames(data)=="Outcome")
    d.A.pos <- which(colnames(data)=="Trt")
    d.PY.pos <- which( substr(colnames(data),1,2)=="PY" )
    d.PW.pos <- which( substr(colnames(data),1,2)=="PW" )
    
    A <- as.numeric( data[,d.A.pos] )
    Y <- as.numeric( data[,d.Y.orig.pos] )
    Ym <- as.matrix( data[,d.PY.pos] )
    Wm <- as.matrix( data[,d.PW.pos] )
    
    Pi  <- link( Ym%*%vv )
    IPW <- 1/(1-Pi)
    Res <- matrix( ((1-A)*IPW - 1), length(A), dim(Wm)[2] )*(Wm)
    sum( ( apply(Res,2,mean) )^2 + sum((vv*penalty.vec[2+length(COEF.PS)])^2) )
  }
  
  GRID.OR <- function(vv){
    d.Y.orig.pos <- which(colnames(data)=="Outcome")
    d.A.pos <- which(colnames(data)=="Trt")
    d.OY.pos <- which( substr(colnames(data),1,2)=="OY" )
    d.OW.pos <- which( substr(colnames(data),1,2)=="OW" )
    
    A <- as.numeric( data[,d.A.pos] )
    Y <- as.numeric( data[,d.Y.orig.pos] )
    Ym <- as.matrix( data[,d.OY.pos] )
    Wm <- as.matrix( data[,d.OW.pos] )
    
    Res <- matrix( (1-A)*( Wm%*%vv - Y ), length(Y), dim(Ym)[2] )*Ym
    sum( ( apply(Res,2,mean) )^2 )
  }
  
  OPTIM1 <- matrix(0,Num.NN+1,length(COEF.PS)+1)
  bb <- Num.NN+1
  OPT <- optim(par=COEF.PS,
               fn=GRID.PS,
               method=Optim.Method)
  OPTIM1[bb,] <- c(OPT$value,OPT$par)
  
  if(Num.NN>0){
    for(bb in 1:Num.NN){
      set.seed(bb)
      OPT <- optim(par=COEF.PS+runif(length(COEF.PS))*radius-0.5*radius,
                   fn=GRID.PS,
                   method=Optim.Method)
      OPTIM1[bb,] <- c(OPT$value,OPT$par)
    }
  } else {
    OPTIM1 <- rbind(OPTIM1,
                    OPTIM1)
  }
  PS.PS <- OPTIM1[which.min(OPTIM1[,1]),-1]
  
  
  OPTIM2 <- matrix(0,Num.NN+1,length(COEF.OR)+1)
  bb <- Num.NN+1
  OPT <- optim(par=COEF.OR,
               fn=GRID.OR,
               method=Optim.Method)
  OPTIM2[bb,] <- c(OPT$value,OPT$par)
  
  if(Num.NN>0){
    for(bb in 1:Num.NN){
      set.seed(bb)
      OPT <- optim(par=COEF.OR+runif(length(COEF.OR))*radius-0.5*radius,
                   fn=GRID.OR,
                   method=Optim.Method)
      OPTIM2[bb,] <- c(OPT$value,OPT$par)
    }
  } else {
    OPTIM2 <- rbind(OPTIM2,
                    OPTIM2)
  }
  PS.OR <- OPTIM2[which.min(OPTIM1[,2]),-1]
  
  
  d.Y.orig.pos <- which(colnames(data)=="Outcome")
  d.A.pos <- which(colnames(data)=="Trt")
  d.PY.pos <- which( substr(colnames(data),1,2)=="PY" )
  d.PW.pos <- which( substr(colnames(data),1,2)=="PW" )
  d.OY.pos <- which( substr(colnames(data),1,2)=="OY" )
  d.OW.pos <- which( substr(colnames(data),1,2)=="OW" )
  
  Pi  <- link( as.matrix(data[,d.PY.pos])%*%PS.PS )
  IPW <- 1/(1-Pi)
  Res <- ((1-data$Trt)*IPW - 1)*(data[,d.PW.pos])
  PSW <- Pi/(1-Pi)
  
  Bridge <- as.matrix(data[,d.OW.pos])%*%PS.OR
  
  ATT1 <- mean( data$Outcome[data$Trt==1] )
  ATT0 <- mean( ((1-data$Trt)*(data$Outcome - Bridge)*PSW+(data$Trt)*Bridge) )/(mean(data$Trt))
  ATT <- ATT1-ATT0
  
  ParaStart <- c(ATT,ATT1,ATT0,PS.PS,PS.OR)
  PS <- as.numeric(ParaStart[-1])
  OPT <- VMAT.Current(PS, data)
  Sigma <- MEAT.Current(PS, data) 
  
  FN.List <- list()
  FN.List$GMM.Moment.Current <- GMM.Moment.Current
  FN.List$GMM.Grad.Current <- GMM.Grad.Current
  FN.List$VMAT.Current <- VMAT.Current
  FN.List$MEAT.Current <- MEAT.Current
  
  GMM <- list()
  GMM$coefficients <- SOLVE.GMM(FN.List,PS,data,Num.NN,radius,Var.Type=Var.Type,
                                Optim.Method=Optim.Method,penalty.vec = penalty.vec)
  
  GMM$vcov <- VMAT.Current(as.numeric(GMM$coefficients), data)/dim(data)[1]
  
  CONT <- rbind( c(1,-1,rep(0,length(GMM$coefficients)-2)), 
                 diag(rep(1,length(GMM$coefficients))) )
  
  CEE.Full <- c(CONT%*%GMM$coefficients)
  names(CEE.Full) <- c("ATT",
                       "ATT1",
                       "ATT0",
                       sprintf("PY_%0.3d",1:dim(Y.PS)[2]),
                       sprintf("OW_%0.3d",1:dim(W.OR)[2]))
  
  CVar.Full <- CONT%*%as.matrix(GMM$vcov)%*%t(CONT)
  
  
  Result <- list()
  Result$Est <- as.numeric( CEE.Full[1] )
  Result$SE  <- as.numeric( sqrt(CVar.Full[1,1]) )
  Result$Full.Est <- CEE.Full
  Result$Full.Var <- CVar.Full
  
  return(Result)
  
  
}







GMM.Moment.Grad.Merge <-  function(theta.collect,
                                   data.collect,
                                   type="PS"){
  
  GT <- list()
  for(TYPE in 1:3){
    
    theta <- theta.collect[[TYPE]]
    data  <- data.collect[[TYPE]]
    
    eps <- 10^(-8)
    
    if(type=="PS"){
      d.Y.orig.pos <- which(colnames(data)=="Outcome")
      d.A.pos <- which(colnames(data)=="Trt")
      d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
      d.W.pos <- which( substr(colnames(data),1,1)=="W" )
      
      GRAD.True <- matrix(0,length(theta),2+length(d.W.pos))
      for(tt in 1:length(theta)){
        theta.pos <- theta.neg <- theta
        theta.pos[tt] <- theta.pos[tt] + eps
        theta.neg[tt] <- theta.neg[tt] - eps
        GRAD.True[tt,] <- apply((GMM.Moment.PS(theta.pos,data) - 
                                   GMM.Moment.PS(theta.neg,data)),2,mean)/(2*eps)
      }
    } else if (type=="OR"){
      d.Y.orig.pos <- which(colnames(data)=="Outcome")
      d.A.pos <- which(colnames(data)=="Trt")
      d.Y.pos <- which( substr(colnames(data),1,1)=="Y" )
      d.W.pos <- which( substr(colnames(data),1,1)=="W" )
      
      GRAD.True <- matrix(0,length(theta),2+length(d.Y.pos))
      for(tt in 1:length(theta)){
        theta.pos <- theta.neg <- theta
        theta.pos[tt] <- theta.pos[tt] + eps
        theta.neg[tt] <- theta.neg[tt] - eps
        GRAD.True[tt,] <- apply((GMM.Moment.OR(theta.pos,data) - 
                                   GMM.Moment.OR(theta.neg,data)),2,mean)/(2*eps)
      }
    } else if (type=="DR"){
      d.Y.orig.pos <- which(colnames(data)=="Outcome")
      d.A.pos <- which(colnames(data)=="Trt")
      d.PY.pos <- which( substr(colnames(data),1,2)=="PY" )
      d.PW.pos <- which( substr(colnames(data),1,2)=="PW" )
      d.OY.pos <- which( substr(colnames(data),1,2)=="OY" )
      d.OW.pos <- which( substr(colnames(data),1,2)=="OW" )
      
      GRAD.True <- matrix(0,length(theta),2+length(d.PW.pos)+length(d.OY.pos))
      for(tt in 1:length(theta)){
        theta.pos <- theta.neg <- theta
        theta.pos[tt] <- theta.pos[tt] + eps
        theta.neg[tt] <- theta.neg[tt] - eps
        GRAD.True[tt,] <- apply((GMM.Moment.DR(theta.pos,data) - 
                                   GMM.Moment.DR(theta.neg,data)),2,mean)/(2*eps)
      }
    }
    
    
    GT[[TYPE]] <- GRAD.True
  }
  
  GT.Full <- matrix(0, 
                    dim(GT[[1]])[1]+dim(GT[[2]])[1]+dim(GT[[3]])[1],
                    dim(GT[[1]])[2]+dim(GT[[2]])[2]+dim(GT[[3]])[2])
  
  GT.Full[1:dim(GT[[1]])[1], 1:dim(GT[[1]])[2]] <- GT[[1]]
  GT.Full[dim(GT[[1]])[1]+1:dim(GT[[2]])[1], 
          dim(GT[[1]])[2]+1:dim(GT[[2]])[2]] <- GT[[2]]
  GT.Full[dim(GT[[1]])[1]+dim(GT[[2]])[1]+1:dim(GT[[3]])[1], 
          dim(GT[[1]])[2]+dim(GT[[2]])[2]+1:dim(GT[[3]])[2]] <- GT[[3]]
  
  
  return( GT.Full )
  
  # dg_1(O,t)/dt_1 dg_2(O,t)/dt_1 dg_3(O,t)/dt_1 ...
  # dg_1(O,t)/dt_2 dg_2(O,t)/dt_2 dg_3(O,t)/dt_2 ...
  # dg_1(O,t)/dt_3 dg_2(O,t)/dt_3 dg_3(O,t)/dt_3 ...
  # (2+num.inner)*(2+dim.outer)
}

