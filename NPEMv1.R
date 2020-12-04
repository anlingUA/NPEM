NPEM <- function(X,M,Y,method="UV"){
  ### Uses Mutual Information to test for significant alpha
  ### (X->M), beta (M->Y), gamma (X->Y) and mediated (X->M->Y)
  ### effects in a 3 variable mediation model
  ###
  ### Args:
  ###   X              - exposure variable of interest
  ###   M              - the mediator of interest
  ###   Y              - the outcome of interest
  ###   Method         - determines which algorithm is used. Currently accepts
  ###                    UV, UVS, BV, BVS. See README for details.
  ### Return:
  ###   NPEM.Results   - List of p-values for different steps of mediation
  ###
  ### Details:
  ###   The returned list contains four features:
  ###   1) gamma.c     - contributed information score(s) from each exposure to response
  ###   2) gamma.p     - p-values for each exposure's total effect on response
  ###   3) alpha.c     - contributed information score(s) from each exposure to mediator
  ###   4) alpha.p     - p-values for each exposure's effect on each mediator
  ###   5) beta.c      - contributed information score(s) from each mediator to response
  ###   6) beta.p      - p-values for each mediator's effect on response
  ###   7) mediation.p - p-values for that the mediator of interest mediates an exposure's effect
  
  
  ###################### Load Required Packages ######################
  library(np)
  library(outliers)
  
  ######################### NPEM ##########################
  
  ### Save variable names ###
  x.names <- colnames(X)
  if(is.null(x.names)) x.names=sapply(1:ncol(X),function(s) paste0("X",s))
  m.names <- colnames(M)
  if(is.null(m.names)) m.names=sapply(1:ncol(M),function(s) paste0("M",s))
  
  ####### UV #######
  if(method=="UV"){
    ### Gamma ###
    A.gamma=apply(X,2,np.MI,Y)
    R.gamma=GrubbsOrderDF(A.gamma,X,Y,pcut=0.05)
    G.gamma=R.gamma[[1]]
    C.gamma=R.gamma[[2]]
    
    ### Alpha ###
    A.alpha=C.alpha=G.alpha=matrix(0L,ncol=ncol(X),nrow=ncol(M))
    for(j in 1:ncol(M)){
      A.alpha[j,]=apply(X,2,function(s) np.MI(s,M[,j]))
      
      R.alpha=GrubbsOrderDF(A.alpha[j,],X,M[,j],pcut=0.05)
      G.alpha[j,]=R.alpha[[1]]
      C.alpha[j,]=R.alpha[[2]]
    }
    
    ### Beta ###
    A.beta=apply(M,2,np.MI,Y)
    R.beta=GrubbsOrderDF(A.beta,M,Y,pcut=0.05)
    G.beta=R.beta[[1]]
    C.beta=R.beta[[2]]
    
    #######Mediation#######
    G.Med=apply(cbind(G.beta,apply(G.alpha,1,min)),1,max)
    
    #Return Outcomes#
    names(C.gamma)=names(G.gamma)=colnames(C.alpha)=colnames(G.alpha)=x.names
    names(C.beta)=names(G.beta)=rownames(C.alpha)=rownames(G.alpha)=names(G.Med)=m.names
    NPEM.Results=list("gamma.c"=C.gamma,
                      "gamma.p"=G.gamma,
                      "alpha.c"=C.alpha,
                      "alpha.p"=G.alpha,
                      "beta.c"=C.beta,
                      "beta.p"=G.beta,
                      "mediation.p"=G.Med)
  }
  
  if(method=="UVS"){
    ### Gamma ###
    A.gamma=apply(X,2,np.MI,Y)
    R.gamma=GrubbsOrderDF(A.gamma,X,Y,pcut=0.05)
    G.gamma=GrubbsOrderDFS(R.gamma[[2]],pcut=0.05)
    C.gamma=R.gamma[[2]]
    
    ### Alpha ###
    A.alpha=C.alpha=G.alpha=matrix(0L,ncol=ncol(X),nrow=ncol(M))
    for(j in 1:ncol(M)){
      A.alpha[j,]=apply(X,2,function(s) np.MI(s,M[,j]))
      
      R.alpha=GrubbsOrderDF(A.alpha[j,],X,M[,j],pcut=0.05)
      G.alpha[j,]=GrubbsOrderDFS(R.alpha[[2]],pcut=0.05)
      C.alpha[j,]=R.alpha[[2]]
    }
    
    ### Beta ###
    A.beta=apply(M,2,np.MI,Y)
    R.beta=GrubbsOrderDF(A.beta,M,Y,pcut=0.05)
    G.beta=GrubbsOrderDFS(R.beta[[2]],pcut=0.05)
    C.beta=R.beta[[2]]
    
    #######Mediation#######
    G.Med=apply(cbind(G.beta,apply(G.alpha,1,min)),1,max)
    
    #Return Outcomes#
    names(C.gamma)=names(G.gamma)=colnames(C.alpha)=colnames(G.alpha)=x.names
    names(C.beta)=names(G.beta)=rownames(C.alpha)=rownames(G.alpha)=names(G.Med)=m.names
    NPEM.Results=list("gamma.c"=C.gamma,
                      "gamma.p"=G.gamma,
                      "alpha.c"=C.alpha,
                      "alpha.p"=G.alpha,
                      "beta.c"=C.beta,
                      "beta.p"=G.beta,
                      "mediation.p"=G.Med)
  }
  
  if(method=="BV"){
    ### Gamma ###
    A.gamma=apply(X,2,np.MI,Y)
    R.gamma=GrubbsOrderDF(A.gamma,X,Y,pcut=0.05)
    G.gamma=R.gamma[[1]]
    C.gamma=R.gamma[[2]]
    
    ### Alpha ###
    MZ=as.factor(1*(M>0))
    dim(MZ)=c(dim(M)[1],dim(M)[2])
    nz=apply(M,2,function(s) which(s>0))
    A.alpha1=A.alpha2=C1.alpha=C2.alpha=G.alpha=matrix(0L,ncol=ncol(X),nrow=ncol(M))
    for(j in 1:ncol(M)){
      nzi=nz[[j]]
      A.alpha1[j,]=apply(X,2,function(s) np.MI(s[nzi],M[nzi,j]))
      A.alpha2[j,]=apply(X,2,function(s) np.MI(s,MZ[,j]))
      
      R.alpha=GrubbsOrderDF2(A.alpha1[j,],A.alpha2[j,],X[nzi,],X,M[nzi,j],MZ[,j],pcut=0.05)
      G.alpha[j,]=R.alpha[[1]]
      C1.alpha[j,]=R.alpha[[2]][,1]
      C2.alpha[j,]=R.alpha[[2]][,2]
    }
    
    ### Beta ###
    A.beta1=apply(as.matrix(1:ncol(M)),1,function(s) np.MI(M[nz[[s]],s],Y[nz[[s]],]))
    A.beta2=apply(MZ,2,np.MI,Y)
    R.beta=GrubbsOrderDF2B(A.beta1,A.beta2,M,MZ,Y,Y,nz,pcut=0.05)
    G.beta=R.beta[[1]]
    C1.beta=R.beta[[2]][,1]
    C2.beta=R.beta[[2]][,2]
    
    #######Mediation#######
    G.Med=apply(cbind(G.beta,apply(G.alpha,1,min)),1,max)
    
    #Return Outcomes
    names(C.gamma)=names(G.gamma)=colnames(C1.alpha)=colnames(C2.alpha)=colnames(G.alpha)=x.names
    names(C1.beta)=names(C2.beta)=names(G.beta)=rownames(C1.alpha)=rownames(C2.alpha)=rownames(G.alpha)=names(G.Med)=m.names
    NPEM.Results=list("gamma.c"=C.gamma,
                      "gamma.p"=G.gamma,
                      "alpha.c"=list(C1.alpha,C2.alpha),
                      "alpha.p"=G.alpha,
                      "beta.c"=list(C1.beta,C2.beta),
                      "beta.p"=G.beta,
                      "mediation.p"=G.Med)
  }
  
  if(method=="BVS"){
    ### Gamma ###
    A.gamma=apply(X,2,np.MI,Y)
    R.gamma=GrubbsOrderDF(A.gamma,X,Y,pcut=0.05)
    G.gamma=GrubbsOrderDFS(R.gamma[[2]],pcut=0.05)
    C.gamma=R.gamma[[2]]
    
    ### Alpha ###
    MZ=as.factor(1*(M>0))
    dim(MZ)=c(dim(M)[1],dim(M)[2])
    nz=apply(M,2,function(s) which(s>0))
    A.alpha1=A.alpha2=C1.alpha=C2.alpha=G.alpha=matrix(0L,ncol=ncol(X),nrow=ncol(M))
    for(j in 1:ncol(M)){
      nzi=nz[[j]]
      A.alpha1[j,]=apply(X,2,function(s) np.MI(s[nzi],M[nzi,j]))
      A.alpha2[j,]=apply(X,2,function(s) np.MI(s,MZ[,j]))
      
      R.alpha=GrubbsOrderDF2(A.alpha1[j,],A.alpha2[j,],X[nzi,],X,M[nzi,j],MZ[,j],pcut=0.05)
      G.alpha[j,]=GrubbsOrderDF2S(R.alpha[[2]],pcut=0.05)
      C1.alpha[j,]=R.alpha[[2]][,1]
      C2.alpha[j,]=R.alpha[[2]][,2]
    }
    
    ### Beta ###
    A.beta1=apply(as.matrix(1:ncol(M)),1,function(s) np.MI(M[nz[[s]],s],Y[nz[[s]],]))
    A.beta2=apply(MZ,2,np.MI,Y)
    R.beta=GrubbsOrderDF2B(A.beta1,A.beta2,M,MZ,Y,Y,nz,pcut=0.05)
    G.beta=GrubbsOrderDF2S(R.beta[[2]],pcut=0.05)
    C1.beta=R.beta[[2]][,1]
    C2.beta=R.beta[[2]][,2]
    
    #######Mediation#######
    G.Med=apply(cbind(G.beta,apply(G.alpha,1,min)),1,max)
    
    #Return Outcomes#
    names(C.gamma)=names(G.gamma)=colnames(C1.alpha)=colnames(C2.alpha)=colnames(G.alpha)=x.names
    names(C1.beta)=names(C2.beta)=names(G.beta)=rownames(C1.alpha)=rownames(C2.alpha)=rownames(G.alpha)=names(G.Med)=m.names
    NPEM.Results=list("gamma.c"=C.gamma,
                      "gamma.p"=G.gamma,
                      "alpha.c"=list(C1.alpha,C2.alpha),
                      "alpha.p"=G.alpha,
                      "beta.c"=list(C1.beta,C2.beta),
                      "beta.p"=G.beta,
                      "mediation.p"=G.Med)
  }
  return(NPEM.Results)
}

np.MI <- function(data.x,data.y){
  require('np')
  bw.x <- npudensbw(data.x, bwmethod='normal-reference')
  bw.y <- npudensbw(data.y, bwmethod='normal-reference')
  data=data.frame(data.x,data.y)
  bw.xy<-npudensbw(data,bws=c(bw.x$bw,bw.y$bw),bandwidth.compute=F,bwmethod='normal-reference')
  kd.x <- npudens(bw.x)
  kd.y <- npudens(bw.y)
  kd.xy <- npudens(bw.xy)
  
  kd.x.fit<-fitted(kd.x)
  kd.y.fit<-fitted(kd.y)
  kd.xy.fit<-fitted(kd.xy)
  summand <- log(kd.xy.fit/(kd.x.fit*kd.y.fit))
  return(mean(summand))
}

GrubbsOrderDF <- function(A,X,Y,pcut=0.5){
  require(outliers)
  require(np)
  names(A)=1:length(A)
  S=NULL
  PS=array(1L,length(A))
  PenStore=array(0,length(A))
  for(index_i in 1:length(A)){
    if(is.null(S)){
      C=A
    } else if(length(A)-length(S)<3){
      S=c(S,as.numeric(names(A)[-S]))
      break
    } else{
      FS=X[,-S]
      PenNew=apply(as.matrix(1:ncol(FS)),1,function(s) np.MI(cbind(FS[,s],X[,fni]),Y)-A[fni])
      PenStore[-S]=(PenStore[-S]*((index_i-2)**2)+PenNew)/((index_i-1)**2)
      C=A[-S]-PenStore[-S]
      names(C)=names(A[-S])
    }
    fni=as.numeric(names(C)[which.max(C)])
    S=c(S,fni)
    if(is.na(sd(C))||sd(C)==0) G=0 else G=(max(C)-mean(C))/sd(C)
    P.val=1-pgrubbs(G,length(C))
    PS[index_i]=P.val
    if(PS[index_i]>pcut) {
      S=c(S,as.numeric(names(A)[-S]))
      break
    }
  }
  outcome=PS[order(S)]
  CIs=A-PenStore
  return(list(outcome,CIs))
}

GrubbsOrderDF2<-function(A1,A2,X1,X2,Y1,Y2,pcut=0.5){
  require(MASS) 
  names(A1)=names(A2)=1:length(A1)
  S=NULL
  PS=array(1L,length(A1))
  FS=data.frame(A1,A2)
  L=FS1=A1
  FS2=A2
  Pen1Store=Pen2Store=array(0,length(A1))
  for(index_i in 1:length(A1)){
    if(is.null(S)){
      C1=A1
      C2=A2
    } else if(length(A1)-length(S)<3){
      S=c(S,as.numeric(names(A1)[-S]))
      break
    } else{
      FS1=X1[,-S]
      FS2=X2[,-S]
      Pen1New=apply(as.matrix(1:ncol(FS1)),1,function(s) np.MI(cbind(FS1[,s],X1[,fni]),Y1)-A1[fni])
      Pen2New=apply(as.matrix(1:ncol(FS2)),1,function(s) np.MI(cbind(FS2[,s],X2[,fni]),Y2)-A2[fni])
      Pen1Store[-S]=(Pen1Store[-S]*((index_i-2)**2)+Pen1New)/((index_i-1)**2)
      Pen2Store[-S]=(Pen2Store[-S]*((index_i-2)**2)+Pen2New)/((index_i-1)**2)
      C1=A1[-S]-Pen1Store[-S]
      C2=A2[-S]-Pen2Store[-S]
      names(C1)=names(A1[-S])
      names(C2)=names(A1[-S])
    }
    C=data.frame(C1,C2)
    # Choose highest contributed information
    covC=NA
    covC=tryCatch(solve(cov(C)),error=function(e){NA})
    oldw <- getOption("warn")
    options(warn = -1)
    if(is.na(covC)){
      S=c(S,as.numeric(names(A1)[-S]))
      break
    }
    L[1:length(L)]=mahalanobis(C,colMeans(C),cov(C))
    fni=as.numeric(names(L)[which.max(L)])
    S=c(S,fni)
    # Caclulate P-value
    P.val=1-pchisq(max(L),df=ncol(C))
    
    # Save P-values
    PS[index_i]=P.val
    PA=p.adjust(PS,"BH")
    
    # Check if stopping conditions are met
    if(PA[index_i]>pcut) {
      S=c(S,as.numeric(names(A1)[-S]))
      break
    }
    options(warn = oldw)
    L=A1[-S]
  }
  outcome=PA[order(S)]
  CIs<-data.frame(A1-Pen1Store,A2-Pen2Store)
  return(list(outcome,CIs))
}

GrubbsOrderDF2B<-function(A1,A2,X1,X2,Y1,Y2,nz,pcut=0.5){
  require(MASS) 
  names(A1)=names(A2)=1:length(A1)
  S=NULL
  PS=array(1L,length(A1))
  FS=data.frame(A1,A2)
  L=FS1=A1
  FS2=A2
  Pen1Store=Pen2Store=array(0,length(A1))
  for(index_i in 1:length(A1)){
    if(is.null(S)){
      C1=A1
      C2=A2
    } else if(length(A1)-length(S)<3){
      S=c(S,as.numeric(names(A1)[-S]))
      break
    } else{
      FS1=X1[,-S]
      FS2=X2[,-S]
      Pen1New=apply(as.matrix(1:ncol(FS1)),1,function(s) np.MI(cbind(FS1[nz[[s]],s],X1[nz[[s]],fni]),Y1[nz[[s]],])-A1[fni])
      Pen2New=apply(as.matrix(1:ncol(FS2)),1,function(s) np.MI(cbind(FS2[,s],X2[,fni]),Y2)-A2[fni])
      Pen1Store[-S]=(Pen1Store[-S]*((index_i-2)**2)+Pen1New)/((index_i-1)**2)
      Pen2Store[-S]=(Pen2Store[-S]*((index_i-2)**2)+Pen2New)/((index_i-1)**2)
      C1=A1[-S]-Pen1Store[-S]
      C2=A2[-S]-Pen2Store[-S]
      names(C1)=names(A1[-S])
      names(C2)=names(A1[-S])
    }
    C=data.frame(C1,C2)
    # Choose highest contributed information
    covC=NA
    covC=tryCatch(solve(cov(C)),error=function(e){NA})
    oldw <- getOption("warn")
    options(warn = -1)
    if(is.na(covC)){
      S=c(S,as.numeric(names(A1)[-S]))
      break
    }
    L[1:length(L)]=mahalanobis(C,colMeans(C),cov(C))
    fni=as.numeric(names(L)[which.max(L)])
    S=c(S,fni)
    # Caclulate P-value
    P.val=1-pchisq(max(L),df=ncol(C))
    
    # Save P-values
    PS[index_i]=P.val
    PA=p.adjust(PS,"BH")
    
    # Check if stopping conditions are met
    if(PA[index_i]>pcut) {
      S=c(S,as.numeric(names(A1)[-S]))
      break
    }
    options(warn = oldw)
    L=A1[-S]
  }
  outcome=PA[order(S)]
  CIs<-data.frame(A1-Pen1Store,A2-Pen2Store)
  return(list(outcome,CIs))
  return(outcome)
}

GrubbsOrderDFS<-function(D,pcut=0.5){
  names(D)=1:length(D)
  S=NULL
  C=D
  PS=array(1L,length(D))
  for(index_i in 1:length(D)){
    if(!is.null(S)){
      C=D[-S]
    }
    fni=as.numeric(names(C)[which.max(C)])
    S=c(S,fni)
    if(is.na(sd(C))||sd(C)==0) G=0 else G=(max(C)-mean(C))/sd(C)
    P.val=1-pgrubbs(G,length(C))
    PS[index_i]=P.val
    if(PS[index_i]>pcut) {
      S=c(S,as.numeric(names(D)[-S]))
      break
    }
  }
  outcome=PS[order(S)]
  return(outcome)
}

GrubbsOrderDF2S<-function(D,pcut=0.5){
  require(MASS) 
  D1=D[,1]
  D2=D[,2]
  names(D1)=names(D2)=1:length(D1)
  S=NULL
  PS=array(1L,length(D1))
  FS=data.frame(D1,D2)
  L=FS1=D1
  FS2=D2
  C=data.frame(D1,D2)
  for(index_i in 1:length(D1)){
    if(!is.null(S)){
      C=data.frame(D1[-S],D2[-S])
    }
    # Choose highest contributed information
    covC=NA
    covC=tryCatch(solve(cov(C)),error=function(e){NA})
    oldw <- getOption("warn")
    options(warn = -1)
    if(is.na(covC)){
      S=c(S,as.numeric(names(D1)[-S]))
      break
    }
    L[1:length(L)]=mahalanobis(C,colMeans(C),cov(C))
    fni=as.numeric(names(L)[which.max(L)])
    S=c(S,fni)
    # Caclulate P-value
    P.val=1-pchisq(max(L),df=ncol(C))
    
    # Save P-values
    PS[index_i]=P.val
    PA=p.adjust(PS,"BH")
    
    # Check if stopping conditions are met
    if(PA[index_i]>pcut) {
      S=c(S,as.numeric(names(D1)[-S]))
      break
    }
    options(warn = oldw)
    L=D1[-S]
  }
  outcome=PA[order(S)]
  return(outcome)
}