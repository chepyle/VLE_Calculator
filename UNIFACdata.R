# collection of UNIFAC structural groups and parameters
# found here: http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
# used data from these sources:
# H. K. Hansen, P. Rasmussen, A. Fredenslund, M. Schiller, and J Gmehling (1991) Ind. Eng. Chem. Res. 30, 2352-2355.
# R. Wittig, J. Lohmann, and J. Gmehling (2003) Ind. Eng. Chem. Res. 42, 183-188.
# K. Balslev and J. Abildskov (2002) Ind. Eng. Chem. Res. 41, 2047-205.
# C. Peng, M. N. Chan, and C. K. Chan (2001) Environ. Sci. Technol. 35, 4495-4501


#checked on example 8-14 in Properties of Gases and Liquids, Fifth Edition (Bruce E. Poling, John M. Prausnitz, John P. O’Connell)


library(openxlsx)
library(minpack.lm)
library(ggplot2)
library(signal)

gr<-read.xlsx('UNIFAC_params.xlsx',sheet=1) # Q and R
ex<-read.xlsx('UNIFAC_params.xlsx',sheet=2) # examples
ex[is.na(ex)]<-0 # remove nas
ajk_sparse<-read.xlsx('UNIFAC_params.xlsx',sheet=3)# bips
Pvap_coef <- read.xlsx('UNIFAC_params.xlsx',sheet=4)  # antoine equations

# setup some global variables for UNIFAC
Q<-gr[,'Qk']
R<-gr[,'Rk']
ID<-gr[,'Group']

amn<-matrix(rep(0,max(unique(gr$Group))^2),ncol=max(unique(gr$Group)))
amn_extend<-matrix(rep(0,length(ID)^2),ncol=length(ID))  # expanded bips for secondary groups
for (n in unique(gr$Group)){
  # unsparse amn with an ambitious command, not sure if rows/columns should be swapped
  amn[ajk_sparse$I[ajk_sparse$J==n],n] <- ajk_sparse$A[ajk_sparse$J==n]
  for (m in ajk_sparse$I[ajk_sparse$J==n]){
    amn_extend[ID==m,
               ID==n]<- ajk_sparse$A[ajk_sparse$J==n & ajk_sparse$I ==m]
  }
}
tau <- exp(-amn_extend/(T+273.15)) # convert  Amn to tau

unifac_gamma <- function(v,x=matrix(c(0.3,0.7),nrow=1),T=25){
  ncomp<-ncol(v)
  
  # UNIFAC equations coded from Tester and Modell (Ch 11 and 13)
  z=10 #usually.. table 11.2?
  
  # some matrix algebra to move stuff out of the loop:
  vmat<-matrix(as.numeric(v),ncol=ncomp)
  
  q <- diag(t(matrix(rep(Q,ncomp),ncol=ncomp)) %*% vmat)
  r <- diag(t(matrix(rep(R,ncomp),ncol=ncomp)) %*% vmat)
  L <- z/2*(r-q)-(r-1)# Table 11.2
  phi <- r*x/(sum(r* x))# volume fraction of compound
  X <- (vmat %*% t(x))/sum( vmat %*% t(x)) # Eq 13-26 mole fraction of groups
  X_pure <- apply(vmat,2,function(x){x/sum(x)}) # mole fraction of groups for pure component
  
  thetaR <- X*Q/(sum(X*Q)) # FG weighted area 
  thetaR_pure <- apply(X_pure*Q,2,function(x){x/sum(x)}) # FG weighted area , pure component, more matrix challenges from R... 
  thetaC <- q*x/(sum(q*x))
  
  log_gamC<-log(phi/x)+z/2*q*log(thetaC/phi)+L-phi/x*sum(x*L)# Table 11.2
  log_gamC[is.nan(log_gamC)]<-0  # correct for log(0 errors)
  log_GAM <- matrix(nrow=length(ID),ncol=length(x))  # duplicate logGAM
  log_GAMpure <-matrix(nrow=length(ID),ncol=length(x))
  
  
  for (i in seq(length(x))) {
    for (k in which(thetaR_pure[,i]>0)) {
      tempsums<-thetaR_pure[,i]*tau[k,]/as.numeric((thetaR_pure[,i]) %*% tau)  
      log_GAMpure[k,i] <- Q[k]*(1-log(sum(thetaR_pure[,i]*tau[,k]))-
                                  sum(tempsums*is.finite(tempsums),na.rm=TRUE))
    }
  }
  
  fa_ix<-which(apply(vmat,1,function(x){any(x>0)}))
  for (k in fa_ix){
    # print(paste0(k,' : ',sum(thetaR*tau[k,])))
    tempsums<-thetaR*tau[k,]/t(t(thetaR) %*% tau)  # use temporary arrays to identify and exclude Inf
    log_GAM[k,] <- Q[k]*(1-log(sum(thetaR*tau[,k]))-
                           sum(tempsums*is.finite(tempsums),na.rm=TRUE)) # Eq 13-26
  }
  tempsum<-vmat*(log_GAM-log_GAMpure) # use temporary arrays to identify and exclude Inf
  log_gamR<- colSums(tempsum*is.finite(tempsum),na.rm=TRUE) # Eq 13-25
  log_gam <- log_gamC + log_gamR # Table 11.2
}

TPvap <- function(T=NULL,P=NULL,mol=c('water','propanol-2')){
  #antoine equation for vapor pressure
  # T in celcius, P in bar, one must be NULL
  ix<-match(mol,Pvap_coef$Name)
  A<-numeric()
  B<- numeric()
  C<-numeric()
  for (i in seq(ix)){
    if (Pvap_coef$Equation.Form[ix[i]]=='log10(P[mmHg]) = A-(B/(T[C]+C))'){ 
      # mmhg/C
      A[i]=Pvap_coef$A[ix[i]]-2.875119 # convert to bar
      B[i]=Pvap_coef$B[ix[i]]
      C[i]=Pvap_coef$C[ix[i]]  
    }else{  # bar/K
      A[i]=Pvap_coef$A[ix[i]]
      B[i]=Pvap_coef$B[ix[i]]
      C[i]=Pvap_coef$C[ix[i]]+273.15 # convert to C
    }
  }
  
  if (is.null(T)){
    Tvap<-B/(A-log10(P))-C
    return(Tvap)
  }else if((is.null(P))){
    Pvap<-10^(A-(B/(T+C)))  # antoine equation bar vs C
    return(Pvap)
  }
}

partialP <- function(v,x,T){
  
  log_gam <- unifac_gamma(v,x,T)
  #print(log_gam)
  return(x*TPvap(T=T,mol=dimnames(v)[[2]])*exp(log_gam))
}


Txy <- function(P,molecules){
  # function to calculate Txy diagrams, currently configured for only binary mixtures, but this can change
  x1range<-c(10^seq(-5,-0.3,by=0.1),1-10^seq(-0.3,-5,by=-0.1))
  v <- sapply(molecules,FUN=function(x){ex[ex[,'Name']==x,-1]})  # chemical columns, functional group rows
  fitT0<-TPvap(P=P,mol=molecules)[2]
  Tmax <- max(TPvap(P=P,mol=molecules))
  
  y1range<-numeric()
  fitT<-numeric()
  fitTy<-numeric()
  log_gam<-matrix(ncol=2,nrow=length(x1range))
  for (ix in seq(x1range)){
    x<-matrix(c(x1range[ix],1-x1range[ix]),nrow=1)
    fitT[ix] <-nls.lm(fitT0,lower=fitT0-10,upper=Tmax,
                      function(t){P-sum(partialP(v,x,t))})$par
    fitT0<-fitT[ix]
    log_gam[ix,]<-unifac_gamma(v,x,fitT[ix])
  }
  
  
  # if gam1*x1=gam2*x2 D(gam1*x1-gam2*x2)/D(x1) is positive, flag a phase split
  spl_ix<-which(diff(sign((x1range*exp(log_gam[,1])-(1-x1range)*exp(log_gam[,2]))))==2)
  
  num_phases<-rep(1,length(x1range))
  num_phases[seq(x1range)>spl_ix[1] & seq(x1range)<=spl_ix[2]] <-2
  
  plot(x1range*exp(log_gam[,1]),x1range)
  plot(x1range*exp(log_gam[,1]),fitT)
  plot(x1range,x1range*exp(log_gam[,1]))
  for (ix in seq(x1range)){
    if (num_phases[ix]==1){  # only plot Txy data from single phase mixtures, otherwise it is a mess!
      x<-matrix(c(x1range[ix],1-x1range[ix]),nrow=1)
      y1range[ix]<-partialP(v,x,fitT[ix])[1]/P
    }else{
      fitT[ix]<-min(fitT[num_phases==1])
    }
  }
  
  return(data.frame(x=x1range,y=y1range,T=fitT))
}


library(reshape2)
#test case for UNIFAC calculation:
P<-1.022 # bar
molecules=c('water','methanol')

## txy plot:
txy<-melt(Txy(P,molecules),id.vars=c('T'))
p<-ggplot(txy,aes(x=value,y=T,color=variable))+geom_line()+geom_point()
print(p)

## xyplot
xy<-melt(txy,id.vars=c('X','Y'))
p<-ggplot(xy,aes(x=X,y=Y,color=value))+geom_line()+geom_point()+geom_abline(slope=1,intercept=0)
print(p)
