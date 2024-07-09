## R functions in support of UV analysis.

expit<-function(x){exp(x)/(1 + exp(x))}
logit<-function(x){log(x) - log(1-x)}


qqplot2<-function(x,tdf=1000,titl=""){
  ## t distribution assumed; reverts to normal if df>100
  qq.samp.mat<-matrix(nrow=10000,ncol=length(x))
  mu.estr<-mean(x)
  sd.estr<-sd(x)
  sdt.estr<-sqrt(tdf/(tdf-2))
  if (tdf>100){
    subtext<-"Normal Model"
  }
  if (tdf<=100){
    subtext<-paste(tdf," df T Model",sep="")
  }
  for (i in 1:10000){
    if (tdf>100){
      qq.samp.mat[i,]<-sort(rnorm(length(x),mean=mu.estr,sd=sd.estr))
    }
    if (tdf<=100){
      qq.samp.mat[i,]<-sort(rt(length(x),df=tdf))
    }
  }
  qq.samp<-apply(qq.samp.mat,2,median)
  qq.samp.ll<-apply(qq.samp.mat,2,quantile,prob=(0.025/nrow(uvDB)))  #bonf adjusted
  qq.samp.ul<-apply(qq.samp.mat,2,quantile,prob=(1-(0.025/nrow(uvDB))))  #bonf adjusted
  rm(qq.samp.mat)
  if (tdf>100){
    q.theo<-qnorm(p=seq(0.5,length(x)-0.5,by=1)/length(x),
                  mean=mu.estr,sd=sd.estr)
  }
  if (tdf<=100){
    q.theo<-qt(p=seq(0.5,length(x)-0.5,by=1)/length(x),df=tdf)
  }
  plot(q.theo,sort(x),xlab="Theoretical Quantiles",ylab="Sample Quantiles",
       main=paste("Posterior Expected Standardized Residuals\n",titl,sep=" "),
       sub=subtext,
       las=1,pch=16,cex=0.5,cex.main=1.05)
  if (tdf>100){
    abline(a=0,b=1,lwd=2,col=4)
    lines(q.theo,qq.samp.ll,col=3,lwd=2)
    lines(q.theo,qq.samp.ul,col=3,lwd=2)
  }
  if (tdf<=100){
    sd.ratio<-(sd.estr/sdt.estr)
    abline(a=0,b=sd.ratio,lwd=2,col=4)
    lines(q.theo,sd.ratio*qq.samp.ll,col=3,lwd=2)
    lines(q.theo,sd.ratio*qq.samp.ul,col=3,lwd=2)
  }
  abline(h=c(2.5,-2.5),lwd=2,col=2,lty=2)
  labs<-paste("B=",as.character(uvDB$batch),":V=",as.character(uvDB$variant),sep="")
  labs<-labs[order.er]
  return(NULL)
}


## "Rao-Blackwellized" estimates of Pr(Dv=1|Data):
fdv<-function(etav,m.s,v.s,lpp.cap,prior.aa,prior.bb){
  cap<-(lpp.cap/log10(exp(1))) ## corresponds to 10^lpp.cap
  nc<-ncol(etav)
  prdv<-matrix(NA,nrow=nrow(etav),ncol=nc)
  lBF<-matrix(NA,nrow=nrow(etav),ncol=nc)
  colnames(prdv)<-colnames(etav)
  colnames(lBF)<-colnames(etav)
  mu0<-m.s[,"alpha[1]"]
  mu1<-m.s[,"alpha[2]"]
  expit<-function(x){exp(x)/(1 + exp(x))}
  logit<-function(x){log(x) - log(1-x)}
  lPriorOdds<-log(prior.aa/prior.bb)  ##m.s[,"a.pp"]
  var0<-v.s[,"sigma2.eta[1]"]
  var1<-v.s[,"sigma2.eta[2]"]
  for (i in 1:nc){
    lBF[,i]<-(dnorm(etav[,i],mean=mu1,sd=sqrt(var1),log=TRUE) - dnorm(etav[,i],mean=mu0,sd=sqrt(var0),log=TRUE))
  }
  lPostOdds<-sweep(lBF,MAR=1,FUN="+",STATS=lPriorOdds)
  prdv<-expit(lPostOdds)
  prdv[(!is.na(lPostOdds))&(lPostOdds>cap)]<-1.0
  prdv[(!is.na(lPostOdds))&(lPostOdds<(-cap))]<-0.0
  ##lPostOdds[(!is.na(lPostOdds))&(lPostOdds>cap)]<-cap
  ##lPostOdds[(!is.na(lPostOdds))&(lPostOdds<(-cap))]<-(-cap)
  prdv.mu<-apply(prdv,2,mean)
  lPostOdds.mu<-logit(prdv.mu)
  lPriorOdds.mu<-mean(lPriorOdds)
  lBF.mu<-(lPostOdds.mu - lPriorOdds.mu)
  return(list(prdv=prdv.mu,
              lPostOdds=lPostOdds.mu,
              lBF=lBF.mu))
}


