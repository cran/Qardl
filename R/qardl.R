#'Qardl function
#'
#'@param formula y~z1+z2
#'@param data the dataframe
#'@param maxlag maximum lag number
#'@param tau the quantile(s) to be estimated, this is generally a number strictly between 0 and 1
#'
#'@importFrom stats lm pchisq BIC as.formula model.frame model.matrix model.response na.omit sd update vcov residuals coef nobs pf pnorm df.residual formula na.exclude
#'@importFrom dplyr sapply
#'@importFrom pbapply pblapply
#'@importFrom quantreg rq
#'@importFrom MASS ginv
#'
#'
#'@examples
#'
#' # Quantile ARDL regression
#' # load data
#' data(exampledata)
#' # Fit the model
#' reg=qardl(y~z1+z2,exampledata,maxlag=7, tau=0.5)
#' reg
#'
#'@return the short-run and the long-run estimated coefficients of the QARDL model
#'
#'
#'@export

qardl<-function(formula,data,maxlag=4,tau=NULL){
  #***************************************
  #* Dr.Taha Zaghdoudi
  #* Research Assistant, University of Jendouba, Tunisa
  #* *************************************
  #* This is an R code for quantile ARDL of
  #* Cho, J. S., Kim, T. H., & Shin, Y. (2015).
  #* Quantile cointegration in the autoregressive distributed-lag modeling framework.
  #* Journal of Econometrics, 188(1), 281-300.
  #***********************************
  #* Matlab verison of Sangwoo Park (2020)
  #* ************************************
  #formula= ltotren~lgeop+lcpu
  #data= dusa
  #maxlag=7
  #tau=0.25
  aa<-unlist(as.character(formula))
  lhs   <- aa[[2]]
  core  <- aa[[3]]
  fm <- paste(lhs,"~",core)
  ffm<-as.formula(fm)
  fffm<-update(ffm, ~ . -1)
  mf <- model.frame(formula=fffm, data=data)
  xx <- model.matrix(attr(mf, "terms"), data=mf)
  k<-ncol(mf)
  y <- model.response(mf)
  y<-as.matrix(y)
  colnames(y)<-lhs
  dy<-diff(y)
  colnames(dy)<-paste("D",lhs,sep=".")
  dx=diff(xx)
  colnames(dx)<-paste("D",colnames(dx),sep=".")
  df=nrow(y)-ncol(xx)

  #**************************************
  #* grid search
  #* ************************************
  grid0=expand.grid(rep(list(1:maxlag),k-1),KEEP.OUT.ATTRS = FALSE)
  obs=nrow(y)
  dx=rbind(c(NA),dx)
  pd=matrix(c(0),nrow(grid0),1)
  for(i in 1:nrow(grid0)){
    pd[i]= grid0[i,2]-1
  }
  grid=cbind(grid0,pd)

  ss=pblapply(1:nrow(grid0) ,function(i) lm(y~nlagm1(y,grid[i,1])+lags(xx,grid[i,2])+mlags(dx,grid[i,3])))
  ak=sapply(1:nrow(grid0),function(x) BICC(ss[[x]]))

  cii=as.matrix(ak)
  np=which.min(cii)
  pg=grid[np,]

  lay1=nlagm1(y,pg[,1])
  lagx=as.matrix(lags(xx,pg[,2]))
  colnames(lagx)=paste("L",colnames(xx),sep = ".")
  lagdx=mlags(dx,pg[,3])
  rhnams<-c("(Intercept)",colnames(lay1),colnames(lagx), colnames(lagdx))
  fits= rq(y~lay1+lagx+lagdx,tau=tau )

  #fits=rq(y~nlagm1(y,2)+lags(xx,1)+mlags(dx,0), tau=tau)

  #*******************************************
  #* tout ce qui est au dessus est correct rest le long run
  names(fits$coefficients)<-rhnams
  sels<-summary(fits,se ="iid" )#"boot", bsmethod= "xy")
  #for(i in 1:length(tau)){
 # rownames(sels[[i]]$coefficients)<-rhnams
 # }
  rownames(sels$coefficients)<-rhnams


  errors=sels$residuals
  p=pg[[1]]
  q=pg[[2]]

  #**************************
  #long run estimation
  #**************************
  coeff<-sels$coefficients
  nlvars<-length(coeff[,1])
  lvars<-coeff[(p+2):nlvars,1]
  coof<-lvars/(1-sum(coeff[2:(p+1)]))
  cof<- matrix(coof, length(lvars),1)
  #**************************
  #SE by delta method
  #**************************
  #*
  braw=dx[-1,]
  tw=cbind(matrix(c(1),nrow(braw),1), braw)
  mm=(crossprod(xx[(q+1):obs,])-crossprod(xx[(q+1):obs,],tw[q:(obs-1),])%*%ginv(crossprod(tw[q:(obs-1),]))%*%crossprod(tw[q:(obs-1),],xx[(q+1):obs,]))/(obs-q)^2
  hb = (4.5*dnorm(qnorm(tau))^4/(obs*(2*qnorm(tau)^2+1)^2))^0.2
  fh = mean(dnorm(-errors/hb))/hb
  bb = 1/((1-sum(coeff[2:(p+1)]))*fh)
  psu=matrix(c(0),2,1)
  psu[1,1]=tau
  psu[2,1]=tau
  qq=(min(psu,1)-tau^2)*bb^2
  lrvcov=kronecker(qq, ginv(mm))
  lrse=sqrt(diag(lrvcov/(obs-1)^2))
  lrt<-coof/lrse
  lrpv<-2 * pnorm(-abs(lrt))
  lres<-cbind(cof,lrse,lrt,lrpv)
  #lres
  rownames(lres)<-names(lvars)
  colnames(lres)<-c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )
  #***********************************
  #* testing short run parameters
  #* *********************************
  if(p>q){
    xxi=xx[(p+1):obs,]
    wi=dx[(p+1):obs,]
    yi=as.matrix(lay1[(p+1):obs,])
    fits2=NULL
    kk=matrix(c(0), (obs-p),p)
    for(i in 1: ncol(yi)){
      fits2= rq(yi[,i]~xxi+wi,tau=tau )
      uu=fits2$residuals
      kk[,i]=uu

    }
    tilw=tw[p:(obs-1),]
    #lll = (crossprod(kk)-crossprod(kk,tilw)%*%ginv(crossprod(tilw))%*%crossprod(tilw,kk))/(obs-p)
  }else{
    xxi=xx[(q+1):obs,]
    wi=dx[(q+1):obs,]
    yi=as.matrix(lay1[(q+1):obs,])
    fits2=NULL
    kk=matrix(c(0), (obs-q),p)
    for(i in 1: ncol(yi)){
      fits2= rq(yi[,i]~xxi+wi,tau=tau )
      uu=fits2$residuals
      kk[,i]=uu
    }

    tilw=tw[(q:obs-1),]
    #lll = (crossprod(kk) -crossprod(kk,tilw)%*%ginv(crossprod(tilw))%*%crossprod(tilw,kk))/(obs-q)
  }


  out=list(sels=sels,lres=lres,fits=fits, lrvcov=lrvcov, mm=mm, qq=qq, hb=hb, fh=fh,
           bb=bb,p=p,q=q, kk=kk, tw=tw, tilw=tilw, obs=obs)
  class(out) <- "qardl"

  return(out)


}

