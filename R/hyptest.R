#'hyptest function
#'
#'@param formula y~z1+z2
#'@param data the dataframe
#'@param maxlag maximum lag number
#'@param tau the quantile(s) to be estimated, this is generally a number strictly between 0 and 1
#'
#'@importFrom Matrix bdiag
#'
#'
#'@examples
#'
#' # Quantile ARDL regression
#' # load data
#' data(exampledata)
#' # Fit the model
#' hyp=hyptest(y~z1+z2,exampledata,maxlag=7, tau=c(0.2,0.5,0.75))
#' summary(hyp)
#'
#'@return the short-run phi and gamma wald test and the long-run beata wald test
#'
#'
#'@export


hyptest=function(formula,data,maxlag=7,tau=NULL ){

  if(length(tau)<=1){stop("tau length must be more than 2")}

  #formula= ltotren~lgeop+lcpu
  #data= dusa
  #maxlag=7
  #tau=c(0.25,0.5,0.75)
  aa<-unlist(as.character(formula))
  lhs   <- aa[[2]]
  core  <- aa[[3]]
  #get regressors count
  nr=length(all.vars(formula))-1
 #run multiple quardl models
 mrq=lapply(1:length(tau), function(i) qardl(formula = formula,data=data,maxlag=maxlag,tau=tau[i]))

#*******************************************
#* estimate the long-run covariance matrix
#* ****************************************

bb=lapply(1:length(tau), function(i) mrq[[i]]$bb )
bb=do.call(rbind,bb)

ss=length(tau)
qq = matrix(0,ss,ss)

for(jj in 1:ss){
  for(ii in 1:ss){
    psu = matrix(0,2,1)
    psu[1,1] = tau[jj]
    psu[2,1] = tau[ii]
    qq[jj,ii] = (min(psu,1) - tau[jj]*tau[ii])*bb[jj]*bb[ii]
  }
}

lrvcov3= kronecker(qq, ginv(mrq[[1]]$mm))
#******************************************
#* construct the hypothisis
#* ***************************************

#hyps
ca1 = matrix(0,nr,2*length(tau))
ca1[1,1] = 1
ca1[1,2+1] = -1
ca1[2,2+1] = 1
ca1[2,2*2+1]= -1
sm1 = matrix(0,2,1)

# betas from long-run
bigbt1=lapply(1:length(tau), function(i) mrq[[i]]$lres[,1] )
bigbt2= data.frame(do.call(rbind, bigbt1))
bigbt3=as.matrix(bigbt2 %>% dplyr::select(starts_with('L.')))
bigbt=matrix(t(bigbt3),nrow(bigbt3)*ncol(bigbt3),1)
#*************************************
#* wald test for the beata long-run
#* *********************************
wlr=wtestlr(bigbt,lrvcov3, ca1,sm1,data)

#********************************************
#* short-run wald test phi
#* *****************************************
p=mrq[[1]]$p
kk=lapply(1:length(tau), function(i) mrq[[i]]$kk )
kk=do.call(cbind,kk)

srlrv=(crossprod(kk) -crossprod(kk,mrq[[1]]$tilw)%*%ginv(crossprod(mrq[[1]]$tilw))%*%crossprod(mrq[[1]]$tilw,kk))/(mrq[[1]]$obs-mrq[[1]]$q)

fh=lapply(1:length(tau), function(i) mrq[[i]]$fh )
fh=do.call(rbind,fh)

cc=matrix(0,ss,ss)

for(jj in 1:ss){
  for(ii in 1:ss){
    psu = matrix(0,2,1)
    psu[1,1] = tau[jj]
    psu[2,1] = tau[ii]
    cc[jj,ii] = (min(psu,1) - tau[jj]*tau[ii])/(fh[ii]*fh[jj])
  }
}
lll=srlrv
psu3 = ginv(lll[1:p,1:p])%*%lll[1:p,1:p]%*%ginv(lll[1:p,1:p])
bigpi = cc%x%psu3
bigphi=lapply(1:length(tau), function(i) mrq[[i]]$sels$coefficients[2:(p+1),1])# 3=regressor count+1 for intercept
bigphi=do.call(rbind,bigphi)
bigphi=matrix(t(bigphi),nrow(bigphi)*ncol(bigphi),1)

# hyp wald test phi short-run

ca2 = matrix(0,2,p*length(tau))
ca2[1,1] = 1
ca2[1,p+1] = -1
ca2[2,p+1] = 1
ca2[2,2*p+1] = -1
sm2 = sm1
wsrphi=wtestlsr(bigphi, bigpi, ca2, sm2, data)

#*****************************************
#* short run gamma wald test
#*****************************************
bigdam1=lapply(1:length(tau), function(i) mrq[[i]]$sels$coefficients[,1])
bigdam2=data.frame(do.call(rbind,bigdam1))
bigdam3=as.matrix(dplyr::select(bigdam2, starts_with("L.")))
bigdam=matrix(bigdam3,nrow(bigdam3)*ncol(bigdam3),1)
bigdam4=t(bigbt3%x%matrix(1,p,1))
bilam=as.matrix(bdiag(lapply(1:length(tau), function(i) matsplitter(bigdam4,2,p)[,,i]))) # nember of regressors
bigff = bilam%*%bigpi%*%t(bilam)
wsrg=wtestlsr(bigdam, bigff, ca1, sm1, data)

out= list(wlr=wlr, wsrphi=wsrphi, wsrg=wsrg, bigdam=bigdam,
          bigff=bigff,bigphi=bigphi,bigpi=bigpi, bigbt=bigbt, lrvcov3=lrvcov3, ca1=ca1,sm1=sm1, ca2=ca2)

class(out) <- "hyptest"

return(out)
}


