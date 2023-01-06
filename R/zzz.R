



#' Function lagm: lagged matrix
#'@param x is a matrix
#'@param nLags the number of lags
#'
#'#load data
#'data(exampledata)
#'m=as.matrix(exampledata$y)
#'lagm(m,2)
#'
#'@keywords internal
#'@export

lagm <- function(m, nLags) {
  # JM: This code is redundant. More than 2 arguments will cause an error
  # anyway, less than 2 arguments will cause an error about missing
  # arguments that's more informative than yours.

  # nargin <- length(as.list(match.call())) - 1
  # if (nargin != 2) {
  #   stop('Check function inputs')
  # }

  if(!is.matrix(m))
    stop("Trying to lag something that's not a matrix")

  d <- dim(m)

  #Add column names if they don't exist yet
  if(is.null(colnames(m)))
    colnames(m) <- as.character(seq_len(d[2]))

  #Check
  if(nLags > d[1])
    stop(sprintf("You try to create %d lags but there's only %d rows in m.",
                 nLags, d[1]))

  lagM <- matrix(NA,nrow=d[1], ncol = d[2]*nLags)

  for(i in seq_len(nLags)){
    #Make ids for the columns in result
    cid <- seq(1:d[2]) + d[2]*(i-1)

    lagM[(i+1):d[1],cid] <- m[1:(d[1]-i),]
  }

  #mnames=paste("L",colnames(m),sep = ".")

  cnames <- outer(colnames(m),seq_len(nLags), FUN = paste, sep = "_")

  colnames(lagM) <- c(cnames)

  return(lagM)

}


#' Function nlagm: lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'
#'#load data
#'data(exampledata)
#'m=as.matrix(exampledata$y)
#'nlagm(m,2)
#'
#'@keywords internal
#'@export
nlagm=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    m=lagm(x,nlags)
    mlag=cbind(x,m)
  }

  mlag
}
#' Function nlagm1: lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'
#'#load data
#'data(exampledata)
#'m=as.matrix(exampledata$y)
#'nlagm1(m,2)
#'
#'@keywords internal
#'@export
nlagm1=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    mlag=as.matrix(lagm(x,nlags))
  }

  mlag
}


#' Function lags: lagged matrix
#'@param x is a matrix
#'@param k the number of lags
#'
#'#load data
#'data(exampledata)
#'m=as.matrix(exampledata$y)
#'lags(m,k=2)
#'
#'@keywords internal
#'@export
lags <- function(x, k=1) {
  i<-is.vector(x)
  if(is.vector(x)) x<-matrix(x) else x<-matrix(x,nrow(x))
  if(k>0) {
    x <- rbind(matrix(rep(NA, k*ncol(x)),ncol=ncol(x)), matrix(x[1:(nrow(x)-k),], ncol=ncol(x)))
  }
  else {
    x <- rbind(matrix(x[(-k+1):(nrow(x)),], ncol=ncol(x)),matrix(rep(NA, -k*ncol(x)),ncol=ncol(x)))
  }
  if(i) x[1:length(x)] else x
}

#' Function mlags: lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'
#'#load data
#'data(exampledata)
#'m=as.matrix(exampledata$y)
#'mlags(m,2)
#'
#'@keywords internal
#'@export
mlags=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    m=lags(x,nlags)
    mlag=cbind(x,m)
  }

  mlag
}

#' Function BIC
#'@param f lm object
#'
#'
#'@keywords internal
#'@export
BICC <- function(f) {
  sample.size<-f$df + length(f$coeff)
  bic<-log(sum(residuals(f)^2)/sample.size)+(length(f$coeff)/sample.size)*log(sample.size)
  bic

}


#' Function matsplitter
#'@param M matrix
#'@param r matrix rows
#'@param c matrix columns
#'
#'aaa=matrix(c(seq(50)),25,2)
#'matsplitter(aaa,5,2)
#'
#'@keywords internal
#'@export
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}

#' Function wtelslr beta
#'@param cof the estimated long-run parameter by qardl
#'@param covar estimated covariance matrix
#'@param rca R matrix in the null hypothesis
#'@param rsm r matrix in the null hypothesis
#'@param data the same data set used for qardl
#'
#'#'data(exampledata)
#'hyp=hyptest(y~z1+z2,exampledata,maxlag=7,tau=c(0.25,0.5,0.75))
#'wlr=wtestlr(hyp$bigbt,hyp$lrvcov3, hyp$ca1,hyp$sm1,exampledata)
#'wlr
#'@keywords internal
#'@export
wtestlr=function(cof,covar,rca,rsm,data){
  nn = nrow(data)
  wtlrb1 = (nn-1)^2*t(rca%*%cof-rsm)%*%ginv(rca%*%covar%*%t(rca))%*%(rca%*%cof-rsm)
  rnkl = ncol(rca)
  pval=pchisq(wtlrb1,rnkl)
  out= list(wtlrb1=wtlrb1,pval=pval)
  return(out)
}


#' Function wtelslsr gamma and phi
#'@param cof the estimated short-run parameter by qardl
#'@param covar estimated covariance matrix
#'@param rca R matrix in the null hypothesis
#'@param rsm r matrix in the null hypothesis
#'@param data the same data set used for qardl
#'
#'data(exampledata)
#'hyp=hyptest(y~z1+z2,exampledata,maxlag=7,tau=c(0.25,0.5,0.75))
#'wsrg=wtestlsr(hyp$bigdam, hyp$bigff, hyp$ca1, hyp$sm1, exampledata)
#'wsrg
#'@keywords internal
#'@export
wtestlsr=function(cof,covar,rca,rsm,data){
  nn = nrow(data)
  wtlrb1 = (nn-1)*t(rca%*%cof-rsm)%*%ginv(rca%*%covar%*%t(rca))%*%(rca%*%cof-rsm)
  rnkl = nrow(rca)
  pval=pchisq(wtlrb1,rnkl)
  out= list(wtlrb1=wtlrb1,pval=pval)
  return(out)
}
