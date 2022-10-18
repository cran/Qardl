#' Function lagm which return a lagged matrix
#'@param m is a matrix
#'@param nLags the number of lags
#'
#'@examples
#'
#'  data(exampledata)
#'  lagm(as.matrix(exampledata$y),2)
#'@return A lagged matrix
#'
#'@keywords internal
#'

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


#' Function nlagm return a  lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'@examples
#'
#'  data(exampledata)
#'  nlagm(as.matrix(exampledata$y),2)
#'
#'@return A lagged matrix
#'@keywords internal
nlagm=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    m=lagm(x,nlags)
    mlag=cbind(x,m)
  }

  mlag
}
#' Function nlagm1 which return a lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'@examples
#'
#' data(exampledata)
#' nlagm1(as.matrix(exampledata$y),2)
#'
#'@return A lagged matrix
#'@keywords internal
nlagm1=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    mlag=as.matrix(lagm(x,nlags))
  }

  mlag
}


#' Function lags which return a lagged matrix
#'@param x is a matrix
#'@param k the number of lags
#'@examples
#'
#'  data(exampledata)
#' lags(exampledata$y,2)
#'
#'@return A lagged matrix
#'@keywords internal
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

#' Function mlags which return a lagged matrix
#'@param x is a matrix
#'@param nlags the number of lags
#'@examples
#'
#'  data(exampledata)
#'  mlags(as.matrix(exampledata$y),2)
#'
#'@return  A lagged matrix
#'@keywords internal

mlags=function(x,nlags){
  if(nlags==0){
    mlag=as.matrix(x)
  }else{
    m=lags(x,nlags)
    mlag=cbind(x,m)
  }

  mlag
}


