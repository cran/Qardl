% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{lags}
\alias{lags}
\title{Function lags: lagged matrix}
\usage{
lags(x, k = 1)
}
\arguments{
\item{x}{is a matrix}

\item{k}{the number of lags

#load data
data(exampledata)
m=as.matrix(exampledata$y)
lags(m,k=2)}
}
\description{
Function lags: lagged matrix
}
\keyword{internal}
