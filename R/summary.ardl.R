#' Summary of a ardl model
#'
#' \code{summary} method for a \code{\link{qardl}} model.
#'
#' @param object is the object of the function
#' @param ... not used
#'
#' @return an object of the S3 class \code{summary.qardl} with the
#' following components:
#'
#' @importFrom stats printCoefmat
#' @rdname summary.qardl
#' @name summary.qardl
#' @export
summary.qardl<-function(object,...)
{
  cat("==============================================================\n")
  print(object$sels)
  cat("==============================================================\n")
  cat("\nLong-run coefficients\n")
  printCoefmat(object$lres,has.Pvalue = TRUE,signif.stars = TRUE)
  cat("==============================================================\n")
}
