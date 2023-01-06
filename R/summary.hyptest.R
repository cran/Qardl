#' Summary of a hyptest
#'
#' \code{summary} method for a \code{\link{hyptest}}.
#'
#' @param object is the object of the function
#' @param ... not used
#'
#' @return an object of the S3 class \code{summary.hyptest} with the
#' following components:
#'
#' @rdname summary.hyptest
#' @name summary.hyptest
#' @export
summary.hyptest<-function(object,...)
{
  cat("==============================================================\n")
  cat(" Wald test (Beta):",object$wlr$wtlrb1, "Prob:",object$wlr$pval, "\n" )
  cat(" Wald test (phi):",object$wsrphi$wtlrb1, "Prob:",object$wsrphi$pval, "\n" )
  cat(" Wald test (gamma):",object$wsrg$wtlrb1, "Prob:",object$wsrg$pval, "\n" )
  cat("==============================================================\n")
}
