
#' Internal function to generate a matrix of coefficients that will be multiplied to the reference mean to generate an experimental model.
#'
#'
#'
#' @noRd

build_mean_mat <- function(fAvec, fBvec, iA, a, b, label_list = label_list, bystart=TRUE)
{
    if(all(diff(diff(fAvec))<1e-4) | !bystart)
    {
      Bdelta <- diff(fAvec)[1]
    }
    if(iA!=1)
    {
      Bincrements <- seq(0, Bdelta*(a-1), Bdelta)
      mmat <- outer(Bincrements, fBvec, FUN = "+")
    }
    else if (iA==1 & bystart)
    {
      mmat <- rep(fBvec, a)
      mmat <- matrix(mmat, a, b, byrow = TRUE)
    } else if (iA==1 & !bystart)
    {
      mmat <- rep(fBvec, each=a)
      mmat <- matrix(mmat, a, b)
    }
  dimnames(mmat) <- label_list
  mmat
}
