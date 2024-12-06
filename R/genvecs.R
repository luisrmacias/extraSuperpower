
#' Internal function to generate a vector of coefficients to be applied to the reference mean to generate an experimental model.
#'
#' The vector's length is length equal to the number of levels of the experimental factor, defined in reps.
#' he coefficients are estimated from a user defined change. Default is to generate linear vector, this can be modified with the
#' @noRd
genvecs <- function(change, reps, bystart=FALSE, scaler=refmean)
{
  if(reps%%1>0 | reps<=0)
  {stop("Number of factor levels (nlfA and nlfB) must be positive integers")}
  if(bystart)
  {
    reps <- reps-1
    increments <- scaler*change - scaler
    increments <- increments*(1:reps)
    increments <- c(1, 1+ (increments/scaler))
  } else if (!bystart)
  {
    increments <- seq(1, change, length.out=reps)
  }
  increments
}
