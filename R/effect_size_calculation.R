#' Effect size calculation 
#'
#' Calculate effect sizes for two-way factorial designs from matrices of expected mean and standard deviation values
#' for each combination of factor levels. The output given is Cohen's f. Calculations are done as exemplified in the 
#' G*Power 3.1 manual.
#'
#' @param matrices_obj List of 2 matrices, named mean.mat and sd.mat. It is the output of the mean_sd_matrix function
#' @examples 
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' effects_treat_time <- mean_sd_matrix(refmean = refmean, f1effect = treateff, f2effect = timeeff, nlf1 = treatgroups, nlf2 = timepoints, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
#' anoveff(effects_treat_time)
#'
#' @export
anoveff <- function(matrices_obj)
  {
  mean_matrix <- matrices_obj$mean.mat
  factors <- names(dimnames(mean_matrix))
  sd_matrix <- matrices_obj$sd.mat
  sim_mean_mat <- mean_matrix - mean(mean_matrix)
  overallvariance <- mean(sd_matrix^2)
  fac1var <- mean(rowMeans(sim_mean_mat)^2)
  fac1f <- sqrt(fac1var/overallvariance)
  fac2var <- mean(colMeans(sim_mean_mat)^2)
  fac2f <- sqrt(fac2var/overallvariance)
  resids <- sim_mean_mat - outer(rowMeans(sim_mean_mat), colMeans(sim_mean_mat), "+")
  interactf <- sqrt(mean(resids^2)/overallvariance)
  ess <- c(fac1f, fac2f, interactf)
  names(ess) <- c(factors, paste(factors, collapse =":"))
  list("Cohen's f effect sizes for the input means and standard deviations are:", ess)
}
