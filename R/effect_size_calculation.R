#' Effect size calculation
#'
#' Calculate effect sizes for two-way factorial designs from matrices of expected mean and standard deviation values
#' at each combination of factor levels. The output given is Cohen's f. Calculations are done as exemplified in the
#' G*Power 3.1 manual.
#'
#' @param matrices_obj List of 2 matrices, named mean.mat and sd.mat. This is the minimal output of the `calculate_mean_matrix` function. The full output from `calculate_mean_matrix` is also valid.
#'
#' @return Vector of length 3. The first two elements are the effect sizes for the main effects factor A and factor B,
#' respectively. The third element is the interaction effect size.
#'
#' @examples
#'
#' # no interaction effect expected
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#'
#' factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
#'
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean,
#'                                             nlfA = treatgroups, nlfB = timepoints,
#'                                             fAeffect = treateff, fBeffect = timeeff,
#'                                             label_list = factors_levels_names)
#' effsize(effects_treat_time)
#'
#' # we add cell specific interaction effect keeping design and main effect coefficients
#' cellswithinteraction <- matrix(c(rep(2,3), 3:5), 3,2)
#' #second level of factor A interacts with 3rd, 4th and 5th level of factor B
#'
#' effects_treat_time_interact <- calculate_mean_matrix(refmean = refmean,
#'                                                      nlfA = treatgroups, nlfB = timepoints,
#'                                                      fAeffect = treateff, fBeffect = timeeff,
#'                                                      label_list = factors_levels_names,
#'                                                      groupswinteraction = cellswithinteraction,
#'                                                      interact=1.3)
#' effsize(effects_treat_time_interact)
#'
#'
#' @export
effsize <- function(matrices_obj)
  {
  if(length(matrices_obj)==2 & any(unlist(sapply(matrices_obj, class))=="gg"))
  {
    matrices_obj <- matrices_obj$matrices_obj
  }
  if(any(!c("mean.mat", "sd.mat") %in% names(matrices_obj)) | !is.list(matrices_obj))
  {
   stop("A list with matrices called 'mean.mat' and 'sd.mat' is required")
  }
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
  ess <- t(data.frame(ess))
  rownames(ess) <- "Cohen's f"
  return(ess)
}
