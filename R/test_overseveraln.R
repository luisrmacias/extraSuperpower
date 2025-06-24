#' Test simulated two-way factorial design experiments over different sample sizes.
#'
#' Wrapper to test data simulated under independent or repeated measurements and under different outcome distributions
#' with different sample sizes.
#'
#' @param data data.frame - data.frame with modeled outcome values, factor level labels, iteration number and sample size.
#' @param test character - Statistical test to be applied, possible values are 'ANOVA', 'rank' and 'permutation'.
#' @param plot logical - Should the power curve be plotted. Default is TRUE.
#' @param target_power Desired power to be attained. Accepts values between 0 and 1, defaults to 0.8.
#' @param title Title for the graph. Defaults to 'Power curve from exact ANOVA test'
#' @param target_line Set to TRUE. If FALSE no target line will be drawn. Overrides target_power.
#' @param alpha_line - logical Should a line at the set type I error be plotted
#' @param alpha - numeric Type I error probability
#'
#' @return Data frame with power and confidence intervals for the main effects and interaction for each of the sample sizes.
#' Also presented in graphical form if 'plot=TRUE'.
#'
#' @examples
#'
#'
#' ## In this example we simulate an independent sample design with skewed outcome
#' ## Model was specified with the 'calculate_mean_matrix function' (see ?calculate_mean_matrix)
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.25
#' timeeff <- 0.85
#' factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
#'
#' indep_matrix <- calculate_mean_matrix(refmean = refmean,
#'                                       fAeffect = treateff, fBeffect = timeeff,
#'                                       nlfA = treatgroups, nlfB = timepoints,
#'                                       label_list = factors_levels_names)
#'
#' indep_skewsim <- simulate_twoway_nrange(indep_matrix, seq(6, 12, 2),
#'                             distribution = "skewed", skewness = 1.8, nsims=5)
#' ##used low number of iterations to reduce computation time
#'
#' test_power_overkn(indep_skewsim, test="rank")
#'
#' @export
test_power_overkn <- function(data, test="ANOVA", plot=TRUE, target_power = NULL, title = NULL, target_line=TRUE, alpha_line=TRUE, alpha=0.05)
{
  res <- lapply(data, twoway_simulation_testing, test=test, alpha=alpha)
  res <- do.call(rbind, res)
  res$effect <- rownames(res)[1:3]
  if (plot)
  {
    p <- plot_powercurves(res, target_power=target_power, title=title, target_line=target_line, alpha_line=alpha_line, alpha=alpha)
    list(power_table = res, power_curve = p)
  } else if (!plot)
    res
}

