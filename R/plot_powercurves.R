#' Plots the output of test_twoway_nrange
#'
#' Internal function, called by test_twoway_nrange, to plot power against sample size.
#'
#' @param power_over_nrange data.frame with sample sizes and corresponding powers to be plotted
#' @param target_power Numeric. Desired power to be attained. Accepts values between 0 and 1, defaults to 0.8.
#' @param title Character. Title for the graph. Defaults to 'Power curve from exact ANOVA test'
#' @param target_line Logical. If FALSE no target line will be drawn. Overrides target_power. Default is TRUE.
#' @param alpha Numeric. Type I error rate.
#' @param alpha_line Logical. Should a dashed line at the set alpha level be drawn. Default is TRUE.
#'
#' @return Plot with power curves.
#'
#' @examples
#'
#' ## 'cornorm_model' is created with the calculate_mean_matrix function
#' refmean <- 10
#' treateff <- 1.2
#' timeeff <- 0.75
#'
#' treatgroups <- 3
#' treatgroups_names <- c("wt", "DrugA", "DrugB")
#'
#' timepoints <- 4
#' timepoints_names <- paste0("t", 1:timepoints)
#'
#' nameslist <- list(treatment=treatgroups_names, time=timepoints_names)
#'
#' rho = 0.7
#'
#' cornorm_model <- calculate_mean_matrix(refmean = refmean, fAeffect = treateff, fBeffect = timeeff,
#' nlfA = treatgroups, nlfB = timepoints,
#' rho = rho, withinf = "fB", label_list = nameslist)
#'
#' nset <- seq(7, 14, 2)
#' cornorm_sim <- simulate_twoway_nrange(cornorm_model, nset, repeated_measurements=TRUE, nsims=5)
#'
#' ##used small number of iterations to reduce computation time
#'
#' power_results <- test_power_overkn(cornorm_sim, test="rank", plot=TRUE)
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#'
#' @importFrom scales pretty_breaks
#'
#' @export
plot_powercurves <- function(power_over_nrange, target_power = NULL, title = NULL, target_line=TRUE, alpha_line=TRUE, alpha=0.05)
{
  if(!is.data.frame(power_over_nrange))
  {
    stop("The 'power_over_nrange' object should be a 'data.frame'")
  }
  if(is.null(target_power))
  {
    target_power = 0.8
  }
  if(is.null(title))
  {
    title = "Power curve from two-way factorial simulation"
  }
  if(target_power<0 | target_power>1)
  {
    stop("The target power parameter should be between 0 and 1")
  }
  if(!is.character(title))
  {
    stop("title should be a character string")
  }
  if("n" %in% names(power_over_nrange) & ncol(power_over_nrange)<=5)
  {
    p <- ggplot2::ggplot(power_over_nrange, ggplot2::aes(x=n, y=power, group=effect, color=effect))
    xlabel <- expression(paste(italic(n), " per group"))
  } else if ("mean.group.size" %in% names(power_over_nrange) & ncol(power_over_nrange)==7)
  {
    p <- ggplot2::ggplot(power_over_nrange, ggplot2::aes(x=mean.group.size, y=power, group=effect, color=effect))
    xlabel <- expression(paste("Mean ", italic(n), " per group"))
  }
  p <- p + ggplot2::geom_line(linewidth=1.5, position = ggplot2::position_dodge(0.2)) + ggplot2::geom_point(size=2.4, position = ggplot2::position_dodge(0.2))
  if(all(c("lower.bound.ci", "upper.bound.ci") %in% names(power_over_nrange)))
  {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=lower.bound.ci, ymax=upper.bound.ci), linewidth=1.2, width=0.2, position = ggplot2::position_dodge(0.2))
    ylab <- "Power (95% confidence interval)"
  } else if (!"lower.bound.ci"%in% names(power_over_nrange))
  {
    ylab <- "Power"
  }
  p <- p + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(length(unique(power_over_nrange$n)))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(col="Effect", title=title, y= ylab, x=xlabel)
  if(target_line)
  {
    p <- p + ggplot2::geom_hline(yintercept = target_power, linetype="dashed", color = "red", linewidth=1.1)
  }
  if(alpha_line)
  {
    p <- p + ggplot2::geom_hline(yintercept = alpha, linetype="dashed", color = "yellowgreen", linewidth=1.1)
  }
  p <- p + ggplot2::theme_minimal()
  return(p)
}
