#' Plots the output of test_twoway_nrange
#'
#' Internal function, called by test_twoway_nrange, to plot power against sample size.
#'
#' @param power_over_nrange data.frame with sample sizes and corresponding powers to be plotted
#' @param target_power Desired power to be attained. Accepts values between 0 and 1, defaults to 0.8.
#' @param title Title for the graph. Defaults to 'Power curve from exact ANOVA test'
#' @param target_line Set to TRUE. If FALSE no target line will be drawn. Overrides target_power.
#'
#' @return Plot with power curves.
#'
#' @examples
#'
#' ## 'cornorm_model' is created with the calculate_mean_matrix function
#'
#' nrange <- 12:20
#' cornorm_sim <- simulate_twoway_nrange(nrange, cornorm_model)
#'
#' power_results <- test_twoway_nrange(cornorm_sim, test="rank", plot=TRUE)
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom scales pretty_breaks
#'
#' @export
plot_powercurves <- function(power_over_nrange, target_power = NULL, title = NULL, target_line=TRUE, alpha_line=TRUE, alpha=0.05)
{
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

  p <- ggplot2::ggplot(power_over_nrange, ggplot2::aes(x=n, y=power, group=effect, color=effect)) + ggplot2::geom_line(linewidth=1.5)
  p <- p + ggplot2::geom_point(size=2.4, position = position_dodge(0.2)) + ggplot2::geom_errorbar(ggplot2::aes(ymin=lower.bound.ci, ymax=upper.bound.ci), linewidth=1.2, width=0.2, , position = position_dodge(0.2))
  p <- p + ggplot2::scale_x_continuous(
    breaks = scales::pretty_breaks(length(unique(power_over_nrange$n)))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(col="Effect", title=title, y="Power (95% confidence interval)", x=expression(paste(italic(n), " per group")))
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
