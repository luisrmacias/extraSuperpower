#' Plots the output of exact_twoway_anova_power
#'
#' Creates a power curve from the output of the exact_twoway_anova_power function.
#'
#' @param exact_anova_power List that contains a data.frame with sample sizes and corresponding powers to be plotted.
#' @param target_power Desired power to be attained. Accepts values between 0 and 1, defaults to 0.8.
#' @param title Title for the graph. Defaults to 'Power curve from exact ANOVA test'
#' @param target_line Set to TRUE. If FALSE no target line will be drawn. Overrides target_power.
#'
#' @return Plot with power curves.
#'
#' @examples
#'
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' cellswithinteraction <- matrix(c(rep(2,3), 3:5), 3,2) #second level of factor A interacts with 3rd, 4th and 5th level of factor B
#' effects_treat_time_interact <- mean_sd_matrix(refmean = refmean, nlfA = treatgroups, nlfB = timepoints, fAeffect = treateff, fBeffect = timeeff, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints),
#'                                               groupswinteraction = cellswithinteraction, interact=1.3)
#' fxs <- anoveff(effects_treat_time_interact)

#' powerresult <- exact_twoway_anova_power(treatgroups, timepoints, fxs, 5:20)

#' plot_powercurve(powerresult, target_line = 0.9)
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
plot_powercurve <- function(exact_anova_power, target_power = NULL, title = NULL, target_line=TRUE)
{
  if(is.null(target_power))
  {
    target_power = 0.8
  }
  if(is.null(title))
  {
    title = "Power curve from exact ANOVA test"
  }
  if(target_power<0 | target_power>1)
  {
    stop("The target power parameter should be between 0 and 1")
  }
  if(!is.character(title))
  {
    stop("title should be a character string")
  }
  plotdata <- reshape2::melt(exact_anova_power$powercurve, id.vars="n", variable.name="Effect")
  names(plotdata)[3] <- "Power"

  p <- ggplot2::ggplot(plotdata, ggplot2::aes(x=n, y=Power, group=Effect, color=Effect)) + ggplot2::geom_line(linewidth=1.5) + ggplot2::geom_point(size=2.4)
  p <- p + ggplot2::scale_x_continuous(
    breaks = scales::pretty_breaks(length(unique(plotdata$n)))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(col="Effect", title=title, y="Power", x=expression(paste(italic(n), " per group")))
  if(target_line)
  {
    p <- p + ggplot2::geom_hline(yintercept = target_power, linetype="dashed", color = "red", linewidth=1.1)
  }
  return(p)
}
