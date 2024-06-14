#' Calculate power for global main effects and interaction from two-way factorial simulated data
#'
#' This functions takes the output of either the 'twoway_simulation_independent' or the 'twoway_simulation_correlated' functions
#' and calculates the power of the sample size used in the simulation under parametric analysis of variance, rank based analysis of variance or
#' permutation testing.
#'
#' Permutation testing under a repeated measurements design may be very time consuming.
#'
#' @param data - Simulation obtained from the 'twoway_simulation_independent' or the 'twoway_simulation_correlated' functions.
#' @param test - The test to be applied. Possible values are "ANOVA" (default), "rank" and "permutation".
#' @param alpha - Type I error rate. Default is 0.05.
#'
#' @return A data.frame with the power and 95% confidence interval for each of the main effects and their interaction.
#'
#' @examples ## After creating a 'matrices_obj' with the 'calculate_mean_matrix' function.
#'
#' n <- 7
#' correlated_sim <- twoway_simulation_correlated(group_size=n, matrices_obj=matrices_obj)
#' ## defaults to 500 iterations
#'
#' twoway_simulation_testing(correlated_sim)
#' ## defaults to parametric analysis of variance
#'
#' twoway_simulation_testing(correlated_sim, test="rank")
#' ## rank based analysis of variance
#'
#' twoway_simulation_testing(correlated_sim, test="permutation")
#' ## warning, permutation testing taked 25-40 seconds per iteration
#'
#' @export
twoway_simulation_testing <- function(data, test="ANOVA", alpha=0.05)
  {
  require(lmPerm)
  if(is.list(data) & is.null(dim(data)))
  {
    withinf <- data$withinf
    simulation <- data$simulated_data
    indep_vars <- names(simulation)[4:5]
    names(simulation)[4:5] <- c("indep_var1", "indep_var2")
    simulation <- split(simulation, simulation$iteration)
    if(withinf=="fA")
    {
      if(test=="ANOVA")
        {
        pvec <- sapply(simulation, function(i)
          suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = "indep_var1", between = "indep_var2",  data = i))$anova_table$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = "indep_var1", between = "indep_var2",  data = simulation[[1]])$anova_table))
      } else if(test=="permutation")
      {
        pvec <- sapply(simulation, function(i) ez::ezPerm(wid=subject, dv = y, within = indep_var1, between = indep_var2, data = i)$p)
        pvecnames <- ez::ezPerm(wid=subject, dv = y, within = indep_var1,  between = indep_var2, data = simulation[[1]])$Effect
      } else if(test=="rank")
      {
        pvec <- sapply(simulation, function(i) nparLD::f1.ld.f1(y=i$y, time = i$indep_var1, group = i$indep_var2, subject = i$subject,
                                                                plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE)$ANOVA.test[,3])
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
      }
    } else if (withinf=="fB")
    {
      if(test=="ANOVA")
      {
        pvec <- sapply(simulation, function(i)
          suppressMessages(afex::aov_ez(id = "subject", dv = "y", between = "indep_var1", within = "indep_var2", data = i))$anova_table$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", between = "indep_var1", within = "indep_var2",  data = simulation[[1]])$anova_table))

      } else if(test=="permutation")
      {
        pvec <- sapply(simulation,
                       function(i) ez::ezPerm(wid=subject, dv = y, between = indep_var1, within = indep_var2, data=i)$p)
        pvecnames <- ez::ezPerm(wid=subject, dv = y, between = indep_var1, within = indep_var2, data = simulation[[1]])$Effect
      } else if(test=="rank")
      {
        pvec <- sapply(simulation, function(i) nparLD::f1.ld.f1(y=i$y, time = i$indep_var2, group = i$indep_var1, subject = i$subject,
                                                                 plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE)$ANOVA.test[,3])
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
      }
    } else if (withinf=="both")
    {
      if(test=="ANOVA")
      {
        pvec <- sapply(simulation, function(i)
          suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = c("indep_var1", "indep_var2"), data = i))$anova_table$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = c("indep_var1", "indep_var2"),  data = simulation[[1]])$anova_table))
      }else if(test=="permutation")
      {
        pvec <- sapply(simulation,
                       function(i) ez::ezPerm(wid=subject, dv = y, within = .(indep_var1, indep_var2),  data = i)$p)
        pvecnames <- ez::ezPerm(wid=subject, dv = y, within = .(indep_var1, indep_var2), data = simulation[[1]])$Effect
      } else if(test=="rank")
      {
        pvec <- sapply(simulation, function(i) suppressWarnings(nparLD::ld.f2(y=i$y, time1 = i$indep_var1, time2 = i$indep_var2, subject = i$subject,
                                                             plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE))$ANOVA.test[,3])
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
      }
    }
    pvecnames <- gsub("indep_var1", indep_vars[1], pvecnames)
    pvecnames <- gsub("indep_var2", indep_vars[2], pvecnames)
    } else if (is.data.frame(data))
  {
    indep_vars <- names(data)[4:5]
    simulation <- split(data, data$iteration)

    if(test=="ANOVA")
      {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2], "+ Error(1|subject)"))
      pvec <- sapply(simulation,
                     function(i) suppressMessages(afex::aov_car(frml, i)$anova_table$`Pr(>F)`))
      pvecnames <- rownames(suppressMessages(afex::aov_car(frml, simulation[[1]])$anova_table))
    } else if(test=="rank")
    {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2]))
      pvec <- sapply(simulation,
                     function(i) Rfit::raov(frml, i)$table[,5])
    } else if(test=="permutation")
    {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2]))
      pvec <- sapply(simulation,
                     function(i) {
                       coefs <- summary(lmPerm::aovp(frml, data = i))
                       coefs[[1]][1:3,5]
                       })
      coefs <- summary(lmPerm::aovp(frml, data=simulation[[1]]))
      pvecnames <- stringr::str_trim(rownames(coefs[[1]])[1:3])
    }
  }
  pprops <- rowSums(pvec<alpha)/ncol(pvec)
  lb <- round(pprops - qnorm(1-(0.05/2))*sqrt((pprops*(1-pprops))/ncol(pvec)), 4)
  lb[lb<0] <- 0.0000
  ub <- round(pprops + qnorm(1-(0.05/2))*sqrt((pprops*(1-pprops))/ncol(pvec)), 4)
  ub[ub<0] <- 0.0000
  if(test=="rank")
  {
    data.frame(power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
  } else if (test!="rank")
  {
    names(pprops) <- pvecnames
    data.frame(power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
  }
  }
