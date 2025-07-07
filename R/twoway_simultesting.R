#' Calculate power for global main effects and interaction from two-way factorial simulated data
#'
#' This functions takes the output of either the 'twoway_simulation_independent' or the 'twoway_simulation_correlated' functions
#' and calculates the power of the sample size used in the simulation under parametric analysis of variance, rank based analysis of variance or
#' permutation testing.
#'
#'
#' @param data - Simulation obtained from the 'twoway_simulation_independent' or the 'twoway_simulation_correlated' functions.
#' @param test - The test to be applied. Possible values are "ANOVA" (default), "rank" and "permutation".
#' @param alpha - Type I error rate. Default is 0.05.
#'
#' @importFrom utils capture.output
#' @importFrom stats as.formula
#' @importFrom stats qnorm
#'
#' @return A data.frame with the power and 95% confidence interval for each of the main effects and their interaction.
#'
#' @examples
#'
#' ## After creating a 'matrices_obj' with the 'calculate_mean_matrix' function.
#'
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' rho <- 0.8
#' withinf <- "fB"
#' factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
#'
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean,
#'                                             fAeffect = treateff, fBeffect = timeeff,
#'                                             nlfA = treatgroups,  nlfB = timepoints,
#'                                             rho = rho, withinf = withinf,
#'                                             label_list = factors_levels_names)
#'
#' n <- 7
#' correlated_sim <- twoway_simulation_correlated(group_size=n, matrices_obj=effects_treat_time,
#'                   nsims=20)
#' ##used smaller number of iterations to reduce computation time
#'
#' twoway_simulation_testing(correlated_sim)
#' ## defaults to parametric analysis of variance
#'
#' twoway_simulation_testing(correlated_sim, test="rank")
#' ## rank based analysis of variance
#'
#' ## permutation test is another option
#'
#'
#' @export
twoway_simulation_testing <- function(data, test="ANOVA", alpha=0.05)
{
  if(is.list(data) & is.null(dim(data)))
  {
    # checkFunction <- function() {
    #   # check if called in testthat or interactive
    #   tb <- .traceback(x = 0)
    #   if(!any(unlist(lapply(tb, function(x) any(grepl("test_env", x))))) && interactive())
    #   {
    #     cat("Permutation testing for repeated measurement designs can take several minutes\n
    #     depending on sample size and number of groups.")
    #     user_input <- readline("Do you wish to proceed? (y/n)  ")
    #   if(user_input != 'y') stop('Exiting')
    #   }
    #}

    withinf <- data$withinf
    simulation <- data$simulated_data
    message(paste("Testing power on a repeated observations design experiment.\n Sample size =", unique(simulation$n),"\n"))
    indep_vars <- names(simulation)[4:5]
    names(simulation)[4:5] <- c("indep_var1", "indep_var2")
    simulation <- split(simulation, simulation$iteration)
    if(withinf=="fA")
    {
      if(test=="ANOVA")
      {
        pvec <- sapply(simulation, function(i)
          suppressMessages(suppressWarnings(afex::aov_ez(id="subject", dv="y", within="indep_var1", between ="indep_var2", data=i)$anova_table))$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = "indep_var1", between = "indep_var2",  data = simulation[[1]])$anova_table))
        pvec <- pvec[c(2,1,3),]
        pvecnames <- pvecnames[c(2,1,3)]
        message("Performing ANOVA testing on simulated data")
      } else if(test=="permutation")
      {
        ##cat(paste('Permutation testing with n=', mean(group_size), 'starts'))
        fmla <- as.formula("y ~ indep_var1*indep_var2+ Error(subject/indep_var1)")
        pvec <- sapply(simulation,
                       function(i) permuco::aovperm(fmla, data = i)$table$`resampled P(>F)`)
        pvec <- pvec[c(2,1,3),]
        pvecnames <- rownames(permuco::aovperm(fmla, data = simulation[[1]])$table)
        pvecnames <- pvecnames[c(2,1,3)]
        message("Performing permutation testing on simulated data")
      } else if(test=="rank")
      {
        pvec <- NULL
        for (i in seq(simulation))
        {
          capture.output(res <- nparLD::f1.ld.f1(y=simulation[[i]]$y, time = simulation[[i]]$indep_var1, group = simulation[[i]]$indep_var2, subject = simulation[[i]]$subject,
                                                 plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE)$ANOVA.test[,3],
                         file = nullfile())
          pvec <- cbind(pvec, res)
        }
        pvec <- pvec[c(2,1,3),]
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
        message("Performing rank testing on simulated data")
      }
    } else if (withinf=="fB")
    {
      if(test=="ANOVA")
      {
        pvec <- sapply(simulation, function(i)
          suppressMessages(suppressWarnings(afex::aov_ez(id="subject", dv="y", between="indep_var1", within="indep_var2", data=i)$anova_table))$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", between = "indep_var1", within = "indep_var2",  data = simulation[[1]])$anova_table))
        message("Performing ANOVA testing on simulated data")
      } else if(test=="permutation")
      {
        ##checkFunction()
        ##cat(paste('Permutation testing with n=', mean(group_size), 'starts'))
        fmla <- as.formula("y ~ indep_var1*indep_var2+ Error(subject/indep_var2)")
        pvec <- sapply(simulation,
                       function(i) permuco::aovperm(fmla, data = i)$table$`resampled P(>F)`)
        pvecnames <- rownames(permuco::aovperm(fmla, data = simulation[[1]])$table)
        message("Performing permutation testing on simulated data")
      } else if(test=="rank")
      {
        pvec <- NULL
        for (i in seq(simulation))
        {
          capture.output(res <- nparLD::f1.ld.f1(y=simulation[[i]]$y, time = simulation[[i]]$indep_var2, group = simulation[[i]]$indep_var1, subject = simulation[[i]]$subject,
                                                 plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE)$ANOVA.test[,3],
                         file = nullfile())
          pvec <- cbind(pvec, res)
        }
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
        message("Performing rank testing on simulated data")
      }
    } else if (withinf=="both")
    {
      if(test=="ANOVA")
      {
        pvec <- sapply(simulation, function(i)
          suppressMessages(suppressWarnings(afex::aov_ez(id="subject", dv="y", within=c("indep_var1", "indep_var2"), data=i)$anova_table))$`Pr(>F)`)
        pvecnames <- rownames(suppressMessages(afex::aov_ez(id = "subject", dv = "y", within = c("indep_var1", "indep_var2"),  data = simulation[[1]])$anova_table))
        message("Performing ANOVA testing on simulated data")
      }else if(test=="permutation")
      {
        ##checkFunction()
        ##cat(paste('Permutation testing with n=', mean(group_size), 'starts'))
        fmla <- as.formula("y ~ indep_var1*indep_var2+ Error(subject/indep_var1 + indep_var2)")
        pvec <- sapply(simulation, function(i) permuco::aovperm(fmla, data = i)$table$`resampled P(>F)`)
        pvecnames <- rownames(permuco::aovperm(fmla, data = simulation[[1]])$table)
        message("Performing permutation testing on simulated data")
      } else if(test=="rank")
      {
        pvec <- NULL
        for (i in seq(simulation))
        {
          capture.output(res <- nparLD::ld.f2(y=simulation[[i]]$y, time1 = simulation[[i]]$indep_var1, time2 = simulation[[i]]$indep_var2, subject = simulation[[i]]$subject,
                                                 plot.RTE = FALSE, order.warning = FALSE, description = FALSE, show.covariance = FALSE)$ANOVA.test[,3],
                         file = nullfile())
          pvec <- cbind(pvec, res)
        }
        pvecnames <- c("indep_var1", "indep_var2", "indep_var1:indep_var2" )
        message("Performing rank testing on simulated data")
      }
    }
    pvecnames <- gsub("indep_var1", indep_vars[1], pvecnames)
    pvecnames <- gsub("indep_var2", indep_vars[2], pvecnames)
  } else if (is.data.frame(data))
  {
    indep_vars <- names(data)[4:5]
    simulation <- split(data, data$iteration)
    if(length(unique(simulation[[1]]$n))==1)
    {
      message(paste("Testing power on an independent observations design experiment.\nSample size per group =", unique(simulation[[1]]$n), "\n"))
    } else if (length(unique(simulation[[1]]$n))>1)
    {
      message(paste("Testing power on an independent observations design experiment.\nMean sample size per group =", round(mean(simulation[[1]]$n),1), "\n"))
    }

    if(test=="ANOVA")
    {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2], "+ Error(1|subject)"))
      pvec <- sapply(simulation,
                     function(i) suppressMessages(afex::aov_car(frml, i)$anova_table$`Pr(>F)`))
      pvecnames <- rownames(suppressMessages(afex::aov_car(frml, simulation[[1]])$anova_table))
      message("Performing ANOVA testing on simulated data")
    } else if(test=="rank")
    {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2]))
      pvec <- sapply(simulation,
                     function(i) Rfit::raov(frml, i)$table[,5])
      pvecnames <- rownames(pvec)
      message("Performing rank testing on simulated data")
    } else if(test=="permutation")
    {
      frml <- as.formula(paste("y ~ ", indep_vars[1], "*", indep_vars[2]))
      pvec <- sapply(simulation, function(i) {
                       res <- permuco::aovperm(frml, data = i)
                       res <- res$table[-nrow(res$table),]
                       res <- res$`resampled P(>F)`
                       })
      pvecnames <- rownames(permuco::aovperm(frml, simulation[[1]])$table)
      pvecnames <- pvecnames[-grep("Residuals", pvecnames)]
      message("Performing permutation testing on simulated data")
    }
  }
  pprops <- rowSums(pvec<alpha)/ncol(pvec)
  lb <- round(pprops - qnorm(1-(0.05/2))*sqrt((pprops*(1-pprops))/ncol(pvec)), 4)
  lb[lb<0] <- 0.0000
  ub <- round(pprops + qnorm(1-(0.05/2))*sqrt((pprops*(1-pprops))/ncol(pvec)), 4)
  ub[ub<0] <- 0.0000

  if(length(unique(simulation[[1]]$n))==1) {
    n <- unique(simulation[[1]]$n)
  } else if (length(unique(simulation[[1]]$n))>1)
  {
    n <- range(simulation[[1]]$n)
  }

  # if(test=="rank")
  # {
  #   if(length(unique(simulation[[1]]$n))==1)
  #   {
  #     data.frame(n = n, power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
  #   }
  #   else if (length(unique(simulation[[1]]$n))>1)
  #   {
  #     data.frame("smallest group" = n[1], "largest group" = n[2], "mean group size" = mean(n), power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
  #   }
  # } else if (test!="rank")
  # {
    names(pprops) <- pvecnames
    if(length(unique(simulation[[1]]$n))==1)
    {
      data.frame(n = n, power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
    }
    else if (length(unique(simulation[[1]]$n))>1)
    {
      data.frame("smallest group" = n[1], "largest group" = n[2], "mean group size" = mean(n), power=pprops, "lower bound ci" = lb, "upper bound ci" = ub)
    }
  #}
}
