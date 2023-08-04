#' Modification of the ANOVA_power function from the Superpower package in order to obtain the simulated values
#' from which power calculations are performed.
#'
#' @param design_result Output from the ANOVA_design function
#' @param nsims Number of simulations to perform  
#' @export
anova_dat <- function (design_result, alpha_level = Superpower_options("alpha_level"), 
          correction = Superpower_options("correction"), p_adjust = "none", 
          nsims = 1000, seed = NULL, verbose = Superpower_options("verbose"), 
          emm = Superpower_options("emm"), emm_model = Superpower_options("emm_model"), 
          contrast_type = Superpower_options("contrast_type"), emm_p_adjust = "none", 
          emm_comp = NULL) 
{
  cohen_f <- partial_eta_squared <- non_centrality <- NULL
  if (is.element(p_adjust, c("holm", "hochberg", "hommel", 
                             "bonferroni", "BH", "BY", "fdr", "none")) == FALSE) {
    stop("p_adjust must be of an acceptable adjustment method: see ?p.adjust")
  }
  if (is.element(correction, c("none", "GG", "HF")) == FALSE) {
    stop("Correction for sphericity can only be none, GG, or HF")
  }
  if (nsims < 10) {
    stop("The number of repetitions in simulation must be at least 10; suggested at least 1000 for accurate results")
  }
  if (!is.null(seed)) {
    sysSeed <- .GlobalEnv$.Random.seed
    on.exit({
      if (!is.null(sysSeed)) {
        .GlobalEnv$.Random.seed <- sysSeed
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  }
  if (missing(alpha_level)) {
    alpha_level <- 0.05
  }
  if (alpha_level >= 1 | alpha_level <= 0) {
    stop("alpha_level must be less than 1 and greater than zero")
  }
  design <- design_result$design
  factornames <- design_result$factornames
  n <- design_result$n
  mu = design_result$mu
  sd <- design_result$sd
  r <- design_result$r
  factors <- design_result$factors
  design_factors <- design_result$design_factors
  sigmatrix <- design_result$sigmatrix
  dataframe <- design_result$dataframe
  design_list <- design_result$design_list
  if (grepl("w", design_result$design) == TRUE && length(unique(design_result$n)) > 
      1) {
    stop("Unequal group sizes are not possible when the design contains within factors")
  }
  n_vec <- n
  n <- max(n)
  frml1 <- design_result$frml1
  frml2 <- design_result$frml2
  aov_result <- suppressMessages({
    afex::aov_car(frml1, data = dataframe, include_aov = FALSE, 
            anova_table = list(es = "pes", p_adjust_method = p_adjust))
  })
  pairs_result_df = NULL
  possible_pc <- (((prod(as.numeric(strsplit(design, "\\D+")[[1]])))^2) - 
                    prod(as.numeric(strsplit(design, "\\D+")[[1]])))/2
    sim_data <- as.data.frame(matrix(ncol = 2 * (2^factors - 
                                                   1) + 2 * possible_pc, nrow = nsims))
  paired_tests <- combn(unique(dataframe$cond), 2)
  paired_p <- numeric(possible_pc)
  paired_d <- numeric(possible_pc)
  within_between <- sigmatrix[lower.tri(sigmatrix)]
  names(sim_data) = c(paste("anova_", rownames(aov_result$anova_table), 
                              sep = ""), paste("anova_es_", rownames(aov_result$anova_table), 
                                               sep = ""), paste("p_", paste(paired_tests[1, ], paired_tests[2, 
                                               ], sep = "_"), sep = ""), paste("d_", paste(paired_tests[1, 
                                               ], paired_tests[2, ], sep = "_"), sep = ""))
  
  for (i in 1:nsims) {
    dataframe <- design_result$dataframe
    dataframe$y <- suppressMessages({
      melt(as.data.frame(MASS::mvrnorm(n = n, mu = mu, Sigma = as.matrix(sigmatrix))))$value
    })
    if (i==1) 
    {
      data_sim <- dataframe
      data_sim$sim <- i
    }else if (i>1)
    {
      dataframe$sim <- i
      data_sim <- rbind(data_sim, dataframe)
    }
    if (length(n_vec) > 1) {
      for (k in 1:length(unique(dataframe$cond))) {
        if ((n - n_vec[k]) > 0) {
          dataframe <- dataframe[-sample(which(dataframe$cond == 
                                                 unique(dataframe$cond)[k]), (n - n_vec[k])), 
          ]
        }
      }
    }
    aov_result <- suppressMessages({
      aov_car(frml1, data = dataframe, include_aov = FALSE, 
              anova_table = list(es = "pes", p_adjust_method = p_adjust, 
                                 correction = correction))
    })
     sim_data[i, ] <- c(aov_result$anova_table[[6]], aov_result$anova_table[[5]], 
                        p.adjust(paired_p, method = p_adjust), paired_d)
  }
  plotData <- suppressMessages(melt(sim_data[1:(2^factors - 
                                                  1)], value.name = "p"))
  SalientLineColor <- "#535353"
  LineColor <- "Black"
  BackgroundColor <- "White"
  p <- plotData$p
  swr = function(string, nwrap = 10) {
    paste(strwrap(string, width = 10), collapse = "\n")
  }
  swr = Vectorize(swr)
  plotData$variable = swr(chartr("_:", "  ", plotData$variable))
  plt1 = ggplot(plotData, aes(x = p)) + scale_x_continuous(breaks = seq(0, 
                                                                        1, by = 0.1), labels = seq(0, 1, by = 0.1)) + geom_histogram(colour = "black", 
                                                                                                                                     fill = "white", breaks = seq(0, 1, by = 0.01)) + geom_vline(xintercept = alpha_level, 
                                                                                                                                                                                                 colour = "red") + facet_grid(variable ~ .) + labs(x = "p") + 
    theme_bw()
  p_paired <- sim_data[(2 * (2^factors - 1) + 1):(2 * (2^factors - 
                                                         1) + possible_pc)]
  plotData <- suppressMessages(melt(p_paired, value.name = "p"))
  p <- plotData$p
  plotData$variable = swr(chartr("_:", "  ", plotData$variable))
  plt2 = ggplot(plotData, aes(x = p)) + scale_x_continuous(breaks = seq(0, 
                                                                        1, by = 0.1), labels = seq(0, 1, by = 0.1)) + geom_histogram(colour = "black", 
                                                                                                                                     fill = "white", breaks = seq(0, 1, by = 0.01)) + geom_vline(xintercept = alpha_level, 
                                                                                                                                                                                                 colour = "red") + facet_grid(variable ~ .) + labs(x = expression(p)) + 
    theme_bw()
  power = as.data.frame(apply(as.matrix(sim_data[(1:(2^factors - 
                                                       1))]), 2, function(x) mean(ifelse(x < alpha_level, 1, 
                                                                                         0) * 100)))
  es = as.data.frame(apply(as.matrix(sim_data[((2^factors):(2 * 
                                                              (2^factors - 1)))]), 2, function(x) mean(x)))
  main_results <- data.frame(power, es)
  names(main_results) = c("power", "effect_size")
  power_paired = as.data.frame(apply(as.matrix(sim_data[(2 * 
                                                           (2^factors - 1) + 1):(2 * (2^factors - 1) + possible_pc)]), 
                                     2, function(x) mean(ifelse(x < alpha_level, 1, 0) * 100)))
  es_paired = as.data.frame(apply(as.matrix(sim_data[(2 * (2^factors - 
                                                             1) + possible_pc + 1):(2 * (2^factors - 1) + 2 * possible_pc)]), 
                                  2, function(x) mean(x)))
  pc_results <- data.frame(power_paired, es_paired)
  names(pc_results) = c("power", "effect_size")
    emm_results = NULL
  manova_result = NULL
  data_sim
}
