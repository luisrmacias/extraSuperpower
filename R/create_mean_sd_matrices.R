#' Create input for simulation based two-way factorial experiments
#'
#' This function will generate a matrix of expected mean values for *ab* factor level combinations of a two-way
#' factorial design by assuming linear effects with possible departure from linearity by interaction. It will also
#' provide a standard deviation matrix for these *ab* combinations of factor levels. If the design has repeated measures,
#' it will additionally provide correlation and covariance matrices calculated depending on which factor has repeated
#' measurements or is the 'within' factor.
#'
#' @param refmean Numeric - expected mean for first level of both factors A and B
#' @param nlfA Integer - number of levels of factor A
#' @param nlfB Integer - number of levels of factor B
#' @param fAeffect Numeric - multiple that defines how cell means are modified by factor A. With the default `endincrement` (`TRUE`), determines the last level of factor A with respect to its first level. When `endincrement=FALSE` this multiple applies from one level to the next.
#' @param fBeffect Numeric - multiple that defines how cell means are modified by factor B. With the default `endincrement` (`TRUE`), determines the last level of factor B with respect to its first level. When `endincrement=FALSE` this multiple applies from one level to the next.
#' @param groupswinteraction Vector length 2 or n*2 matrix - Combination of levels from factors A and B in which interaction is expected
#' @param interact Numeric - value by which the mean from cell or cells indicated in `groupswinteraction` is multiplied after it has been calculated accordingly to `fAeffect` and `fBeffect`
#' @param label_list List length 2 - vectors with the names of the factor levels. The objects in this list should be named as the factors. The use of this option is encouraged as these names are used for plotting and inherited to downstream functions.
#' @param sdproportional Logical - whether the standard deviation for each combination of factor levels is a proportion of the respective factor level combination mean, defaults to `TRUE`
#' @param sdratio Numeric - value by which the expected mean value of a factor level combination is multiplied to obtain the respective standard deviation, defaults to 0.2.
#' @param endincrement Logical - determines if the multiples provided in `fAeffect` and `fBeffect` refer to change between first and last levels (default) or level to level changes.
#' @param rho Vector length 1 or 2, or 2 by 2 matrix - Controls how the correlation and hence de covariance matrix is built. See 'details' and `?gencorrelationmat` examples.
#' @param withinf Character - Names the factor with repeated measures. Possible values are NULL, "fA", "fB" or "both"
#' @param plot Logical - Should a line plot with the modeled mean and standard deviations be part of the output. Default is `TRUE`
#'
#' @return If `rho` and `withinf` are left at their default values of 0 and NULL, respectively, a cell mean matrix, a cell standard deviation matrix and
#' optionally a graph that represents both.
#' @return If `rho` is between -1 and 1 but different to 0 and `withinf` is either "fA", "fB" or "both", correlation and covariance matrices are generated
#' along with the aforementioned output.
#'
#' @details
#' The user must provide a reference mean (usually mean in control or untreated group), the expected change for each factor
#' from first to last level (or from one level to the next) and the number of levels in each factor.
#'
#' The user can also specify factor level combinations in which interaction is assumed and its magnitude with respect to the
#' reference mean. The cell mean matrix will be modified accordingly and this can also have an effect of the standard deviation
#' matrix.
#'
#' We were motivated by sample size calculation for two-way factorial designs with *1,2,...,a* levels of factor *A* and
#' *1,2,...,b* levels of factor *B* in which the mean outcome value for replicates of cell *A=1, B=1* are known.
#' Furthermore, there is an expected change in level mean for each of the factors. Finally, interaction can be explicitly
#' introduced to level combinations in which it is expected to occur.
#'
#' If a repeated measures experiment is intended `withinf` must be set to "fA", "fB" or "both", depending on which is the 'within' factor.
#' If `rho` is a vector length 1, the within subject correlation will be constant for the factor defined in `withinf`. If `rho` is a vector
#' length 2 and `withinf` is either "fA" or "fB" a correlation gradient will be created from the first to second value of `rho`. If `rho` is
#' a vector length 2 and `withinf="both"`, the first element of `rho` will be the correlation within factor A, while the second element will
#' be the correlation within factor B. If `rho` is a 2*2 matrix, only possible if `withinf="both"`, a correlation gradient will be created
#' across rows of `rho` for each of the factors.
#'
#'
#' @examples
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
#' ## Independent design
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean,
#'                                             fAeffect = treateff, fBeffect = timeeff,
#'                                             nlfA = treatgroups, nlfB = timepoints,
#'                                             label_list = factors_levels_names)
#' ## Inspect plot to check if matrices correspond to design
#' effects_treat_time$meansplot
#' n <- 20
#' independent_experiment <- twoway_simulation_independent(group_size = n,
#'                                       matrices_obj = effects_treat_time)
#'
#' head(independent_experiment, 10)
#'
#' ## Repeated measures design, suppose subjects from 4 independent treatment groups measured
#' ## at 5 different timepoints.
#' ## We use the same parameters as the independent design example, except we add within factor level
#' ## correlation and we specify that factor B is the within factor.
#'
#'
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' rho <- 0.8
#' withinf <- "fB"
#'
#'factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
#'
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean, fAeffect = treateff,
#'                       fBeffect = timeeff, nlfA = treatgroups,  nlfB = timepoints,
#'                       rho = rho, withinf = withinf, label_list = factors_levels_names)
#'
#' ## Plot should look the same, structure within data can be checked once simulated
#' effects_treat_time$meansplot
#'
#' n <- 20
#' repeatedmeasures_experiment <- twoway_simulation_correlated(group_size = n,
#'                                           matrices_obj = effects_treat_time)
#' head(repeatedmeasures_experiment, 10)
#'
#' @export
calculate_mean_matrix <- function(refmean, nlfA, nlfB, fAeffect, fBeffect, groupswinteraction=NULL, interact=1, label_list = NULL,
                                  sdproportional = TRUE, sdratio=0.2, endincrement=TRUE,  rho=0, withinf=NULL, plot=TRUE)
{
  if(is.null(withinf))
  {
    withinf <- "none"
  }
  ##Check if labels correspond to group sizes. If no labels are given, assign them.
  if(!is.null(label_list))
  {
    if(length(label_list)!=2)
    {
      stop("Label names must be a list of length 2")
    } else if (length(label_list[[1]]) != nlfA | length(label_list[[2]]) != nlfB)
    {
      stop("\nNumber of labels must match number of levels for each factor")
    }
  }
  if(is.null(label_list))
  {
    label_list <- list(fA = LETTERS[1:nlfA], fB = letters[1:nlfB])
  }
  if(any(abs(rho)>1))
  {
    stop("\nRho must be a number between -1 and 1.")
  }
  if(length(rho)==2 & withinf=="both")
  {
    message("The first element of 'rho' will be the correlation for factor A, the second element of 'rho' the correlation for factor B")
  }
  if(is.matrix(rho) & withinf!="both")
  {
    stop("'rho' can only be a matrix if both factor A and factor B are within factors. In that case 'within' should be set to 'both'")
  }
  if((fAeffect==0|fBeffect==0) & isTRUE(sdproportional))
  {
    if(withinf=="none")
    {
      warning("\nBy setting '0' effects and proportional SD you will obtain groups with standard deviation of 0.")
    }
    else if(withinf=="fA" | withinf=="fB" | withinf=="both")
    {
      warning("\nBy setting '0' effects and proportional SD you will obtain groups with standard deviation\nof 0 and individuals with 0 covariance.")
    }
  }

  # Bincrements <- fAvec[-1] - fAvec[1]
  # effmat <- t(sapply(1:length(Bincrements), function(x) fBvec[-1] + Bincrements[x]))
  # if(nlfB>2)
  # {
  #   effmat <- rbind(fBvec[-1], effmat)
  # } else if (nlfB==2)
  # {
  #   effmat <- c(fBvec[-1], effmat)
  # }
  # effmat <- cbind(fAvec, effmat)
  # dimnames(effmat) <- label_list

  ## Generation of mean matrix
  fAvec <- genvecs(change = fAeffect, reps = nlfA, bystart = endincrement, scaler = refmean)
  fBvec <- genvecs(change = fBeffect, reps = nlfB, bystart = endincrement, scaler = refmean)

  effmat <- build_mean_mat(fAvec = fAvec, fBvec = fBvec, iA = fAeffect, a = nlfA, b = nlfB,
                           label_list = label_list, bystart=endincrement)

  ## Modify mean matrix depending on interaction terms
  if(interact!=1)
  {
    if (!(is.matrix(groupswinteraction) | length(groupswinteraction)==2))
    {
      stop("\nPlease provide a vector length 2 or a n*2 matrix")
    }
    if (length(groupswinteraction==2) & is.vector(groupswinteraction))
    {
      groupswinteraction <- t(as.matrix(groupswinteraction))
    }
    for(x in seq(nrow(groupswinteraction)))
    {
      i <- groupswinteraction[x,1]
      j <- groupswinteraction[x,2]
      effmat[i,j] <- effmat[i,j]*interact
    }
  }
  mean_matrix <- effmat*refmean
  ## Generation of standard deviation matrix
  if (sdproportional)
  {
    sd_matrix <- abs(mean_matrix*sdratio)
  }else if (!sdproportional)
  {
    sd_matrix <- abs(mean(mean_matrix)*sdratio)
  }
  ## Function ends here if measurements are independent
  if(withinf=="none")
  {
    matrices_obj <- list(mean.mat = mean_matrix, sd.mat= sd_matrix)
  } else if(withinf!="none")
  {
    if(withinf!="fA" & withinf!="fB" & withinf!="both")
    {
      stop("Possible values for the 'withinf' parameter are 'fA', 'fB' or 'both'")
    }
    ## Generation of correlation and covariance matrices
    cormat <- suppressMessages(gencorrelationmat(mean_matrix = mean_matrix, rho = rho, withinf = withinf,
                                                   nlfA = nlfA, nlfB = nlfB))
    sigmat <- suppressMessages(gencovariancemat(correlation_matrix = cormat, sd_matrix = sd_matrix, withinf = withinf,
                               nlfA = nlfA, nlfB = nlfB))
    matrices_obj <- list(within.factor = withinf, mean.mat = mean_matrix, sd.mat= sd_matrix, cormat = cormat, sigmat = sigmat)
  }
  if(plot)
  {
    meansplot <- graph_twoway_assumptions(matrices_obj = matrices_obj)
    list(matrices_obj=matrices_obj, meansplot=meansplot)
  }
  else if(!plot)
  {
    matrices_obj
  }
}

