#' Create input for simulation based two-way factorial experiments
#'
#' The calculate_mean_matrix will generate a matrix of expected mean values for i*j group combinations of a two-way
#' factorial design, as well as a standard deviation matrix for these i*j groups.
#' If the design has repeated measures it will additionally provide correlation and covariance matrices calculated
#' depending on which factors are 'within' factors in the design.
#'
#' The user must provide a reference mean (usually mean in control or untreated group), the expected change for each factor
#' from one level to the next (or from the first to last level) and the number of levels in each factor.
#'
#' Also, if the user specifies factor level combinations which are expected to present interaction and its
#' magnitude with respect to the reference mean, the expected change in the cell means will be incorporated to the
#' aforementioned matrices.
#'
#' We were motivated by sample size calculation for two-way factorial designs with a,b,...,i levels of factor A and
#' a,b,...,j levels of factor B in which the mean outcome value for replicates of cell A=a and B=a are known.
#' Furthermore, there is an expected change in level mean for each of the factors. Finally, interaction can be explicitly
#' introduced to level combinations in which it is expected to occur.
#'
#' @param refmean Numeric - expected mean for first level of both factors A and B
#' @param nlfA Integer - number of levels of factor A
#' @param nlfB integer - number of levels of factor B
#' @param fAeffect Numeric - multiple by which the refmean is modified when going from one level to the next of factor A when endincrement is FALSE (default), or multiple by which the last level of factor A is modified with respect to refmean when endincrement is TRUE
#' @param fBeffect Numeric - multiple by which the refmean is modified when going from one level to the next of factor B when endincrement is FALSE (default), or multiple by which the last level of factor B is modified with respect to refmean when endincrement is TRUE
#' @param groupswinteraction vector length 2 or n*2 matrix - Combination of levels from factors A and B in which interaction is expected
#' @param interact Numeric - value by which the mean from cell or cells indicated in groupswinteraction is multiplied after it has been calculated accordingly to fAeffect and fBeffect
#' @param label_list List length 2 - vectors with the names of the factor levels. The objects in this list should be named as the factors. The use of this option is encouraged as these names are inherited to ANOVA_design.
#' @param sdproportional Logical - whether the standard deviation for each combination of factor levels is a proportion of the respective factor level combination mean, defaults to TRUE
#' @param sdratio Numeric - value by which the expected mean value of a factor level combination is multiplied to obtain the respective standard deviation, defaults to 0.2.
#' @param endincrement Logical - determines if the multiples provided in fAeffect and fBeffect refer to level to level changes (default) or change between first and last levels.
#' @param rho Numeric between -1 and 1 - correlation in outcome within subjects
#' @param whithinf Character - Names the factor with repeated measures. Possible values are NULL, "fA", "fB" or "both"
#' @param plot Logical - Should a line plot with the modeled mean and standard deviations be part of the output. Default=TRUE
#'
#' @return If rho and whithinf are left at their default values of 0 and NULL, respectively, a list with two objects.
#' @return The first consist of two matrices, one of expected means for each cell of the two-way factorial experiment, one of expected standard deviations for said cells.
#' @return If rho is between -1 and 1 but different to 0 and whithinf is either "fA", "fB" or "both", along with the above mentioned output the output will include correlation and covariance matrices.
#'
#' @examples
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' ## Independent design
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean, fAeffect = treateff, fBeffect = timeeff, nlfA = treatgroups, nlfB = timepoints, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
#' ## Inspect plot to check if matrices correspond to design
#' n <- 20
#' independent_experiment <- twoway_simulation_independent(group_size = n, matrices_obj = effects_treat_time)
#' head(independent_experiment, 10)
#'
#' ## Repeated measures design, suppose subjects from 4 independent treatment groups measured at 5 different timepoints.
#' ## We use the same parameters as the independent design example, except we add within factor level correlation and we specify that factor B is the within factor.
#'
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' rho <- 0.8
#' withinf <- "fB"
#'
#' effects_treat_time <- calculate_mean_matrix(refmean = refmean, fAeffect = treateff, fBeffect = timeeff, nlfA = treatgroups,  nlfB = timepoints,
#' rho = rho, withinf = withinf, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
#'
#' ## Inspect plot to check if matrices correspond to design
#' n <- 20
#' repeatedmeasures_experiment <- twoway_simulation_correlated(group_size = n, matrices_obj = effects_treat_time)
#' head(repeatedmeasures_experiment, 10)
#'
#' @export
calculate_mean_matrix <- function(refmean, nlfA, nlfB, fAeffect, fBeffect, groupswinteraction=NULL, interact=1, label_list = NULL,
                                  sdproportional = TRUE, sdratio=0.2, endincrement=FALSE,  rho=0, withinf=NULL, plot=TRUE)
{
  ##Check if labels correspond to group sizes. If no labels are given, assign them.
  if(!is.null(label_list))
  {
    if(!all.equal(as.vector(sapply(label_list, length)), c(nlfA, nlfB)))
    {
      stop("\nNumber of labels must match number of levels for each factor")}
  }
  if(is.null(label_list))
  {
    label_list <- list(fA = letters[1:nlfA], fB = letters[1:nlfB])
  }
  if(rho < -1 | rho>1)
  {
    stop("\nRho must be a number between -1 and 1.")
  }
  if((fAeffect==0|fBeffect==0) & isTRUE(sdproportional))
  {
    if(is.null(withinf))
    {
      warning("\nBy setting '0' effects and proportional SD you will obtain groups with standard deviation of 0.")
    }
    else if(withinf=="fA" | withinf=="fB" | withinf=="both")
    {
      warning("\nBy setting '0' effects and proportional SD you will obtain groups with standard deviation\nof 0 and individuals with 0 covariance.")
    }
  }
  ## Generation of mean matrix
  fAvec <- genvecs(change = fAeffect, reps = nlfA, bystart = endincrement)
  fBvec <- genvecs(change = fBeffect, reps = nlfB, bystart = endincrement)
  effects <- rowSums(expand.grid(fAvec, fBvec-1))
  effmat <- matrix(effects, nlfA, nlfB, dimnames = label_list)
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
  if(is.null(withinf))
  {
    matrices_obj <- list(mean.mat = mean_matrix, sd.mat= sd_matrix)
  } else if(!is.null(withinf))
  {
    if(withinf!="fA" & withinf!="fB" & withinf!="both")
    {
      stop("Possible values for the 'withinf' parameter are 'fA', 'fB' or 'both'")
    }
    ## Generation of correlation and covariance matrices
    partialoutput <- gencovmat(mean_matrix = mean_matrix, sd_matrix = sd_matrix, rho = rho,
                               label_list = label_list, withinf = withinf, nlfA = nlfA, nlfB = nlfB)
    matrices_obj <- list(within.factor = withinf, mean.mat = mean_matrix, sd.mat= sd_matrix, cormat = partialoutput[[1]], sigmat = partialoutput[[2]])
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

