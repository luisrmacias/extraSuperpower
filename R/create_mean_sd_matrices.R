#' Create input for simulation based two-way factorial experiments
#'
#' The mean_sd_matrix function is intended as a companion function to the Superpower package ANOVA_design function,
#' followed by the Superpower ANOVA_power function.
#' We were motivated by sample size calculation for two-way factorial designs with a,b,...,i levels of factor A and
#' a,b,...,j levels of factor B in which the mean outcome value for replicates of cell A=a and B=a are known.
#' Furthermore, there is an estimate of the change in level mean for each of the factors. Finally, some level
#' combinations are expected to present interaction.
#'
#' With user provided reference mean (mean for A=a and B=a), the expected change for each factor from one level to
#' the next (or from the first to last level) and the number of levels in each factor, the function creates a
#' matrix of expected means and a matrix of expected standard deviations with a cell for each combination of levels
#' i and j.
#'
#' Aditionally, if the user specifies factor level combinations which are expected to present interaction and its
#' magnitude with respect to the reference mean, the expected change in the cell means will be incorporated to the
#' aforementioned matrices.
#'
#' @param refmean Expected mean for lowest level of both factors A and B
#' @param nlfA Number of levels of factor A
#' @param nlfB Number of levels of factor B
#' @param fAeffect Multiple by which the refmean is modified when going from one level to the next of factor A when endincrement is FALSE (default), or multiple by which the last level of factor A is modified with respect to refmean when endincrement is TRUE
#' @param fBeffect Multiple by which the refmean is modified when going from one level to the next of factor B when endincrement is FALSE (default), or multiple by which the last level of factor B is modified with respect to refmean when endincrement is TRUE
#' @param groupswinteraction Combination of levels from factors A and B in which interaction is expected
#' @param interact Value by which the mean from cell or cells indicated in groupswinteraction is multiplied after it has been calculated accordingly to fAeffect and fBeffect
#' @param label_list List that contains vectors with the names of the factor levels. The objects in this list should be named as the factors. The use of this option is encouraged as these names are inherited to ANOVA_design.
#' @param sdproportional Whether the standard deviation for each combination of factor levels is a proportion of the respective factor level combination mean, defaults to TRUE
#' @param sdratio Value by which the expected mean value of a factor level combination is multiplied to obtain the respective standard deviation, defaults to 0.2.
#' @param endincrement Determines if the multiples provided in fAeffect and fBeffect refer to level to level changes (default) or change between first and last levels.
#' @examples
#' refmean <- 1
#' treatgroups <- 4
#' timepoints <- 5
#' treateff <- 1.5
#' timeeff <- 0.85
#' effects_treat_time <- mean_sd_matrix(refmean = refmean, fAeffect = treateff, fBeffect = timeeff, nlfA = treatgroups, nlfB = timepoints, label_list = list(treatment=letters[1:treatgroups], time=1:timepoints))
#' muvec <- as.vector(effects_treat_time$mean.mat)
#' sdvec <- as.vector(effects_treat_time$sd.mat)
#' bwdes <- Superpower::ANOVA_design("5w*4b", n=15, mu =  muvec, sd = sdvec, label_list = rev(dimnames(effects_treat_time$mean.mat)))
#' @export
mean_sd_matrix <- function(refmean, nlfA, nlfB, fAeffect, fBeffect, groupswinteraction=NULL, interact=1, label_list = NULL, sdproportional = TRUE, sdratio=0.2, endincrement=FALSE)
  {
  if(!is.null(label_list))
  {
  if(!all.equal(as.vector(sapply(label_list, length)), c(nlfA, nlfB)))
  {
    stop("Number of labels must match number of levels for each factor")}
  }
  if(is.null(label_list))
  {
    label_list <- list(f1 = letters[1:nlfA], f2 = letters[1:nlfB])
  }
  f1vec <- genvecs(change = fAeffect, reps = nlfA, bystart = endincrement)
  f2vec <- genvecs(change = fAeffect, reps = nlfB, bystart = endincrement)
  effects <- rowSums(expand.grid(f1vec, f2vec-1))
  effmat <- matrix(effects, nlfA, nlfB, dimnames = label_list)
  if(interact!=1)
  {
    if (!(is.matrix(groupswinteraction) | length(groupswinteraction)==2))
    {
      stop("Please provide a vector length 2 or a n*2 matrix")
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
  if (sdproportional)
  {
    sd_matrix <- abs(mean_matrix*sdratio)
  }else if (!sdproportional)
  {
    sd_matrix <- abs(mean(mean_matrix)*sdratio)
  }
  list(mean.mat = mean_matrix, sd.mat= sd_matrix)
}

