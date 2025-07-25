% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mean_sd_matrices.R
\name{calculate_mean_matrix}
\alias{calculate_mean_matrix}
\title{Create input for simulation based two-way factorial experiments}
\usage{
calculate_mean_matrix(
  refmean,
  nlfA,
  nlfB,
  fAeffect,
  fBeffect,
  groupswinteraction = NULL,
  interact = 1,
  label_list = NULL,
  sdproportional = TRUE,
  sdratio = 0.2,
  endincrement = TRUE,
  rho = 0,
  withinf = NULL,
  plot = TRUE
)
}
\arguments{
\item{refmean}{Numeric - expected mean for first level of both factors A and B}

\item{nlfA}{Integer - number of levels of factor A}

\item{nlfB}{integer - number of levels of factor B}

\item{fAeffect}{Numeric - multiple that defines how cell means are modified by factor A. With the default \code{endincrement} (\code{TRUE}), determines the last level of factor A with respect to its first level. When \code{endincrement=FALSE} this multiple applies from one level to the next.}

\item{fBeffect}{Numeric - multiple that defines how cell means are modified by factor B. With the default \code{endincrement} (\code{TRUE}), determines the last level of factor B with respect to its first level. When \code{endincrement=FALSE} this multiple applies from one level to the next.}

\item{groupswinteraction}{vector length 2 or n*2 matrix - Combination of levels from factors A and B in which interaction is expected}

\item{interact}{Numeric - value by which the mean from cell or cells indicated in \code{groupswinteraction} is multiplied after it has been calculated accordingly to \code{fAeffect} and \code{fBeffect}}

\item{label_list}{List length 2 - vectors with the names of the factor levels. The objects in this list should be named as the factors. The use of this option is encouraged as these names are used for plotting and inherited to downstream functions.}

\item{sdproportional}{Logical - whether the standard deviation for each combination of factor levels is a proportion of the respective factor level combination mean, defaults to \code{TRUE}}

\item{sdratio}{Numeric - value by which the expected mean value of a factor level combination is multiplied to obtain the respective standard deviation, defaults to 0.2.}

\item{endincrement}{Logical - determines if the multiples provided in \code{fAeffect} and \code{fBeffect} refer to change between first and last levels (default) or level to level changes.}

\item{rho}{Vector length 1 or 2, or 2 by 2 matrix - Controls how the correlation and hence de covariance matrix is built. See 'details' and \code{?gencorrelationmat} examples.}

\item{withinf}{Character - Names the factor with repeated measures. Possible values are NULL, "fA", "fB" or "both"}

\item{plot}{Logical - Should a line plot with the modeled mean and standard deviations be part of the output. Default is \code{TRUE}}
}
\value{
If \code{rho} and \code{withinf} are left at their default values of 0 and NULL, respectively, a cell mean matrix, a cell standard deviation matrix and
optionally a graph that represents both.

If \code{rho} is between -1 and 1 but different to 0 and \code{withinf} is either "fA", "fB" or "both", correlation and covariance matrices are generated
along with the aforementioned output.
}
\description{
This function will generate a matrix of expected mean values for \emph{ab} factor level combinations of a two-way
factorial design by assuming linear effects with possible departure from linearity by interaction. It will also
provide a standard deviation matrix for these \emph{ab} combinations of factor levels. If the design has repeated measures,
it will additionally provide correlation and covariance matrices calculated depending on which factor has repeated
measurements or is the 'within' factor.
}
\details{
The user must provide a reference mean (usually mean in control or untreated group), the expected change for each factor
from first to last level (or from one level to the next) and the number of levels in each factor.

The user can also specify factor level combinations in which interaction is assumed and its magnitude with respect to the
reference mean. The cell mean matrix will be modified accordingly and this can also have an effect of the standard deviation
matrix.

We were motivated by sample size calculation for two-way factorial designs with \emph{1,2,...,a} levels of factor \emph{A} and
\emph{1,2,...,b} levels of factor \emph{B} in which the mean outcome value for replicates of cell \emph{A=1, B=1} are known.
Furthermore, there is an expected change in level mean for each of the factors. Finally, interaction can be explicitly
introduced to level combinations in which it is expected to occur.

If a repeated measures experiment is intended \code{withinf} must be set to "fA", "fB" or "both", depending on which is the 'within' factor.
If \code{rho} is a vector length 1, the within subject correlation will be constant for the factor defined in \code{withinf}. If \code{rho} is a vector
length 2 and \code{withinf} is either "fA" or "fB" a correlation gradient will be created from the first to second value of \code{rho}. If \code{rho} is
a vector length 2 and \code{withinf="both"}, the first element of \code{rho} will be the correlation within factor A, while the second element will
be the correlation within factor B. If \code{rho} is a 2*2 matrix, only possible if \code{withinf="both"}, a correlation gradient will be created
across rows of \code{rho} for each of the factors.
}
\examples{
refmean <- 1
treatgroups <- 4
timepoints <- 5
treateff <- 1.5
timeeff <- 0.85
factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)
## Independent design
effects_treat_time <- calculate_mean_matrix(refmean = refmean,
                                            fAeffect = treateff, fBeffect = timeeff,
                                            nlfA = treatgroups, nlfB = timepoints,
                                            label_list = factors_levels_names)
## Inspect plot to check if matrices correspond to design
effects_treat_time$meansplot
n <- 20
independent_experiment <- twoway_simulation_independent(group_size = n,
                                      matrices_obj = effects_treat_time)

head(independent_experiment, 10)

## Repeated measures design, suppose subjects from 4 independent treatment groups measured
## at 5 different timepoints.
## We use the same parameters as the independent design example, except we add within factor level
## correlation and we specify that factor B is the within factor.


refmean <- 1
treatgroups <- 4
timepoints <- 5
treateff <- 1.5
timeeff <- 0.85
rho <- 0.8
withinf <- "fB"

factors_levels_names <- list(treatment=letters[1:treatgroups], time=1:timepoints)

effects_treat_time <- calculate_mean_matrix(refmean = refmean, fAeffect = treateff,
                      fBeffect = timeeff, nlfA = treatgroups,  nlfB = timepoints,
                      rho = rho, withinf = withinf, label_list = factors_levels_names)

## Plot should look the same, structure within data can be checked once simulated
effects_treat_time$meansplot

n <- 20
repeatedmeasures_experiment <- twoway_simulation_correlated(group_size = n,
                                          matrices_obj = effects_treat_time)
head(repeatedmeasures_experiment, 10)

}
