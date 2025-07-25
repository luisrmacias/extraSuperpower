#' Simulate measurements repeated over either or both factors of a two-way design
#'
#' Both regular and internal function.
#' As regular function takes input generated by the `calculate_mean_matrix` function and iteratively simulates repeated measures two-way factorial experiments.
#' Data are sampled from a normal, skewed normal or truncated normal distribution.
#'
#' As internal function runs with a single iteration inside `graph_twoway_assumptions`, which in itself is inside  'calculate_mean_matrix' to generate data for the cell mean and standard deviation plot.
#'
#' @param group_size Integer or matrix - Sample size for each group (combination of factor levels). If `balanced=TRUE` (default) `group_size` must be an integer. If `balanced=FALSE` `group_size` must be a matrix.
#' @param matrices_obj List - Output generated by `calculate_mean_matrix` that include cell mean and covariance matrices
#' @param distribution Character - Type of distribution from which to sample, possible values are "normal", "skewed" and "truncated"
#' @param shape Vector - Degree of skewness in the distribution. May be a single value, have a length equal to the number of levels of any one of the factors or a length equal to the product of the length of each factor.
#' @param inferior_limit Numeric - Value for which the distribution is truncated on the left. Only valid if `distribution="truncated.normal"`
#' @param superior_limit Numeric - Value for which the distribution is truncated on the right. Only valid if `distribution="truncated.normal"`
#' @param balanced Logical - Whether the study will be performed with the same number of subjects in all groups. Default is `TRUE`. See 'details'.
#' @param loss Character - Type of selection of subjects in groups that have less observations than `max(group_size)`. Possible values are 'random' and 'sequential'. Ignored if `balanced=TRUE`. See 'details'.
#' @param nsims Integer - Number of iterations
#'
#' @details
#' For unbalanced repeated measures designs, this function generates a simulation with `max(group_size)` for all combinations of factors and then eliminates observations.
#' If `loss="random"` elimination of in those factor combinations that have less participants or study subjects will occur at random. If `loss="sequential"` the participants or subjects
#' from the groups with less observations will be a subset of participants or subjects of groups with more observations. This may not sound like the most efficient way to proceed, is quite fast anyhow.
#'
#' The 'n' column in the output will reflect how many observations each factor combination has. This should match the input matrix.
#'
#' @return Dataframe with simulated outcome values, factor level labels and iteration number.
#'
#' @examples
#' ## Repeated measures design, suppose subjects from 4 independent treatment groups
#' ## measured at 5 different timepoints.
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
#' ## Inspect plot to check if matrices correspond to design
#' effects_treat_time$meansplot
#'
#' n <- 20
#' repeatedmeasures_experiment <- twoway_simulation_correlated(group_size = n,
#'                                 matrices_obj = effects_treat_time)
#'
#' head(repeatedmeasures_experiment, 10)
#'
#' @export
twoway_simulation_correlated <- function(group_size, matrices_obj, distribution="normal", shape=0, inferior_limit= -Inf, superior_limit=Inf, balanced=TRUE, loss=NULL, nsims=200)
{
  if(!is.character(matrices_obj[[1]][[1]]))
  {
    stop("Your 'matrices_obj' does not specify which factor is within participants.\nIf you have a repeated measures design you need to specify the within factor and its correlation when using 'calculate_mean_matrix'.")
  }
  if(!all(sapply(matrices_obj[-1], is.matrix)) & !is.numeric(matrices_obj[[2]][1]))
  {
    matrices_obj <- matrices_obj$matrices_obj
  }
  label_list <- dimnames(matrices_obj$mean.mat)
  factor_levels <- dim(matrices_obj$mean.mat)
  nlevels <- prod(factor_levels)
  mean_matrix <- as.vector(t(matrices_obj$mean.mat))
  if (!balanced & length(group_size)==1)
  {
    stop("Are you sure you want to set balanced to false?")
  }
  if (balanced & length(group_size)>1)
  {
    stop("If you wish a balanced design a single integer is required as 'group_size' argument.")
  }
  if(balanced & length(group_size==1))
  {
    if(as.integer(group_size)!=group_size)
    {
      stop("For balanced designs 'group_size' must be a single integer")
    }
  }

  if((distribution=="normal" | distribution=="skewed") & (superior_limit<Inf | inferior_limit>-Inf))
  {
    warning(paste("Superior and inferior limits are ignored when distribution is", distribution))
  }
  if((distribution=="normal" | distribution=="truncated.normal") & (any(shape>0)))
  {
    warning(paste("Skewness is ignored when distribution is", distribution))
  }

  sigmatrix <- matrices_obj$sigmat
  withinf <- matrices_obj$within.factor
  sampn <- max(group_size)
  #generating factor data frame
  if(withinf=="both")
  {
    subject <- rep(1:sampn, nlevels)
  }
  else if (withinf=="fA")
  {
    subject <- rep(1:(sampn*factor_levels[2]), factor_levels[1])
  }
  else if (withinf=="fB")
  {
    subject <- as.vector(sapply(1:factor_levels[1], function(x)
    {y <- x-1; rep((1+(sampn*y)):(sampn*x), factor_levels[2])}))
  }
  fdata <- as.data.frame(MASS::mvrnorm(sampn, mean_matrix, sigmatrix))
  fdata$subject <- 1:sampn
  fdata <- reshape2::melt(fdata, id.vars = "subject", variable.name = "cond",
                          value.name = "y")
  fdata$subject <- factor(subject)
  for (j in 1:2)
  {
    assigned_factor_levels <- rep(as.list(paste(names(label_list)[j], label_list[[j]], sep = "_")),
                                  each = sampn * nlevels/prod(factor_levels[1:j]),
                                  times = nlevels/prod(factor_levels[j:2]))
    assigned_factor_levels <- unlist(assigned_factor_levels)
    assigned_factor_levels <- factor(assigned_factor_levels, levels = unique(assigned_factor_levels))
    fdata <- cbind(fdata, assigned_factor_levels)
  }
  names(fdata)[4:5] <- names(label_list)
  ##generating outcomes
  if(distribution=="normal")
  {
    y <- suppressMessages(replicate(nsims, reshape2::melt(as.data.frame(MASS::mvrnorm(n = sampn, mu = mean_matrix, Sigma = sigmatrix)))$value))
  }
  else if (distribution=="skewed")
  {
    if(!is.numeric(shape))
    {stop("The shape must be numeric")}
    shapelen <- length(shape)
    if(shapelen>1 & (shapelen!=factor_levels[1] & shapelen!=factor_levels[2] & shapelen!=nlevels))
    {stop("Length of shape should be either 1, the number of levels of 1 of the factors of the product of number of levels of both factors")}
    if(shapelen==1)
    {
      gam <- rep(shape, nlevels)
    } else if (shapelen==factor_levels[1])
    {
      gam <- rep(shape, factor_levels[2])
    } else if (shapelen==factor_levels[2])
    {
      gam <- rep(shape, factor_levels[1])
    } else if (shapelen==nlevels)
    {
      gam <- shape
    }
    # include parameter shrinkage

    cp <- list(mean=mean_matrix, var.cov=sigmatrix, gamma1=gam)
    dp <- try(sn::cp2dp(cp, "SN"), silent = TRUE)
    while(class(dp)=="try-error")
    {
      gam <- gam/3
      cp <- list(mean=mean_matrix, var.cov=sigmatrix, gamma1=gam)
      dp <- try(sn::cp2dp(cp, "SN"), silent = TRUE)
    }
    # if (class(dp)=="try-error")
    # {stop("Please use a smaller shape value or values. Shape is restricted by correlation.")}
    # y <- suppressMessages(replicate(nsims, reshape2::melt(as.data.frame(CensMFM::rMSN(n = sampn, mu=as.vector(mean_matrix), Sigma = sigmatrix, shape = gam)))$value))
    y <- suppressMessages(replicate(nsims, reshape2::melt(as.data.frame(sn::rmsn(n = sampn, dp=dp)))$value))
    # op <- list(xi=refs$mean.mat, Psi=refs$sigmat, lambda=gam)
    # dp <- sn::op2dp(op, "SN")
    # y <- suppressMessages(replicate(nsims, reshape2::melt(as.data.frame(sn::rmsn(n = sampn, dp=dp)))$value))
  }
  else if (distribution=="truncated.normal")
  {
    if(inferior_limit>min(mean_matrix)*1.2 | superior_limit<max(mean_matrix)*0.85)
    {
      warning("The lower or upper bound for the truncated distribution is too extreme for efficient sampling.")
    }
    y <- suppressWarnings(replicate(nsims, reshape2::melt(tmvtnorm::rtmvnorm2(sampn, mean=mean_matrix, sigma = sigmatrix, lower = rep(inferior_limit, nlevels), upper = rep(superior_limit, nlevels),  algorithm ="gibbs"))$value))
  }

  sim <- lapply(seq(nsims),
                function(x)
                {
                  fdata$y <- y[,x]
                  fdata$iteration <- x
                  fdata
                })
  sim <- do.call(rbind, sim)
  if(!balanced & length(group_size)>1)
  {
    if(length(loss)==0)
    {
      stop("Value of 'loss' must be either \"random\" or \"sequential\" in unbalaced repeated measures designs.")
    }
    tosample <- sim[sim$iteration==1,]
    if(loss=="random")
    {
      keep <- sapply(1:length(levels(tosample$cond)),
                     function(x) sample(tosample$subject[tosample$cond==levels(tosample$cond)[x]], t(group_size)[x]))
      names(keep) <- levels(tosample$cond)
      sim <- lapply(levels(sim$cond), function(x)
      {
        selection <- sim$subject[sim$cond==x] %in% keep[[grep(x, names(keep))]]
        sim[sim$cond==x & selection,]
      })
    } else if(loss=="sequential")
    {
      if(withinf=="fA")
      {
        withcol <- 4
        betwcol <- 5
        group_size <- t(group_size)
      }
      if(withinf=="fB")
      {
        withcol <- 5
        betwcol <- 4
      }
      if(withinf=="fA" | withinf=="fB")
      {
        initial <- lapply(1:length(levels(tosample[,betwcol])),
               function(x)(sample(unique(tosample$subject[tosample[,betwcol]==levels(tosample[,betwcol])[x]]),
                                group_size[x,1])))
        gather <- NULL
        for(i in 1:length(initial))
        {
          ord <- order(group_size[i,], decreasing = TRUE)
          avail <- initial[[i]]
          j <- sample(avail, group_size[i,2])
          keep <- list(avail, j)
          avail <- j

          for(j in ord[-(1:2)])
          {
            j <- sample(avail, group_size[i,j])
            avail <- j
            keep <- rlist::list.append(keep, j)
          }
          names(keep) <- levels(tosample[,withcol])[ord]
          gather <- c(gather, keep)
        }
      } else if (withinf=="both")
        {
        topn <- max(group_size)
        nvec <- topn - as.vector(t(group_size))
        initial <- sample(1:topn, topn-unique(nvec)[1])
        dup <- initial
        for(i in seq(length(unique(nvec))))
        {
          tmp <- sample(dup, max(group_size)-unique(nvec)[i])
          keep <- c(dup, tmp)
          dup <- tmp
        }
        keep <- c(initial, keep)
        keep <- order(table(keep), decreasing = TRUE)
        gather <- sapply(t(group_size), function(x) keep[1:x])
        }
      if(withinf=="fA")
      {
        namesec <- order(sapply(strsplit(levels(tosample$cond), "_"), "[", 2))
        names(gather) <- levels(tosample$cond)[namesec]
      } else if(withinf=="fB" | withinf=="both")
      {
        names(gather) <- levels(tosample$cond)
      }
      sim <- lapply(levels(sim$cond), function(x)
      {
        selection <- sim$subject[sim$cond==x] %in% gather[[grep(x, names(gather))]]
        sim[sim$cond==x & selection,]
      })
    }
    sim <- do.call(rbind, sim)
    sim$subject <- droplevels(sim$subject)
    sim$n <- rep(unlist(mapply(rep, t(group_size), t(group_size))), nsims)
  }
  if(balanced & length(group_size)==1)
  {
    sim$n <- group_size
  }
  list(simulated_data = sim, withinf = withinf)
}
