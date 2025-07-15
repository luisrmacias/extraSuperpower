#' Function that generates a correlation matrix taking as input number of factors for each level, factor or factors that present correlation and rho value or values.
#' Additionally, a mean matrix is required to check consistency.
#'
#' May be run independently or internally as part of `calculate_mean_matrix`.
#'
#' @param mean_matrix Matrix - cell mean value matrix
#' @param rho Vector length 1 or 2, or 2 by 2 matrix - Controls how the correlation and hence de covariance matrix is built. See details.
#' @param withinf Character- Factor for which measurements are repeated, options are NULL, "fA", "fB" and "both". If NULL (default) independent measurements will be considered.
#' @param label_list List length 2 - Names of factor levels
#' @param nlfA Integer - number of levels of factor A
#' @param nlfB Integer - number of levels of factor B
#'
#' @details
#' For a repeated measures experiment `withinf` must be set to "fA", "fB" or "both", depending on which is the 'within' factor.
#' If `rho` is a vector length 1, the within subject correlation will be constant for the factor defined in `withinf`. If `rho` is a vector
#' length 2 and `withinf` is either "fA" or "fB" a correlation gradient will be created from the first to second value of `rho`. If `rho` is
#' a vector length 2 and `withinf="both"`, the first element of `rho` will be the correlation within factor A, while the second element will
#' be the correlation within factor B. If `rho` is a 2*2 matrix, only possible if `withinf="both"`, a correlation gradient will be created
#' across rows of `rho` for each of the factors.
#'
#' @importFrom methods is
#'
#' @return Correlation matrix
#'
#' @examples
#'
#' meanvals <- c(seq(3,9,2),seq(2,8,2),seq(1,7,2))
#' mean_matrix <- matrix(meanvals, 3, 4, byrow = TRUE,
#'                    dimnames = list(A=LETTERS[1:3], B=letters[1:4]))
#'
#' mean_matrix
#'
#' gencorrelationmat(mean_matrix = mean_matrix, rho = 0.7, withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' ##correlation gradient over levels of factor B
#' gencorrelationmat(mean_matrix = mean_matrix, rho = c(0.7, 0.4), withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' ##gradient both factors
#'
#' rhovals <- matrix(c(0.7, 0.4), 2, 2, byrow = TRUE)
#' gencorrelationmat(mean_matrix = mean_matrix, rho = rhovals, withinf = "both", nlfA = 3, nlfB = 4)
#'
#' @export
gencorrelationmat <- function(mean_matrix, rho, label_list=NULL, withinf, nlfA, nlfB)
  {
  if(nlfA!=nrow(mean_matrix))
  {
    stop(paste0("Number of rows in 'mean_matrix' must be equal to number of levels of '", names(label_list)[1], "'"))
  }

  if(nlfB!=ncol(mean_matrix))
  {
    stop(paste0("Number of columns in 'mean_matrix' must be equal to number of levels of '", names(label_list)[2], "'"))
  }

  if((length(rho)>2 & is.vector(rho)) | (is.matrix(rho) & length(rho)!=4))
  {
    stop("'rho' must be a single value, a vector length 2 or a 2 by 2 matrix")
  }

  if(is(rho, "matrix") & withinf!="both")
  {
    stop("'rho' can only be a matrix if 'withinf' is set to 'both'")
  }
  generic_labels <- list(fA = LETTERS[1:nlfA], fB = letters[1:nlfB])

    if(is.null(label_list) & identical(dimnames(mean_matrix), generic_labels))
  {
    label_list <-  generic_labels
  } else if (is.null(label_list) & !identical(dimnames(mean_matrix), generic_labels))
  {
    label_list <- dimnames(mean_matrix)
    message("Correlation matrix names assigned based on names from the mean matrix")
  } else if (!is.null(label_list) & identical(dimnames(mean_matrix), generic_labels) & all.equal(sapply(label_list, length), c(nlfA, nlfB), check.attributes=FALSE))
  {
    warning("The correlation matrix will be generated with user provided names although the mean matrix has generic names")
  } else if (!is.null(label_list) & !identical(dimnames(mean_matrix), label_list))
  {
    stop("Provided label list and mean matrix names do not match")
  }

  cnames <- expand.grid(label_list[[2]], label_list[[1]])
  cnames <- paste(cnames$Var2, cnames$Var1, sep = "_")

  cormat <- diag(1, prod(nlfA, nlfB))
  rownames(cormat) <- colnames(cormat) <- cnames

  if(withinf=="fA")
  {
    if(length(rho)==1)
    {
      tmpmat <- matrix(rho, nlfA, nlfA)
      diag(tmpmat) <- 1
    } else if (length(rho)==2)
    {
      rho <- seq(rho[1], rho[2], length.out=nlfA-1)
      tmpmat <- diag(nlfA)
      for(i in seq(nlfA-1))
      {
        j <- i-1
        tmpmat[-(1:i),i] <- rho[1:(length(rho)-j)]
      }
      tmpmat[upper.tri(tmpmat)] <- t(tmpmat)[upper.tri(tmpmat)]
    }
    rowpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), rownames(cormat)))
    colpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), colnames(cormat)))
    if(all.equal(rowpos, colpos))
    {
      for(i in seq(nlfB))
      {
        cormat[rowpos[,i], colpos[,i]] <- tmpmat
      }
    }
  }
  if(withinf=="fB")
  {
    if(length(rho)==1)
    {
      tmpmat <- matrix(rho, nlfB, nlfB)
      diag(tmpmat) <- 1
    } else if (length(rho)==2)
    {
      rho <- seq(rho[1], rho[2], length.out=nlfB-1)
      tmpmat <- diag(nlfB)
      for(i in seq(nlfB-1))
      {
        j <- i-1
        tmpmat[-(1:i),i] <- rho[1:(length(rho)-j)]
      }
      tmpmat[upper.tri(tmpmat)] <- t(tmpmat)[upper.tri(tmpmat)]
    }
    rowpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), rownames(cormat)))
    colpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), colnames(cormat)))
    if(all.equal(rowpos, colpos))
    {
      for(i in seq(nlfA))
      {
        cormat[rowpos[,i], colpos[,i]] <- tmpmat
      }
    }
  }
  if(withinf=="both")
  {
    if(length(rho)==1)
    {
      cormat[cormat!=1] <- rho
    } else if(length(rho)>1)
    {
    if (length(rho)==2)
    {
      tmpmatA <- matrix(rho[1], nlfA, nlfA)
      diag(tmpmatA) <- 1
      tmpmatB <- matrix(rho[2], nlfB, nlfB)
      diag(tmpmatB) <- 1
    } else if (is.matrix(rho) & length(rho)==4)
    {
      rhoA <- seq(rho[1,1], rho[1,2], length.out=nlfA-1)
      tmpmatA <- diag(nlfA)
      rhoB <- seq(rho[2,1], rho[2,2], length.out=nlfB-1)
      tmpmatB <- diag(nlfB)
      for(i in seq(nlfA-1))
      {
        j <- i-1
        tmpmatA[-(1:i),i] <- rhoA[1:(length(rhoA)-j)]
      }
      tmpmatA[upper.tri(tmpmatA)] <- rev(tmpmatA[lower.tri(tmpmatA)])
      for(i in seq(nlfB-1))
      {
        j <- i-1
        tmpmatB[-(1:i),i] <- rhoB[1:(length(rhoB)-j)]
      }
      tmpmatB[upper.tri(tmpmatB)] <- rev(tmpmatB[lower.tri(tmpmatB)])
    }
      rowpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), rownames(cormat)))
      colpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), colnames(cormat)))
      if(all.equal(rowpos, colpos))
      {
        for(x in seq(nlfB))
        {
          cormat[rowpos[,x], colpos[,x]] <- tmpmatA
        }
      }

      rowpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), rownames(cormat)))
      colpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), colnames(cormat)))
      if(all.equal(rowpos, colpos))
      {
        for(i in seq(nlfA))
        {
          cormat[rowpos[,i], colpos[,i]] <- tmpmatB
        }
      }
    }

  }
  cormat
}

