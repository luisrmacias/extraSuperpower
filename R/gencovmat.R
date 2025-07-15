#' Function that generates a covariance matrix taking as input a correlation matrix and a standard deviation matrix or value.
#'
#' May be run independently or internally as part of 'calculate_mean_matrix'.
#'
#' @param correlation_matrix Matrix - Expected correlation between combinations of factor levels
#' @param sd_matrix Numeric or matrix - Standard deviation value or matrix of standard deviation values for combinations of factor levels.
#' @param withinf Character- Factor for which measurements are repeated, options are NULL, "fA", "fB" and "both". If NULL (default) independent measurements will be considered.
#' @param label_list List length 2 - Names of factor levels
#' @param nlfA Integer - number of levels of factor A
#' @param nlfB Integer - number of levels of factor B
#'
#' @return Covariance matrix
#'
#' @examples
#'
#' meanvals <- c(seq(3,9,2),seq(2,8,2),seq(1,7,2))
#' mean_matrix <- matrix(meanvals, 3, 4, byrow = TRUE,
#'                    dimnames = list(A=LETTERS[1:3], B=letters[1:4]))
#'
#' mean_matrix
#' sd_matrix <- mean_matrix*0.2
#'
#' cor_matrix <- gencorrelationmat(mean_matrix = mean_matrix,
#'                                 rho = 0.7, withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' gencovariancemat(cor_matrix, sd_matrix, withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' ##correlation gradient over levels of factor B
#' cor_matrix <- gencorrelationmat(mean_matrix = mean_matrix,
#'                                 rho = c(0.7, 0.4), withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' gencovariancemat(cor_matrix, sd_matrix, withinf = "fB", nlfA = 3, nlfB = 4)
#'
#' @export
gencovariancemat <- function(correlation_matrix, sd_matrix, withinf, label_list=NULL, nlfA, nlfB)
  {

  generic_labels <- list(fA = LETTERS[1:nlfA], fB = letters[1:nlfB])

  if(length(sd_matrix)==1)
  {
    sd_matrix <- matrix(sd_matrix, nrow = nlfA, ncol = nlfB)
    if(is.null(label_list))
    {
      label_list <- generic_labels
    }
  } else if(length(sd_matrix)>1)
  {
    if(nlfA!=nrow(sd_matrix))
    {
      stop(paste("Number of rows in sd_matrix must be equal to number of levels of ", names(label_list)[1]))
    }

    if(nlfB!=ncol(sd_matrix))
    {
      stop(paste("Number of columns in sd_matrix must be equal to number of levels of ", names(label_list)[2]))
    }

    if(is.null(label_list) & identical(dimnames(sd_matrix), generic_labels))
    {
      label_list <-  generic_labels
    } else if (is.null(label_list) & !identical(dimnames(sd_matrix), generic_labels))
      {
        label_list <- dimnames(sd_matrix)
        cnames <- expand.grid(label_list[[2]], label_list[[1]])
        cnames <- paste(cnames$Var2, cnames$Var1, sep = "_")
        rownames(correlation_matrix) <- colnames(correlation_matrix) <- cnames
        message("Covariance matrix names assigned based on names from the standard deviation matrix")
      } else if (!is.null(label_list) & identical(dimnames(sd_matrix), generic_labels) & all.equal(sapply(label_list, length), c(nlfA, nlfB), check.attributes=FALSE))
      {
        warning("The covariance matrix will be generated with user provided names although the standard deviation matrix has generic names")
      } else if (!is.null(label_list) & !identical(dimnames(sd_matrix), label_list))
      {
      stop("Provided label list and standard deviation matrix names do not match")
    }
  }

  cnames <- expand.grid(label_list[[2]], label_list[[1]])
  cnames <- paste(cnames$Var2, cnames$Var1, sep = "_")
  sigmat <- diag(0, prod(nlfA, nlfB))
  rownames(sigmat) <- colnames(sigmat) <- cnames


  if(withinf=="fA")
    {
    rowpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), rownames(correlation_matrix)))
    colpos <- sapply(1:nlfB, function(x) grep(paste0("_", label_list[[2]][x], "$"), colnames(correlation_matrix)))
    if(all.equal(rowpos, colpos))
    {
      for(i in seq(nlfB))
      {
        sigmat[rowpos[,i], colpos[,i]] <- correlation_matrix[rowpos[,i], colpos[,i]]*tcrossprod(sd_matrix[,i])
      }
    }
  }
  if(withinf=="fB")
  {
    rowpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), rownames(correlation_matrix)))
    colpos <- sapply(1:nlfA, function(x) grep(paste0("^", label_list[[1]][x], "_"), colnames(correlation_matrix)))
    if(all.equal(rowpos, colpos))
    {
      for(i in seq(nlfA))
      {
        i <- i-1
        sigmat[(i*nlfB)+(1:nlfB), (i*nlfB)+(1:nlfB)] <- correlation_matrix[(i*nlfB)+(1:nlfB), (i*nlfB)+(1:nlfB)]*tcrossprod(sd_matrix[i+1,])
      }
    }
  }
  if(withinf=="both")
  {
    sigmat <- correlation_matrix*tcrossprod(as.vector(sd_matrix))
  }
  if(!identical(dim(correlation_matrix), dim(sigmat)) & identical(colnames(correlation_matrix), colnames(sigmat)) & identical(rownames(correlation_matrix), rownames(sigmat)))
  {stop("Factors specified for correlation matrix are different from factor or factors specified for covariance matrix")}
  if(!(identical(correlation_matrix>=0, sigmat>=0) & identical(correlation_matrix<=0, sigmat<=0)) & all(sd_matrix>0))
  {stop("The within factors of the correlation and covariance matrices are inconsistent")}
  rhokind <- unique(as.vector(correlation_matrix))
  rhokind <- rhokind[-which(rhokind==0|rhokind==1)]
  if(length(rhokind)>1 | any(rhokind<0))
  {
    sigmat <- Matrix::nearPD(sigmat)$mat
    sigmat <- as.matrix(sigmat)
  }
  sigmat
}
