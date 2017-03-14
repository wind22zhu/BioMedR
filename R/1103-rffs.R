#' Random Forest Cross-Valdidation for feature selection 
#' 
#' Random Forest Cross-Valdidation for feature selection 
#' 
#' This function shows the cross-validated prediction performance of models with 
#' sequentially reduced number of predictors (ranked by variable importance) via 
#' a nested cross-validation procedure.
#' 
#' @param trainx matrix or data frame containing columns of predictor variables
#' 
#' @param trainy vector of response, must have length equal to the number of rows in 
#'            \code{trainx}
#'            
#' @param cv.fold The fold, the defalut is 5.
#' 
#' @param scale If \code{"log"}, reduce a fixed proportion (\code{step}) of variables 
#'              at each step, otherwise reduce \code{step} variables at a time
#' 
#' @param step If \code{log=TRUE}, the fraction of variables to remove at each step, 
#'             else remove this many variables at a time
#' 
#' @param mtry A function of number of remaining predictor variables to use as the 
#'             \code{mtry} parameter in the \code{randomForest} call
#'              
#' @param recursive Whether variable importance is (re-)assessed at each step of variable 
#'                  reduction                             
#'
#' @return  A list with the following three components::
#' \itemize{
#' \item \code{n.var} - vector of number of variables used at each step
#' \item \code{error.cv} - corresponding vector of error rates or MSEs at each step
#' \item \code{res} - list of n.var components, each containing the feature importance values from 
#' the cross-validation
#' }

#' @keywords  rf.fs
#'
#' @aliases rf.fs
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{rf.cv}} for the Cross-Validation of Classification and 
#' Regression models using Random Forest 
#' 
#' @export rf.fs
#' 
#' @importFrom randomForest randomForest
#' 
#' @references
#' Svetnik, V., Liaw, A., Tong, C. and Wang, T., Application of Breiman's Random 
#' Forest to Modeling Structure-Activity Relationships of Pharmaceutical Molecules, 
#' MCS 2004, Roli, F. and Windeatt, T. (Eds.) pp. 334-343.
#' 
#' @examples
#' training = read.csv(system.file('sysdata/training1.csv', package = 'BioMedR'), header = TRUE)
#' y = training[, 1]
#' x = training[, -1]
#' result = rf.fs(x, y)
#'
rf.fs <- function(trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,
                 mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE) {
    classRF <- is.factor(trainy)
    n <- nrow(trainx)
    p <- ncol(trainx)
    if (scale == "log") {
        k <- floor(log(p, base=1/step))
        n.var <- round(p * step^(0:(k-1)))
        same <- diff(n.var) == 0
        if (any(same)) n.var <- n.var[-which(same)]
        if (! 1 %in% n.var) n.var <- c(n.var, 1)
    } else {
        n.var <- seq(from=p, to=1, by=step)
    }
    k <- length(n.var)
    cv.pred <- vector(k, mode="list")
    imp.var.num <- vector(k, mode="list")
    imp.var.idx <- vector(k, mode="list")
    for (i in 1:k) cv.pred[[i]] <- trainy
    for (i in 1:k) imp.var.num[[i]] <- matrix(NA, cv.fold, n.var[i])
    for (i in 1:k) imp.var.idx[[i]] <- matrix(NA, cv.fold, n.var[i])
    ## Generate the indices of the splits
    ## Stratify the classes for classification problem.
    ## For regression, bin the response into 5 bins and stratify.
    if(classRF) {
        f <- trainy
    } else {
        ##f <- cut(trainy, c(-Inf, quantile(trainy, 1:4/5), Inf))
		f <- factor(rep(1:5, length=length(trainy))[order(order(trainy))])
    }
    nlvl <- table(f)
    idx <- numeric(n)
    for (i in 1:length(nlvl)) {
        idx[which(f == levels(f)[i])] <-  sample(rep(1:cv.fold, length=nlvl[i]))
    }

    for (i in 1:cv.fold) {
        ## cat(".")
        all.rf <- randomForest(trainx[idx != i, , drop=FALSE],
                               trainy[idx != i],
                               trainx[idx == i, , drop=FALSE],
                               trainy[idx == i],
                               mtry=mtry(p), importance=TRUE)
        cv.pred[[1]][idx == i] <- all.rf$test$predicted
        impvar <- (1:p)[order(all.rf$importance[,1], decreasing=TRUE)]
        imp.var.num[[1]][i,] <- as.numeric(all.rf$importance[,1])
        imp.var.idx[[1]][i,] <- impvar
        impidx <- impvar
        for (j in 2:k) {
            imp.idx <- impvar[1:n.var[j]]
            sub.rf <-
                randomForest(trainx[idx != i, imp.idx, drop=FALSE],
                             trainy[idx != i],
                             trainx[idx == i, imp.idx, drop=FALSE],
                             trainy[idx == i],
                             mtry=mtry(n.var[j]), importance=recursive)
            cv.pred[[j]][idx == i] <- sub.rf$test$predicted
            imp.var.num[[j]][i,] <- as.numeric(sub.rf$importance[,1])
            
            ## For recursive selection, use importance measures from the sub-model.
            
            imp <-(1:length(imp.idx))[order(sub.rf$importance[,1], decreasing=TRUE)]
              
            impvar.idx <- c()
            for (l in 1L:n.var[j]){
              impvar.idx[l] = impidx[imp[l]]
            } 
            impidx <- impvar.idx 
            if (recursive) {
              impvar = impidx
            }
            imp.var.idx[[j]][i,] <- impvar.idx
            
        }
    NULL
    }
  
    imp.var = vector(k, mode="list")
    for (i in 1:k) imp.var[[i]] <- matrix(0, cv.fold, p)
    imp.var[[1]] = imp.var.num[[1]]
    for (i in 2L:k) {
        for (j in 1L:cv.fold) {
            for (l in 1L:n.var[i]){
                imp.var[[i]][j, imp.var.idx[[i]][j,l]] = imp.var.num[[i]][j,l]
            }
        }
    }
    ## cat("\n")
    res <- vector(k, mode="list")
    for(i in 1:k) {
        res[[i]] <- apply(imp.var[[i]], 2L, sum)
    }
  
  
    if(classRF) {
        error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
    } else {
        error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
    }
    names(error.cv) <- names(cv.pred) <- n.var
    list(n.var = n.var, error.cv = error.cv, res = res)
}
