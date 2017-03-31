.result = function(Pred,cl){
  res = table(Pred,cl)
  TP = as.double(res[2, 2])   # maybe not [2,2] cause 2 means "malignant"--positive
  TN = as.double(res[1, 1])  # maybe not [1,1] cause 1 means "benign"--negative
  FP = as.double(res[2, 1])
  FN = as.double(res[1, 2])
  
  ACC = (TP + TN) / (TP + TN + FP + FN)
  SE = TP / (TP + FN)
  SP = TN/(FP + TN)
  F1 = 2 * TP / (2 * TP + TP + FN) # F1 score is the harmonic mean of precision and sensitivity
  MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN+FN))
  result = list()
  result$pred_label = Pred
  result$table = res
  result$ACC = ACC
  result$SE = SE
  result$SP = SP
  result$F1 = F1
  result$MCC = MCC
  result
}

#' The Cross-Validation of Classification and Regression models using Partial Least Squares  
#' 
#' The Cross-Validation of Classification and Regression models using Partial Least Squares  
#' 
#' This function performs k-fold cross validation for partial least squares regression and 
#' classification.
#' 
#' @param xtr A data frame or a matrix of predictors.
#' 
#' @param ytr A response vector. If a factor, classification is assumed, otherwise 
#'            regression is assumed.
#'            
#' @param cv.fold The fold, the defalut is 5.
#' 
#' @param maxcomp Maximum number of components included within the models, 
#'                if not specified, default is the variable (column) numbers in x.                                                        
#'
#' @return the retrun a list containing four components:
#' \itemize{
#' \item \code{plspred} - the predicted values of the input data based on cross-validation
#' \item \code{Error} - error for all samples
#' \item \code{RMSECV} - Root Mean Square Error for cross-validation
#' \item \code{Q2} - R2 for cross-validation
#' }
#' 
#' @keywords  pls.cv
#'
#' @aliases pls.cv
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{rf.cv}} for the Cross-Validation of Classification and 
#' Regression models using Random Forest 
#' 
#' @importFrom pls plsr
#' 
#' @export pls.cv
#' 
#' @examples
#' training = read.csv(system.file('sysdata/training2.csv', package = 'BioMedR'), header = TRUE)
#' y = training[, 1]
#' x = training[, -1]
#' pls.tr <- pls.cv(x, y)
#' 

pls.cv <- function(xtr, ytr,
                   cv.fold = 5, maxcomp = NULL) {
  
  mx = dim(xtr)[1]
  plspred <- ytr
  index = rep(1:cv.fold, nrow(xtr))
  ind = index[1:nrow(xtr)]
  for (k in 1:cv.fold) {
    cat(".")
    xcal <- xtr[ind != k, ] 
    ycal <- ytr[ind != k]
    xtest <- xtr[ind == k, ] 
    ytest <- ytr[ind == k]  
      xcal <- scale(xcal, center = TRUE, scale = TRUE)
      xcen <- attributes(xcal)$'scaled:center'
      xsca <- attributes(xcal)$'scaled:scale'
      xtest <- scale(xtest, xcen, xsca)
  
      ycal <- scale(ycal, center = TRUE, scale = FALSE)
      ycen <- attributes(ycal)$'scaled:center'
    
    mvrout <- pls::plsr(ycal ~ ., ncomp = min(c(maxcomp, dim(xcal))), 
                 data = data.frame(xcal, ycal), scale = FALSE,
                 method = 'simpls')
      plspred[ind == k]<- predict(mvrout, comps = 1:min(c(maxcomp, dim(xcal))), 
                                  xtest) + ycen     
  }
    cat("\n")
    RMSECV = sqrt(t(ytr - plspred) %*% (ytr - plspred) / mx)
    q2 = 1 - t(ytr - plspred) %*% (ytr - plspred) / (t(ytr - mean(ytr)) %*% (ytr - mean(ytr))) 
    err = ytr - plspred
    ret <- list(plspred = plspred, Error = err, RMSECV = RMSECV, Q2 = q2)
    return(ret)
}
