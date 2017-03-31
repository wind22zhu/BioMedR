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

#' The Cross-Validation of Classification and Regression models using Random Forest 
#' 
#' The Cross-Validation of Classification and Regression models using Random Forest 
#' 
#' \code{rf.cv} implements Breiman's random forest algorithm for classification and 
#' regression. here we use it to make a k-fold cross-validation
#' 
#' @param xtr A data frame or a matrix of predictors.
#' 
#' @param ytr A response vector. If a factor, classification is assumed, otherwise 
#'            regression is assumed.
#'            
#' @param cv.fold The fold, the defalut is 5.
#' 
#' @param type method type.
#' 
#' @param trees Number of trees to grow.  This should not be set to too small a 
#'              number, to ensure that every input row gets predicted at least 
#'              a few times. 
#'              
#' @param mtrysize Number of variables randomly sampled as candidates at each 
#'                 split.  Note that the default values are different for 
#'                 classification (sqrt(p) where p is number of variables in 
#'                 \code{xtr}) and regression (p/3)                             
#'
#' @return if type is regression, the retrun a list containing four components:
#' \itemize{
#' \item \code{RFpred} - the predicted values of the input data based on cross-validation
#' \item \code{Error} - error for all samples
#' \item \code{RMSECV} - Root Mean Square Error for cross-validation
#' \item \code{Q2} - R2 for cross-validation
#' }
#' if type is classification, the retrun a list containing four components:
#' \itemize{
#' \item \code{table} - confusion matrix
#' \item \code{ACC} - accuracy
#' \item \code{SE} - sensitivity
#' \item \code{SP} - specifivity
#' \item \code{F1} - a measure of a test's accuracy.
#' \item \code{MCC} - Mathews correlation coefficient
#' \item \code{RFPred} - the predicted values
#' \item \code{prob} - the predicted probability values
#' }
#' @keywords  rf.cv
#'
#' @aliases rf.cv
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{pls.cv}} for the Cross-Validation of Classification and 
#' Regression models using PLS 
#' 
#' @export rf.cv
#' 
#' @importFrom randomForest randomForest
#' 
#' @references
#' Breiman, L. (2001), \emph{Random Forests}, Machine Learning 45(1), 5-32.
#' 
#' @examples
#' training = read.csv(system.file('sysdata/training2.csv', package = 'BioMedR'), header = TRUE)
#' y = training[, 1]
#' x = training[, -1]
#' rf.tr <- rf.cv(x, y)


rf.cv <- function (xtr, ytr, cv.fold = 5, type = 'regression', 
                      trees = 500, mtrysize = 10) { 
  
  mx = dim(xtr)[1]
  rfpred <- ytr
  prob = matrix(nrow = length(ytr), ncol = 2)
  index = rep(1:cv.fold, nrow(xtr))
  ind = index[1:nrow(xtr)]
  for (k in 1:cv.fold) {
    cat(".")
    xcal <- xtr[ind != k, ] 
    ycal <- ytr[ind != k]
    xtest <- xtr[ind == k, ] 
    ytest <- ytr[ind == k]   
    rfout <- randomForest::randomForest(ycal~., data = data.frame(xcal, ycal),  
                          ntrees = trees, mtry = mtrysize,         
                          importance = FALSE)
    if (type == 'regression') {
      rfpred[ind==k] <- predict(rfout, xtest)
    } else if (type == 'classification') {
      rfpred[ind == k] = predict(rfout, xtest, type = "response")
      prob[ind == k, ] = predict(rfout, xtest, type = "prob")
    }    
  }
  cat("\n")
  if (type == 'regression') {
    RMSECV = sqrt(t(ytr - rfpred) %*% (ytr - rfpred) / mx)
    q2 = 1 - t(ytr - rfpred) %*% (ytr - rfpred) / (t(ytr - mean(ytr)) %*% (ytr - mean(ytr))) 
    err = ytr - rfpred
    ret <- list(RFpred = rfpred, Error = err, RMSECV = RMSECV, Q2 = q2)
  } else if (type == 'classification') {
    ret = list()
    r = .result(rfpred, ytr)
    ret$table = r$table
    ret$ACC = r$ACC
    ret$SE = r$SE
    ret$SP = r$SP
    ret$F1 = r$F1
    ret$MCC = r$MCC
    ret$RFPred = rfpred
    ret$prob = prob
  }
  return(ret)
}
