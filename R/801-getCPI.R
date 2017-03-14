.inner = function (a, b, f, ...) {

    # For computing column-by-column (pseudo)-tensor product type interactions

    f = match.fun(f)
    apply(b, 2, function (x) apply(a, 1, function (y) f(x, y, ...)))

}

.getCPICombine = function (drugmat, protmat) {

    return(cbind(drugmat, protmat))

}

.getCPITensor = function (drugmat, protmat, row, dcol, pcol) {

    return(array(.inner(t(drugmat), protmat, '*'), c(row, dcol * pcol)))

}

.getPPICombine = function (protmat1, protmat2) {
  
  return(cbind(protmat1, protmat2))
  
}

.getPPITensor = function (protmat1, protmat2, row, col) {
  
  return(array(.inner(t(protmat1), protmat2, '*'), c(row, col^2)))
  
}

.getPPIEntry = function (protmat1, protmat2) {
  
  return(cbind((protmat1 * protmat2), (protmat1 + protmat2)))
  
}

.getDDICombine = function (DNAmat1, DNAmat2) {
  
  return(cbind(DNAmat1, DNAmat2))
  
}

.getDDITensor = function (DNAmat1, DNAmat2, row, col) {
  
  return(array(.inner(t(DNAmat1), DNAmat2, '*'), c(row, col^2)))
  
}

.getDDIEntry = function (DNAmat1, DNAmat2) {
  
  return(cbind((DNAmat1 * DNAmat2), (DNAmat1 + DNAmat2)))
  
}

.getDPIinner = function (a, b, f, ...) {
  
  # For computing column-by-column (pseudo)-tensor product type interactions
  
  f = match.fun(f)
  apply(b, 2, function (x) apply(a, 1, function (y) f(x, y, ...)))
  
}

.getDPICombine = function (DNAmat, protmat) {
  
  return(cbind(DNAmat, protmat))
  
}

.getDPITensor = function (DNAmat, protmat, row, dcol, pcol) {
  
  return(array(.getDPIinner(t(DNAmat), protmat, '*'), c(row, dcol * pcol)))
  
}

.getCCICombine = function (drugmat1, drugmat2) {
  
  return(cbind(drugmat1, drugmat2))
  
}

.getCCITensor = function (drugmat1, drugmat2, row, col) {
  
  return(array(.inner(t(drugmat1), drugmat2, '*'), c(row, col^2)))
  
}

.getCCIEntry = function (drugmat1, drugmat2) {
  
  return(cbind((drugmat1 * drugmat2), (drugmat1 + drugmat2)))
  
}

.getCDIinner = function (a, b, f, ...) {
  
  # For computing column-by-column (pseudo)-tensor product type interactions
  
  f = match.fun(f)
  apply(b, 2, function (x) apply(a, 1, function (y) f(x, y, ...)))
  
}

.getCDICombine = function (drugmat, DNAmat) {
  
  return(cbind(drugmat, DNAmat))
  
}

.getCDITensor = function (drugmat, DNAmat, row, dcol, pcol) {
  
  return(array(.getCDIinner(t(drugmat), DNAmat, '*'), c(row, dcol * pcol)))
  
}

#' @title Generating Interaction Descriptors
#'
#' @description Generating Interaction Descriptors
#' 
#' @details This function calculates the interaction descriptors
#' by three types of interaction:
#' \itemize{
#' \item \code{combine} - combine the two descriptor matrix, 
#' result has \code{(p1 + p2)} columns
#' \item \code{tensorprod} - calculate column-by-column 
#' (pseudo)-tensor product type interactions, result has 
#' \code{(p1 * p2)} columns
#' \item \code{entrywise} - calculate entrywise product and 
#' entrywise sum of the two matrices, then combine them, 
#' result has \code{(p + p)} columns
#' }
#' 
#' @param drugmat The compound descriptor matrix.
#' @param protmat The protein descriptor matrix.
#' @param type The interaction type, one or more of  
#' \code{"combine"}, \code{"tensorprod"}, and \code{"entrywise"}.
#' 
#' @return A matrix containing the interaction descriptors
#' 
#' @keywords getCPI getCDI getDPI getPPI getDDI getCCI interaction 
#'
#' @aliases getCPI
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export getCPI
#' 
#' @examples
#' x = matrix(1:10, ncol = 2)
#' y = matrix(1:15, ncol = 3)
#' # getCPI 
#' getCPI(x, y, 'combine')
#' # getCDI
#' getCDI(x, y, 'tensorprod')
#' # getDPI
#' getDPI(x, y, type = c('combine', 'tensorprod'))
#' getDPI(x, y, type = c('tensorprod', 'combine'))
#' 

getCPI = function (drugmat, protmat, type = c('combine', 'tensorprod')) {

    if (!is.matrix(drugmat)) drugmat = as.matrix(drugmat)
    if (!is.matrix(protmat)) protmat = as.matrix(protmat)

    drugrow = nrow(drugmat)
    protrow = nrow(protmat)

    if (drugrow != protrow) stop('Matrix row count must match')

    drugcol = ncol(drugmat)
    protcol = ncol(protmat)

    if (missing(type)) stop('Must provide at least one interaction type')

    if (all(type == 'combine')) {

        result = .getCPICombine(drugmat, protmat)

        } else if (all(type == 'tensorprod')) {

            result = .getCPITensor(drugmat, protmat, row = drugrow, 
                                   dcol = drugcol, pcol = protcol)

            } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {

                result = cbind(.getCPICombine(drugmat, protmat), 
                               .getCPITensor(drugmat, protmat, row = drugrow, 
                                             dcol = drugcol, pcol = protcol))

                } else {
                    
                    stop('Interaction type must be in "tensorprod" and "combine" or both')

                }

    return(result)

}


#' @rdname getCPI
#'
#' @aliases getPPI
#' 
#' @param protmat1 The first protein descriptor matrix, 
#' must have the same ncol with \code{protmat2}.
#' @param protmat2 The second protein descriptor matrix, 
#' must have the same ncol with \code{protmat1}.
#' 
#' @export getPPI
#' 
#' @examples
#' x = matrix(1:10, ncol = 2)
#' y = matrix(5:14, ncol = 2)
#' 
#' # getPPI 
#' getPPI(x, y, type = 'combine')
#' getPPI(x, y, type = 'tensorprod')
#' # getDDI
#' getDDI(x, y, type = 'entrywise')
#' getDDI(x, y, type = c('combine', 'tensorprod'))
#' # getCCI
#' getCCI(x, y, type = c('combine', 'entrywise'))
#' getCCI(x, y, type = c('entrywise', 'tensorprod'))
#' getCCI(x, y, type = c('combine', 'entrywise', 'tensorprod'))

getPPI = function (protmat1, protmat2, type = c('combine', 'tensorprod', 
                                                'entrywise')) {
  
  if (!is.matrix(protmat1)) protmat1 = as.matrix(protmat1)
  if (!is.matrix(protmat2)) protmat2 = as.matrix(protmat2)
  
  protrow1 = nrow(protmat1)
  protrow2 = nrow(protmat2)
  
  if (protrow1 != protrow2) stop('Matrix row count must match')
  
  protcol1 = ncol(protmat1)
  protcol2 = ncol(protmat2)
  
  if (protcol1 != protcol2) stop('Matrix column count must match')
  
  if (missing(type)) stop('Must provide at least one interaction type')
  
  if (all(type == 'combine')) {
    
    result = .getPPICombine(protmat1, protmat2)
    
  } else if (all(type == 'tensorprod')) {
    
    result = .getPPITensor(protmat1, protmat2, 
                           row = protrow1, col = protcol1)
    
  } else if (all(type == 'entrywise')) {
    
    result = .getPPIEntry(protmat1, protmat2)
    
  } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {
    
    result = cbind(.getPPICombine(protmat1, protmat2), 
                   .getPPITensor(protmat1, protmat2, 
                                 row = protrow1, col = protcol1))
    
  } else if (length(setdiff(type, c('tensorprod', 'entrywise'))) == 0L) {
    
    result = cbind(.getPPITensor(protmat1, protmat2, 
                                 row = protrow1, 
                                 col = protcol1), 
                   .getPPIEntry(protmat1, protmat2))
    
  } else if (length(setdiff(type, c('combine', 'entrywise'))) == 0L) {
    
    result = cbind(.getPPICombine(protmat1, protmat2), 
                   .getPPIEntry(protmat1, protmat2))
    
  } else if (length(setdiff(type, c('tensorprod', 
                                    'combine', 
                                    'entrywise'))) == 0L) {
    
    result = cbind(.getPPICombine(protmat1, protmat2), 
                   .getPPITensor(protmat1, protmat2, 
                                 row = protrow1, 
                                 col = protcol1), 
                   .getPPIEntry(protmat1, protmat2))
    
  } else {
    
    stop('Interaction type must be in "tensorprod", "combine" and "entrywise"')
    
  }
  
  return(result)
  
}

#' @rdname getCPI
#' 
#' @aliases getDDI
#' 
#' @param DNAmat1 The first DNA descriptor matrix, 
#' must have the same ncol with \code{DNAmat2}.
#' @param DNAmat2 The second DNA descriptor matrix, 
#' must have the same ncol with \code{DNAmat1}.
#' 
#' @export getDDI
#' 

getDDI = function (DNAmat1, DNAmat2, type = c('combine', 'tensorprod', 
                                              'entrywise')) {
  
  if (!is.matrix(DNAmat1)) DNAmat1 = as.matrix(DNAmat1)
  if (!is.matrix(DNAmat2)) DNAmat2 = as.matrix(DNAmat2)
  
  DNArow1 = nrow(DNAmat1)
  DNArow2 = nrow(DNAmat2)
  
  if (DNArow1 != DNArow2) stop('Matrix row count must match')
  
  DNAcol1 = ncol(DNAmat1)
  DNAcol2 = ncol(DNAmat2)
  
  if (DNAcol1 != DNAcol2) stop('Matrix column count must match')
  
  if (missing(type)) stop('Must provide at least one interaction type')
  
  if (all(type == 'combine')) {
    
    result = .getDDICombine(DNAmat1, DNAmat2)
    
  } else if (all(type == 'tensorprod')) {
    
    result = .getDDITensor(DNAmat1, DNAmat2, 
                           row = DNArow1, col = DNAcol1)
    
  } else if (all(type == 'entrywise')) {
    
    result = .getDDIEntry(DNAmat1, DNAmat2)
    
  } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {
    
    result = cbind(.getDDICombine(DNAmat1, DNAmat2), 
                   .getDDITensor(DNAmat1, DNAmat2, 
                                 row = DNArow1, col = DNAcol1))
    
  } else if (length(setdiff(type, c('tensorprod', 'entrywise'))) == 0L) {
    
    result = cbind(.getDDITensor(DNAmat1, DNAmat2, 
                                 row = DNArow1, 
                                 col = DNAcol1), 
                   .getDDIEntry(DNAmat1, DNAmat2))
    
  } else if (length(setdiff(type, c('combine', 'entrywise'))) == 0L) {
    
    result = cbind(.getDDICombine(DNAmat1, DNAmat2), 
                   .getDDIEntry(DNAmat1, DNAmat2))
    
  } else if (length(setdiff(type, c('tensorprod', 
                                    'combine', 
                                    'entrywise'))) == 0L) {
    
    result = cbind(.getDDICombine(DNAmat1, DNAmat2), 
                   .getDDITensor(DNAmat1, DNAmat2, 
                                 row = DNArow1, 
                                 col = DNAcol1), 
                   .getDDIEntry(DNAmat1, DNAmat2))
    
  } else {
    
    stop('Interaction type must be in "tensorprod", "combine" and "entrywise"')
    
  }
  
  return(result)
  
}
#' @rdname getCPI
#' 
#' @aliases getDPI
#' 
#' @param DNAmat The DNA descriptor matrix.
#' 
#' @export getDPI
#' 

getDPI = function (DNAmat, protmat, type = c('combine', 'tensorprod')) {
  
  if (!is.matrix(DNAmat)) DNAmat = as.matrix(DNAmat)
  if (!is.matrix(protmat)) protmat = as.matrix(protmat)
  
  DNArow = nrow(DNAmat)
  protrow = nrow(protmat)
  
  if (DNArow != protrow) stop('Matrix row count must match')
  
  DNAcol = ncol(DNAmat)
  protcol = ncol(protmat)
  
  if (missing(type)) stop('Must provide at least one interaction type')
  
  if (all(type == 'combine')) {
    
    result = .getDPICombine(DNAmat, protmat)
    
  } else if (all(type == 'tensorprod')) {
    
    result = .getDPITensor(DNAmat, protmat, row = DNArow, 
                           dcol = DNAcol, pcol = protcol)
    
  } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {
    
    result = cbind(.getDPICombine(DNAmat, protmat), 
                   .getDPITensor(DNAmat, protmat, row = DNArow, 
                                 dcol = DNAcol, pcol = protcol))
    
  } else {
    
    stop('Interaction type must be in "tensorprod" and "combine" or both')
    
  }
  
  return(result)
  
}

#' @rdname getCPI
#' 
#' @aliases getCCI
#' 
#' @param drugmat1 The first compound descriptor matrix, 
#' must have the same ncol with \code{drugmat2}.
#' @param drugmat2 The second compound descriptor matrix, 
#' must have the same ncol with \code{drugmat1}.
#' 
#' @export getCCI
#' 


getCCI = function (drugmat1, drugmat2, type = c('combine', 'tensorprod', 
                                                'entrywise')) {
  
  if (!is.matrix(drugmat1)) drugmat1 = as.matrix(drugmat1)
  if (!is.matrix(drugmat2)) drugmat2 = as.matrix(drugmat2)
  
  drugrow1 = nrow(drugmat1)
  drugrow2 = nrow(drugmat2)
  
  if (drugrow1 != drugrow2) stop('Matrix row count must match')
  
  drugcol1 = ncol(drugmat1)
  drugcol2 = ncol(drugmat2)
  
  if (drugcol1 != drugcol2) stop('Matrix column count must match')
  
  if (missing(type)) stop('Must provide at least one interaction type')
  
  if (all(type == 'combine')) {
    
    result = .getCCICombine(drugmat1, drugmat2)
    
  } else if (all(type == 'tensorprod')) {
    
    result = .getCCITensor(drugmat1, drugmat2, 
                           row = drugrow1, col = drugcol1)
    
  } else if (all(type == 'entrywise')) {
    
    result = .getCCIEntry(drugmat1, drugmat2)
    
  } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {
    
    result = cbind(.getCCICombine(drugmat1, drugmat2), 
                   .getCCITensor(drugmat1, drugmat2, 
                                 row = drugrow1, col = drugcol1))
    
  } else if (length(setdiff(type, c('tensorprod', 'entrywise'))) == 0L) {
    
    result = cbind(.getCCITensor(drugmat1, drugmat2, 
                                 row = drugrow1, 
                                 col = drugcol1), 
                   .getCCIEntry(drugmat1, drugmat2))
    
  } else if (length(setdiff(type, c('combine', 'entrywise'))) == 0L) {
    
    result = cbind(.getCCICombine(drugmat1, drugmat2), 
                   .getCCIEntry(drugmat1, drugmat2))
    
  } else if (length(setdiff(type, c('tensorprod', 
                                    'combine', 
                                    'entrywise'))) == 0L) {
    
    result = cbind(.getCCICombine(drugmat1, drugmat2), 
                   .getCCITensor(drugmat1, drugmat2, 
                                 row = drugrow1, 
                                 col = drugcol1), 
                   .getCCIEntry(drugmat1, drugmat2))
    
  } else {
    
    stop('Interaction type must be in "tensorprod", "combine" and "entrywise"')
    
  }
  
  return(result)
  
}

#' @rdname getCPI
#' 
#' @aliases getCDI
#' 
#' @export getCDI
#' 


getCDI = function (drugmat, DNAmat, type = c('combine', 'tensorprod')) {
  
  if (!is.matrix(drugmat)) drugmat = as.matrix(drugmat)
  if (!is.matrix(DNAmat)) DNAmat = as.matrix(DNAmat)
  
  drugrow = nrow(drugmat)
  DNArow = nrow(DNAmat)
  
  if (drugrow != DNArow) stop('Matrix row count must match')
  
  drugcol = ncol(drugmat)
  DNAcol = ncol(DNAmat)
  
  if (missing(type)) stop('Must provide at least one interaction type')
  
  if (all(type == 'combine')) {
    
    result = .getCDICombine(drugmat, DNAmat)
    
  } else if (all(type == 'tensorprod')) {
    
    result = .getCDITensor(drugmat, DNAmat, row = drugrow, 
                           dcol = drugcol, pcol = DNAcol)
    
  } else if (length(setdiff(type, c('tensorprod', 'combine'))) == 0L) {
    
    result = cbind(.getCDICombine(drugmat, DNAmat), 
                   .getCDITensor(drugmat, DNAmat, row = drugrow, 
                                 dcol = drugcol, pcol = DNAcol))
    
  } else {
    
    stop('Interaction type must be in "tensorprod" and "combine" or both')
    
  }
  
  return(result)
  
}

