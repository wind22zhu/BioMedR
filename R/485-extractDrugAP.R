#' Calculate the Atom Pair Fingerprints
#' 
#' Calculate the Atom Pair Fingerprints
#' 
#' Generates fingerprints from descriptor vectors such as atom pairs 
#' stored in \code{APset} or \code{list} containers. The obtained 
#' fingerprints can be used for structure similarity comparisons, 
#' searching and clustering. Due to their compact size, computations 
#' on fingerprints are often more time and memory efficient than on 
#' their much more complex atom pair counterparts.
#' 
#' @param x Object of classe \code{APset} or \code{list} of vectors
#'
#' @param descnames Descriptor set to consider for fingerprint encoding. If a single value 
#'                  from 1-4096 is provided then the function uses the corresponding number 
#'                  of the most frequent atom pairs stored in the \code{apfp} data set provided 
#'                  by the package. Alternatively, one can provide here any custom atom pair 
#'                  selection in form of a \code{character} vector.
#'  
#' @return matrix or character vectors
#' 
#' @keywords extrDrugAP
#'
#' @aliases extrDrugAP
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#' 
#' @export extrDrugAP
#' 
#' @examples
#' data(sdfbcl)
#'
#' apbcl <- convSDFtoAP(sdfbcl)
#' mol <- extrDrugAP(x = apbcl, descnames = 1024)
#' 
extrDrugAP <- function (x, descnames = 1024) {
  
    if (length(descnames) == 1) {
        data(apfp)
        descnames <- as.character(apfp$AP)[1:descnames]
    }	
    if (class(x)=="APset") { 
        apfp <- matrix(0, nrow = length(x), ncol = length(descnames), 
                                     dimnames = list(cid(x), descnames))
        apsetlist <- ap(x)
        for (i in cid(x)) apfp[i, descnames %in% as.character(apsetlist[[i]])] <- 1
    } else if (class(x) == "list") {
        apfp <- matrix(0, nrow = length(x), ncol = length(descnames), 
                                     dimnames = list(names(x), descnames))
        for (i in names(x)) apfp[i, descnames %in% as.character(x[[i]])] <- 1
    } else {
        stop ("x needs to be of class APset or list")
    }
    return(apfp)
}
