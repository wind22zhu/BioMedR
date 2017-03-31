#' Fingerprints from descriptor vectors
#'
#' Fingerprints from descriptor vectors
#'
#' Generates fingerprints from descriptor vectors such as atom pairs stored in \code{APset} or \code{list} 
#' containers. The obtained fingerprints can be used for structure similarity comparisons, searching and 
#' clustering. Due to their compact size, computations on fingerprints are often more time and memory 
#' efficient than on their much more complex atom pair counterparts.
#' 
#' @param x Object of classe \code{APset} or \code{list} of vectors
#' 
#' @param descnames Descriptor set to consider for fingerprint encoding. If a single value from 1-4096 is 
#'                  provided then the function uses the corresponding number of the most frequent atom pairs 
#'                  stored in the \code{apfp} data set provided by the package. Alternatively, one can provide 
#'                  here any custom atom pair selection in form of a \code{character} vector.
#'                                  
#' @return A FPset
#' 
#' @keywords convAPtoFP
#'
#' @aliases convAPtoFP
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{convSDFtoAP}} for Atom pair library.
#' 
#' @export convAPtoFP
#' 
#' @examples
#' data(sdfbcl)
#' apbcl = convSDFtoAP(sdfbcl)
#' fpbcl = convAPtoFP(apbcl)
#' 
convAPtoFP <- function (x, descnames = 1024) {
  fpset = ChemmineR::desc2fp(x, descnames = descnames, type = 'FPset')
}