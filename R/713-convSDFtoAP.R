#' Atom pair library
#' 
#' Atom pair library
#' 
#' Creates from a SDFset a searchable atom pair library that is stored in a 
#' container of class APset.
#' 
#' @param sdfset Objects of classes \code{SDFset} or \code{SDF}
#'
#' @param type if \code{type="AP"}, the function returns \code{APset}/\code{AP} 
#'             objects; if \code{type="character"}, it returns the result as a 
#'             \code{character} vector of length one. The latter is useful for 
#'             storing AP data in tabular files.
#' 
#' @param uniquePairs When the same atom pair occurs more than once in a single compound, should the 
#'                    names be unique or not? Setting this to true will take slightly longer to compute.
#' 
#' @return \item{APset}{ if input is \code{SDFset}}
#'         \item{AP}{ if input is \code{SDF}}
#' 
#' @keywords convSDFtoAP 
#'
#' @aliases convSDFtoAP
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#' 
#' @export convSDFtoAP
#' 
#' @examples
#' data(sdfbcl)
#'
#' apset <- convSDFtoAP(sdfbcl)
#' 
convSDFtoAP <- function(sdfset, type="AP",uniquePairs=TRUE) {
  apset = ChemmineR::sdf2ap(sdfset)
}