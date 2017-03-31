#' Tripeptide Composition Descriptor
#'
#' Tripeptide Composition Descriptor
#' 
#' This function calculates the Tripeptide Composition descriptor (Dim: 8000).
#' 
#' @param x A character vector, as the input protein sequence. 
#'
#' @return A length 8000 named vector
#' 
#' @keywords extract TC extrProtTC Tripeptide Composition
#'
#' @aliases extrProtTC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See  \code{\link{extrProtDC}} for  dipeptide composition descriptors.
#' 
#' @export extrProtTC
#' 
#' @references
#' M. Bhasin, G. P. S. Raghava.
#' Classification of Nuclear Receptors Based on 
#' Amino Acid Composition and Dipeptide Composition. 
#' \emph{Journal of Biological Chemistry}, 2004, 279, 23262.
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'BioMedR'))[[1]]
#' extrProtTC(x)
#' 

extrProtTC = function (x) {

    if (checkProt(x) == FALSE) stop('x has unrecognized amino acid type')

    AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    DCDict = as.vector((outer(AADict, AADict, paste, sep = '')))
    TCDict = as.vector((outer(DCDict, AADict, paste, sep = '')))

    xSplitted = strsplit(x, split = '')[[1]]
    n  = nchar(x)
    TC = summary(factor(paste(paste(xSplitted[-c(n, n-1)], 
                                    xSplitted[-c(1, n)], sep = ''), 
                              xSplitted[-c(1, 2)], sep = ''), 
                        levels = TCDict), maxsum = 8001)/(n-2)

    return(TC)

}
