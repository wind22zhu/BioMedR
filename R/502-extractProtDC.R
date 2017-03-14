#' Dipeptide Composition Descriptor
#'
#' Dipeptide Composition Descriptor
#' 
#' This function calculates the Dipeptide Composition descriptor (Dim: 400).
#' 
#' @param x A character vector, as the input protein sequence. 
#'
#' @return A length 400 named vector
#' 
#' @keywords extract DC extrProtDC Dipeptide Composition
#'
#' @aliases extrProtDC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{extrProtTC}} for tripeptide composition descriptors.
#' 
#' @export extrProtDC
#' 
#' @references
#' M. Bhasin, G. P. S. Raghava.
#' Classification of Nuclear Receptors Based on 
#' Amino Acid Composition and Dipeptide Composition. 
#' \emph{Journal of Biological Chemistry}, 2004, 279, 23262.
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'BioMedR'))[[1]]
#' extrProtDC(x)
#' 

extrProtDC = function (x) {

    if (checkProt(x) == FALSE) stop('x has unrecognized amino acid type')

    AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    DCDict = as.vector((outer(AADict, AADict, paste, sep = '')))

    xSplitted = strsplit(x, split = '')[[1]]
    n  = nchar(x)
    DC = summary(factor(paste(xSplitted[-n], xSplitted[-1], sep = ''), 
                        levels = DCDict), maxsum = 401)/(n-1)

    return(DC)

}
