#' Amino Acid Composition Descriptor
#'
#' Amino Acid Composition Descriptor
#' 
#' This function calculates the Amino Acid Composition descriptor (Dim: 20).
#' 
#' @param x A character vector, as the input protein sequence. 
#'
#' @return A length 20 named vector
#' 
#' @keywords extract AAC extrProtACC Amino Acid Composition
#'
#' @aliases extrProtAAC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>,
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{extrProtDC}} and \code{\link{extrProtTC}} 
#'          for Dipeptide Composition and Tripeptide Composition descriptors.
#' 
#' @export extrProtAAC
#' 
#' @references
#' M. Bhasin, G. P. S. Raghava.
#' Classification of Nuclear Receptors Based on 
#' Amino Acid Composition and Dipeptide Composition. 
#' \emph{Journal of Biological Chemistry}, 2004, 279, 23262.
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'BioMedR'))[[1]]
#' extrProtAAC(x)
#' 

extrProtAAC = function (x) {

    if (checkProt(x) == FALSE) stop('x has unrecognized amino acid type')

    AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

    AAC = summary(factor(strsplit(x, split = '')[[1]], levels = AADict), 
                  maxsum = 21)/nchar(x)

    return(AAC)

}
