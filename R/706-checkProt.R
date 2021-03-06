#' Check if the protein sequence's amino acid types are the 20 default types
#'
#' Check if the protein sequence's amino acid types are the 20 default types
#' 
#' This function checks if the protein sequence's amino acid types are 
#' the 20 default types.
#' 
#' @param x A character vector, as the input protein sequence.
#'
#' @return Logical. \code{TRUE} if all of the amino acid types of the sequence
#'         are within the 20 default types.
#' 
#' @keywords protein sequence amino acid type check
#'
#' @aliases checkProt
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export checkProt
#' 
#' @examples
#' x = readFASTA(system.file('protseq/P00750.fasta', package = 'BioMedR'))[[1]]
#' checkProt(x)  # TRUE
#' checkProt(paste(x, 'Z', sep = ''))  # FALSE
#' 

checkProt = function (x) {

    AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

    return(all(strsplit(x, split = '')[[1]] %in% AADict))

}
