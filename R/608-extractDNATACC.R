#' The Trinucleotide-based Auto-cross Covariance Descriptor
#' 
#' The Trinucleotide-based Auto-cross Covariance Descriptor
#'
#' This function calculates the trinucleotide-based auto-cross covariance descriptor
#' 
#' @param x the input data, which should be a list or file type.
#' 
#' @param index the physicochemical indices, it should be a list and there are 12
#'              different physicochemical indices (Table 2), which the users can choose.
#'
#' @param nlag an integer larger than or equal to 0 and less than or equal to L-2 (L means the length 
#'             of the shortest DNA sequence in the dataset). It represents the distance between two dinucleotides.
#'
#' @param normaliztion with this option, the final feature vector will be normalized based
#'                  on the total occurrences of all kmers. Therefore, the elements in the feature vectors 
#'                  represent the frequencies of kmers. The default value of this parameter is False.
#' 

#' @param allprop all the 12 physicochemical indices will be
#'                employed to generate the feature vector. Its default value is False.
#'       
#' @param customprops the users can use their own indices to generate the feature vector. It should be a dict, 
#'                    the key is dinucleotide (string), and its corresponding value is a list type.                       
#' 
#' @return A vector 
#' 
#' @keywords extract TACC
#' 
#' @aliases extrDNATACC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNATACC
#' 
#' @seealso See \code{\link{extrDNATAC}} and \code{\link{extrDNATCC}}
#'          
#' @note if the user defined physicochemical indices have not been normalized, it should be normalized.
#' 
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNATACC(x)
#'

extrDNATACC = function (x, index = c('Dnase I', 'Nucleosome'), nlag = 2, normaliztion = FALSE, 
                     customprops = NULL, allprop = FALSE) {
  tac = extrDNATAC(x, index, nlag, customprops = customprops, 
                allprop = allprop, normaliztion = normaliztion)
  tcc = extrDNATCC(x, index, nlag, customprops = customprops,
                normaliztion = normaliztion)  
  tacc = c(tac, tcc)
  return(tacc) 
}