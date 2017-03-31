#' The Dinucleotide-based Auto-cross Covariance Descriptor
#' 
#' The Dinucleotide-based Auto-cross Covariance Descriptor
#'
#' This function calculates the dinucleotide-based auto-cross covariance descriptor
#' 
#' @param x the input data, which should be a list or file type.
#' 
#' @param index the physicochemical indices, it should be a list and there are 38
#'              different physicochemical indices (Table 1), which the users can choose.
#'
#' @param nlag an integer larger than or equal to 0 and less than or equal to L-2 (L means the length 
#'             of the shortest DNA sequence in the dataset). It represents the distance between two dinucleotides.
#'
#'  @param normaliztion with this option, the final feature vector will be normalized based
#'                  on the total occurrences of all kmers. Therefore, the elements in the feature vectors 
#'                  represent the frequencies of kmers. The default value of this parameter is False.
#' 

#' @param allprop all the 38 physicochemical indices will be
#'                employed to generate the feature vector. Its default value is False.
#'       
#' @param customprops the users can use their own indices to generate the feature vector. It should be a dict, 
#'                    the key is dinucleotide (string), and its corresponding value is a list type.                       
#' 
#' @return A vector 
#' 
#' @keywords extract DACC
#' 
#' @aliases extrDNADACC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNADACC
#' 
#' @seealso See \code{\link{extrDNADAC}} and \code{\link{extrDNADCC}}
#'          
#' @note if the user defined physicochemical indices have not been normalized, it should be normalized.
#'
#' @references 
#' Dong Q, Zhou S, Guan J. A new taxonomy-based protein fold recognition approach based on 
#' autocross-covariance transformation. \emph{Bioinformatics}, 2009, 25(20): 2655-2662.
#'
#' 
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNADACC(x)

extrDNADACC = function (x, index = c('Twist', 'Tilt'), nlag = 2, normaliztion = FALSE, 
                     customprops = NULL, allprop = FALSE) {
  dac = extrDNADAC(x, index, nlag, customprops = customprops, 
                allprop = allprop, normaliztion = normaliztion)
  dcc = extrDNADCC(x, index, nlag, customprops = customprops, 
                allprop = allprop, normaliztion = normaliztion)  
  dacc = c(dac, dcc)
  return(dacc) 
}