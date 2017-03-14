#' The Reverse chars
#' 
#' The Reverse chars
#'
#' This function calculates Reverse chars
#' 
#' @param x the input data, which should be a string.                      
#' 
#' @return A vector 
#' 
#' @keywords extract reverse_chars
#' 
#' @aliases revchars
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export revchars
#'          
#' @note if the user defined physicochemical indices have not been normalized, it should be normalized.
#'
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' revchars(x)
#' 

revchars = function (x) {
  reversed_split = rev(strsplit(as.character(x), split = "")[[1]])
  reverse_chars = paste(unlist(reversed_split), collapse = "")
  return(reverse_chars)
}