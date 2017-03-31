#' Calculate The Basic Kmer Feature Vector
#'
#' Calculate The Basic Kmer Feature Vector
#'
#' This function calculate the basic kmer feature vector.
#'
#' @param k the k value of kmer, it should be an integer larger than 0. 
#' 
#' @param alphabet the 
#' 
#' @return The result character vector
#'
#' @keywords kmer index
#'
#' @aliases make_kmer_index
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#'
#' @seealso See \code{\link{extrDNAkmer}} 
#'
#' @export make_kmer_index
#'
#'
#' @examples
#' make_kmer_index(2, alphabet = "ACGT")

make_kmer_index = function (k, alphabet = 'ACGT') {
  
  dict = unlist(strsplit(alphabet, ''))
  make_index = list()
  temp = rep(dict, k)
  make_index = sort(unique(apply(expand.grid(split(temp,
                                                   rep(1:k,each=4))), 1, paste, collapse = '')))
  return(make_index)
  
}


