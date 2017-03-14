#' The Basic Kmer Descriptor
#' 
#' The Basic Kmer Descriptor
#'
#' This function calculates the basic kmer descriptor 
#' 
#' @param x the input data, which should be a list or file type.
#' 
#' @param k the k value of kmer, it should be an integer larger than 0.
#'
#' @param normalize with this option, the final feature vector will be normalized based
#'                  on the total occurrences of all kmers. Therefore, the elements in the feature vectors 
#'                  represent the frequencies of kmers. The default value of this parameter is False.
#' 
#' @param upto generate all the kmers: 1mer, 2mer, ..., kmer. The output feature vector is 
#'             the combination of all these kmers. The default value of this parameter is False.
#'
#' @param reverse make reverse complements into a single feature, The default value of this parameter is False.
#'                if reverse is True, this method returns the reverse compliment kmer feature vector.
#' 
#' @return A vector 
#' 
#' @keywords extract kmer
#' 
#' @aliases extrDNAkmer
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNAkmer
#' 
#' @seealso See \code{\link{make_kmer_index}}
#'          
#' @note if the parameters normalize and upto are both True, and then the feature vector is 
#'       the combination of all these normalized kmers, e.g. the combination of normalized 1-kmer 
#'       and normalized 2-kmer when k=2, normalize=True, upto=True.
#'
#' @references 
#' Noble W S, Kuehn S, Thurman R, et al. Predicting the in vivo signature of human gene regulatory sequences. 
#' \emph{Bioinformatics}, 2005, 21 Suppl 1, i338-343.
#' Lee D, Karchin R, Beer M A. Discriminative prediction of mammalian enhancers from 
#' DNA sequence. \emph{Genome research}. 2005, 21, 2167-2180.
#'
#' 
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNAkmer(x)
#'

extrDNAkmer = function (x, k = 2, upto = FALSE, 
                 normalize = FALSE, reverse = FALSE) {
  
  dict = c("A", "C", "G", "T")
  make_index = list()
  kmer = c()
  if (upto) l = 1 else l = k
  
  for (j in l:k) {
    make_index = make_kmer_index(j, alphabet = 'ACGT')
    xPaste = c()
    n = nchar(x)
    xSplit = strsplit(x, split = "")[[1]]
    for (i in 1:j) {
      temp = xSplit[i:(n-j+i)]
      xPaste = paste0(xPaste, temp)
    }
    temp_kmer = summary(factor(xPaste, levels = make_index), maxsum = 4^j)
    if (reverse) {
      reverse_index = sapply(make_index, revchars)
      reverse_index = chartr('ACGT', 'TGCA', reverse_index)
      temp_kmer_reverse = summary(factor(xPaste, levels = reverse_index), maxsum = 4^j)
      cmp = which((make_index > reverse_index) - (make_index < reverse_index) < 0)
      temp_kmer[cmp] = temp_kmer[cmp] + temp_kmer_reverse[cmp]
      cmp = which((make_index > reverse_index) - (make_index < reverse_index) <= 0)
      temp_kmer = temp_kmer[cmp]
    }
    if (normalize) temp_kmer = temp_kmer/sum(temp_kmer)
    kmer = append(kmer, temp_kmer)
  }
  return(kmer)
}

