#' The Increment Of Diversity Descriptors
#' 
#' The Increment Of Diversity Descriptors
#'
#' This function calculates the The Basic Kmer Descriptor 
#' 
#' @param k the k value of kmer, it should be an integer larger than 0,the default value is 6.
#'
#' @param x the input data, which should be a list or file type.
#' 
#' @param pos the positive source data, which should be a or type.
#' 
#' @param neg the negative source data, which should be or type.
#' 
#' @param upto generate all the kmers: 1mer, 2mer, ..., kmer. The output feature vector is 
#'             the combination of all these kmers. The default value of this parameter is True
#' 
#' @return if upto is True, A length \code{k * 2} named vector, \code{k} is the k value of kmer;
#'         if upto is False, A length 2 named vector
#' 
#' @keywords extract the increment of diversity 
#' 
#' @aliases IncDiv
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNAIncDiv
#' 
#' @seealso See \code{\link{extrDNAkmer}} 
#'
#' @references 
#' Chen W, Luo L, Zhang L. The organization of nucleosomes around splice sites. 
#' \emph{Nucleic acids research}, 2010, 38(9): 2788-2798.
#' Liu G, Liu J, Cui X, et al. Sequence-dependent prediction of recombination hotspots in
#' Saccharomyces cerevisiae. \emph{Journal of theoretical biology}, 2012, 293: 49-54.
#' 
#' @examples
#'  
#' pos = readFASTA(system.file('dnaseq/pos.fasta', package = 'BioMedR'))
#' neg = readFASTA(system.file('dnaseq/neg.fasta', package = 'BioMedR'))
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNAIncDiv(k = 6, x, pos, neg)



extrDNAIncDiv = function (k = 6, x, pos, neg,
                                 upto = TRUE) {
  
  if (upto) l = 1 
    else l = k
  ID = c()
  for (j in l:k){
    if (upto == FALSE || j == 1){
      temp_pos_s_vec = data.frame(lapply(pos, extrDNAkmer, k=j))
      temp_neg_s_vec = data.frame(lapply(neg, extrDNAkmer, k=j))
      temp_pos_s_vec = rowSums(temp_pos_s_vec)
      temp_neg_s_vec = rowSums(temp_neg_s_vec)
    }
    temp_vec = lapply(x, extrDNAkmer, k = j) 
    temp_id.pos = lapply(temp_vec, id_x_s, vec_s = temp_pos_s_vec, diversity_s = diversity(temp_pos_s_vec))
    temp_id.neg = lapply(temp_vec, id_x_s, vec_s = temp_neg_s_vec, diversity_s = diversity(temp_neg_s_vec))
    temp_id = mapply(c, temp_id.pos, temp_id.neg, SIMPLIFY = FALSE)
    if (length(ID) == 0) {
      ID = temp_id
    } else {
      ID = mapply(c, ID, temp_id, SIMPLIFY = FALSE)
    }
  }
  return(ID)
}

#' diversity function for make_kmer_vec
#'
#' diversity
#'
#' @return diversity
#' 
#' @keywords internal

diversity = function (vec) {
  vec = vec[which(vec != 0)]
  diversity = sum(vec) * log2(sum(vec)) - sum(vec * log2(vec))
  return(diversity)
}

#' id_x_s function for make_kmer_vec
#'
#' id_x_s
#'
#' @return id
#' 
#' @keywords internal
#' 

id_x_s = function (vec_x, vec_s, diversity_s) {
  vec_x_s = rowSums(data.frame(vec_x[1:length(vec_s)], vec_s))
  id_x_s = diversity(vec_x_s) - diversity(vec_x) - diversity_s
  return(id_x_s)
}
