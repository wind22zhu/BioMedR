#' The Pseudo K-tupler Composition Descriptor
#' 
#' The Pseudo K-tupler Composition Descriptor
#'
#' This function calculates the pseudo k-tupler composition Descriptor
#' 
#' @param x the input data, which should be a list or file type.
#' 
#' @param normalize with this option, the final feature vector will be normalized based
#'                  on the total occurrences of all kmers. Therefore, the elements in the feature vectors 
#'                  represent the frequencies of kmers. The default value of this parameter is False.
#'                  
#' @param lambda an integer larger than or equal to 0 and less than or equal to L-2 (L means the length of the shortest 
#'               sequence in the dataset). It represents the highest counted rank (or tier) of the correlation along a 
#'               DNA sequence. Its default value is 3. 
#' 
#' @param k an integer larger than 0 represents the k-tuple. Its default value is 3.
#' 
#' @param w the weight factor ranged from 0 to 1. Its default value is 0.05.                                       
#' 
#' @param customprops the users can use their own indices to generate the feature vector. It should be a dict, 
#'                    the key is dinucleotide (string), and its corresponding value is a list type.                       
#' 
#' @return A vector 
#' 
#' @keywords extract PseKNC
#' 
#' @aliases extrDNAPseKNC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNAPseKNC
#' 
#' @seealso See \code{\link{extrDNAPseDNC}}
#'          
#' @note if the user defined physicochemical indices have not been normalized, it should be normalized.
#'
#' @references 
#' Guo S H, Deng E Z, Xu L Q, et al. iNuc-PseKNC: a sequence-based predictor for predicting nucleosome positioning 
#' in genomes with pseudo k-tuple nucleotide composition. \emph{Bioinformatics}, 2014: btu083.
#' 
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNAPseKNC(x)
#' 

extrDNAPseKNC = function (x, lambda = 1, k = 3, normalize = FALSE,
                                  w = 0.5, customprops = NULL) {
  
  if (checkDNA(x) == FALSE) 
    stop ("x has unrecognized  type !")
  
  Didx = data.frame(AA = c(0.06, 0.5, 0.09, 1.59, 0.11, -0.11),
                    AC = c(1.50, 0.50, 1.19, 0.13, 1.29, 1.04),
                    AG = c(0.78, 0.36, -0.28, 0.68, -0.24, -0.62),
                    AT = c(1.07, 0.22, 0.83, -1.02, 2.51, 1.17),
                    CA = c(-1.38, -1.36, -1.01, -0.86, -0.62, -1.25),
                    CC = c(0.06, 1.08, -0.28, 0.56, -0.82, 0.24),
                    CG = c(-1.66, -1.22, -1.38, -0.82, -0.29, -1.39),
                    CT = c(0.78, 0.36, -0.28, 0.68, -0.24, -0.62),
                    GA = c(-0.08, 0.5, 0.09, 0.13, -0.39, 0.71),
                    GC = c(-0.08, 0.22, 2.3, -0.35, 0.65, 1.59),
                    GG = c(0.06, 1.08, -0.28, 0.56, -0.82, 0.24),
                    GT = c(1.50, 0.50, 1.19, 0.13, 1.29, 1.04),
                    TA = c(-1.23, -2.37, -1.38, -2.24, -1.51, -1.39),
                    TC = c(-0.08, 0.5, 0.09, 0.13, -0.39, 0.71),
                    TG = c(-1.38, -1.36, -1.01, -0.86, -0.62, -1.25),
                    TT = c(0.06, 0.5, 0.09, 1.59, 0.11, -0.11))
 
  Ddict = make_kmer_index(k = 2)
  if (!is.null(customprops)) {
    if (normalize) {
      n0 = dim(customprops)[1]
      H0 = as.matrix (customprops)
      H = matrix(ncol = 16, nrow = n0)
      for (i in 1:n0) H[i, ] = (H0[i, ] - mean(H0[i, ]))/(sqrt(sum((H0[i,] - 
                                                                     mean(H0[i, ])) ^ 2) / 16))
      colnames(H) = Ddict
    }
    Didx = rbind(Didx, H)
  }
  n = dim(Didx)[1]  
  
  Theta = vector("list", lambda)
  for (i in 1:lambda) Theta[[i]] = vector("list", n)
  xSplit = strsplit(x, split = "")[[1]]
  xPaste = paste0(xSplit[1:length(xSplit) - 1], xSplit[2:length(xSplit)])
  N = length(xPaste)
  for (i in 1:lambda) {
    temp = c()
    for (j in 1:(N - lambda)) {
     temp = append(temp,mean((Didx[, xPaste[j]] - Didx[, xPaste[j + i]])^2))
    }
    Theta[[i]] = sum(temp)/(N-i)
  }
  theta = unlist(Theta)
  xPaste2 = c()
  for(i in 1:k){
    temp = xSplit[i:(nchar(x) - k + i)]
    xPaste2 = paste0(xPaste2, temp)
  }
  
  fc = summary(factor(xPaste2, levels = make_kmer_index(k)), maxsum = 4 ^ k)
  fc = fc/sum(fc)
  
  Xc1 = fc/(1 + (w * sum(theta)))
  names(Xc1) = paste("Xc1.", names(Xc1), sep = "")
  Xc2 = (w * theta)/(1 + (w * sum(theta)))
  names(Xc2) = paste("Xc2.lambda.", 1:lambda, sep = "")
  Xc = c(Xc1, Xc2)
  
  Xc = round(Xc, 3)
  return(Xc)
}