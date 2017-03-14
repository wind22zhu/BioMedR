#' The Pseudo Dinucleotide Composition Descriptor
#' 
#' The Pseudo Dinucleotide Composition Descriptor
#'
#' This function calculates the pseudo dinucleotide composition Descriptor
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
#' @param w the weight factor ranged from 0 to 1. Its default value is 0.05.                                       
#' 
#' @param customprops the users can use their own indices to generate the feature vector. It should be a dict, 
#'                    the key is dinucleotide (string), and its corresponding value is a list type.                       
#' 
#' @return A vector 
#' 
#' @keywords extract PseDNC
#' 
#' @aliases extrDNAPseDNC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNAPseDNC
#' 
#' @seealso See \code{\link{extrDNAPseKNC}} 
#'          
#' @note if the user defined physicochemical indices have not been normalized, it should be normalized.
#'
#' @references 
#' Chen W, Feng P M, Lin H, et al. iRSpot-PseDNC: identify recombination spots with pseudo dinucleotide composition. 
#' \emph{Nucleic acids research}, 2013: gks1450.
#' 
#' @examples
#'  
#' x = 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'
#' extrDNAPseDNC(x)
#' 

extrDNAPseDNC = function (x, lambda = 3, w = 0.05, 
                       normalize = FALSE, customprops = NULL) {
  
  if (checkDNA(x) == FALSE) 
    stop("x has unrecognized  type !")
  
  Didx = data.frame(AA = c(0.06, 0.5, 0.27, 1.59, 0.11, -0.11),
                    AC = c(1.50, 0.50, 0.80, 0.13, 1.29, 1.04),
                    AG = c(0.78, 0.36, 0.09, 0.68, -0.24, -0.62),
                    AT = c(1.07, 0.22, 0.62, -1.02, 2.51, 1.17),
                    CA = c(-1.38, -1.36, -0.27, -0.86, -0.62, -1.25),
                    CC = c(0.06, 1.08, 0.09, 0.56, -0.82, 0.24),
                    CG = c(-1.66, -1.22, -0.44, -0.82, -0.29, -1.39),
                    CT = c(0.78, 0.36, 0.09, 0.68, -0.24, -0.62),
                    GA = c(-0.08, 0.5, 0.27, 0.13, -0.39, 0.71),
                    GC = c(-0.08, 0.22, 1.33, -0.35, 0.65, 1.59),
                    GG = c(0.06, 1.08, 0.09, 0.56, -0.82, 0.24),
                    GT = c(1.50, 0.50, 0.80, 0.13, 1.29, 1.04),
                    TA = c(-1.23, -2.37, -0.44, -2.24, -1.51, -1.39),
                    TC = c(-0.08, 0.5, 0.27, 0.13, -0.39, 0.71),
                    TG = c(-1.38, -1.36, -0.27, -0.86, -0.62, -1.25),
                    TT = c(0.06, 0.5, 0.27, 1.59, 0.11, -0.11))
  
  if (!is.null(customprops))
    Didx = rbind(Didx, customprops)
  
  n = dim(Didx)[1]
  Ddict = make_kmer_index(k = 2) 
  if (normalize) { 
    H = matrix(ncol = 16, nrow = n)
    Didx = as.matrix(Didx)
    for (i in 1:n) H[i, ] = (Didx[i, ] - mean(Didx[i, ]))/(sqrt(sum((Didx[i,] - 
                                      mean(Didx[i, ])) ^ 2) / 16))
    H = round(H, 3)
    colnames(H) = Ddict
  }
  else H = Didx
  
  Theta = vector("list", lambda)
  for (i in 1:lambda) Theta[[i]] = vector("list", n)
  xSplit = strsplit(x, split = "")[[1]]
  xPaste = paste0(xSplit[1:length(xSplit) - 1], xSplit[2:length(xSplit)])
  N = length(xPaste)
  for (i in 1:lambda) {
    temp=c()
    for (j in 1:(N - i)) {
      temp = append(temp, mean((H[, xPaste[j]] - H[, xPaste[j + i]]) ^ 2))
    }
    Theta[[i]] = temp
  }
  theta = sapply(Theta, mean)
  
  fc = summary(factor(xPaste, levels = Ddict), maxsum = 16)
  fc = fc/sum(fc)
  
  Xc1 = fc/(1 + (w * sum(theta)))
  names(Xc1) = paste("Xc1.", names(Xc1), sep = "")
  Xc2 = (w * theta)/(1 + (w * sum(theta)))
  names(Xc2) = paste("Xc2.lambda.", 1:lambda, sep = "")
  Xc = c(Xc1, Xc2)
  
  Xc = round(Xc, 3)
  return (Xc)
}