#' The Dinucleotide-based Auto Covariance Descriptor
#' 
#' The Dinucleotide-based Auto Covariance Descriptor
#'
#' This function calculates the dinucleotide-based auto covariance descriptor
#' 
#' @param x the input data, which should be a list or file type.
#' 
#' @param index the physicochemical indices, it should be a list and there are 38
#'              different physicochemical indices (Table 1), which the users can choose.
#'
#' @param nlag an integer larger than or equal to 0 and less than or equal to L-2 (L means the length 
#'             of the shortest DNA sequence in the dataset). It represents the distance between two dinucleotides.
#'
#' @param normaliztion with this option, the final feature vector will be normalized based
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
#' @keywords extract DAC
#' 
#' @aliases extrDNADAC
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @export extrDNADAC
#' 
#' @seealso See \code{\link{extrDNADCC}} and \code{\link{extrDNADACC}}
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
#' extrDNADAC(x)


extrDNADAC = function (x, index = c('Twist', 'Tilt'), nlag = 2, normaliztion = FALSE, 
                    customprops = NULL, allprop = FALSE) {
  
  if (!checkDNA(x)) 
    stop ('x has character except A, G, C, T' )
  
   Didx = read.csv(system.file('sysdata/38_dinucleotide_physicochemical_indices.csv', package = 'BioMedR'), header = TRUE, sep = ";")
   didx = Didx[, -1]
   row.names(didx) = Didx[, 1]
  
   if (allprop) index = Didx[,1]
   if (!is.null(customprops)) {
     didx = rbind(didx, customprops)
     index = c(as.character(index), rownames(customprops))
   } 
  
  n = length(index)
  
  #### normaliztion ####
  Ddict = make_kmer_index(k = 2)
  if (normaliztion) {
    indmean = rowMeans(didx[index, ])
    indsd   = apply(didx[index, ], 1, sd) 
    
    Pr = data.frame(matrix(ncol = 16, nrow = n))
    for (i in 1:n) Pr[i, ] = (didx[index[i], ] - indmean[i])/indsd[i]
    names(Pr) = Ddict 
  }
  else Pr = didx[index, ]
 
  xSplit = strsplit(x, split = "")[[1]]
  xPaste = paste0(xSplit[1:length(xSplit) - 1], xSplit[2:length(xSplit)])
  P = vector('list', n)
  for (i in 1:n) P[[i]] = xPaste
  for (i in 1:n) {
    for (j in Ddict) {
     try(P[[i]][which(P[[i]] == j)] <- Pr[i, j], silent = TRUE)
     }
  }
  P = lapply(P, as.numeric)

  DAC = vector('list', n)
  N  = length(xPaste)



  for (i in 1:n) {
     for (j in 1:nlag) {
       Pbar = lapply(seq_along(P), function(i) P[[i]][1:(N-j)])
       Pbar = sapply(Pbar, sum)/nchar(x)
       DAC[[i]][j] = (sum((P[[i]][1:(N - j)] - Pbar[i]) * (P[[i]][(1:(N - j)) + j] - Pbar[i])))/(N-j)
     }
   }

   DAC = unlist(DAC)

   names(DAC) = as.vector(t(outer(index, paste('.lag', 1:nlag, sep = ''), 
                               paste, sep = '')))
   
   DAC = round(DAC, 3)
   return(DAC)


}