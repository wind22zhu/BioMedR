#' Calculate the MACCS Molecular Fingerprints (in Compact Format)
#'
#' Calculate the MACCS Molecular Fingerprints (in Compact Format)
#' 
#' The popular 166 bit MACCS keys described by MDL.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A list, each component represents one of the molecules, each element 
#' in the component represents the index of which element in the fingerprint is 1.
#' Each component's name is the length of the fingerprints.
#' 
#' @keywords extrDrugMACCS
#'
#' @aliases extrDrugMACCS
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugMACCS
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugMACCSComplete}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugMACCS(mol)
#' head(fp)
#' 

extrDrugMACCS = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'maccs', verbose = !silent)

        fp = vector('list', 1)
        fp[[1]] = x@bits
        names(fp) = x@nbit

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'maccs', verbose = !silent)

            fp = vector('list', length(molecules))

            for (i in 1:length(molecules)) {

                fp[[i]] = x[[i]]@bits
                names(fp)[i] = x[[i]]@nbit

            }

        }

    return(fp)

}

#' Calculate the MACCS Molecular Fingerprints (in Complete Format)
#'
#' Calculate the MACCS Molecular Fingerprints (in Complete Format)
#' 
#' The popular 166 bit MACCS keys described by MDL.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return An integer vector or a matrix. Each row represents one molecule, 
#' the columns represent the fingerprints.
#' 
#' @keywords extrDrugMACCSComplete
#'
#' @aliases extrDrugMACCSComplete
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugMACCSComplete
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugMACCS}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugMACCSComplete(mol)
#' dim(fp)
#' 

extrDrugMACCSComplete = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'maccs', verbose = !silent)

        fp = integer(166)
        fp[x@bits] = 1L

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'maccs', verbose = !silent)

            fp = matrix(0L, nrow = length(molecules), ncol = 166)

            for (i in 1:length(molecules)) fp[ i, x[[i]]@bits ] = 1L

        }

    return(fp)

}
