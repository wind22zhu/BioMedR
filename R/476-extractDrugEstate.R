#' Calculate the E-State Molecular Fingerprints (in Compact Format)
#'
#' Calculate the E-State Molecular Fingerprints (in Compact Format)
#' 
#' 79 bit fingerprints corresponding to the E-State atom types 
#' described by Hall and Kier.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A list, each component represents one of the molecules, each element 
#' in the component represents the index of which element in the fingerprint is 1.
#' Each component's name is the length of the fingerprints.
#' 
#' @keywords extrDrugEstate
#'
#' @aliases extrDrugEstate
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugEstate
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugEstateComplete}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugEstate(mol)
#' head(fp)
#' 

extrDrugEstate = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'estate', verbose = !silent)

        fp = vector('list', 1)
        fp[[1]] = x@bits
        names(fp) = x@nbit

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'estate', verbose = !silent)

            fp = vector('list', length(molecules))

            for (i in 1:length(molecules)) {

                fp[[i]] = x[[i]]@bits
                names(fp)[i] = x[[i]]@nbit

            }

        }

    return(fp)

}

#' Calculate the E-State Molecular Fingerprints (in Complete Format)
#'
#' Calculate the E-State Molecular Fingerprints (in Complete Format)
#' 
#' 79 bit fingerprints corresponding to the E-State atom types 
#' described by Hall and Kier.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return An integer vector or a matrix. Each row represents one molecule, 
#' the columns represent the fingerprints.
#' 
#' @keywords extrDrugEstateComplete
#'
#' @aliases extrDrugEstateComplete
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugEstateComplete
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugEstate}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugEstateComplete(mol)
#' dim(fp)
#' 

extrDrugEstateComplete = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'estate', verbose = !silent)

        fp = integer(79)
        fp[x@bits] = 1L

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'estate', verbose = !silent)

            fp = matrix(0L, nrow = length(molecules), ncol = 79)

            for (i in 1:length(molecules)) fp[ i, x[[i]]@bits ] = 1L

        }

    return(fp)

}
