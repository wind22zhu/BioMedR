#' Calculate the PubChem Molecular Fingerprints (in Compact Format)
#'
#' Calculate the PubChem Molecular Fingerprints (in Compact Format)
#' 
#' Calculate the 881 bit fingerprints defined by PubChem.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A list, each component represents one of the molecules, each element 
#' in the component represents the index of which element in the fingerprint is 1.
#' Each component's name is the length of the fingerprints.
#' 
#' @keywords extrDrugPubChem
#'
#' @aliases extrDrugPubChem
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugPubChem
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugPubChemComplete}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugPubChem(mol)
#' head(fp)
#' 

extrDrugPubChem = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'pubchem', verbose = !silent)

        fp = vector('list', 1)
        fp[[1]] = x@bits
        names(fp) = x@nbit

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'pubchem', verbose = !silent)

            fp = vector('list', length(molecules))

            for (i in 1:length(molecules)) {

                fp[[i]] = x[[i]]@bits
                names(fp)[i] = x[[i]]@nbit

            }

        }

    return(fp)

}

#' Calculate the PubChem Molecular Fingerprints (in Complete Format)
#'
#' Calculate the PubChem Molecular Fingerprints (in Complete Format)
#' 
#' Calculate the 881 bit fingerprints defined by PubChem.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return An integer vector or a matrix. Each row represents one molecule, 
#' the columns represent the fingerprints.
#' 
#' @keywords extrDrugPubChemComplete
#'
#' @aliases extrDrugPubChemComplete
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugPubChemComplete
#' 
#' @importFrom rcdk get.fingerprint
#' 
#' @seealso \link{extrDrugPubChem}
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' fp  = extrDrugPubChemComplete(mol)
#' dim(fp)
#' 

extrDrugPubChemComplete = function (molecules, silent = TRUE) {

    if (length(molecules) == 1) {

        x = get.fingerprint(molecules, type = 'pubchem', verbose = !silent)

        fp = integer(881)
        fp[x@bits] = 1L

        } else {

            x = lapply(molecules, get.fingerprint, 
                       type = 'pubchem', verbose = !silent)

            fp = matrix(0L, nrow = length(molecules), ncol = 881)

            for (i in 1:length(molecules)) fp[ i, x[[i]]@bits ] = 1L

        }

    return(fp)

}
