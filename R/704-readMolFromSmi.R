#' Read Molecules from SMILES Files and Return Parsed Java 
#' Molecular Object or Plain Text List
#'
#' Read Molecules from SMILES Files and Return Parsed Java 
#' Molecular Object or Plain Text List
#' 
#' This function reads molecules from SMILES strings and return 
#' parsed Java molecular object or plain text list 
#' needed by \code{extrDrug...()} functions.
#' 
#' @param smifile Character vector, containing SMILES file location(s).
#' @param type \code{'mol'} or \code{'text'}. \code{'mol'} returns parsed 
#'        Java molecular object, used for \code{'text'} returns (plain-text) 
#'        character string list. For common molecular descriptors and 
#'        fingerprints, use \code{'mol'}. 
#'        For descriptors and fingerprints calculated by OpenBabel,
#'        i.e. functions named \code{extrDrugOB...()} , use \code{'text'}.
#' 
#' @return A list, containing parsed Java molecular object or character strings.
#' 
#' @keywords readMolFromSmi MOL SMILES
#'
#' @aliases readMolFromSmi
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{readMolFromSDF}} for reading molecules 
#' from SDF files and returning parsed Java molecular object.
#' 
#' @export readMolFromSmi
#' 
#' @importFrom rcdk parse.smiles
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol1 = readMolFromSmi(smi, type = 'mol')
#' mol2 = readMolFromSmi(smi, type = 'text')
#' 

readMolFromSmi = function (smifile, type = c('mol', 'text')) {

    if (type == 'mol') {

        txt = scan(smifile, what = 'complex', quiet = TRUE)
        smi = as.character(txt)
        mol = parse.smiles(smi)

        } else if (type == 'text') {

            mol = scan(smifile, what = 'complex', quiet = TRUE)

            } else {

                stop('type must be one of "mol" or "text"')

            }

    return(mol)

}
