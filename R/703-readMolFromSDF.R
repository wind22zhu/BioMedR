#' Read Molecules from SDF Files and Return Parsed Java Molecular Object
#'
#' Read Molecules from SDF Files and Return Parsed Java Molecular Object
#' 
#' This function reads molecules from SDF files and return 
#' parsed Java molecular object needed by \code{extrDrug...} functions.
#' 
#' @param sdffile Character vector, containing SDF file location(s).
#' 
#' @return A list, containing parsed Java molecular object.
#' 
#' @keywords readMolFromSDF MOL SDF
#'
#' @aliases readMolFromSDF
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @seealso See \code{\link{readMolFromSmi}} for reading molecules by SMILES 
#' string and returning parsed Java molecular object.
#' 
#' @export readMolFromSDF
#' 
#' @importFrom rcdk load.molecules
#' 
#' @examples
#' mol  = readMolFromSDF(system.file('compseq/DB00859.sdf', package = 'BioMedR'))
#' mols = readMolFromSDF(c(system.file('compseq/DB00859.sdf', package = 'BioMedR'), 
#'                         system.file('compseq/DB00860.sdf', package = 'BioMedR')))
#' 

readMolFromSDF = function (sdffile) {

    mol = load.molecules(sdffile)

    return(mol)

}
