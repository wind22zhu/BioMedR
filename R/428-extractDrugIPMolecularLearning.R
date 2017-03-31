#' Calculates the Descriptor that Evaluates the Ionization Potential
#'
#' Calculates the Descriptor that Evaluates the Ionization Potential
#'
#' Calculate the ionization potential of a molecule. 
#' The descriptor assumes that explicit hydrogens have been added to the molecules.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns one column named \code{MolIP}.
#' 
#' @keywords extrDrugIPMolecularLearning Ionization Potential
#'
#' @aliases extrDrugIPMolecularLearning
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugIPMolecularLearning
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugIPMolecularLearning(mol)
#' head(dat)
#' 

extrDrugIPMolecularLearning = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.IPMolecularLearningDescriptor', 
                  verbose = !silent)

    return(x)

}
