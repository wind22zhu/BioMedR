#' @title Calculates the Moreau-Broto Autocorrelation Descriptors using Partial Charges
#'
#' @description Calculates the Moreau-Broto Autocorrelation Descriptors using Partial Charges
#'
#' @details Calculates the ATS autocorrelation descriptor, 
#' where the weight equal to the charges.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns 5 columns named 
#' \code{ATSc1}, \code{ATSc2}, \code{ATSc3}, \code{ATSc4}, \code{ATSc5}.
#' 
#' @keywords extrDrugAutocorrelationcharge Autocorrelation Charge
#'
#' @aliases extrDrugAutocorrelationcharge
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugAutocorrelationcharge
#' 
#' @importFrom rcdk eval.desc
#' 
#' @name Autocorrelation
#' 
#' @examples
#' # Calculates the Moreau-Broto Autocorrelation Descriptors using Partial Charges
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugAutocorrelationcharge(mol)
#' head(dat)
#' 

extrDrugAutocorrelationcharge = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge', 
                  verbose = !silent)

    return(x)

}

#' @rdname Autocorrelation
#' 
#' @title Calculates the Moreau-Broto Autocorrelation Descriptors using Atomic Weight
#'
#' @description Calculates the Moreau-Broto Autocorrelation Descriptors using Atomic Weight
#'
#' @details Calculates the ATS autocorrelation descriptor, 
#' where the weight equal to the scaled atomic mass.
#'
#' @return extrDrugAutocorrelationMass: This function returns 5 columns named 
#' \code{ATSm1}, \code{ATSm2}, \code{ATSm3}, \code{ATSm4}, \code{ATSm5}.
#' 
#' @keywords extrDrugAutocorrelationMass Autocorrelation Mass
#'
#' @aliases extrDrugAutocorrelationMass
#' 
#' @export extrDrugAutocorrelationMass
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Moreau, Gilles, and Pierre Broto. 
#' The autocorrelation of a topological structure: a new molecular descriptor.
#' Nouv. J. Chim 4 (1980): 359-360.
#' 
#' @examples
#' # Calculates the Moreau-Broto Autocorrelation Descriptors using Atomic Weight
#' dat = extrDrugAutocorrelationMass(mol)
#' head(dat)
#' 

extrDrugAutocorrelationMass = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname Autocorrelation
#' 
#' @title Calculates the Moreau-Broto Autocorrelation Descriptors using Polarizability
#'
#' @description Calculates the Moreau-Broto Autocorrelation Descriptors using Polarizability
#'
#' @details Calculates the ATS autocorrelation descriptor using polarizability.
#'
#' @return extrDrugAutocorrelationPolarizability: This function returns 5 columns named 
#' \code{ATSp1}, \code{ATSp2}, \code{ATSp3}, \code{ATSp4}, \code{ATSp5}.
#' 
#' @keywords extrDrugAutocorrelationPolarizability 
#' Autocorrelation Polarizability
#' 
#' @aliases extrDrugAutocorrelationPolarizability
#' 
#' @export extrDrugAutocorrelationPolarizability
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Moreau-Broto Autocorrelation Descriptors using Polarizability
#' dat = extrDrugAutocorrelationPolarizability(mol)
#' head(dat)
#' 

extrDrugAutocorrelationPolarizability = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability', 
                verbose = !silent)
  
  return(x)
  
}


