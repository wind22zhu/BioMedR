#' Descriptor that Characterizing Molecular Complexity 
#' in Terms of Carbon Hybridization States
#'
#' Descriptor that Characterizing Molecular Complexity 
#' in Terms of Carbon Hybridization States
#' 
#' This descriptor calculates the fraction of sp3 carbons to sp2 carbons.
#' Note that it only considers carbon atoms and rather than use a simple 
#' ratio it reports the value of Nsp3/(Nsp3 + Nsp2). 
#' The original form of the descriptor (i.e., simple ratio) 
#' has been used to characterize molecular complexity, 
#' especially in the are of natural products, 
#' which usually have a high value of the sp3 to sp2 ratio.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns one column named \code{HybRatio}.
#' 
#' @keywords extrDrugHybridizationRatio Hybridization Ratio
#'
#' @aliases extrDrugHybridizationRatio
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugHybridizationRatio
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugHybridizationRatio(mol)
#' head(dat)
#' 

extrDrugHybridizationRatio = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor', 
                  verbose = !silent)

    return(x)

}
