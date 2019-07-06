#' @title Calculates Atom Additive logP and Molar Refractivity Values Descriptor
#'
#' @description Calculates Atom Additive logP and Molar Refractivity Values Descriptor
#' 
#' @details Calculates ALOGP (Ghose-Crippen LogKow) and the Ghose-Crippen molar 
#' refractivity as described by Ghose, A.K. and Crippen, G.M. 
#' Note the underlying code in CDK assumes that aromaticity 
#' has been detected before evaluating this descriptor. 
#' The code also expects that the molecule will have 
#' hydrogens explicitly set. For SD files, this is 
#' usually not a problem since hydrogens are explicit. 
#' But for the case of molecules obtained from SMILES, 
#' hydrogens must be made explicit.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. This function returns three columns 
#' named \code{ALogP}, \code{ALogp2} and \code{AMR}.
#' 
#' @keywords extrDrugALOGP ALOGP
#'
#' @aliases extrDrugALOGP
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugALOGP
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Ghose, A.K. and Crippen, G.M. , 
#' Atomic physicochemical parameters for three-dimensional structure-directed 
#' quantitative structure-activity relationships. 
#' I. Partition coefficients as a measure of hydrophobicity, 
#' Journal of Computational Chemistry, 1986, 7:565-577.
#' 
#' Ghose, A.K. and Crippen, G.M. , 
#' Atomic physicochemical parameters for three-dimensional-structure-directed 
#' quantitative structure-activity relationships. 
#' 2. Modeling dispersive and hydrophobic interactions, 
#' Journal of Chemical Information and Computer Science, 1987, 27:21-35.
#' 
#' @name property
#' 
#' @examples
#' # Calculates Atom Additive logP and Molar Refractivity Values Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugALOGP(mol)
#' head(dat)
#' 

extrDrugALOGP = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor', 
                  verbose = !silent)

    return(x)

}

#' @rdname property
#' 
#' @title Calculates the Sum of the Atomic Polarizabilities Descriptor
#'
#' @description Calculates the Sum of the Atomic Polarizabilities Descriptor
#'
#' @details Calculates the sum of the atomic polarizabilities 
#' (including implicit hydrogens) descriptor. 
#' Polarizabilities are taken from 
#' \url{https://doi.org/10.1063/1.459444}.
#' 
#' @keywords extrDrugApol Apol
#'
#' @aliases extrDrugApol
#' 
#' @export extrDrugApol
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Sum of the Atomic Polarizabilities Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugApol(mol)
#' head(dat)
#' 

extrDrugApol = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname property
#' 
#' @title Calculates the Descriptor that Describes the Sum of the Absolute 
#' Value of the Difference between Atomic Polarizabilities of 
#' All Bonded Atoms in the Molecule
#'
#' @description Calculates the Descriptor that Describes the Sum of the Absolute 
#' Value of the Difference between Atomic Polarizabilities of 
#' All Bonded Atoms in the Molecule
#'
#' @details This descriptor calculates the sum of the absolute value of the 
#' difference between atomic polarizabilities of all bonded atoms 
#' in the molecule (including implicit hydrogens) with polarizabilities 
#' taken from \url{https://doi.org/10.1063/1.459444}.
#' This descriptor assumes 2-centered bonds.
#' 
#' @keywords extrDrugBPol BPol Polarizability
#'
#' @aliases extrDrugBPol
#' 
#' @export extrDrugBPol
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # the Sum of the Absolute Value of the Difference between Atomic 
#' # Polarizabilities of All Bonded Atoms in the Molecule
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugBPol(mol)
#' head(dat)
#' 

extrDrugBPol = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname property
#' 
#' @title Descriptor that Calculates the Number of Hydrogen Bond Acceptors
#' 
#' @description Descriptor that Calculates the Number of Hydrogen Bond Acceptors
#' 
#' @details This descriptor calculates the number of hydrogen bond acceptors 
#' using a slightly simplified version of the PHACIR atom types. 
#' The following groups are counted as hydrogen bond acceptors:
#' any oxygen where the formal charge of the oxygen is 
#' non-positive (i.e. formal charge <= 0) except
#' \enumerate{
#' \item an aromatic ether oxygen (i.e. an ether oxygen that is 
#' adjacent to at least one aromatic carbon)
#' \item an oxygen that is adjacent to a nitrogen
#' }
#' and any nitrogen where the formal charge of the nitrogen is 
#' non-positive (i.e. formal charge <= 0) except a nitrogen 
#' that is adjacent to an oxygen.
#' 
#' @keywords extrDrugHBondAcceptorCount HBond Acceptor Count
#'
#' @aliases extrDrugHBondAcceptorCount
#' 
#' @export extrDrugHBondAcceptorCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number of Hydrogen Bond Acceptors
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugHBondAcceptorCount(mol)
#' head(dat)
#' 

extrDrugHBondAcceptorCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}
#' @rdname property
#' 
#' @title Descriptor that Calculates the Number of Hydrogen Bond Donors
#' 
#' @description Descriptor that Calculates the Number of Hydrogen Bond Donors
#' 
#' @details This descriptor calculates the number of hydrogen bond donors using 
#' a slightly simplified version of the PHACIR atom types 
#' (\url{https://www.ncbi.nlm.nih.gov/pubmed/12653513}). 
#' The following groups are counted as hydrogen bond donors:
#' \itemize{
#' \item Any-OH where the formal charge of the oxygen is non-negative 
#' (i.e. formal charge >= 0)
#' \item Any-NH where the formal charge of the nitrogen is non-negative 
#' (i.e. formal charge >= 0)
#' }
#' 
#' @keywords extrDrugHBondDonorCount Bond Donor Count
#'
#' @aliases extrDrugHBondDonorCount
#' 
#' @export extrDrugHBondDonorCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number of Hydrogen Bond Donors
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugHBondDonorCount(mol)
#' head(dat)
#' 

extrDrugHBondDonorCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname property
#' 
#' @title Descriptor that Calculates the Number Failures of the Lipinski's Rule Of Five
#'
#' @description Descriptor that Calculates the Number Failures of the Lipinski's Rule Of Five
#'  
#' @details This descriptor calculates the number failures of the Lipinski's Rule Of Five:
#' \url{https://www.ncbi.nlm.nih.gov/pubmed/24981612}.
#' 
#' @keywords extrDrugRuleOfFive Lipinski Rule Five
#'
#' @aliases extrDrugRuleOfFive
#' 
#' @export extrDrugRuleOfFive
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number Failures of the Lipinski's Rule Of Five
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugRuleOfFive(mol)
#' head(dat)
#' 

extrDrugRuleOfFive = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname property
#' 
#' @title Descriptor of Topological Polar Surface Area Based on 
#' Fragment Contributions (TPSA)
#'
#' @description Descriptor of Topological Polar Surface Area Based on 
#' Fragment Contributions (TPSA)
#'
#' @details Calculate the descriptor of topological polar surface area 
#' based on fragment contributions (TPSA).
#' 
#' @keywords extrDrugTPSA Topological Polar Surface Area
#'
#' @aliases extrDrugTPSA
#' 
#' @export extrDrugTPSA
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Ertl, P., Rohde, B., & Selzer, P. (2000). 
#' Fast calculation of molecular polar surface area as a sum of 
#' fragment-based contributions and its application to the prediction 
#' of drug transport properties. 
#' Journal of medicinal chemistry, 43(20), 3714-3717.
#' 
#' @examples
#' # Descriptor of Topological Polar Surface Area Based on 
#' # Fragment Contributions (TPSA)
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugTPSA(mol)
#' head(dat)
#' 

extrDrugTPSA = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname property
#' 
#' @title Descriptor that Calculates the Total Weight of Atoms
#'
#' @description Descriptor that Calculates the Total Weight of Atoms
#' 
#' @details This descriptor calculates the molecular weight.
#' 
#' @keywords extrDrugWeight Weight
#'
#' @aliases extrDrugWeight
#' 
#' @export extrDrugWeight
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Total Weight of Atoms
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugWeight(mol)
#' head(dat)
#' 

extrDrugWeight = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname property
#' 
#' @title Descriptor that Calculates the Prediction of logP 
#' Based on the Atom-Type Method Called XLogP
#'
#' @description Descriptor that Calculates the Prediction of logP 
#' Based on the Atom-Type Method Called XLogP
#' 
#' @details Prediction of logP based on the atom-type method called XLogP.
#' 
#' @keywords extrDrugLogP XLogP
#'
#' @aliases extrDrugLogP
#' 
#' @export extrDrugLogP
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Wang, R., Fu, Y., and Lai, L., 
#' A New Atom-Additive Method for Calculating Partition Coefficients, 
#' Journal of Chemical Information and Computer Sciences, 1997, 37:615-621.
#' 
#' Wang, R., Gao, Y., and Lai, L., 
#' Calculating partition coefficient by atom-additive method, 
#' Perspectives in Drug Discovery and Design, 2000, 19:47-66.
#' 
#' @examples
#' # Calculates the Prediction of logP 
#' # Based on the Atom-Type Method Called XLogP
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugLogP(mol)
#' head(dat)
#' 

extrDrugLogP = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor', 
                verbose = !silent)
  
  return(x)
  
}
