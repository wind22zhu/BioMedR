#' @title Calculates the Number of Amino Acids Descriptor
#'
#' @description Calculates the Number of Amino Acids Descriptor
#'
#' @details Calculates the number of each amino acids (total 20 types) 
#' found in the molecues.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns 20 columns named 
#' \code{nA}, \code{nR}, \code{nN}, \code{nD}, \code{nC},
#' \code{nF}, \code{nQ}, \code{nE}, \code{nG}, \code{nH}, 
#' \code{nI}, \code{nP}, \code{nL} \code{nK}, \code{nM}, 
#' \code{nS}, \code{nT}, \code{nY} \code{nV}, \code{nW}.
#' 
#' @keywords extrDrugAminoAcidCount Amino Acid Count
#'
#' @aliases extrDrugAminoAcidCount
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugAminoAcidCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @name Constitutional 
#' 
#' @examples
#' # Calculates the Number of Amino Acids Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugAminoAcidCount(mol)
#' head(dat)
#' 

extrDrugAminoAcidCount = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor', 
                  verbose = !silent)

    return(x)

}

#' @rdname Constitutional 
#' 
#' @title Calculates the Number of Aromatic Atoms Descriptor
#'
#' @description Calculates the Number of Aromatic Atoms Descriptor
#'
#' @details Calculates the number of aromatic atoms of a molecule.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns one column named \code{naAromAtom}.
#' 
#' @keywords extrDrugAromaticAtomsCount Aromatic Atoms Count
#'
#' @aliases extrDrugAromaticAtomsCount
#' 
#' @export extrDrugAromaticAtomsCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number of Aromatic Atoms Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugAromaticAtomsCount(mol)
#' head(dat)
#' 

extrDrugAromaticAtomsCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname Constitutional
#' 
#' @title Calculates the Number of Aromatic Bonds Descriptor
#'
#' @description Calculates the Number of Aromatic Bonds Descriptor
#'
#' @details Calculates the number of aromatic bonds of a molecule.
#' 
#' @keywords extrDrugAromaticBondsCount Aromatic Bond Count
#'
#' @aliases extrDrugAromaticBondsCount
#' 
#' @export extrDrugAromaticBondsCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number of Aromatic Bonds Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugAromaticBondsCount(mol)
#' head(dat)
#' 

extrDrugAromaticBondsCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname Constitutional
#' 
#' @title Calculates the Number of Atom Descriptor
#'
#' @description Calculates the Number of Atom Descriptor
#'
#' @details Calculates the number of atoms of a certain element type in a molecule. 
#' By default it returns the count of all atoms.
#' 
#' @keywords extrDrugAtomCount Atom Count
#'
#' @aliases extrDrugAtomCount
#' 
#' @export extrDrugAtomCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Number of Atom Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugAtomCount(mol)
#' head(dat)
#' 

extrDrugAtomCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname Constitutional
#' 
#' @title Calculates the Descriptor Based on the Number of Bonds of a 
#' Certain Bond Order
#'
#' @description Calculates the Descriptor Based on the Number of Bonds of a 
#' Certain Bond Order
#' 
#' @details Calculates the descriptor based on the number of bonds of a 
#' certain bond order.
#' 
#' @keywords extrDrugBondCount Bond Count
#'
#' @aliases extrDrugBondCount
#' 
#' @export extrDrugBondCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Descriptor Based on the Number of Bonds of a 
#' # Certain Bond Order
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugBondCount(mol)
#' head(dat)
#' 

extrDrugBondCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname Constitutional
#' 
#' @title Descriptor that Calculates the Number of Atoms in the Largest Chain
#' 
#' @description Descriptor that Calculates the Number of Atoms in the Largest Chain
#' 
#' @details This descriptor calculates the number of atoms in the largest chain. 
#' Note that a chain exists if there are two or more atoms. 
#' Thus single atom molecules will return \code{0}.
#' 
#' @keywords extrDrugLargestChain Largest Chain
#'
#' @aliases extrDrugLargestChain
#' 
#' @export extrDrugLargestChain
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Descriptor that Calculates the Number of Atoms in the Largest Chain
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugLargestChain(mol)
#' head(dat)
#' 

extrDrugLargestChain = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname Constitutional
#' 
#' @title Descriptor that Calculates the Number of Atoms in the Largest Pi Chain
#'
#' @description Descriptor that Calculates the Number of Atoms in the Largest Pi Chain
#' 
#' @details This descriptor calculates the number of atoms in the largest pi chain.
#' 
#' @keywords extrDrugLargestPiSystem Largest Pi Chain
#'
#' @aliases extrDrugLargestPiSystem
#' 
#' @export extrDrugLargestPiSystem
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Descriptor that Calculates the Number of Atoms in the Largest Pi Chain
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugLargestPiSystem(mol)
#' head(dat)
#' 

extrDrugLargestPiSystem = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname Constitutional
#' 
#' @title Descriptor that Calculates the Number of Atoms in the Longest Aliphatic Chain
#'
#' @description Descriptor that Calculates the Number of Atoms in the Longest Aliphatic Chain
#' 
#' @details This descriptor calculates the number of atoms in the longest aliphatic chain.
#' 
#' @keywords extrDrugLongestAliphaticChain Longest Aliphatic Chain
#'
#' @aliases extrDrugLongestAliphaticChain
#' 
#' @export extrDrugLongestAliphaticChain
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Descriptor that Calculates the Number of Atoms in the Longest Aliphatic Chain
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugLongestAliphaticChain(mol)
#' head(dat)
#' 

extrDrugLongestAliphaticChain = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname Constitutional
#' 
#' @title Descriptor that Calculates the Number of Nonrotatable Bonds on A Molecule
#' 
#' @description Descriptor that Calculates the Number of Nonrotatable Bonds on A Molecule
#' 
#' @details The number of rotatable bonds is given by the SMARTS specified by 
#' Daylight on SMARTS tutorial
#' (\url{http://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html#EXMPL})
#' 
#' @keywords extrDrugRotatableBondsCount Rotatable Bonds Count
#'
#' @aliases extrDrugRotatableBondsCount
#' 
#' @export extrDrugRotatableBondsCount
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Descriptor that Calculates the Number of Nonrotatable Bonds on A Molecule
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugRotatableBondsCount(mol)
#' head(dat)
#' 

extrDrugRotatableBondsCount = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

