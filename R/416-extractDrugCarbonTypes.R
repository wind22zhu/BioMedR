#' @title Topological Descriptor Characterizing the Carbon Connectivity 
#' in Terms of Hybridization
#'
#' @description Topological Descriptor Characterizing the Carbon Connectivity 
#' in Terms of Hybridization
#'
#' @details Calculates the carbon connectivity in terms of hybridization. 
#' The function calculates 9 descriptors in the following order:
#' 
#' \itemize{
#' \item \code{C1SP1} - triply hound carbon bound to one other carbon
#' \item \code{C2SP1} - triply bound carbon bound to two other carbons
#' \item \code{C1SP2} - doubly hound carbon bound to one other carbon
#' \item \code{C2SP2} - doubly bound carbon bound to two other carbons
#' \item \code{C3SP2} - doubly bound carbon bound to three other carbons
#' \item \code{C1SP3} - singly bound carbon bound to one other carbon
#' \item \code{C2SP3} - singly bound carbon bound to two other carbons
#' \item \code{C3SP3} - singly bound carbon bound to three other carbons
#' \item \code{C4SP3} - singly bound carbon bound to four other carbons
#' }
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns 9 columns named 
#' \code{C1SP1}, \code{C2SP1}, \code{C1SP2}, \code{C2SP2}, \code{C3SP2}, 
#' \code{C1SP3}, \code{C2SP3}, \code{C3SP3} and \code{C4SP3}.
#' 
#' @keywords extrDrugCarbonTypes Carbon Types
#'
#' @aliases extrDrugCarbonTypes
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugCarbonTypes
#' 
#' @importFrom rcdk eval.desc
#' 
#' @name topology
#' 
#' @examples
#' # Topological Descriptor Characterizing the Carbon Connectivity 
#' # in Terms of Hybridization
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugCarbonTypes(mol)
#' head(dat)
#' 

extrDrugCarbonTypes = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor', 
                  verbose = !silent)

    return(x)

}

#' @rdname topology
#' 
#' @title Calculates the Eccentric Connectivity Index Descriptor 
#'
#' @description Calculates the Eccentric Connectivity Index Descriptor 
#' 
#' @details Eccentric Connectivity Index (ECI) is a topological descriptor combining 
#' distance and adjacency information. This descriptor is described by Sharma et al.
#' and has been shown to correlate well with a number of physical properties. 
#' The descriptor is also reported to have good discriminatory ability.
#' The eccentric connectivity index for a hydrogen supressed molecular graph 
#' is given by
#' \deqn{x_i^c = \sum_{i = 1}^{n} E(i) V(i)}
#' where E(i) is the eccentricity of the i-th atom 
#' (path length from the i-th atom to the atom farthest from it) 
#' and V(i) is the vertex degree of the i-th atom.
#' 
#' @keywords extrDrugECI Eccentric Connectivity Index
#'
#' @aliases extrDrugECI
#' 
#' @export extrDrugECI
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Sharma, V. and Goswami, R. and Madan, A.K. (1997), 
#' Eccentric Connectivity Index: A Novel Highly Discriminating 
#' Topological Descriptor for Structure-Property and Structure-Activity Studies, 
#' Journal of Chemical Information and Computer Sciences, 37:273-282
#' 
#' @examples
#' # Calculates the Eccentric Connectivity Index Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugECI(mol)
#' head(dat)
#' 

extrDrugECI = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname topology
#' 
#' @title Calculates the FMF Descriptor
#'
#' @description Calculates the FMF Descriptor
#'
#' @details Calculates the FMF descriptor characterizing molecular complexity 
#' in terms of its Murcko framework. This descriptor is the ratio of 
#' heavy atoms in the framework to the total number of heavy atoms 
#' in the molecule. By definition, acyclic molecules which have no frameworks, 
#' will have a value of 0. Note that the authors consider an isolated ring 
#' system to be a framework (even though there is no linker). 
#' 
#' @keywords extrDrugFMF FMF
#'
#' @aliases extrDrugFMF
#' 
#' @export extrDrugFMF
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Yang, Y., Chen, H., Nilsson, I., Muresan, S., & Engkvist, O. (2010). 
#' Investigation of the relationship between topology and selectivity 
#' for druglike molecules. Journal of medicinal chemistry, 
#' 53(21), 7709-7714.
#' 
#' @examples
#' # Calculates the FMF Descriptor
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugFMF(mol)
#' head(dat)
#' 

extrDrugFMF = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname topology
#' 
#' @title Calculate Complexity of a System
#' 
#' @description Calculate Complexity of a System
#' 
#' @details This descriptor calculates the complexity of a system. 
#' The complexity is defined in Nilakantan, R. et al. as:
#' \deqn{C = abs(B^2 - A^2 + A) + \frac{H}{100}}
#' where C is complexity, A is the number of non-hydrogen atoms, 
#' B is the number of bonds and H is the number of heteroatoms.
#' 
#' @keywords extrDrugFragmentComplexity Fragment Complexity
#'
#' @aliases extrDrugFragmentComplexity
#' 
#' @export extrDrugFragmentComplexity
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Nilakantan, R. and Nunn, D.S. and Greenblatt, 
#' L. and Walker, G. and Haraki, K. and Mobilio, D., 
#' A family of ring system-based structural fragments 
#' for use in structure-activity studies: 
#' database mining and recursive partitioning., 
#' Journal of chemical information and modeling, 2006, 46:1069-1077
#' 
#' @examples
#' # Calculate Complexity of a System
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugFragmentComplexity(mol)
#' head(dat)
#' 

extrDrugFragmentComplexity = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname topology
#' 
#' @title Calculate Molecular Distance Edge (MDE) Descriptors for C, N and O
#'
#' @description Calculate Molecular Distance Edge (MDE) Descriptors for C, N and O
#'
#' @details This descriptor calculates the 10 molecular distance edge (MDE) descriptor 
#' described in Liu, S., Cao, C., & Li, Z, and in addition it calculates 
#' variants where O and N are considered.
#' 
#' @keywords extrDrugMDE MDE Molecular Distance Edge
#'
#' @aliases extrDrugMDE
#' 
#' @export extrDrugMDE
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Liu, S., Cao, C., & Li, Z. (1998). 
#' Approach to estimation and prediction for normal boiling point (NBP) 
#' of alkanes based on a novel molecular distance-edge (MDE) vector, lambda. 
#' Journal of chemical information and computer sciences, 38(3), 387-394.
#' 
#' @examples
#' # Calculate Molecular Distance Edge (MDE) Descriptors for C, N and O
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugMDE(mol)
#' head(dat)
#' 

extrDrugMDE = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname topology
#' 
#' @title Descriptor that Calculates the Petitjean Number of a Molecule
#'
#' @description Descriptor that Calculates the Petitjean Number of a Molecule
#' 
#' @details This descriptor calculates the Petitjean number of a molecule.
#' According to the Petitjean definition, the eccentricity of a vertex 
#' corresponds to the distance from that vertex to the most remote vertex 
#' in the graph. 
#' 
#' The distance is obtained from the distance matrix as the count of edges 
#' between the two vertices. If \code{r(i)} is the largest matrix entry 
#' in row \code{i} of the distance matrix \code{D}, then the radius is defined 
#' as the smallest of the \code{r(i)}. The graph diameter \code{D} is defined as 
#' the largest vertex eccentricity in the graph. 
#' (\url{http://www.edusoft-lc.com/molconn/manuals/400/chaptwo.html})
#' 
#' @keywords extrDrugPetitjeanNumber Petitjean
#'
#' @aliases extrDrugPetitjeanNumber
#' 
#' @export extrDrugPetitjeanNumber
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Petitjean Number of a Molecule
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugPetitjeanNumber(mol)
#' head(dat)
#' 

extrDrugPetitjeanNumber = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname topology
#' 
#' @title Descriptor that Calculates the Petitjean Shape Indices
#' 
#' @description Descriptor that Calculates the Petitjean Shape Indices
#' 
#' @details The topological and geometric shape indices described Petitjean 
#' and Bath et al. respectively. Both measure the anisotropy in a molecule.
#' 
#' @keywords extrDrugPetitjeanShapeIndex Petitjean Geometric Shape Index
#'
#' @aliases extrDrugPetitjeanShapeIndex
#' 
#' @export extrDrugPetitjeanShapeIndex
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Petitjean, M., 
#' Applications of the radius-diameter diagram to the classification of 
#' topological and geometrical shapes of chemical compounds, 
#' Journal of Chemical Information and Computer Science, 
#' 1992, 32:331-337
#' 
#' Bath, P.A. and Poirette, A.R. and Willet, P. and Allen, F.H. , 
#' The Extent of the Relationship between the Graph-Theoretical 
#' and the Geometrical Shape Coefficients of Chemical Compounds, 
#' Journal of Chemical Information and Computer Science, 1995, 35:714-716.
#' 
#' @examples
#' # Calculates the Petitjean Shape Indices
#' sdf = system.file('sysdata/test.sdf', package = 'BioMedR')
#' mol = readMolFromSDF(sdf)
#' dat = extrDrugPetitjeanShapeIndex(mol)
#' head(dat)
#' 

extrDrugPetitjeanShapeIndex = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname topology
#' 
#' @title Descriptor that Calculates the Volume of A Molecule
#'
#' @description Descriptor that Calculates the Volume of A Molecule
#'  
#' @details This descriptor calculates the volume of a molecule.
#' 
#' @keywords extrDrugVABC Volume VABC
#'
#' @aliases extrDrugVABC
#' 
#' @export extrDrugVABC
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Volume of A Molecule
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugVABC(mol)
#' head(dat)
#' 

extrDrugVABC = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname topology
#' 
#' @title Descriptor that Calculates the Vertex Adjacency Information of A Molecule
#'
#' @description Descriptor that Calculates the Vertex Adjacency Information of A Molecule
#' 
#' @details Vertex adjacency information (magnitude): 
#' \eqn{1 + \log_2^m} where \eqn{m} is the number of heavy-heavy bonds. 
#' If \eqn{m} is zero, then \code{0} is returned.
#' 
#' @keywords extrDrugVAdjMa Vertex Adjacency Magnitude
#'
#' @aliases extrDrugVAdjMa
#' 
#' @export extrDrugVAdjMa
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Vertex Adjacency Information of A Molecule
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugVAdjMa(mol)
#' head(dat)
#' 

extrDrugVAdjMa = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname topology
#' 
#' @title Descriptor that Calculates the Weighted Path (Molecular ID)
#'
#' @description Descriptor that Calculates the Weighted Path (Molecular ID)
#' 
#' @details This descriptor calculates the weighted path (molecular ID) 
#' described by Randic, characterizing molecular branching. 
#' Five descriptors are calculated, based on the implementation in the ADAPT 
#' software package. Note that the descriptor is based on identifying all paths 
#' between pairs of atoms and so is NP-hard. 
#' This means that it can take some time for large, complex molecules.
#'
#' @return WTPT 
#' \code{WTPT.1}, \code{WTPT.2}, \code{WTPT.3}, \code{WTPT.4}, \code{WTPT.5}:
#' \itemize{
#' \item \code{WTPT.1} - molecular ID
#' \item \code{WTPT.2} - molecular ID / number of atoms
#' \item \code{WTPT.3} - sum of path lengths starting from heteroatoms
#' \item \code{WTPT.4} - sum of path lengths starting from oxygens
#' \item \code{WTPT.5} - sum of path lengths starting from nitrogens
#' }
#' 
#' @keywords extrDrugWeightedPath Weighted Path
#'
#' @aliases extrDrugWeightedPath
#' 
#' @export extrDrugWeightedPath
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Randic, M., On molecular identification numbers (1984).
#' Journal of Chemical Information and Computer Science, 24:164-175.
#' 
#' @examples
#' # Calculates the Weighted Path (Molecular ID)
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugWeightedPath(mol)
#' head(dat)
#' 

extrDrugWeightedPath = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor', 
                verbose = !silent)
  
  return(x)
  
}


#' @rdname topology
#' 
#' @title Descriptor that Calculates Wiener Path Number and Wiener Polarity Number
#'
#' @description Descriptor that Calculates Wiener Path Number and Wiener Polarity Number
#' 
#' @details This descriptor calculates the Wiener numbers, including the 
#' Wiener Path number and the Wiener Polarity Number. Wiener path number: 
#' half the sum of all the distance matrix entries; Wiener polarity number: 
#' half the sum of all the distance matrix entries with a value of 3. 
#' 
#' @keywords extrDrugWienerNumbers Wiener Numbers
#'
#' @aliases extrDrugWienerNumbers
#' 
#' @export extrDrugWienerNumbers
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Wiener, H. (1947). 
#' Structural determination of paraffin boiling points. 
#' Journal of the American Chemical Society, 69(1), 17-20.
#' 
#' @examples
#' # Calculates Wiener Path Number and Wiener Polarity Number
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugWienerNumbers(mol)
#' head(dat)
#' 

extrDrugWienerNumbers = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname topology
#' 
#' @title Descriptor that Calculates the Sum of the Squared Atom Degrees 
#' of All Heavy Atoms
#'
#' @description Descriptor that Calculates the Sum of the Squared Atom Degrees 
#' of All Heavy Atoms
#' 
#' @details Zagreb index: the sum of the squares of atom degree over 
#' all heavy atoms \code{i}.
#' 
#' @keywords extrDrugZagrebIndex Zagreb
#'
#' @aliases extrDrugZagrebIndex
#' 
#' @export extrDrugZagrebIndex
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Calculates the Sum of the Squared Atom Degrees 
#' # of All Heavy Atoms
#' smi = system.file('vignettedata/test.smi', package = 'BioMedR')
#' mol = readMolFromSmi(smi, type = 'mol')
#' dat = extrDrugZagrebIndex(mol)
#' head(dat)
#' 

extrDrugZagrebIndex = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

