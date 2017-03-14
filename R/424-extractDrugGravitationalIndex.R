#' @title Descriptor Characterizing the Mass Distribution of the Molecule.
#'
#' @description Descriptor Characterizing the Mass Distribution of the Molecule.
#' 
#' @details Descriptor characterizing the mass distribution of the molecule described by 
#' Katritzky et al. For modelling purposes the value of the descriptor is 
#' calculated both with and without H atoms. 
#' Furthermore the square and cube roots of the descriptor 
#' are also generated as described by Wessel et al.
#' 
#' @param molecules Parsed molucule object.
#' @param silent Logical. Whether the calculating process 
#' should be shown or not, default is \code{TRUE}.
#'
#' @return A data frame, each row represents one of the molecules, 
#' each column represents one feature. 
#' This function returns 9 columns:
#' \itemize{
#' \item \code{GRAV.1} - gravitational index of heavy atoms
#' \item \code{GRAV.2} - square root of gravitational index of heavy atoms
#' \item \code{GRAV.3} - cube root of gravitational index of heavy atoms
#' \item \code{GRAVH.1} - gravitational index - hydrogens included
#' \item \code{GRAVH.2} - square root of hydrogen-included gravitational index
#' \item \code{GRAVH.3} - cube root of hydrogen-included gravitational index
#' \item \code{GRAV.4} - grav1 for all pairs of atoms (not just bonded pairs)
#' \item \code{GRAV.5} - grav2 for all pairs of atoms (not just bonded pairs)
#' \item \code{GRAV.6} - grav3 for all pairs of atoms (not just bonded pairs)
#' }
#' 
#' @keywords extrDrugGravitationalIndex Gravitational Index
#'
#' @aliases extrDrugGravitationalIndex
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>, 
#'         Nan Xiao <\url{http://r2s.name}>
#' 
#' @export extrDrugGravitationalIndex
#' 
#' @importFrom rcdk eval.desc
#' 
#' @references
#' Katritzky, A.R. and Mu, L. and Lobanov, V.S. and Karelson, M., 
#' Correlation of Boiling Points With Molecular Structure. 
#' 1. A Training Set of 298 Diverse Organics and a 
#' Test Set of 9 Simple Inorganics, 
#' J. Phys. Chem., 1996, 100:10400-10407.
#' 
#' Wessel, M.D. and Jurs, P.C. and Tolan, J.W. and Muskal, S.M. , 
#' Prediction of Human Intestinal Absorption of Drug Compounds 
#' From Molecular Structure, 
#' Journal of Chemical Information and Computer Sciences, 1998, 38:726-735.
#' 
#' @name geometric
#' 
#' @examples
#' sdf = system.file('sysdata/test.sdf', package = 'BioMedR')
#' mol = readMolFromSDF(sdf)
#' # Descriptor Characterizing the Mass Distribution of the Molecule
#' dat = extrDrugGravitationalIndex(mol)
#' head(dat)
#' 

extrDrugGravitationalIndex = function (molecules, silent = TRUE) {

    x = eval.desc(molecules, 
                  'org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor', 
                  verbose = !silent)

    return(x)

}

#' @rdname geometric
#' 
#' @title Calculates the Ratio of Length to Breadth Descriptor
#' 
#' @description Calculates the Ratio of Length to Breadth Descriptor
#' 
#' @details Calculates the Ratio of Length to Breadth, as a result ti does not perform 
#' any orientation and only considers the X & Y 
#' extents for a series of rotations about the Z axis 
#' (in 10 degree increments).
#'
#' @return extrDrugLengthOverBreadth: 
#' This function returns two columns named \code{LOBMAX} and \code{LOBMIN}:
#' \itemize{
#' \item \code{LOBMAX} - The maximum L/B ratio;
#' \item \code{LOBMIN} - The L/B ratio for the rotation that results in the 
#' minimum area (defined by the product of the X & Y extents for that orientation).
#' }
#' 
#' @keywords extrDrugLengthOverBreadth Length Breadth
#'
#' @aliases extrDrugLengthOverBreadth
#' 
#' @export extrDrugLengthOverBreadth
#' 
#' @importFrom rcdk eval.desc
#' 
#' @note extrDrugLengthOverBreadth : The descriptor assumes that the 
#' atoms have been configured.
#' 
#' @examples
#' # Calculates the Ratio of Length to Breadth Descriptor
#' dat = extrDrugLengthOverBreadth(mol)
#' head(dat)
#' 

extrDrugLengthOverBreadth = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

#' @rdname geometric
#' 
#' @title Descriptor that Calculates the Principal Moments of 
#' Inertia and Ratios of the Principal Moments
#'
#' @description Descriptor that Calculates the Principal Moments of 
#' Inertia and Ratios of the Principal Moments
#' 
#' @details A descriptor that calculates the moment of inertia and radius of gyration. 
#' Moment of inertia (MI) values characterize the mass distribution of a molecule. 
#' Related to the MI values, ratios of the MI values along the three principal 
#' axes are also well know modeling variables. 
#' This descriptor calculates the MI values along the 
#' X, Y and Z axes as well as the ratio's X/Y, X/Z and Y/Z. 
#' Finally it also calculates the radius of 
#' gyration of the molecule.
#'
#' @return extrDrugMomentOfInertia:  
#' This function returns 7 columns named 
#' \code{MOMI.X}, \code{MOMI.Y}, \code{MOMI.Z}, 
#' \code{MOMI.XY}, \code{MOMI.XZ}, \code{MOMI.YZ}, \code{MOMI.R}:
#' \itemize{
#' \item \code{MOMI.X} - MI along X axis
#' \item \code{MOMI.Y} - MI along Y axis
#' \item \code{MOMI.Z} - MI along Z axis
#' \item \code{MOMI.XY} - X/Y
#' \item \code{MOMI.XZ} - X/Z
#' \item \code{MOMI.YZ} - Y/Z
#' \item \code{MOMI.R} - Radius of gyration
#' }
#' One important aspect of the algorithm is that if the eigenvalues 
#' of the MI tensor are below \code{1e-3}, 
#' then the ratio's are set to a default of 1000.
#' 
#' @keywords extrDrugMomentOfInertia Moment Inertia
#' 
#' @aliases extrDrugMomentOfInertia
#' 
#' @export extrDrugMomentOfInertia
#' 
#' @importFrom rcdk eval.desc
#' 
#' @examples
#' # Descriptor that Calculates the Principal Moments of 
#' # Inertia and Ratios of the Principal Moments
#' dat = extrDrugMomentOfInertia(mol)
#' head(dat)
#' 

extrDrugMomentOfInertia = function (molecules, silent = TRUE) {
  
  x = eval.desc(molecules, 
                'org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor', 
                verbose = !silent)
  
  return(x)
  
}

