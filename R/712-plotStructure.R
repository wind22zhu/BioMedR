#' Plots compound structure(s) for molecules stored in SDF and SDFset containers
#' 
#' Plots compound structure(s) for molecules stored in SDF and SDFset containers.
#' 
#' The function \code{plotStructure} depicts a single 2D compound structure based on 
#' the XY-coordinates specified in the atom block of an \code{SDF}. The functions depend on 
#' the availability of the XY-coordinates in the source SD file and only 2D (not 3D) 
#' representations are plotted correctly.
#' 
#' @param sdf Object of class \code{SDF}
#' 
#' @param atomcex Font size for atom labels
#' 
#' @param atomnum If \code{TRUE}, then the atom numbers are included in the plot. 
#'                They are the position numbers of each atom in the atom 
#'                block of an \code{SDF}.
#'                
#' @param no_print_atoms  Excludes specified atoms from being plotted.
#'
#' @param noHbonds If \code{TRUE}, then the C-hydrogens and their bonds - explicitly 
#'                 defined in an SDF - are excluded from the plot.
#'                 
#' @param bondspacer Numeric value specifying the plotting distance for 
#'                   double/triple bonds.                
#'  
#' @param colbonds Highlighting of subgraphs in main structure by providing a numeric 
#'                 vector of atom numbers, here position index in atom block. 
#'                 The bonds of connected atoms will be plotted in the color provided 
#'                 under \code{bondcol}.     
#'                 
#' @param bondcol A character or numeric vector of length one to specify the color to 
#'                use for substructure highlighting under \code{colbonds}.      
#'                
#' @param ...  Arguments to be passed to/from other methods.                                                  
#' 
#' @return Prints summary of SDF/SDFset to screen and plots their structures to 
#'         graphics device.
#'                                                                                                                                                                           
#' @keywords plot Structure
#'
#' @aliases plotStructure
#'
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#'
#' @export plotStructure
#'
#' @references 
#' ...
#' 
#' @examples
#' 
#' data(sdfbcl)
#' plotStructure(sdfbcl[[1]])
#' plotStructure(sdf = sdfbcl[[2]], atomcex = 1.2, atomnum = FALSE, 
#'               no_print_atoms = c("C"), noHbonds = TRUE, bondspacer = 0.08)

plotStructure <- function (sdf, atomcex = 1.2, atomnum = FALSE,
                           no_print_atoms = c("C"), noHbonds = TRUE,
                           bondspacer = 0.12, colbonds = NULL,
                           bondcol = "red", ...) {
  
    ChemmineR::plotStruc(sdf, atomcex = atomcex, atomnum = atomnum, 
                         no_print_atoms = no_print_atoms, 
                         noHbonds = noHbonds, bondspacer = bondspacer, 
                         colbonds = colbonds,
                         bondcol = bondcol, ...)
}


