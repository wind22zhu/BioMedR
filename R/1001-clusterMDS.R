#' visualize clustering result using multi-dimensional scaling
#' 
#' 'clusterMDS' takes clustering result returned by 'clusterCMP' and 
#' generate multi-dimensional scaling plot for visualization purpose.
#' 
#' 'clusterMDS' internally calls the 'cmdscale' function to generate a 
#' set of points in 2-D for the compounds in selected clusters.Note that for 
#' compounds in clusters smaller than the cutoff size, they will not be 
#' considered in this calculation - their entries in 'distmat' will be 
#' discarded if 'distmat' is provided, and distances involving them will 
#' not be computed if 'distmat' is not provided.
#' To determine the value for 'size.cutoff', you can use 'cluster.sizestat' 
#' to see the size distribution of clusters.
#' Because 'clusterCMP' function allows you to perform multiple clustering 
#' processes simultaneously with different cutoff values, the 'cls' parameter 
#' may point to a data frame containing multiple clustering results. The user 
#' can use 'cluster.result' to specify which result to use. By default, this 
#' is set to 1, and the first clustering result will be used in visualization. 
#' Whatever the value is, in interactive mode (described below), all clustering 
#' result will be displayed when a compound is selected in the interactive plot.
#' If the colors provided in 'color.vector' are not enough to distinguish 
#' clusters by colors, the function will silently reuse the colors, resulting 
#' multiple clusters colored in the same color. 
#' By default, 'dimensions' is set to 2, and the built-in 'plot' function will 
#' be used for plotting. If you need to do 3-Dimensional plotting, set 
#' 'dimensions' to 3, and pass the returned value to 3D plot utilities, such 
#' as 'scatterplot3d' or 'rggobi'. This package does not perform 3D plot 
#' on its own.
#' 
#' @param db  The desciptor database
#' 
#' @param cls The clustering result returned by 'clusterCMP'.
#' 
#' @param size.cutoff The cutoff size for clusters considered in this 
#' visualization. Clusters of size smaller than the cutoff will not be considered.
#'
#' @param distmat A distance matrix that corresponds to the 'db'. If not provided, 
#' it will be computed on-the-fly in an efficient manner.
#'                
#' @param color.vector Colors to be used in the plot. If the number of colors in 
#' the vector is not enough for the plot, colors will be reused. If not provided, 
#' color will be generated and randomly sampled from 'rainbow'.             
#'                         
#' @param cluster.result Used to select the clustering result if multiple 
#' clustering results are present in 'cls'.                     
#'  
#' @param dimensions Dimensionality to be used in visualization. See details.
#' 
#' @param quiet Whether to supress the progress bar.
#' 
#' @param highlight.compounds A vector of compound IDs, corresponding to 
#' compounds to be highlighted in the plot. A highlighted compound is 
#' represented as a filled circle.
#' 
#' @param highlight.color Color used for highlighted compounds. If not set, 
#' a highlighted compounds will have the same color as that used for other 
#' compounds in the same cluster.                       
#'                                                                                                                                                                                    
#' @return This function returns a data frame of MDS coordinates and 
#' clustering result. This value can be passed to 3D plot utilities such as
#' 'scatterplot3d' and 'rggobi'.
#' 
#' The last column of the output gives whether the compounds have been 
#' clicked in the interactive mode.
#' 
#' @keywords cluster visualize clusterMDS
#'
#' @aliases clusterMDS
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{clusterCMP}} for cluster compounds using a 
#' descriptor database.
#' 
#' @export clusterMDS
#' 
#' @references
#' ...
#' 
#' @examples
#' data(sdfbcl) 
#' apbcl = convSDFtoAP(sdfbcl)
#' clusters <- clusterCMP(apbcl, cutoff = c(0.5, 0.4))
#' clusterMDS(apbcl, clusters, size.cutoff = 2, quiet = TRUE)
#' 
clusterMDS <- function (db, cls, size.cutoff, distmat = NULL, 
                        color.vector = NULL, cluster.result = 1,
                        dimensions = 2, quiet = FALSE, highlight.compounds = NULL, 
                        highlight.color = NULL) {
  
    coord = ChemmineR::cluster.visualize(db, cls, size.cutoff, distmat=distmat, color.vector=color.vector, 
                                         cluster.result=cluster.result, dimensions=dimensions, quiet=quiet, 
                                         highlight.compounds=highlight.compounds, 
                                         highlight.color=highlight.color)
    return(coord)
}
