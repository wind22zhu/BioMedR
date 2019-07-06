#' cluster compounds using a descriptor database
#'
#' 'clusterCMP' uses structural compound descriptors and clusters the compounds 
#' based on their pairwise distances. \code{clusterCMP} uses single linkage to 
#' measure distance between clusters when it merges clusters. It accepts both a 
#' single cutoff and a cutoff vector. By using a cutoff vector, it can generate 
#' results similar to hierarchical clustering after tree cutting.
#' 
#' \code{clusterCMP} will compute distances on the fly if \code{use.distances} 
#' is not set. Furthermore, if \code{save.distances} is not set, the distance 
#' values computed will never be stored and any distance between two compounds is 
#' guaranteed not to be computed twice. Using this method, \code{clusterCMP} can 
#' deal with large databases when a distance matrix in memory is not feasible. The 
#' speed of the clustering function should be slowed when using a transient distance 
#' calculation. 
#' When \code{save.distances} is set, \code{clusterCMP} will be forced to compute 
#' the distance matrix and save it in memory before the clustering. This is useful 
#' when additional clusterings are required in the future without re-computed the 
#' distance matrix. Set \code{save.distances} to TRUE if you only want to force the 
#' clustering to use this 2-step approach; otherwise, set it to the filename under 
#' which you want the distance matrix to be saved. After you save it, when you need 
#' to reuse the distance matrix, you can 'load' it, and supply it to \code{clusterCMP} 
#' via the \code{use.distances} argument.
#' \code{clusterCMP} supports a vector of several cutoffs. When you have multiple 
#' cutoffs, \code{clusterCMP} still guarantees that pairwise distances will never be
#' recomputed, and no copy of distances is kept in memory. It is guaranteed to be as 
#' fast as calling \code{clusterCMP} with a single cutoff that results in the longest 
#' processing time, plus some small overhead linear in processing time.
#' 
#' @param db  The desciptor database
#' 
#' @param cutoff The clustering cutoff. Can be a single value or a vector. The cutoff 
#'               gives the maximum distance between two compounds in order to group 
#'               them in the same cluster.
#' 
#' @param is.similarity Set when the cutoff supplied is a similarity cutoff. This cutoff 
#'                      is the minimum similarity value between two compounds such that 
#'                      they will be grouped in the same cluster.
#'
#' @param save.distances whether to save distance for future clustering. See details below.
#'                
#' @param use.distances Supply pre-computed distance matrix.
#' 
#' @param quiet Whether to suppress the progress information.
#' 
#' @param ... Further arguments to be passed to similarity.
#'                                                                                                                                                                                                            
#' @return 
#' Returns a \code{data.frame}. Besides a variable giving compound ID, each of the other 
#' variables in the data frame will either give the cluster IDs of compounds under some 
#' clustering cutoff, or the size of clusters that the compounds belong to. When N cutoffs 
#' are given, in total 2*N+1 variables will be generated, with N of them giving the cluster 
#' ID of each compound under each of the N cutoffs, and the other N of them giving the 
#' cluster size under each of the N cutoffs. The rows are sorted by cluster sizes.
#' 
#' @keywords cluster compounds clusterCMP
#'
#' @aliases clusterCMP
#' 
#' @author Min-feng Zhu <\email{wind2zhu@@163.com}>
#' 
#' @seealso See \code{\link{clusterStat}} for generate statistics on sizes of clusters. 
#' 
#' @export clusterCMP
#' 
#' @references
#' ...
#' 
#' @examples
#' data(sdfbcl) 
#' apbcl <- convSDFtoAP(sdfbcl)
#' fpbcl <- convAPtoFP(apbcl)
#' clusters <- clusterCMP(db = apbcl, cutoff = c(0.5, 0.85))
#' clusters2 <- clusterCMP(fpbcl, cutoff = c(0.5, 0.7), method = "Tversky") 
#' # clusters <- clusterCMP(apbcl, cutoff = 0.65, save.distances = "distmat.rda")
#' # load("distmat.rda")
#' # clusters <- clusterCMP(apbcl, cutoff = 0.60, use.distances = distmat)
#' 
clusterCMP <- function (db, cutoff, is.similarity = TRUE, save.distances = FALSE,
                        use.distances = NULL, quiet = FALSE, ...) {
  
  cluster_id = ChemmineR::cmp.cluster(db = db, cutoff = cutoff, is.similarity = is.similarity, 
                                        save.distances = save.distances, use.distances = use.distances,
                                      quiet = quiet, ...)
  return(cluster_id)
}