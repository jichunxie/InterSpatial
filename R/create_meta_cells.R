#' Hierarchical Graph Clustering for scRNA-seq
#'
#' Applies PCA and hierarchical clustering (HGC) to scRNA-seq data using Seurat and ClusteringTree.
#'
#' @param dat_scRNA Seurat object with SCT assay.
#' @param k Number of clusters to cut from hierarchical tree.
#' @param num_pca Number of principal components.
#' @param p_hvg Number of highly variable genes.
#'
#' @return A data frame with cell names and assigned cluster IDs.
#'
#' @importFrom Seurat VariableFeatures ScaleData RunPCA FindNeighbors
#' @importFrom HGC FindClusteringTree
#' @importFrom stats cutree
#' @export

find_cluster_scRNA_HGC <- function(dat_scRNA,k, num_pca=100,p_hvg=2000)
{
  
  hvg_sc <- VariableFeatures(dat_scRNA)[1:p_hvg]
  dat_scRNA <- ScaleData(dat_scRNA, features = hvg_sc)
  dat_scRNA <- RunPCA(dat_scRNA, assay="SCT", npcs=num_pca, features = hvg_sc)
  
  
  
  set.seed(12345)
  dat_scRNA_cluster <- FindNeighbors(dat_scRNA,reduction = "pca", dims = 1:num_pca)
  
  dat_scRNA_cluster <- FindClusteringTree(dat_scRNA_cluster)
  
  
  clusters <- cutree(dat_scRNA_cluster@graphs$ClusteringTree, k=k)
  
  
  clusters_dat_scRNA <- data.frame("cell"=colnames(dat_scRNA_cluster),
                                   "cluster"=clusters)
  
  return (clusters_dat_scRNA)
}





#' Compute Probability Vector for a Single Cluster
#'
#' Computes mean probabilities for all spatial regions given one cluster of scRNA-seq data.
#'
#' @param prob_data Probability matrix (cells x spatial spots).
#' @param cluster_labels A vector of cluster labels for each cell.
#' @param cluster_id The cluster ID to compute the average for.
#'
#' @return A numeric vector of average probabilities for the selected cluster.
#'
#' @importFrom Matrix rowSums
#' @export

compute_single_cluster_prob <- function(prob_data,cluster_labels, cluster_id)
{
  nnz_cell <- which(Matrix::rowSums(prob_data)>0)
  cells <- which(cluster_labels==cluster_id)
  nnz_cluster_cells <- intersect(cells, nnz_cell)
  
  if (length(nnz_cluster_cells) > 1) {
    res <- colMeans(as.matrix(prob_data[as.numeric(nnz_cluster_cells), ]))
  } else {
    if (length(nnz_cluster_cells) == 1) {
      res <- as.matrix(prob_data[nnz_cluster_cells, ])
    } else {
      res <- rep(0, ncol(prob_data))
    }
  }
  res
}

#' Aggregate Probabilities Across All Clusters
#'
#' Computes a new matrix of average probabilities for each cluster using parallelization.
#'
#' @param prob_data Probability matrix (cells x spatial spots).
#' @param cluster_labels A vector of cluster labels for each cell.
#' @param nCores Optional number of CPU cores to use.
#'
#' @return A sparse matrix (clusters x spatial spots).
#'
#' @importFrom Matrix Matrix
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export


compute_all_clusters_prob_matrix <- function(prob_data,cluster_labels,nCores=NULL)
{
  unique_clusters <- unique(cluster_labels)
  if(is.null(nCores)==TRUE)
  {
    cores= min(20,detectCores()-2)
  }else{
    cores= min(nCores,detectCores()-2)
  }
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  prob_data_hier_cluster_list <- foreach(cluster_id = unique_clusters,
                                         .export = c('compute_single_cluster_prob'),
                                         .packages = c('Matrix')
  ) %dopar% {
    prob_cluster <- compute_single_cluster_prob(prob_data=prob_data,cluster_labels =cluster_labels ,
                                                cluster_id=cluster_id)
    prob_cluster
  }
  stopCluster(cl)
  
  
  k  <- length(unique_clusters)
  prob_data_hier_cluster <- matrix(0,nrow=k,ncol=ncol(prob_data))
  for(i in 1:k)
  {
    prob_data_hier_cluster[i,]=prob_data_hier_cluster_list[[i]]
    
  }
  
  prob_data_hier_cluster <- Matrix(prob_data_hier_cluster,sparse=TRUE)
  return(prob_data_hier_cluster)
  
}

