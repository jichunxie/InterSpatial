#' Cell clustering using HGC
#'
#' Divide the single cell data into required number of clusters/meta-cells
#'
#' @param dat_scRNA A Seurat object containing scRNA data, with SCT assay
#' @param k Number of desired cell clusters or meta-cells
#' @param num_pca Number of PCA
#' @param p_hvg Number of highly variable genes
#'
#' @return Returns the cluster labels of each cell
#' @export
#'
#' @importFrom HGC FindNeighbors FindClusteringTree cutree
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




#' Avg. Spatial probability of a cluster/meta-cell
#'
#' @param prob_data The cell by spot probability matrix
#' @param cluster_labels Assigned cluster labels of the cells. Length should be equal to the nrow(prob_data), and the i^th entry should correspond to the i^th row of prob_data
#' @param cluster_id The input cluster id
#'
#' @return The avg. probability distribution of the selected cluster
#'
#' @importFrom Matrix rowSums
#'
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


#' Avg. Spatial probability of all clusters/meta-cells
#'
#' @param prob_data The cell by spot probability matrix
#' @param cluster_labels Assigned cluster labels of the cells. Length should be equal to the nrow(prob_data), and the i^th entry should correspond to the i^th row of prob_data
#' @param nCores Number of cores for parallelization
#'
#' @return The spatial probability matrix of the meta-cells. Each row corresponds to one such cluster.
#'
#' @export
compute_all_clusters_prob_matrix <- function(prob_data,cluster_labels,nCores=NULL)
{
  unique_clusters <- unique(cluster_labels)
  if(is.null(nCores)==TRUE)
  {
    cores= min(20,detectCores()-2)
  }else{
    cores=nCores
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
