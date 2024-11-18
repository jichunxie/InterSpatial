#' Cell by Spot Probability Matrix
#'
#' Computes mapping probabilities of the single cells to the spatial spots.
#'
#' @param dat_scRNA A Seurat object or a matrix containing scRNA data.
#' @param dat_visium A Seurat object containing spatial transcriptomic data with a spatial assay.
#' @param p_hvg Number of highly variable genes to use from each data.
#' @param num_pca Number of principle components
#' @param alpha Level of significance of FDR control
#' @param nCores Number of cores to use for parallelization
#'
#' @return Returns the cell by spot probability matrix, each row corresponds to the spatial distribution of one cell.
#' @export
#' @importFrom parallel makeCluster stopCluster
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel

create_probability_matrix <- function(dat_scRNA, dat_visium,p_hvg=2000, num_pca=100,alpha=0.05,nCores=NULL)
{

  gene_chr_scRNA <- check_ensembel(dat_scRNA)
  gene_chr_visium <- check_ensembel(dat_visium)

  all.genes.sc <- rownames(dat_scRNA@assays[["SCT"]]@counts)
  all.genes.spatial <- rownames(dat_visium)

  hvg_sc <- VariableFeatures(dat_scRNA)

  dat_visium <- FindVariableFeatures(dat_visium,nfeatures=5000)
  hvg_dat_visium <- VariableFeatures(dat_visium)

  if(gene_chr_scRNA!=gene_chr_visium)
  {
    temp <- convert_gene_names_scRNA(dat_scRNA = dat_scRNA)
    dat_scRNA <- temp$dat_scRNA
    annots <- temp$gene_conversion
    hvg_sc <- annots[which(annots[,1] %in% hvg_sc),2]
  }


  common_genes <- Reduce(intersect,list(all.genes.sc,all.genes.spatial,c(hvg_dat_visium[1:p_hvg],hvg_sc[1:p_hvg])))

  ### Dimensionality reduction- dat_visium ###

  dat_visium <- ScaleData(dat_visium,assay = "SCT", verbose = FALSE, features = common_genes)
  dat_visium <- RunPCA(dat_visium, assay = "SCT", npcs=num_pca,verbose = FALSE, features = common_genes)

  dat_visium_cluster <- FindNeighbors(dat_visium,reduction = "pca", dims = 1:num_pca)
  dat_visium_cluster <- FindClusters(dat_visium_cluster)
  # Dimensionality reduction- scrna

  dat_scRNA <- ScaleData(dat_scRNA, features = common_genes)
  dat_scRNA <- RunPCA(dat_scRNA, assay="SCT", npcs=num_pca, features = common_genes)

  dat_scRNA <- FindNeighbors(dat_scRNA,reduction = "pca", dims = 1:num_pca)
  dat_scRNA <- FindClusters(dat_scRNA,resolution = 0.3)

  dat_visium_cell_loading <- dat_visium@reductions[["pca"]]@cell.embeddings # m * num_pca
  allen_cell_loading <- dat_scRNA@reductions[["pca"]]@cell.embeddings  # n * num_pca


  allen_loading <- as.matrix(dat_scRNA@reductions[["pca"]]@feature.loadings)  #p * num_pca

  dat_visium_Assay <- as.matrix(GetAssayData(dat_visium)[c(common_genes),])  #p * m

  allen_assay <- as.matrix(GetAssayData(dat_scRNA)[common_genes,])  #p * n
  dat_visium_loading <- as.matrix(dat_visium@reductions[["pca"]]@feature.loadings)  #p * num_pca

  common_genes <- intersect(intersect(rownames(allen_loading),rownames(allen_assay)),
                            intersect(rownames(dat_visium_loading), rownames(dat_visium_Assay)))


  allen_loading <- allen_loading[common_genes,]
  allen_assay <- allen_assay[common_genes,]
  dat_visium_loading <- dat_visium_loading[common_genes,]
  dat_visium_Assay <- dat_visium_Assay[common_genes,]

  M1 <- t(dat_visium_Assay) %*% allen_loading  #m * num_pca
  M2 <- t(allen_assay) %*% dat_visium_loading  #n * num_pca

  dat_visium_M <- cbind(dat_visium_cell_loading, M1)  #m * 2num_pca
  allen_M <- cbind(M2, allen_cell_loading)  #n * 2num_pca


  dat_visium_M_norm <- t(scale(t(dat_visium_M), center = TRUE, scale = TRUE))
  allen_M_norm <- t(scale(t(allen_M), center = TRUE, scale = TRUE))

  ## Correlation Matrix ##

  cor_M <- (dat_visium_M_norm %*% t(allen_M_norm)) / (2*num_pca-1)
  dim(cor_M)

  cor_M <- Matrix(cor_M,sparse=TRUE)




  #### FDR: Spatial Probability Distribution for Cells ####

  spot_present_ind <- which(rownames(cor_M)%in%rownames(dat_visium@images$slice@coordinates))


  spot_absent_ind <- setdiff(c(1:nrow(cor_M)),spot_present_ind)


  if(length(spot_absent_ind)>0)
  {
    cor_M[spot_absent_ind,] <- 0
  }

  # cor_M_sparse <- matrix(0,nrow=nrow(cor_M),ncol=ncol(cor_M))
  # rownames(cor_M_sparse)=rownames(cor_M)
  # colnames(cor_M_sparse)=colnames(cor_M)



  cluster_dat_visium <- dat_visium_cluster$SCT_snn_res.0.8
  cluster_sc <- dat_scRNA$SCT_snn_res.0.3
  nclust_dat_visium <- length(unique(cluster_dat_visium))
  nclust_sc <- length(unique(cluster_sc))



  unique_cluster_dat_visium <- unique(cluster_dat_visium)
  unique_cluster_sc <- unique(cluster_sc)

# return(list("corr_matrix"=cor_M,
#             "single_cell_clusters"=unique_cluster_sc,
#             "spatial_clusters"=unique_cluster_dat_visium))

if(is.null(nCores)==TRUE)
{
  cores= min(20,detectCores()-2)
}else{
  cores=nCores
}
cl <- makeCluster(cores)
registerDoParallel(cl)

# Create a grid of all combinations of 'ind1' and 'ind2'
task_grid <- expand.grid(ind1 = 1:length(unique_cluster_dat_visium),
                         ind2 = 1:length(unique_cluster_sc))

# Run computations in parallel
results_list <- foreach(
  task = iter(task_grid, by = 'row'),
  .packages = c('Matrix'),  # Include any other necessary packages
  .combine = 'list',
  .multicombine = TRUE,
  .export = c('spot_present_ind', 'cluster_dat_visium', 'cluster_sc',
              'cor_M', 'alpha', 'unique_cluster_dat_visium', 'unique_cluster_sc')
) %dopar% {
  ind1 <- task$ind1
  ind2 <- task$ind2

  # Extract indices for the current clusters
  dat_visium_ind <- intersect(
    spot_present_ind,
    which(cluster_dat_visium == unique(cluster_dat_visium)[ind1])
  )
  sc_ind <- which(cluster_sc == unique(cluster_sc)[ind2])

  # Proceed only if indices are not empty
  if (length(dat_visium_ind) > 0 && length(sc_ind) > 0) {
    # Extract the submatrix
    vec <- cor_M[dat_visium_ind, sc_ind]

    # Separate positive and negative values
    neg_values <- as.numeric(vec[which(vec < 0)])
    pos_values <- as.numeric(vec[which(vec > 0)])

    if (length(pos_values) == 0) {
      thres <- 0
    } else {
      T_seq <- sort(pos_values)
      T_seq <- seq(T_seq[1], T_seq[length(T_seq)], length.out = min(100, length(T_seq)))

      FDR <- sapply(T_seq, function(t) {
        num_neg <- length(which(neg_values < -t))
        num_pos <- length(which(pos_values > t))
        num_neg / max(1, num_pos)
      })

      if (all(FDR >= alpha)) {
        thres <- Inf
      } else {
        thres <- T_seq[min(which(FDR < alpha))]
      }
    }

    vec[which(vec <= thres)] <- 0

    list(
      dat_visium_ind = dat_visium_ind,
      sc_ind = sc_ind,
      vec = vec
    )
  } else {
    # Return NULL if no indices to process
    NULL
  }
}

# Stop the cluster after computations are done
stopCluster(cl)

# Remove NULL results (if any)
results_list <- Filter(Negate(is.null), results_list)

cor_M_sparse <- Matrix(0,nrow(cor_M),ncol(cor_M),sparse=TRUE)
# Update 'cor_M_sparse' with the results
for (res in results_list) {
  # Get indices
  dat_visium_ind <- res$dat_visium_ind
  sc_ind <- res$sc_ind
  vec_matrix <- res$vec

  # Assign the processed submatrix to 'cor_M_sparse'
  cor_M_sparse[dat_visium_ind, sc_ind] <- as.numeric(vec_matrix)
}

# Normalize cor_M_sparse_new
col_sums <- colSums(cor_M_sparse)
col_sums[col_sums == 0] <- 1  # Avoid division by zero
cor_M_sparse_new <- cor_M_sparse
cor_M_sparse_new@x <- cor_M_sparse_new@x / rep.int(col_sums, diff(cor_M_sparse_new@p))

prob_data <- t(cor_M_sparse_new)

return(prob_data)


}


#' Parallel Sum of Squared Numbers
#'
#' Demonstrates using `doParallel` inside an R package function to compute the sum of squares in parallel.
#'
#' @param numbers A numeric vector.
#' @param nCores Number of cores for parallelization (default: auto-detect).
#' @return Sum of squared numbers.
#' @examples
#' parallel_sum(1:10)
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @export
parallel_sum <- function(numbers, nCores = NULL) {
  if (!is.numeric(numbers)) stop("'numbers' must be numeric.")

  if (is.null(nCores)) {
    nCores <- max(1, parallel::detectCores() - 1)
  }

  cl <- parallel::makeCluster(nCores)
  on.exit(parallel::stopCluster(cl))  # Ensure cleanup
  doParallel::registerDoParallel(cl)

  result <- foreach::foreach(num = numbers, .combine = sum) %dopar% {
    num^2
  }

  return(result)
}
