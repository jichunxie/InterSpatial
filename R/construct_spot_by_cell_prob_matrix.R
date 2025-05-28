#' Create Cross-Modality Probability Matrix
#'
#' Constructs a sparse probability matrix between scRNA and Visium spatial data using PCA projections and a parallel FDR filtering step.
#'
#' @param dat_scRNA Seurat object with scRNA-seq data.
#' @param dat_visium Seurat object with spatial data.
#' @param p_hvg Number of highly variable genes to use.
#' @param num_pca Number of principal components.
#' @param alpha FDR significance threshold.
#' @param nCores Optional number of CPU cores to use.
#'
#' @return A sparse probability matrix (spatial clusters x scRNA clusters).
#'
#' @importFrom Seurat VariableFeatures FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters GetAssayData
#' @importFrom Matrix sparseMatrix
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export


create_probability_matrix <- function(dat_scRNA, dat_visium, p_hvg=2000, num_pca=100, alpha=0.05, nCores=NULL) {
  
  print("Checking gene ensemble information for scRNA and Visium data...")
  gene_chr_scRNA <- check_ensembel(dat_scRNA)
  gene_chr_visium <- check_ensembel(dat_visium)
  
  print("Extracting gene names and variable features...")
  all.genes.sc <- rownames(dat_scRNA@assays[["SCT"]]@counts)
  all.genes.spatial <- rownames(dat_visium)
  hvg_sc <- VariableFeatures(dat_scRNA)
  
  print("Finding variable features for Visium data...")
  dat_visium <- FindVariableFeatures(dat_visium, nfeatures=5000)
  hvg_dat_visium <- VariableFeatures(dat_visium)
  
  if (gene_chr_scRNA != gene_chr_visium) {
    print("Gene ensemble mismatch detected. Converting gene names...")
    temp <- convert_gene_names_scRNA(dat_scRNA = dat_scRNA)
    dat_scRNA <- temp$dat_scRNA
    annots <- temp$gene_conversion
    hvg_sc <- annots[which(annots[,1] %in% hvg_sc),2]
    all.genes.sc <- rownames(dat_scRNA@assays[["SCT"]]@counts)
  }
  
  print("Identifying common genes for integration...")
  common_genes <- Reduce(intersect, list(all.genes.sc, all.genes.spatial, c(hvg_dat_visium[1:p_hvg], hvg_sc[1:p_hvg])))
  
  print("Performing dimensionality reduction on Visium data...")
  dat_visium <- ScaleData(dat_visium, assay = "SCT", verbose = FALSE, features = common_genes)
  dat_visium <- RunPCA(dat_visium, assay = "SCT", npcs=num_pca, verbose = FALSE, features = common_genes)
  dat_visium_cluster <- FindNeighbors(dat_visium, reduction = "pca", dims = 1:num_pca)
  dat_visium_cluster <- FindClusters(dat_visium_cluster)
  
  print("Performing dimensionality reduction on scRNA data...")
  dat_scRNA <- ScaleData(dat_scRNA, features = common_genes)
  dat_scRNA <- RunPCA(dat_scRNA, assay="SCT", npcs=num_pca, features = common_genes)
  dat_scRNA <- FindNeighbors(dat_scRNA, reduction = "pca", dims = 1:num_pca)
  dat_scRNA <- FindClusters(dat_scRNA, resolution = 0.3)
  
  print("Extracting PCA embeddings and feature loadings...")
  dat_visium_cell_loading <- dat_visium@reductions[["pca"]]@cell.embeddings
  allen_cell_loading <- dat_scRNA@reductions[["pca"]]@cell.embeddings
  allen_loading <- as.matrix(dat_scRNA@reductions[["pca"]]@feature.loadings)
  dat_visium_Assay <- as.matrix(GetAssayData(dat_visium)[common_genes,])
  allen_assay <- as.matrix(GetAssayData(dat_scRNA)[common_genes,])
  dat_visium_loading <- as.matrix(dat_visium@reductions[["pca"]]@feature.loadings)
  
  print("Filtering to common genes across datasets...")
  common_genes <- intersect(intersect(rownames(allen_loading), rownames(allen_assay)),
                            intersect(rownames(dat_visium_loading), rownames(dat_visium_Assay)))
  
  allen_loading <- allen_loading[common_genes,]
  allen_assay <- allen_assay[common_genes,]
  dat_visium_loading <- dat_visium_loading[common_genes,]
  dat_visium_Assay <- dat_visium_Assay[common_genes,]
  
  print("Computing projection matrices...")
  M1 <- t(dat_visium_Assay) %*% allen_loading
  M2 <- t(allen_assay) %*% dat_visium_loading
  
  print("Combining data for correlation analysis...")
  dat_visium_M <- cbind(dat_visium_cell_loading, M1)
  allen_M <- cbind(M2, allen_cell_loading)
  
  print("Normalizing correlation matrices...")
  dat_visium_M_norm <- t(scale(t(dat_visium_M), center = TRUE, scale = TRUE))
  allen_M_norm <- t(scale(t(allen_M), center = TRUE, scale = TRUE))
  
  print("Computing correlation matrix...")
  cor_M <- (dat_visium_M_norm %*% t(allen_M_norm)) / (2 * num_pca - 1)
  
  
  print("Adjusting correlation matrix based on spatial information...")
  spot_present_ind <- which(rownames(cor_M) %in% rownames(dat_visium@images$slice@coordinates))
  spot_absent_ind <- setdiff(seq_len(nrow(cor_M)), spot_present_ind)
  
  cor_M_present <- cor_M[spot_present_ind, , drop = FALSE]
  cor_M_absent <- Matrix(0, nrow = length(spot_absent_ind), ncol = ncol(cor_M), sparse = TRUE)
  rownames(cor_M_absent) <- rownames(cor_M)[spot_absent_ind]
  colnames(cor_M_absent) <- colnames(cor_M)
  cor_M <- rbind(cor_M_present, cor_M_absent)
  
  # Reorder to original row order
  row_order <- order(c(spot_present_ind, spot_absent_ind))
  
  cor_M <- cor_M[row_order, , drop = FALSE]
  
  print("Setting up clustering variables...")
  cluster_dat_visium <- dat_visium_cluster$SCT_snn_res.0.8
  cluster_sc <- dat_scRNA$SCT_snn_res.0.3
  unique_cluster_dat_visium <- unique(cluster_dat_visium)
  unique_cluster_sc <- unique(cluster_sc)
  
  print("Determining available CPU cores for parallel processing...")
  if (is.null(nCores)) {
    cores = min(20, detectCores() - 2)
  } else {
    cores = min(nCores, detectCores() - 2)
  }
  
  
  # --- Build task list instead of full submatrix_list --- #
  task_data <- list()
  
  for (ind1 in seq_along(unique_cluster_dat_visium)) {
    for (ind2 in seq_along(unique_cluster_sc)) {
      dat_visium_ind <- intersect(spot_present_ind, which(cluster_dat_visium == unique_cluster_dat_visium[ind1]))
      sc_ind <- which(cluster_sc == unique_cluster_sc[ind2])
      
      if (length(dat_visium_ind) > 0 && length(sc_ind) > 0) {
        vec <- as.vector(cor_M[dat_visium_ind, sc_ind])
        task_data[[length(task_data) + 1]] <- list(
          dat_visium_ind = dat_visium_ind,
          sc_ind = sc_ind,
          vec = vec
        )
      }
    }
  }
  
  print(paste("Using", cores, "CPU cores for parallel processing..."))
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  print("Running FDR calculations in parallel...")
  # task_grid <- expand.grid(ind1 = seq_along(unique_cluster_dat_visium),
  #                          ind2 = seq_along(unique_cluster_sc))
  # 
  # submatrix_list <- list()
  # for (i in seq_len(nrow(task_grid))) {
  #   ind1 <- task_grid[i, "ind1"]
  #   ind2 <- task_grid[i, "ind2"]
  #   
  #   dat_visium_ind <- intersect(spot_present_ind, which(cluster_dat_visium == unique_cluster_dat_visium[ind1]))
  #   sc_ind <- which(cluster_sc == unique_cluster_sc[ind2])
  #   
  #   if (length(dat_visium_ind) > 0 && length(sc_ind) > 0) {
  #     sub_mat <- as.matrix(cor_M[dat_visium_ind, sc_ind])
  #     submatrix_list[[i]] <- list(dat_visium_ind = dat_visium_ind,
  #                                 sc_ind = sc_ind,
  #                                 sub_mat = sub_mat)
  #   } else {
  #     submatrix_list[[i]] <- NULL
  #   }
  # }
  # 
  # valid_subs <- which(!sapply(submatrix_list, is.null))
  # task_grid <- task_grid[valid_subs, ]
  # submatrix_list <- submatrix_list[valid_subs]
  # 
  # # 4. Parallel FDR using only small submatrices
  # results_list <- foreach(i = seq_along(submatrix_list),
  #                         .packages = c("Matrix"),
  #                         .combine = function(x, y) c(x, list(y)),
  #                         .init = list()) %dopar% {
  #                           sm <- submatrix_list[[i]]
  #                           vec <- as.numeric(sm$sub_mat)
  #                           neg_values <- vec[vec < 0]
  #                           pos_values <- vec[vec > 0]
  #                           
  #                           thres <- if (length(pos_values) == 0) 0 else {
  #                             T_seq <- seq(min(pos_values), max(pos_values), length.out = min(100, length(pos_values)))
  #                             FDR <- sapply(T_seq, function(t) length(neg_values[neg_values < -t]) / max(1, length(pos_values[pos_values > t])))
  #                             if (all(FDR >= alpha)) Inf else T_seq[min(which(FDR < alpha))]
  #                           }
  #                           
  #                           vec[vec <= thres] <- 0
  #                           list(dat_visium_ind = sm$dat_visium_ind,
  #                                sc_ind = sm$sc_ind,
  #                                vec = vec)
  #                         }
  
  
  
  
  results_list <- foreach(task = task_data,
                          .packages = "Matrix",
                          .combine = function(x, y) c(x, list(y)),
                          .init = list()) %dopar% {
                            vec <- task$vec
                            pos_values <- vec[vec > 0]
                            neg_values <- vec[vec < 0]
                            
                            pos_values <- sort(pos_values[!is.na(pos_values)])
                            neg_values <- sort(neg_values[!is.na(neg_values)])
                            
                            if (length(pos_values) == 0) {
                              return(NULL)
                            }
                            
                            T_seq <- seq(min(pos_values), max(pos_values), length.out = min(100, length(pos_values)))
                            n_pos <- length(pos_values)
                            n_neg <- length(neg_values)
                            
                            pos_counts <- findInterval(T_seq, pos_values, left.open = TRUE, rightmost.closed = TRUE)
                            neg_counts <- findInterval(-T_seq, neg_values, left.open = TRUE, rightmost.closed = TRUE)
                            
                            FDR <- neg_counts / pmax(n_pos - pos_counts, 1)
                            
                            if (all(FDR >= alpha)) {
                              thres = Inf
                            } else {
                              thres= T_seq[min(which(FDR < alpha))]
                            }
                            vec[vec <= thres] <- 0
                            
                            list(
                              dat_visium_ind = task$dat_visium_ind,
                              sc_ind = task$sc_ind,
                              vec = vec
                            )
                          }
  
  
  
  stopCluster(cl)
  print("Parallel computation complete. Compiling results...")
  
  
  
  # Remove NULL elements
  results_list <- Filter(Negate(is.null), results_list)
  
  all_rows <- unlist(lapply(results_list, function(res) rep(res$dat_visium_ind, length(res$sc_ind))))
  all_cols <- unlist(lapply(results_list, function(res) rep(res$sc_ind, each = length(res$dat_visium_ind))))
  all_values <- unlist(lapply(results_list, function(res) as.numeric(res$vec)))
  
  print("Constructing sparse matrix efficiently...")
  
  # ✅ Filter out zero values (ensures only nonzero values remain)
  nonzero_mask <- all_values != 0
  all_rows <- all_rows[nonzero_mask]
  all_cols <- all_cols[nonzero_mask]
  all_values <- all_values[nonzero_mask]
  # Create sparse matrix in one go
  cor_M_sparse <- sparseMatrix(
    i = all_rows,
    j = all_cols,
    x = all_values,
    dims = c(nrow(cor_M), ncol(cor_M))
  )
  
  colnames(cor_M_sparse) <- colnames(cor_M)
  rownames(cor_M_sparse) <- rownames(cor_M)
  #cor_M_sparse <- Matrix(cor_M_sparse, sparse = TRUE)
  
  # for (res_ind in 1:length(results_list)) {
  #   res <- results_list[[res_ind]]
  #   cor_M_sparse[res$dat_visium_ind, res$sc_ind] <- as.numeric(res$vec)
  #   print(res_ind)
  # }
  
  print("Normalizing final probability matrix...")
  col_sums <- colSums(cor_M_sparse)
  col_sums[col_sums == 0] <- 1
  cor_M_sparse_new <- cor_M_sparse
  cor_M_sparse_new@x <- cor_M_sparse_new@x / rep.int(col_sums, diff(cor_M_sparse_new@p))
  
  print("Probability matrix computation complete!")
  return(t(cor_M_sparse_new))
}

