#' Hierarchical Graph Clustering for scRNA-seq
#'
#' Applies PCA and hierarchical clustering (HGC) to scRNA-seq data using Seurat and ClusteringTree.
#'
#' @param dat_scRNA Seurat object with SCT assay.
#' @param k Number of clusters to cut from hierarchical tree.
#' @param num_pca Number of principal components.
#' @param p_hvg Number of highly variable genes.
#' @param cell_type_column Optional. Column name in Seurat object metadata to assign cell types for cluster summary.
#'
#' @return A data frame with cell names and assigned cluster IDs.
#'
#' @importFrom Seurat VariableFeatures ScaleData RunPCA FindNeighbors
#' @importFrom HGC FindClusteringTree
#' @importFrom stats cutree
#' @importFrom dplyr %>% n
#' @export

find_cluster_scRNA_HGC <- function(dat_scRNA,k, num_pca=100,p_hvg=2000, cell_type_column = NULL)
{

  hvg_sc <- VariableFeatures(dat_scRNA)[1:p_hvg]
  dat_scRNA <- ScaleData(dat_scRNA, features = hvg_sc)
  dat_scRNA <- RunPCA(dat_scRNA, assay="SCT", npcs=num_pca, features = hvg_sc)



  set.seed(12345)
  dat_scRNA_cluster <- FindNeighbors(dat_scRNA,reduction = "pca", dims = 1:num_pca)

  dat_scRNA_cluster <- FindClusteringTree(dat_scRNA_cluster)


  clusters <- cutree(dat_scRNA_cluster@graphs$ClusteringTree, k=k)


  clusters_df <- data.frame("cell" = colnames(dat_scRNA_cluster),
                            "cluster" = clusters)

  if (!is.null(cell_type_column) && cell_type_column %in% colnames(dat_scRNA@meta.data)) {
    clusters_df$cell_type <- dat_scRNA@meta.data[clusters_df$cell, cell_type_column]

    clusters_summary <- clusters_df %>%
      group_by(cluster, cell_type) %>%
      summarize(group_size = n(), .groups = "drop")
    count_matrix_all_meta_cells <- as.matrix(xtabs(group_size ~ cluster + cell_type, data = clusters_summary))
    count_matrix_all_meta_cells  <- prop.table(count_matrix_all_meta_cells ,margin=1)

    attr(clusters_df, "cell_type_proportion_meta_cells") <- count_matrix_all_meta_cells
  }

  return(clusters_df)
}


#' Assign Dominant Cell Type Labels to All Meta-Cells
#'
#' For each cluster, computes the dominant cell type and its proportion from the cluster_labels data frame.
#'
#' @param cluster_labels A data frame with columns: 'cell', 'cluster', and 'cell_type'.
#'
#' @return A data frame with columns: cluster, assigned_cell_type, and cell_type_proportion.
#'
#' @importFrom dplyr group_by summarize %>%
#' @export
assign_cell_type_labels_to_meta_cells <- function(cluster_labels) {
  if (!all(c("cell", "cluster", "cell_type") %in% colnames(cluster_labels))) {
    stop("cluster_labels must have columns named 'cell', 'cluster', and 'cell_type'")
  }

  # Count cell type occurrences per cluster
  clusters_summary <- cluster_labels %>%
    group_by(cluster, cell_type) %>%
    summarize(group_size = n(), .groups = "drop")

  # Compute proportion matrix
  prop_matrix <- as.matrix(xtabs(group_size ~ cluster + cell_type, data = clusters_summary))
  prop_matrix <- prop.table(prop_matrix, margin = 1)

  # Determine max cell type and its proportion for each cluster
  assigned_types <- apply(prop_matrix, 1, function(row) {
    cell_type <- names(row)[which.max(row)]
    proportion <- max(row)
    c(assigned_cell_type = cell_type, cell_type_proportion = proportion)
  })

  assigned_df <- data.frame(
    meta_cell_id = as.integer(rownames(prop_matrix)),
    cell_type = assigned_types["assigned_cell_type", ],
    cell_type_proportion = as.numeric(assigned_types["cell_type_proportion", ]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  return(assigned_df)
}


#' Compute Probability Vector for a Single Cluster
#'
#' Computes the mean probabilities for all spatial regions given one cluster of scRNA-seq data.
#'
#' @param prob_data Probability matrix (cells x spatial spots), rownames are cell IDs.
#' @param cluster_labels A data frame with columns 'cell' and 'cluster'.
#' @param cluster_id The cluster ID to compute the average for.
#'
#' @return A numeric vector of average probabilities for the selected cluster.
#'
#' @importFrom Matrix rowSums
#' @export
compute_single_cluster_prob <- function(prob_data, cluster_labels, cluster_id) {
  if (!all(c("cell", "cluster") %in% colnames(cluster_labels))) {
    stop("cluster_labels must have columns named 'cell' and 'cluster'")
  }

  matched_idx <- match(rownames(prob_data), cluster_labels$cell)
  if (any(is.na(matched_idx))) {
    stop("Some rownames in prob_data do not match cell IDs in cluster_labels")
  }

  cluster_vec <- cluster_labels$cluster[matched_idx]

  nnz_cell <- which(Matrix::rowSums(prob_data) > 0)
  cluster_cells <- which(cluster_vec == cluster_id)
  nnz_cluster_cells <- intersect(cluster_cells, nnz_cell)

  if (length(nnz_cluster_cells) > 1) {
    res <- colMeans(as.matrix(prob_data[nnz_cluster_cells, , drop = FALSE]))
  } else if (length(nnz_cluster_cells) == 1) {
    res <- as.numeric(prob_data[nnz_cluster_cells, ])
  } else {
    res <- rep(0, ncol(prob_data))
  }

  return(res)
}

#' Aggregate Probabilities Across All Clusters (Non-Parallel, Zero-filled)
#'
#' Computes a matrix of average probabilities per cluster using matrix operations.
#' Returns 0 rows for clusters with no valid (nonzero) cells.
#'
#' @param prob_data Probability matrix (cells x spatial spots), rownames are cell IDs.
#' @param cluster_labels A data frame with columns 'cell' and 'cluster'.
#'
#' @return A sparse matrix (clusters x spatial spots).
#'
#' @importFrom Matrix Matrix rowSums
#' @export
compute_all_clusters_prob_matrix <- function(prob_data, cluster_labels) {
  if (!all(c("cell", "cluster") %in% colnames(cluster_labels))) {
    stop("cluster_labels must have columns named 'cell' and 'cluster'")
  }

  matched_idx <- match(rownames(prob_data), cluster_labels$cell)
  if (any(is.na(matched_idx))) {
    stop("Some rownames in prob_data do not match cell IDs in cluster_labels")
  }

  cluster_vec_all <- cluster_labels$cluster[matched_idx]
  cluster_levels <- sort(unique(cluster_vec_all))

  nonzero_rows <- which(Matrix::rowSums(prob_data) > 0)
  prob_data_nz <- prob_data[nonzero_rows, , drop = FALSE]
  cluster_vec_nz <- cluster_vec_all[nonzero_rows]

  # Sparse one-hot encoding of cluster membership
  cluster_factor <- factor(cluster_vec_nz, levels = cluster_levels)
  cluster_design <- Matrix::sparseMatrix(
    i = seq_along(cluster_factor),
    j = as.integer(cluster_factor),
    x = 1,
    dims = c(length(cluster_factor), length(cluster_levels)),
    dimnames = list(NULL, as.character(cluster_levels))
  )

  # Matrix multiplication: sum expression per cluster
  raw_sum <- t(cluster_design) %*% prob_data_nz

  rs <- Matrix::rowSums(raw_sum)

  # Precompute number of nonzeros in each row
  row_nnz <- tabulate(raw_sum@i+1 , nbins = nrow(raw_sum))

  # Normalize values only where row sum > 0
  nz_rows <- which(rs > 0)
  raw_sum_t <- t(raw_sum)
  raw_sum_t@x <- raw_sum_t@x / rep(rs[nz_rows], row_nnz[nz_rows])
  raw_sum <- t(raw_sum_t)

  return(raw_sum)

}




#' Create Meta-Cell Gene Expression Matrix
#'
#' Aggregates log-normalized gene expression from individual cells into meta-cell (cluster-level) profiles.
#'
#' @param dat_scRNA A Seurat object containing raw RNA assay counts.
#' @param cluster_labels A data frame with columns \code{"cell"} and \code{"cluster"}, mapping each cell to a cluster.
#'
#' @return A sparse matrix of dimensions (clusters x genes) with log-normalized expression.
#'         An attribute \code{"gene_symbol_type"} is attached indicating whether gene names are ENSEMBL or SYMBOL.
#'
#' @details This function performs log-normalization using CPM scaling (counts per 10k),
#'          aggregates expression per cluster using matrix operations, and returns a cluster-by-gene sparse matrix.
#'          The function also auto-detects whether the gene identifiers are in ENSEMBL format.
#'
#' @importFrom Matrix Matrix Diagonal sparseMatrix
#' @export

create_meta_cell_gene_count_matrix <- function(dat_scRNA, cluster_labels)
{
  ## Log counts Per Cells ##

  count_mx <- dat_scRNA@assays$RNA@counts
  count_mx@x <- count_mx@x / rep.int(colSums(count_mx), diff(count_mx@p))
  gc()
  count_mx <- 1e4*count_mx
  gc()
  count_mx@x <- log(1+count_mx@x,base=2)

  # set cluster name
  cell_clusters <- setNames(cluster_labels$cluster,cluster_labels$cell)

  # Match order of columns (cells) in count matrix
  cell_clusters <- cell_clusters[colnames(count_mx)]

  # Get unique cluster IDs
  clusters <- sort(unique(cell_clusters))
  k <- length(clusters)

  # Preallocate result matrix
  avg_expr_sparse <- matrix(0, nrow = nrow(count_mx), ncol = k)
  rownames(avg_expr_sparse) <- rownames(count_mx)
  colnames(avg_expr_sparse) <- paste0("Cluster_", clusters)




  # Step 1: Ensure cell order matches
  cell_clusters <- cell_clusters[colnames(count_mx)]  # align to count_mx cols
  clusters <- sort(unique(cell_clusters))
  k <- length(clusters)
  n_cells <- length(cell_clusters)

  # Step 2: Create sparse one-hot matrix (C x K)
  Z <- sparseMatrix(
    i = seq_len(n_cells),
    j = match(cell_clusters, clusters),
    x = 1,
    dims = c(n_cells, k)
  )


  # Step 3: Matrix multiplication: (G x C) %*% (C x K) = (G x K)
  cluster_sums <- count_mx %*% Z  # gene x cluster matrix of summed expression

  # Step 4: Normalize by cluster sizes
  cluster_sizes <- as.vector(colSums(Z))  # number of cells per cluster
  avg_expr_sparse <- cluster_sums %*% Diagonal(x = 1 / cluster_sizes)

  # Step 5: Add row/col names
  rownames(avg_expr_sparse) <- rownames(count_mx)
  colnames(avg_expr_sparse) <- paste0("Cluster_", clusters)



  avg_expr_sparse <- t(Matrix(avg_expr_sparse,sparse = TRUE))

  ## Convert ENS to SYMBOL if needed
  gene_chr_scRNA <- check_ensembel(dat_scRNA)

  attr(avg_expr_sparse, "gene_symbol_type") <- gene_chr_scRNA
  return (avg_expr_sparse)

}


#' Convert Gene Names from ENSEMBL to SYMBOL in Meta-Cell Matrix
#'
#' Renames gene columns in a meta-cell matrix from ENSEMBL IDs to gene symbols.
#'
#' @param count_matrix_meta_cells A sparse matrix or regular matrix with meta-cells (rows) and genes (columns).
#' @param gene_symbol_type Character. Type of current gene names. Use \code{"ENSEMBL"} to trigger conversion. Default is \code{"SYMBOL"}.
#' @param var_metadata Optional. A named vector or data frame providing a custom mapping from ENSEMBL to SYMBOL.
#'        If a data frame, it must contain columns \code{ENSEMBL} and \code{SYMBOL}.
#'
#' @return Matrix with gene columns renamed to SYMBOL. Unmapped ENSEMBL IDs are replaced with \code{"NA_<index>"} placeholders.
#'
#' @details Uses \code{AnnotationDbi::mapIds} to map ENSEMBL to SYMBOL. Handles NA values gracefully by assigning placeholders.
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @export

convert_meta_cell_gene_count_matrix_gene_symbol <- function(count_matrix_meta_cells, gene_symbol_type="SYMBOL", var_metadata=NULL)
{
  if (gene_symbol_type != "ENSEMBL") {
    message("No conversion needed. Returning original matrix.")
    return(count_matrix_meta_cells)
  }

  gene_ids <- colnames(count_matrix_meta_cells)

  # Step 1: Create mapping
  if (is.null(var_metadata)) {
    symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                        keys = gene_ids,
                                        column = "SYMBOL",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")
  } else if (is.data.frame(var_metadata)) {
    symbol_map <- setNames(var_metadata$SYMBOL, var_metadata$ENSEMBL)[gene_ids]
  } else if (is.vector(var_metadata)) {
    symbol_map <- var_metadata[gene_ids]
  } else {
    stop("var_metadata must be either NULL, a data frame, or a named vector.")
  }

  # Step 2: Replace NA values with placeholder
  na_indices <- which(is.na(symbol_map))
  symbol_map[na_indices] <- paste0("NA_", na_indices)

  # Step 3: Update column names
  colnames(count_matrix_meta_cells) <- symbol_map

  return(count_matrix_meta_cells)
}
