#' Create Distance Matrix for Visium Spots
#'
#' Computes the pairwise Euclidean distance matrix for Visium spatial spots.
#'
#' @param dat_visium Seurat object with Visium spatial data.
#' @param prob_data Matrix of probabilities (cells x spots).
#' @param visium.data.dir Path to the folder containing Visium data (the one used in Load10X_Spatial).
#' @param csv_path_spot_positions Optional path to the tissue positions CSV file. If NULL, will try to auto-detect from visium.data.dir.
#'
#' @return A numeric distance matrix (spots x spots).
#'
#' @importFrom utils read.csv
#' @export
create_distance_matrix <- function( prob_data,dat_visium, visium.data.dir, csv_path_spot_positions = NULL) {

  # Fallback: auto-detect tissue_positions CSV from visium.data.dir
  if (is.null(csv_path_spot_positions)) {
    spatial_dir <- file.path(visium.data.dir, "spatial")
    csv_files <- list.files(spatial_dir, full.names = TRUE)
    csv_matches <- grep("tissue_positions.*\\.csv$", csv_files, value = TRUE)

    if (length(csv_matches) == 0) {
      stop("No CSV file containing 'tissue_positions' found in spatial directory. Please provide csv_path_spot_positions manually.")
    }

    csv_path_spot_positions <- csv_matches[1]
  }

  # Read coordinates
  positions_csv <- read.csv(csv_path_spot_positions, header = TRUE)

  # Normalize column names if needed
  if (!"pxl_row_in_fullres" %in% colnames(positions_csv)) {
    colnames(positions_csv) <- c("barcode", "in_tissue", "array_row", "array_col",
                                 "pxl_row_in_fullres", "pxl_col_in_fullres")
  }

  coefficient <- dat_visium@images$slice@scale.factors$lowres

  row_coordinates_csv <- positions_csv$pxl_row_in_fullres * coefficient
  col_coordinates_csv <- positions_csv$pxl_col_in_fullres * coefficient

  full_coordinates <- matrix(NA, nrow = ncol(prob_data), ncol = 2)
  for (i in seq_len(ncol(dat_visium))) {
    ind_csv <- which(positions_csv$barcode == colnames(prob_data)[i])
    full_coordinates[i, 1] <- row_coordinates_csv[ind_csv]
    full_coordinates[i, 2] <- col_coordinates_csv[ind_csv]
  }

  D <- as.matrix(dist(full_coordinates))
  return(D)
}



#' Compute Wasserstein Distance Between Two Cluster Distributions
#'
#' Computes the Wasserstein distance between two probability distributions using a precomputed cost matrix.
#'
#' @param index Index of the row in `ind_pair`.
#' @param ind_pair A matrix with two columns indicating index pairs (cluster1, cluster2).
#' @param Prob_Matrix A probability matrix (clusters x features).
#' @param Distance_Matrix A cost matrix used to compute Wasserstein distance.
#'
#' @return A numeric value representing the Wasserstein distance.
#'
#' @importFrom transport wasserstein
#' @export

wasserstein_pair <- function(index,ind_pair,Prob_Matrix,Distance_Matrix){
  ind1 <- ind_pair[index,1]
  ind2 <- ind_pair[index,2]
  prob_data1 <- as.vector(Prob_Matrix[ind1, ])
  prob_data2 <- as.vector(Prob_Matrix[ind2, ])

  wasserstein(prob_data1, prob_data2, costm = Distance_Matrix)
}


#' Compute Pairwise Wasserstein Distances Between Cluster Probabilities
#'
#' Computes Wasserstein distances between all source and receiver meta-cell pairs
#' based on their probability distributions across spatial spots.
#'
#' @param prob_data A matrix (meta-cells × spatial spots) of normalized probabilities.
#' @param source_cell_id Integer vector of row indices in `prob_data` indicating source meta-cells.
#' @param receiver_cell_id Integer vector of row indices in `prob_data` indicating receiver meta-cells.
#' @param Distance_Matrix A numeric distance matrix (spots × spots), typically from `create_distance_matrix()`.
#' @param out_file_prefix Prefix for naming output RDS files.
#' @param chunk_size Integer. Number of meta-cell pairs to compute per batch (default: 10,000).
#' @param out_dir Output directory for saving distance files (default: current directory).
#' @param nCores Integer. Number of CPU cores to use (default: all available minus 2, up to max 20).
#' @param parallelize Logical. Whether to use parallel computing (default: TRUE).
#'
#' @return A data frame with Wasserstein distances between meta-cell pairs.
#'
#' @importFrom transport wasserstein
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
compute_wasserstein_distance <- function(prob_data,
                                         source_cell_id,
                                         receiver_cell_id,
                                         source_name = NULL,
                                         receiver_name = NULL,
                                         Distance_Matrix,
                                         out_file_prefix = "wass_distance_meta_cells",
                                         chunk_size = 10000,
                                         out_dir = ".",
                                         nCores = NULL,
                                         parallelize = TRUE) {
  neighbour_ind <- expand.grid(ind1 = source_cell_id, ind2 = receiver_cell_id)
  non_zero_spots <- which(colSums(prob_data > 0) > 0)

  Prob_Matrix <- as.matrix(prob_data[, non_zero_spots])
  D_mat <- as.matrix(Distance_Matrix[non_zero_spots, non_zero_spots])

  total_pairs <- nrow(neighbour_ind)
  n_chunks <- ceiling(total_pairs / chunk_size)
  all_chunks <- list()

  # Setup parallel cluster outside loop (if applicable)
  if (parallelize) {
    cores <- if (is.null(nCores)) min(20, detectCores() - 2) else min(nCores, detectCores() - 2)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)  # Stop only once
  }

  for (chunk_idx in seq_len(n_chunks)) {
    cat("\n⏳ Running chunk", chunk_idx, "of", n_chunks, "\n")
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, total_pairs)

    chunk_inds <- neighbour_ind[start_idx:end_idx, ]
    time <- proc.time()

    if (parallelize) {
      clusterExport(cl, varlist = c("Prob_Matrix", "D_mat", "chunk_inds"), envir = environment())

      results <- foreach(i = seq_len(nrow(chunk_inds)),
                         .combine = 'c',
                         .packages = 'transport',
                         .errorhandling = 'pass') %dopar% {
                           p1 <- as.numeric(Prob_Matrix[chunk_inds[i, 1], ])
                           p2 <- as.numeric(Prob_Matrix[chunk_inds[i, 2], ])
                           wasserstein(p1, p2, costm = D_mat)
                         }
    } else {
      results <- numeric(nrow(chunk_inds))
      for (i in seq_len(nrow(chunk_inds))) {
        p1 <- as.numeric(Prob_Matrix[chunk_inds[i, 1], ])
        p2 <- as.numeric(Prob_Matrix[chunk_inds[i, 2], ])
        results[i] <- wasserstein(p1, p2, costm = D_mat)
      }
    }

    distance_chunk <- chunk_inds
    distance_chunk$wass_dist <- results
    all_chunks[[chunk_idx]] <- distance_chunk

    # Save individual chunk
    chunk_file <- file.path(out_dir, paste0(out_file_prefix, "_chunk", chunk_idx, ".rds"))
    saveRDS(distance_chunk, chunk_file)

    cat("✅ Saved chunk", chunk_idx, "with", nrow(chunk_inds), "pairs\n")
    print(proc.time() - time)
  }

  # Combine and save full output
  full_df <- do.call(rbind, all_chunks)
  colnames(full_df)[1:2] <- c("clust_ind_source", "clust_ind_receptor")
  if (is.null(source_name)) {
    source_name <- "source"
  }
  if (is.null(receiver_name)) {
    receiver_name <- "receiver"
  }

  saveRDS(full_df, file.path(out_dir, paste0(out_file_prefix, "_FULL.rds")))

  return(list(
    wasserstein_distance = full_df,
    source_name = source_name,
    receiver_name = receiver_name
  ))
}


#' Compute Wasserstein Distances Between Meta-Cells of Two Cell Types
#'
#' Computes Wasserstein distances between spatial probability vectors of source and receiver cell-type-specific meta-cells.
#' If `meta_cell_labels` is not provided, it derives meta-cell labels by assigning dominant cell types from `cluster_labels`.
#'
#' @param prob_data A numeric matrix (meta-cells × spatial spots) representing probability distributions.
#' @param dat_visium A Seurat Visium object containing the spatial structure of spots.
#' @param visium.data.dir Path to the Visium data directory (should contain the `spatial/` folder).
#' @param source_cell_type Character. Name of the source cell type.
#' @param receiver_cell_type Character. Name of the receiver cell type.
#' @param meta_cell_labels Optional. A data frame with columns \code{meta_cell_id}, \code{cell_type}, and optionally \code{cell_type_proportion}.
#' @param cluster_labels Optional. A data frame with columns \code{cell}, \code{cluster}, and \code{cell_type} (used if \code{meta_cell_labels} is not supplied).
#' @param thres_source_proportion Minimum required cell type proportion in a source meta-cell (default: 0.5).
#' @param thres_receiver_proportion Minimum required cell type proportion in a receiver meta-cell (default: 0.5).
#' @param out_file_prefix Prefix for saving chunked RDS outputs.
#' @param chunk_size Number of meta-cell pairs to compute per batch (default: 10,000).
#' @param out_dir Output directory for saving distance files (default: current directory).
#' @param nCores Number of CPU cores to use. If NULL, chooses \code{min(20, detectCores()-2)}.
#' @param parallelize Logical. Whether to parallelize computation (default: TRUE).
#'
#' @return A data frame with Wasserstein distances between each source-receiver meta-cell pair.
#'
#' @export
compute_wasserstein_distance_between_cell_type <- function(prob_data,
                                                           dat_visium,
                                                           visium.data.dir,
                                                           source_cell_type,
                                                           receiver_cell_type,
                                                           meta_cell_labels = NULL,
                                                           cluster_labels = NULL,
                                                           thres_source_proportion = 0.5,
                                                           thres_receiver_proportion = 0.5,
                                                           out_file_prefix = "wass_distance_meta_cells_to_meta_cells",
                                                           chunk_size = 10000,
                                                           out_dir = ".",
                                                           nCores = NULL,
                                                           parallelize = TRUE) {
  # Create spatial distance matrix
  D <- create_distance_matrix(prob_data = prob_data,
                              dat_visium = dat_visium,
                              visium.data.dir = visium.data.dir)

  nnz_meta_cells <- which(rowSums(prob_data) > 0)

  if (is.null(meta_cell_labels)) {
    if (is.null(cluster_labels)) {
      stop("Either 'meta_cell_labels' or 'cluster_labels' must be provided.")
    }
    if (!all(c("cell", "cluster", "cell_type") %in% colnames(cluster_labels))) {
      stop("cluster_labels must contain 'cell', 'cluster', and 'cell_type' columns.")
    }

    meta_cell_labels <- assign_cell_type_labels_to_meta_cells(cluster_labels)

    meta_cell_df_source <- meta_cell_labels %>%
      dplyr::filter(cell_type == source_cell_type, cell_type_proportion > thres_source_proportion)
    meta_cell_df_receiver <- meta_cell_labels %>%
      dplyr::filter(cell_type == receiver_cell_type, cell_type_proportion > thres_receiver_proportion)

    source_meta_cells <- intersect(meta_cell_df_source$meta_cell_id, nnz_meta_cells)
    receiver_meta_cells <- intersect(meta_cell_df_receiver$meta_cell_id, nnz_meta_cells)

  } else {
    if (!all(c("meta_cell_id", "cell_type") %in% colnames(meta_cell_labels))) {
      stop("meta_cell_labels must contain 'meta_cell_id' and 'cell_type' columns.")
    }

    meta_cell_labels <- meta_cell_labels[meta_cell_labels$meta_cell_id %in% nnz_meta_cells, , drop = FALSE]

    if ("cell_type_proportion" %in% colnames(meta_cell_labels)) {
      meta_cell_df_source <- meta_cell_labels %>%
        dplyr::filter(cell_type == source_cell_type, cell_type_proportion > thres_source_proportion)
      meta_cell_df_receiver <- meta_cell_labels %>%
        dplyr::filter(cell_type == receiver_cell_type, cell_type_proportion > thres_receiver_proportion)
    } else {
      meta_cell_df_source <- meta_cell_labels %>%
        dplyr::filter(cell_type == source_cell_type)
      meta_cell_df_receiver <- meta_cell_labels %>%
        dplyr::filter(cell_type == receiver_cell_type)
    }

    source_meta_cells <- intersect(meta_cell_df_source$meta_cell_id, nnz_meta_cells)
    receiver_meta_cells <- intersect(meta_cell_df_receiver$meta_cell_id, nnz_meta_cells)
  }

  message("✅ Selected ", length(source_meta_cells), " source meta-cells (", source_cell_type, ")")
  message("✅ Selected ", length(receiver_meta_cells), " receiver meta-cells (", receiver_cell_type, ")")

  if (length(source_meta_cells) == 0 || length(receiver_meta_cells) == 0) {
    warning("❌ No valid meta-cells found for one or both cell types. Returning NULL.")
    return(NULL)
  }

  distance_all_source_vs_all_receiver <- compute_wasserstein_distance(
    prob_data = prob_data,
    source_cell_id = source_meta_cells,
    receiver_cell_id = receiver_meta_cells,
    source_name = source_cell_type,
    receiver_name = receiver_cell_type,
    Distance_Matrix = D,
    out_file_prefix = out_file_prefix,
    chunk_size = chunk_size,
    out_dir = out_dir,
    nCores = nCores,
    parallelize = parallelize
  )

  return(distance_all_source_vs_all_receiver)
}
