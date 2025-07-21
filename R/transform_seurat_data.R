#' SCTransform Normalization for scRNA-seq Data
#'
#' Applies SCTransform normalization to scRNA-seq data, accepting either a Seurat object or a raw matrix.
#'
#' @param dat_scRNA A Seurat object or raw gene expression matrix (genes x cells).
#' @param thres Threshold for minimum number of cells a gene must be expressed in. If NULL, it defaults to 0.5% of total cells.
#' @param metadata Optional metadata to use if a matrix is provided.
#'
#' @return A normalized Seurat object with SCT assay.
#'
#' @importFrom Seurat CreateSeuratObject SCTransform
#' @importFrom Matrix Matrix
#' @export
transform_scRNA_data <- function(dat_scRNA, thres=NULL, metadata=NULL) {
  if (!inherits(dat_scRNA, "Seurat")) {
    validate_matrix <- is.matrix(dat_scRNA) || inherits(dat_scRNA, "Matrix")
    if(validate_matrix==FALSE) {
      stop("Error: The input is neither a matrix nor a Seurat object.")
    } else {
      input_mx <- Matrix(dat_scRNA, sparse = TRUE)
      dat_scRNA <- CreateSeuratObject(input_mx, meta.data = metadata)
    }
  }
  
  if(is.null(thres)) {
    gene_thres <- floor(ncol(dat_scRNA) * 0.005)
  } else {
    gene_thres <- thres
  }
  
  if(ncol(dat_scRNA) < 2e5) {
    dat_transformed <- SCTransform(dat_scRNA, min_cells = gene_thres, verbose = FALSE)
  } else {
    dat_transformed <- SCTransform(dat_scRNA, min_cells = gene_thres, conserve.memory = TRUE)
  }
  
  return(dat_transformed)
}

#' SCTransform Normalization for Visium Spatial Transcriptomics Data
#'
#' Applies SCTransform to the Spatial assay of a Seurat object from Visium data.
#'
#' @param dat_visium A Seurat object containing a "Spatial" assay.
#' @param thres Optional threshold for minimum cell expression per gene.
#'
#' @return A Seurat object with normalized Spatial assay.
#'
#' @importFrom Seurat SCTransform
#' @export
transform_spatial_data <- function(dat_visium, thres=NULL) {
  if (!"Spatial" %in% names(dat_visium@assays)) {
    stop("Error: The input does not have a Spatial Assay.")
  }
  
  if(length(which(colSums(dat_visium@assays$Spatial@counts)==0))>0)
  {
    dat_visium <- dat_visium[,-which(colSums(dat_visium@assays$Spatial@counts)==0)]
  }
  
  n_spatial <- ncol(dat_visium)
  
  if (is.null(thres)) {
    gene_thres <- floor(n_spatial * 0.005)
  } else {
    gene_thres <- thres
  }
  
  dat_visium <- SCTransform(dat_visium, assay = "Spatial", verbose = FALSE, min_cells = gene_thres)
  return(dat_visium)
}
