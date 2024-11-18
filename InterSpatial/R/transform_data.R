#' Transform scRNA Data
#'
#' Transforms scRNA data using SCTransform from the Seurat package.
#'
#' @param dat_scRNA A Seurat object or a matrix containing scRNA data.
#' @param thres Optional threshold for gene filtering. If NULL, defaults to 0.5% of the number of genes.
#' @param metadata Optional metadata to include in the Seurat object.
#' @return A Seurat object with transformed data.

#' @export
#' @importFrom Seurat CreateSeuratObject SCTransform
#' @importFrom Matrix Matrix
transform_scRNA_data <- function(dat_scRNA, thres=NULL,metadata=NULL)
{

  if (!inherits(dat_scRNA, "Seurat")) {
    validate_matrix <- is.matrix(dat_scRNA) || inherits(dat_scRNA, "Matrix")
    if(validate_matrix==FALSE)
    {
      stop("Error: The input is neither a matrix nor a Seurat object.")
    }else{
      input_mx <- Matrix(dat_scRNA,sparse = TRUE)
      dat_scRNA <- CreateSeuratObject(input_mx,meta.data = metadata)
    }

  }

  if(is.null(thres)==TRUE)
  {
    gene_thres <- floor(nrow(dat_scRNA)*0.005)
  }else{
    gene_thres <- thres
  }

  if(nrow(dat_scRNA)<1e5)
  {
    dat_transformed <- SCTransform(dat_scRNA, min_cells=gene_thres,
                                   verbose = FALSE)
  }else{
    dat_transformed <- SCTransform(dat_transformed, min_cells=gene_thres,conserve.memory=TRUE)
  }

  return (dat_transformed)
}



#' Transform Spatial Data
#'
#' Transforms Spatial Transcriptomic data using SCTransform from the Seurat package.
#'
#' @param dat_visium A Seurat object containing spatial transcriptomic data with a spatial assay.
#' @param thres Optional threshold for gene filtering. If NULL, defaults to 0.5% of the number of genes.
#' @return A Seurat object with transformed data.

#' @export
#' @importFrom Seurat CreateSeuratObject SCTransform

transform_spatial_data <- function(dat_visium, thres=NULL)
{
  if(!c("Spatial")%in% names(dat_visium@assays))
  {
    stop("Error: The input does not have a Spatial Assay.")
  }else{
    # counts_matrix <- GetAssayData(dat_visium, assay = "Spatial", slot = "counts")
    # zero_col_indices <- which(colSums(counts_matrix) == 0)
    #
    # if (length(zero_col_indices) > 0) {
    #   dat_visium <- dat_visium[, -zero_col_indices]
    # }


    p_spatial <- nrow(dat_visium)
    n_spatial <- ncol(dat_visium)


    if(is.null(thres)==TRUE)
    {
      gene_thres <- floor(n_spatial*0.005)
    }else{
      gene_thres <- thres
    }

    dat_visium <- SCTransform(dat_visium, assay = "Spatial", verbose = FALSE,min_cells=gene_thres)
    return(dat_visium)
  }


}



#' Check Gene Name Format
#'
#' Checks the gene name format in scRNA/spatial data
#'
#' @param input A Seurat object containing scRNA/spatial transcriptomic data .
#' @return Returns the format of the gene names, ENSEMBEL or SYMBOL

#' @export

check_ensembel <- function(input)
{
  if(sum(grepl("ENS",rownames(input)))/nrow(input)>0.1)
  {
    return ("ENSEMBEL")
  }else{
    return ("SYMBOL")
  }
}


#' Convert Gene Name Types
#'
#' Converts the gene name types of the scRNA data. This is needed if the gene name types doesn't match in scRNA and spatial data.
#'
#' @param dat_SCRNA A Seurat object containing scRNA data after SCTransform.
#' @param assay_name The assay object of which we want to change the gene name type. Default value is "SCT".
#' @return Returns scRNA Seurat object with converted gene names.

#' @export

convert_gene_names_scRNA <- function(dat_scRNA, assay_name="SCT")
{

  gene_chr_scRNA <- check_ensembel(dat_scRNA)

  if(gene_chr_scRNA=="ENSEMBEL")
  {
    all.genes.sc <- rownames(dat_scRNA@assays[[assay_name]]@counts)
    annots <- AnnotationDbi::select(org.Hs.eg.db, keys=all.genes.sc,
                                    columns="SYMBOL", keytype="ENSEMBL")

    annots <- annots  %>%filter(duplicated(ENSEMBL) == FALSE)
    ind_na <- which(is.na(annots[,2])==TRUE)
    annots[ind_na,2] <- paste0("NA_",ind_na)


    all.genes.sc.symbol <- annots[,2]
    all.genes.sc.ens <- annots[,1]
    all.genes.sc <- all.genes.sc.symbol
    rownames(dat_scRNA@assays[[assay_name]]@counts) = all.genes.sc.symbol
    rownames(dat_scRNA@assays[[assay_name]]@data)= all.genes.sc.symbol
    rownames(dat_scRNA@assays[[assay_name]]@scale.data)= annots[which(annots[,1] %in% rownames(dat_scRNA@assays[[assay_name]]@scale.data)),2]

  }else{
    all.genes.sc <- rownames(dat_scRNA@assays[[assay_name]]@counts)
    annots <- AnnotationDbi::select(org.Hs.eg.db, keys=all.genes.sc,
                                    columns="ENSEMBL", keytype="SYMBOL")

    annots <- annots  %>%filter(duplicated(SYMBOL) == FALSE)
    ind_na <- which(is.na(annots[,2])==TRUE)
    annots[ind_na,2] <- paste0("NA_",ind_na)


    all.genes.sc.ens <- annots[,2]
    all.genes.sc.symbol <- annots[,1]
    all.genes.sc <- all.genes.sc.ens
    rownames(dat_scRNA@assays[[assay_name]]@counts) = all.genes.sc.ens
    rownames(dat_scRNA@assays[[assay_name]]@data)= all.genes.sc.ens
    rownames(dat_scRNA@assays[[assay_name]]@scale.data)= annots[which(annots[,1] %in% rownames(dat_scRNA@assays$SCT@scale.data)),2]

  }

  return (list("dat_scRNA"=dat_scRNA,
               "gene_conversion"=annots))
}

