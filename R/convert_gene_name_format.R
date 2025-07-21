#' Check if Gene Names are ENSEMBL or SYMBOL
#'
#' Determines whether the input matrix contains ENSEMBL IDs or gene symbols.
#'
#' @param input A matrix or Seurat object with gene names as rownames.
#'
#' @return A character string: either "ENSEMBEL" or "SYMBOL".
#'
#' @export
check_ensembel <- function(input)
{
  if (sum(grepl("ENS", rownames(input))) / nrow(input) > 0.8) {
    return ("ENSEMBEL")
  } else {
    return ("SYMBOL")
  }
}

#' Convert Gene Names in scRNA-seq Data
#'
#' Converts between ENSEMBL IDs and gene symbols in the given Seurat object.
#'
#' @param dat_scRNA A Seurat object to be modified in place.
#' @param assay_name Assay name (default: "SCT").
#'
#' @return A list with updated Seurat object and gene conversion data frame.
#'
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @importFrom dplyr filter
#' @export
convert_gene_names_scRNA <- function(dat_scRNA, assay_name="SCT")
{
  gene_chr_scRNA <- check_ensembel(dat_scRNA)

  if (gene_chr_scRNA == "ENSEMBEL") {
    all.genes.sc <- rownames(dat_scRNA@assays[[assay_name]]@counts)
    annots <- AnnotationDbi::select(org.Hs.eg.db, keys = all.genes.sc,
                                    columns = "SYMBOL", keytype = "ENSEMBL")

    annots <- annots %>% filter(duplicated(ENSEMBL) == FALSE)
    ind_na <- which(is.na(annots[, 2]) == TRUE)
    annots[ind_na, 2] <- paste0("NA_", ind_na)

    all.genes.sc.symbol <- annots[, 2]
    all.genes.sc.ens <- annots[, 1]
    all.genes.sc <- all.genes.sc.symbol

    rownames(dat_scRNA@assays[[assay_name]]@counts) <- all.genes.sc.symbol
    rownames(dat_scRNA@assays[[assay_name]]@data) <- all.genes.sc.symbol
    rownames(dat_scRNA@assays[[assay_name]]@scale.data) <- annots[which(annots[, 1] %in%
                                                                          rownames(dat_scRNA@assays[[assay_name]]@scale.data)), 2]

  } else {
    all.genes.sc <- rownames(dat_scRNA@assays[[assay_name]]@counts)
    annots <- AnnotationDbi::select(org.Hs.eg.db, keys = all.genes.sc,
                                    columns = "ENSEMBL", keytype = "SYMBOL")

    annots <- annots %>% filter(duplicated(SYMBOL) == FALSE)
    ind_na <- which(is.na(annots[, 2]) == TRUE)
    annots[ind_na, 2] <- paste0("NA_", ind_na)

    all.genes.sc.ens <- annots[, 2]
    all.genes.sc.symbol <- annots[, 1]
    all.genes.sc <- all.genes.sc.ens

    rownames(dat_scRNA@assays[[assay_name]]@counts) <- all.genes.sc.ens
    rownames(dat_scRNA@assays[[assay_name]]@data) <- all.genes.sc.ens
    rownames(dat_scRNA@assays[[assay_name]]@scale.data) <- annots[which(annots[, 1] %in%
                                                                          rownames(dat_scRNA@assays[[assay_name]]@scale.data)), 2]
  }

  return(list("dat_scRNA" = dat_scRNA,
              "gene_conversion" = annots))
}
