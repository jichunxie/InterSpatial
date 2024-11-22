% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_data.R
\name{convert_gene_names_scRNA}
\alias{convert_gene_names_scRNA}
\title{Convert Gene Name Types}
\usage{
convert_gene_names_scRNA(dat_scRNA, assay_name = "SCT")
}
\arguments{
\item{assay_name}{The assay object of which we want to change the gene name type. Default value is "SCT".}

\item{dat_SCRNA}{A Seurat object containing scRNA data after SCTransform.}
}
\value{
Returns scRNA Seurat object with converted gene names.
}
\description{
Converts the gene name types of the scRNA data. This is needed if the gene name types doesn't match in scRNA and spatial data.
}
