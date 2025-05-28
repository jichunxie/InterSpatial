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