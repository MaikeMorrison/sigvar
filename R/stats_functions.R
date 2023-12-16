# sigvar
#' Compute the within-sample diversity and across-sample heterogeneity of mutational signature activity in one or multiple populations of samples.
#'
#' Across-sample variability is quantified by computing the population-genetic statistic Fst across the rows of the signature activity table. This method is based on the R package Fst, which implements an Fst-based assessment of Variability across vectors of relative Abundances. Mean within-sample variability is quantified by computing the mean Gini-Simpson index across all samples. Specify \code{group} if your table of signature activities contains mutiple populations of samples which you wish to compare. Specifying \code{S} allows \code{sigvar} to account for the pairwise cosine similarity among signatures. Specifying \code{w} or \code{time} allows for uneven weighting of the samples.
#'
#' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)}.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to. Use if \code{sig_activity} is a single matrix containing multiple groups of samples you wish to compare.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' @param w Optional; a vector of length \code{nrow(sig_activity)} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight the variability measures by the distance between adjacent samples in the matrix.
#' @param normalized Optional; should normalized Fst be used as the measure of across-sample variability? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized Fst. Fst can only be normalized if it is not weighted.
#'
#' @return A data frame with columns corresponding to the across-sample variability and the mean within-sample variability for the provided samples.
#' @export
#'
#' @examples
#' # Compute the across-sample and mean within-sample variability of  mutational signatures
#' # in ESCC samples grouped by country
#' # We provide a cosine similarity matrix in order to account for cosine similarity among signatures
#'  sigvar(sig_activity = smoker_sigs_chen, K = 3, group = "Smoker", S = smoker_sigs_chen_cossim)
#'
sigvar <- function(sig_activity,
                   K = NULL,
                   group = NULL,
                   S = NULL,
                   w = NULL,
                   time = NULL,
                   normalized = FALSE){

  if(is.null(group)){
    var_table = data.frame(across_sample_heterogeneity = fst(relab_matrix = sig_activity, K = K, S = S, w = w,
                                                             time = time, group = group, normalized = normalized),
                           mean_within_sample_diversity = het_mean(relab_matrix = sig_activity, K = K, S = S, w = w,
                                                                   time = time, group = group))
  }else{
    var_table = dplyr::full_join(fst(relab_matrix = sig_activity, K = K, S = S, w = w,
                                     time = time, group = group, normalized = normalized),
                                 het_mean(relab_matrix = sig_activity, K = K, S = S, w = w,
                                          time = time, group = group))
    cols = ncol(var_table)
    colnames(var_table)[(cols-1):cols] = c("across_sample_heterogeneity", "mean_within_sample_diversity")
  }


  return(var_table)
}


# cossim
#' Compute the cosine similarity matrix for a set of mutational signatures.
#'
#' @param ref_sigs A matrix with K columns containing non-negative entries that sum to 1. Each column represents a signature (e.g., SBS1, SBS2), each row represents a mutation type (e.g., \code{A[C>A]A}), and each entry represents the abundance of that type in the signature.
#'
#' @return A K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' @export
#'
#' @examples
#' # Compute the cosine similarity matrix for the lung cancer
#' # in never smoker datasets from Zhang et al. 2021
#'
#'  cossim(ref_sigs = as.matrix(Sherlock_LCINS_SBS.refs[,1:14]))
#'
cossim <- function(ref_sigs){
  res = t(ref_sigs)%*%ref_sigs/matrix(sqrt(colSums(ref_sigs**2)),
                                      nrow = ncol(ref_sigs),ncol=ncol(ref_sigs))/matrix(sqrt(colSums(ref_sigs**2)),
                                                                                        nrow = ncol(ref_sigs),
                                                                                        ncol=ncol(ref_sigs),byrow = T)
  diag(res)=1
  return(res)
}
