#' sigvar
#' Compute the within- and across-sample variability of a set of
#'
#' #' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(sig_activity))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{sig_activity} is a single matrix containing multiple groups of samples you wish to compare.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#'
#' @return A data frame with columns corresponding to the across-sample variability and the mean within-sample variability for the provided samples.
#' @export
#'
#' @examples
#' #Import Q matrix of carcinogens in mice from Riva et al. Nature Genetics 2020
#' mice_sigs = as.data.frame(mutsig_carcinogens_mice_SBS[,c(25,2:12)])
#'
#' sigvar(sig_activity = mice_sigs, K = 11, group = "chemical")
#' @import FAVA
sigvar <- function(sig_activity,
                    K = NULL,
                    S = NULL,
                    w = NULL,
                    time = NULL,
                    group = NULL,
                    normalized = FALSE){
  var_table = full_join(fava(sig_activity = sig_activity, K = K, S = S, w = w,
                 time = time, group = group, normalized = normalized),
            gini_simpson_mean(sig_activity = sig_activity, K = K, S = S, w = w,
                 time = time, group = group))
  cols = ncol(var_table)
  colnames(var_table)[(cols-1):cols] = c("Across_sample_variability", "Mean_within_sample_variability")
  return(var_table)
}





#' sigFAVA
#' Main function for FAVA. Can test the statistical association between the variability of a Q matrix and a covariable
#' using a moving average with a Gaussian kernel
#'
#' @param Q a matrix of mutational signatures in multiple individuals
#' @param w the bandwidth of the Gaussian kernel
#'
#' @return Test statistics
#' @export
#'
#' @examples
#' #Import Q matrix of carcinogens in mice from Riva et al. Nature Genetics 2020
#' Qmice = as.data.frame(mutsig_carcinogens_mice_SBS[,c(2:12,25)])
#' Qmice[Qmice<0.1] = 0
#' plot_dots(Qmice, group = "chemical" )
#' #test
#' test_covariable(Qlist[[1]],1,1)
#' @importFrom readr read_tsv
sigFAVA <- function(){

}


