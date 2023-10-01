#' sigvar
#' Compute the within- and across-sample variability of mutational signature activity in one or multiple populations of samples.
#'
#' Across-sample variability is quantified by computing the population-genetic statistic Fst across the rows of the signature activity table. This method is based on the R package FAVA, which implements an Fst-based assessment of Variability across vectors of relative Abundances. Mean within-sample variability is quantified by computing the mean Gini-Simpson index across all samples. Specify \code{group} if your table of signature activities contains mutiple populations of samples which you wish to compare. Specifying \code{S} allows \code{sigvar} to account for the pairwise cosine similarity among signatures. Specifying \code{w} or \code{time} allows for uneven weighting of the samples.
#'
#' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)}.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to. Use if \code{sig_activity} is a single matrix containing multiple groups of samples you wish to compare.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' @param w Optional; a vector of length \code{nrow(sig_activity)} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight the variability measures by the distance between adjacent samples in the matrix.
#' @param normalized Optional; should normalized Fst be used as the measure of across-sample variability? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
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
#' @import FAVA
sigvar <- function(sig_activity,
                    K = NULL,
                   group = NULL,
                    S = NULL,
                    w = NULL,
                    time = NULL,
                    normalized = FALSE){
  var_table = dplyr::full_join(FAVA::fava(relab_matrix = sig_activity, K = K, S = S, w = w,
                 time = time, group = group, normalized = normalized),
                 FAVA::gini_simpson_mean(relab_matrix = sig_activity, K = K, S = S, w = w,
                 time = time, group = group))
  cols = ncol(var_table)
  colnames(var_table)[(cols-1):cols] = c("Across_sample_variability", "Mean_within_sample_variability")
  return(var_table)
}

#' #' sigboot
#' #' Estimate uncertainty in the within- and mean across-sample variability of mutational signature activity in one or multiple populations of samples.
#' #'
#' #' Generate bootstrap replicates of one or more tables of mutational signature activity. Compute across-sample and mean within-sample variability of each replicate table. Visualize the resulting distributions of variability values and compare populations with statistical tests. \code{sigboot} takes the same options as \code{sigvar}, so, as with \code{sigvar}, you can separately analyze multiple populations or groups of samples (specify \code{group}), account for cosine similarity among signatures (specify \code{S}), and incorporate non-uniform weightings of samples (specify \code{w} or \code{time}).
#' #'
#' #' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' #' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' #' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)}.
#' #' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to. Use if \code{sig_activity} is a matrix containing multiple groups of samples you wish to compare.
#' #' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' #' @param w Optional; a vector of length \code{nrow(sig_activity)} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' #' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight the variability measures by the distance between adjacent samples in the matrix.
#' #' @param normalized Optional; should Fst normalized by its upper bound conditional on the mean activity of the most abundant signature be used as the measure of across-sample variability? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' #' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE};  \code{save_replicates = FALSE} savea memory when analyzing large datasets.
#' #'
#' #' @return A named list containing the following entries:
#' #' \itemize{
#' #' \item \code{bootstrap_replicates}: If \code{save_replicates = TRUE}, a list of lists. Each element is named for a group provided in \code{sig_activity} and contains a list of \code{n_replicates} bootstrap replicates of the provided matrix. E.g., if \code{n_replicates = 100} and the first group in \code{sig_activity} is named \code{A}, then the first element of \code{bootstrap_replicates} is itself a list of 100 matrices, each representing a bootstrap replicate of matrix A.
#' #' \item \code{statistics}: A dataframe containing \code{Across_sample_variability} and \code{Mean_within_sample_variability} computed for each bootstrap replicate matrix in \code{bootstrap_replicates}. If \code{group} is specified, the preceding columns indicate which group, or combination of groups, the row corresponds to.
#' #' \item \code{plot_variability}: A ggplot2 scatter plot depicting the bootstrap distributions of mean within-sample variability (x-axis) and across-sample variability (y-axis) for each group in \code{sig_activity}. Ellipses provide 95% confidence regions.
#' #' \item \code{plot_boxplot}: A ggplot2 boxplot separately depicting the bootstrap distribution of each type of variability for each group.
#' #' \item \code{test_kruskal_wallis_across}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of across-sample variability. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' #' \item \code{test_pairwise_wilcox_across}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of across-sample variability. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed. The result is a matrix of p-values whose entries correspond to each pair of relative abundance matrices.
#' #' #' \item \code{test_kruskal_wallis_within}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of mean within-sample variability. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' #' \item \code{test_pairwise_wilcox_within}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of mean within-sample variability. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed. The result is a matrix of p-values whose entries correspond to each pair of relative abundance matrices.
#' #' @export
#' #'
#' #' @examples
#' #' # Estimate the uncertainty in the across-sample and mean within-sample variability of
#' #' mutational signatures in ESCC samples grouped by country
#' #' # We provide a cosine similarity matrix in order to account for cosine similarity among signatures
#' #' smoker_boot = sigboot(sig_activity = smoker_sigs_chen, K = 3, n_replicates = 500,
#' #'                     group = "Smoker", S = smoker_sigs_chen_cossim,
#' #'                     seed = 1)
#' #'
#' #' # The statistics entry of the output gives bootstrap replicates of the two
#' #' # variability measures for each group
#' #' head(smoker_boot$statistics)
#' #'
#' #' # sigboot visualizes these bootstrap replicates in two types of plots:
#' #' smoker_boot$plot_variability
#' #' smoker_boot$plot_boxplot
#' #'
#' #' # sigboot also computes statistics comparing the two types of variability for each group:
#' #' smoker_boot$test_pairwise_wilcox_across
#' #' smoker_boot$test_pairwise_wilcox_within
#' #'
#' #' @import FAVA
#' sigboot <- function(sig_activity, n_replicates,
#'                    K = NULL,
#'                    group = NULL,
#'                    seed = NULL,
#'                    S = NULL,
#'                    w = NULL,
#'                    time = NULL,
#'                    normalized = FALSE,
#'                    save_replicates = FALSE){
#'
#'   boot_out = bootstrap_fava(sig_activity,
#'                             K = K, S = S,
#'                  n_replicates = n_replicates, seed = seed, group = group,
#'                  time = time, w = w, normalized = normalized, save_replicates = save_replicates)
#'
#'   statistics = boot_out$statistics %>% dplyr::select(-gini_simpson_pooled)
#'   colnames(statistics)[(ncol(statistics)-2):ncol(statistics)] = c("group",
#'                                                                   "Across_sample_variability",
#'                                                                   "Mean_within_sample_variability")
#'
#'   plot_ellipses = ggplot2::ggplot(statistics, ggplot2::aes(x = Mean_within_sample_variability,
#'                                            y = Across_sample_variability,
#'                                            color = group)) +
#'     ggplot2::geom_point(alpha = 0.5) +
#'     ggplot2::stat_ellipse(size = 1) +
#'     ggplot2::theme_bw()
#'
#'
#'   plot_boxplot = ggplot2::ggplot(tidyr::pivot_longer(statistics,
#'                                                      cols = c(Across_sample_variability,
#'                                                               Mean_within_sample_variability),
#'                                                      names_to = "Statistic", values_to = "Variability"),
#'                                  ggplot2::aes(x = group,
#'                                               y = Variability)) +
#'     ggplot2::geom_violin(ggplot2::aes(fill = group), alpha = 0.5, width = 1, color = "white") +
#'     ggplot2::geom_boxplot(width = 0.2, alpha = 0) +
#'     ggplot2::facet_wrap(~ Statistic, ncol = 1, scales = "free_y") +
#'     ggplot2::theme_bw() +
#'     ggplot2::theme(legend.position = "none")
#'
#'   if(!is.null(group)){
#'     test_kruskal_wallis_across = kruskal.test(x = statistics$Across_sample_variability,
#'                                               g = statistics$group)
#'     test_pairwise_wilcox_across = pairwise.wilcox.test(x = statistics$Across_sample_variability,
#'                                                        g = statistics$group)
#'
#'     test_kruskal_wallis_within = kruskal.test(x = statistics$Mean_within_sample_variability,
#'                                               g = statistics$group)
#'     test_pairwise_wilcox_within = pairwise.wilcox.test(x = statistics$Mean_within_sample_variability,
#'                                                        g = statistics$group)
#'   }else{
#'     test_kruskal_wallis_across =
#'       "This statistical test can only be performed if a grouping variable or list of matrices is provided."
#'     test_pairwise_wilcox_across =
#'       "This statistical test can only be performed if a grouping variable or list of matrices is provided."
#'     test_kruskal_wallis_within =
#'       "This statistical test can only be performed if a grouping variable or list of matrices is provided."
#'     test_pairwise_wilcox_within =
#'       "This statistical test can only be performed if a grouping variable or list of matrices is provided."
#'   }
#'
#'
#'   return(list(bootstrap_replicates =
#'                 ifelse(save_replicates,
#'                        boot_out$bootstrap_replicates,
#'                        "Set save_replicates = TRUE if you would like to save the bootstrap replicate matrices."),
#'               statistics = statistics,
#'               plot_variability = plot_ellipses,
#'               plot_boxplot = plot_boxplot,
#'               test_kruskal_wallis_across = test_kruskal_wallis_across,
#'               test_pairwise_wilcox_across = test_pairwise_wilcox_across,
#'               test_kruskal_wallis_within = test_kruskal_wallis_within,
#'               test_pairwise_wilcox_within = test_pairwise_wilcox_within))
#' }



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


