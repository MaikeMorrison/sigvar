#' sigvar
#' Compute the within- and across-sample variability of mutational signature activity in one or multiple populations of samples.
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

#' sigboot
#' Estimate uncertainty in the within- and across-sample variability of mutational signature activity in one or multiple populations of samples.
#'
#' Generate bootstrap replicates of one or more tables of mutational signature activity. Compute within-and mean across-sample variability of each replicate table. Visualize the resulting distributions of variability values and compare populations with statistical tests.
#'
#' #' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)}.
#' @param seed Optional; a number to set as the random seed. Use if reproducibility of random results is desired.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(sig_activity))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{sig_activity} is a single matrix containing multiple groups of samples you wish to compare.
#' @param normalized Optional; should normalized FAVA be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized FAVA. FAVA can only be normalized if it is not weighted.
#' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE}; use \code{save_replicates = FALSE} to save memory when analyzing large datasets.
#'
#' @return A named list containing the following entries:
#' \itemize{
#' \item \code{bootstrap_replicates}: If \code{save_replicates = TRUE}, a named list of lists. Each element is named for a group provided in \code{sig_activity} and contains a list of \code{n_replicates} bootstrap replicates of the provided matrix. E.g., if \code{n_replicates = 100} and the first group in \code{sig_activity} is named \code{A}, then the first element of \code{bootstrap_replicates} is itself a list of 100 matrices, each representing a bootstrap replicate of matrix A.
#' \item \code{statistics}: A dataframe containing \code{Across_sample_variability} and \code{Mean_within_sample_variability} computed for each bootstrap replicate matrix in \code{bootstrap_replicates}. If \code{group} is specified, the preceding columns indicate which group, or combination of groups, the row corresponds to.
#' \item \code{plot_variability}: A ggplot2 scatter plot depicting the bootstrap distributions of mean within-sample variability (x-axis) and across-sample variability (y-axis) for each group in \code{sig_activity}. Ellipses provide 95% confidence
#' \item \code{plot_violin}: A ggplot2 violin plot depicting the bootstrap distribution of FAVA for each matrix in \code{sig_activity}.
#' \item \code{plot_ecdf}: A ggplot2 empirical cumulative distribution function plot depicting the bootstrap distribution of FAVA for each matrix in \code{sig_activity}.
#' \item \code{test_kruskal_wallis}: Results of a Kruskal-Wallis test performed on the bootstrap distributions of FAVA. This test is a non-parametric statistical test of whether all provided bootstrap distributions are identically distributed.
#' \item \code{test_pairwise_wilcox}: Results of a Wilcoxon rank-sum test performed on the bootstrap distributions of FAVA. This test is a non-parameteric statistical test of whether \emph{each pairwise combination} of provided bootstrap distributions is identically distributed. The result is a matrix of p-values whose entries correspond to each pair of relative abundance matrices.
#' }
#' @examples
#'
#'
#'
#' @import FAVA
sigboot <- function(sig_activity, n_replicates,
                   K = NULL,
                   seed = NULL,
                   S = NULL,
                   w = NULL,
                   time = NULL,
                   group = NULL,
                   normalized = FALSE,
                   save_replicates = FALSE){

  boot_out = bootstrap_fava(sig_activity,
                            K = K, S = S,
                 n_replicates = n_replicates, seed = seed, group = group,
                 time = time, w = w, normalized = normalized, save_replicates = save_replicates)

  statistics = boot_out$statistics %>% dplyr::select(-gini_simpson_pooled)
  colnames(statistics)[(ncol(statistics)-2):ncol(statistics)] = c("group",
                                                                  "Across_sample_variability",
                                                                  "Mean_within_sample_variability")

  plot_ellipses = ggplot2::ggplot(statistics, ggplot2::aes(x = Mean_within_sample_variability,
                                           y = Across_sample_variability,
                                           color = group)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::stat_ellipse(size = 1) +
    ggplot2::theme_bw()


  plot_boxplot = ggplot2::ggplot(tidyr::pivot_longer(statistics,
                                                     cols = c(Across_sample_variability,
                                                              Mean_within_sample_variability),
                                                     names_to = "Statistic", values_to = "Variability"),
                                 ggplot2::aes(x = group,
                                              y = Variability)) +
    ggplot2::geom_violin(ggplot2::aes(fill = group), alpha = 0.5, width = 1, color = "white") +
    ggplot2::geom_boxplot(width = 0.2, alpha = 0) +
    ggplot2::facet_wrap(~ Statistic, ncol = 1, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  if(!is.null(group)){
    test_kruskal_wallis_across = kruskal.test(x = statistics$Across_sample_variability,
                                              g = statistics$group)
    test_pairwise_wilcox_across = pairwise.wilcox.test(x = statistics$Across_sample_variability,
                                                       g = statistics$group)

    test_kruskal_wallis_within = kruskal.test(x = statistics$Mean_within_sample_variability,
                                              g = statistics$group)
    test_pairwise_wilcox_within = pairwise.wilcox.test(x = statistics$Mean_within_sample_variability,
                                                       g = statistics$group)
  }else{
    test_kruskal_wallis_across =
      "This statistical test can only be performed if a grouping variable or list of matrices is provided."
    test_pairwise_wilcox_across =
      "This statistical test can only be performed if a grouping variable or list of matrices is provided."
    test_kruskal_wallis_within =
      "This statistical test can only be performed if a grouping variable or list of matrices is provided."
    test_pairwise_wilcox_within =
      "This statistical test can only be performed if a grouping variable or list of matrices is provided."
  }


  return(list(bootstrap_replicates =
                ifelse(save_replicates,
                       boot_out$bootstrap_replicates,
                       "Set save_replicates = TRUE if you would like to save the bootstrap replicate matrices."),
              statistics = statistics,
              plot_variability = plot_ellipses,
              plot_boxplot = plot_boxplot,
              test_kruskal_wallis_across = test_kruskal_wallis_across,
              test_pairwise_wilcox_across = test_pairwise_wilcox_across,
              test_kruskal_wallis_within = test_kruskal_wallis_within,
              test_pairwise_wilcox_within = test_pairwise_wilcox_within))
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


