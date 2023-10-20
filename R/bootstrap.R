# sigboot
#' Statistically compare the within-sample diversity and across-sample heterogeneity of the mutational signature activity in two or more groups of tumor samples.
#'
#' \code{sigboot} uses bootstrapping to evaluate differences in within-sample diversity and across-sample heterogeneity of mutational signature activity  between pairs of groups of tumor samples. \code{sigboot} takes the same options as \code{sigvar}, so, as with \code{sigvar}, you can separately analyze multiple populations or groups of samples (specify \code{group}), and account for cosine similarity among signatures (specify \code{S}). \code{sigboot} follows the bootstrapping procedure defined by Efron and Tibshirani (1993). Details on the bootstrapping procedure are available in the Methods section of the accompanying paper.
#'
#' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' @param group A string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)-length(group)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' @param normalized Optional; should the across-sample variability (Fst) be normalized by its upper bound conditional on the mean activity of the most abundant signature be used as the measure of across-sample variability? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized Fst. Fst can only be normalized if it is not weighted.
#' @param seed Optional; an integer to be used as a random seed for the simulations.
#' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE};  \code{save_replicates = FALSE} savea memory when analyzing large datasets.
#' @param alternative Optional; do you want to do a one- or two.sided test? Default is \code{alternative = "two.sided"}. If you wish to do a one-sided test, specify either \code{alternative = "less"} or \code{alternative = "greater"}.
#'
#' @return A named list containing the following entries:
#' \itemize{
#' \item \code{p_values}: The probability of observing the observed difference in variability between each pair of groups if there were no difference between groups. Computed as the fraction of bootstrap differences greater than or equal to the observed difference. Depends on what \code{alternative} is specified ("greater", "lesser", or "two.sided").
#' \item \code{bootstrap_distribution_plot}: The distribution of bootstrap replicate differences in each variability value. The observed differences are shown in red. The further the red points are from 0, the more significant the statistical difference between groups.
#' \item \code{observed_stats}: The observed diversity statistics for the groups.
#' \item \code{bootstrap_stats}: The bootstrap replicate diversity statistics for the groups.
#' \item \code{bootstrap_replicates}: The bootstrap replicate matrices, reported only if  \code{save_replicates = TRUE}.
#'
#' @export
#'
#' @examples
#' # Estimate the uncertainty in the across-sample and mean within-sample variability of
#' mutational signatures in ESCC samples grouped by country
#' # We provide a cosine similarity matrix in order to account for cosine similarity among signatures
#' smoker_boot = sigboot(sig_activity = smoker_sigs_chen, K = 3, n_replicates = 500,
#'                     group = "Smoker", S = smoker_sigs_chen_cossim,
#'                     seed = 1)
#'
#' # We find that there is a strong, significant difference between smokers and non-smokers
#' # in their mean within-sample diversity, but that the difference in across-sample heterogeneity
#' # is less strong
#' smoker_boot$P_values
#'
#' # We can visualize this result by looking at the observed differences between
#' # smokers and non-smokers and the bootstrap replicates, which assume no difference
#' # between these two groups
#' smoker_boot$bootstrap_distribution_plot
#'
#'
#' @import FAVA
sigboot <- function(sig_activity,
                    n_replicates,
                    K = ncol(sig_activity) - length(group),
                    group,
                    S = NULL, normalized = FALSE,
                    seed = NULL,
                    save_replicates = FALSE,
                    alternative = "two.sided"){

  # Set random seed (optional)
  if(!is.null(seed)){
    set.seed(seed)
  }

  # If multiple groups are provided, make a new grouping column
  multiple_groups = FALSE
  if(length(group)>1){
    multiple_groups = TRUE
    sig_activity = dplyr::mutate(sig_activity,
                                 group =  apply( sig_activity[ , group ] , 1 , paste , collapse = "_" ),
                                 .before = 1)
    group_multiple = group
    group = "group"

    group_table = dplyr::distinct(dplyr::select(sig_activity, dplyr::all_of(c("group", group_multiple))))
  }

  # How many groups are there in the data? Do we need to do multiple pairwise comparisons?
  groups = as.character(unique(sig_activity[[group]]))
  if(length(groups) == 2){
    out = pairwise_comparison(sig_activity = sig_activity, n_replicates = n_replicates, K = K,
                              group = group, S = S,normalized = normalized,
                              save_replicates = save_replicates, alternative =  alternative)

    out$P_values = out$P_values %>% data.frame() %>% t() %>%
      cbind(group_1 = groups[[1]], group_2 = groups[[2]],.)
    rownames( out$P_values ) = NULL

    if(multiple_groups){
      out$observed_stats = dplyr::left_join(group_table, out$observed_stats)
      out$bootstrap_stats = dplyr::left_join(group_table, out$bootstrap_stats)
    }
    return(out)
  }else{
    # Make a list of all unique pairs of groups
    group_pairs = t(combn(groups, 2))

    # Do the bootstrap comparison procedure for each group:
    bootstrap_list = list()
    for(pair in 1:nrow(group_pairs)){
      group_pair = group_pairs[pair,]
      sig_pair = sig_activity[sig_activity[[group]] %in% group_pair, ]

      bootstrap_list[[pair]] = pairwise_comparison(sig_activity = sig_pair,
                                                   n_replicates = n_replicates,
                                                   K = K, group = group,
                                                   S = S, normalized = normalized,
                                                   save_replicates = save_replicates,
                                                   alternative = alternative)
    }

    # Combine all elements from each category

    # 1 - P-values
    p_values = cbind(data.frame(t(combn(groups, 2))) %>% `colnames<-`(c("group_1", "group_2")),
                     lapply(bootstrap_list, function(list) list$P_values) %>%
                       do.call(rbind, .))

    # 2 - bootstrap_distribution_plot
    bootstrap_distribution_plots = lapply(bootstrap_list, function(list)
      list$bootstrap_distribution_plot) %>% `names<-`(apply(group_pairs, 1, paste0, collapse = "--"))

    # 3 - observed_stats
    observed_stats = lapply(bootstrap_list, function(list) list$observed_stats) %>%
      do.call(rbind, .) %>% dplyr::distinct() %>%
      dplyr::arrange(group)

    if(multiple_groups){ observed_stats = left_join(group_table, observed_stats) }


    # 4 - bootstrap_stats
    bootstrap_stats = lapply(bootstrap_list, function(list) list$bootstrap_stats) %>%
      do.call(rbind, .)

    if(multiple_groups){ bootstrap_stats = left_join(group_table, bootstrap_stats) }

    # 5 - bootstrap_replicates

    bootstrap_replicates = lapply(bootstrap_list, function(list) list$bootstrap_replicates) %>%
      do.call(rbind, .)

    return(list(P_values = p_values,
                bootstrap_distribution_plot = bootstrap_distribution_plots,
                observed_stats = observed_stats,
                bootstrap_stats = bootstrap_stats,
                bootstrap_replicates = bootstrap_replicates))

  }
}



pairwise_comparison <- function(sig_activity,
                                n_replicates,
                                K = ncol(sig_activity),
                                group,
                                S = NULL, normalized = FALSE,
                                save_replicates = FALSE,
                                alternative){

  # Confirm there are only two groups provided
  groups = unique(sig_activity[[group]]) %>% as.character
  if(length(groups) != 2){
    stop(paste0("There must be exactly 2 groups. There are ", length(groups), " groups in the provided sig_activity matrix."))
  }

  # Split the data into the two groups
  A = sig_activity[sig_activity[[group]] == groups[1], (ncol(sig_activity) - K + 1):ncol(sig_activity)]
  B = sig_activity[sig_activity[[group]] == groups[2], (ncol(sig_activity) - K + 1):ncol(sig_activity)]
  pooled = rbind(A, B)

  m = nrow(A)
  n = nrow(B)
  N = m+n

  # Generate bootstrap replicates of the pooled groups
  rep_list = list()
  for(rep in 1:n_replicates){
    A_rep = pooled[sample(1:N, m, replace = TRUE),]
    B_rep = pooled[sample(1:N, n, replace = TRUE),]

    rep_list[[rep]] = list(A = A_rep, B = B_rep)
  }

  # Compute statistics for each replicate
  bootstrap_stats = lapply(rep_list,
                           function(rep){
                             lapply(rep,
                                    function(matrix){
                                      c("across_sample_heterogeneity" =
                                          fst(relab_matrix = matrix, K = K, S = S, normalized = normalized),
                                        "mean_within_sample_diversity" =
                                          het_mean(relab_matrix = matrix, K = K, S = S),
                                        "pooled_diversity" =
                                          het_pooled(relab_matrix = matrix, K = K, S = S))
                                    }) %>%
                               do.call(rbind, .) %>%
                               data.frame %>%
                               tibble::rownames_to_column(var = "group")
                           }) %>%
    do.call(rbind, .)


  bootstrap_stats_A = bootstrap_stats %>% dplyr::filter(group == "A") %>% dplyr::select(-group)
  bootstrap_stats_B = bootstrap_stats %>% dplyr::filter(group == "B") %>% dplyr::select(-group)


  observed_stats_A = c("across_sample_heterogeneity" =
                         fst(relab_matrix = A, K = K, S = S, normalized = normalized),
                       "mean_within_sample_diversity" =
                         het_mean(relab_matrix = A, K = K, S = S),
                       "pooled_diversity" =
                         het_pooled(relab_matrix = A, K = K, S = S))
  observed_stats_B = c("across_sample_heterogeneity" =
                         fst(relab_matrix = B, K = K, S = S, normalized = normalized),
                       "mean_within_sample_diversity" =
                         het_mean(relab_matrix = B, K = K, S = S),
                       "pooled_diversity" =
                         het_pooled(relab_matrix = B, K = K, S = S))

  # Compute the difference between the two scrambled bootstrap populations
  bootstrap_difference = bootstrap_stats_A - bootstrap_stats_B
  observed_difference = observed_stats_A - observed_stats_B

  diff_label = paste0("Difference (", groups[[1]], " - ", groups[[2]], ")")

  # COMPUTE P-VALUES
  if(alternative == "greater"){
    p_values = sapply(colnames(bootstrap_difference), function(stat) mean((bootstrap_difference[[stat]] >= observed_difference[[stat]])) )
  } else if(alternative == "less"){
    p_values = sapply(colnames(bootstrap_difference), function(stat) mean((bootstrap_difference[[stat]] <= observed_difference[[stat]])) )
  } else if(alternative == "two.sided"){
    # two.sided p-value is from pg 212, eqn 15.26 of Efron and Tibshirani book
    p_values = sapply(colnames(bootstrap_difference),
                      function(stat) mean((abs(bootstrap_difference[[stat]]) >= abs(observed_difference[[stat]]))) )
  } else {
    stop("Valid options for alternative are 'greater', 'less', and 'two.sided'.")
  }

  # Replace P-values of 0 with less than 1 over the number of replicates
  p_values[p_values ==  0] = paste0("<", round(1/n_replicates, 10))



  # PLOT
  boot_diff_long <- bootstrap_difference %>%
    tidyr::pivot_longer(cols = 1:3, names_to = "Statistic", values_to = "Difference")

  obs_diff_long <- data.frame(Statistic = names(observed_difference),
                              Difference = observed_difference)

  plot = ggplot2::ggplot(boot_diff_long, ggplot2::aes(x = Statistic, y = Difference)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_violin(fill = "grey", color = NA, width = 1, alpha = 0.5) +
    ggplot2::geom_boxplot(width = 0.3, alpha = 0) +
    # ggbeeswarm::geom_beeswarm(alpha = 0.4) +
    ggplot2::geom_jitter(alpha = 0.4, width = 0.05) +
    ggplot2::geom_point(data = obs_diff_long, color = "red", size = 5) +
    ggplot2::ylab(diff_label) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(labels = c("FAVA" = "Across-sample\nheterogeneity",
                                         "gini_simpson_mean" = "Mean within-\nsample diversity",
                                         "gini_simpson_pooled" =  "Diversity of\npooled samples"))
  # plot
  bootstrap_replicates = rep_list

  if(!save_replicates){
    bootstrap_replicates = NULL
  }

  bootstrap_stats = cbind(group = c(rep(groups[[1]], n_replicates),
                                    rep(groups[[2]], n_replicates)),
                          rbind(bootstrap_stats_A, bootstrap_stats_B))

  observed_stats = cbind(group = c(groups[[1]], groups[[2]]),
                         rbind(observed_stats_A, observed_stats_B) %>%
                           apply(c(1,2), as.numeric) %>% `row.names<-`(NULL)) %>%
    data.frame()

  return(list(P_values = p_values,
              bootstrap_distribution_plot = plot,
              observed_stats = observed_stats,
              bootstrap_stats = bootstrap_stats,
              bootstrap_replicates = bootstrap_replicates))
}
