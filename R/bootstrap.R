# sigboot
#' Statistically compare the within-sample diversity and across-sample heterogeneity of the mutational signature activity in two or more groups of tumor samples.
#'
#' \code{sigboot} uses bootstrapping to evaluate differences in within-sample diversity and across-sample heterogeneity of mutational signature activity  between pairs of groups of tumor samples. \code{sigboot} takes the same options as \code{sigvar}, so, as with \code{sigvar}, you can separately analyze multiple populations or groups of samples (specify \code{group}), and account for cosine similarity among signatures (specify \code{S}). \code{sigboot} follows the bootstrapping procedure defined by Efron and Tibshirani (1993). Details on the bootstrapping procedure are available in the Methods section of the accompanying paper.
#'
#' @param sig_activity A matrix or data frame with \code{I} rows. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param n_replicates The number of bootstrap replicate matrices to generate for each provided relative abundance matrix.
#' @param group A string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to.
#' @param K Optional; an integer specifying the number of mutational signatures included in \code{sig_activity}. Default is \code{K=ncol(sig_activity)-length(group)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[j,k]} is the similarity between signature \code{j} and signature \code{k}, equaling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(K)}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} of \code{sig_activity} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(sig_activity), nrow(sig_activity))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row (or another continuous, increasing variable). Include if you wish to weight SVA by the temporal distance between samples.
#' @param normalized Optional; should the across-sample variability (Fst) be normalized by its upper bound conditional on the mean activity of the most abundant signature be used as the measure of across-sample variability? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized Fst. Fst can only be normalized if it is not weighted.
#' @param seed Optional; an integer to be used as a random seed for the simulations.
# #' @param save_replicates Optional; should all of the bootstrap replicate matrices be included in the output? Default is \code{save_replicates = FALSE};  \code{save_replicates = FALSE} savea memory when analyzing large datasets.
#' @param alternative Optional; do you want to do a one- or two.sided test? Default is \code{alternative = "two.sided"}. If you wish to do a one-sided test, specify either \code{alternative = "less"} or \code{alternative = "greater"}.
#'
#' @return A named list containing the following entries:
#' \itemize{
#' \item \code{p_values}: The probability of observing the observed difference in variability between each pair of groups if there were no difference between groups. Computed as the fraction of bootstrap differences greater than or equal to the observed difference. Depends on what \code{alternative} is specified ("greater", "lesser", or "two.sided").
#' \item \code{bootstrap_distribution_plot}: The distribution of bootstrap replicate differences in each variability value. The observed differences are shown in red. The further the red points are from 0, the more significant the statistical difference between groups.
#' \item \code{observed_stats}: The observed diversity statistics for the groups.
#' \item \code{bootstrap_stats}: The bootstrap replicate diversity statistics for the groups.
#' \item \code{bootstrap_replicates}: The bootstrap replicate matrices, reported only if  \code{save_replicates = TRUE}.}
#'
#' @import utils
#' @export
#'
#' @examples
#' # Estimate the uncertainty in the across-sample and mean within-sample variability of
#' # mutational signatures in ESCC samples grouped by country
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
sigboot <- function(sig_activity,
                           n_replicates,
                           group,
                           K = NULL,
                           S = NULL,
                           w = NULL,
                           time = NULL,
                           normalized = FALSE,
                           seed = NULL,
                           # save_replicates = FALSE,
                           alternative = "two.sided"){

  # To appease R cmd check
  P_value <- P_value_numeric <- Comparison <- Difference <- combn <- Statistic <- . <- NULL

  if((!is.null(time)) && (!is.null(w))){
    stop("Please specify either time or w, but not both.")
  }

  # Set random seed (optional)
  if(!is.null(seed)){
    set.seed(seed)
  }

  # If multiple grouping variables are provided, make a new grouping column
  multiple_groups = FALSE
  if(length(group)>1){
    multiple_groups = TRUE
    sig_activity = dplyr::mutate(sig_activity,
                                 group =  apply( sig_activity[ , group ] , 1 , paste , collapse = "_" ),
                                 .before = 1)
    group_multiple = group
    group = "group"

    group_table = dplyr::distinct(dplyr::select(sig_activity, dplyr::all_of(c("group", group_multiple))))

    sig_activity = sig_activity %>%
      dplyr::select(-dplyr::all_of(group_multiple))
  }


  relab_check_out = relab_checker(relab = sig_activity, K = K, group = group, time = time)

  relab_clean = relab_check_out$relab_matrix
  groups = relab_check_out$group
  times = relab_check_out$time
  K = ncol(relab_clean)
  if(is.null(w) & !is.null(time)){
    w = time_weights(times = times, group = groups)
  }

  # How many groups are there in the data? Do we need to do multiple pairwise comparisons?
  if(length(unique(groups)) < 2){
    stop(paste0("bootstrap_fava must be provided with multiple groups to compare. The grouping column '",
                group,
                "' contains only the following group: '",
                groups, "'\n" ))
  }
  if(length(unique(groups)) == 2){
    bootstrap_list = pairwise_comparison(group_pair = unique(groups), relab_clean = relab_clean,
                                         n_replicates = n_replicates, groups = groups, K = K, S = S, w = w,
                                         normalized = normalized, alternative = alternative)


    p_values = rbind(bootstrap_list$P_value) %>% data.frame %>%
      mutate(Comparison = bootstrap_list$comparison, .before = 1)

    # data.frame(Comparison = bootstrap_list$comparison, P_value = bootstrap_list$P_value) %>%
    # dplyr::mutate(P_value_numeric = as.numeric(P_value),
    #               P_value = ifelse(P_value_numeric==0, paste0("<", 1/n_replicates), P_value_numeric))

    # 2 - observed_difference
    observed_difference = rbind(bootstrap_list$observed_difference) %>% data.frame %>%
      mutate(Comparison = bootstrap_list$comparison, .before = 1)


    # 3 - bootstrap_difference
    bootstrap_difference = rbind(bootstrap_list$bootstrap_difference) %>% data.frame %>%
      mutate(Comparison = bootstrap_list$comparison, .before = 1)

    # 4 - bootstrap_plot
    bootstrap_plot = ggplot2::ggplot(data = bootstrap_difference %>%
                                       tidyr::pivot_longer(cols = c("across_sample_heterogeneity",
                                                                    "mean_within_sample_diversity"),
                                                           names_to = "Statistic", values_to = "Difference"),
                                     mapping = ggplot2::aes(x = Statistic, y = Difference)) +
      ggplot2::geom_violin(fill = "grey", color = NA, alpha = 0.6) +
      ggplot2::geom_boxplot(alpha = 0, width = 0.3) +
      ggplot2::geom_jitter(alpha = 0.3, width = 0.05, size = 2) +
      ggplot2::theme_bw() +
      ggplot2::ylab(paste0("Difference in statistic\n(", unique(bootstrap_difference$Comparison), ")")) +
      ggplot2::geom_point(data = observed_difference %>%
                            tidyr::pivot_longer(cols = c("across_sample_heterogeneity",
                                                         "mean_within_sample_diversity"),
                                                names_to = "Statistic", values_to = "Difference"), color = "red", size = 5) +
      ggplot2::scale_x_discrete(labels = c("across_sample_heterogeneity" = "Across-sample\nheterogeneity",
                                           "mean_within_sample_diversity" = "Mean within-\nsample diversity"))



    return(list(P_values = p_values,
                bootstrap_distribution_plot = bootstrap_plot,
                observed_difference = observed_difference,
                bootstrap_difference = bootstrap_difference))

  }else{
    # Make a list of all unique pairs of groups
    group_pairs = t(combn(unique(groups), 2))

    # Do the bootstrap comparison procedure for each group:
    bootstrap_list = list()
    for(pair in 1:nrow(group_pairs)){
      group_pair = group_pairs[pair,]

      bootstrap_list[[pair]] = pairwise_comparison(group_pair = group_pair, relab_clean = relab_clean, n_replicates = n_replicates,
                                                   groups = groups, K = K, S = S, w = w,
                                                   normalized = normalized, alternative = alternative)


    }

    # Combine all elements from each category

    # 1 - P-values
    p_values = lapply(bootstrap_list, function(pair) c(Comparison = pair$comparison, pair$P_value)) %>%
      do.call(rbind, .) %>%
      data.frame() %>%
      dplyr::mutate(across(-Comparison, as.numeric))
                    # P_value = ifelse(P_value_numeric==0, paste0("<", 1/n_replicates), P_value_numeric))

    # cbind(data.frame(t(combn(groups, 2))) %>% `colnames<-`(c("group_1", "group_2")),
    #                lapply(bootstrap_list, function(list) list$P_values) %>%
    #                  do.call(rbind, .))

    # # 2 - bootstrap_distribution_plot
    # bootstrap_distribution_plots = lapply(bootstrap_list, function(list)
    #   list$bootstrap_distribution_plot) %>% `names<-`(apply(group_pairs, 1, paste0, collapse = "--"))

    # 2 - observed_difference
    observed_difference = lapply(bootstrap_list, function(list) c(Comparison = list$comparison, list$observed_difference)) %>%
      do.call(rbind, .) %>% data.frame   %>%
      dplyr::mutate(across(-Comparison, as.numeric),
                    Comparison= unlist(Comparison))
    # observed_difference$Comparison = stringr::str_replace_all(observed_difference$Comparison, ' - ', " -\n")
    # if(multiple_groups){ observed_difference = left_join(group_table, observed_difference) }


    # 3 - bootstrap_difference
    bootstrap_difference = lapply(bootstrap_list, function(list) mutate(list$bootstrap_difference, Comparison =list$comparison, .before = 1)) %>%
      do.call(rbind, .) %>%
      data.frame()
      # tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Comparison", values_to = "Difference")
    # bootstrap_difference$Comparison = stringr::str_replace_all(bootstrap_difference$Comparison, '\\.\\.\\.', " -\n")

    # if(multiple_groups){ bootstrap_difference = left_join(group_table, bootstrap_difference) }

    # 4 - bootstrap_plot

    bootstrap_plot = ggplot2::ggplot(data = bootstrap_difference %>%
                      tidyr::pivot_longer(cols = c("across_sample_heterogeneity",
                                                   "mean_within_sample_diversity"),
                                          names_to = "Statistic", values_to = "Difference"),
                    mapping = ggplot2::aes(x = Statistic, y = Difference)) +
      ggplot2::geom_violin(fill = "grey", color = NA, alpha = 0.6) +
      ggplot2::geom_boxplot(alpha = 0, width = 0.3) +
      ggplot2::geom_jitter(alpha = 0.3, width = 0.05, size = 2) +
      ggplot2::theme_bw() +
      ggplot2::ylab(paste0("Difference in statistic")) +
      ggplot2::facet_wrap(~ Comparison) +
      ggplot2::geom_point(data = observed_difference %>%
                            tidyr::pivot_longer(cols = c("across_sample_heterogeneity",
                                                         "mean_within_sample_diversity"),
                                                names_to = "Statistic", values_to = "Difference"), color = "red", size = 5) +
          ggplot2::scale_x_discrete(labels = c("across_sample_heterogeneity" = "Across-sample\nheterogeneity",
                                               "mean_within_sample_diversity" = "Mean within-\nsample diversity"))

    # bootstrap_plot = ggplot2::ggplot(data = bootstrap_difference, mapping = ggplot2::aes(x = Comparison, y = Difference)) +
    #   ggplot2::geom_violin(fill = "grey", color = NA, alpha = 0.6) +
    #   ggplot2::geom_boxplot(alpha = 0, width = 0.3) +
    #   ggplot2::geom_jitter(alpha = 0.3, width = 0.05, size = 2) +
    #   ggplot2::theme_bw() +
    #   ggplot2::ylab("Difference in FAVA") +
    #   ggplot2::geom_point(data = observed_difference, color = "red", size = 5)

    # # 5 - bootstrap_replicates
    #
    # bootstrap_replicates = lapply(bootstrap_list, function(list) list$bootstrap_replicates) %>%
    #   do.call(rbind, .)

    return(list(P_values = p_values,
                bootstrap_distribution_plot = bootstrap_plot,
                observed_difference = observed_difference,
                bootstrap_difference = bootstrap_difference))

  }
}



pairwise_comparison <- function(group_pair,
                                relab_clean,
                                n_replicates,
                                groups,

                                K,
                                S,
                                w,
                                # time,
                                normalized,

                                # save_replicates,
                                alternative){
  # To appease R cmd check
  . <- NULL

  # Confirm there are only two groups provided
  if(length(group_pair) != 2){
    stop(paste0("There must be exactly 2 groups. There are ", length(groups), " groups in the provided relab_pair matrix."))
  }



  # Split the data into the two groups
  A = relab_clean[groups == group_pair[[1]],]
  B = relab_clean[groups == group_pair[[2]],]
  pooled = rbind(A, B)


  m = nrow(A)
  n = nrow(B)
  N = m+n
  if(!is.null(w)){
    wA = w[groups == group_pair[[1]]]
    wB = w[groups == group_pair[[2]]]
    wPooled = c(wA, wB)
  }else{wPooled = NULL; wA = NULL; wB = NULL}

  # Generate bootstrap replicates of the pooled groups
  rep_list = list()
  for(rep in 1:n_replicates){
    a_samp = sample(1:N, m, replace = TRUE)
    b_samp = sample(1:N, n, replace = TRUE)

    rep_list[[rep]] = list(A = a_samp, B = b_samp)
  }

  # Compute statistics for each replicate
  bootstrap_differences = lapply(rep_list,
                                 function(rep){
                                   # Compute for A
                                   sigvar(sig_activity = pooled[rep$A,], K = K, S = S, normalized = normalized,
                                          w = `if`(is.null(NULL), NULL, wPooled[rep$A]/sum(wPooled[rep$A]))) -
                                     # Compute for B
                                     sigvar(sig_activity = pooled[rep$B,], K = K, S = S, normalized = normalized,
                                            w = `if`(is.null(NULL), NULL, wPooled[rep$B]/sum(wPooled[rep$B])))
                                 }

  ) %>%
    do.call(rbind, .) %>%
    data.frame()

  # bootstrap_stats$Difference = bootstrap_stats[,group_pair[[1]]] - bootstrap_stats[,group_pair[[2]]]

  # Compute the difference between the two original populations
  observed_difference = sigvar(sig_activity = A, K = K, S = S, normalized = normalized, w = wA) -
    sigvar(sig_activity = B, K = K, S = S, normalized = normalized, w = wB)

  diff_label = paste0( #"Difference (",
    group_pair[[1]], " - ", group_pair[[2]]#, ")"
  )

  # COMPUTE P-VALUES
  if(alternative == "greater"){
    p_val = mapply(col = bootstrap_differences, true = observed_difference,
                   FUN = function(col, true){mean(col >= true)})
  } else if(alternative == "less"){
    p_val =  mapply(col = bootstrap_differences, true = observed_difference,
                    FUN = function(col, true){mean(col <= true)})
  } else if(alternative == "two.sided"){
    # two.sided p-value is from pg 212, eqn 15.26 of Efron and Tibshirani book
    p_val =  mapply(col = bootstrap_differences, true = observed_difference,
                    FUN = function(col, true){mean(abs(col) >= abs(true))})
  } else {
    stop("Valid options for alternative are 'greater', 'less', and 'two.sided'.")
  }

  return(list(comparison = diff_label,
              P_value = p_val,
              observed_difference = observed_difference,
              bootstrap_difference = bootstrap_differences
  ))
}





















# sigboot <- function(sig_activity,
#                     n_replicates,
#                     K = ncol(sig_activity) - length(group),
#                     group,
#                     S = NULL,
#                     w = NULL,
#                     time = NULL,
#                       normalized = FALSE,
#                     seed = NULL,
#                     save_replicates = FALSE,
#                     alternative = "two.sided"){
#
#   # Satisfy R cmd check
#   . <- NULL
#
#
#   if((!is.null(time)) && (!is.null(w))){
#     stop("Please specify either time or w, but not both.")
#   }
#
#
#   # Set random seed (optional)
#   if(!is.null(seed)){
#     set.seed(seed)
#   }
#
#   # If multiple groups are provided, make a new grouping column
#   multiple_groups = FALSE
#   if(length(group)>1){
#     multiple_groups = TRUE
#     sig_activity = dplyr::mutate(sig_activity,
#                                  group =  apply( sig_activity[ , group ] , 1 , paste , collapse = "_" ),
#                                  .before = 1)
#     group_multiple = group
#     group = "group"
#
#     group_table = dplyr::distinct(dplyr::select(sig_activity, dplyr::all_of(c("group", group_multiple))))
#   }
#
#
#
#   relab_check_out = relab_checker(relab = sig_activity, K = K, group = group, time = time)
#
#   relab_clean = relab_check_out$relab_matrix
#   groups = relab_check_out$group
#   times = relab_check_out$time
#   K = ncol(relab_clean)
#   if(is.null(w) & !is.null(time)){
#     w = time_weights(times = times, group = groups)
#   }
#
#
#
#   # How many groups are there in the data? Do we need to do multiple pairwise comparisons?
#   if(length(unique(groups)) < 2){
#     stop(paste0("sigboot must be provided with multiple groups to compare. The grouping column '",
#                 group,
#                 "' contains only the following group: '",
#                 groups, "'\n" ))
#   }
#   groups = as.character(unique(sig_activity[[group]]))
#   if(length(unique(groups)) == 2){
#     out = pairwise_comparison(sig_activity = sig_activity, n_replicates = n_replicates, K = K,
#                               groups = groups, S = S, w = w, time = time,
#                               normalized = normalized,
#                               save_replicates = save_replicates, alternative =  alternative)
#
#     out$P_values = out$P_values %>% data.frame() %>% t() %>%
#       cbind(group_1 = groups[[1]], group_2 = groups[[2]],.)
#     rownames( out$P_values ) = NULL
#
#     if(multiple_groups){
#       out$observed_stats = dplyr::left_join(group_table, out$observed_stats)
#       out$bootstrap_stats = dplyr::left_join(group_table, out$bootstrap_stats)
#     }
#     return(out)
#   }else{
#     # Make a list of all unique pairs of groups
#     group_pairs = t(utils::combn(groups, 2))
#
#     # Do the bootstrap comparison procedure for each group:
#     bootstrap_list = list()
#     for(pair in 1:nrow(group_pairs)){
#       group_pair = group_pairs[pair,]
#       sig_pair = sig_activity[sig_activity[[group]] %in% group_pair, ]
#
#       bootstrap_list[[pair]] = pairwise_comparison(sig_activity = sig_pair,
#                                                    n_replicates = n_replicates,
#                                                    K = K, group = group,
#                                                    S = S, w = w, time = time,
#                                                    normalized = normalized,
#                                                    save_replicates = save_replicates,
#                                                    alternative = alternative)
#     }
#
#     # Combine all elements from each category
#
#     # 1 - P-values
#     p_values = cbind(data.frame(t(utils::combn(groups, 2))) %>% `colnames<-`(c("group_1", "group_2")),
#                      lapply(bootstrap_list, function(list) list$P_values) %>%
#                        do.call(rbind, .))
#
#     # 2 - bootstrap_distribution_plot
#     bootstrap_distribution_plots = lapply(bootstrap_list, function(list)
#       list$bootstrap_distribution_plot) %>% `names<-`(apply(group_pairs, 1, paste0, collapse = "--"))
#
#     # 3 - observed_stats
#     observed_stats = lapply(bootstrap_list, function(list) list$observed_stats) %>%
#       do.call(rbind, .) %>% dplyr::distinct() %>%
#       dplyr::arrange(group)
#
#     if(multiple_groups){ observed_stats = dplyr::left_join(group_table, observed_stats) }
#
#
#     # 4 - bootstrap_stats
#     bootstrap_stats = lapply(bootstrap_list, function(list) list$bootstrap_stats) %>%
#       do.call(rbind, .)
#
#     if(multiple_groups){ bootstrap_stats = dplyr::left_join(group_table, bootstrap_stats) }
#
#     # 5 - bootstrap_replicates
#
#     bootstrap_replicates = lapply(bootstrap_list, function(list) list$bootstrap_replicates) %>%
#       do.call(rbind, .)
#
#     return(list(P_values = p_values,
#                 bootstrap_distribution_plot = bootstrap_distribution_plots,
#                 observed_stats = observed_stats,
#                 bootstrap_stats = bootstrap_stats,
#                 bootstrap_replicates = bootstrap_replicates))
#
#   }
# }
#
#
#
# pairwise_comparison <- function(sig_activity,
#                                 n_replicates,
#                                 K = ncol(sig_activity),
#                                 group,
#                                 S = NULL,
#                                 w = NULL,
#                                 time = NULL,
#                                 normalized = FALSE,
#                                 save_replicates = FALSE,
#                                 alternative){
#   Statistic <- Difference <- . <- NULL
#
#   # Confirm there are only two groups provided
#   groups = unique(sig_activity[[group]]) %>% as.character
#   if(length(groups) != 2){
#     stop(paste0("There must be exactly 2 groups. There are ", length(groups), " groups in the provided sig_activity matrix."))
#   }
#
#   # Split the data into the two groups
#   A = sig_activity[sig_activity[[group]] == groups[1], (ncol(sig_activity) - K + 1):ncol(sig_activity)]
#   B = sig_activity[sig_activity[[group]] == groups[2], (ncol(sig_activity) - K + 1):ncol(sig_activity)]
#   pooled = rbind(A, B)
#
#   m = nrow(A)
#   n = nrow(B)
#   N = m+n
#
#   # Generate bootstrap replicates of the pooled groups
#   rep_list = list()
#   for(rep in 1:n_replicates){
#     A_rep = pooled[sample(1:N, m, replace = TRUE),]
#     B_rep = pooled[sample(1:N, n, replace = TRUE),]
#
#     rep_list[[rep]] = list(A = A_rep, B = B_rep)
#   }
#
#   # Compute statistics for each replicate
#   bootstrap_stats = lapply(rep_list,
#                            function(rep){
#                              lapply(rep,
#                                     function(matrix){
#                                       c("across_sample_heterogeneity" =
#                                           fst(relab_matrix = matrix, K = K, S = S, normalized = normalized),
#                                         "mean_within_sample_diversity" =
#                                           het_mean(relab_matrix = matrix, K = K, S = S)#,
#                                         # "pooled_diversity" =
#                                         #   het_pooled(relab_matrix = matrix, K = K, S = S)
#                                         )
#                                     }) %>%
#                                do.call(rbind, .) %>%
#                                data.frame %>%
#                                tibble::rownames_to_column(var = "group")
#                            }) %>%
#     do.call(rbind, .)
#
#
#   bootstrap_stats_A = bootstrap_stats %>% dplyr::filter(group == "A") %>% dplyr::select(-group)
#   bootstrap_stats_B = bootstrap_stats %>% dplyr::filter(group == "B") %>% dplyr::select(-group)
#
#
#   observed_stats_A = c("across_sample_heterogeneity" =
#                          fst(relab_matrix = A, K = K, S = S, normalized = normalized),
#                        "mean_within_sample_diversity" =
#                          het_mean(relab_matrix = A, K = K, S = S)#,
#                        # "pooled_diversity" =
#                        #   het_pooled(relab_matrix = A, K = K, S = S)
#                        )
#   observed_stats_B = c("across_sample_heterogeneity" =
#                          fst(relab_matrix = B, K = K, S = S, normalized = normalized),
#                        "mean_within_sample_diversity" =
#                          het_mean(relab_matrix = B, K = K, S = S)#,
#                        # "pooled_diversity" =
#                        #   het_pooled(relab_matrix = B, K = K, S = S)
#                        )
#
#   # Compute the difference between the two scrambled bootstrap populations
#   bootstrap_difference = bootstrap_stats_A - bootstrap_stats_B
#   observed_difference = observed_stats_A - observed_stats_B
#
#   diff_label = paste0("Difference (", groups[[1]], " - ", groups[[2]], ")")
#
#   # COMPUTE P-VALUES
#   if(alternative == "greater"){
#     p_values = sapply(colnames(bootstrap_difference), function(stat) mean((bootstrap_difference[[stat]] >= observed_difference[[stat]])) )
#   } else if(alternative == "less"){
#     p_values = sapply(colnames(bootstrap_difference), function(stat) mean((bootstrap_difference[[stat]] <= observed_difference[[stat]])) )
#   } else if(alternative == "two.sided"){
#     # two.sided p-value is from pg 212, eqn 15.26 of Efron and Tibshirani book
#     p_values = sapply(colnames(bootstrap_difference),
#                       function(stat) mean((abs(bootstrap_difference[[stat]]) >= abs(observed_difference[[stat]]))) )
#   } else {
#     stop("Valid options for alternative are 'greater', 'less', and 'two.sided'.")
#   }
#
#   # Replace P-values of 0 with less than 1 over the number of replicates
#   p_values[p_values ==  0] = paste0("<", round(1/n_replicates, 10))
#
#
#
#   # PLOT
#   boot_diff_long <- bootstrap_difference %>%
#     tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Statistic", values_to = "Difference")
#
#   obs_diff_long <- data.frame(Statistic = names(observed_difference),
#                               Difference = observed_difference)
#
#   plot = ggplot2::ggplot(boot_diff_long, ggplot2::aes(x = Statistic, y = Difference)) +
#     ggplot2::geom_hline(yintercept = 0) +
#     ggplot2::geom_violin(fill = "grey", color = NA, width = 1, alpha = 0.5) +
#     ggplot2::geom_boxplot(width = 0.3, alpha = 0) +
#     # ggbeeswarm::geom_beeswarm(alpha = 0.4) +
#     ggplot2::geom_jitter(alpha = 0.4, width = 0.05) +
#     ggplot2::geom_point(data = obs_diff_long, color = "red", size = 5) +
#     ggplot2::ylab(diff_label) +
#     ggplot2::theme_bw() +
#     ggplot2::scale_x_discrete(labels = c("FAVA" = "Across-sample\nheterogeneity",
#                                          "gini_simpson_mean" = "Mean within-\nsample diversity",
#                                          "gini_simpson_pooled" =  "Diversity of\npooled samples"))
#   # plot
#   bootstrap_replicates = rep_list
#
#   if(!save_replicates){
#     bootstrap_replicates = NULL
#   }
#
#   bootstrap_stats = cbind(group = c(rep(groups[[1]], n_replicates),
#                                     rep(groups[[2]], n_replicates)),
#                           rbind(bootstrap_stats_A, bootstrap_stats_B))
#
#   observed_stats = cbind(group = c(groups[[1]], groups[[2]]),
#                          rbind(observed_stats_A, observed_stats_B) %>%
#                            apply(c(1,2), as.numeric) %>% `row.names<-`(NULL)) %>%
#     data.frame()
#
#   return(list(P_values = p_values,
#               bootstrap_distribution_plot = plot,
#               observed_stats = observed_stats,
#               bootstrap_stats = bootstrap_stats,
#               bootstrap_replicates = bootstrap_replicates))
# }
