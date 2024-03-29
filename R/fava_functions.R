
# Function to process relab matrix and return correct grouping columns, w, K, S, etc.

process_relab <- function(relab_matrix,
                          K = NULL,
                          S = NULL,
                          w = NULL,
                          time = NULL,
                          group = NULL){

  # Satisfy R cmd check
  grouping_var_multiple <- NULL

  # Arrange data by time if not already
  if(!is.null(time)){
    if(!is.null(group)){
      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[group[[1]]]], .data[[time]])
    } else{
      relab_matrix = dplyr::arrange(.data = relab_matrix, .data[[time]])
    }}

  # Define K if not provided
  if(is.null(K)){
    K = ncol(relab_matrix) - (length(group)) - (!is.null(time))
  }

  # Define S if not provided
  if(is.null(S)){
    S = diag(K)
  }

  # Modify relab_matrix if there are multiple groups
  multiple_groups = FALSE
  relab_grouping_vars = NULL
  if(length(group) > 1){
    relab_grouping_vars = dplyr::select(relab_matrix, dplyr::all_of(group))

    relab_grouping_vars$grouping_var_multiple = apply(relab_grouping_vars, 1, paste, collapse = "_")

    relab_matrix = dplyr::mutate(relab_matrix,
                                 grouping_var_multiple = relab_grouping_vars$grouping_var_multiple,
                                 .before = 1)

    if(any(table(relab_matrix$grouping_var_multiple)<2)){
      ignore = names(which(table(relab_matrix$grouping_var_multiple)<2))
      warning("Only analyzing combinations of grouping variables with at least two samples. Ignoring the following combinations of grouping variables: ", paste(ignore, collapse = "  "))
      relab_matrix = dplyr::filter(relab_matrix,
                                   grouping_var_multiple %in%
                                     names(which(table(relab_matrix$grouping_var_multiple) >= 2)))
    }

    group = "grouping_var_multiple"

    multiple_groups = TRUE
  }

  if(length(group) == 1){
    if(any(table(relab_matrix[[group]])<2)){
      ignore = names(which(table(relab_matrix[[group]])<2))
      warning("Only analyzing groups with at least two samples. Ignoring the following groups: ", paste(ignore, collapse = "  "))
      relab_matrix = relab_matrix[relab_matrix[[group]] %in% names(which(table(relab_matrix[[group]]) >= 2)),]
    }
  }


  relab_matrix_clean = relab_checker(relab_matrix,
                                     K = K,
                                     group = group,
                                     time = time)

  if(any(is.na(relab_matrix_clean$group))){
    stop("There is at least one NA in a grouping column.")
  }


  if(!is.null(time)){
    w = time_weights(times = relab_matrix_clean$time, group = relab_matrix_clean$group)
  }

  w_default = FALSE
  if((!is.null(group)) & (is.null(w))){
    w_default = TRUE
    w = rep(1/table(relab_matrix_clean$group), table(relab_matrix_clean$group))
  }

  if(is.null(w)){
    w = rep(1/nrow(relab_matrix), nrow(relab_matrix))
  }

  S = as.matrix(S)
  S_checker(S = S, K = K, relab_matrix = relab_matrix)

  if(!missing(w) & length(w) != nrow(relab_matrix)){
    stop("Length of w must equal number of rows of relab_matrix.")
  }

  return(list(K = K,
              S = S,
              w = w,
              w_default = w_default,
              time = time,
              group = group,
              multiple_groups = multiple_groups,
              relab_matrix_clean = relab_matrix_clean,
              relab_grouping_vars = relab_grouping_vars))

}



# het -----------------------------------------------------------------
#' Compute the Gini-Simpson index of a compositional vector
#'
#' This function computes the Gini-Simpson index, a statistical measure of variability also known as the Gini-Simpson index, of vector of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the variability of the vector. Values of 0 are achieved when the vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param q A vector with \code{K=length(q)} non-negative entries that sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=length(q)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(q))}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute unweighted Gini-Simpson index:
#' het(q = c(0.4, 0.3, 0.3))
#'
#' # Compute Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(3)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' het(q = c(0.4, 0.3, 0.3), S = similarity_matrix)
#' @export
het <- function(q, K = length(q), S = diag(K)){
  S = as.matrix(S)
  S_checker(S = S, K = K)
  q = unlist(q)[(length(q)-K+1):length(q)]
  q = as.numeric(q)

  if(round(sum(q), 4) != 1){stop("Vector does not sum to 1")}
  if(length(q) == 1){q = c(1,0)}

  1 - sum(q * c(S %*% q))
}

# het_mean -----------------------------------------------------------------
#' Compute the mean Gini-Simpson index of the rows in a matrix of compositional vectors
#'
#' This function computes the mean Gini-Simpson index, a statistical measure of variability also known as the Gini-Simpson index, of a set of vectors of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the mean variability of the vectors. Values of 0 are achieved when each vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vectors are equal to (1/K, 1/K, ..., 1/K).
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight Fst by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the mean Gini-Simpson index of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' mean(sapply(list(q1, q2, q3, q4), het))
#'
#' # Or we could use het_mean:
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' het_mean(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute mean Gini-Simpson index ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' het_mean(relative_abundances, w = row_weights)
#'
#' # Compute mean Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' het_mean(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' het_mean(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
het_mean <- function(relab_matrix,
                              K = NULL,
                              S = NULL,
                              w = NULL,
                              time = NULL,
                              group = NULL
                    # K = ncol(relab_matrix),
                    # w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                    # S = diag(ncol(relab_matrix))
                    ){

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars

  if(is.null(group)){
    het_mean_fast(relab_matrix_clean$relab_matrix, K, w, S)
  } else{
    gs_list = c()
    for(subgroup in unique(relab_matrix_clean$group)){

      include = relab_matrix_clean$group == subgroup

      relab_sub = relab_matrix_clean$relab_matrix[include,]

      if(w_default){
        w_sub = w[names(w) == subgroup]
      } else{
        w_sub = w[include]
      }

      gs_list = c(gs_list,
                  het_mean_fast(relab_sub, K, w_sub, S))
    }
    gs_df = data.frame(unique(relab_matrix_clean$group), gs_list)
    colnames(gs_df) = c(group, "het_mean")
    # names(gs_list) = unique(relab_matrix_clean$group)

    if(multiple_groups){
      gs_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), gs_df)
    }
    return(gs_df)
  }

  # OLD VERSION BELOW
  # # K, S, w, time, and group are optional arguments
  #
  # # S = as.matrix(S)
  # I = nrow(relab_matrix)
  # # S_checker(S = S, K = K)
  #
  #
  # if(!missing(w) && length(w) != nrow(relab_matrix)){
  #   stop("Length of w must equal number of rows of relab_matrix.")
  # }
  #
  # relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  #
  # # Average Gini-Simpson index of each of the I subpopulations
  # sum(w *sapply(1:I, function(i){ het(q = unlist(relab_matrix[i,]), S = S) }))

}

# het_pooled -----------------------------------------------------------------
#' Compute the pooled Gini-Simpson index of the rows in a matrix of compositional vectors
#'
#' This function computes the Gini-Simpson index of a "pooled" vector equal to \code{colMeans(relab_matrix)}. Values of 0 are achieved when this pooled vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when this pooled vector is equal to (1/K, 1/K, ..., 1/K).
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight Fst by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # To compute the pooled Gini-Simpson index of
#' # the following compositional vectors...
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#'
#' # we could compute the mean manually:
#' qPooled = (q1 + q2 + q3 + q4)/4
#' het(qPooled)
#'
#' # Or we could use het_pooled:
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' het_pooled(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute pooled Gini-Simpson index ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' het_pooled(relative_abundances, w = row_weights)
#'
#' # Compute pooled Gini-Simpson index assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' het_pooled(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' het_pooled(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
het_pooled <- function(relab_matrix,
                                K = NULL,
                                S = NULL,
                                w = NULL,
                                time = NULL,
                                group = NULL
                      # K = ncol(relab_matrix),
                      # w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                      # S = diag(ncol(relab_matrix))
                      ){

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars

  if(is.null(group)){
    het_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S)
  } else{
    gs_list = c()
    for(subgroup in unique(relab_matrix_clean$group)){

      include = relab_matrix_clean$group == subgroup

      relab_sub = relab_matrix_clean$relab_matrix[include,]

      if(w_default){
        w_sub = w[names(w) == subgroup]
      } else{
        w_sub = w[include]
      }

      gs_list = c(gs_list,
                  het_pooled_fast(relab_sub, K, w_sub, S))
    }
    gs_df = data.frame(unique(relab_matrix_clean$group), gs_list)
    colnames(gs_df) = c(group, "het_pooled")
    # names(gs_list) = unique(relab_matrix_clean$group)

    if(multiple_groups){
      gs_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), gs_df)
    }
    return(gs_df)
  }

  # OLD VERSION BELOW
  # # w and S are optional arguments
  #
  # I = nrow(relab_matrix)
  #
  # # S = as.matrix(S)
  # # S_checker(S = S, K = K)
  #
  # if(missing(S)){
  #   S = diag(K)
  # }
  #
  # if(missing(w)){
  #   w = rep(1, I)/I
  # }
  #
  #
  # if(!missing(w) && length(w) != nrow(relab_matrix)){
  #   stop("Length of w must equal number of rows of relab_matrix.")
  # }
  #
  # relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  #
  # het(q = colSums(sweep(x = relab_matrix, MARGIN = 1, w, `*`)), S = S)

}



# fst_norm -----------------------------------------------------------------
#' Compute the normalized Fst of a matrix of compositional vectors
#'
#' This function computes the normalized Fst given the number of rows and the mean abundance of the most abundant category.
#' We employ the normalization employed in the [FSTruct package](https://github.com/MaikeMorrison/FSTruct) by
#' [Morrison, Alcala, and Rosenberg (2020)](https://doi.org/10.1111/1755-0998.13647).
#'
#' @param relab_matrix A matrix with \code{I=nrow(relab_matrix)} rows, each containing \code{K=ncol(relab_matrix)} non-negative entries that sum to 1.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix and the number of entries
#' that sum to 1 (\code{K}) must be specified.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the weighted fst of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fst_norm(relative_abundances)
#' @export
fst_norm <- function(relab_matrix, K = ncol(relab_matrix)){
  relab_matrix = relab_matrix[,(ncol(relab_matrix)-K+1):ncol(relab_matrix)]
  I <- nrow(relab_matrix) # Here, I is the number of individuals (number of subpopulations)
  p <- colSums(relab_matrix) # yields vector of summed allele frequencies across populations
  sig1 <- max(p)
  J <- ceiling(1 / sig1)
  sig1.frac <- sig1 - floor(sig1)

  if (sig1 == I) {
    fstMax <- 0
    fst <- 0
    ratio <- 0
  } else {
    if (sig1 <= 1) {
      fstMax <- ((I - 1) * (1 - sig1 * (J - 1) * (2 - J * sig1))) /
        (I - (1 - sig1 * (J - 1) * (2 - J * sig1)))
    } else {
      fstMax <- (I * (I - 1) - sig1^2 + floor(sig1) - 2 * (I - 1) * sig1.frac + (2 * I - 1) * sig1.frac^2) / (I * (I - 1) - sig1^2 - floor(sig1) + 2 * sig1 - sig1.frac^2)
    }

    fst <-
      (sum(relab_matrix^2) / I - sum(colSums(relab_matrix / I)^2)) /
      (1 - sum(colSums(relab_matrix / I)^2))

    fst / fstMax
  }
}



# FUNCTIONS FOR WEIGHTINGS

# S_checker ------------------------------------------------------------------
# An internal function to check if similarity matrices are up to spec and fix any issues automatically.
# S - a similarity matrix
# K - the number of taxa
S_checker <- function(S, K, relab_matrix = NULL) {

  if(!isSymmetric(S)){
    warning("S is not symmetric.")
  }
  if(any(diag(S)!=1)){
    stop("All diagonal elements of S must be equal to 1.")
  }
  if(nrow(S) != K | ncol(S) != K){
    stop("S must have K rows and K columns.")
  }
  if(any(S>1)){
    stop("All elements of S must be less than or equal to 1.")
  }
  if(any(S<0)){
    stop("All elements of S must be positive.")
  }
  if(!is.null(relab_matrix)){
    taxa_names = colnames(relab_matrix)[(ncol(relab_matrix)-K + 1):ncol(relab_matrix)]
    if(any(taxa_names!=colnames(S))){
      warning("The column names of the similarity matrix S do not match the names of the K categories in relab_matrix.")
    }
    if(any(taxa_names!=rownames(S))){
      warning("The row names of the similarity matrix S do not match the names of the K categories in relab_matrix.")
    }
  }
}


# time_weights -----------------------------------------------------------------
#' Compute a normalized weighting vector based on a vector of sampling times.
#'
#' This function takes a vector of sampling times, \eqn{t = (t_1, t_2, \ldots, t_I)}{latex}
#' and computes a normalized vector which can be used to weight each sample based on
#' the time between the subsequent and the preceding samples. The weighting vector \eqn{w}
#' is defined such that each entry, \eqn{w_i = d_i / 2T}, where \eqn{T=t_I - t_1} and
#' \eqn{d_i = t_{i+1} - t_{i-1}} for \eqn{i} not equal to 1 or I. \eqn{d_1 = t_2-t_1} and \eqn{d_I = t_I-t_{I-1}}.
#'
#' @param times A numeric vector of sampling times. Each entry must be non-negative and
#' greater than the previous entry.
#' @param group Optional; a character vector specifying the group identity of each
#' sampling time. Use if there are samples from multiple replicates or subjects
#' in one data set.
#' @examples

#' time_vector = c(1, 8, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
#'                 32, 33, 34, 35, 36, 37, 38, 39, 44, 50, 57, 64)
#'
#' time_weights(times = time_vector)
#' @export
time_weights <- function(times, group = NULL){

  if(is.null(group)){
    I = length(times)
    T = times[I] - times[1]

    if(I<2){
      stop("times must have length greater than 1.")
    }
    if(any(sapply(2:I, function(i) times[i] - times[i-1]) <= 0)){
      stop("times must be increasing. Each entry must be greater than the previous entry.")
    }

    di <- c(times[2] - times[1])
    for(i in 2:(I-1)){
      di[i] = times[i+1] - times[i-1]
    }
    di[I] = times[I] - times[I-1]
    return(di/(2*T))
  }else{

    wi = c()
    for(name in unique(group)){
      time_name = times[which(group == name)]

      I = length(time_name)
      T = time_name[I] - time_name[1]

      if(I<2){
        stop("Within each group, times must have length greater than 1.")
      }
      if(any(sapply(2:I, function(i) time_name[i] - time_name[i-1]) <= 0)){
        stop("Within each group, times must be increasing. Each entry must be greater than the previous entry.")
      }

      wi <- c(wi, (time_name[2] - time_name[1])/(2*T))
      for(i in 2:(I-1)){
        wi = c(wi, (time_name[i+1] - time_name[i-1])/(2*T))
      }
      wi = c(wi, (time_name[I] - time_name[I-1])/(2*T))
    }
    return(wi)
  }

}



# fst -----------------------------------------------------------------
#' Compute the Fst of a matrix of compositional vectors
#'
#' This function computes the population genetic statistic Fst on any matrix with rows that sum to 1. Values of 0 are achieved when each row is a permutation of (1,0,..., 0) and at least two categories have non-zero abundance across all rows. The value equals 1 when each row is identical.
#'
#' @param relab_matrix A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
#' a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
#' the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
#' each row must sum to 1.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param S Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight Fst by the distance between samples.
#' @param group Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.
#' @param normalized Optional; should normalized Fst be used? Default is \code{normalized = FALSE}; use \code{normalized = TRUE} to compute normalized Fst. Fst can only be normalized if it is not weighted.
#' @returns A numeric value between 0 and 1.
#' @examples
#' # Compute the Fst of
#' # the following compositional vectors:
#' q1 = c(1,   0,   0,   0)
#' q2 = c(0.5, 0.5, 0,   0)
#' q3 = c(1/4, 1/4, 1/4, 1/4)
#' q4 = c(0,   0,   1,   0)
#
#' relative_abundances = matrix(c(q1, q2, q3, q4),
#'                   byrow = TRUE, nrow = 4)
#'
#' fst(relative_abundances)
#'
#'# Incoporating weights:
#'
#' # Compute fst ignoring
#' # rows 2 and 3
#' row_weights = c(0.5, 0, 0, 0.5)
#' fst(relative_abundances, w = row_weights)
#'
#' # Compute fst assuming that
#' # categories 1 and 2 are identical:
#' similarity_matrix = diag(4)
#' similarity_matrix[1,2] = 1
#' similarity_matrix[2,1] = 1
#' fst(relative_abundances, S = similarity_matrix)
#'
#' # Assume categories 1 and 2 are identical AND
#' # ignore rows 2 and 4:
#' row_weights = c(0.5, 0, 0.5, 0)
#' fst(relative_abundances, w = row_weights, S = similarity_matrix)
#' @export
fst <- function(relab_matrix,
                 K = NULL,
                 S = NULL,
                 w = NULL,
                 time = NULL,
                 group = NULL,
                 normalized = FALSE){

  if(normalized == TRUE && any(!sapply(list(time, w, S), is.null))){
    stop("Fst can be either normalized or weighted, but not both. Please specify `normalized = TRUE` if you wish to compute normalized Fst OR provide the weighting parameters w and/or S.")
  }

  process_out = process_relab(relab_matrix = relab_matrix, K = K, S = S, w = w, time = time, group = group)

  K = process_out$K
  S = process_out$S
  w = process_out$w
  time = process_out$time
  group = process_out$group

  w_default = process_out$w_default
  multiple_groups = process_out$multiple_groups
  relab_matrix_clean = process_out$relab_matrix_clean
  relab_grouping_vars = process_out$relab_grouping_vars



  if(is.null(group)){
    if(normalized){
      fst_norm(relab_matrix = relab_matrix_clean$relab_matrix)
    }else{
      (het_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S) -
         het_mean_fast(relab_matrix_clean$relab_matrix, K, w, S))/
        het_pooled_fast(relab_matrix_clean$relab_matrix, K, w, S)
    }} else{
      fst_list = c()
      for(subgroup in unique(relab_matrix_clean$group)){

        include = relab_matrix_clean$group == subgroup

        relab_sub = relab_matrix_clean$relab_matrix[include,]

        if(w_default){
          w_sub = w[names(w) == subgroup]
        } else{
          w_sub = w[include]
        }

        fst_list = c(fst_list,
                      ifelse(normalized,
                             fst_norm(relab_matrix = relab_sub),
                             (het_pooled_fast(relab_sub, K, w_sub, S) -
                                het_mean_fast(relab_sub, K, w_sub, S))/
                               het_pooled_fast(relab_sub, K, w_sub, S)))
      }
      fst_df = data.frame(unique(relab_matrix_clean$group), fst_list)
      colnames(fst_df) = c(group, "Fst")
      # names(fst_list) = unique(relab_matrix_clean$group)

      if(multiple_groups){
        fst_df = dplyr::left_join(dplyr::distinct(relab_grouping_vars), fst_df)
      }
      return(fst_df)
    }
}




# fast versions of functions to use in fst function:
het_fast <- function(q, S = diag(length(q)), K = length(q)){

  1 - sum(q * c(S %*% q))
}

het_mean_fast <- function(relab_matrix,
                              K = ncol(relab_matrix),
                              w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                              S = diag(ncol(relab_matrix))){
  I = nrow(relab_matrix)

  # Average Gini-Simpson index of each of the I subpopulations
  sum(w *sapply(1:I, function(i){ het_fast(q = unlist(relab_matrix[i,]), S = S) }))

}

het_pooled_fast <- function(relab_matrix,
                                K = ncol(relab_matrix),
                                w = rep(1/nrow(relab_matrix), nrow(relab_matrix)),
                                S = diag(ncol(relab_matrix))){

  I = nrow(relab_matrix)

  het_fast(q = colSums(sweep(x = relab_matrix, MARGIN = 1, w, `*`)), S = S)

}
