
# Q_checker ------------------------------------------------------------------
# An internal function to check if Q matrices are up to spec and fix any issues automatically.
# Q - a Q matrix provided as input to another function
# K - the number of ancestral clusters
# rep - an optional parameter used if Q matrix is one of a list, in order to provide more useful warning messages
Q_checker <- function(Q, K, rep) {
  # Check if Q matrix is within a list, and extract if needed
  if (is.list(Q) && !is.data.frame(Q) && !is.array(Q)) {
    Q <- Q[[1]]
  }
  # Check if Q matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last K columns.
  if (ncol(Q) > K) {
    Q <- Q[, (ncol(Q) - K + 1):ncol(Q)]
  }

  # convert Q matrix entries to numbers
  Q <- data.matrix(Q)

  # Name Q matrix columns q1, q2, ..., qK
  colnames(Q) <- paste0("q",1:K)

  # Check if Q matrix has any missing values, and give warning if necessary
  if(any(is.na(Q))){
    # Identify location of missing entries
    na.pos <- sapply(which(is.na(Q)),
                     function(index) c(index %% nrow(Q), ceiling(index/nrow(Q)))) %>%
      t()
    # Format missing entries as a string
    na.pos.format <- list()
    for(row in 1:nrow(na.pos)){
      na.pos.format[row] <- paste0("(", na.pos[row,1], ", ", na.pos[row,2], ")")
    }
    na.pos.format.string <- as.character(na.pos.format) %>% paste(collapse = ", ")

    stop(paste0("There is at least one NA value in your Q matrix. The missing entries are found in the following positions: ",
                na.pos.format.string))
  }

  # check if matrix rows sum to 1, and give useful warnings if rounding is necessary
  sums <- rowSums(Q) %>% round(5)
  if (any(sums != 1)) {
    if (missing(rep)) {
      warning("At least one Q matrix has rows which do not sum to exactly 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row.")
    } else {
      warning(paste0(
        "At least one of the rows of Q matrix number ", rep,
        " (restricted to the last K columns) does not sum to 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row."
      ))
    }
    # Normalize each row of the matrix by dividing by the rowsums
    Q <- Q / sums
  }
  return(Q)
}


# plot_dots -----------------------------------------------------------------
#' Plot the mean mutational signature attributions of each group of samples
#'
#' This function generates a dot plot with signatures along the y-axis,
#' group/populations along the x-axis, dot size corresponding to the fraction
#' of samples in each population containing a signature, and dot color
#' corresponding to the mean relative attribution of mutations to that
#' signature.
#'
#' @param Q A dataframe, matrix, or array representing mutatational signature
#'   relative attributions.
#'   Each row represents a sample. The first \code{ncol(Q)-K} columns may
#'   contain other information about the sample, and must contain the grouping
#'   variable. When restricted to the last \code{K} columns, the rows of this
#'   matrix must sum to 1.
#' @param K The number of signatures in the matrix. Each vector must
#'   must have \code{K} categories. Default is the number of columns in \code{Q}
#'   minus 1 (since one column must specify the group).
#' @param group A string specifying the name of the column that describes which
#'   group/population each sample belongs to. Default is the first column name.
#' @return A ggplot object containing a dor plot visualization of the mean
#'   mutational signature attributions
#' @examples
#'   # Make an example matrix.
#'   # Each row is a sample. Rows sum to 1.
#' Q = matrix(c(
#'     "A", .4, .2, .4,
#'     "A", .5, .3, .2,
#'     "A", .5, .4, .1,
#'     "B", .6, .1, .3,
#'     "B", .6, .3, .1
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   )
#'
#'   colnames(Q) = c("pop", "SBS1", "SBS2", "SBS3")
#'   Q = dplyr::mutate(data.frame(Q), dplyr::across(c(SBS1, SBS2, SBS3),
#'                                                  as.numeric))
#'
#' plot_dots(Q,
#'  K = 3, # How many categories per vector?
#'   group = "pop"
#' )
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_dots <- function(Q, K=ncol(Q)-1, group = colnames(Q)[1]) {

  signatures = colnames(Q)[(ncol(Q)-K + 1):ncol(Q)]

  Q_sigs = cbind(data.frame("group" = Q[[group]]),
                 Q_checker(Q, K)) %>%
    `colnames<-`(c("group", signatures))

  Q_present = Q_sigs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(signatures),
                                   function(x) sum(x > 0)/dplyr::n())) %>%
    tidyr::pivot_longer(cols = signatures, names_to = "Signature",
                        values_to = "Proportion\nof tumors\nwith sig.")


  Q_means = Q_sigs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(signatures), mean)) %>%
    tidyr::pivot_longer(cols = signatures, names_to = "Signature",
                        values_to = "Mean\nrelative\ncontribution")

  dplyr::inner_join(Q_present, Q_means) %>%
  ggplot2::ggplot(ggplot2::aes(y = Signature, x = group,
             color = `Mean\nrelative\ncontribution`,
             size = `Proportion\nof tumors\nwith sig.`)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_area(max_size=8,
                    limits = c(0,1), breaks = c(0, 0.5, 1)) +
    ggplot2::scale_color_gradientn(colours = c("#E81F27", "#881F92", "#2419F9"),
                          limits = c(0,1), breaks = c(0, 0.5, 1)) +
    ggplot2::theme_bw()  +
    # theme(axis.title = element_blank(),
    #       axis.text.y = ggtext::element_markdown(size = 12),
    #       axis.text.x = element_text(size = 12, color = c("#199CC9", "#D86600")))+
    ggplot2::guides(color = ggplot2::guide_colourbar(barheight = 3)) +
    ggplot2::xlab(group)

}
