
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
  # Check if K is too large
  if(K>ncol(Q)){stop("K is larger than the number of columns in the provided matrix.")}
  # Check if Q matrix is in STRUCTURE/ADMIXTURE output form, or if it only contains columns of ancestry coefficients
  # If it is in STRUCTURE/ADMIXTURE output form, extract the ancestry coefficients--the last K columns.
  if (ncol(Q) > K) {
    Q <- Q[, (ncol(Q) - K + 1):ncol(Q)]
  }

  # convert Q matrix entries to numbers
  Q <- data.matrix(Q)

  # Name Q matrix columns q1, q2, ..., qK
  #colnames(Q) <- paste0("q",1:K)

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


# plot_signature_prop -----------------------------------------------------------------
#' Plot the signature proportion for a set of samples
#'
#'
#' @param Q A dataframe, matrix, or array representing mutational signature
#'   relative contributions.
#'   Each row represents a sample.
plot_signature_prop <- function(Q){
  # convert to long format
  Q_long <- Q %>% mutate(Sample=paste0("P",1:nrow(Q))) %>%
    tidyr::pivot_longer(-Sample,names_to = "Signature",values_to = "Proportion") %>%
    mutate(Signature = factor(Signature, ordered = TRUE, levels=colnames(Q)[order(colMeans(Q),decreasing = T)]) )

  # order samples
  Sample_order = Q_long %>% dplyr::group_by(Signature) %>% dplyr::arrange(Proportion) %>% dplyr::pull(Sample) %>% unique()
  Q_long = Q_long %>% dplyr::mutate(Sample = factor(Sample,levels=Sample_order))

  # create color palette
  sig_palette = sigFAVA::sbs_palette[levels(Q_long$Signature)[levels(Q_long$Signature)%in%names(sigFAVA::sbs_palette)]]
  sig_palette = c(sig_palette,
                  PNWColors::pnw_palette("Sailboat",length(levels(Q_long$Signature))-length(sig_palette) ) %>%
                    `names<-`(sort(levels(Q_long$Signature)[!levels(Q_long$Signature)%in%names(sigFAVA::sbs_palette)]) ))
  sig_palette = sig_palette[order(names(sig_palette))]

  ggplot2::ggplot(Q_long,
       ggplot2::aes(x = Sample, y = Proportion,
           fill = Signature, color = Signature)) +
  ggplot2::geom_bar(stat = "identity") +
  #facet_wrap(~ Group, scales = "free_x", ncol = 4)+
  scale_color_manual(values = sig_palette) + #, breaks = rev(SBS_sigs)[1:10]) +
  scale_fill_manual( values = sig_palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(#legend.position = "bottom",
    panel.grid = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 11),
    axis.text = ggplot2::element_blank()) +
    ggplot2::ylab("") + ggplot2::xlab("Sample")
}

# plot_dots -----------------------------------------------------------------
#' Plot the mean mutational signature contributions of each group of samples
#'
#' This function generates a dot plot with signatures along the y-axis,
#' group/populations along the x-axis, dot size corresponding to the fraction
#' of samples in each population containing a signature, and dot color
#' corresponding to the mean relative contribution of mutations to that
#' signature.
#'
#' @param Q A dataframe, matrix, or array representing mutational signature
#'   relative contributions.
#'   Each row represents a sample. The first \code{ncol(Q)-K} columns may
#'   contain other information about the sample, and must contain the grouping
#'   variable. When restricted to the last \code{K} columns, the rows of this
#'   matrix must sum to 1.
#' @param group A string specifying the name of the column that describes which
#'   group/population each sample belongs to. Default is the first column name.
#' @param K Optional; the number of signatures in the matrix. Default is
#'   \code{K=ncol(Q)-1-as.numeric(!missing(facet))}, the number of columns in
#'    \code{Q} minus either 1 (since one column must specify the group) or 2 if
#'   \code{facet} is provided, since one column must specify the variable on
#'   which to facet.
#' @param facet Optional; a string specifying the name of the column by which
#'   you would like to facet your plot.
#' @param pivot Optional; set \code{pivot=TRUE} if you would like to plot groups
#'   on the y-axis and signatures on the x-axis. Default is \code{pivot=FALSE}.
#' @param max_dotsize Optional; a number specifying the maximum size for each dot.
#' @return A ggplot object containing a dor plot visualization of the mean
#'   mutational signature contributions
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
plot_dots <- function(Q, group = colnames(Q)[1],
                      K=ncol(Q)-1-as.numeric(!missing(facet)),
                      max_dotsize = 5,
                      pivot = FALSE,
                      facet) {

  facet_true = !missing(facet)
  # will there be a few facet panels, or many? used to determine # of columns
  facets_few = ifelse(facet_true, length(unique(unlist(Q[facet])))<4, FALSE)

  if(length(K)>0){ if(K>(ncol(Q)-1) ) warning(paste0("K too large, not enough columns in K; K reduced to ncol(Q)-1=",ncol(Q)-1))}
  signatures = colnames(Q)[colnames(Q)!=group][1:min(K,ncol(Q)-1)]

  Q_sigs = cbind(data.frame("group" = Q[[group]]),
                 Q_checker(Q %>% dplyr::select(-group), K)) %>%
    `colnames<-`(c("group", signatures))

  Q_present = Q_sigs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(signatures),
                                   function(x) sum(x > 0)/dplyr::n())) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(signatures), names_to = "Signature",
                        values_to = "Proportion_present")

  Q_means = Q_sigs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(signatures), function(col) mean(col[col>0]))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(signatures),
                        names_to = "Signature",
                        values_to = "Mean_contribution")


  plot_data = dplyr::inner_join(Q_present, Q_means) %>%
    dplyr::mutate(Signature = factor(Signature, ordered = TRUE,
                                     levels = (signatures))) %>%
    {if(facet_true){dplyr::left_join(.,
                                     Q %>%
                                       dplyr::transmute(group = get(group),
                                                        facet = get(facet)) %>%
                                       dplyr::distinct())
    }else{.}}

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(y = Signature, x = group,
                               color = Mean_contribution,
                               size = Proportion_present)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(guide = ggplot2::guide_legend(title.position = "top",
                                                        direction = "horizontal"),
                                   #max_size = max_dotsize,
                                   limits = c(0,1), range = c(-1, max_dotsize),
                                   breaks = c(0.5, 1),
                                   name = "Proportion of\ntumors with\nsignature") +
    ggplot2::scale_color_gradientn(guide = ggplot2::guide_colorbar(title.position = "top",
                                                          barwidth = 4,
                                                          direction = "horizontal"),
                                   colours = c("#E81F27", "#881F92", "#2419F9"),
                                   limits = c(0,1), breaks = c(0, 0.5, 1),
                                   name = "Mean relative\ncontribution in\ntumors with\nsignature") +
    ggplot2::theme_bw()  +
    # ggplot2::guides(color = ggplot2::guide_colourbar(barheight = 3)) +
    # {if(!pivot)ggplot2::guides(color = ggplot2::guide_colourbar(barheight = 3))}  +

    {if(pivot)ggplot2::coord_flip()} +
    #{if(pivot)
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #}+

    # {if(pivot & !facet_true)ggplot2::theme(legend.position = "top")}+
    # {if(pivot & !facet_true)ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 3))}+
    # {if(pivot & facet_true)ggplot2::guides(color = ggplot2::guide_colourbar(barheight = 3))}+

    {if(facet_true & pivot & facets_few)ggplot2::facet_wrap(~facet, scales = "free_y", ncol = 1)}+
    {if(facet_true & pivot & !facets_few)ggplot2::facet_wrap(~facet, scales = "free_y")}+
    {if(facet_true & !pivot)ggplot2::facet_wrap(~facet, scales = "free_x")}+

    {if(facet_true)ggplot2::theme(strip.background = ggplot2::element_blank(),
                                  strip.text = ggplot2::element_text(size = 12))}+

    ggplot2::xlab(group)
}
