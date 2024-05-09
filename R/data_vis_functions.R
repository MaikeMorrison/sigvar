
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

    stop(paste0("There is at least one NA value in your signature activity matrix. The missing entries are found in the following positions: ",
                na.pos.format.string))
  }

  # check if matrix rows sum to 1, and give useful warnings if rounding is necessary
  sums <- rowSums(Q) %>% round(5)
  if (any(sums != 1)) {
    if (missing(rep)) {
      warning("At least one signature activity matrix has rows which do not sum to exactly 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row.")
    } else {
      warning(paste0(
        "At least one of the rows of signature activity matrix number ", rep,
        " (restricted to the last K columns) does not sum to 1. Rounding the sum of each row to 1 by dividing all entries by the sum of the row."
      ))
    }
    # Normalize each row of the matrix by dividing by the rowsums
    Q <- Q / sums
  }
  return(Q)
}




# # plot_signature_prop ----------------------------------------
#' Visualize a relative abundance matrix as a stacked bar plot.
#'
#' This function enables graphical visualization of a matrix of compostional data.
#'In the output plot, each vertical bar represents a single vector;
#' the height of each color in the bar corresponds to the abundance of each category
#' in that vector. Because this function produces a
#' ggplot object, its output can be modified using
#' standard ggplot2 syntax.
#'
#' @param relab_matrix A matrix with \code{I=nrow(relab_matrix)} rows, each containing \code{K=ncol(relab_matrix)} non-negative entries that sum to 1.
#' If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix and the number of entries
#' that sum to 1 (\code{K}) must be specified.
#' @param group Optional; a string specifying the name of the column that describes which group each row (sample) belongs to. Use if \code{matrices} is a single matrix containing multiple groups of samples you wish to compare.
#' @param time Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight FAVA by the distance between samples.
#' @param w Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.
#' @param K Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.
#' @param arrange Optional; controls horizontal ordering of samples and vertical ordering of categories.
#'   If \code{arrange = TRUE} or \code{arrange = "both"}, samples are ordered by the categories of greatest
#'   abundance and categories are ordered in decreasing abundance. If \code{arrange = "vertical"}, sample
#'   order is unchanged but categories are ordered in decreasing abundance. If \code{arrange = "horizontal"},
#'   samples are ordered by the most abundant categories, but category order is unchanged. If \code{arrange} is missing
#'   or \code{arrange = FALSE}, neither order is changed.
#' @return A ggplot object containing a bar plot visualization of the relative abundance matrix.
#' @examples
#'
#' # Make an example matrix of compositional data
#' # Each row is an individual. Rows sum to 1.
#' population_A = matrix(c(
#'     .5, .3, .2,
#'     .4, .2, .4,
#'     .5, .4, .1,
#'     .6, .1, .3,
#'     .2, 0, .8
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   )
#'
#'   plot_signature_prop(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = FALSE
#'               )
#'   plot_signature_prop(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = "horizontal"
#'               )
#'   plot_signature_prop(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = "vertical"
#'               )
#'    plot_signature_prop(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = TRUE  # could also be "both"
#'               )
#'
#'
#' # You can modify the plot as you would any ggplot2 object
#' plot_signature_prop(relab_matrix = population_A,
#'               K = 3, # How many categories per vector?
#'               arrange = TRUE
#'               ) +
#'   # Below are example, optional modifications to the default plot
#'   ggplot2::ggtitle("Population A") +
#'   ggplot2::scale_fill_brewer("Blues") +
#'   ggplot2::scale_color_brewer("Blues") +
#'   ggplot2::xlab("Individuals")
#'   # Note that both scale_fill and scale_color are needed to change the color of the bars.
#'
#'
#'   # Plot a dataset which has 2 populations
#'
#'   population_B = matrix(c(
#'     .9, 0, .1,
#'     .6, .4, 0,
#'     .7, 0, .3,
#'     .3, .4, .3,
#'     .5, .3, .2
#'   ),
#'   nrow = 5,
#'   byrow = TRUE
#'   )
#'
#'
#'   populations_AB = cbind(data.frame(c("A", "A", "A", "A", "A",
#'                                      "B", "B", "B", "B", "B")),
#'                          rbind(population_A, population_B))
#'   colnames(populations_AB) = c("population", "category_1", "category_2", "category_3")
#'
#'
#'  plot_signature_prop(relab_matrix = populations_AB, group = "population")
#'  plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "vertical")
#'  plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "horizontal")
#'  plot_signature_prop(relab_matrix = populations_AB, group = "population", arrange = "both")
#'
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @export
plot_signature_prop <- function(relab_matrix, group = NULL, time = NULL, w = NULL, K = NULL, arrange = FALSE) {

  relab_checker_out = relab_checker(relab = relab_matrix, K = K, group = group, time = time)

  K = ncol(relab_checker_out$relab_matrix)

  # Repeat rows to account for time or weight (w) if provided,
  # otherwise return relab_matrix unaltered
  relab_edited = relab_sample_weighter(relab = relab_matrix, K = K, time = time, w = w, group = group)

  # Re-arrange rows or columns as specified by arrange
  # otherwise return relab_edited unaltered
  relab_edited = arrange_categories(relab_matrix = relab_edited,
                                    arrange = arrange,
                                    K = K, group = group, time = time)


  # Generate the data to plot
  relab_plot = dplyr::mutate(relab_edited, ID = 1:nrow(relab_edited), .before = 1)


  start =  2 + (!is.null(group)) + (!is.null(time))

  relab_plot_long = tidyr::pivot_longer(relab_plot, cols = start:ncol(relab_plot))


  relab_plot_long$ID = factor(relab_plot_long$ID , levels = unique(relab_plot$ID), ordered = TRUE)
  relab_plot_long$name = factor(relab_plot_long$name, levels = colnames(relab_plot)[start:ncol(relab_plot)] %>% rev, ordered = TRUE)

  if(!is.null(group)){

    if(is.null(relab_checker_out$group)){
      stop("The group provided is not a column name in relab_matrix. Please provide a valid group.")
    }

    ggplot2::ggplot(data = relab_plot_long, ggplot2::aes(fill = .data$name,
                                                         color = .data$name,
                                                         y = .data$value,
                                                         x = .data$ID)) +
      ggplot2::geom_bar(position = "stack", stat = "identity",
                        width = 1) + ggplot2::theme_void() + ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 12)) +
      ggplot2::facet_wrap(~ group, scales = "free_x")
  }else{

    ggplot2::ggplot(data = relab_plot_long, ggplot2::aes(fill = .data$name,
                                                         color = .data$name,
                                                         y = .data$value,
                                                         x = .data$ID)) +
      ggplot2::geom_bar(position = "stack", stat = "identity",
                        width = 1) + ggplot2::theme_void() + ggplot2::ylab("") +
      ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(size = 12))

  }
}





# plot_dots -----------------------------------------------------------------
#' Plot the mean mutational signature activities of each group of samples
#'
#' This function generates a dot plot with signatures along the y-axis,
#' group/populations along the x-axis, dot size corresponding to the fraction
#' of samples in each population containing a signature, and dot color
#' corresponding to the mean relative activity of mutations to that
#' signature.
#'
#' @param sig_activity A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents a sample, each column represents a mutational signature, and each entry represents the abundance of that signature in the sample. If \code{sig_activity} contains any metadata, it must be on the left-hand side of the matrix, the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of each row must sum to 1.
#' @param group A string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each sample belongs to. Default is the first column of \code{sig_activity}: \code{colnames(sig_activity)[1]}.
#' @param K Optional; the number of signatures in the matrix. Default is
#'   \code{K=ncol(sig_activity)-1-as.numeric(!missing(facet))}, the number of columns in
#'    \code{sig_activity} minus either 1 (since one column must specify the group) or 2 if
#'   \code{facet} is provided, since one column must specify the variable on
#'   which to facet.
#' @param median Optional; specify `median = TRUE` if you would like to color points by the median signature activity across samples, rather than the mean. Default value is `median = FALSE`.
#' @param normalized Optional; specify `normalized = FALSE` if you would like to plot absolute signature activities rather than relative signature activities. I.e., if your `sig_activity` matrix contains rows representing absolute signature activities (i.e., they do not necessarily sum to 1) and you wish that these activities are not automatically normalized. This determines the color scale for the resulting plot. Default is `normalized = TRUE`.
#' @param facet Optional; a string specifying the name of the column by which
#'   you would like to facet your plot.
#' @param pivot Optional; set \code{pivot=TRUE} if you would like to plot groups
#'   on the y-axis and signatures on the x-axis. Default is \code{pivot=FALSE}.
#' @param max_dotsize Optional; a number specifying the maximum size for each dot.
#' @param threshold Optional; a number between 0 and 1. For each signature, the size of the associated dot corresponds to the proportion of tumors with that signature at an abundance exceeding this threshold.  Default is \code{threshold = 0}.
#' @return A ggplot object containing a dot plot visualization of the mean
#'   mutational signature activities
#' @examples
#'   # Make an example matrix.
#'   # Each row is a sample. Rows sum to 1.
#' sig_activity = matrix(c(
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
#'   colnames(sig_activity) = c("pop", "SBS1", "SBS2", "SBS3")
#'   sig_activity = dplyr::mutate(data.frame(sig_activity), dplyr::across(c(SBS1, SBS2, SBS3),
#'                                                  as.numeric))
#'
#' plot_dots(sig_activity,
#'  K = 3, # How many categories per vector?
#'   group = "pop"
#' )
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_dots <- function(sig_activity, group = colnames(sig_activity)[1],
                      K=ncol(sig_activity)-1-as.numeric(!missing(facet)),
                      max_dotsize = 5,
                      pivot = FALSE,
                      median = FALSE,
                      normalized = TRUE,
                      facet, threshold = 0) {

  # Satisfy R cmd check
  Signature <- . <- Mean_activity <- Proportion_present <- NULL


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

  facet_true = !missing(facet)
  # will there be a few facet panels, or many? used to determine # of columns
  facets_few = ifelse(facet_true, length(unique(unlist(sig_activity[facet])))<=4, FALSE)

  if(length(K)>0){ if(K>(ncol(sig_activity)-1) ) warning(paste0("K too large, not enough columns in K; K reduced to ncol(sig_activity)-1=",ncol(sig_activity)-1))}

  signatures = colnames(sig_activity)[(ncol(sig_activity) - K + 1):ncol(sig_activity)]

  #  if(facet_true){
  #   signatures = colnames(sig_activity)[colnames(sig_activity)!=group & colnames(sig_activity)!=facet ][1:min(K,ncol(sig_activity)-1)]
  # }else{
  #   signatures = colnames(sig_activity)[colnames(sig_activity)!=group][1:min(K,ncol(sig_activity)-1)]
  # }


  if(facet_true){
    sig_activity_sigs = cbind(data.frame("group" = sig_activity[[group]],
                                         "facet" = sig_activity[[facet]]),
                              (if(normalized){Q_checker(sig_activity %>% dplyr::select(-c(group, facet)), K)}else{
                                sig_activity[,(ncol(sig_activity) - K + 1):ncol(sig_activity)]})) %>%
      `colnames<-`(c("group", "facet", signatures))

    sig_activity_present = sig_activity_sigs %>%
      dplyr::group_by(group, facet) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(signatures),
                                     function(x) mean(x > threshold))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(signatures), names_to = "Signature",
                          values_to = "Proportion_present")


  }else{
    sig_activity_sigs = cbind(data.frame("group" = sig_activity[[group]]),
                              (if(normalized){Q_checker(sig_activity %>% dplyr::select(-group), K)}else{
                                sig_activity[,(ncol(sig_activity) - K + 1):ncol(sig_activity)]})) %>%
      `colnames<-`(c("group", signatures))

    sig_activity_present = sig_activity_sigs %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(signatures),
                                     function(x) mean(x > threshold))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(signatures), names_to = "Signature",
                          values_to = "Proportion_present")
  }





  if(median){
    sig_activity_means = sig_activity_sigs %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(signatures), function(col) median(col[col>0]))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(signatures),
                          names_to = "Signature",
                          values_to = "Mean_activity")
  }else{
    sig_activity_means = sig_activity_sigs %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(signatures), function(col) mean(col[col>0]))) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(signatures),
                          names_to = "Signature",
                          values_to = "Mean_activity")
  }




  plot_data = dplyr::inner_join(sig_activity_present, sig_activity_means) %>%
    dplyr::mutate(Signature = factor(Signature, ordered = TRUE,
                                     levels = (signatures))) %>%
    {if(facet_true){dplyr::left_join(.,
                                     sig_activity %>%
                                       dplyr::transmute(group = get(group),
                                                        facet = get(facet)) %>%
                                       dplyr::distinct())
    }else{.}}




  ggplot2::ggplot(plot_data,
                  ggplot2::aes(y = Signature, x = group,
                               color = Mean_activity,
                               size = Proportion_present)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(guide = ggplot2::guide_legend(title.position = "top",
                                                                 direction = "horizontal"),
                                   #max_size = max_dotsize,
                                   limits = c(threshold,1), range = c(-1, max_dotsize),
                                   breaks = c(0.5, 1),
                                   name = "Proportion of\ntumors with\nsignature") +
    {if(median)ggplot2::scale_color_gradientn(guide = ggplot2::guide_colorbar(title.position = "top",
                                                                              barwidth = 4,
                                                                              direction = "horizontal"),
                                              colours = c("#E81F27", "#881F92", "#2419F9"),
                                              limits = c(0,ifelse(normalized, 1, max(plot_data$Mean_activity, na.rm = TRUE))),
                                              breaks = signif(seq(from = 0,
                                                                  to = ifelse(normalized, 1,
                                                                              max(plot_data$Mean_activity, na.rm = TRUE)),
                                                                  length.out = 3), 2),
                                              name = ifelse(normalized,
                                                            "Median relative\nactivity in\ntumors with\nsignature",
                                                            "Median activity\nin tumors with\nsignature"))} +
    {if(!median)ggplot2::scale_color_gradientn(guide = ggplot2::guide_colorbar(title.position = "top",
                                                                               barwidth = 4,
                                                                               direction = "horizontal"),
                                               colours = c("#E81F27", "#881F92", "#2419F9"),
                                               limits = signif(c(0, ifelse(normalized, 1,
                                                                               max(plot_data$Mean_activity, na.rm = TRUE))), 2),
                                               breaks = signif(seq(from = 0,
                                                            to = ifelse(normalized, 1,
                                                                        max(plot_data$Mean_activity, na.rm = TRUE)),
                                                            length.out = 3), 2),
                                               name = ifelse(normalized,
                                                             "Mean relative\nactivity in\ntumors with\nsignature",
                                                             "Mean activity\nin tumors with\nsignature"))} +
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

    {if(facet_true & pivot & facets_few)ggforce::facet_col(dplyr::vars(facet), scales = "free_y", space = "free")}+ # This used to have ncol = 1 when it was facet_wrap
    {if(facet_true & pivot & !facets_few)ggplot2::facet_grid(.~facet, scales = "free_y", space = "free_y")}+
    {if(facet_true & !pivot)ggplot2::facet_grid(.~facet, scales = "free_x", space = "free_x")}+

    {if(facet_true)ggplot2::theme(strip.background = ggplot2::element_blank(),
                                  strip.text = ggplot2::element_text(size = 12))}+

    ggplot2::xlab(group)
}




# plot_SBS_spectrum -----------------------------------------------------------------
#' Plot an SBS mutational spectrum
#'
#' @param SBS_table A matrix with rows corresponding to the 96 single-base substitutions and columns corresponding to distinct mutational spectra you wish to plot.
#' @return A ggplot2 bar chart depicting SBS mutational spectra
#' @examples
#' SBS_table = dplyr::select(COSMIC3.3.1_SBS, SBS1, SBS5)
#' plot_SBS_spectrum(SBS_table)
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
plot_SBS_spectrum <- function(SBS_table) {

  # Satisfy R cmd check
  Type <- name <- Relative_abundance <- Sub <- NULL


  # CHECKS ON SBS_table ------------------------------------------------------

  # does SBS_table have 96 rows?
  if(nrow(SBS_table) != 96){stop("SBS_table must have exactly 96 rows, one for each single-base substitution.")}


  # does SBS_table have only numeric columns?
  col_count = ncol(SBS_table)
  SBS_table = dplyr::select_if(SBS_table, is.numeric)

  if(col_count != ncol(SBS_table)){warning(paste0("The ",
                                                  col_count - ncol(SBS_table),
                                                  " column(s) containing non-numeric values were ommitted."))}

  # does SBS_table have only columns that sum to 1?
  if(any(round(colSums(SBS_table),5) != 1)){
    SBS_table = apply(SBS_table, 2, function(col) col/sum(col))

    warning("At least one column did not sum to 1. The columns have each been divided by their sum so that they now sum to 1.")
  }


  sbs <- sigvar::COSMIC3.3.1_SBS %>%
    dplyr::select(Type)
  sbs$Sub = stringr::str_split(sbs$Type, "\\[|\\]", simplify = TRUE)[,2]

  sbs$Context =  stringr::str_split(sbs$Type, "\\[|\\]|>") %>% lapply(function(x) paste0(x[1], x[2], x[4])) %>% unlist

  sub_pal = c("#02BCED", "#010101", "#E22926",
              "#CAC8C9", "#A0CE62", "#ECC6C5") %>%
    `names<-`(unique(sort(sbs$Sub)))

  sbs = dplyr::left_join(sbs,
                         data.frame(Sub = names(sub_pal),
                                    color = sub_pal))


  strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = sub_pal, color = sub_pal),
                               text_x = ggh4x::elem_list_text(color = c("black", "white", "white", rep("black", 3))))

  plot_data_wide = cbind(sbs, SBS_table) %>%
    dplyr::mutate(name = glue::glue("<b style='color:#BEBEBE'>{stringr::str_sub(Context, 1,1)}<b style='color:{color}'>{stringr::str_sub(Context, 2, 2)}<b style='color:#BEBEBE'>{stringr::str_sub(Context, 3,3)}"), .before = 5)

  plot_data_long = tidyr::pivot_longer(plot_data_wide, cols = colnames(SBS_table),
                                       names_to = "Spectrum", values_to = "Relative_abundance")

  ggplot2::ggplot(plot_data_long, ggplot2::aes(x = name, y = Relative_abundance, fill = Sub)) +
    ggplot2::geom_bar(stat = "identity") +
    ggh4x::facet_grid2(Spectrum ~ Sub, strip = strip, scales = "free")  +
    ggplot2::theme_minimal() +
    # scale_x_discrete(expand = c(0, 0))+
    # scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0,
      # size = 6, family = "mono"),
      axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 0, size = 4,
                                             family = "mono", margin=ggplot2::margin(0)),
      # axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 8),
      strip.text.y = ggplot2::element_text(size = 12),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank())+
    ggplot2::scale_fill_manual(values = sub_pal) +
    ggplot2::ylab("Relative abundance") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3)


}

