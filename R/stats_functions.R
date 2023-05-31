#' test_covariable
#' Test the statistical association between the variability of a Q matrix and a covariable
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
#' test_covariable(Qlist[[1]],1)
#' @importFrom readr read_tsv
test_covariable <- function(Q,w){
  test = NA
  return(test)
}
