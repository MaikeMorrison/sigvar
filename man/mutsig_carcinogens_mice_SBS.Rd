% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mutsig_carcinogens_mice_SBS}
\alias{mutsig_carcinogens_mice_SBS}
\title{Single Base Substitution mutational signatures for 181 mice exposed to known or suspected carcinogens}
\format{
\subsection{\code{mutsig_carcinogens_mice_SBS}}{

A tibble with 181 rows and 46 columns:
\describe{
\item{Sample}{Identifier of the tumor}
\item{mSBS1, mSBS5, ..., mSBS_N3}{Number of mutations for each signature}
\item{Tissue}{Tissue of origin}
\item{DIAGNOSIS}{Type of cancer}
\item{dose}{Treatment dose}
\item{chemical}{Known or suspected carcinogen chemical administered}
}
}
}
\source{
\url{https://github.com/team113sanger/mouse-mutatation-signatures/blob/6a00d910df40d178c373ac4d57849918e9814951/figure1/mexposure.rds}
}
\usage{
mutsig_carcinogens_mice_SBS
}
\description{
Single Base Substitution mutational signature attributions
from Riva et al. Nat Genet 52, 1189–1197 (2020). https://doi.org/10.1038/s41588-020-0692-4
}
\keyword{datasets}
