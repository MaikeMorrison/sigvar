% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mutsig_carcinogens_mice_bootstrap_p_vals}
\alias{mutsig_carcinogens_mice_bootstrap_p_vals}
\title{Bootstrap p-values comparing SVA results for carcinogen-exposed murine tumors to spontaneous tumors}
\format{
\subsection{\code{mutsig_carcinogens_mice_bootstrap_p_vals}}{

A tibble with 29 rows and 7 columns:
\describe{
\item{group}{Organ of tumor and exposure status}
\item{chemical}{Exposure substance}
\item{Tissue}{Tissue of origin}
\item{group_2}{Tumor/exposure reference each tumor/exposure group was compared to (either "Spontaneous_liver"  or"Spontaneous_lung" )}
\item{across_sample_heterogeneity}{P-values comparing values of across-sample heterogeneity between group and group_2}
\item{mean_within_sample_diversity}{P-values comparing values of within-sample diversity between group and group_2}
\item{pooled_diversity}{P-values comparing the pooled diversity between group and group_2}
}
}
}
\usage{
mutsig_carcinogens_mice_bootstrap_p_vals
}
\description{
Data on single Base Substitution mutational signature attributions
from Riva et al. Nat Genet 52, 1189–1197 (2020). https://doi.org/10.1038/s41588-020-0692-4
}
\keyword{datasets}
