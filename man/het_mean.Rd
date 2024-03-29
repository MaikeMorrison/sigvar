% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fava_functions.R
\name{het_mean}
\alias{het_mean}
\title{Compute the mean Gini-Simpson index of the rows in a matrix of compositional vectors}
\usage{
het_mean(relab_matrix, K = NULL, S = NULL, w = NULL, time = NULL, group = NULL)
}
\arguments{
\item{relab_matrix}{A matrix or data frame with rows containing non-negative entries that sum to 1. Each row represents
a sample, each column represents a category, and each entry represents the abundance of that category in the sample.
If \code{relab_matrix} contains any metadata, it must be on the left-hand side of the matrix,
the right \code{K} entries of each row must sum to 1, and \code{K} must be specified. Otherwise, all entries of
each row must sum to 1.}

\item{K}{Optional; an integer specifying the number of categories in the data. Default is \code{K=ncol(relab_matrix)}.}

\item{S}{Optional; a K x K similarity matrix with diagonal elements equal to 1 and off-diagonal elements between 0 and 1. Entry \code{S[i,k]} for \code{i!=k} is the similarity between category and \code{i} and category \code{k}, equalling 1 if the categories are to be treated as identical and equaling 0 if they are to be treated as totally dissimilar. The default value is \code{S = diag(ncol(relab_matrix))}.}

\item{w}{Optional; a vector of length \code{I} with non-negative entries that sum to 1. Entry \code{w[i]} represents the weight placed on row \code{i} in the computation of the mean abundance of each category across rows. The default value is \code{w = rep(1/nrow(relab_matrix), nrow(relab_matrix))}.}

\item{time}{Optional; a string specifying the name of the column that describes the sampling time for each row. Include if you wish to weight Fst by the distance between samples.}

\item{group}{Optional; a string (or vector of strings) specifying the name(s) of the column(s) that describes which group(s) each row (sample) belongs to. Use if \code{relab_matrix} is a single matrix containing multiple groups of samples you wish to compare.}
}
\value{
A numeric value between 0 and 1.
}
\description{
This function computes the mean Gini-Simpson index, a statistical measure of variability also known as the Gini-Simpson index, of a set of vectors of non-negative entries which sum to 1. The function returns a number between 0 and 1 which quantifies the mean variability of the vectors. Values of 0 are achieved when each vector is a permutation of (1,0,..., 0). The value approaches 1 as the number of categories K increases when the vectors are equal to (1/K, 1/K, ..., 1/K).
}
\examples{
# To compute the mean Gini-Simpson index of
# the following compositional vectors...
q1 = c(1,   0,   0,   0)
q2 = c(0.5, 0.5, 0,   0)
q3 = c(1/4, 1/4, 1/4, 1/4)
q4 = c(0,   0,   1,   0)

# we could compute the mean manually:
mean(sapply(list(q1, q2, q3, q4), het))

# Or we could use het_mean:
relative_abundances = matrix(c(q1, q2, q3, q4),
                  byrow = TRUE, nrow = 4)

het_mean(relative_abundances)

# Incoporating weights:

# Compute mean Gini-Simpson index ignoring
# rows 2 and 3
row_weights = c(0.5, 0, 0, 0.5)
het_mean(relative_abundances, w = row_weights)

# Compute mean Gini-Simpson index assuming that
# categories 1 and 2 are identical:
similarity_matrix = diag(4)
similarity_matrix[1,2] = 1
similarity_matrix[2,1] = 1
het_mean(relative_abundances, S = similarity_matrix)

# Assume categories 1 and 2 are identical AND
# ignore rows 2 and 4:
row_weights = c(0.5, 0, 0.5, 0)
het_mean(relative_abundances, w = row_weights, S = similarity_matrix)
}
