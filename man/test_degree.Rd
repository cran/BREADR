% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_degree.R
\name{test_degree}
\alias{test_degree}
\title{test_degree}
\usage{
test_degree(in_tibble, row, degree, verbose = TRUE)
}
\arguments{
\item{in_tibble}{a tibble that is the output of the callRelatedness() function.}

\item{row}{either the row number or pair name for which the posterior distribution is to be plotted.}

\item{degree}{the degree of relatedness to be tested.}

\item{verbose}{a logical (boolean) for whether all test output should be printed to screen.}
}
\value{
the associated p-value for the test
}
\description{
Test if a degree of relatedness is consistent with an observed PMR
}
\examples{
test_degree(relatedness_example, 1, 1)
}
