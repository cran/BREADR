#' relatedness_example
#'
#' this is an example of the tibble made by callRelatedness()
#'
#' @format ## `relatedness_example`
#' A data frame with 15 rows and 13 columns:
#' \describe{
#' \item{row}{The row number}
#' \item{pair}{the pair of individuals that are compared.}
#' \item{relationship}{the highest posterior probability estimate of the degree of relatedness.}
#' \item{pmr}{the pairwise mismatch rate (mismatch/nsnps).}
#' \item{sd}{the estimated standard deviation of the pmr.}
#' \item{mismatch}{the number of sites which did not match for each pair.}
#' \item{nsnps}{the number of overlapping snps that were compared for each pair.}
#' \item{ave_re}{ the value for the background relatedness used for normalisation.}
#' \item{Same_Twins}{the posterior probability associated with a same individual/twins classification.}
#' \item{First_Degree}{the posterior probability associated with a first-degree classification.}
#' \item{Second_Degree}{the posterior probability associated with a second-degree classification.}
#' \item{Unrelated}{the posterior probability associated with an unrelated classification.}
#' \item{BF}{A strength of confidence in the Bayes Factor associated with the highest posterior probability classification compared to the 2nd highest.}
#' }
"relatedness_example"
