## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BREADR)

## -----------------------------------------------------------------------------
ind_path <- "path/to/eigenstrat/indfile"
snp_path <- "path/to/eigenstrat/snpfile"
geno_path <- "path/to/eigenstrat/genofile"

## ----eval=FALSE---------------------------------------------------------------
# counts_example <- processEigenstrat(
#   indfile = ind_path,
#   snpfile = snp_path,
#   genofile = geno_path
# )

## ----eval=FALSE---------------------------------------------------------------
# counts_example <- processEigenstrat(
#   indfile = ind_path,
#   snpfile = snp_path,
#   genofile = geno_path,
#   outfile = "path_to_save_tsv"
# )

## ----example, eval = FALSE----------------------------------------------------
# library(BREADR)
# counts_example

## ----echo = FALSE-------------------------------------------------------------
counts_example

## ----eval = FALSE-------------------------------------------------------------
# relatedness_example <- callRelatedness(counts_example)
# relatedness_example

## ----echo = FALSE-------------------------------------------------------------
relatedness_example <- callRelatedness(counts_example)
relatedness_example

## ----fig.dim=c(7.2,5), warning=FALSE, message=FALSE, fig.cap = '*__Figure 1__: plotDOUGH results for the example data set, with the Between Groups population kept.*'----
indfile1 <- fs::path_package(package = "BREADR", "extdata/example.ind.txt")
plotDOUGH(
  counts_example,
  indfile = indfile1,
  nsnps_cutoff = 500
)

## ----echo=FALSE, fig.dim=c(7.2,5), warning=FALSE, message=FALSE, fig.cap = '*__Figure 2__: plotDOUGH results for the example data set where the PMR values have been artifically modified to yield different baseline (median) PMR values.*'----
indfile2 <- fs::path_package(package = "BREADR", "extdata/example2.ind.txt")
counts_example %>%
  dplyr::mutate(
    mismatch = dplyr::case_when(
      stringr::str_count(pair, 'Ind[123]') == 2 ~ floor(0.22 * nsnps),
      stringr::str_count(pair, 'Ind[456]') == 2 ~ floor(0.16 * nsnps),
      T ~ floor(0.19 * nsnps)
    ) +
      sample(10:20, dplyr::n(), replace = TRUE),
    pmr = mismatch / nsnps
  ) %>%
  plotDOUGH(
    indfile = indfile2,
    nsnps_cutoff = 500
  )

## ----warning=FALSE, message=FALSE, fig.dim=c(7.2,5), fig.cap = '*__Figure 3__: Pairwise diagnostic plot for relatedness call for Ind1 and Ind 2 called using the number option for row. (A) the distribution of PMR values given the number of overlapping SNPs (n=1518) for each degree of relatedness (filled colour), and the observed PMR with 95% confidence intervals. (B) relative posterior probabilities for all four possible assignments of degrees of relatedness. Note we added panel labels A and B using the labels parameter.*'----
plotLOAF(relatedness_example)

## ----fig.dim=c(7.2,5), warning=FALSE, fig.cap = '*__Figure 4__: Pairwise diagnostic plot for relatedness call for Ind1 and Ind 2 called using the text option for row. (A) the distribution of PMR values given the number of overlapping SNPs (n=1518) for each degree of relatedness (filled colour), and the observed PMR with 95% confidence intervals. (B) relative posterior probabilities for all four possible assignments of degrees of relatedness. Note we added panel labels A and B using the labels parameter.*'----
plotSLICE(relatedness_example, row = 1, labels = c('A', 'B'))

## ----fig.dim=c(7.2,5), warning=FALSE, fig.cap = '*__Figure 5__: Pairwise diagnostic plot for relatedness call for Ind1 and Ind 2 called using the text option for row. (A) the distribution of PMR values given the number of overlapping SNPs (n=1518) for each degree of relatedness (filled colour), and the observed PMR with 95% confidence intervals. (B) relative posterior probabilities for all four possible assignments of degrees of relatedness. Note we added panel labels A and B using the labels parameter.*'----
plotSLICE(relatedness_example, row = "Ind1 - Ind2", labels = c('A', 'B'))

## -----------------------------------------------------------------------------
test_degree(relatedness_example, row = 1, degree = 3, verbose = TRUE)

## ----fig.dim=c(7.2,5), warning=FALSE, message=FALSE, fig.cap = '*__Figure 6__: sensitivity analysis for the relationship between Ind1 and Ind2. The prior probabilities for all degrees of relatedness are tested at values between 0.05 and 0.5. Bars indicate the proportion of times a degree of relatedness was the maximum posterior relationship. Here we see that higher prior probabilities for Second-Degree and Unrelated yield highest posterior probability assignments for those degrees of relatedness, respectively. Hence, we must consider these prior probabilities carefully when assigning this relationship.*'----
BREADR::priorSensitivity(
  in_BREADR = relatedness_example,
  row = "Ind1 - Ind2",
  degrees = c(3, 4),
  grid_space = 0.05,
  maxPrior = 0.5
)

## ----fig.dim=c(7.2,3.5), warning=FALSE, message=FALSE,fig.cap = '*__Figure 7__: sensitivity analysis for the relationship between Ind3 and Ind4. The prior probabilities for all degrees of relatedness are tested at values between 0.05 and 0.5. Bars indicate the proportion of times a degree of relatedness was the maximum posterior relationship, and the dashed black line indicates 50%. This clearly shows that no matter the prior probability values for any of the degrees of relatedness, the highest posterior probability assignment is always Same/Twins.*'----
BREADR::priorSensitivity(
  in_BREADR = relatedness_example,
  row = "Ind3 - Ind4",
  degrees = c(1),
  grid_space = 0.025,
  maxPrior = 0.5
)

