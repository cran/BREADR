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

## ----warning=FALSE, message=FALSE---------------------------------------------
plotLOAF(relatedness_example)

## -----------------------------------------------------------------------------
plotSLICE(relatedness_example, row = 1)

## -----------------------------------------------------------------------------
plotSLICE(relatedness_example, row = "Ind1 - Ind2")

## -----------------------------------------------------------------------------
test_degree(relatedness_example, 1, 3, verbose = TRUE)

