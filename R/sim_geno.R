#' sim_geno
#'
#' Simulated geno file of eigenstrat format
#' 
#' @param n_ind number of individuals
#' @param n_snp number of SNPs
#' @param filename filename of export 
#'
#' @return NULL exports a file
#' @export
#'
#' @examples
#' \dontrun{
#' sim_geno(10, 5, "geno.txt")
#' }
sim_geno <- function(n_ind, n_snp, filename){
  N <- n_ind * n_snp
  M <- matrix(
    sample(c(0, 1, 2, 9), size = N, replace = TRUE), 
    nrow = n_snp, ncol = n_ind
  )
  MASS::write.matrix(M, file = filename, sep = "")
}
# pacman::p_load(conflicted, tidyverse, targets)
# sim_geno(10, 5, "geno.txt") |> print()