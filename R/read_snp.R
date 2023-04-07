#' split line
#'
#' takes a line for a SNP file and splits into parts.
#' @param x line from SNP file
#'
#' @return tibble with 6 columns.
#' @export
#'
#' @examples
#' split_line("1_14.570829090394763     1        0.000000              14 A X")
#' split_line("rs3094315	1	0.0	752566	G	A")
split_line <- function(x){
  x <- stringr::str_trim(x)
  df <- stringr::str_split(x, " +|\t") %>% purrr::pluck(1)
  stopifnot(length(df) == 6)
  df <- tibble::tibble(
    snp = df[1],
    chr = as.numeric(df[2]),
    pos = as.numeric(df[3]),
    site = as.numeric(df[4]),
    anc = df[5],
    der = df[6]
    )
  return(df)
}
#' read_snp
#'
#' @param filename a SNP text file.
#'
#' @return tibble with column headings: snp (CHR), chr (DBL),
#' pos (DBL), site (DBL), anc (CHR), and der (CHR).
#' @export
#'
#' @examples
#' std_snpfile <- system.file("extdata", "example.snp.txt", package = "BREADR")
#' broken_snpfile <- system.file("extdata", "broken.snp.txt", package = "BREADR")
#' read_snp(std_snpfile)
#' read_snp(broken_snpfile)
read_snp <- function(filename){
    df <- data.table::fread(
      filename,
      col.names=c("snp", "chr", "pos", "site", "anc", "der"),
      colClasses=list(character=c(1,2,5,6))
    )
    df <- tibble::as_tibble(df)
    return(df)
}
# pacman::p_load(tidyverse, targets)
# read_snp("inst/extdata/broken.snp.txt")
# read_snp("inst/extdata/example.snp.txt")
