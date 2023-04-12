#' read_ind
#'
#' @param filename a IND text file.
#'
#' @return tibble with column headings: ind (CHR), sex (CHR), pop (CHR)
#' @export
#'
#' @examples
#' ind_snpfile <- system.file("extdata", "example.ind.txt", package = "BREADR")
#' read_ind(ind_snpfile)
read_ind <- function(filename){
    df <- data.table::fread(
      filename, header = FALSE,
      col.names=c("ind", "sex", "pop"),
      colClasses=list(character=1:3)
    )
    df <- tibble::as_tibble(df)
    return(df)
}

