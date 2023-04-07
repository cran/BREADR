#' test_degree
#'
#' Test if a degree of relatedness is consistent with an observed PMR
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param row either the row number or pair name for which the posterior distribution is to be plotted.
#' @param degree the degree of relatedness to be tested.
#' @param verbose a logical (boolean) for whether all test output should be printed to screen.
#'
#' @return the associated p-value for the test
#' @export
#'
#' @examples
#' test_degree(relatedness_example, 1, 1)
test_degree <- function(in_tibble,row,degree,verbose=TRUE){
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((ncol(in_tibble)!=12)|(!all(names(in_tibble)==col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }

  # check row is a numeric or char variable, and if char, find the row.
  if(is.character(row)){
    rowtmp <- row
    row <- which(row==in_tibble$pair)
    if(length(row)==0){
      stop(paste0('There is no pair called: ',rowtmp))
    }
  }else if(is.numeric(row)&(row>nrow(in_tibble))){
    stop('row number requested was higher than the number of rows in input tibble!')
  }else if(!is.numeric(row)){
    stop('row must be a number or a pair name.')
  }

  if(!is.numeric(degree)){
    stop('Degree must be numeric.')
  }
  # if(degree>10|degree<0|(degree%%1!=0)){ # Removed constraint that integer
  if(degree>10|degree<0){
    stop('Degree must be a value k such that 0<k<10.')
  }

  M <- in_tibble$ave_rel[row]
  pk <- (1-0.5^(degree+1))*M
  N <- in_tibble$nsnps[row]
  x <- in_tibble$mismatch[row]
  pair <- in_tibble$pair[row]
  # pval <- 2*stats::pbinom(x,N,pk,ifelse(x>pk*N,F,T))
  pval <- min(
    c(
      2*stats::pbinom(x,N,pk,T),
      2*stats::pbinom(x,N,pk,F)
    )
  )
  # CI <- x/N+c(-1,1)*1.96*sqrt(x/N*(1-x/N)/N)
  pvaltext <- ifelse(pval<0.0001,
                     sprintf('p-value          : %s\n',format(pval,scientific=T)),
                     sprintf('p-value          : %.4f\n',pval))


  degree_string <- dplyr::case_when(degree==1 ~ "1st",
                             degree==2 ~ "2nd",
                             degree==3 ~ "3rd",
                             T ~ paste0(degree,"th"))
  if(verbose){
    cat(sprintf('Testing H0       : "%s" are %s-degree relatives.\n',pair,degree_string))
    cat(sprintf('Expected PMR     : %.4f\n',pk))
    cat(sprintf('Observed PMR     : %.4f\n',x/N))
    cat(sprintf('Estimated degree : %.4f\n',convertP(x/N,M)))
    cat(pvaltext)
    cat(ifelse(pval>=0.05,'Decision         : Retain H0\n','Decision         : Reject H0\n'))
  }
  return(pval)
}
