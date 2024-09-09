#' callRelatedness
#'
#'
#' A function that takes PMR observations, and (given a prior distribution for
#' degrees of relatedness) returns the posterior probabilities of all pairs of
#' individuals being (a) the same individual/twins, (b) first-degree related,
#' (c) second-degree related or (d) "unrelated" (third-degree or higher).
#' The highest posterior probability degree of relatedness is also returned as
#' a hard classification.
#' Options include setting the background relatedness (or using the sample
#' median), a minimum number of overlapping SNPs if one uses the sample median
#' for background relatedness, and a minimum number of overlapping SNPs
#' for including pairs in the analysis.
#'
#' @param pmr_tibble a tibble that is the output of the processEigenstrat
#' function.
#' @param class_prior the prior probabilities for same/twin, 1st-degree,
#' 2nd-degree, unrelated, respectively.
#' @param average_relatedness a single numeric value, or a vector of numeric
#' values, to use as the average background relatedness. If NULL, the sample
#' median is used.
#' @param median_co if average_relatedness is left NULL, then the minimum cutoff
#' for the number of overlapping snps to be included in the median calculation
#' is 500.
#' @param filter_n the minimum number of overlapping SNPs for which pairs are
#' removed from the entire analysis. If NULL, default is 1.
#'
#' @return   results_tibble: A tibble containing 13 columns:
#' - row: The row number
#' - pair: the pair of individuals that are compared.
#' - relationship: the highest posterior probability estimate of the degree of relatedness.
#' - pmr: the pairwise mismatch rate (mismatch/nsnps).
#' - sd: the estimated standard deviation of the pmr.
#' - mismatch: the number of sites which did not match for each pair.
#' - nsnps: the number of overlapping snps that were compared for each pair.
#' - ave_re;: the value for the background relatedness used for normalisation.
#' - Same_Twins: the posterior probability associated with a same individual/twins classification.
#' - First_Degree: the posterior probability associated with a first-degree classification.
#' - Second_Degree: the posterior probability associated with a second-degree classification.
#' - Unrelated: the posterior probability associated with an unrelated classification.
#' - BF: A strength of confidence in the Bayes Factor associated with the highest posterior probability classification compared to the 2nd highest. (No longer included)
#' @export
#'
#' @examples
#' callRelatedness(counts_example,
#'   class_prior=rep(0.25,4),
#'   average_relatedness=NULL,
#'   median_co=5e2,filter_n=1
#' )
callRelatedness <- function(pmr_tibble,
                    class_prior=rep(0.25,4),
                    average_relatedness=NULL,
                    median_co=5e2,filter_n=1
){

  # Check pmr_tibble is of the form required:
  cols_check <- c('pair','nsnps','mismatch','pmr')
  if((ncol(pmr_tibble)!=4)|!(all(names(pmr_tibble)==cols_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(cols_check,collapse=', '),'.'))
  }

  # Check that the prior distribution makes sense:
  if(any(class_prior<0)|any(class_prior>1)){
    stop('Posterior probabilities for degrees of relatedness must be values between 0 and 1.')
  }else if(sum(class_prior)!=1){
    stop('Posterior probabilities for degrees of relatedness must sum to 1.')
  }

  # Check that the average relatedness makes sense:
  if(!is.null(average_relatedness)){
    # if((average_relatedness<=0)|average_relatedness>=1){ OLD VERSION
    if(any(average_relatedness<=0) | any(average_relatedness>=1)){
      stop('The average relatedness must be a value between 0 and 1.')
    }
  }

  # Check median estimator cut off:
  if((median_co<=0)){
    stop('The median cut off number of overlapping snps must be greater than 0.')
  }else if(median_co>=max(pmr_tibble$nsnps,na.rm=T)){
    stop('The median cut off number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }

  # Check median estimator cut off:
  if((filter_n<=0)){
    stop('The cut off for the number of overlapping snps must be greater than 0.')
  }else if(filter_n>=max(pmr_tibble$nsnps,na.rm=T)){
    stop('The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }

  # local variables
  class_vec <- c('Same_Twins','First_Degree','Second_Degree','Unrelated')

  # If no user-defined value for average_relatedness is given, use the median (above nsnps cut off)
  if(is.null(average_relatedness)){
    M <- pmr_tibble %>%
      dplyr::filter(nsnps>median_co) %>%
      dplyr::pull(pmr) %>%
      stats::median(na.rm=T)
    filter_tibble <- pmr_tibble %>%
      dplyr::filter((nsnps>=filter_n)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(ave_rel=M)
  }else{
    M <- average_relatedness
    filter_tibble <- pmr_tibble %>%
      dplyr::mutate(ave_rel=M) %>%
      dplyr::filter((nsnps>=filter_n))
  }

  results_tibble <- filter_tibble %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Same_Twins=(stats::dbinom(mismatch,nsnps,0.5*ave_rel,log=T)+log(class_prior[1])),
                  First_Degree=(stats::dbinom(mismatch,nsnps,(1-0.5^2)*ave_rel,log=T)+log(class_prior[2])),
                  Second_Degree=(stats::dbinom(mismatch,nsnps,(1-0.5^3)*ave_rel,log=T)+log(class_prior[3])),
                  Unrelated=weightedBinom(mismatch,nsnps,ave_rel,log=T)+log(class_prior[4])) %>%
    # dplyr::mutate(BF=diff(sort(c(Same_Twins,First_Degree,Second_Degree,Unrelated))[c(3,4)])) %>%
    dplyr::mutate(normConst=matrixStats::logSumExp(c(Same_Twins,First_Degree,Second_Degree,Unrelated))%>%exp()) %>%
    dplyr::mutate(relationship=class_vec[which.max(c(Same_Twins,First_Degree,Second_Degree,Unrelated))]) %>%
    dplyr::ungroup() %>%
    # dplyr::mutate(BF=dplyr::case_when(BF <= log(10^0.5) ~ 'Weak Evidence',
    #                            BF <= log(10^1) ~ 'Substantial Evidence',
    #                            BF <= log(10^1.5) ~ 'Strong Evidence',
    #                            BF <= log(10^2) ~ 'Very Strong Evidence',
    #                            T ~ 'Decisive')) %>%
    dplyr::mutate(Same_Twins=exp(Same_Twins)/normConst,
                  First_Degree=exp(First_Degree)/normConst,
                  Second_Degree=exp(Second_Degree)/normConst,
                  Unrelated=exp(Unrelated)/normConst) %>%
    dplyr::mutate(sd=sqrt(pmr*(1-pmr)/nsnps),
                  relationship=factor(relationship,levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated')),
                  row=1:dplyr::n()) %>%
    dplyr::select(row,pair,relationship,pmr,sd,mismatch,nsnps,ave_rel,everything(),-normConst)
    return(results_tibble)
}

