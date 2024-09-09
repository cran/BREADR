#' process Eigenstrat data
#'
#'A function that takes paths to an eigenstrat trio (ind, snp and geno file) and
#'returns the pairwise mismatch rate for all pairs on a thinned set of SNPs.
#'Options include choosing thinning parameter, subsetting by population names,
#'and filtering out SNPs for which deamination is possible.
#'
#' @param indfile path to eigenstrat ind file
#' @param genofile path to eigenstrat geno file.
#' @param snpfile path to eigenstrat snp file.
#' @param filter_length the minimum distance between sites to be compared (to reduce the effect of LD).
#' @param pop_pattern a character vector of population names to filter the ind file if only some populations are to compared.
#' @param filter_deam a TRUE/FALSE for if C->T and G->A sites should be ignored.
#' @param outfile (OPTIONAL) a path and filename to which we can save the output of the function as a TSV,
#' if NULL, no back up saved. If no outfile, then a tibble is returned.
#' @param chromosomes the chromosome to filter the data on.
#' @param verbose controls printing of messages to console
#'
#' @return out_tibble: A tibble containing four columns:
# - pair [char]: the pair of individuals that are compared.
# - nsnps [numeric]: the number of overlapping snps that were compared for each pair.
# - mismatch [numeric]: the number of sites which did not match for each pair.
# - pmr [numeric]: the pairwise mismatch rate (mismatch/nsnps).
#' @export
#'
#' @examples
#' # Use internal files to the package as an example
#' indfile <- system.file("extdata", "example.ind.txt", package = "BREADR")
#' genofile <- system.file("extdata", "example.geno.txt", package = "BREADR")
#' snpfile <- system.file("extdata", "example.snp.txt", package = "BREADR")
#' processEigenstrat_old(
#' indfile, genofile, snpfile,
#' filter_length=1e5,
#' pop_pattern=NULL,
#' filter_deam=FALSE
#' )
processEigenstrat_old <- function(indfile, genofile, snpfile,
                              filter_length=NULL, pop_pattern=NULL, filter_deam=FALSE,
                              outfile=NULL, chromosomes=NULL, verbose = TRUE){
  # Check to see if eigenstrat files exist
  if(!file.exists(indfile)){
    stop(paste0('indfile ',indfile,' does not exist!'))
  }
  if(!file.exists(genofile)){
    stop(paste0('genofile ',genofile,' does not exist!'))
  }
  if(!file.exists(snpfile)){
    stop(paste0('snpfile ',snpfile,' does not exist!'))
  }

  # Just check to see if a length filter was included and is reasonable, and if not,
  # let the user know of the default and/or problem.
  if(is.null(filter_length)){
    filter_length <- 1e5
    if(verbose){
      cat(sprintf('No site distance filter.\nUsing default minimum of 1e5.\n'))
    }
  }
  if(!is.numeric(filter_length)|(filter_length<=0)){
    stop('filter_length must be a positive number.')
  }

  # Check if all pop_names are in the possible populations.
  # Change to delim so works in different countries.
  # ind_raw <- readr::read_delim(
  #   indfile,
  #   col_names=c('ind','sex','pop'),
  #   col_types=c('ccc'),
  #   delim = "\t"
  # )
  ind_raw <- read_ind(indfile)
  if(!all(pop_pattern%in%ind_raw$pop)){
    paste0('Following populations in pop_pattern do not exist in indfile: ',
           paste0(setdiff(pop_pattern,ind_raw$pop),collapse=', ')) %>%
      stop()
  }

  # Check that filter deam is a logical variable.
  if(!is.logical(filter_deam)){
    stop('filter_deam must be a logical (boolean) variable, i.e., TRUE or FALSE.')
  }

  # Check that the folder in which we plan to save the outfile actually exists.
  if(!is.null(outfile)){
    if(!file.exists(dirname(outfile))){
      stop(paste0('Output folder ',dirname(outfile),' does not exist!'))
    }
  }

  # Collect SNP details
  if(verbose){
    cat('Reading in SNP data.\n')
  }
  ## Read in SNP file
  snps <- read_snp(snpfile)
  ## Add relative position
  snps <-
    snps %>%
    dplyr::mutate(
      relative_pos=1:dplyr::n()
    )
  ## Mutate and filter for DEAM
  if(filter_deam){
    snps <-
      snps %>%
      dplyr::mutate(
        deam=dplyr::case_when(
          anc=='C'&der=='T' ~ T,
          anc=='G'&der=='A' ~ T,
          T ~ F)
      ) %>%
      dplyr::filter(!deam) %>%
      dplyr::select(-deam)
  }
  # if(filter_deam){
  #   snps <- readr::read_delim(
  #     snpfile,
  #     col_names=c('snp','chr','pos','site','anc','der'),
  #     col_types=c('cdddcc'),
  #     delim = "\t"
  #   ) %>%
  #     dplyr::mutate(relative_pos=1:dplyr::n()) %>%
  #     dplyr::mutate(deam=dplyr::case_when(anc=='C'&der=='T' ~ T,
  #                                  anc=='G'&der=='A' ~ T,
  #                                  T ~ F)) %>%
  #     dplyr::filter(!deam) %>%
  #     dplyr::select(-deam)
  # }else{
  #   snps <- readr::read_delim(
  #     snpfile,
  #     col_names=c('snp','chr','pos','site','anc','der'),
  #     col_types=c('cdddcc'),
  #     delim = "\t"
  #   ) %>%
  #     dplyr::mutate(relative_pos=1:dplyr::n())
  # }
  # Filter for requested chromosomes
  available_chr <- snps$chr %>%
    unique()
  if(!is.null(chromosomes)){
    if(!any(chromosomes%in%available_chr)){
      stop(paste0('Chromosome names do not match SNP file:\n',
                  'Requested: ',paste0(chromosomes,collapse=', '),'\n',
                  'Available: ',paste0(available_chr,collapse=', '),'\n'))
    }else{
      snps <- snps %>%
        dplyr::filter(chr%in%chromosomes)
    }
  }
  # Print information about chromosomes
  if(verbose){
    cat(paste0('Analysing chromosomes:\n',paste0(unique(snps$chr),collapse=', '),'\n'))
  }
  # Collect (potentially filtered) ind details
  if(is.null(pop_pattern)){
    ind <-  ind_raw %>%
      dplyr::mutate(row_number=1:dplyr::n())
  }else{
    ind <- ind_raw %>%
      dplyr::mutate(row_number=1:dplyr::n()) %>%
      dplyr::filter(pop%in%pop_pattern)
  }

  # Throw error if ind has no rows.
  n_ind <- nrow(ind)
  if(n_ind<=1){
    if(is.null(pop_pattern)){
      stop(paste0(n_ind,' individuals in final, filtered ind file. Check ind file.'))
    }else{
      stop(paste0(n_ind,' individuals in final, filtered ind file. Check ind file, and check pop_pattern.'))
    }
  }
  geno_list <- vector(mode='list',length=n_ind)

  if(verbose){
    cat('Starting to read in genotype data.\n')
  }
  pb1 = utils::txtProgressBar(min=1,max=length(geno_list),initial=1,style=3)
  for(i in 1:n_ind){
    geno_list[[i]] <- (system(paste0('cut -c ',ind$row_number[i],' ',genofile),intern=T) %>% dplyr::na_if("9"))[snps$relative_pos]
    utils::setTxtProgressBar(pb1,i)

  }
  if(verbose){
    cat('  Complete.\n\n')

    cat('Starting to compare genotypes and calculate PMR.\n')
  }
  counter <- 1
  out_tibble <- tibble::tibble(pair=rep('x',choose(n_ind,2)),nsnps=0,mismatch=0,pmr=0)
  pb2 = utils::txtProgressBar(min=1,max=choose(n_ind,2),initial=1,style=3)
  for(i in 1:(n_ind-1)){
    ind_i <- ind$ind[i]
    vi <- geno_list[[i]]
    for(j in (i+1):n_ind){
      ind_j <- ind$ind[j]
      vj <- geno_list[[j]]

      snps_filtered_prefiltered <- snps  %>%
        dplyr::mutate(row_num=1:dplyr::n()) %>%
        dplyr::filter(!is.na(vi),!is.na(vj))

      if(nrow(snps_filtered_prefiltered)>0){
        snps_filtered <- snps_filtered_prefiltered %>%
          dplyr::group_by(chr) %>%
          dplyr::mutate(keep=getfilter(site,filter_length)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(keep) %>%
          dplyr::select(-keep)

        out_tibble$pair[counter] <- paste0(ind_i,' - ',ind_j)
        out_tibble$nsnps[counter] <- (!is.na(vi[snps_filtered$row_num])&!is.na(vj[snps_filtered$row_num])) %>% sum()
        out_tibble$mismatch[counter] <- (vi[snps_filtered$row_num]!=vj[snps_filtered$row_num]) %>% sum(na.rm=T)
        out_tibble$pmr[counter] <-  out_tibble$mismatch[counter]/out_tibble$nsnps[counter]
      }else{
        out_tibble$pair[counter] <- paste0(ind_i,' - ',ind_j)
        out_tibble$nsnps[counter] <- 0
        out_tibble$mismatch[counter] <- 0
        out_tibble$pmr[counter] <-  NA
      }
      counter <- counter+1
      utils::setTxtProgressBar(pb2,counter)
    }
  }

  if(!is.null(outfile)){
    if(verbose){
      cat(sprintf('\n\nSaving processed data to %s.\n',outfile))
    }
    readr::write_delim(out_tibble,file=outfile,delim='\t')
  }
  if(verbose){
    cat('\nComplete.\n\n')
  }

  out_tibble %>% return()
}

