#' plotLOAF
#'
#' Plots all (sorted by increasing value) observed PMR values with
#' maximum posterior probability classifications represented by colour
#' and shape.
#' Options include a cut off for the minimum number of overlapping SNPs,
#' the max number of pairs to plot and x-axis font size.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param nsnps_cutoff the minimum number of overlapping SNPs for which
#' pairs are removed from the plot. If NULL, default is 500.
#' @param N the number of (sorted by increasing PMR) pairs to plot.
#' Avoids plotting all pairs (many of which are unrelated).
#' @param fntsize the fontsize for the x-axis names.
#' @param verbose if TRUE, then information about the plotting process is sent
#' to the console
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' relatedness_example
#' plotLOAF(relatedness_example)
plotLOAF <- function(in_tibble,nsnps_cutoff=NULL,N=NULL,fntsize=7, verbose = TRUE){

  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((ncol(in_tibble)!=12)|(!all(names(in_tibble)==col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }

  # Test that the overlapping snp cut off is sensible.
  if(is.null(nsnps_cutoff)){
    nsnps_cutoff <- 5e2
    if(verbose){
      cat(sprintf('No minimum number of overlapping SNPs given.\nUsing default minimum of 500.\n'))
    }
  }
  if((nsnps_cutoff<=0)){
    stop('The cut off for the number of overlapping snps must be greater than 0.')
  }else if(nsnps_cutoff>=max(in_tibble$nsnps,na.rm=T)){
    stop('The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }

  # test that the number of pairs to plot is sensible.
  if(is.null(N)){
    N <- min(nrow(in_tibble),50)
    if(verbose){
      cat(sprintf('No upper limit on number of pairs to plot given.\nPlotting first %i pairs.\n',N))
    }
  }else if(N>sum(in_tibble$nsnps>nsnps_cutoff)){
    N <- sum(in_tibble$nsnps>nsnps_cutoff)
    if(verbose){
      cat(sprintf('Number of pairs to plot was greater than the number of available pairs.\nPlotting first %i pairs.\n',N))
    }
  }

  # Define class breaks
  class_breaks <- tibble::tibble(relationship=factor(c('Same_Twins','First_Degree','Second_Degree','Unrelated'),
                                                     levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated')),
                                 meanPMR=((1-0.5^c(1:3,Inf)))*in_tibble$ave_rel[1],
                                 col=ggcolorhue(4))
  # Produce plot
  gg2 <- in_tibble %>%
    dplyr::arrange(pmr) %>%
    dplyr::mutate(relationship=factor(relationship,levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated'))) %>%
    dplyr::mutate(ymin=pmr-2*sd,ymax=pmr+2*sd) %>%
    dplyr::filter((nsnps>=nsnps_cutoff)) %>%
    dplyr::slice(1:N) %>%
    ggplot2::ggplot(ggplot2::aes(x=forcats::fct_reorder(pair,pmr),y=pmr))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=meanPMR,col=relationship),linetype='dashed',data=class_breaks,linewidth=1)+
    ggplot2::theme_bw()+
    ggplot2::geom_errorbar(ggplot2::aes(forcats::fct_reorder(pair,pmr),ymin=ymin,ymax=ymax))+
    ggplot2::geom_point(ggplot2::aes(fill=relationship,shape=relationship),size=2, show.legend = TRUE)+
    ggplot2::ylab("PMR")+
    ggplot2::xlab(" ")+
    # scale_y_continuous(limits=c(0,0.25))+
    ggplot2::scale_fill_discrete(drop=FALSE,name=NULL)+
    ggplot2::scale_colour_discrete(drop=FALSE,name=NULL,guide='none')+
    ggplot2::scale_shape_manual(drop=FALSE,name=NULL,values=c(21:24))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1, size=fntsize),legend.position='bottom')+
    ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(size=5)))

  gg2 %>%
    return()
}
