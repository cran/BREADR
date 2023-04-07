#' saveSLICES
#'
#' Plots all pairwise diagnostic plots (in a tibble as output by callRelatedness),
#' as produced by plotSLICE, to a folder.
#' Options include the width and height of the output files, and the units in
#' which these dimensions are measured.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param outFolder the folder into which all diagnostic plots will be saved
#' @param width the width of the output PDFs.
#' @param height the height of the output PDFs.
#' @param units the units for the height and width of the output PDFs.
#' @param verbose Controls the printing of progress to console.
#'
#' @return nothing
#' @export
#'
#' @examples
#' \donttest{
#' saveSLICES(relatedness_example[1:3, ], outFolder = tempdir())
#' }
saveSLICES <- function(in_tibble,outFolder=NULL,width=297,height=210,units='mm', verbose = TRUE){


  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((ncol(in_tibble)!=12)|(!all(names(in_tibble)==col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }

  # Test to see if outfolder is set and exists
  if(is.null(outFolder)){
    stop('No output folder given.')
  }else{
    if(!dir.exists(outFolder)){
      stop('Output folder does not exist.')
    }
  }

  # Make sure output folder string ends in forward slash
  if(stringr::str_sub(outFolder,-1)!='/'){
    outFolder <- paste0(outFolder,'/')
  }

  # Test if plotting dimensions are numeric, and if the units are valid.
  if(any(!sapply(c(width,height),is.numeric))){
    stop('Height and width must be numeric values')
  }
  if(!units%in%c(c("in","cm","mm","px"))){
    stop("Plot dimension units must be 'in', 'cm', 'mm' or 'px'")
  }

  # Start plotting!
  if(verbose){
    cat('Starting pairwise plot creation.\n')
  }


  pb = utils::txtProgressBar(min=1,max=nrow(in_tibble),initial=1,style=3)
  for(i in 1:nrow(in_tibble)){
    v.sm <- plotSLICE(in_tibble,i,showPlot=F)
    ggplot2::ggsave(v.sm,
                    filename=paste0(outFolder,gsub(' - ','_',in_tibble$pair[i]),'.pdf'),
                    width=width,height=height,units=units)
    utils::setTxtProgressBar(pb,i)
  }
  if(verbose){
    cat(paste0('\nCompleted\nAll plots in: ',outFolder,'\n'))
  }
  return(TRUE)
}
