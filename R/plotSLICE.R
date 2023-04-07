#' plotSLICE
#'
#' A function for plotting the diagnostic information when classifying a
#' specific pair (defined by the row number or pair name) of individuals.
#' Output includes the PDFs for each degree of relatedness (given the number of
#' overlapping SNPs) in panel A, and the normalised posterior probabilities
#' for each possible degree of relatedness.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param row either the row number or pair name for which the
#' posterior distribution is to be plotted.
#' @param title an optional title for the plot. If NULL, the pair
#' from the user-defined row is used.
#' @param class_prior the prior probabilities for same/twin, 1st-degree,
#' 2nd-degree, unrelated, respectively.
#' @param showPlot If TRUE, display plot. If FALSE, just pass plot as a variable.
#' @param which_plot if 1, returns just the plot of the posterior distributions,
#' if 2 returns just the normalised posterior values. Anything else returns both
#' plots.
#' @param labels a length two character vector of labels for plots. Default is
#' no labels.
#'
#' @return a two-panel diagnostic ggplot object
#' @export
#'
#' @examples
#' plotSLICE(relatedness_example, row = 1)
plotSLICE <- function(
    in_tibble,
    row,title=NULL,
    class_prior=rep(1/4,4),
    showPlot=TRUE,
    which_plot = 0,
    labels = NULL){

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

  # Check that the prior distribution makes sense:
  if(any(class_prior<0)|any(class_prior>1)){
    stop('Posterior probabilities for degrees of relatedness must be values between 0 and 1.')
  }else if(sum(class_prior)!=1){
    stop('Posterior probabilities for degrees of relatedness must sum to 1.')
  }

  # Check that showPlot is a logical variable
  if(!is.logical(showPlot)){
    stop('showPlot must be a logical (TRUE/FALSE) variable.')
  }

  # Check that labels is length two is present
  if(!is.null(labels)){
    stopifnot(length(labels) == 2)
  }

  # grab PMR info
  x <- in_tibble$mismatch[row]
  N <- in_tibble$nsnps[row]
  M <- in_tibble$ave_rel[row]

  # Set title if NULL
  if(is.null(title)){
    title <- in_tibble$pair[row]
  }

  # Set max number of points which are plotted at 5000 (uniformly). Set for efficiency.
  if(N>5e3){
    support <- seq(0,N,length.out=5e3)%>%round()
  }else{
    support <-0:N
  }

  # Find limit of x-axis.
  XLIM <- (M)+5*sqrt((M)*(1-M)/N)

  # Define possible classes.
  kclasses <- c('Same_Twins','First_Degree','Second_Degree','Unrelated')

  # find pdf values across "support" for each relatedness category
  distribution.tib <- tibble::tibble(x=rep(support,4),
                                     p.sim=rep((1-0.5^c(1:3,Inf))*M,each=length(support)),
                                     model=rep(kclasses,each=length(support))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y=(model=='Same_Twins'|model=='Second_Degree')*stats::dbinom(x,N,p.sim)+
                    (model=='First_Degree')*stats::dbinom(x,N,p.sim)+#dbbinom(x,N,alpha,beta)+
                    (model=='Unrelated')*weightedBinom(x,N,p.sim),
                  model=factor(model,levels=kclasses))

  # Define plot height for observed PMR below plot
  h <- max(distribution.tib$y)*0.05

  # Find highest probability points for dashed lines
  midpoint.tib <- distribution.tib %>%
    dplyr::group_by(model) %>%
    dplyr::filter(y==max(y))

  # find critical value
  zstar <- stats::qnorm(0.975)

  # make plot of probability curves
  distribution.plot <- distribution.tib %>%
    ggplot2::ggplot(ggplot2::aes(x=x/N,y=y))+
    ggplot2::scale_x_continuous(limits=c(0,XLIM),expand=c(0,0))+
    ggplot2::scale_y_continuous(limits=c(-h,max(distribution.tib$y)),expand=c(0,0))+
    ggplot2::annotate(geom='point',x=x/N,y=-h/2,pch=19,size=2)+
    ggplot2::geom_segment(x=x/N-zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N+zstar*sqrt((x/N)*(1-x/N)/N),y=-h/2,yend=-h/2)+
    ggplot2::geom_segment(x=x/N-zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N-zstar*sqrt((x/N)*(1-x/N)/N),y=-h/4,yend=-3*h/4)+
    ggplot2::geom_segment(x=x/N+zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N+zstar*sqrt((x/N)*(1-x/N)/N),y=-h/4,yend=-3*h/4)+
    ggplot2::theme_bw()+
    ggplot2::xlab(paste0('Pairwise Mismatch Rate (n=',N,')'))+
    ggplot2::ylab(NULL)+
    ggplot2::theme(legend.position='bottom')+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=0,ymax=y,fill=model),col='black',alpha=0.5)+
    ggplot2::geom_segment(data=midpoint.tib,ggplot2::aes(x=p.sim,xend=p.sim,y=0,yend=y*0.99),linetype='dashed',col='black')

  # return just distribution plot is asked for
  if(which_plot == 1){
    return(distribution.plot)
  }
  # Calculate and plot posterior probability values
  posterior.tibble <- tibble::tibble(posterior=c(stats::dbinom(x,N,(1-0.5^c(1:3))*M),weightedBinom(x,N,M))*class_prior/
                                       sum(c(stats::dbinom(x,N,(1-0.5^c(1:3))*M),weightedBinom(x,N,M))*class_prior),
                                     model=factor(kclasses,levels=kclasses))
  posterior.tibble <- posterior.tibble %>%
    dplyr::mutate(
      label =label_formatter(posterior)
      )
  # return(posterior.tibble)
  posterior.plot <- posterior.tibble %>%
    ggplot2::ggplot(ggplot2::aes(x=model,y=posterior))+
    ggplot2::theme_bw()+
    ggplot2::geom_bar(stat='identity',ggplot2::aes(fill=model),col='black',alpha=0.8)+
    ggplot2::geom_label(
      ggplot2::aes(
        label=label,
        colour=model
      ), parse = TRUE
    )+
    ggplot2::ylab("Posterior Probability")+
    ggplot2::xlab('Model')+
    ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(n.dodge=2))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  # Return just posterior plot is asked for
  if(which_plot == 2){
    return(posterior.plot)
  }
  # Combine plots
  # Old code left in as archaeology. Why would you not in a paper for relationship between samples.
  # gg1 <- gridExtra::arrangeGrob(distribution.plot+ggplot2::theme(legend.position='none'),
  #                               posterior.plot+ggplot2::theme(legend.position='none'),
  #                               nrow=1,
  #                               top = grid::textGrob(title,gp=grid::gpar(fontsize=20,font=3)))
  plot_grid <- ggpubr::ggarrange(
    distribution.plot,
    posterior.plot,
    nrow = 1,
    legend = "none",
    labels = labels
  )
  plot_grid <- ggpubr::annotate_figure(
    plot_grid,
    top = ggpubr::text_grob(title,size = 20)
  )
  # Plot if required
  if(showPlot){
    return(plot_grid)
  }
  invisible(plot_grid)
}
