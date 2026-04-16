utils::globalVariables(
  c("sex", "pop1", "pop2", "group")
)
#' plotDOUGH
#'
#' @param in_tibble A tibble containing the output from processEigenstrat().
#' @param indfile Path to an ind file with columns for the individuals IDs, genetic sex (not used) and population.
#' @param nsnps_cutoff Minimum number of overlapping sites for a pair to be included.
#' @param fntsize Font size for x axis tick labels.
#' @param verbose Controls the printing of progress to console.
#' @param removeBetween TRUE/FALSE for whether to include the “between population” PMR values as their own group (FALSE), or to exclude them (TRUE).
#'
#' @returns plots
#'
#' @export
plotDOUGH <- function(
  in_tibble,
  indfile,
  nsnps_cutoff = NULL,
  fntsize = 7,
  verbose = TRUE,
  removeBetween = FALSE
) {
  # Check in_tibble is real
  # Test that the in_tibble is of the correct form
  if (nrow(in_tibble) == 0) {
    stop('The input tibble is empty')
  }
  col_names_check <- c('pair', 'nsnps', 'mismatch', 'pmr')
  if (!all(col_names_check %in% names(in_tibble))) {
    stop(paste0(
      'Input tibble/data.frame must have 4 columns named: ',
      paste0(col_names_check, collapse = ', '),
      '.'
    ))
  }

  # Test that the overlapping snp cut off is sensible.
  if (is.null(nsnps_cutoff)) {
    nsnps_cutoff <- 5e2
    if (verbose) {
      cat(sprintf(
        'No minimum number of overlapping SNPs given.\nUsing default minimum of 500.\n'
      ))
    }
  }
  if ((nsnps_cutoff <= 0)) {
    stop(
      'The cut off for the number of overlapping snps must be greater than 0.'
    )
  } else if (nsnps_cutoff >= max(in_tibble$nsnps, na.rm = T)) {
    stop(
      'The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!'
    )
  }
  # Check indfile exists
  if (!file.exists(indfile)) {
    stop(paste0('indfile ', indfile, ' does not exist!'))
  }

  IND <- readr::read_tsv(
    indfile,
    show_col_types = FALSE,
    col_names = c('ind', 'sex', 'pop')
  ) %>%
    dplyr::select(-sex)

  plot_data <- in_tibble %>%
    dplyr::filter(nsnps > nsnps_cutoff) %>%
    tidyr::separate(
      col = pair,
      into = c('ind1', 'ind2'),
      sep = ' - ',
      remove = FALSE
    ) %>%
    dplyr::left_join(IND, by = c('ind1' = 'ind')) %>%
    dplyr::left_join(IND, by = c('ind2' = 'ind'), suffix = c('1', '2')) %>%
    dplyr::mutate(group = ifelse(pop1 == pop2, pop1, 'Between Groups'))

  # Check that all individuals are accounted for
  fileInds <- c(plot_data$ind1, plot_data$ind2) %>%
    base::unique() %>%
    base::sort()

  missingInds <- base::setdiff(x = fileInds, y = IND$ind)

  if (length(missingInds) > 0) {
    missing_str <- missingInds %>%
      paste0(collapse = ', ') %>%
      sprintf(
        'The following individuals are in the BREADR analysis, but not the ind file: %s',
        .
      )
    stop(missing_str)
  }

  LEVELS <- plot_data$group %>%
    base::unique() %>%
    base::sort() %>%
    base::setdiff('Between Groups') %>%
    c(., 'Between Groups')

  outplot <- ggstatsplot::ggbetweenstats(
    data = plot_data %>%
      dplyr::mutate(group = factor(group, levels = LEVELS)) %>%
      dplyr::filter(
        group != ifelse(removeBetween, 'Between Groups', 'Can be Ignored Here')
      ),
    x = group,
    y = pmr,
    point.args = list(
      position = ggplot2::position_jitterdodge(dodge.width = 0.6),
      alpha = 0.5,
      size = 2,
      shape = 16
    ),
    type = 'nonparametric',
    plot.type = "box",
    drop = FALSE,
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    centrality.plotting = TRUE,
    bf.message = FALSE,
    ggtheme = ggplot2::theme_bw(),
    xlab = 'Population',
    ylab = 'PMR'
  )
  return(outplot)
}
