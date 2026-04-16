utils::globalVariables(
  c(
    "prior_Same_Twins",
    "prior",
    "degree",
    "p",
    "prior_First_Degree",
    "prior_Second_Degree",
    "case",
    "prior_Unrelated"
  )
)
#' priorSensitivity
#'
#' @param in_BREADR A tibble containing the output of callRelatedness().
#' @param row Either the row number or pair name for which the sensitivity analysis should be run.
#' @param degrees A vector of integers identifying which degrees of relatedness degrees to plot. Same/Twins (1), First-Degree (2), Second-Degree (3) and Unrelated (4).
#' @param grid_space The space between prior probability values that are tested (smaller values mean a finer grid size).
#' @param maxPrior Maximum value a prior probability can take (the right-hand limit of the x-axis).
#'
#' @returns plot
#'
#' @export
priorSensitivity <- function(
  in_BREADR,
  row,
  degrees = NULL,
  grid_space = 0.05,
  maxPrior = 0.5
) {
  # Test that the in_tibble is of the correct form
  if (base::nrow(in_BREADR) == 0) {
    stop('The input tibble is empty')
  }
  col_names_check <- c(
    'row',
    'pair',
    'relationship',
    'pmr',
    'sd',
    'mismatch',
    'nsnps',
    'ave_rel',
    'Same_Twins',
    'First_Degree',
    'Second_Degree',
    'Unrelated'
  )
  if ((!base::all(col_names_check %in% names(in_BREADR)))) {
    stop(base::paste0(
      'Input tibble/data.frame must have columns named: ',
      paste0(col_names_check, collapse = ', '),
      '.'
    ))
  }
  if (
    !(base::all(
      degrees %in% c('Same_Twins', 'First_Degree', 'Second_Degree', 'Unrelated')
    ) |
      all(degrees %in% (1:4)))
  ) {
    stop(
      'degrees must be a vector of either:\n1. Numbers from c(1,2,3,4)\n2.Strings from c(Same_Twins,First_Degree,Second_Degree,Unrelated)'
    )
  }
  if (!base::is.numeric(grid_space) | (grid_space <= 0) | (grid_space > 0.5)) {
    stop('grid_space must be a number between 0 and 0.5')
  }
  if (
    !base::is.numeric(maxPrior) |
      (maxPrior <= 0) |
      (maxPrior < grid_space) |
      maxPrior <= 0.25
  ) {
    stop('maxPrior must be a number between 0.25 and grid_space')
  }
  # check row is a numeric or char variable, and if char, find the row.
  if (base::is.character(row)) {
    rowtmp <- row
    row <- base::which(row == in_BREADR$pair)
    if (base::length(row) == 0) {
      stop(base::paste0('There is no pair called: ', rowtmp))
    }
  } else if (base::is.numeric(row) & (row > base::nrow(in_BREADR))) {
    stop(
      'row number requested was higher than the number of rows in input tibble!'
    )
  } else if (!base::is.numeric(row)) {
    stop('row must be a number or a pair name.')
  }

  hues <- base::seq(15, 375, length = 5)
  PAL <- grDevices::hcl(h = hues, l = 65, c = 100)[1:4]
  CALLS <- c('Same_Twins', 'First_Degree', 'Second_Degree', 'Unrelated')
  if (base::is.null(degrees)) {
    degrees <- CALLS
  }
  if (base::is.numeric(degrees)) {
    degrees <- CALLS[degrees]
  }

  grid_in <- base::seq(
    from = grid_space,
    to = min(1 - 3 * grid_space, maxPrior),
    by = grid_space
  )
  grid_check <- base::list(
    Same_Twins = grid_in,
    First_Degree = grid_in,
    Second_Degree = grid_in,
    Unrelated = grid_in
  ) %>%
    base::expand.grid() %>%
    tidyr::as_tibble() %>%
    dplyr::filter(
      abs(Same_Twins + First_Degree + Second_Degree + Unrelated - 1) < 1e-6
    ) %>%
    dplyr::arrange(Same_Twins, First_Degree, Second_Degree, Unrelated) %>%
    dplyr::rename_with(function(x) {
      return(paste0('prior_', x))
    })

  PAIR <- in_BREADR$pair[row]

  base::sprintf(
    'Starting sensitivity analysis (%i separate analyses)....\n"',
    nrow(grid_check)
  ) %>%
    cat()
  pb <- utils::txtProgressBar(
    min = 1,
    max = nrow(grid_check),
    initial = 1,
    style = 3
  )
  for (i in 1:base::nrow(grid_check)) {
    if (i == 1) {
      res_list <- base::vector(mode = 'list', length = base::nrow(grid_check))
    }
    cur_prior <- grid_check %>%
      dplyr::slice(i) %>%
      base::unlist() %>%
      base::as.vector()
    res_list[[i]] <- callRelatedness(
      in_BREADR %>%
        dplyr::select(pair, nsnps, mismatch, pmr),
      average_relatedness = in_BREADR$ave_rel,
      class_prior = cur_prior
    ) %>%
      dplyr::select(Same_Twins, First_Degree, Second_Degree, Unrelated) %>%
      dplyr::slice(row)
    utils::setTxtProgressBar(pb, i)
  }
  cat("\nCollating results....\n")
  ST <- res_list %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(grid_check, .) %>%
    dplyr::select(
      prior_Same_Twins,
      Same_Twins,
      First_Degree,
      Second_Degree,
      Unrelated
    ) %>%
    dplyr::rename(prior = prior_Same_Twins) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      call = CALLS[which.max(c(
        Same_Twins,
        First_Degree,
        Second_Degree,
        Unrelated
      ))]
    ) %>%
    dplyr::group_by(prior) %>%
    dplyr::summarise(
      Same_Twins = sum(call == 'Same_Twins') / dplyr::n(),
      First_Degree = sum(call == 'First_Degree') / dplyr::n(),
      Second_Degree = sum(call == 'Second_Degree') / dplyr::n(),
      Unrelated = sum(call == 'Unrelated') / dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::gather(degree, p, 2:5) %>%
    dplyr::mutate(case = 'Same_Twins')
  FD <- res_list %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(grid_check, .) %>%
    dplyr::select(
      prior_First_Degree,
      Same_Twins,
      First_Degree,
      Second_Degree,
      Unrelated
    ) %>%
    dplyr::rename(prior = prior_First_Degree) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      call = CALLS[which.max(c(
        Same_Twins,
        First_Degree,
        Second_Degree,
        Unrelated
      ))]
    ) %>%
    dplyr::group_by(prior) %>%
    dplyr::summarise(
      Same_Twins = sum(call == 'Same_Twins') / dplyr::n(),
      First_Degree = sum(call == 'First_Degree') / dplyr::n(),
      Second_Degree = sum(call == 'Second_Degree') / dplyr::n(),
      Unrelated = sum(call == 'Unrelated') / dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::gather(degree, p, 2:5) %>%
    dplyr::mutate(case = 'First_Degree')
  SD <- res_list %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(grid_check, .) %>%
    dplyr::select(
      prior_Second_Degree,
      Same_Twins,
      First_Degree,
      Second_Degree,
      Unrelated
    ) %>%
    dplyr::rename(prior = prior_Second_Degree) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      call = CALLS[which.max(c(
        Same_Twins,
        First_Degree,
        Second_Degree,
        Unrelated
      ))]
    ) %>%
    dplyr::group_by(prior) %>%
    dplyr::summarise(
      Same_Twins = sum(call == 'Same_Twins') / dplyr::n(),
      First_Degree = sum(call == 'First_Degree') / dplyr::n(),
      Second_Degree = sum(call == 'Second_Degree') / dplyr::n(),
      Unrelated = sum(call == 'Unrelated') / dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::gather(degree, p, 2:5) %>%
    dplyr::mutate(case = 'Second_Degree')
  UR <- res_list %>%
    dplyr::bind_rows() %>%
    dplyr::bind_cols(grid_check, .) %>%
    dplyr::select(
      prior_Unrelated,
      Same_Twins,
      First_Degree,
      Second_Degree,
      Unrelated
    ) %>%
    dplyr::rename(prior = prior_Unrelated) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      call = CALLS[which.max(c(
        Same_Twins,
        First_Degree,
        Second_Degree,
        Unrelated
      ))]
    ) %>%
    dplyr::group_by(prior) %>%
    dplyr::summarise(
      Same_Twins = sum(call == 'Same_Twins') / dplyr::n(),
      First_Degree = sum(call == 'First_Degree') / dplyr::n(),
      Second_Degree = sum(call == 'Second_Degree') / dplyr::n(),
      Unrelated = sum(call == 'Unrelated') / dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::gather(degree, p, 2:5) %>%
    dplyr::mutate(case = 'Unrelated')

  out.gg.dat <- list(ST, FD, SD, UR) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(case %in% degrees) %>%
    dplyr::mutate(
      case = case %>%
        factor(
          levels = c(
            'Same_Twins',
            'First_Degree',
            'Second_Degree',
            'Unrelated'
          ),
          labels = c('Same/Twins', '1st-Degree', '2nd-Degree', 'Unrelated')
        ),
      degree = degree %>%
        factor(
          levels = c(
            'Same_Twins',
            'First_Degree',
            'Second_Degree',
            'Unrelated'
          ),
          labels = c('Same/Twins', '1st-Degree', '2nd-Degree', 'Unrelated')
        )
    )
  out.gg <- out.gg.dat %>%
    ggplot2::ggplot(ggplot2::aes(
      x = prior,
      y = p * 100,
      alpha = p,
      fill = degree
    )) +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = 'identity', col = 'black') +
    ggplot2::geom_hline(yintercept = 50, linetype = 'dashed') +
    ggplot2::facet_wrap(ggplot2::vars(case), nrow = base::length(degrees)) +
    ggplot2::xlab('Prior Probability') +
    ggplot2::ylab('Highest Posterior Calls (%)') +
    ggplot2::scale_y_continuous(breaks = base::seq(0, 100, by = 20)) +
    ggplot2::scale_fill_discrete(name = 'Degree') +
    ggplot2::scale_alpha_continuous(guide = 'none') +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ggtitle(PAIR) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  base::cat('Done.\n')
  return(out.gg)
}
