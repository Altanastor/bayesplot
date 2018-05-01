#' Old version of \code{mcmc_areas}.
#'
#' @export
#' @template args-mcmc-x
#' @template args-pars
#' @template args-regex_pars
#' @template args-transformations
#' @param ... Currently unused.
#' @param prob The probability mass to include in the inner interval (for
#'   \code{mcmc_intervals}) or in the shaded region (for \code{mcmc_areas}). The
#'   default is \code{0.5} (50\% interval).
#' @param prob_outer The probability mass to include in the outer interval. The
#'   default is \code{0.9} for \code{mcmc_intervals} (90\% interval)
#'   and \code{1} for \code{mcmc_areas}.
#' @param point_est The point estimate to show. Either \code{"median"} (the
#'   default), \code{"mean"}, or \code{"none"}.
#' @param rhat An optional numeric vector of \eqn{\hat{R}}{Rhat} estimates, with
#'   one element per parameter included in \code{x}. If \code{rhat} is provided,
#'   the intervals/areas and point estimates in the resulting plot are colored
#'   based on \eqn{\hat{R}}{Rhat} value. See \code{\link{rhat}} for methods for
#'   extracting \eqn{\hat{R}}{Rhat} estimates.
#' @param bw,adjust,kernel For \code{mcmc_areas}, optional arguments passed to
#'   \code{\link[stats]{density}} to override default kernel density estimation
#'   parameters.
#'
#' @template return-ggplot
#'
mcmc_areas_old <- function(x,
                       pars = character(),
                       regex_pars = character(),
                       transformations = list(),
                       ...,
                       prob = 0.5,
                       prob_outer = 1,
                       point_est = c("median", "mean", "none"),
                       rhat = numeric(),
                       bw = NULL,
                       adjust = NULL,
                       kernel = NULL) {
  check_ignored_arguments(...)
  stopifnot(prob_outer >= prob)
  x <- prepare_mcmc_array(x, pars, regex_pars, transformations)
  .mcmc_intervals(
    x = merge_chains(x),
    prob_inner = prob,
    prob_outer = prob_outer,
    point_est = point_est,
    rhat = rhat,
    show_density = TRUE,
    bw = bw,
    adjust = adjust,
    kernel = kernel
  )
}


# internal ----------------------------------------------------------------

# @param x A matrix (not a 3-D array) created by merge_chains()
.mcmc_intervals <- function(x,
                            prob_inner = 0.5,
                            prob_outer = 0.95,
                            point_est = c("median", "mean", "none"),
                            rhat = numeric(),
                            show_density = FALSE,
                            bw = NULL,
                            adjust = NULL,
                            kernel = NULL) {
  n_param <- ncol(x)
  parnames <- colnames(x)

  probs <- c(0.5 - prob_outer / 2,
             0.5 - prob_inner / 2,
             0.5,
             0.5 + prob_inner / 2,
             0.5 + prob_outer / 2)

  quantiles <- t(apply(x, 2, quantile, probs = probs))
  y <- as.numeric(seq(n_param, 1, by = -1))
  x_lim <- c(min(quantiles[, 1]), max(quantiles[, 5]))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range

  data <- data.frame(parnames, y, quantiles)
  colnames(data) <- c("parameter", "y", "ll", "l", "m", "h", "hh")
  point_est <- match.arg(point_est)
  no_point_est <- identical(point_est, "none")
  if (point_est == "mean")
    data$m <- unname(colMeans(x))

  color_by_rhat <- isTRUE(length(rhat) > 0)
  if (color_by_rhat) {
    rhat <- factor(factor_rhat(rhat), levels = c("high", "ok", "low"))
    if (length(rhat) != nrow(data))
      stop(
        "'rhat' has length ", length(rhat),
        " but 'x' has ", nrow(data), " parameters.",
        call. = FALSE
      )

    data$rhat <- rhat
  }

  graph <- ggplot(data)

  # faint vertical line at zero if zero is within x_lim
  if (0 > x_lim[1] && 0 < x_lim[2])
    graph <- graph + vline_0(color = "gray90", size = 0.5)

  if (show_density) {
    # density outline
    n_dens_pts <- 512
    y_dens <- matrix(0, nrow = n_dens_pts, ncol = n_param)
    x_dens <- matrix(0, nrow = n_dens_pts, ncol = n_param)
    for (i in seq_len(n_param)) {
      d_temp <- compute_dens_i(
        x = x[, i],
        from = quantiles[i, 1],
        to = quantiles[i, 5],
        n = n_dens_pts,
        bw = bw,
        adjust = adjust,
        kernel = kernel
      )
      x_dens[, i] <- d_temp$x
      y_max <- max(d_temp$y)
      y_dens[, i] <- d_temp$y / y_max * 0.8 + y[i]
    }
    df_dens <- data.frame(
      x = as.vector(x_dens),
      y = as.vector(y_dens),
      name = rep(parnames, each = n_dens_pts)
    )
    if (color_by_rhat)
      df_dens$rhat <- rep(rhat, each = n_dens_pts)

    dens_args <- list(
      data = df_dens,
      mapping = aes_(
        x = ~ x,
        y = ~ y,
        group = ~ name,
        color = if (!color_by_rhat) NULL else ~ rhat
      )
    )
    if (!color_by_rhat)
      dens_args$color <- get_color("d")
    g_dens <- do.call("geom_line", dens_args)

    #shaded interval
    y_poly <- matrix(0, nrow = n_dens_pts + 2, ncol = n_param)
    x_poly <- matrix(0, nrow = n_dens_pts + 2, ncol = n_param)
    for (i in seq_len(n_param)) {
      d_temp <- compute_dens_i(
        x = x[, i],
        from = quantiles[i, 2],
        to = quantiles[i, 4],
        n = n_dens_pts,
        bw = bw,
        adjust = adjust,
        kernel = kernel
      )
      x_poly[, i] <-
        c(d_temp$x[1], as.vector(d_temp$x), d_temp$x[n_dens_pts])
      y_max <- max(d_temp$y)
      y_poly[, i] <-
        as.vector(c(0, as.vector(d_temp$y) / y_max * 0.8, 0) + y[i])
    }
    df_poly <-
      data.frame(
        x = as.vector(x_poly),
        y = as.vector(y_poly),
        name = rep(parnames, each = n_dens_pts + 2)
      )
    if (color_by_rhat)
      df_poly$rhat <- rep(rhat, each = n_dens_pts + 2)
    g_poly <-
      geom_polygon(data = df_poly, aes_(
        x = ~ x,
        y = ~ y,
        group = ~ name,
        fill = if (color_by_rhat) ~ rhat else ~ y
      ))

    # point estimate
    df_dens$parameter <- df_dens$name
    pt_data <- dplyr::summarise_(
      # find y value at which to stop vertical pt est segment
      dplyr::group_by_(
        dplyr::left_join(df_dens, data[, c("parameter", "m")],
                         by = "parameter"),
        .dots = list(~ parameter)
      ),
      .dots = list(maxy = ~ y[which.min(abs(x - m))])
    )
    segment_args <- list(
      data = dplyr::left_join(data, pt_data, by = "parameter"),
      mapping = aes_(
        x = ~ m,
        xend = ~ m,
        y = ~ y,
        yend = ~ maxy,
        color = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 1
    )
    if (!color_by_rhat)
      segment_args$color <- get_color("m")
    g_point <- do.call("geom_segment", segment_args)

    # bottom line
    bottom_args <- list(
      mapping = aes_(
        x = ~ ll,
        xend = ~ hh,
        y = ~ y,
        yend = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat
      )
    )
    if (!color_by_rhat)
      bottom_args$color <- get_color("d")
    g_bottom <- do.call("geom_segment", bottom_args)

    graph <- graph + g_poly
    if (!no_point_est)
      graph <- graph + g_point
    graph <- graph + g_bottom + g_dens

    if (color_by_rhat) {
      graph <- graph +
        scale_fill_diagnostic("rhat") +
        scale_color_diagnostic("rhat")
    } else {
      graph <- graph + scale_fill_gradient(low = get_color("l"),
                                           high = get_color("l"),
                                           guide = "none")
    }

  } else { # No densities

    # outer interval
    graph <-
      graph + geom_segment(aes_(
        x = ~ ll,
        xend = ~ hh,
        y = ~ y,
        yend = ~ y
      ),
      colour = get_color("m")
      )

    # inner interval
    segment_args <- list(
      mapping = aes_(
        x = ~ l,
        xend = ~ h,
        y = ~ y,
        yend = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 2,
      show.legend = FALSE
    )
    if (!color_by_rhat)
      segment_args$color <- get_color("d")
    graph <- graph + do.call("geom_segment", segment_args)

    # point estimate
    point_args <- list(
      mapping = aes_(
        x = ~ m,
        y = ~ y,
        color = if (!color_by_rhat) NULL else ~ rhat,
        fill = if (!color_by_rhat) NULL else ~ rhat
      ),
      size = 4,
      shape = 21
    )
    if (!color_by_rhat) {
      point_args$color <- get_color("dh")
      point_args$fill <- get_color("l")
    }

    if (!no_point_est)
      graph <- graph + do.call("geom_point", point_args)

    if (color_by_rhat)
      graph <- graph +
      scale_color_diagnostic("rhat") +
      scale_fill_diagnostic("rhat")
  }

  graph +
    scale_y_continuous(
      breaks = y,
      labels = parnames,
      limits = c(0.5, n_param + 1)
    ) +
    xlim(x_lim) +
    legend_move(ifelse(color_by_rhat, "top", "none")) +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    yaxis_ticks(size = 1) +
    xaxis_title(FALSE)
}


# compute kernel density estimates
# all arguments are passed to stats::density
compute_dens_i <- function(x, bw, adjust, kernel, n, from, to) {
  args <- c(
    # can't be null
    list(
      x = x,
      from = from,
      to = to,
      n = n
    ),
    # might be null
    bw = bw,
    adjust = adjust,
    kernel = kernel
  )
  do.call("density", args)
}
