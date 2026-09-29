# Count photobleaching steps in one fluorescence-intensity trajectory.
#
# Pipeline adapted by Bingzhi Wang (2023) from `PB_mainPUB.py`
# (Copyright (C) 2015 Konstantinos Tsekouras, GNU GPL v3 or later); ported to R (2026).

#' Priors for one trajectory
#'
#' What the user supplies for each trajectory, all in raw intensity units. Usually
#' created by [click_priors()] or read from a file with [read_priors()].
#'
#' @param t_last frame (0-based) of the last bleaching step (last fluorophore ->
#'   background).
#' @param mB,vB mean and variance of the background (after the last step).
#' @param mF mean intensity added by one fluorophore.
#' @param vF variance of the single-fluorophore level.
#' @return A list of class `photobleach_priors`.
#' @examples
#' p <- priors(t_last = 120, mB = 1000, vB = 62500, mF = 1500, vF = 185000)
#' @export
priors <- function(t_last, mB, vB, mF, vF) {
  structure(list(t_last = as.integer(t_last), mB = as.numeric(mB), vB = as.numeric(vB),
                 mF = as.numeric(mF), vF = as.numeric(vF)),
            class = "photobleach_priors")
}

#' @export
print.photobleach_priors <- function(x, ...) {
  cat(sprintf("<priors> t_last = %d, mB = %g, vB = %g, mF = %g, vF = %g\n",
              x$t_last, x$mB, x$vB, x$mF, x$vF))
  invisible(x)
}

# Pad the trace to a whole number of windows and lift it above zero.
.prepare <- function(trace, window) {
  x <- as.numeric(trace)
  rem <- length(x) %% window
  if (rem) {
    # repeat the last (window - rem) frames: these are background, so harmless
    x <- c(x, rep_len(utils::tail(x, window - rem), window - rem))
  }
  lift <- if (min(x) <= 0) abs(min(x)) + 1 else 0   # the model needs positive intensities
  list(x = x + lift, lift = lift)
}

#' Count photobleaching steps in one trajectory
#'
#' Runs the full algorithm: time reversal, photobleaching-rate prior and a windowed
#' mSIC step search. The fluorophore count at the end of one window is carried into
#' the next.
#'
#' @param trace numeric vector of intensities, one value per frame.
#' @param priors a [priors()] object for this trajectory.
#' @param window window size in frames (the most important tuning parameter).
#' @param progress optional function `function(k, n)` called after each window.
#' @return A list of class `photobleach_result` with
#'   \describe{
#'     \item{n_fluorophores}{maximum number of simultaneously active fluorophores}
#'     \item{n_steps}{number of frames where the fluorophore count changes}
#'     \item{active}{active fluorophores per frame (real time)}
#'     \item{steps}{per frame: >0 = that many fluorophores bleached, <0 = blinked on}
#'     \item{fit}{idealised intensity `mB + active * mF` (raw units)}
#'   }
#' @examples
#' traces <- read_traces(system.file("extdata", "synthetic_traces.txt",
#'                                   package = "photobleach"))
#' pri <- read_priors(system.file("extdata", "synthetic_priors.csv",
#'                                package = "photobleach"))
#' res <- count_steps(traces[, 1], pri[["0"]])
#' res$n_fluorophores
#' @export
count_steps <- function(trace, priors, window = 100, progress = NULL) {
  n <- length(trace)
  prep <- .prepare(trace, window)
  x <- prep$x
  stats <- c(priors$mB + prep$lift, priors$vB, priors$mF, priors$vF)

  # Work in reversed time: fluorophores now appear one by one on top of the background.
  r <- rev(x)
  t0 <- length(r) - priors$t_last        # position of the last bleaching step in `r`
  lam <- estimate_lambda_bar(r, stats, t0)
  lam_sum <- rate_sum(r, lam, t0)

  jumps  <- numeric(length(r))
  active <- numeric(length(r))
  level  <- 0                            # the reversed trace starts at background
  n_windows <- length(r) %/% window
  for (k in seq_len(n_windows)) {
    sl <- ((k - 1) * window + 1):(k * window)
    fw <- fit_window(r[sl], level, stats, length(r), lam_sum)
    jumps[sl]  <- fw$jumps
    active[sl] <- level + c(0, cumsum(fw$jumps)[-window])
    level <- fw$level                    # carry the fluorophore count into the next window
    if (!is.null(progress)) progress(k, n_windows)
  }

  # back to real time, drop the padding
  active <- rev(active)[seq_len(n)]
  steps  <- rev(jumps)[seq_len(n)]
  structure(list(
    n_fluorophores = as.integer(max(active)),
    n_steps = sum(diff(active) != 0),
    active = active,
    steps = steps,
    fit = priors$mB + active * priors$mF
  ), class = "photobleach_result")
}

#' @export
print.photobleach_result <- function(x, ...) {
  cat(sprintf("<photobleach result> %d fluorophores, %d steps over %d frames\n",
              x$n_fluorophores, x$n_steps, length(x$active)))
  invisible(x)
}
