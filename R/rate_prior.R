# Crude estimate of the photobleaching rate used in the mSIC prior.
#
# Ported from the original `LeffFinderPUB.py` (classes `LbarFind`, `PriorSlicer`)
# Copyright (C) 2015 Konstantinos Tsekouras, GNU GPL v3 or later.
# Simplified and annotated by Bingzhi Wang (2023); ported to R (2026).
#
# In the reversed trace, fluorophores switch on as an inhomogeneous Poisson process:
# the more fluorophores are already on, the faster the next one appears
# (lambda_eff = k * lambda_bar after the k-th event). The mSIC only needs a rough,
# trace-wide value, so we return the sum of lambda_eff over all points.

#' Per-fluorophore bleaching rate lambda_bar
#'
#' @param r reversed (and lifted) trace.
#' @param stats `c(mB, vB, mF, vF)`.
#' @param t0 0-based index in `r` of the first "switch-on" (= last bleaching step in
#'   real time).
#' @param n_edge points at the end of `r` (= start of the real trace) used to guess the
#'   initial number of fluorophores.
#' @return lambda_bar (per frame).
#' @keywords internal
estimate_lambda_bar <- function(r, stats, t0, n_edge = 25) {
  mB <- stats[[1]]; mF <- stats[[3]]
  n_max <- ceiling((mean(utils::tail(r, n_edge)) - mB) / mF)   # rough fluorophores at t = 0
  lit_frames <- length(r) - t0                                  # frames with >= 1 fluorophore
  if (n_max <= 1) return(10 / (max(n_max, 1) * lit_frames))
  harmonic <- sum(1 / seq_len(n_max - 1))                       # 1 + 1/2 + ... + 1/(n_max-1)
  10 / (harmonic * lit_frames)
}

#' Sum of lambda_eff over the trace
#'
#' lambda_eff(t) = k(t) * lambda_bar, where k(t) is the expected number of events so
#' far: the k-th event is expected ~1/(k * lambda_bar) frames after the previous one,
#' starting at `t0`.
#'
#' @inheritParams estimate_lambda_bar
#' @param lambda_bar output of [estimate_lambda_bar()].
#' @return A single number.
#' @keywords internal
rate_sum <- function(r, lambda_bar, t0) {
  lit_frames <- length(r) - t0
  marks <- 0; t <- 0; k <- 1
  repeat {
    t <- t + ceiling(1 / (k * lambda_bar))
    if (k == 1 && t > lit_frames)
      stop("lambda_bar too small for this trace; check the last-step position", call. = FALSE)
    if (t >= lit_frames) break
    marks <- c(marks, t)
    k <- k + 1
  }
  marks <- t0 + marks
  k_of_t <- findInterval(seq_along(r) - 1, marks)   # number of marks <= t
  sum(k_of_t * lambda_bar)
}
