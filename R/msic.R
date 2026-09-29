# Step finding with the modified Schwarz Information Criterion (mSIC).
#
# This is "Step 3" of the method of Tsekouras et al., Mol. Biol. Cell 27:3601 (2016),
# ported from the original `SeekerPUB.py` (class `mSICer`)
# Copyright (C) 2015 Konstantinos Tsekouras, GNU GPL v3 or later.
# Simplified, vectorised and annotated by Bingzhi Wang (2023); ported to R (2026).
#
# The trace is analysed *reversed in time*, so fluorophores appear to switch ON one
# after another starting from background. In one window of the trace the model is:
#
#     n active fluorophores  ->  intensity ~ Normal(mean = mB + n*mF, var = vB + n*vF)
#
# A candidate "step arrangement" is a vector as long as the window, zero everywhere
# except at step positions, where it holds the change in n (+1..+3, or -1..-3 for
# blinking). Its score is
#
#     mSIC = -2 log(likelihood)  +  prior penalty on the number/size of steps
#
# and the arrangement with the lowest score wins. The search is:
#   1. exhaustive over every arrangement with one or two steps in the window;
#   2. greedy: repeatedly add one more step wherever it lowers the score the most
#      (at most `max_extra_steps` times), stopping as soon as no addition helps.
#
# Positions inside a window are 0-based in the comments and helper arguments below,
# to stay line-by-line comparable with the Python version.

#' Allowed changes in fluorophore number per step
#' @keywords internal
JUMPS <- c(1L, 2L, 3L, -3L, -2L, -1L)

.prior_cache <- new.env(parent = emptyenv())

#' Prior part of the mSIC (Eq. 11 in Tsekouras et al. 2016)
#'
#' Penalises (i) many steps and (ii) many events squeezed into few steps
#' (overlapping bleaching), weighted by combinatorial degeneracy factors.
#'
#' @param abs_jumps integer vector of |jump| for every step, e.g. `c(1, 1, 2)`.
#' @param n_total number of points in the whole trace (N).
#' @param rate_sum crude overall photobleaching rate (sum of lambda_eff over the trace).
#' @param gamma cut-off gamma_0 of the prior that penalises overlapping events.
#' @return The prior penalty (a single number).
#' @keywords internal
step_prior <- function(abs_jumps, n_total, rate_sum,
                       gamma = getOption("photobleach.gamma0", 0.5)) {
  abs_jumps <- sort(as.integer(abs_jumps))
  key <- paste(paste(abs_jumps, collapse = ","), n_total,
               sprintf("%.17g", rate_sum), gamma, sep = "|")
  hit <- .prior_cache[[key]]
  if (!is.null(hit)) return(hit)

  K <- length(abs_jumps)
  m <- sum(abs_jumps)
  degeneracy <- sum(lfactorial(tabulate(abs_jumps)))
  p <- 2 * (degeneracy + lfactorial(m - 1) + K * log(n_total) - K * log(rate_sum) -
              lfactorial(m - K) - lfactorial(K))
  p <- p + 2 * gamma * (m - K + 1) / K + log(m - K + 2) + log(m - K + 1)
  p <- p - 2 * log(m - K + 2 - (m - K + 1) * exp(-gamma / K))

  assign(key, p, envir = .prior_cache)
  p
}

# Fast likelihood of any piecewise-constant fluorophore number on one window.
# Prefix sums give sum(x) and sum(x^2) of any segment in O(1).
.make_window <- function(x, stats) {
  x <- as.numeric(x)
  list(w = length(x), mB = stats[[1]], vB = stats[[2]], mF = stats[[3]], vF = stats[[4]],
       S1 = c(0, cumsum(x)), S2 = c(0, cumsum(x^2)))
}

# -2 log-likelihood (up to a constant) of x[start..end] (0-based, inclusive) at `level`
# fluorophores. Vectorised over start/end/level; an empty segment (start = end + 1) costs 0.
.seg <- function(win, start, end, level) {
  n   <- end - start + 1
  var <- level * win$vF + win$vB
  mu  <- level * win$mF + win$mB
  sx  <- win$S1[end + 2] - win$S1[start + 1]
  sxx <- win$S2[end + 2] - win$S2[start + 1]
  n * log(var) + (sxx - 2 * mu * sx + n * mu * mu) / var
}

# mSIC of a full step arrangement. Returns list(score, level_out); score = Inf if the
# arrangement ever needs a negative number of fluorophores.
.score <- function(win, jumps, level_in, n_total, rate_sum) {
  idx  <- which(jumps != 0)
  vals <- jumps[idx]
  if (any(level_in + cumsum(vals) < 0)) return(list(score = Inf, level = NA))
  pos <- idx - 1L                           # 0-based positions
  total <- 0; start <- 0; level <- level_in
  for (k in seq_along(pos)) {
    if (pos[k] == win$w - 1) break          # a step on the very last point has no data after
    total <- total + .seg(win, start, pos[k], level)   # it (kept from the original code)
    start <- pos[k] + 1
    level <- level + vals[k]
  }
  total <- total + .seg(win, start, win$w - 1, level)
  total <- total + step_prior(abs(vals), n_total, rate_sum)
  list(score = total, level = level)
}

#' Find the best step arrangement in one window
#'
#' @param x intensities of this window (reversed time).
#' @param level_in number of active fluorophores at the start of the window.
#' @param stats numeric vector `c(mB, vB, mF, vF)`: background / single-fluorophore
#'   mean and variance.
#' @param n_total number of points in the whole (padded) trace.
#' @param rate_sum output of [rate_sum()].
#' @return A list with `jumps` (as long as `x`, non-zero at steps) and `level`
#'   (number of active fluorophores at the end of the window).
#' @keywords internal
fit_window <- function(x, level_in, stats, n_total, rate_sum) {
  win <- .make_window(x, stats)
  w <- win$w
  L <- level_in
  max_extra <- getOption("photobleach.max_extra_steps", 9L)

  # --- baseline: no step at all -------------------------------------------------
  best <- .seg(win, 0, w - 1, L)
  best_jumps <- numeric(w)
  best_level <- L

  # --- 1. exhaustive search: exactly one step at position a ---------------------
  a <- 0:(w - 1)
  for (j in JUMPS) {
    if (L + j < 0) next
    cost <- .seg(win, 0, a, L) + .seg(win, a + 1, w - 1, L + j) +
      step_prior(abs(j), n_total, rate_sum)
    i <- which.min(cost)
    if (cost[i] < best) {
      best <- cost[i]
      best_jumps <- numeric(w); best_jumps[i] <- j
      best_level <- if (i < w) L + j else L
    }
  }

  # --- 1b. exhaustive search: two steps at positions a < b ----------------------
  if (w > 1) {
    A <- rep(0:(w - 2), times = (w - 1):1)                    # same order as numpy's
    B <- unlist(lapply(0:(w - 2), function(k) (k + 1):(w - 1)))  # triu_indices(w, k=1)
    head <- .seg(win, 0, A, L)
    for (j1 in JUMPS) for (j2 in JUMPS) {
      if (L + j1 < 0 || L + j1 + j2 < 0) next
      cost <- head + .seg(win, A + 1, B, L + j1) + .seg(win, B + 1, w - 1, L + j1 + j2) +
        step_prior(c(abs(j1), abs(j2)), n_total, rate_sum)
      i <- which.min(cost)
      if (cost[i] < best) {
        best <- cost[i]
        best_jumps <- numeric(w)
        best_jumps[A[i] + 1] <- j1
        best_jumps[B[i] + 1] <- j2
        best_level <- L + j1 + (if (B[i] < w - 1) j2 else 0)
      }
    }
  }

  # --- 2. greedy: keep adding the single most helpful extra step ----------------
  for (round in seq_len(max_extra)) {
    base <- best_jumps
    improved <- FALSE
    for (p in which(base == 0)) {
      for (j in JUMPS) {
        cand <- base; cand[p] <- j
        s <- .score(win, cand, L, n_total, rate_sum)
        if (s$score < best) {
          best <- s$score; best_jumps <- cand; best_level <- s$level; improved <- TRUE
        }
      }
    }
    if (!improved) break
  }

  list(jumps = best_jumps, level = best_level)
}
