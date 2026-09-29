#' Simulate photobleaching trajectories with a known number of fluorophores
#'
#' Useful to try out and benchmark the step counter. Each trajectory holds 1 to
#' `max_fluor` fluorophores with exponential lifetimes; the noise variance grows with
#' the number of active fluorophores (`sdB^2 + n * sdF^2`).
#'
#' R's random number generator differs from NumPy's, so these traces are not the same
#' as the ones in `inst/extdata` (which were made with the Python version and are
#' shipped for exact comparison).
#'
#' @param n_traj number of trajectories.
#' @param n_frames frames per trajectory.
#' @param mB,mF background and single-fluorophore mean intensity.
#' @param sdB,sdF background and per-fluorophore noise standard deviation.
#' @param mean_lifetime mean bleaching time in frames.
#' @param max_fluor maximum fluorophores per trajectory.
#' @param seed random seed.
#' @param outdir if given, write `synthetic_traces.txt`, `synthetic_truth.csv` and
#'   `synthetic_priors.csv` (the true priors, so no clicking is needed) there.
#' @return A list with `traces` (frames x trajectories matrix), `truth` (data frame)
#'   and `priors` (named list of [priors()]).
#' @examples
#' sim <- simulate_traces(n_traj = 3)
#' res <- count_steps(sim$traces[, 1], sim$priors[["0"]])
#' c(true = sim$truth$n_fluorophores[1], found = res$n_fluorophores)
#' @export
simulate_traces <- function(n_traj = 20, n_frames = 400, mB = 1000, mF = 1500,
                            sdB = 250, sdF = 350, mean_lifetime = 60, max_fluor = 6,
                            seed = 1, outdir = NULL) {
  set.seed(seed)
  traces <- matrix(0, n_frames, n_traj)
  truth <- data.frame(trajectory = seq_len(n_traj) - 1L, n_fluorophores = 0L, t_last = 0L)
  for (k in seq_len(n_traj)) {
    n <- sample.int(max_fluor, 1)
    t_bleach <- sort(as.integer(stats::rexp(n, 1 / mean_lifetime))) + 5L
    t_bleach <- pmin(t_bleach, n_frames - 60L)          # leave some background at the end
    active <- vapply(0:(n_frames - 1), function(t) sum(t_bleach > t), numeric(1))
    traces[, k] <- mB + active * mF + stats::rnorm(n_frames, 0, sqrt(sdB^2 + active * sdF^2))
    truth$n_fluorophores[k] <- n
    truth$t_last[k] <- max(t_bleach)
  }
  pri <- lapply(truth$t_last, function(t) priors(t, mB, sdB^2, mF, sdB^2 + sdF^2))
  names(pri) <- as.character(truth$trajectory)

  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    utils::write.table(format(round(traces, 2), nsmall = 2, trim = TRUE),
                       file.path(outdir, "synthetic_traces.txt"),
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.csv(truth, file.path(outdir, "synthetic_truth.csv"),
                     row.names = FALSE, quote = FALSE)
    write_priors(file.path(outdir, "synthetic_priors.csv"), pri)
  }
  list(traces = traces, truth = truth, priors = pri)
}
