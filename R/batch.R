# Batch photobleaching step counting for single-molecule fluorescence trajectories.
#
# Author: Bingzhi Wang (2023), ported to R (2026). Built on the step-finding code of
# Konstantinos Tsekouras and Sina Jazani (Copyright (C) 2015, GNU GPL v3 or later).
#
# Input : a whitespace-separated text file, one trajectory per COLUMN, one frame per row
#         (e.g. traces exported from Igor Pro).
# Output: one CSV per trajectory, summary.csv, priors.csv and two histograms.
#
# Trajectory indices and frame numbers are 0-based throughout, so priors.csv files are
# interchangeable with the Python version of this tool.

#' Read a file of trajectories
#'
#' @param path whitespace-separated text file with one trajectory per column and one
#'   frame per row, no header.
#' @return A numeric matrix, frames x trajectories.
#' @export
read_traces <- function(path) {
  m <- as.matrix(utils::read.table(path, header = FALSE, colClasses = "numeric"))
  dimnames(m) <- NULL
  m
}

# ---------------------------------------------------------------------------------
# 1. Priors: set interactively by clicking on each trace
# ---------------------------------------------------------------------------------

.pvar <- function(v) mean((v - mean(v))^2)          # population variance (numpy's np.var)

# Read one line from the user, also when running under Rscript.
.ask <- function(prompt) {
  if (interactive()) return(readline(prompt))
  cat(prompt)
  con <- file("stdin")
  on.exit(close(con))
  ans <- readLines(con, n = 1)
  if (length(ans)) ans else ""
}

# Make sure there is an on-screen graphics device that supports locator().
.ensure_interactive_device <- function() {
  if (grDevices::dev.cur() > 1 && grDevices::dev.interactive()) return(invisible())
  if (interactive() && grDevices::dev.interactive(orNone = TRUE)) {
    grDevices::dev.new()
    return(invisible())
  }
  opened <- tryCatch({
    sys <- Sys.info()[["sysname"]]
    if (sys == "Darwin") grDevices::quartz(width = 11, height = 4.5)
    else if (.Platform$OS.type == "windows") grDevices::windows(width = 11, height = 4.5)
    else grDevices::x11(width = 11, height = 4.5)
    TRUE
  }, error = function(e) FALSE)
  if (!opened)
    stop("No interactive graphics device is available for clicking the priors. ",
         "Run from an interactive R session, or supply a priors file.", call. = FALSE)
  invisible()
}

#' Set the priors of one trajectory by clicking on its trace
#'
#' Plots the trace and waits for 5 left-clicks:
#' 1. the last bleaching step (where the trace drops to background for good);
#' 2. and 3. the start and end of a background stretch;
#' 4. and 5. the start and end of a stretch where exactly one fluorophore is on.
#'
#' Right-click (or Esc / "Finish" in RStudio) skips the trajectory. After the fifth
#' click you are asked in the console to accept, redo or skip.
#'
#' @param trace numeric vector of intensities.
#' @param i trajectory index, shown in the title.
#' @return A [priors()] object, or `NULL` if the trajectory is skipped.
#' @export
click_priors <- function(trace, i = 0) {
  .ensure_interactive_device()
  n <- length(trace)
  frames <- 0:(n - 1)
  labels <- c("last step", "bg start", "bg end", "1F start", "1F end")
  cols <- c("#C0392B", "#2471A3", "#2471A3", "#1E8449", "#1E8449")

  repeat {
    graphics::par(mar = c(4, 4, 4.5, 1))
    graphics::plot(frames, trace, type = "l", lwd = 0.8, col = "grey30",
                   xlab = "frame", ylab = "intensity")
    graphics::title(sprintf(
      "Trajectory %d: left-click 5 points (right-click / Esc to skip)\n%s", i,
      "1) last bleaching step   2-3) start/end of background   4-5) start/end of ONE fluorophore"),
      cex.main = 0.85, font.main = 1)

    pts <- numeric(0)
    for (k in 1:5) {
      p <- graphics::locator(1)
      if (is.null(p)) return(NULL)                     # right-click / Esc: skip
      pts <- c(pts, p$x)
      graphics::abline(v = p$x, col = cols[k], lty = 2)
      graphics::mtext(labels[k], side = 3, at = p$x, line = 0.1, cex = 0.7, col = cols[k])
    }

    x <- pmin(pmax(round(pts), 0), n - 1)
    bg_idx  <- sort(x[2:3]); one_idx <- sort(x[4:5])
    if (diff(bg_idx) < 2 || diff(one_idx) < 2) {
      message("Each stretch must span at least 2 frames; please click again.")
      next
    }
    ans <- tolower(trimws(.ask(sprintf(
      "Trajectory %d: [Enter] accept, [r] redo, [s] skip: ", i))))
    if (ans == "r") next
    if (ans == "s") return(NULL)
    break
  }

  # same half-open slices as trace[a:b] in Python
  bg  <- trace[(bg_idx[1] + 1):bg_idx[2]]
  one <- trace[(one_idx[1] + 1):one_idx[2]]
  priors(
    t_last = x[1],
    mB = mean(bg),
    vB = .pvar(bg),
    mF = mean(one) - mean(bg),   # intensity added by one fluorophore
    vF = .pvar(one)              # as in the original code: total variance of the
  )                              # 1-fluorophore level (model: vB + vF)
}

#' Write and read priors files
#'
#' `priors.csv` has columns `trajectory,t_last,mB,vB,mF,vF`; a skipped trajectory has
#' `skip` in the `t_last` column. The format is the same as in the Python version.
#'
#' @param path file path.
#' @param priors named list of [priors()] objects (or `NULL` for skipped trajectories),
#'   named by 0-based trajectory index.
#' @return `read_priors()` returns such a named list.
#' @export
write_priors <- function(path, priors) {
  num <- function(v) format(v, digits = 15, scientific = FALSE, trim = TRUE)
  lines <- vapply(names(priors), function(i) {
    p <- priors[[i]]
    if (is.null(p)) paste0(i, ",skip,,,,")
    else paste(i, p$t_last, num(p$mB), num(p$vB), num(p$mF), num(p$vF), sep = ",")
  }, character(1))
  writeLines(c("trajectory,t_last,mB,vB,mF,vF", lines), path)
  invisible(path)
}

#' @rdname write_priors
#' @export
read_priors <- function(path) {
  d <- utils::read.csv(path, colClasses = "character", strip.white = TRUE)
  out <- vector("list", nrow(d))
  names(out) <- as.character(as.integer(d$trajectory))
  for (k in seq_len(nrow(d))) {
    if (d$t_last[k] == "skip") next
    out[[k]] <- priors(as.integer(as.numeric(d$t_last[k])), as.numeric(d$mB[k]),
                       as.numeric(d$vB[k]), as.numeric(d$mF[k]), as.numeric(d$vF[k]))
  }
  out
}

# ---------------------------------------------------------------------------------
# 2. Outputs
# ---------------------------------------------------------------------------------

.write_trajectory <- function(path, trace, res) {
  body <- sprintf("%d,%.4f,%d,%d,%.4f", seq_along(trace) - 1L, trace,
                  as.integer(res$active), as.integer(res$steps), res$fit)
  writeLines(c("frame,intensity,active_fluorophores,step,fit", body), path)
}

.save_histograms <- function(path, counts, mF) {
  grDevices::png(path, width = 1500, height = 600, res = 150)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
  graphics::hist(counts, breaks = seq(0.5, max(counts) + 0.5, by = 1),
                 col = "#4C72B0", border = "white",
                 xlab = "fluorophores per trajectory", ylab = "trajectories",
                 main = "Photobleaching step count")
  graphics::hist(mF, col = "#4C72B0", border = "white",
                 xlab = "single-fluorophore intensity (mF)", ylab = "trajectories",
                 main = "Single-fluorophore brightness")
}

# ---------------------------------------------------------------------------------
# 3. Batch run
# ---------------------------------------------------------------------------------

#' Count photobleaching steps in a whole file of trajectories
#'
#' Loads all trajectories, gets the priors (by clicking, or from an earlier
#' `priors.csv`), runs [count_steps()] on each and writes the results.
#'
#' Files written to `outdir`:
#' \describe{
#'   \item{traj_XXX.csv}{per frame: `intensity`, `active_fluorophores`, `step`, `fit`}
#'   \item{summary.csv}{per trajectory: `n_fluorophores`, `n_steps`, `mB`, `vB`, `mF`, `vF`}
#'   \item{priors.csv}{your clicks (only when clicking), reusable via `priors`}
#'   \item{histograms.png}{fluorophore counts and single-fluorophore brightness}
#' }
#'
#' @param traces path to a text file (one trajectory per column) or a numeric matrix
#'   (frames x trajectories).
#' @param outdir output directory (created if needed).
#' @param window window size in frames (default 100).
#' @param priors `NULL` to click the priors on each trace, a path to a `priors.csv`
#'   from an earlier run, or a list as returned by [read_priors()].
#' @param only optional vector of 0-based trajectory indices to analyse.
#' @param verbose print progress.
#' @return The summary as a data frame (invisibly).
#' @examples
#' \donttest{
#' out <- file.path(tempdir(), "results")
#' run_batch(system.file("extdata", "synthetic_traces.txt", package = "photobleach"),
#'           outdir = out,
#'           priors = system.file("extdata", "synthetic_priors.csv", package = "photobleach"))
#' }
#' @export
run_batch <- function(traces, outdir = "results", window = 100, priors = NULL,
                      only = NULL, verbose = TRUE) {
  data <- if (is.character(traces)) read_traces(traces) else as.matrix(traces)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  todo <- if (is.null(only)) seq_len(ncol(data)) - 1L else as.integer(only)

  if (is.null(priors)) {
    pri <- vector("list", length(todo))
    names(pri) <- as.character(todo)
    for (k in seq_along(todo)) {
      p <- click_priors(data[, todo[k] + 1], todo[k])
      if (!is.null(p)) pri[[k]] <- p
    }
    write_priors(file.path(outdir, "priors.csv"), pri)
  } else if (is.character(priors)) {
    pri <- read_priors(priors)
  } else {
    pri <- priors
  }

  rows <- list()
  for (i in todo) {
    p <- pri[[as.character(i)]]
    if (is.null(p)) {
      if (verbose) message(sprintf("trajectory %d: skipped", i))
      next
    }
    prog <- if (verbose) function(k, n) cat(sprintf("\rtrajectory %d: window %d/%d", i, k, n))
    res <- count_steps(data[, i + 1], p, window = window, progress = prog)
    if (verbose) cat(sprintf("  -> %d fluorophores, %d steps\n", res$n_fluorophores, res$n_steps))
    .write_trajectory(file.path(outdir, sprintf("traj_%03d.csv", i)), data[, i + 1], res)
    rows[[length(rows) + 1]] <- data.frame(
      trajectory = i, n_fluorophores = res$n_fluorophores, n_steps = res$n_steps,
      mB = p$mB, vB = p$vB, mF = p$mF, vF = p$vF)
  }

  if (!length(rows)) return(invisible(NULL))
  summary <- do.call(rbind, rows)
  utils::write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE, quote = FALSE)
  .save_histograms(file.path(outdir, "histograms.png"), summary$n_fluorophores, summary$mF)
  if (verbose) cat(sprintf("Done: %d trajectories analysed, results in %s/\n",
                           nrow(summary), outdir))
  invisible(summary)
}
