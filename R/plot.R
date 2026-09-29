#' Plot a trajectory together with its fitted photobleaching steps
#'
#' @param x a `traj_XXX.csv` file written by [run_batch()], or a data frame with
#'   columns `frame`, `intensity`, `active_fluorophores` and `fit`.
#' @param file optional PNG file to write instead of drawing on the current device.
#' @param ... further arguments passed to [graphics::plot()].
#' @return The data (invisibly).
#' @examples
#' out <- file.path(tempdir(), "results")
#' run_batch(system.file("extdata", "synthetic_traces.txt", package = "photobleach"),
#'           outdir = out, only = 0, verbose = FALSE,
#'           priors = system.file("extdata", "synthetic_priors.csv", package = "photobleach"))
#' plot_fit(file.path(out, "traj_000.csv"))
#' @export
plot_fit <- function(x, file = NULL, ...) {
  d <- if (is.character(x)) utils::read.csv(x) else x
  if (!is.null(file)) {
    grDevices::png(file, width = 1350, height = 480, res = 150)
    on.exit(grDevices::dev.off())
  }
  op <- graphics::par(mar = c(4, 5, 1, 1), las = 1)
  on.exit(graphics::par(op), add = TRUE, after = FALSE)
  graphics::plot(d$frame, d$intensity, type = "l", lwd = 0.7, col = "grey60",
                 xlab = "frame", ylab = "intensity", ...)
  graphics::lines(d$frame, d$fit, type = "s", lwd = 2, col = "#D62728")
  graphics::legend("topright", bty = "n", lwd = c(0.7, 2), col = c("grey60", "#D62728"),
                   legend = c("data", sprintf("fit: %d fluorophores",
                                              as.integer(max(d$active_fluorophores)))))
  invisible(d)
}
