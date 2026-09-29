ext <- function(f) system.file("extdata", f, package = "photobleach")

test_that("counts on the synthetic data match the Python reference implementation", {
  traces <- read_traces(ext("synthetic_traces.txt"))
  pri <- read_priors(ext("synthetic_priors.csv"))
  # n_fluorophores / n_steps produced by the Python version (window = 100)
  ref_n     <- c(3, 4, 1, 2, 4, 5, 2, 2, 5, 1, 4, 6, 1, 5, 3, 6, 1, 4, 4, 5)
  ref_steps <- c(7, 3, 1, 2, 3, 5, 3, 8, 4, 3, 3, 5, 1, 4, 3, 4, 1, 4, 2, 4)
  for (i in c(0, 5, 11, 15)) {
    res <- count_steps(traces[, i + 1], pri[[as.character(i)]])
    expect_equal(res$n_fluorophores, ref_n[i + 1])
    expect_equal(res$n_steps, ref_steps[i + 1])
  }
})

test_that("counts are close to the ground truth", {
  traces <- read_traces(ext("synthetic_traces.txt"))
  pri <- read_priors(ext("synthetic_priors.csv"))
  truth <- utils::read.csv(ext("synthetic_truth.csv"))
  found <- vapply(0:3, function(i) count_steps(traces[, i + 1], pri[[as.character(i)]])$n_fluorophores,
                  numeric(1))
  expect_true(all(abs(found - truth$n_fluorophores[1:4]) <= 1))
})

test_that("priors round-trip through priors.csv, including skipped trajectories", {
  pri <- list("0" = priors(10, 1, 2, 3, 4), "1" = NULL, "2" = priors(20, 1.5, 2.5, 3.5, 4.5))
  pri["1"] <- list(NULL)
  f <- tempfile(fileext = ".csv")
  write_priors(f, pri)
  back <- read_priors(f)
  expect_equal(names(back), c("0", "1", "2"))
  expect_null(back[["1"]])
  expect_equal(back[["2"]]$mF, 3.5)
  expect_equal(back[["0"]]$t_last, 10L)
})

test_that("step_prior matches a hand-computed value", {
  # K = 1, m = 1: 2*(log N - log R) + 2*gamma + log 2 - 2*log(2 - exp(-gamma))
  N <- 400; R <- 3; g <- 0.5
  expect_equal(photobleach:::step_prior(1, N, R, g),
               2 * (log(N) - log(R)) + 2 * g + log(2) - 2 * log(2 - exp(-g)))
})

test_that("run_batch writes all outputs", {
  out <- file.path(tempdir(), "pb_test")
  s <- run_batch(ext("synthetic_traces.txt"), outdir = out, only = c(0, 2),
                 priors = ext("synthetic_priors.csv"), verbose = FALSE)
  expect_equal(nrow(s), 2)
  expect_true(all(file.exists(file.path(out, c("traj_000.csv", "traj_002.csv",
                                               "summary.csv", "histograms.png")))))
})
