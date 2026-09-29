#' photobleach: count photobleaching steps in single-molecule trajectories
#'
#' Batch implementation of the mSIC step-finding method of Tsekouras et al. (2016).
#' Start with [run_batch()] for a whole file of trajectories, or [count_steps()] for
#' a single trace.
#'
#' @section Options:
#' \describe{
#'   \item{`photobleach.gamma0`}{cut-off gamma_0 of the prior that penalises
#'     overlapping bleaching events (default `0.5`).}
#'   \item{`photobleach.max_extra_steps`}{number of greedy rounds after the exhaustive
#'     one- and two-step search in each window (default `9`).}
#' }
#' Set them with e.g. `options(photobleach.gamma0 = 0.3)`.
#'
#' @references Tsekouras K, Custer TC, Jashnsaz H, Walter NG, Pressé S (2016).
#'   A novel method to accurately locate and count large numbers of steps by
#'   photobleaching. Mol Biol Cell 27(22):3601-3615. \doi{10.1091/mbc.E16-06-0404}
#' @keywords internal
"_PACKAGE"
