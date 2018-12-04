#' Helper function for Westfall Young Single Step
#'
#' This helper function compares permutated test statistics values under H0 with sample test statistics
#' under H1
#'
#' @param abs.Zs.H0.1row A vector of permutated test statistics values under H0
#' @param abs.Zs.H1.1samp One sample of H1 values
#' @param oo Order matrix of test statistics in descending order
#'
#' @return returns a vector of 1s and 0s with length of M outcomes
#' @export
comp.rawt.SS <- function(abs.Zs.H0.1row, abs.Zs.H1.1samp, oo) {

  M<-length(abs.Zs.H0.1row)
  maxt <- rep(NA, M)
  for (m in 1:M) {maxt[m] <- max(abs.Zs.H0.1row[m]) > abs.Zs.H1.1samp[m]}
  return(as.integer(maxt))
}
