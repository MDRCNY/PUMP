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

#' Helper Functions for WestFallYoung Step down
#'
#' @param abs.Zs.H0.1row blah blah
#' @param abs.Zs.H1.1samp blah blah
#' @param oo blah blah
#'
#' @return blah blah
comp.rawt.SD <- function(abs.Zs.H0.1row, abs.Zs.H1.1samp, oo) {

  M<-length(abs.Zs.H0.1row)
  maxt <- rep(NA, M)
  nullt.oo<-abs.Zs.H0.1row[oo]
  rawt.oo<-abs.Zs.H1.1samp[oo]
  maxt[1] <- max(nullt.oo) > rawt.oo[1]
  for (h in 2:M) {maxt[h] <- max(nullt.oo[-(1:(h-1))]) > rawt.oo[h]}
  return(as.integer(maxt))
}

#' Adjust Single Step WestFallYoung
#'
#' @param snum blah blah
#' @param abs.Zs.H0 blah blah
#' @param abs.Zs.H1 blah blah
#'
#' @return blah blah
adjust.allsamps.WYSS<-function(snum,abs.Zs.H0,abs.Zs.H1) {

  adjp.WY<-matrix(NA,snum,ncol(abs.Zs.H0))
  doWY<-for (s in 1:snum) {
    ind.B<-t(apply(abs.Zs.H0, 1, comp.rawt.SS, abs.Zs.H1.1samp=abs.Zs.H1[s,]))
    adjp.WY[s,]<-colMeans(ind.B)
  }
  return(adjp.WY)
}

#' Adjust allsamps WYSD
#'
#' @param snum blah blah
#' @param abs.Zs.H0 blah blah
#' @param abs.Zs.H1 blah blah
#' @param order.matrix blah blah
#' @param ncl blah blah
#'
#' @return blah blah
#'
adjust.allsamps.WYSD<-function(snum,abs.Zs.H0,abs.Zs.H1,order.matrix,ncl) {

  cl <- snow::makeCluster(ncl)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl=cl, list("comp.rawt.SD"))
  M<-ncol(abs.Zs.H0)
  adjp.WY<-matrix(NA,snum,M)
   s <- 1:snum
   `%dopar%` <- `%dopar%`
  doWY <- foreach::foreach(s, .combine=rbind) %dopar% {
    ind.B<-t(apply(abs.Zs.H0, 1, comp.rawt.SD, abs.Zs.H1.1samp=abs.Zs.H1[s,], oo=order.matrix[s,]))
    pi.p.m <- colMeans(ind.B)
    # enforcing monotonicity
    adjp.minp <- numeric(M)
    adjp.minp[1] <- pi.p.m[1]
    for (h in 2:M) {adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])}
    adjp.WY[s,] <- adjp.minp[order.matrix[s,]]
  }
  return(doWY)
  parallel::stopCluster(cl)
}


#' t.mean.h1 helper function
#'
#' @param MDES blah blah
#' @param J blah blah
#' @param n.j blah blah
#' @param R2.1 blah blah
#' @param p blah blah
#'
#' @return blah blah
t.mean.H1<-function(MDES,J,n.j,R2.1,p) {

  MDES * sqrt(p*(1-p)*J*n.j) / sqrt(1-R2.1)
}

#' Degrees of Freedom
#'
#' @param J blah blah
#' @param n.j blah blah
#' @param numCovar.1 blah
#'
#' @return blah blah

df<-function(J,n.j,numCovar.1) {

  J*n.j - J - numCovar.1 - 1

}






















