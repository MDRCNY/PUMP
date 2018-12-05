if(!exists("LinInterpolate", mode = "function")) source("R/utils.R")
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
#'
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
#'
#'
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
#'
#'
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
#'
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

#' Block RCT2 power function
#'
#' @param M the number of hypothesis tests (outcomes)
#' @param MDES vector of MDES's of length M - can be zero if not assuming all nulls are false
#' @param J number of blocks
#' @param n.j units per block
#' @param p the proportion of samples that are assigned the treatment
#' @param alpha significance level
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#' @param numCovar.2 number of Level 2 baseline covariates
#' @param R2.1 a vector of length M corresponding to R^2 for M outcomes of Level 1(R^2 = variation in the data explained by the model)
#' @param R2.2 a vector of length M corresponding to R^2 for M outcomes of Level 2
#' @param ICC intraclass correlation
#' @param mod.type "c" for constant effects, "f" for fixed effects, "r" for random effects
#' @param sigma correlation matrix for correlations between test statistics (this need to be flexible across multiple test statistics. Now, it's set to be at 1)
#' @param omega NULL
#' @param tnum number of test statistics (samples) for all procedures other than WY & number of permutations for WY (i.e permutation samples) WY has permutation which is reassignment of treatment variables. We do not see the treatment reassignment for permutation explicitly here.
#' @param snum number of samples for WY (i.e are they less than the other test. If so, why?)
#' @param ncl the number of clusters to use for parallel processing. It has a default of 2.
#'
#' @return power results
#' @export
#'
power.blockedRCT.2<-function(M, MDES, J, n.j,
                             p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC,
                             mod.type, sigma, omega,
                             tnum = 10000, snum=1000, ncl=2) {

  # number of false nulls
  numfalse<-sum(1*MDES>0)

  # compute Q(m) for all false nulls
  t.shift<-t.mean.H1(MDES,J,n.j,R2.1,p)
  t.df<-df(J,n.j,numCovar.1)
  t.shift.mat<-t(matrix(rep(t.shift,tnum),M,tnum)) # repeating shift.beta on every row

  # generate test statistics and p-values under null and alternative $s=\frac{1}{2}$
  Zs.H0<- mvtnorm::rmvt(tnum, sigma = sigma, df = t.df, delta = rep(0,M),type = c("shifted", "Kshirsagar"))
  Zs.H1 <- Zs.H0 + t.shift.mat
  pvals.H0<- stats::pt(-abs(Zs.H0),df=t.df) * 2
  pvals.H1<- stats::pt(-abs(Zs.H1),df=t.df) * 2
  abs.Zs.H0 <- abs(Zs.H0)
  abs.Zs.H1 <- abs(Zs.H1)

  # adjust p-values for all but Westfall-Young
  # mt.rawp2adjp <- multtest::mt.rawp2adjp
  adjp<-apply(pvals.H1,1,multtest::mt.rawp2adjp,proc=c("Bonferroni","Holm","BH"),alpha=alpha)
  rawp<-do.call(rbind,lapply(adjp,grab.pval,proc="rawp"))
  adjp.BF<-do.call(rbind,lapply(adjp,grab.pval,proc="Bonferroni"))
  adjp.HO<-do.call(rbind,lapply(adjp,grab.pval,proc="Holm"))
  adjp.BH<-do.call(rbind,lapply(adjp,grab.pval,proc="BH"))

  # adjust p-values for Westfall-Young (single-step and step-down)
  order.matrix<-t(apply(abs.Zs.H1,1,order,decreasing=TRUE))
  adjp.SS<-adjust.allsamps.WYSS(snum,abs.Zs.H0,abs.Zs.H1)
  adjp.WY<-adjust.allsamps.WYSD(snum,abs.Zs.H0,abs.Zs.H1,order.matrix,ncl)
  # combine all adjusted p-values in list (each entry is matrix for given MTP)
  adjp.all<-list(rawp,adjp.BF,adjp.HO,adjp.BH,adjp.SS,adjp.WY)

  # for each MTP, get matrix of indicators of whether adjusted p-value is less than alpha
  reject<-function(x) {as.matrix(1*(x<alpha))}
  reject.all<-lapply(adjp.all,reject)
  # in each row for each MTP matrix, count number of p-values less than 0.05, in rows corresponding to false nulls
  lt.alpha<-function(x) {apply(as.matrix(x[,MDES>0]),1,sum)}
  lt.alpha.all<-lapply(reject.all,lt.alpha)
  # indiv power for WY, BH, and HO is mean of columns of dummies of whether adjusted pvalues were less than alpha
  power.ind.fun<-function(x) {apply(x,2,mean)}
  power.ind.all<-lapply(reject.all,power.ind.fun)
  power.ind.all.mat<-do.call(rbind,power.ind.all)
  # m-min powers for all procs (including complete power when m=M)
  power.min.fun <- function(x,M) {
    power.min<-numeric(M)
    cnt<-0
    for (m in 1:M) {
      power.min[m]<-mean(x>cnt)
      cnt<-cnt+1
    }
    return(power.min)
  }
  power.min<-lapply(lt.alpha.all,power.min.fun,M=M)
  power.min.mat<-do.call(rbind,power.min)

  # complete power is probability all false nulls rejected when p-values not adjusted
  # this is row 1 and column number = numfalse
  power.cmp<-rep(power.min.mat[1,M],length(power.min)) # should it be numfalse or M?

  # combine all power for all definitions
  all.power.results<-cbind(power.ind.all.mat,power.min.mat[,-M],power.cmp)
  # take mean of all individual power estimates
  mean.ind.power <- apply(as.matrix(all.power.results[,1:M][,MDES>0]),1,mean)
  # revise final matrix to report this mean individual power and return results
  all.power.results<-cbind(mean.ind.power,all.power.results)
  colnames(all.power.results)<-c("indiv",paste0("indiv",1:M),paste0("min",1:(M-1)),"complete")
  rownames(all.power.results)<-c("rawp","BF","HO","BH","WY-SS","WY-SD")
  return(all.power.results)

}






















