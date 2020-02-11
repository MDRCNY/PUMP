# the if statement checks if we have a grab.pval function. If not, pull it from utils.R file.
if(!exists("grab.pval", mode = "function")) source("R/utils.R")

#' Helper function for Westfall Young Single Step
#'
#' The  function  comp.rawt.SS is  needed  to  implement  the  Westfall-Young single-step multiple
#' testing procedure (MTP). It operates on one row of null test statistics.
#'
#' @param abs.Zs.H0.1row A vector of permutated test statistics values under H0
#' @param abs.Zs.H1.1samp One sample of raw statistics
#' @param oo Order matrix of test statistics in descending order (Only used in Step Down)
#'
#' @return returns a vector of 1s and 0s with length of M outcomes
#'
#'
comp.rawt.SS <- function(abs.Zs.H0.1row, abs.Zs.H1.1samp, oo) {

  # getting the number of M outcomes from 1 row of H0
  M <- length(abs.Zs.H0.1row)
  # creating an empty vector of length M to save boolean values
  maxt <- rep(NA, M)

  for (m in 1:M) {

    # comparing the maximum of null test values of M outcomes to each of the alternative test raw sample values
    # saving each M boolean in maxt vector
    maxt[m] <- max(abs.Zs.H0.1row) > abs.Zs.H1.1samp[m]

  }
  return(as.integer(maxt))
}

#' Helper Functions for WestFallYoung Step down
#'
#' @param abs.Zs.H0.1row A vector of permutated test statistics values under H0
#' @param abs.Zs.H1.1samp One sample of raw statistics
#' @param oo Order matrix of test statistics in descending order
#' @return returns a vector of 1s and 0s with lengths of M outcomes
#' @export
#'
comp.rawt.SD <- function(abs.Zs.H0.1row, abs.Zs.H1.1samp, oo) {

  # getting M number of outcomes from 1 row of statistics
  M <- length(abs.Zs.H0.1row)
  # creating an empty vector of length M to save boolean values
  maxt <- rep(NA, M)
  # saving the null test statistics
  nullt.oo <- abs.Zs.H0.1row[oo]
  # saving the raw test statistics under H1
  rawt.oo <- abs.Zs.H1.1samp[oo]
  # saving the first boolean by comparing the max of null values with the first of raw test statistics
  maxt[1] <- max(nullt.oo) > rawt.oo[1]

  # Step-down comparison where the next max of null values is compared to the next raw test statistics
  for (h in 2:M) {
    maxt[h] <- max(nullt.oo[-(1:(h-1))]) > rawt.oo[h]
  } # end of for loop

  return(as.integer(maxt))
}

#' WestFallYoung Single Step Adjustment Function
#'
#' This adjustment function utilizes the comp.rawt.SS helper function to compare
#' each row of the matrix sample test statistics under
#' alternative hypothesis to all the rows in the matrix of the test statistics under the complete null (i.e think a distribution).
#' Furthermore, it carries out the comparison for all the samples of raw test statistics under the alternative.
#'
#' @param snum the number of samples for which test statistics under the alternative hypothesis
#' are compared with the distribution (matrix) of test statistics under the complete null (this distribution
#' is obtained through drawing test values under H0 with a default of 10,000)
#' @param abs.Zs.H0 a matrix of test statistics under the complete null
#' @param abs.Zs.H1 a matrix of raw test statistics under the alternative
#'
#' @return a matrix of adjusted test statistics values

adjust.allsamps.WYSS<-function(snum,abs.Zs.H0,abs.Zs.H1) {

  # creating the matrix to store the adjusted test values with the number of samples &
  # number of M outcomes
  adjp.WY<-matrix(NA,snum,ncol(abs.Zs.H0))
  # looping through all the samples of raw test statistics under the alternative hypothesis
  doWY<-for (s in 1:snum) {

    # using apply to compare the distribution of test statistics under H0 with 1 sample of the raw statistics under H1
    ind.B<-t(apply(abs.Zs.H0, 1, comp.rawt.SS, abs.Zs.H1.1samp=abs.Zs.H1[s,]))
    # calculating the p-value for each sample
    adjp.WY[s,]<-colMeans(ind.B)

  }
  return(adjp.WY)
}

#' Westfall Young Step Down Function
#'
#' This adjustment function utilizes the comp.rawt.SD helper function to compare
#' each row of the matrix sample test statistics under
#' alternative hypothesis to all the rows in the matrix of the test statistics under the complete null (i.e think a distribution).
#' Furthermore, it carries out the comparison for all the samples of raw test statistics under the alternative.
#'
#' @param snum the number of samples for which test statistics under the alternative hypothesis
#' are compared with the distribution (matrix) of test statistics under the complete null (this distribution
#' is obtained through drawing test values under H0 with a default of 10,000)
#' @param abs.Zs.H0 a matrix of test statistics under the complete null
#' @param abs.Zs.H1 a matrix of raw test statistics under the alternative
#' @param order.matrix Order matrix of test statistics in descending order
#' @param ncl number of clusters to be made for parallelization. The default is 2.
#'
#' @return a matrix of adjusted test statistics values

adjust.allsamps.WYSD<-function(snum,abs.Zs.H0,abs.Zs.H1,order.matrix,ncl) {

  # creates clusters to run parallelization on
  cl <- snow::makeCluster(ncl)
  # leveraging snow to run multiple cores for foreach loops
  doParallel::registerDoParallel(cl)
  # registering the comp.rawt.SD function in global enivronment of each node
  parallel::clusterExport(cl=cl, list("comp.rawt.SD"))
  # getting M number of outcomes vector
  M <- ncol(abs.Zs.H0)
  # setting up the matrix to save the adjusted p values
  adjp.WY<-matrix(NA,snum,M)
  # dopar is a special function that has to be explicitly called from the foreach package
  # dopar accepts only 2 parameters. The number of times to execute the parallelization and the
  # series of steps to execute
  `%dopar%` <- foreach::`%dopar%`
  # making s a local variable to perpetuate across (created to bypass a package requirement)
  s = 1:snum
    doWY <- foreach::foreach(s= 1:snum, .combine=rbind) %dopar% {
    # using apply to compare the distribution of test statistics under H0 with 1 sample of the raw statistics under H1
    ind.B <- t(apply(abs.Zs.H0, 1, comp.rawt.SD, abs.Zs.H1.1samp=abs.Zs.H1[s,], oo=order.matrix[s,]))
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

#' t.mean.h1 function for generating the mean of test statistics under the joint alternative hypothesis
#'
#' The function t.mean.H1 computes the means of the test statistics under the joint alternative hypothesis.
#'
#' @param MDES  a vector of length M corresponding to the minimum detectable effect sizes (MDESs) for the M outcomes
#' @param J the number of blocks
#' @param n.j the harmonic means of the number of units per block
#' @param R2.1 a vector of length M corresponding to R^2 for Level-1 covariates for M outcomes
#' @param p the proportion of test statistics assigned to treatment within each block group
#'
#' @return mean of the test statistics under the joint alternative hypothesis

t.mean.H1<-function(MDES,J,n.j,R2.1,p) {

  MDES * sqrt(p*(1-p)*J*n.j) / sqrt(1-R2.1)
}

#' Degrees of Freedom
#'
#' This function calculates the degree of freedom for 2 level RCT block design
#'
#' @param J number of blocks
#' @param n.j units per block
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#'
#' @return the degree of freedom

df<-function(J,n.j,numCovar.1) {

  J*n.j - J - numCovar.1 - 1

}

#' Block RCT2 power function
#'
#' This functions calculates power for all definitions of power (individual, d-minimal, complete) for all the different MTPs
#' (Bonferroni, Holms, Bejamini-Hocheberg, Westfall-Young Single Step, Westfall-Young Step Down). The function works for
#' 2-Level Block RCT
#'
#' @param M the number of hypothesis tests (outcomes)
#' @param MDES a vector of length M corresponding to the MDESs for the M outcomes
#' @param Ai a single entry vector specifying the estimated number of outcomes with a non-zero effect
#' @param J the number of blocks
#' @param n.j the harmonic mean of the number of units per block
#' @param p the proportion of samples that are assigned to the treatment
#' @param alpha the family wise error rate (FWER)
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#' @param numCovar.2 number of Level 2 baseline covariates (set to 0 for this design)
#' @param R2.1 a vector of length M corresponding to R^2 for M outcomes of Level 1 (R^2 = variation in the data explained by the model)
#' @param R2.2 a vector of length M corresponding to R^2 for M outcomes of Level 2 (R^2 = variation in the data explained by the model)
#' @param ICC intraclass correlation
#' @param mod.type "c" for constant effects, "f" for fixed effects, "r" for random effects (parameter not in use at the moment)
#' @param sigma correlation matrix for correlations between test statistics (parameter not in use at the moment.Default is set to 0.99)
#' @param omega NULL (parameter in development)
#' @param tnum the number of test statistics (samples) for all procedures other than Westfall-Young & number of permutations for WY. The default is set at 10,000
#' @param snum the number of samples for Westfall-Young. The default is set at 1,000.
#' @param ncl the number of clusters to use for parallel processing. The default is set at 2.
#' @param rho correlation between outcomes
#' @param updateProgress the callback function to update the progress bar (User does not have to input anything)
#'
#' @importFrom multtest mt.rawp2adjp
#' @return power results across all definitions of power and MTP
#' @export

power.blockedRCT.2<-function(M, MDES, Ai, J, n.j,
                             p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2 = NULL, ICC,
                             mod.type, sigma = 0,rho = 0.99, omega = NULL,
                             tnum = 10000, snum=1000, ncl=2, updateProgress = NULL) {

  # Error handling when user put in actual effect number that is greater than the total number of outcomes
  if( Ai > M){

    stop('The number of outcomes with actual effects cannot be greater than the total number of outcomes of an experiment.
         Please adjust your inputs.')

  } # end of if statement

  # MDES must be the length of Actual Impacts
  MDES <- rep(MDES,Ai)

  # the difference between the length of M and Ai would be zero as they do not have any impacts
  noeffect <- rep(0, M-Ai)
  MDES <- c(MDES, noeffect)

  # Setting a default Sigma up
  sigma <- matrix(rho, M, M)
  diag(sigma) <- 1

  # number of false nulls (i.e they are really not nulls)
  numfalse <- sum(1*MDES>0)

  # compute Q(m) for all false nulls. We are calculating the test statistics for when the alternative hypothesis is true.
  t.shift<-t.mean.H1(MDES,J,n.j,R2.1,p)
  t.df<-df(J,n.j,numCovar.1)
  t.shift.mat<-t(matrix(rep(t.shift,tnum),M,tnum)) # repeating shift.beta on every row

  # generate test statistics and p-values under null and alternative $s=\frac{1}{2}$
  # rmvt draws from a multivariate t-distribution
  Zs.H0<- mvtnorm::rmvt(tnum, sigma = sigma, df = t.df, delta = rep(0,M),type = c("shifted", "Kshirsagar"))
  Zs.H1 <- Zs.H0 + t.shift.mat

  # calculates p values from quantiles and degrees of freedom
  pvals.H0<- stats::pt(-abs(Zs.H0),df=t.df) * 2
  pvals.H1<- stats::pt(-abs(Zs.H1),df=t.df) * 2

  # getting the absolute values of the test statistics
  abs.Zs.H0 <- abs(Zs.H0)
  abs.Zs.H1 <- abs(Zs.H1)

  # 1st call back to progress bar on progress of calculation: P values generation
  if (is.function(updateProgress) & !is.null(abs.Zs.H0)) {
    msg  <- paste0("P-values have been generated!") # Priamry text we want to display
    updateProgress(message = msg) # Passing back the progress messages onto the callback function
  } # if the function is being called, run the progress bar


  # explicit invocation of the function and saving it so that we can pass it along in apply without the
  # explicit calling inside apply
  mt.rawp2adjp <- multtest::mt.rawp2adjp

  #adjp<-multtest::mt.rawp2adjp(pvals.H1,proc=c("Bonferroni","Holm","BH"),alpha=alpha)
  #rawp <- grab.pval(adjp, proc = "rawp")
  #adjp.BF <- grab.pval(adjp, proc = "Bonferroni")
  #adjp.HO <- grab.pval(adjp, proc = "Holm")
  #adjp.BH <- grab.pval(adjp, proc = "BH")

  # adjust p-values for all but Westfall-Young using multtest's adjustment function
  # via Bonferroni, Holm and Bejamini-Hocheberg
  adjp<-apply(pvals.H1,1,mt.rawp2adjp,proc=c("Bonferroni","Holm","BH"),alpha=alpha)

  # seperating out p values that are adjusted by Bonferroni, Holm and Benjamini-Hocheberg
  grab.pval<-function(...,proc) {return(...$adjp[order(...$index),proc])}
  rawp<-do.call(rbind,lapply(adjp,grab.pval,proc="rawp"))
  adjp.BF<-do.call(rbind,lapply(adjp,grab.pval,proc="Bonferroni"))
  adjp.HO<-do.call(rbind,lapply(adjp,grab.pval,proc="Holm"))
  adjp.BH<-do.call(rbind,lapply(adjp,grab.pval,proc="BH"))

  # 2nd call back to progress bar on adjustments except WestFall Young is done
  if (is.function(updateProgress) & !is.null(adjp.BH)) {
    msg  <- paste0("Multiple adjustments done except for WestFall Young.") # Priamry text we want to display
    updateProgress(message = msg) # Passing back the progress messages onto the callback function
  } # if the function is being called, run the progress bar

  # adjust p-values for Westfall-Young (single-step and step-down)
  order.matrix<-t(apply(abs.Zs.H1,1,order,decreasing=TRUE))
  adjp.SS<-adjust.allsamps.WYSS(snum,abs.Zs.H0,abs.Zs.H1)
  adjp.WY<-adjust.allsamps.WYSD(snum,abs.Zs.H0,abs.Zs.H1,order.matrix,ncl)

  # combine all adjusted p-values in list (each entry is a matrix for given MTP)
  adjp.all<-list(rawp,adjp.BF,adjp.HO,adjp.BH,adjp.SS,adjp.WY)

  # 3rd call back to progress bar on adjustments: Westfall Young Single Step and Step Down is done
  if (is.function(updateProgress) & !is.null(adjp.all)) {
    msg  <- paste0("Multiple adjustments done for WestFall Young SS & SD.") # Priamry text we want to display
    updateProgress(message = msg) # Passing back the progress messages onto the callback function
  } # if the function is being called, run the progress bar

  # for each MTP, get matrix of indicators for whether the adjusted p-value is less than alpha
  reject<-function(x) {as.matrix(1*(x<alpha))}
  reject.all<-lapply(adjp.all,reject)

  # Helper function: In each row for each MTP matrix, count number of p-values less than 0.05,
  # in rows corresponding to false nulls
  lt.alpha<-function(x) {apply(as.matrix(x[,MDES>0]),1,sum)}
  lt.alpha.all<-lapply(reject.all,lt.alpha)

  # indiv power for WY-SS, WY-SD, BH, HO, BF is mean of columns of booleans of whether adjusted pvalues were less than alpha
  # in other words, the null has been rejected
  power.ind.fun<-function(x) {apply(x,2,mean)}
  power.ind.all<-lapply(reject.all,power.ind.fun)
  power.ind.all.mat<-do.call(rbind,power.ind.all)

  # 4th call back to progress bar: Individual power calculations are done
  if (is.function(updateProgress) & !is.null(power.ind.all.mat)) {
    msg  <- paste0("Individual power calculation is done.") # Priamry text we want to display
    updateProgress(message = msg) # Passing back the progress messages onto the callback function
  } # if the function is being called, run the progress bar

  # Helper function: m-min powers for all MTPs (including complete power when m=M)
  power.min.fun <- function(x,M) {
    power.min<-numeric(M)
    cnt<-0
    for (m in 1:M) {
      power.min[m]<-mean(x>cnt)
      cnt<-cnt+1
    }
    return(power.min)
  } # end of calculating d-minimal power

  # calculating d-minimal power
  power.min<-lapply(lt.alpha.all,power.min.fun,M=M)
  power.min.mat<-do.call(rbind,power.min)

  # complete power is the power to detect outcomes at least as large as the MDES on all outcomes
  # separating out complete power from d-minimal power by taking the last entry
  power.cmp<-rep(power.min.mat[1,M],length(power.min)) # should it be numfalse or M?

  # combine all power for all definitions
  all.power.results<-cbind(power.ind.all.mat,power.min.mat[,-M],power.cmp)

  # calculating average individual power
  mean.ind.power <- apply(as.matrix(all.power.results[,1:M][,MDES>0]),1,mean)

  # revise final matrix to report this mean individual power and return results
  all.power.results<-cbind(mean.ind.power,all.power.results)

  # setting the col and row names for all power results table
  colnames(all.power.results)<-c("indiv",paste0("indiv",1:M),paste0("min",1:(M-1)),"complete")
  rownames(all.power.results)<-c("rawp","BF","HO","BH","WY-SS","WY-SD")

  #5th call back to progress bar: All power calculation done
  if (is.function(updateProgress) & !is.null(all.power.results)) {
    msg  <- paste0("All definitions of power calculation are done.") # Priamry text we want to display
    updateProgress(message = msg) # Passing back the progress messages onto the callback function
  } # if the function is being called, run the progress bar

  return(all.power.results)

}

#' Midpoint function
#'
#' Calculating the midpoint between the lower and upper bound by calculating half the distance between the two
#' and adding the lower bound to it. The function is a helper function in determining the MDES that falls within
#' acceptable power range.
#'
#' @param lower lower bound
#' @param upper upper bound
#' @importFrom stats dist
#' @return returns midpoint value

midpoint<-function(lower,upper) {

  lower+(dist(c(lower,upper))/2)
} # midpoint function to calculate the right power

#' MDES (minimum detectable effect size) function
#'
#' The minimum detectable effect size function calculates the most feasible minimum detectable effect size
#' for a given MTP, power and power definition. The goal is to find the MDES value that satisfies the margin of error
#' set in the parameter in the power value.
#'
#'
#' @param M the number of hypothesis tests (outcomes)
#' @param numFalse the number of false nulls. This parameter is used for non-Shiny calculations. For Shiny calculations, refer to Ai_mdes
#' @param J the number of blocks
#' @param n.j the harmonic mean of the number of units per block
#' @param power required statistical power for the experiment
#' @param power.definition definition of statistical power from individual, d-minimal to complete power
#' @param MTP type of multiple testing procedure
#' @param marginError the margin of error for MDES estimation based on targeted power value
#' @param p the proportion of samples that are assigned to the treatment
#' @param alpha the family wise error rate (FWER)
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#' @param numCovar.2 number of Level 2 baseline covariates (set to 0 for this design)
#' @param R2.1 a vector of length M corresponding to R^2 for M outcomes of Level 1 (R^2 = variation in the data explained by the model)
#' @param R2.2 a vector of length M corresponding to R^2 for M outcomes of Level 2 (R^2 = variation in the data explained by the model)
#' @param ICC intraclass correlation
#' @param mod.type "c" for constant effects, "f" for fixed effects, "r" for random effects (parameter not in use at the moment)
#' @param sigma correlation matrix for correlations between test statistics (parameter not in use at the moment.Default is set to 0.99)
#' @param omega NULL (parameter in development)
#' @param tnum the number of test statistics (samples) for all procedures other than Westfall-Young & number of permutations for WY. The default is set at 10,000.
#' @param snum the number of samples for Westfall-Young. The default is set at 1,000.
#' @param Ai_mdes a single entry vector specifying the estimated number of outcomes with a non-zero effect
#' @param updateProgress this is the progress bar function that will be passed to the main MDES calculation function
#' @param rho correlation between outcomes. This generates the sigma matrix.
#' @param ncl ncl the number of clusters to use for parallel processing. The default is set at 2.
#' @importFrom stats qt
#' @return mdes results
#' @export

MDES.blockedRCT.2<-function(M, numFalse,Ai_mdes, J, n.j, power, power.definition, MTP, marginError,
                            p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC,
                            mod.type, sigma = 0, rho = 0.99,omega,
                            tnum = 10000, snum=2, ncl=2, updateProgress=NULL) {

  # Setting up Sigma values
  sigma <- matrix(rho, M, M)
  diag(sigma) <- 1

  # Checks on what we are estimating, sample size
  print(paste("Estimating MDES for target ",power.definition,"power of ",round(power,4)))

  # Check to see if the MTP is Westfall Young and it has enough samples. Otherwise, enforce the requirement.
  if (MTP=="WY-SD" & snum < 1000){
    print("For the step-down Westfall-Young procedure, it is recommended that sample (snum) be at least 1000.")
    snum <- 1000
  } # end of if

  if (MTP!="WY-SD"){
    snum <- 2
  } # end of if

  # Compute Q(m)
  Q.m <- sqrt( (1-R2.1) / (p*(1-p)*J*n.j) )
  t.df <- df(J,n.j,numCovar.1)

  # For raw and BF, compute critical values
  crit.alpha <- qt(p=(1-alpha/2),df=t.df)
  crit.alphaxM <- qt(p=(1-alpha/M/2),df=t.df)

  # Compute raw and BF MDES for individual power
  crit.beta <- ifelse(power > 0.5, qt(power,df=t.df), qt(1-power,df=t.df))
  MDES.raw <- ifelse(power > 0.5, Q.m * (crit.alpha + crit.beta), Q.m * (crit.alpha - crit.beta))
  MDES.BF <- ifelse(power > 0.5, Q.m * (crit.alphaxM + crit.beta), Q.m * (crit.alphaxM - crit.beta))

  # SETTING THE MDES BOUNDS FOR INDIVIDUAL AND OTHER TYPES OF POWER from using raw and bf mdes bounds #

  ### INDIVIDUAL POWER ###
  if (power.definition =="indiv") {

    if (MTP == "raw"){

      # Attaching the MDES result to power results for tabular output
      mdes.results <- t(data.frame(c(MDES.raw,power))) # transpose the MDES raw and power to have the results columnwise
      colnames(mdes.results) <- c("MDES without adjustment", paste0(power.definition, " power"))

      return (mdes.results)

    } #Raw MDES if anybody ever asked for it

    if (MTP == "BF"){

      # Attaching the MDES result to power results for tabular output
      mdes.results <- t(data.frame(c(MDES.BF,power))) #transpose the MDES raw and power to have the results columnwise
      colnames(mdes.results) <- c(paste0( MTP, " adjusted MDES"), paste0(power.definition, " power"))

      return(mdes.results)

    } # Bonferroni adjusted MDES for Individual Power

  } # if we are doing power for raw (i.e. unadjusted) and Bonferroni

  # For individual power, other MDES's will be between MDES.raw and MDES.BF, so make starting value the midpoint!
  if (MTP %in% c("HO","BH","WY-SS","WY-SD") & power.definition == "indiv") {

    lowhigh <- c(MDES.raw,MDES.BF)
    try.MDES <- midpoint(MDES.raw,MDES.BF)

  } # MDES for MTP that is not Bonferroni and for individual power

  ### NOT INDIVIDUAL POWER ###

  # For other scenarios, set low-high intervals and compute midpoint
  # For cases where the power definition is not individual power, restrict it between 0 and 1.
  ifelse (power.definition =="indiv", lowhigh<-c(MDES.raw,1), lowhigh<-c(0,1))

  # Searching for the right MDES through a while loop with 20 iterations as the limit

  try.MDES <- midpoint(lowhigh[1],lowhigh[2]) # Initializing MDES for first attempt
  ii <- 0 # Iteration counter
  target.power <- 0 # Initializing a target power

  # While loop through until the iteration is past 20 or we have met the target.power as we search for the right MDES
  # within the margin of error we have specified.

  while (ii < 20 & (target.power < power - marginError | target.power > power + marginError)) {

    # Passing our callback function
    if (is.function(updateProgress)) {
      text <- paste0("Optiomal MDES is currently in the interval between ",round(lowhigh[1],4)," and ",round(lowhigh[2],4),". ") # Secondary text we want to display
      msg  <- paste0("Trying MDES of ",round(try.MDES,4)," ... ") # Priamry text we want to display
      updateProgress(message = msg, detail = text) # Passing back the progress messages onto the callback function
    } # if the function is being called, run the progress bar

    # Function to calculate the target power to check in with the pre-specified power in the loop
    runpower <- power.blockedRCT.2(M = M, MDES = try.MDES, Ai = Ai_mdes, J = J, n.j = n.j,
                                   p = p, alpha = alpha, numCovar.1 = numCovar.1,numCovar.2=0, R2.1 = R2.1, R2.2 = R2.2, ICC = ICC,
                                   mod.type = mod.type, sigma = sigma, omega = omega,
                                   tnum = tnum, snum = snum, ncl = ncl)

    # Pull out the power value corresponding to the MTP and definition of power
    target.power <- runpower[MTP,power.definition]

    # Displaying the progress of mdes calculation via target power
    if (is.function(updateProgress)) {

      msg <- paste("Estimated power for this MDES is",round(target.power,4)) # Text for estimating power
      updateProgress(message = msg)

    } # checking on Progress Update for MDES

    # If the calculated target.power is within the margin of error of the prescribed power, break and return the results
    if(target.power > power - marginError & target.power < power + marginError){

      mdes.results <- data.frame(try.MDES[1], target.power)

      colnames(mdes.results) <- c(paste0(MTP, " adjusted MDES"),paste0(power.definition, " power")) # Giving the proper colnames

      return(mdes.results)

    } # Return results if our targeted power is within a margin of error of the specified power

    # Check if the calculated target power is greater than the prescribed power
    is.over <- target.power > power

    # if we are overpowered, we can detect EVEN SMALLER effect size so we would shrink the effect range with the
    # high end of the bound being the current MDES. Else it would be the opposite.

    if(!is.over) {
      lowhigh[1] <- try.MDES
    }
    if(is.over) {
      lowhigh[2] <- try.MDES
    }

    # re-establish the midpoint
    try.MDES <- midpoint(lowhigh[1],lowhigh[2])

    # run through another iteration
    ii <- ii + 1

  } # end while

} # MDES blockedRCT 2

#' Calculating Sample for Raw (Unadjusted)
#'
#' This is a Helper function for getting Sample Size when no adjustments has been made to the test statistics.
#' The function starts with PowerUp package function mrss.bira2cl but that function seems to have a bug -
#' Only works if we pass in numeric values and not if we pass in objects that hold those values.
#' Additionally, mrss.bira2cl only computes J, not n.j.
#'
#' @param J the number of blocks
#' @param n.j the harmonic mean of the number of units per block
#' @param J0 starting values for J0 to look for optimal J and n.j
#' @param n.j0 starting values for n.j0 to look for optimal J and n.j
#' @param whichSS which type of sample size to optimize for. J or n.j
#' @param MDES minimum detectable effect size
#' @param power required statistical power for the experiment
#' @param p the proportion of samples that are assigned to the treatment
#' @param alpha the family wise error rate (FWER)
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#' @param numCovar.2 number of Level 2 baseline covariates (set to 0 for this design)
#' @param R2.1 a vector of length M corresponding to R^2 for M outcomes of Level 1 (R^2 = variation in the data explained by the model)
#' @param R2.2 a vector of length M corresponding to R^2 for M outcomes of Level 2 (R^2 = variation in the data explained by the model)
#' @param ICC intraclass correlation
#' @param mod.type "c" for constant effects, "f" for fixed effects, "r" for random effects (parameter not in use at the moment)
#' @param sigma correlation matrix for correlations between test statistics (parameter not in use at the moment.Default is set to 0.99)
#' @param omega NULL (parameter in development)
#' @param two.tailed a boolean value for whether we are looking at two-tailed distribution or not
#' @param num.iter number of iterations to look for sample size. The default is set at 100
#' @param tol tolerance from initializing sample values
#' @return raw sample returns
#' @export

SS.blockedRCT.2.RAW <- function(J, n.j, J0=10, n.j0=10,
                                whichSS, MDES, power, p,
                                alpha, numCovar.1, numCovar.2=0,
                                R2.1, R2.2, ICC, mod.type, sigma,
                                omega, two.tailed = TRUE, num.iter = 100, tol=0.1) {


  i <- 0 # starting the iterator
  conv <- FALSE # boolean value for convergence

  while (i <= num.iter & conv == FALSE) {
    # checking which type of sample we are estimating
    if (whichSS =="J"){
      df <- J0 * (n.j - 1) - numCovar.1 - 1 # degree of freedom calculation
    }
    if (whichSS =="n.j") {
      df <- J * (n.j0 - 1) - numCovar.1 - 1
    }

    if (df < 0 | is.infinite(df)) {
      break
    }

    T1 <- ifelse(two.tailed == TRUE, abs(qt(alpha/2, df)),
                 abs(qt(alpha, df)))

    T2 <- abs(qt(power, df))

    MT <- ifelse(power >= 0.5, T1 + T2, T1 - T2)

    if (whichSS=="J") {

      J1 <- (MT/MDES)^2 * ((1 - R2.1)/(p * (1 - p) * n.j))

      if (abs(J1 - J0) < tol) {

        conv <- TRUE
      }
      J0 <- (J1 + J0)/2
    }
    if (whichSS=="n.j") {
      n.j1 <- (MT/MDES)^2 * ((1 - R2.1)/(p * (1 - p) * J))

      if (abs(n.j1 - n.j0) < tol) {
        conv <- TRUE
      }
      n.j0 <- (n.j1 + n.j0)/2
    }

    i <- i + 1
  }
  if (whichSS=="J") {
    J <- ifelse(df > 0, round(J0), NA)
    return(ceiling(J))
  }
  if (whichSS=="n.j") {
    n.j <- ifelse(df > 0, round(n.j0), NA)
    return(ceiling(n.j))
  }

}

#These currently only work if numFalse = M and if MDES is the same or all outcomes.

#' Sample Function
#'
#' @param M the number of hypothesis tests (outcomes)
#' @param numFalse the number of false nulls. This parameter is used for non-Shiny calculations. For Shiny calculations, refer to Ai_mdes
#' @param typesample the type of the number of sample we would like to estimate: either block J or n.j (harmonic mean within block. For Shiny use)
#' @param J the number of blocks (set to NULL if you do not want to estimate this one)
#' @param n.j the harmonic mean of blocks (set to NULL if you do not want to estimate this one)
#' @param J0 the initial value for the sample number of blocks. The default is set at 10.
#' @param n.j0 the initial value for the harmonic mean for the number of samples within block. The default is set at 10.
#' @param MDES minimum detectable effet size
#' @param power required statistical power for the experiment
#' @param power.definition definition of statistical power from individual, d-minimal to complete power
#' @param MTP type of multiple testing procedure in use from Bonferroni, Benjamini-Hocheberg, Holms, Westfall-Young Single Step, Westfall-Young Step Down
#' @param marginError the margin of error for MDES estimation based on targeted power value
#' @param p the proportion of samples that are assigned to the treatment
#' @param alpha the family wise error rate (FWER)
#' @param numCovar.1 number of Level 1 baseline covariates (not including block dummies)
#' @param numCovar.2 number of Level 2 baseline covariates (set to 0 for this design)
#' @param R2.1 a vector of length M corresponding to R^2 for M outcomes of Level 1 (R^2 = variation in the data explained by the model)
#' @param R2.2 a vector of length M corresponding to R^2 for M outcomes of Level 2 (R^2 = variation in the data explained by the model)
#' @param ICC intraclass correlation
#' @param mod.type "c" for constant effects, "f" for fixed effects, "r" for random effects (parameter not in use at the moment)
#' @param sigma correlation matrix for correlations between test statistics (parameter not in use at the moment.Default is set to 0.99)
#' @param omega NULL (parameter in development)
#' @param tnum the number of test statistics (samples) for all procedures other than Westfall-Young & number of permutations for WY. The default is set at 10,000.
#' @param snum the number of samples for Westfall-Young. The default is set at 1,000.
#' @param ncl ncl the number of clusters to use for parallel processing. The default is set at 2.
#' @param num.iter the number of iterations to look for the optimal sample size. The default is set at 20
#' @param updateProgress a call back function for our internal use in our Shiny application
#' @param rho correlation between outcomes when sigma is generated
#'
#' @return Sample number returns
#' @export

SS.blockedRCT.2 <- function(M, numFalse, typesample, J, n.j,
                            J0, n.j0, MDES, power, power.definition,
                            MTP, marginError,p, alpha, numCovar.1,
                            numCovar.2=0, R2.1, R2.2,ICC,mod.type,
                            sigma = 0, rho = 0.99, omega,tnum = 10000,
                            snum=2, ncl=2, num.iter = 20, updateProgress=NULL) {

  # SET UP #
  sigma <- matrix(rho, M, M)
  diag(sigma) <- 1

  # indicator for which sample to compute. J is for blocks. n.j is for harmonic mean of samples within block
  if(typesample == "J"){

    doJ <- TRUE
    don.j <- FALSE

  } else if (typesample == "n.j") {

    don.j <- TRUE
    doJ <- FALSE

  } # Sample within block

  ifelse(doJ,whichSS<-"J",whichSS<-"n.j")

  # Progress Message for the Type of Sample we are estimating, the type of power and the targeted power value
  if(is.function(updateProgress)){
    msg <- (paste0("Estimating ",whichSS," for target ",power.definition," power of ",round(power,4), " .")) #msg to be displayed in the progress bar
    updateProgress(message = msg)
  } # For printing via update progress function


  # Compute J or n.j for raw and BF SS for INDIVIDUAL POWER. We are estimating bounds like we estimated MDES bounds.
  # for now assuming only two tailed tests
  if (doJ) {
    J.raw <- SS.blockedRCT.2.RAW(J=NULL, n.j, J0=J0, n.j0=n.j0, whichSS, MDES, power, p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC, mod.type, sigma, omega, num.iter = 100, tol=0.1)
    J.BF <- SS.blockedRCT.2.RAW(J=NULL, n.j, J0=J0, n.j0=n.j0, whichSS, MDES, power, p, alpha/M, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC, mod.type, sigma, omega, num.iter = 100, tol=0.1)
  }

  if (don.j) {
    n.j.raw <- SS.blockedRCT.2.RAW(J, n.j=NULL, J0=J0, n.j0=n.j0, whichSS, MDES, power, p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC, mod.type, sigma, omega, num.iter = 100, tol=0.1)
    n.j.BF <- SS.blockedRCT.2.RAW(J, n.j=NULL, J0=J0, n.j0=n.j0, whichSS, MDES, power, p, alpha/M, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC, mod.type, sigma, omega, num.iter = 100, tol=0.1)
  }

  # So below we focus on just one type of sample being estimated: Either block or samples within block
  if (doJ) {
    ss.raw <- J.raw
    ss.BF <- J.BF
  }

  if (don.j) {
    ss.raw <- n.j.raw
    ss.BF <- n.j.BF
  }

  ### INDIVIDUAL POWER for Raw and BF ###
  if (power.definition=="indiv") {

    if (MTP == "raw"){

      # saving the sample estimates for Individual Power with the MTP type
      raw.ss <- data.frame("Raw","Indivdual", ss.raw)
      colnames(raw.ss) <- c("Type of MTP", "Type of Power", "Sample Size")

      # To check if the function is returning a table

      return(raw.ss)

    } #MTP raw

    if (MTP == "BF") {

      # saving the sample estimates for Individual Power with the MTP type
      ss.BF <- data.frame("Bonferroni", "Individual", ss.BF)
      colnames(ss.BF) <- c("Type of MTP", "Type of Power", "Sample Size")

      return(ss.BF)

    } #MTP BF

  } # Individual power

  ### INDIVIDUAL POWER FOR NON BF MTPs ###

  # Like the MDES calculation, the sample size would be between raw and Bonferroni. There is no adjustment and there is very
  # conservative adjustment

  # For individual power, other J's or n.j's will be between raw and BF, so make starting value the midpoint
  if (MTP %in% c("HO","BH","WY-SS","WY-SD") & power.definition == "indiv") {

    lowhigh <- c(ss.raw,ss.BF)
    try.ss <- midpoint(lowhigh[1],lowhigh[2])

  }
  # For minimal powers, makes starting value = raw
  if (power.definition != "indiv")  {

    lowhigh <- c(0,ss.BF)
    try.ss <- midpoint(lowhigh[1],lowhigh[2])

  }
  ii <- 0
  target.power <- 0
  while (ii < num.iter & (target.power < power - marginError | target.power > power + marginError) ) {

    if (is.function(updateProgress)) {

      #if statement to get the more apprpriate sample text for progress bar
      if (whichSS == "J"){
        typeofsample <- "Number of blocks"
      } else {
        typeofsample <- "Samples within blocks"
      }

      # Updating the progress of the search for sample
      text <- paste0(typeofsample, " is in the interval between ",round(lowhigh[1],4)," and ",round(lowhigh[2],4),".")
      msg <- paste0("Trying ",typeofsample," of ",round(try.ss,4), " .")
      updateProgress(message = msg, detail = text)

    }

    # We check the Power of the experimental set up with the given estimated sample size to see if it is
    # within the targeted power

    if (doJ) {

      runpower <- power.blockedRCT.2(M, MDES, Ai = numFalse, J= try.ss,n.j,
                                     p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC,
                                     mod.type, sigma, omega,
                                     tnum, snum, ncl)
    }

    if (don.j) {

      runpower <- power.blockedRCT.2(M, MDES, Ai = numFalse, J, n.j=try.ss,
                                     p, alpha, numCovar.1, numCovar.2=0, R2.1, R2.2, ICC,
                                     mod.type, sigma, omega,
                                     tnum, snum, ncl)

    }

    # Pulling the power result out from the table
    target.power <- runpower[MTP,power.definition]

    # Providing message on current targeted power
    if (is.function(updateProgress)){
      text <- paste0("Estimated power for this ",whichSS," is ",target.power)
      updateProgress(detail = text)
    }

    # checking if the estimation is over or not
    is.over <- target.power > power

    # if the target power is within the margin we have set, we will return the estimated sample
    if(target.power > power - marginError & target.power < power + marginError) {

      # estimated sample for a given MTP, type of power
      try.ss.numeric <- ceiling(as.numeric(try.ss))

      # The estimated sample table with MTP type, Power, Sample Size and the target power
      est.sample <- data.frame(MTP, power.definition, try.ss.numeric, target.power)
      colnames(est.sample) <- c("Type of MTP", "Type of Power", "Sample Size", "Target Power")

      return(est.sample)

    }

    # if the test is underpowered or overpowered, repeat the loop with new bounds as below.

    if(!is.over) {
      p.off <- (power - target.power) / power
      lowhigh[1] <- try.ss
    }

    if(is.over) {
      lowhigh[2] <- try.ss
      p.off <- (target.power - power) / power
    }

    lowhigh.dist <- lowhigh[2]-lowhigh[1]
    try.ss <- ifelse(target.power < power, (try.ss + lowhigh[2])/2, (try.ss + lowhigh[1])/2) # midpoint
    ii <- ii + 1
  } # end while

  if (ii==num.iter & !(target.power > power - marginError & target.power < power + marginError)) {

    text <- paste0("Reached maximum iterations without converging on MDES estimate within margin of error. Try increasing maximum number of iterations (num.iter).")
    updateProgress(detail = text)
  }
} # SS.blockedRCT.2



## indiv, BF, J
# test.SS <- SS.blockedRCT.2(M, numFalse = M, J=NULL, n.j, J0=J0, n.j0=n.j0, MDES = rep(mdes1,M), power=test.power["BF","indiv"], power.definition = "indiv", MTP = "BF", marginError = 0.005,p, alpha, numCovar.1=0, numCovar.2=0, R2.1=r2, R2.2=0, ICC=0, mod.type="constant", sigma=sigma, omega=NULL,  tnum = 10000, snum=2, ncl=4)
# print(test.SS)
# ## indiv, BH, n.j
# test.SS <- SS.blockedRCT.2(M, numFalse = M, J, n.j=NULL, J0=J0, n.j0=n.j0, MDES = mdes1, power=test.power["BH","indiv"], power.definition = "indiv", MTP = "BH", marginError = 0.005,p, alpha, numCovar.1=0, numCovar.2=0, R2.1=r2, R2.2=0, ICC=0, mod.type="constant", sigma=sigma, omega=NULL,  tnum = 10000, snum=2, ncl=4)
# print(test.SS)
# ## min1, BH, J
# test.SS <- SS.blockedRCT.2(M, numFalse = M, J=NULL, n.j, J0=J0, n.j0=n.j0, MDES = mdes1, test.power["BH","min1"], power.definition = "min1", MTP = "BH", marginError = 0.005,p, alpha, numCovar.1=0, numCovar.2=0, R2.1=r2, R2.2=0, ICC=0, mod.type="constant", sigma=sigma, omega=NULL,  tnum = 10000, snum=2, ncl=4)
# print(test.SS)
#
# test.SS <- SS.blockedRCT.2(M, numFalse = M, J=NULL, n.j, power=0.417, power.definition, MTP, marginError = 0.005, p, alpha, numCovar.1=0, numCovar.2=0, R2.1=r2, R2.2=0, ICC=0 mod.type="constant", sigma=sigma, omega=NULL,tnum = 10000, snum=2, ncl=2)
