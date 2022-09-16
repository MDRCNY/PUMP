# ------------------------------#
# check correlation of test statistics
# ------------------------------#


get_rawt <- function(d_m, model.params.list, Tbar, n.sims = 100)
{
    rawt.all <- matrix(NA, nrow = n.sims, ncol = model.params.list$M)
    
    start.time <- Sys.time()
    
    # number of simulations
    for(s in 1:n.sims)
    {
        if (s == 10)
        {
            end.time <- Sys.time()
            iter.time <- difftime(end.time, start.time, 'secs')[[1]]/10
            finish.time <- round((iter.time * n.sims)/60)
            msg <- paste('Estimated time to finish', n.sims,
                         'simulation iterations:',
                         finish.time, 'minutes')
            message(msg)
        }
        if (s %% 100 == 0){ message(paste0("Now processing simulation ", s, " of ", n.sims)) }
        
        sim.data <- gen_sim_data(d_m = d_m, model.params.list, Tbar = Tbar)

        # calculate t statistics
        dat.all <- makelist_samp(sim.data, sim.data$T.x)
        rawpt.out <- get_rawpt(dat.all, d_m = d_m, model.params.list = model.params.list)
        rawt <- sapply(rawpt.out[['rawpt']], function(s){ return(s[['tstat']])})
        rawt.all[s,] <- rawt
    }
    
    return(rawt.all)
}

#' @title Check correlation of test statistics (simulation function)
#' 
#' @description Estimates the pairwise correlations
#' between test statistics for all outcomes.
#' 
#' Takes in two options:
#' - a pumpresult object
#' OR
#' - a list of necessary data-generating parameters
#' - the context (d_m)
#' - Tbar
#' 
#' Note that this function can take several minutes to run.
#'
#' @inheritParams gen_sim_data
#' @param rho.V matrix; correlation matrix of level 3 covariates.
#' @param rho.w0 matrix; correlation matrix of level 3 random effects.
#' @param rho.w1 matrix; correlation matrix of level 3 random impacts.
#' @param rho.X matrix; correlation matrix of level 2 covariates.
#' @param rho.u0 matrix; correlation matrix of level 2 random effects.
#' @param rho.u1 matrix; correlation matrix of level 2 random impacts.
#' @param rho.C matrix; correlation matrix of level 1 covariates.
#' @param rho.r matrix; correlation matrix of level 1 residuals.
#' @param n.sims numeric; Number of simulated datasets to generate.
#' More datasets will achieve a more accurate result
#' but also increase computation time.
#'
#' @return matrix; M x M correlation matrix between test statistics.
#' 
#'
#' @export
#' 
#' @examples
#' pp <- pump_power( d_m = "d3.2_m3ff2rc",
#'                   MTP = "BF",
#'                   MDES = rep( 0.10, 2 ),
#'                   M = 2,
#'                   J = 4, # number of schools/block
#'                   K = 10, # number RA blocks
#'                   nbar = 50,
#'                   Tbar = 0.50, # prop Tx
#'                   alpha = 0.05, # significance level
#'                   numCovar.1 = 5, numCovar.2 = 3,
#'                   R2.1 = 0.1, R2.2 = 0.7,
#'                   ICC.2 = 0.05, ICC.3 = 0.4,
#'                   rho = 0.4, # how correlated test statistics are
#'                   tnum = 200
#' )
#' cor.tstat <- check_cor(
#'     pump.object = pp, n.sims = 4
#' )
#' est.cor <- mean(cor.tstat[lower.tri(cor.tstat)])

check_cor <- function(pump.object = NULL,
                      rho.V = NULL, rho.w0 = NULL, rho.w1 = NULL,
                      rho.X = NULL, rho.u0 = NULL, rho.u1 = NULL,
                      rho.C = NULL, rho.r = NULL,
                      d_m = NULL, model.params.list = NULL, Tbar = 0.5,
                      n.sims = 100)
{
    
    if(is.null(pump.object))
    {
        if(is.null(d_m) | is.null(model.params.list))
        {
            stop("You must provide either a pump object
                or both a d_m string and list of model params.")
        }
    } else
    {
        if(!is.null(d_m) | !is.null(model.params.list))
        {
            stop("You must provide either a pump object
                or both a d_m string and list of model params.")
        }
        model.params.list <- params(pump.object)
        d_m <- d_m(pump.object)
        Tbar <- model.params.list$Tbar
        model.params.list$rho.default <- model.params.list$rho
        # fill in correlations
        if(!is.null(rho.V))
        {
            model.params.list$rho.V <- rho.V
        }
        if(!is.null(rho.w0))
        {
            model.params.list$rho.w0 <- rho.w0
        }
        if(!is.null(rho.w1))
        {
            model.params.list$rho.w1 <- rho.w1
        }
        if(!is.null(rho.X))
        {
            model.params.list$rho.X <- rho.X
        }
        if(!is.null(rho.u0))
        {
            model.params.list$rho.u0 <- rho.u0
        }
        if(!is.null(rho.u1))
        {
            model.params.list$rho.u1 <- rho.u1
        }
        if(!is.null(rho.C))
        {
            model.params.list$rho.C <- rho.C
        }
        if(!is.null(rho.r))
        {
            model.params.list$rho.r <- rho.r
        }
        
    }
    
    rawt.all <- get_rawt(
        d_m = d_m,
        model.params.list = model.params.list,
        Tbar = Tbar, 
        n.sims = n.sims
    )
    
    est.cor.mat <- stats::cor(rawt.all)
    
    return(est.cor.mat)
}
