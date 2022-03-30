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

get_cor <- function(rawt.all)
{
    
    # calculate correlation
    cor.tstat <- stats::cor(rawt.all)
    est.cor <- cor.tstat[lower.tri(cor.tstat)]
    
    return(est.cor)
}

#' @title Check correlation of test statistics (simulation function)
#' 
#' @description Given a design, model, and set of model parameters,
#' estimates the pairwise correlations between test statistics for
#' all outcomes.
#' 
#'
#' @inheritParams gen_sim_data
#' @param n.sims Number of simulated datasets to generate.
#' More datasets will achieve a more accurate result
#' but also increase computation time.
#'
#' @return dataframe; contains estimated pairwise
#' correlations between all outcomes, plus information
#' about input parameters.
#' 
#'
#' @export
check_cor <- function(d_m, model.params.list, Tbar, n.sims = 100)
{
    
    rawt.all <- get_rawt(
        d_m = d_m,
        model.params.list = model.params.list,
        Tbar = Tbar, 
        n.sims = n.sims
    )
    
    est.cor <- get_cor(rawt.all)
    
    return(est.cor)
}
