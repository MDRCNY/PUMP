# ------------------------------#
# generate simulation data
# ------------------------------#



process_and_generate_param_list <- function(
        d_m, param.list, pump.object, Tbar = NULL, 
        include_Tx = TRUE ) {
    
    # If first argument is a pump.object, swap arguments
    if ( !is.null(d_m) && (is.pumpresult(d_m) || is.pumpgridresult(d_m)) ) {
        pump.object <- d_m
        d_m <- NULL
    }
    
    if (is.null(pump.object)) {
        if ( include_Tx && (is.null(d_m) || is.null(param.list)) ) {
            stop(paste("You must provide either a pump object or both",
                    "a design string (d_m) and list of model params."))
        }
    } else {
        if ( !is.null(d_m) || !is.null(param.list) ) {
            stop(paste("You must provide either a pump object",
                       "or a design string (d_m) and list of model params",
                       "pair (not both)."))
        }
        
        param.list <- params(pump.object)
        param.list$rho.default <- param.list$rho
        
        d_m <- d_m(pump.object)

        # Convert pump object to param list
        if ( include_Tx ) {
            Tbar <- param.list$Tbar
        }
        
        if ( pump_type(pump.object) == "sample" ) {
            param.list$nbar = 10 
            param.list[ pump.object$Sample.type ] = pump.object$Sample.size
        }
    }
    
    param.list$d_m = d_m
    if ( include_Tx && !is.null( Tbar ) ) {
        param.list$Tbar = Tbar
    }
    return( param.list )
}




#' @title Generate correlation matrix (simulation function)
#'
#' @param M scalar; dimension of matrix.
#' @param rho.scalar scalar; rho value.
#'
#' @return matrix; M x M correlation matrix
#' with rho.scalar as diagonal.
#' 
#' @export
gen_corr_matrix <- function(M, rho.scalar)
{
    rho.matrix <- diag(M) + 
        rho.scalar * matrix(1, M, M) - rho.scalar * diag(M)
    return(rho.matrix)
}


#' generate covariance matrix between two variables
#'
#' @param D dimension of matrix
#' @param var1.vec vector of variances of first variable
#' @param var2.vec vector of variances of second variable
#' @param rho.matrix matrix of correlations
#'
#' @return Sigma matrix of covariance
#' @keywords internal
gen_cov_matrix <- function(D, var1.vec, var2.vec, rho.matrix) {
    Sigma <- matrix(NA, D, D)
    for (k in 1:D) {
        for (l in 1:D) {
            Sigma[k,l] <- rho.matrix[k,l] * 
                sqrt(var1.vec[k]) * sqrt(var2.vec[l])
        }
    }
    return(Sigma)
}


#' generate a parameterized covariance matrix from the provided 3 blocks
#' @param Sigma.w level 3 covariance matrix
#' @param Sigma.z level 2 covariance matrix
#' @param Sigma.wz covariance between level 2 and level 3
#' 
#' @return 2M x 2M matrix for generating correlated pairs of random effects
#' @keywords internal
gen_RE_cov_matrix <- function(Sigma.w, Sigma.z, Sigma.wz) {
    stopifnot( nrow(Sigma.w) == ncol(Sigma.w) )
    stopifnot( all( dim(Sigma.z) == dim( Sigma.wz ) ) )
    stopifnot( all( dim(Sigma.w) == dim( Sigma.wz ) ) )
    
    M <- nrow(Sigma.w)
    
    # full covariance matrix
    Sigma.wz.full                                    <- matrix(NA, 2 * M, 2 * M)
    Sigma.wz.full[1:M, 1:M]                          <- Sigma.w
    Sigma.wz.full[(M + 1):(2 * M), (M + 1):(2 * M)]  <- Sigma.z
    Sigma.wz.full[1:M, (M + 1):(2 * M)]              <- Sigma.wz
    Sigma.wz.full[(M + 1):(2 * M), 1:M]              <- t(Sigma.wz)
    
    Sigma.wz.full
}




#' @title Generate base simulated multi-level data (simulation
#'   function)
#'
#' @description Generates simulated data for multi-level RCTs for
#'   pump-supported designs and models for both unobserved potential
#'   outcomes. This function does not generate treatment assignments
#'   or observed outcomes--see gen_sim_data() for that.
#'
#'   This method takes in a list of necessary data-generating
#'   parameters, following the rest of the package.
#'
#'   This function is beyond the main scope of calculating power, and
#'   is instead used for simulating data. For more info on use, see
#'   the simulation vignette.
#'
#' @inheritParams gen_sim_data
#' @param dgp.params TRUE means param.list is already converted to DGP
#'   parameters, FALSE means it needs to be converted via
#'   `convert_params()`.
#' @seealso gen_sim_data
#'
#' @return list; potential outcomes given control y0, treatment y1,
#'   covariates V.k, X.jk, C.ijk, or list of dataframes if
#'   return.as.dataframe = TRUE.
#'
#' @export
gen_base_sim_data <- function(param.list, pump.object = NULL,
                              return.as.dataframe = TRUE, 
                              no.list = TRUE,
                              dgp.params = FALSE ) {
    
    if ( !dgp.params ) {
        param.list = process_and_generate_param_list( 
            d_m = NULL, param.list = param.list, Tbar = NULL,
            pump.object = pump.object, include_Tx = FALSE 
        )
    
        param.list = convert_params( param.list )   
    }
    
    # ------------------------------#
    # setup: convert model params.list to variables
    # ------------------------------#
    has.level.three <- param.list[['has.level.three']]
    M        <- param.list[['M']];       J       <- param.list[['J']]
    K        <- param.list[['K']];       nbar    <- param.list[['nbar']]
    S.id     <- param.list[['S.id']];    D.id    <- param.list[['D.id']]
    Xi0      <- param.list[['Xi0']];     Xi1     <- param.list[['Xi1']]
    rho.V    <- param.list[['rho.V']];   xi      <- param.list[['xi']]
    eta0.sq  <- param.list[['eta0.sq']]; eta1.sq <- param.list[['eta1.sq']]
    rho.w0   <- param.list[['rho.w0']];  rho.w1  <- param.list[['rho.w1']]
    kappa.w  <- param.list[['kappa.w']]
    rho.X    <- param.list[['rho.X']];   delta   <- param.list[['delta']]
    tau0.sq  <- param.list[['tau0.sq']]; tau1.sq <- param.list[['tau1.sq']]
    rho.u0   <- param.list[['rho.u0']];  rho.u1  <- param.list[['rho.u1']]
    kappa.u  <- param.list[['kappa.u']]
    rho.C    <- param.list[['rho.C']];   gamma   <- param.list[['gamma']]
    rho.r    <- param.list[['rho.r']]
    
    stopifnot( !is.null( nbar ) )
   
    # ------------------------------#
    # Generate school and district IDs
    # ------------------------------#
    # generates vector of school and district IDs, 
    # assuming equal sizes of everything
    if ( is.null( S.id ) ) {
        assignments <- gen_cluster_ids( nbar, J, K )
        # N-length vector of indiv school assignments i.e. (1,1,2,2,3,3)
        S.id        <- assignments[['S.id']] 
        # N-length vector of indiv district assignments i.e. (1,1,1,2,2,2)
        D.id        <- assignments[['D.id']]  
    }
    N <- length( S.id )
    
    # ------------------------------#
    # Districts: Level 3
    # ------------------------------#
    
    gen.level.three.data <- function(K, M, rho.V, 
                                     eta0.sq, eta1.sq, 
                                     rho.w0, rho.w1, kappa.w)
    {
        # generate district covariates
        Sigma.V  <- gen_cov_matrix(M, rep(1, M), rep(1, M), rho.V)
        V.k      <- matrix(mvtnorm::rmvnorm(K, mean = rep(0, M), 
                                            sigma = Sigma.V), K, M)
        
        # covariance of random district effects
        Sigma.w0       <- gen_cov_matrix(M, eta0.sq, eta0.sq, rho.w0)
        # covariance random district impacts
        Sigma.w1       <- gen_cov_matrix(M, eta1.sq, eta1.sq, rho.w1)
        # covariance between impacts and effects
        Sigma.w        <- gen_cov_matrix(M, eta0.sq, eta1.sq, kappa.w)
        # full covariance matrix
        Sigma.w.full   <- gen_RE_cov_matrix( Sigma.w0, Sigma.w1, Sigma.w )
        
        # generate full vector of district random effects and impacts
        w01.k <- matrix(mvtnorm::rmvnorm(K, mean = rep(0, 2 * M), 
                                         sigma = Sigma.w.full), 
                        nrow = K, ncol = 2 * M)
        w0.k  <- w01.k[,1:M, drop = FALSE]
        w1.k  <- w01.k[,(M + 1):(2 * M), drop = FALSE]
        return(list(V.k = V.k, w0.k = w0.k, w1.k = w1.k))
    }
    
    if ( has.level.three ) {
        stopifnot( !is.null( K ) )
        level.three.data <- gen.level.three.data(K, M, rho.V, 
                                                 eta0.sq, eta1.sq, 
                                                 rho.w0, rho.w1, kappa.w) 
        V.k <- level.three.data[['V.k']]
        w0.k <- level.three.data[['w0.k']]
        w1.k <- level.three.data[['w1.k']]
    } else {
        V.k <- NULL
        w0.k <- NULL
        w1.k <- NULL
    }
    
    # ------------------------------#
    # Schools: Level 2
    # ------------------------------#
    
    gen.level.two.data <- function(J, K, M, rho.X, 
                                   tau0.sq, tau1.sq, 
                                   rho.u0, rho.u1, kappa.u)
    {
        # generate school covariates
        Sigma.X       <- gen_cov_matrix(M, rep(1, M), rep(1, M), rho.X)
        X.jk          <- matrix(
            mvtnorm::rmvnorm(J*K, mean = rep(0, M), sigma = Sigma.X),
            nrow = J*K, ncol = M
        )
        
        # covariance of school random effects
        Sigma.u0       <- gen_cov_matrix(M, tau0.sq, tau0.sq, rho.u0)
        # covariance of school random impacts
        Sigma.u1       <- gen_cov_matrix(M, tau1.sq, tau1.sq, rho.u1)
        # covariance of school random effects and impacts
        Sigma.u        <- gen_cov_matrix(M, tau0.sq, tau1.sq, kappa.u)
        # full covariance matrix
        Sigma.u.full   <- gen_RE_cov_matrix( Sigma.u0, Sigma.u1, Sigma.u )
        
        # generate full vector of school random effects and impacts
        u01.jk <- matrix(
            mvtnorm::rmvnorm(J*K, mean = rep(0, 2*M), sigma = Sigma.u.full),
            nrow = J*K, ncol = 2*M
        )
        u0.jk  <- u01.jk[,1:M, drop = FALSE]
        u1.jk  <- u01.jk[,(M + 1):(2 * M), drop = FALSE]
        
        return(list(X.jk = X.jk, u0.jk = u0.jk, u1.jk = u1.jk))
    }
    
    level.two.data <- gen.level.two.data(J, K, M, rho.X, 
                                         tau0.sq, tau1.sq, 
                                         rho.u0, rho.u1, kappa.u) 
    X.jk  <- level.two.data[['X.jk']]
    u0.jk <- level.two.data[['u0.jk']]
    u1.jk <- level.two.data[['u1.jk']]
    
    # ------------------------------#
    # Individuals: Level 1
    # ------------------------------#
    
    gen.level.one.data <- function(N, M, rho.C, rho.r)
    {
        # generate individual covariates
        Sigma.C <- gen_cov_matrix(M, rep(1, M), rep(1, M), rho.C)
        C.ijk   <- matrix(mvtnorm::rmvnorm(N, mean = rep(0, M), 
                                           sigma = Sigma.C), N, M)
        
        # generate individual residuals
        Sigma.r <- gen_cov_matrix(M, rep(1, M), rep(1, M), rho.r)
        r.ijk   <- matrix(mvtnorm::rmvnorm(N, mean = rep(0,M), 
                                           sigma = Sigma.r), N, M)
        
        return(list(C.ijk = C.ijk, r.ijk = r.ijk))
    }
    
    level.one.data <- gen.level.one.data(N, M, rho.C, rho.r) 
    C.ijk <- level.one.data[['C.ijk']]
    r.ijk <- level.one.data[['r.ijk']]
    
    # ------------------------------#
    # reformat everything into N x M matrices
    # ------------------------------#
    # for example, D is K x M, now I populate V.k, which is N x M,
    # by filling in district information for each individual
    
    V.ijk  <- w0.ijk <- w1.ijk <-
        X.ijk  <- u0.ijk <- u1.ijk <- matrix(NA, N, M)
    
    Xi0.ijk <- matrix(Xi0, nrow = N, ncol = M, byrow = TRUE) 
    Xi1.ijk <- matrix(Xi1, nrow = N, ncol = M, byrow = TRUE) 
    
    # loop through each individual student
    for (i in 1:N)
    {
        if ( has.level.three )
        {
            # fill in values from district level variables
            V.ijk[i,]    <- V.k[D.id[i],]
            w0.ijk[i,]   <- w0.k[D.id[i],]
            w1.ijk[i,]   <- w1.k[D.id[i],]
        }
        
        # fill in values from school level variables
        X.ijk[i,]    <- X.jk[S.id[i],]
        u0.ijk[i,]   <- u0.jk[S.id[i],]
        u1.ijk[i,]   <- u1.jk[S.id[i],]
    }
    
    # ------------------------------#
    # generate potential outcomes
    # ------------------------------#
    
    # district level
    if ( has.level.three ) {
        mu0.ijk <- Xi0.ijk  + xi * V.ijk + w0.ijk
        mu1.ijk <- Xi1.ijk               + w1.ijk
    } else {
        mu0.ijk <- Xi0.ijk
        mu1.ijk <- Xi1.ijk
    }
    # school level
    theta0.ijk <- mu0.ijk + delta * X.ijk + u0.ijk
    
    # treatment impact
    psi1.ijk   <- mu1.ijk                + u1.ijk
    
    # individual level
    Y0.ijk     <- theta0.ijk  + gamma * C.ijk + r.ijk
    Y1.ijk     <- Y0.ijk                      + psi1.ijk
    
    colnames(Y0.ijk) <- colnames(Y1.ijk) <- paste0('m', 1:M)
    Y0.ijk <- data.frame(Y0.ijk)
    Y1.ijk <- data.frame(Y1.ijk)
    
    ID <- data.frame( S.id = S.id, D.id = D.id )
    
    res = NA
    if ( has.level.three ) {
        res <- list(Y0 = Y0.ijk, Y1 = Y1.ijk, 
                    V.k = V.ijk, X.jk = X.ijk, 
                    C.ijk = C.ijk, ID = ID )
    } else {
        ID$D.id <- NULL
        res <- list(Y0 = Y0.ijk, Y1 = Y1.ijk, 
                    V.k = NULL, X.jk = X.ijk, 
                    C.ijk = C.ijk, ID = ID )
    }
    
    if ( return.as.dataframe ) {
        res <- makelist_samp(res)
        if ( no.list && length(res) == 1 ) {
            return( res[[1]] )
        } else {
            return( res )
        }
    } else {
        return( res )
    }
}





#' @title Generate simulated multi-level data (simulation function)
#'
#' @description Generates simulated data for multi-level RCTs for
#'   pump-suppored designs and models for both unobserved and observed
#'   potential outcomes.
#'
#'   Takes in two options: - a pumpresult object OR - a list of
#'   necessary data-generating parameters - the context (d_m) - Tbar
#'   (proportion assigned to treatment)
#'
#'   This function is beyond the main scope of calculating power, and
#'   is instead used for simulating data. For more info on use, see
#'   the simulation vignette.
#'
#' @param pump.object A pumpresult object.
#' @param d_m string; a single context, which is a design and model
#'   code. See pump_info() for list of choices.
#' @param param.list list; model parameters such as ICC, R2, etc. See
#'   simulation vignette for details.
#' @param Tbar scalar; the proportion of samples that are assigned to
#'   the treatment.
#' @param return.as.dataframe TRUE means return list of dataframes,
#'   one for each outcome.  FALSE means return components of the
#'   covariates, etc., in a list.
#' @param no.list Only relevant if return.as.dataframe=TRUE.
#'   no.list=TRUE means if M=1 return the dataframe, not a list of
#'   length 1.  FALSE means return a list of length 1, even if there
#'   is only 1 outcome.
#' @return list; potential outcomes, covariates, observed outcomes,
#'   and treatment assignment.
#'
#'
#' @export
#' @examples
#'
#' pp <- pump_power( d_m = "d3.2_m3ff2rc",
#'                   MTP = "BF",
#'                   MDES = rep( 0.10, 3 ),
#'                   M = 3,
#'                   J = 3, # number of schools/block
#'                   K = 21, # number RA blocks
#'                   nbar = 258,
#'                   Tbar = 0.50, # prop Tx
#'                   alpha = 0.05, # significance level
#'                   numCovar.1 = 5, numCovar.2 = 3,
#'                   R2.1 = 0.1, R2.2 = 0.7,
#'                   ICC.2 = 0.05, ICC.3 = 0.4,
#'                   rho = 0.4,
#'                   tnum = 200
#' )
#' sim.data <- gen_sim_data(pump.object = pp)
#' 
gen_sim_data <- function(
        d_m = NULL, param.list = NULL, Tbar = 0.5,
        pump.object = NULL,
        return.as.dataframe = TRUE,
        no.list = TRUE ) {
    
    param.list = process_and_generate_param_list( d_m = d_m,
                                                  param.list = param.list,
                                                  Tbar = Tbar,
                                                  pump.object = pump.object )
    d_m = param.list$d_m
    Tbar = param.list$Tbar
    
    sim.data <- gen_base_sim_data( param.list, return.as.dataframe = FALSE )
    sim.data$T.x <- gen_T.x(
        d_m = d_m,
        S.id = sim.data$ID$S.id,
        D.id = sim.data$ID$D.id,
        Tbar = Tbar
    )
    sim.data$Yobs <- gen_Yobs(sim.data, T.x = sim.data$T.x)
    
    if ( return.as.dataframe ) {
        res <- makelist_samp(sim.data)
        if ( no.list && length(res) == 1 ) {
            return( res[[1]] )
        } else {
            return( res )
        }
    } else {
        return(sim.data)
    }
}



#' @title Converts model params into DGP params (simulation function)
#'
#' @description Converts user-provided parameters such as ICC and
#'   omega into data-generating parameters for the multilevel random
#'   effects model used to produce simulated data, such as variance
#'   values and covariate coefficients.
#'
#'   This function is beyond the main scope of calculating power, and
#'   is instead used for simulating data. For more info on use, see
#'   the simulation vignette.
#'
#' @param param.list list; model parameters such as ICC, R2, etc.
#'
#' @return list; data-generating parameters.
#'
#' @export
convert_params <- function(param.list) {
    
    # save out useful parameters
    M <- param.list[['M']]
    ICC.2 <- param.list[['ICC.2']]
    ICC.3 <- param.list[['ICC.3']]
    R2.1 <- param.list[['R2.1']]
    R2.2 <- param.list[['R2.2']]
    R2.3 <- param.list[['R2.3']]
    omega.2 <- param.list[['omega.2']]
    omega.3 <- param.list[['omega.3']]
    rho.default <- param.list[['rho.default']]
    
    # If only single outcome, rho not needed.  Give default so code runs.
    if ( M == 1 && is.null( rho.default ) ) {
        rho.default = 0
    }
    
    # If no district info, set district parameters to 0
    has.level.three <- TRUE
    if ( is.null( param.list$K ) ) {
        has.level.three <- FALSE
        ICC.3 <- param.list$ICC.3 <- rep(0, M)
        R2.3 <- param.list$R2.3 <- rep(0, M)
        omega.3 <- param.list$omega.3 <- rep(0, M)
        K <- param.list$K <- 1
    }
    
    if (is.null(rho.default) & any(c(
        is.null(param.list$rho.V), 
        is.null(param.list$rho.w0),
        is.null(param.list$rho.w1),
        is.null(param.list$rho.X), 
        is.null(param.list$rho.u0),
        is.null(param.list$rho.u1),
        is.null(param.list$rho.C), 
        is.null(param.list$rho.r))))
    {
        stop(paste('Please provide either a rho.default',
                   'or ALL necessary correlation matrices.'))
    }
    
    default.rho.matrix <- gen_corr_matrix(M = M, rho.scalar = rho.default)
    
    if (is.null(param.list$rho.V))
    {
        param.list$rho.V  <- default.rho.matrix 
    }
    if (is.null(param.list$rho.w0))
    {
        param.list$rho.w0  <- default.rho.matrix 
    }
    if (is.null(param.list$rho.w1))
    {
        param.list$rho.w1 <- default.rho.matrix 
    }
    
    if (is.null(param.list$rho.X))
    {
        param.list$rho.X  <- default.rho.matrix 
    }
    if (is.null(param.list$rho.u0))
    {
        param.list$rho.u0  <- default.rho.matrix 
    }
    if (is.null(param.list$rho.u1))
    {
        param.list$rho.u1 <- default.rho.matrix 
    }
    
    if (is.null(param.list$rho.C))
    {
        param.list$rho.C  <- default.rho.matrix 
    }
    if (is.null(param.list$rho.r))
    {
        param.list$rho.r  <- default.rho.matrix 
    }
    
    convert.scalar <- function(x, M)
    {
        if ( is.null(x) )
        {
            return( rep(0, M) )
        } else if ( length(x) == 1 )
        {
            return( rep(x, M) )
        } else
        {
            return(x)
        }
    }
    param.list$MDES <- convert.scalar(param.list$MDES, M)
    R2.1    <- convert.scalar(R2.1, M)
    R2.2    <- convert.scalar(R2.2, M)
    R2.3    <- convert.scalar(R2.3, M)
    ICC.2   <- convert.scalar(ICC.2, M)
    ICC.3   <- convert.scalar(ICC.3, M)
    omega.2 <- convert.scalar(omega.2, M)
    omega.3 <- convert.scalar(omega.3, M)
    
    # defaults for kappa matrices
    if (is.null(param.list[['kappa.w']]))
    {
        param.list$kappa.w <- matrix(0, M, M) 
    }
    if (is.null(param.list[['kappa.u']]))
    {
        param.list$kappa.u <- matrix(0, M, M) 
    }
    
    # default for grand mean
    if (is.null(param.list[['Xi0']]))
    {
        param.list$Xi0 <- 0
    }
    
    # check ICC is valid
    if ( any( ICC.2 + ICC.3 >= 1 ) )
    {
        stop(paste('ICC.2 + ICC.3 must be less than 1.
                   ICC.2:', ICC.2, 'ICC3:', ICC.3))
    }
    
    # random intercepts variances
    tau0.sq <- ( (1 - R2.2) / (1 - R2.1) ) * ( ICC.2 / (1 - ICC.2 - ICC.3) )
    eta0.sq <- ( (1 - R2.3) / (1 - R2.1) ) * ( ICC.3 / (1 - ICC.2 - ICC.3) )
    # covariate coefficients
    delta   <- sqrt( (R2.2 / (1 - R2.1)) *  ( ICC.2 / (1 - ICC.2 - ICC.3) ) )
    xi      <- sqrt( (R2.3 / (1 - R2.1)) *  ( ICC.3 / (1 - ICC.2 - ICC.3) ) )
    gamma   <- sqrt( (R2.1 / (1 - R2.1)) )
    
    # random impacts variances
    tau1.sq <- omega.2 * (tau0.sq + delta^2)
    eta1.sq <- omega.3 * (eta0.sq + xi^2)
    
    # grand mean impact
    Xi1 <- param.list[['MDES']] *
        sqrt(xi^2 + gamma^2 + delta^2 + eta0.sq + tau0.sq + 1)
    
    new.param.list <- list(
        has.level.three = has.level.three
        , M = param.list[['M']]                    # number of outcomes
        , J = param.list[['J']]                    # number of schools
        , K = param.list[['K']]                    # number of districts
        , nbar = param.list[['nbar']]              # number of individuals per school
        , Xi0 = param.list[['Xi0']]                # scalar grand mean outcome under no treatment
        , Xi1 = Xi1                                # scalar grand mean impact
    )
    
    if ( has.level.three ) {
        new.param.list <- c( new.param.list, list( 
            # -------------------------------------------- level 3
            xi = xi                             # M-vector of coefficient of district covariates
            , rho.V = param.list[['rho.V']]     # MxM correlation matrix of district covariates
            , eta0.sq = eta0.sq                 # M-vector of variances of district random effects
            , eta1.sq = eta1.sq                 # M-vector of variances of district impacts
            , rho.w0 = param.list[['rho.w0']]   # MxM matrix of correlations for district random effects
            , rho.w1 = param.list[['rho.w1']]   # MxM matrix of correlations for district impacts
            , kappa.w = param.list[['kappa.w']] # MxM matrix of correlations between district random effects and impacts
        ) )
    }
    
    new.param.list <- c( new.param.list, list(
        # -------------------------------------------- level 2
        delta = delta                             # M-vector of coefficients of school covariates
        , rho.X = param.list[['rho.X']]           # MxM correlation matrix of school covariates
        , tau0.sq = tau0.sq                       # M-vector of variances of school random effects
        , tau1.sq = tau1.sq                       # M-vector of variances of school impacts
        , rho.u0 = param.list[['rho.u0']]         # MxM matrix of correlations for school random effects
        , rho.u1 = param.list[['rho.u1']]         # MxM matrix of correlations for school impacts
        , kappa.u = param.list[['kappa.u']]       # MxM matrix of correlations between school random effects and impacts
        # -------------------------------------------- level 1
        , gamma = gamma                           # M-vector of coefficients of individual covariates
        , rho.C = param.list[['rho.C']]           # MxM correlation matrix of individual covariates
        , rho.r = param.list[['rho.r']]           # MxM matrix of correlations for individual residuals
    ) )
    
    nulls <- vapply( new.param.list, is.null, logical(1) )
    new.param.list <- new.param.list[ !nulls ]
    
    return(new.param.list)
}


#' @title Generates school and district assignments (simulation
#'   function)
#'
#' @description Generates simple default schools and districts IDs for
#'   individual students for the purpose of simulations. This assumes
#'   equal sized schools in equal sized districts.
#'
#'   This function is beyond the main scope of calculating power, and
#'   is instead used for simulating data. For more info on use, see
#'   the simulation vignette.
#'
#' @param K scalar; number of districts.
#' @param J scalar; number of schools per district.
#' @param nbar scalar; number of individuals per school.
#'
#' @return list; school and district assignments (S.id, D.id) for each
#'   individual.
#'
#' @export
gen_cluster_ids <- function(nbar, J, K){
    
    J <- ifelse(is.null(J), 1, J)
    K <- ifelse(is.null(K), 1, K)
    N <- nbar * J * K
    
    # vector of assignments to schools
    S.id <- rep(NA, N)
    start.index <- 1
    end.index <- nbar
    for (j in 1:(K*J))
    {
        S.id[start.index:end.index] <- j
        start.index <- end.index + 1
        end.index <- end.index + nbar
    }
    
    D.id <- rep(NA, N)
    start.index <- 1
    n.k <- N/K
    end.index <- n.k
    for (k in 1:K)
    {
        D.id[start.index:end.index] <- k
        start.index <- end.index + 1
        end.index <- end.index + n.k
    }
    stopifnot( all( !is.na( S.id ) ) )
    stopifnot( all( !is.na( D.id ) ) )
    
    return(list(S.id = S.id, D.id = D.id))
}




#' @title Generate treatment assignment vector (simulation function)
#' 
#' @description Given a RCT design and supporting information,
#' generates treatment assignments for each student.
#' 
#' This function is beyond the main scope of calculating power,
#' and is instead used for simulating data.
#' For more info on use, see the simulation vignette.
#'
#' @param d_m string; design and model.
#' @param S.id vector; school assignments.
#' @param D.id vector; district assignments.
#' @param Tbar scalar; probability of treatment assignment.
#'
#' @return vector; treatment assignments for each unit.
#' 
#' @export
gen_T.x <- function(d_m, S.id, D.id, Tbar)
{
    od <- d_m
    d_m <- strsplit(d_m, "_")
    d_m <- d_m[[1]][[1]]
    
    if (d_m == "d1.1") {
        T.x <- randomizr::simple_ra(N = length(S.id), prob = Tbar)
    } else if ( d_m == "d2.1" || d_m == "d3.1" ) {
        T.x <- randomizr::block_ra( S.id, prob = Tbar )
    } else if ( d_m == "d2.2" ) { 
        T.x <- randomizr::cluster_ra( S.id, prob = Tbar )
    } else if ( d_m == "d3.3" ) {
        T.x <- randomizr::cluster_ra( D.id, prob = Tbar )
    } else if ( d_m == "d3.2" ) {
        T.x <- randomizr::block_and_cluster_ra( 
            blocks = D.id, clusters = S.id, prob = Tbar 
        )
    } else {
        stop(print(paste('design', d_m, 'not implemented yet')))
    }
    return(T.x)
}



#' @title Generate observed outcomes (simulation function)
#'
#' @description Takes in a full dataset of both observed and latent
#'   potential outcomes and the treatment assignment vector, and
#'   returns only the observed outcomes.
#'
#'   This function is beyond the main scope of calculating power, and
#'   is instead used for simulating data. For more info on use, see
#'   the simulation vignette.
#'
#' @param full.data data.frame; full dataset of potential outcomes.
#' @param T.x vector; binary assignment to treat/control.
#'
#' @return vector; observed outcomes
#'
#' @export
gen_Yobs <- function(full.data, T.x) {
    Yobs <- full.data$Y0
    Yobs[T.x == 1,] <- full.data$Y1[T.x == 1,]
    return(Yobs)
}


#' Convert multi-outcome data structure to list for each outcome.
#'
#' Given the simulated multi-outcome structure, make a list of
#' complete (tidy) rectangular datasets, one for each outcome.
#'
#' @param samp.obs a single iteration of observed data
#' @param T.x vector of treatment assignments
#' @return List of dataframes.
#' @keywords internal
makelist_samp <- function(samp.obs, T.x = NULL ) {
    
    if ( is.null( T.x ) ) {
        T.x = samp.obs$T.x
    }
    
    M = ncol( samp.obs$C.ijk )
    
    tx_assigned = TRUE
    if ( is.null( samp.obs[['Yobs']] ) ) {
        stopifnot( is.null( T.x ) )
        tx_assigned = FALSE
    } else {
        stopifnot( !is.null( T.x ) )
    }
    
    mdat.rn <- rep( NULL, M )
    for (m in 1:M)
    {
        # level 3
        if (!is.null(samp.obs[['V.k']]))
        {
            mdat.rn[[m]] <- data.frame(
                V.k         = samp.obs[['V.k']][,m],
                X.jk        = samp.obs[['X.jk']][,m],
                C.ijk       = samp.obs[['C.ijk']][,m],
                S.id        = as.factor(samp.obs$ID$S.id),
                D.id        = as.factor(samp.obs$ID$D.id)
            )
        } else
            # level 2
        {
            mdat.rn[[m]] <- data.frame(
                X.jk        = samp.obs[['X.jk']][,m],
                C.ijk       = samp.obs[['C.ijk']][,m],
                S.id        = as.factor(samp.obs$ID$S.id)
            )
        }
        
        if ( tx_assigned ) {
            mdat.rn[[m]]$Yobs = samp.obs[['Yobs']][,m]
            mdat.rn[[m]]$T.x         = T.x
        }
    }
    return(mdat.rn)
}



