test_that("simulation function works", {
    
    M <- 3
    rho.default <- 0.5
    default.rho.matrix <- gen_corr_matrix(M = M, rho.scalar = rho.default)
    default.kappa.matrix <- matrix(0, M, M) 
    model.params.list <- list(
        M = 3                             # number of outcomes
        , J = 30                          # number of schools
        , K = 10                          # number of districts
        # (for two-level model, set K = 1)
        , nbar = 50                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        ################################################## grand mean outcome and impact
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = rep(0.125, M)            # minimum detectable effect size      
        ################################################## level 3: districts
        , numCovar.3 = 1                  # number of district covariates
        , R2.3 = rep(0.1, M)              # percent of district variation
        # explained by district covariates
        , rho.V = default.rho.matrix      # MxM correlation matrix of district covariates
        , ICC.3 = rep(0.2, M)             # district intraclass correlation
        , omega.3 = rep(0.1, M)           # ratio of district effect size variability
        # to random effects variability
        , rho.w0 = default.rho.matrix     # MxM matrix of correlations for district random effects
        , rho.w1 = default.rho.matrix     # MxM matrix of correlations for district impacts
        , kappa.w =  default.kappa.matrix # MxM matrix of correlations between district
        # random effects and impacts
        ################################################## level 2: schools
        , numCovar.2 = 1                  # number of school covariates
        , R2.2 = rep(0.1, M)              # percent of school variation
        # explained by school covariates
        , rho.X = default.rho.matrix      # MxM correlation matrix of school covariates
        , ICC.2 = rep(0.2, M)             # school intraclass correlation	
        , omega.2 = rep(0.1, M)           # ratio of school effect size variability
        # to random effects variability
        , rho.u0 = default.rho.matrix     # MxM matrix of correlations for school random effects
        , rho.u1 = default.rho.matrix     # MxM matrix of correlations for school impacts
        , kappa.u = default.kappa.matrix  # MxM matrix of correlations between school
        # random effects and impacts
        ################################################## level 1: individuals
        , numCovar.1 = 1                  # number of individual covariates
        , R2.1 = rep(0.1, M)              # percent of indiv variation explained
        # by indiv covariates
        , rho.C = default.rho.matrix      # MxM correlation matrix of individual covariates
        , rho.r = default.rho.matrix      # MxM matrix of correlations for individual residuals 
    )
    
    sim.data <- gen_sim_data(d_m = 'd3.3',
                             model.params.list = model.params.list,
                             Tbar = 0.5)
    expect_true( !is.null(sim.data))
    expect_equal(8, length(sim.data))
    expect_equal(3, ncol(sim.data$Yobs))
    expect_equal(15000, nrow(sim.data$Yobs))
                 
    pp <- pump_power( d_m = "d3.2_m3ff2rc",
                      MTP = "BF",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 3, # number of schools/block
                      K = 21, # number RA blocks
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5, numCovar.2 = 3,
                      R2.1 = 0.1, R2.2 = 0.7,
                      ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.4, # how correlated outcomes are
                      tnum = 200
    )
    sim.data <- gen_sim_data(pump.object = pp)
    expect_true( !is.null(sim.data))
    expect_equal(8, length(sim.data))
    expect_equal(3, ncol(sim.data$Yobs))
    expect_equal(16254, nrow(sim.data$Yobs))
})